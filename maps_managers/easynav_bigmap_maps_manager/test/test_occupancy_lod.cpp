// Copyright 2026 José Miguel Guerrero
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// The occupancy block of bigmap (occupancy.*): the default LOD, but the size cap of the
// distance applies only where something is occupied (a hole's border, a sensed obstacle):
//   * meshes without T-junctions, whatever the holes, obstacles and robot;
//   * within the full-resolution radius, every NavCel over something occupied is of the minimum
//     size, and free ground is larger (as large as its height error allows);
//   * beyond it, what is occupied grows with the distance (as in the distance mode);
//   * the height error still splits free ground;
//   * the incremental update of the sensed obstacles gives the same mesh as building it again.

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <utility>
#include <vector>

#include "easynav_bigmap_maps_manager/known_map.hpp"
#include "easynav_bigmap_maps_manager/lod_mesher.hpp"
#include "easynav_bigmap_maps_manager/lod_policy.hpp"
#include "easynav_bigmap_maps_manager/occupancy_lod.hpp"

using namespace bigmap_maps_manager;  // NOLINT(build/namespaces)

namespace
{

DemGrid flat_dem(int side, double res)
{
  DemGrid g;
  g.side = side;
  g.resolution = res;
  g.radius = (side - 1) / 2 * res;
  g.z.assign(static_cast<std::size_t>(side) * side, 0.0f);
  return g;
}

// Edges used by one triangle only with a mesh vertex strictly inside: a T-junction
std::size_t t_junctions(const LodMesh & m, double res)
{
  using Edge = std::pair<uint32_t, uint32_t>;
  std::map<Edge, int> use;
  for (std::size_t t = 0; t < m.triangle_count(); ++t) {
    for (int k = 0; k < 3; ++k) {
      const uint32_t a = m.tri[3 * t + k], b = m.tri[3 * t + (k + 1) % 3];
      ++use[{std::min(a, b), std::max(a, b)}];
    }
  }
  std::set<std::pair<long, long>> nodes;   // NOLINT(runtime/int)
  auto node = [&](uint32_t v) {
      return std::make_pair(std::lround(m.x[v] / res), std::lround(m.y[v] / res));
    };
  for (uint32_t v = 0; v < m.x.size(); ++v) {
    nodes.insert(node(v));
  }
  std::size_t bad = 0;
  for (const auto & [e, n] : use) {
    if (n != 1) {continue;}
    const auto a = node(e.first), b = node(e.second);
    const long dx = b.first - a.first, dy = b.second - a.second;   // NOLINT(runtime/int)
    const long g = std::gcd(std::labs(dx), std::labs(dy));          // NOLINT(runtime/int)
    for (long k = 1; k < g; ++k) {                                  // NOLINT(runtime/int)
      if (nodes.count({a.first + dx / g * k, a.second + dy / g * k})) {++bad;}
    }
  }
  return bad;
}

// The kinds of the cells under a triangle (0 ground, 1 sensed obstacle, 2 no ground), sampled at a
// quarter of a cell inside it
std::set<int> kinds_under(
  const LodMesh & m, std::size_t t, const DemGrid & dem, const OccupancyLod & occ)
{
  const int c = (dem.side - 1) / 2;
  const double res = dem.resolution;
  double px[3], py[3];
  for (int k = 0; k < 3; ++k) {
    px[k] = m.x[m.tri[3 * t + k]] / res + c;
    py[k] = m.y[m.tri[3 * t + k]] / res + c;
  }
  const double x0 = std::min({px[0], px[1], px[2]}), x1 = std::max({px[0], px[1], px[2]});
  const double y0 = std::min({py[0], py[1], py[2]}), y1 = std::max({py[0], py[1], py[2]});
  auto side = [](double ax, double ay, double bx, double by, double x, double y) {
      return (bx - ax) * (y - ay) - (by - ay) * (x - ax);
    };
  std::set<int> out;
  for (double y = y0 + 0.125; y < y1; y += 0.25) {
    for (double x = x0 + 0.125; x < x1; x += 0.25) {
      const double s0 = side(px[0], py[0], px[1], py[1], x, y);
      const double s1 = side(px[1], py[1], px[2], py[2], x, y);
      const double s2 = side(px[2], py[2], px[0], py[0], x, y);
      if ((s0 >= 0 && s1 >= 0 && s2 >= 0) || (s0 <= 0 && s1 <= 0 && s2 <= 0)) {
        out.insert(occ.cell_kind(static_cast<int>(std::floor(x)), static_cast<int>(std::floor(y))));
      }
    }
  }
  return out;
}

bool occupied(const std::set<int> & kinds)
{
  return kinds.count(1) || (kinds.count(0) && kinds.count(2));
}

// distance from the robot to the farthest vertex of a triangle
double farthest(const LodMesh & m, std::size_t t, double rx, double ry)
{
  double far = 0.0;
  for (int k = 0; k < 3; ++k) {
    const auto v = m.tri[3 * t + k];
    far = std::max(far, static_cast<double>(std::hypot(m.x[v] - rx, m.y[v] - ry)));
  }
  return far;
}

LodParams params(double res)
{
  LodParams p;
  p.min_size = res;
  p.max_size = 32 * res;
  p.high_res_radius = 12 * res;
  p.low_res_distance = 60 * res;
  return p;
}

struct World
{
  DemGrid dem;
  std::vector<std::uint8_t> cell_valid;
};

// Rooms: walls of holes (`thick` cells thick) with doors of 4 cells, random pillars of one cell
// (none if the walls are thick: everything there is thick)
World rooms(int side, double res, std::mt19937 & rng, int thick = 1)
{
  World w{flat_dem(side, res), {}};
  const int nc = side - 1;
  w.cell_valid.assign(static_cast<std::size_t>(nc) * nc, 1);
  auto hole = [&](int i, int j) {
      if (i >= 0 && j >= 0 && i < nc && j < nc) {
        w.cell_valid[static_cast<std::size_t>(j) * nc + i] = 0;
      }
    };
  for (int k = nc / 4; k < nc; k += nc / 4) {
    for (int t = 0; t < nc; ++t) {
      if (t % (nc / 4) >= 3 && t % (nc / 4) < 7) {continue;}   // doors of 4 cells
      for (int e = 0; e < thick; ++e) {
        hole(k + e, t);
        hole(t, k + e);
      }
    }
  }
  std::uniform_int_distribution<int> u(0, nc - 1);
  for (int n = 0; thick == 1 && n < nc / 4; ++n) {
    hole(u(rng), u(rng));
  }
  return w;
}

std::vector<std::pair<int, int>> random_obstacles(int side, std::mt19937 & rng, int n)
{
  std::uniform_int_distribution<int> u(0, side - 2);
  std::vector<std::pair<int, int>> out;
  for (int k = 0; k < n; ++k) {
    const int i = u(rng), j = u(rng);
    for (int a = 0; a < 3; ++a) {
      out.push_back({i + a, j});
    }
  }
  return out;
}

}  // namespace

TEST(OccupancyLod, CrackFreeWithHolesAndSensedObstacles)
{
  std::mt19937 rng(20261116);
  for (int trial = 0; trial < 8; ++trial) {
    const int side = (trial % 2) ? 129 : 97;
    const double res = (trial % 3 == 0) ? 1.0 : 0.5;
    auto w = rooms(side, res, rng);
    LodMesher mesher(w.dem, {}, {}, w.cell_valid);
    OccupancyLod occ(mesher.dem(), mesher.error_map().tile_size(), mesher.cell_valid());
    occ.set_obstacles(random_obstacles(side, rng, 20));
    const auto p = params(res);
    std::uniform_real_distribution<double> u(-(side / 2) * res, (side / 2) * res);
    for (int r = 0; r < 4; ++r) {
      const auto mesh = mesher.extract(u(rng), u(rng), p, nullptr, nullptr, &occ);
      EXPECT_EQ(t_junctions(mesh, res), 0u) << "trial " << trial;
    }
  }
}

TEST(OccupancyLod, WithinTheRadiusOnlyWhatIsOccupiedIsFine)
{
  std::mt19937 rng(11);
  const int side = 129;
  const double res = 0.5;
  // thick walls: far away their borders can grow (a thin wall stays fine: see the topology test)
  auto w = rooms(side, res, rng, 8);
  LodMesher mesher(w.dem, {}, {}, w.cell_valid);
  OccupancyLod occ(mesher.dem(), mesher.error_map().tile_size(), mesher.cell_valid());
  // a sensed obstacle on free ground near the robot
  const int c = (side - 1) / 2;
  std::vector<std::pair<int, int>> obs;
  for (int j = c - 3; j <= c - 1; ++j) {
    for (int i = c + 12; i <= c + 14; ++i) {   // beyond the thick wall at c .. c + 7
      if (occ.cell_kind(i, j) == 0) {obs.push_back({i, j});}
    }
  }
  ASSERT_FALSE(obs.empty());
  occ.set_obstacles(obs);
  const auto p = params(res);
  const double rx = 7.0, ry = -3.0;   // in the room, 2 m from the sensed obstacle
  const auto mesh = mesher.extract(rx, ry, p, nullptr, nullptr, &occ);
  std::size_t occupied_inside = 0, occupied_coarse_inside = 0, free_large_inside = 0;
  std::size_t occupied_coarse_far = 0;
  for (std::size_t t = 0; t < mesh.triangle_count(); ++t) {
    const auto kinds = kinds_under(mesh, t, mesher.dem(), occ);
    const double far = farthest(mesh, t, rx, ry);
    const bool small = mesh.leg[t] <= p.min_size + 1e-6;
    if (far < p.high_res_radius) {
      if (occupied(kinds)) {
        ++occupied_inside;
        occupied_coarse_inside += !small;
      } else if (mesh.leg[t] > 2 * p.min_size) {
        ++free_large_inside;
      }
    } else if (far > p.high_res_radius + 10.0 && occupied(kinds) && !small) {
      ++occupied_coarse_far;
    }
  }
  EXPECT_GT(occupied_inside, 0u);
  EXPECT_EQ(occupied_coarse_inside,
    0u) << "occupied NavCels larger than the minimum near the robot";
  EXPECT_GT(free_large_inside, 0u) << "free ground near the robot should merge";
  EXPECT_GT(occupied_coarse_far, 0u) << "far away what is occupied grows, as in the distance mode";
}

TEST(OccupancyLod, TheHeightErrorStillSplitsFreeGround)
{
  const int side = 129;
  const double res = 0.5;
  auto p = params(res);
  // a ridge along x = 3.5 (not on an edge of the coarse levels): 10 degrees down on both sides,
  // no holes, no obstacles
  const double xr = 3.5;
  auto ridge = flat_dem(side, res);
  for (int j = 0; j < side; ++j) {
    for (int i = 0; i < side; ++i) {
      ridge.z[static_cast<std::size_t>(j) * side + i] =
        static_cast<float>(-std::tan(10.0 * M_PI / 180.0) * std::abs(ridge.x(i) - xr));
    }
  }
  LodMesher flat_mesher(flat_dem(side, res));
  LodMesher ridge_mesher(ridge);
  OccupancyLod flat_occ(flat_mesher.dem(), flat_mesher.error_map().tile_size(), {});
  OccupancyLod ridge_occ(ridge_mesher.dem(), ridge_mesher.error_map().tile_size(), {});
  const auto flat = flat_mesher.extract(0.0, 0.0, p, nullptr, nullptr, &flat_occ);
  const auto bumpy = ridge_mesher.extract(0.0, 0.0, p, nullptr, nullptr, &ridge_occ);
  EXPECT_LT(flat.triangle_count(), 64u) << "flat free ground merges up to the largest size";
  EXPECT_GT(bumpy.triangle_count(), flat.triangle_count());
  // no NavCel crosses the ridge near the robot (the error there is the ridge's height)
  for (std::size_t t = 0; t < bumpy.triangle_count(); ++t) {
    float x0 = 1e9f, x1 = -1e9f;
    for (int k = 0; k < 3; ++k) {
      x0 = std::min(x0, bumpy.x[bumpy.tri[3 * t + k]]);
      x1 = std::max(x1, bumpy.x[bumpy.tri[3 * t + k]]);
    }
    if (farthest(bumpy, t, 0.0, 0.0) < p.high_res_radius) {
      EXPECT_FALSE(x0<xr - 1e-3 && x1> xr + 1e-3) << "NavCel across the ridge";
    }
  }
}

TEST(OccupancyLod, IncrementalObstaclesMatchARebuild)
{
  const int side = 129;
  const double res = 0.5;
  std::mt19937 rng(7);
  auto w = rooms(side, res, rng);
  LodMesher mesher(w.dem, {}, {}, w.cell_valid);
  const auto p = params(res);
  OccupancyLod occ(mesher.dem(), mesher.error_map().tile_size(), mesher.cell_valid());
  const double rx = 3.0, ry = -2.0;
  const auto before = mesher.extract(rx, ry, p, nullptr, nullptr, &occ);

  const auto obs = random_obstacles(side, rng, 15);
  EXPECT_TRUE(occ.set_obstacles(obs));
  EXPECT_FALSE(occ.set_obstacles(obs)) << "the same obstacles change nothing";
  const auto with = mesher.extract(rx, ry, p, nullptr, nullptr, &occ);
  EXPECT_EQ(t_junctions(with, res), 0u);

  OccupancyLod fresh(mesher.dem(), mesher.error_map().tile_size(), mesher.cell_valid());
  fresh.set_obstacles(obs);
  const auto again = mesher.extract(rx, ry, p, nullptr, nullptr, &fresh);
  EXPECT_EQ(again.tri, with.tri);
  EXPECT_EQ(again.x, with.x);

  // half of them go: the same as a rebuild with the other half
  const std::vector<std::pair<int, int>> half(obs.begin(), obs.begin() + obs.size() / 2);
  EXPECT_TRUE(occ.set_obstacles(half));
  OccupancyLod fresh_half(mesher.dem(), mesher.error_map().tile_size(), mesher.cell_valid());
  fresh_half.set_obstacles(half);
  EXPECT_EQ(
    mesher.extract(rx, ry, p, nullptr, nullptr, &occ).tri,
    mesher.extract(rx, ry, p, nullptr, nullptr, &fresh_half).tri);

  // all gone: back to the first mesh
  EXPECT_TRUE(occ.set_obstacles({}));
  const auto after = mesher.extract(rx, ry, p, nullptr, nullptr, &occ);
  EXPECT_EQ(after.tri, before.tri);
  EXPECT_EQ(after.x, before.x);
}

// keep_path: the borders of the map's holes the robot has had within full_res (the observed
// KnownMap, never forgotten) stay at the minimum size when the robot is far; sensed obstacles do
// not (they are not marked)
TEST(OccupancyLod, KeptPathBordersStayFineButSensedObstaclesDoNot)
{
  std::mt19937 rng(5);
  const int side = 129;
  const double res = 0.5;
  auto w = rooms(side, res, rng, 8);   // thick walls: without keep_path their borders can grow
  LodMesher mesher(w.dem, {}, {}, w.cell_valid);
  OccupancyLod occ(mesher.dem(), mesher.error_map().tile_size(), mesher.cell_valid());
  auto p = params(res);
  p.observed_margin = 1.0;
  p.known_falloff = 9.0;
  const int c = (side - 1) / 2;
  // a sensed obstacle 5 m from the walls through the centre
  std::vector<std::pair<int, int>> obs;
  for (int i = c + 10; i <= c + 12; ++i) {
    if (occ.cell_kind(i, c - 10) == 0) {obs.push_back({i, c - 10});}
  }
  ASSERT_FALSE(obs.empty());
  occ.set_obstacles(obs);

  // the robot was at A (the centre): the hole borders within full_res are marked
  const double ax = 0.0, ay = 0.0;
  KnownMap kept(side, res);
  const auto border = hole_border_nodes(side, mesher.cell_valid());
  for (int j = 0; j < side; ++j) {
    for (int i = 0; i < side; ++i) {
      if (border[static_cast<std::size_t>(j) * side + i] &&
        std::hypot((i - c) * res - ax, (j - c) * res - ay) <= p.high_res_radius)
      {
        kept.mark_node(i, j);
      }
    }
  }
  kept.update();
  ASSERT_FALSE(kept.empty());

  // now it is far away
  const double bx = 26.0, by = -26.0;
  const auto with = mesher.extract(bx, by, p, nullptr, &kept, &occ);
  const auto without = mesher.extract(bx, by, p, nullptr, nullptr, &occ);
  EXPECT_EQ(t_junctions(with, res), 0u);
  auto coarse_borders_near_a = [&](const LodMesh & m) {
      std::size_t n = 0;
      for (std::size_t t = 0; t < m.triangle_count(); ++t) {
        float cx = 0.0f, cy = 0.0f;
        for (int k = 0; k < 3; ++k) {
          cx += m.x[m.tri[3 * t + k]] / 3.0f;
          cy += m.y[m.tri[3 * t + k]] / 3.0f;
        }
        if (std::hypot(cx - ax, cy - ay) > p.high_res_radius - 1.0) {continue;}
        const auto kinds = kinds_under(m, t, mesher.dem(), occ);
        if (kinds.count(0) && kinds.count(2) && m.leg[t] > p.min_size + 1e-6) {++n;}
      }
      return n;
    };
  EXPECT_EQ(coarse_borders_near_a(with), 0u) << "the kept borders must stay at the minimum size";
  EXPECT_GT(coarse_borders_near_a(without), 0u) << "without keep_path they merge when far";
  // the sensed obstacle near A is not kept: its NavCels are not all of the minimum size
  float largest_on_obstacle = 0.0f;
  for (std::size_t t = 0; t < with.triangle_count(); ++t) {
    if (kinds_under(with, t, mesher.dem(), occ).count(1)) {
      largest_on_obstacle = std::max(largest_on_obstacle, with.leg[t]);
    }
  }
  EXPECT_GT(largest_on_obstacle, p.min_size + 1e-6);
}

// Far from the robot no NavCel opens a wall or closes a door: a door of 4 cells in a wall of one
// cell keeps the two sides connected, and a wall with no door keeps them apart
namespace
{
// whether the NavCels left of x0 reach those right of x1 through NavCels that share an edge
bool connected_across(const LodMesh & m, float x0, float x1)
{
  const std::size_t n = m.triangle_count();
  std::map<std::pair<uint32_t, uint32_t>, std::vector<std::size_t>> by_edge;
  for (std::size_t t = 0; t < n; ++t) {
    for (int k = 0; k < 3; ++k) {
      const uint32_t a = m.tri[3 * t + k], b = m.tri[3 * t + (k + 1) % 3];
      by_edge[{std::min(a, b), std::max(a, b)}].push_back(t);
    }
  }
  auto cx = [&](std::size_t t) {
      return (m.x[m.tri[3 * t]] + m.x[m.tri[3 * t + 1]] + m.x[m.tri[3 * t + 2]]) / 3.0f;
    };
  std::vector<char> seen(n, 0);
  std::vector<std::size_t> stack;
  for (std::size_t t = 0; t < n; ++t) {
    if (cx(t) < x0) {seen[t] = 1; stack.push_back(t);}
  }
  while (!stack.empty()) {
    const auto t = stack.back();
    stack.pop_back();
    if (cx(t) > x1) {return true;}
    for (int k = 0; k < 3; ++k) {
      const uint32_t a = m.tri[3 * t + k], b = m.tri[3 * t + (k + 1) % 3];
      for (const auto u : by_edge[{std::min(a, b), std::max(a, b)}]) {
        if (!seen[u]) {seen[u] = 1; stack.push_back(u);}
      }
    }
  }
  return false;
}
}  // namespace

TEST(OccupancyLod, FarAwayNoNavCelOpensAWallOrClosesADoor)
{
  const int side = 129;
  const double res = 0.5;
  const int nc = side - 1, c = nc / 2;
  auto p = params(res);   // full_res 6 m, low_res_start 30 m, NavCels up to 16 m
  int wrong_without = 0;
  for (const int door : {4, 2, 0}) {
    for (const int wall_i : {c + 33, c + 40, c + 45}) {   // 16.5 .. 22.5 m from the centre
      std::vector<std::uint8_t> cv(static_cast<std::size_t>(nc) * nc, 1);
      const int d0 = 50;
      for (int j = 0; j < nc; ++j) {
        if (j >= d0 && j < d0 + door) {continue;}
        cv[static_cast<std::size_t>(j) * nc + wall_i] = 0;
      }
      LodMesher mesher(flat_dem(side, res), {}, {}, cv);
      OccupancyLod occ(mesher.dem(), mesher.error_map().tile_size(), mesher.cell_valid(),
        static_cast<int>(p.max_size / res));
      // the robot far on the other side: around the wall NavCels could be up to 16 m
      const double rx = -25.0, ry = -25.0;
      const auto mesh = mesher.extract(rx, ry, p, nullptr, nullptr, &occ);
      EXPECT_EQ(t_junctions(mesh, res), 0u);
      const float wx = static_cast<float>((wall_i - c) * res);
      EXPECT_EQ(connected_across(mesh, wx - 0.01f, wx + res + 0.01f), door > 0)
        << "door of " << door << " cells, wall at " << wx << " m";
      // without the topology (no triangle looked at: largest leg of one cell)
      OccupancyLod loose(mesher.dem(), mesher.error_map().tile_size(), mesher.cell_valid(), 1);
      const auto coarse = mesher.extract(rx, ry, p, nullptr, nullptr, &loose);
      wrong_without += connected_across(coarse, wx - 0.01f, wx + res + 0.01f) != (door > 0);
    }
  }
  EXPECT_GT(wrong_without, 0) << "the test should see a door closed or a wall opened without it";
}
