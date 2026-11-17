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

#include "easynav_bigmap_maps_manager/occupancy_lod.hpp"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <stdexcept>

namespace bigmap_maps_manager
{

// The tree of nodes, level by level (s = 2, 4, ..., tile), as RTIN builds it:
//   * B nodes (s): the midpoints of the axis-aligned hypotenuses of length s, i.e. the midpoints
//     of the edges of the squares of side s: (k s + s/2, j s) and (k s, j s + s/2). Their diamond
//     (the two triangles that share that hypotenuse) lies in the four squares of side h = s/2 that
//     have the node as a corner: its mask is made from theirs.
//   * A nodes (s): the midpoints of the diagonal hypotenuses of the triangles with legs s, i.e. the
//     centres of the squares of side s: (k s + h, j s + h). Their diamond is that square, whose
//     four edges are B nodes of the same level: its mask is made from theirs.
// The squares of side 1 are the cells. Each mask is the union of masks that cover at least what it
// covers, so occupied(child) => occupied(parent).

OccupancyLod::OccupancyLod(
  const DemGrid & dem, int tile_size, const std::vector<std::uint8_t> & cell_valid,
  int max_leg_cells)
: max_leg_(max_leg_cells)
{
  if (dem.side < 3 || tile_size < dem.side - 1 || (tile_size & (tile_size - 1)) != 0) {
    throw std::invalid_argument("OccupancyLod: the tile must be a power of two >= side - 1");
  }
  cells_ = dem.side - 1;
  tile_ = tile_size;
  grid_ = tile_ + 1;
  const std::size_t nc = static_cast<std::size_t>(cells_) * cells_;
  if (!cell_valid.empty() && cell_valid.size() != nc) {
    throw std::invalid_argument("OccupancyLod: cell_valid must have one byte per cell");
  }
  ground_.assign(nc, 1);
  for (std::size_t k = 0; k < cell_valid.size(); ++k) {
    ground_[k] = cell_valid[k] ? 1 : 0;
  }
  sensed_.assign(nc, 0);
  masks_.assign(static_cast<std::size_t>(grid_) * grid_, 0);
  topology_.assign(static_cast<std::size_t>(grid_) * grid_, false);
  rebuild(0, 0, cells_ - 1, cells_ - 1);
  // the topology of the holes (static): the masks above tell where there is something to look at,
  // then the flags go into the masks
  build_widths();
  flag_topology(0, 0, tile_, tile_, tile_, 0);
  flag_topology(tile_, tile_, 0, 0, 0, tile_);
  rebuild(0, 0, cells_ - 1, cells_ - 1);
}

bool OccupancyLod::ground(int i, int j) const
{
  return i >= 0 && j >= 0 && i < cells_ && j < cells_ &&
         ground_[static_cast<std::size_t>(j) * cells_ + i] != 0;
}

void OccupancyLod::build_widths()
{
  // per direction (horizontal, vertical, two diagonals), the length of the run of ground through
  // each ground cell, and of hole through each hole cell; the width of its passage (the thickness
  // of its wall) is the shortest of them (a diagonal step is sqrt 2 long)
  const std::size_t nc = static_cast<std::size_t>(cells_) * cells_;
  width_.assign(nc, 255);
  const int dirs[4][2] = {{1, 0}, {0, 1}, {1, 1}, {1, -1}};
  std::vector<std::size_t> run;
  for (const auto & d : dirs) {
    const double step = (d[0] != 0 && d[1] != 0) ? M_SQRT2 : 1.0;
    // the starts of every line in this direction: cells whose previous cell is out of the grid
    for (int j = 0; j < cells_; ++j) {
      for (int i = 0; i < cells_; ++i) {
        const int pi = i - d[0], pj = j - d[1];
        if (pi >= 0 && pj >= 0 && pi < cells_ && pj < cells_) {continue;}
        int a = i, b = j;
        run.clear();
        bool kind = true;        // the kind of the current run
        bool from_edge = true;   // it started at the edge of the grid
        // a run of holes at the edge of the grid goes on beyond it (no ground there): not thin
        auto flush = [&](bool to_edge) {
            if (!run.empty() && (kind || !(from_edge || to_edge))) {
              const auto len = static_cast<int>(std::lround(run.size() * step));
              const auto w = static_cast<std::uint8_t>(std::min(len, 255));
              for (const auto k : run) {
                width_[k] = std::min(width_[k], w);
              }
            }
            run.clear();
          };
        while (a >= 0 && b >= 0 && a < cells_ && b < cells_) {
          const auto k = static_cast<std::size_t>(b) * cells_ + a;
          const bool g = ground_[k] != 0;
          if (g != kind) {
            flush(false);
            from_edge = run.empty() && a == i && b == j;
            kind = g;
          }
          run.push_back(k);
          a += d[0];
          b += d[1];
        }
        flush(true);
      }
    }
  }
}

bool OccupancyLod::breaks(int ax, int ay, int bx, int by, int cx, int cy) const
{
  // what the mesher would make of it: ground if the cell under its centroid has ground
  const double gx = (ax + bx + cx) / 3.0, gy = (ay + by + cy) / 3.0;
  const bool shown = ground(static_cast<int>(std::floor(gx)), static_cast<int>(std::floor(gy)));
  // the cells whose centre is in the triangle (or on its border), on a small local grid
  const int i0 = std::min({ax, bx, cx}), i1 = std::max({ax, bx, cx});
  const int j0 = std::min({ay, by, cy}), j1 = std::max({ay, by, cy});
  const int w = i1 - i0, h = j1 - j0;
  if (w <= 0 || h <= 0) {return false;}
  auto side = [](double px, double py, double qx, double qy, double x, double y) {
      return (qx - px) * (y - py) - (qy - py) * (x - px);
    };
  // 0: not in the triangle, 1: the kind that would be lost or kept apart, 2: the other one
  std::vector<std::uint8_t> cell(static_cast<std::size_t>(w) * h, 0);
  bool any_ground = false, any_hole = false;
  int n_ground = 0, n_hole = 0;
  for (int j = 0; j < h; ++j) {
    for (int i = 0; i < w; ++i) {
      const double x = i0 + i + 0.5, y = j0 + j + 0.5;
      const double s0 = side(ax, ay, bx, by, x, y), s1 = side(bx, by, cx, cy, x, y);
      const double s2 = side(cx, cy, ax, ay, x, y);
      if (!((s0 >= 0 && s1 >= 0 && s2 >= 0) || (s0 <= 0 && s1 <= 0 && s2 <= 0))) {continue;}
      const bool g = ground(i0 + i, j0 + j);
      (g ? any_ground : any_hole) = true;
      ++(g ? n_ground : n_hole);
      // left out over ground of a passage narrower than the triangle, it could close it; shown as
      // ground over a wall thinner than the triangle, it could open it (alone or with a neighbour:
      // it depends on the cells only)
      // (only real cells: the padding beyond the map is not a wall)
      if (g != shown && i0 + i < cells_ && j0 + j < cells_ &&
        width_[static_cast<std::size_t>(j0 + j) * cells_ + (i0 + i)] < std::hypot(ax - cx, ay - cy))
      {
        return true;
      }
      // shown as ground: the pieces of ground count; left out: the pieces of hole
      cell[static_cast<std::size_t>(j) * w + i] = (g == shown) ? 1 : 2;
    }
  }
  if (!any_ground || !any_hole) {return false;}
  // shown as what is not most of it: a mostly free triangle left out cuts the ground around it, a
  // mostly occupied one shown makes an island of ground
  if (shown ? n_hole > n_ground : n_ground > n_hole) {return true;}
  // pieces of the kind it would be shown as: ground 4-connected, holes 8-connected
  const bool eight = !shown;
  int pieces = 0;
  std::vector<int> stack;
  for (int start = 0; start < w * h; ++start) {
    if (cell[start] != 1) {continue;}
    if (++pieces > 1) {return true;}
    cell[start] = 0;
    stack.push_back(start);
    while (!stack.empty()) {
      const int k = stack.back();
      stack.pop_back();
      const int i = k % w, j = k / w;
      for (int dj = -1; dj <= 1; ++dj) {
        for (int di = -1; di <= 1; ++di) {
          if ((di == 0 && dj == 0) || (!eight && di != 0 && dj != 0)) {continue;}
          const int a = i + di, b = j + dj;
          if (a < 0 || b < 0 || a >= w || b >= h) {continue;}
          const int n = b * w + a;
          if (cell[n] == 1) {
            cell[n] = 0;
            stack.push_back(n);
          }
        }
      }
    }
  }
  return false;
}

bool OccupancyLod::flag_topology(int ax, int ay, int bx, int by, int cx, int cy)
{
  if (std::abs(ax - cx) + std::abs(ay - cy) <= 1) {return false;}          // a leaf
  if (std::min({ax, bx, cx}) >= cells_ || std::min({ay, by, cy}) >= cells_) {return false;}
  const int mx = (ax + bx) >> 1, my = (ay + by) >> 1;
  const bool below = flag_topology(cx, cy, ax, ay, mx, my) | flag_topology(bx, by, cx, cy, mx, my);
  if (below) {return true;}   // the masks carry it up
  const auto m = masks_[idx(mx, my)];
  if (!((m & kGround) || (m & kSensed)) || !(m & kNoGround)) {return false;}   // one kind only
  // (the leg as the mesher measures it)
  if (max_leg_ > 0 && std::hypot(ax - cx, ay - cy) > max_leg_ + 1e-9) {return false;}
  if (!breaks(ax, ay, bx, by, cx, cy)) {return false;}
  topology_[idx(mx, my)] = true;
  return true;
}

std::size_t OccupancyLod::topology_count() const
{
  return static_cast<std::size_t>(std::count(topology_.begin(), topology_.end(), true));
}

std::uint8_t OccupancyLod::cell_mask(int i, int j) const
{
  if (i < 0 || j < 0 || i >= tile_ || j >= tile_) {return 0;}          // beyond the tile: nothing
  if (i >= cells_ || j >= cells_) {return kNoGround;}                   // the padding
  const std::size_t k = static_cast<std::size_t>(j) * cells_ + i;
  if (!ground_[k]) {return kNoGround;}
  return sensed_[k] ? kSensed : kGround;
}

int OccupancyLod::cell_kind(int i, int j) const
{
  const auto m = cell_mask(i, j);
  return m == kSensed ? 1 : (m == kGround ? 0 : 2);
}

void OccupancyLod::rebuild(int cx0, int cy0, int cx1, int cy1)
{
  // Region whose masks changed at the level below, in node coordinates (cell centres first)
  double x0 = cx0 + 0.5, y0 = cy0 + 0.5, x1 = cx1 + 0.5, y1 = cy1 + 0.5;
  auto lo = [](double v, int off, int s) {return static_cast<int>(std::ceil((v - off) / s));};
  auto hi = [](double v, int off, int s) {return static_cast<int>(std::floor((v - off) / s));};
  for (int s = 2; s <= tile_; s *= 2) {
    const int h = s / 2, n = tile_ / s;
    // the mask of the square of side h with corner (sx, sy) (A node of level h, or a cell)
    auto square = [&](int sx, int sy) -> std::uint8_t {
        if (sx < 0 || sy < 0 || sx >= tile_ || sy >= tile_) {return 0;}
        return h == 1 ? cell_mask(sx, sy) : masks_[idx(sx + h / 2, sy + h / 2)];
      };
    auto set_b = [&](int ex, int ey) {
        const auto k = idx(ex, ey);
        masks_[k] = static_cast<std::uint8_t>(
          square(ex - h, ey - h) | square(ex, ey - h) | square(ex - h, ey) | square(ex, ey) |
          (topology_[k] ? kTopology : 0));
      };
    // B nodes: those whose four squares reach the changed region
    const double bx0 = x0 - 0.5 * h, by0 = y0 - 0.5 * h, bx1 = x1 + 0.5 * h, by1 = y1 + 0.5 * h;
    for (int j = std::max(0, lo(by0, 0, s)); j <= std::min(n, hi(by1, 0, s)); ++j) {
      for (int k = std::max(0, lo(bx0, h, s)); k <= std::min(n - 1, hi(bx1, h, s)); ++k) {
        set_b(k * s + h, j * s);
      }
    }
    for (int j = std::max(0, lo(by0, h, s)); j <= std::min(n - 1, hi(by1, h, s)); ++j) {
      for (int k = std::max(0, lo(bx0, 0, s)); k <= std::min(n, hi(bx1, 0, s)); ++k) {
        set_b(k * s, j * s + h);
      }
    }
    // A nodes: those with a changed edge
    const double ax0 = bx0 - h, ay0 = by0 - h, ax1 = bx1 + h, ay1 = by1 + h;
    for (int j = std::max(0, lo(ay0, h, s)); j <= std::min(n - 1, hi(ay1, h, s)); ++j) {
      for (int k = std::max(0, lo(ax0, h, s)); k <= std::min(n - 1, hi(ax1, h, s)); ++k) {
        const int cx = k * s + h, cy = j * s + h;
        const auto node = idx(cx, cy);
        masks_[node] = static_cast<std::uint8_t>(
          masks_[idx(cx, cy - h)] | masks_[idx(cx, cy + h)] | masks_[idx(cx - h, cy)] |
          masks_[idx(cx + h, cy)] | (topology_[node] ? kTopology : 0));
      }
    }
    x0 = ax0;
    y0 = ay0;
    x1 = ax1;
    y1 = ay1;
  }
}

bool OccupancyLod::set_obstacles(const std::vector<std::pair<int, int>> & cells)
{
  std::vector<std::uint32_t> now;
  now.reserve(cells.size());
  for (const auto & [i, j] : cells) {
    if (i < 0 || j < 0 || i >= cells_ || j >= cells_) {continue;}
    const std::size_t k = static_cast<std::size_t>(j) * cells_ + i;
    if (!ground_[k]) {continue;}
    now.push_back(static_cast<std::uint32_t>(k));
  }
  std::sort(now.begin(), now.end());
  now.erase(std::unique(now.begin(), now.end()), now.end());
  std::vector<std::uint32_t> changed;
  std::set_symmetric_difference(
    now.begin(), now.end(), obstacles_.begin(), obstacles_.end(), std::back_inserter(changed));
  if (changed.empty()) {return false;}
  for (const auto c : obstacles_) {
    sensed_[c] = 0;
  }
  for (const auto c : now) {
    sensed_[c] = 1;
  }
  obstacles_ = std::move(now);
  int cx0 = cells_, cy0 = cells_, cx1 = -1, cy1 = -1;
  for (const auto c : changed) {
    const int i = static_cast<int>(c % static_cast<std::uint32_t>(cells_));
    const int j = static_cast<int>(c / static_cast<std::uint32_t>(cells_));
    cx0 = std::min(cx0, i);
    cy0 = std::min(cy0, j);
    cx1 = std::max(cx1, i);
    cy1 = std::max(cy1, j);
  }
  rebuild(cx0, cy0, cx1, cy1);
  return true;
}

std::size_t OccupancyLod::memory_bytes() const
{
  return ground_.capacity() + sensed_.capacity() +
         obstacles_.capacity() * sizeof(std::uint32_t) + masks_.capacity() +
         topology_.capacity() / 8 + width_.capacity();
}

}  // namespace bigmap_maps_manager
