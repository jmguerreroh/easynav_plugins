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

#ifndef EASYNAV_BIGMAP_MAPS_MANAGER__OCCUPANCY_LOD_HPP_
#define EASYNAV_BIGMAP_MAPS_MANAGER__OCCUPANCY_LOD_HPP_

#include <cstdint>
#include <utility>
#include <vector>

#include "easynav_bigmap_maps_manager/dem_grid.hpp"

namespace bigmap_maps_manager
{

/// Where the grid is occupied, for the occupancy block of bigmap (`occupancy.*`): which
/// nodes of the RTIN tree have something occupied under them.
///
/// Each cell of the grid is ground, a sensed obstacle (on ground) or no ground (a hole of the map:
/// a wall, a rack). A node is OCCUPIED if the cells under the two triangles that share it as their
/// hypotenuse midpoint (a diamond) hold a sensed obstacle, or hold both ground and no ground (the
/// border of a hole). A diamond with no ground at all is EMPTY: nothing is shown there.
///
/// The cells of a node are kept as a mask of the kinds found under it, made from the masks of its
/// children's nodes (which cover a little more than its diamond: they are diamonds themselves).
/// So occupied(child) => occupied(parent), which is what an RTIN needs to be free of T-junctions
/// (see lod_policy.hpp). The holes are static; the sensed obstacles change (set_obstacles), and
/// only the nodes around what changed are computed again.
///
/// TOPOLOGY. A coarse NavCel stands for ground or for a hole as a whole (the cell under its
/// centroid decides, as LodMesher emits it). Over the border of a hole that can change what is
/// connected: shown as ground over a wall, it opens it; left out over a door or an aisle, it closes
/// it. Such a triangle BREAKS THE TOPOLOGY, and its node is split whatever its distance:
///   * it is shown as what is not most of it (ground over mostly holes: an island; left out over
///     mostly ground: a cut);
///   * shown as ground, the ground under it (4-connected) is in more than one piece; left out, its
///     holes (8-connected) are;
///   * shown as ground, it holds a cell of a wall thinner than its leg; left out, a cell of a
///     passage narrower than its leg. The width of the passage of a ground cell (the thickness
///     of the wall of a hole cell) is the shortest run of its kind through it, horizontal,
///     vertical or diagonal. It depends on the cells only, so a door or a wall shared by two
///     coarse triangles is kept by both.
/// A node is flagged if one of its triangles breaks the topology or one of its children's nodes is
/// flagged (closed upwards too): the zone is refined down to the largest NavCels that still show
/// the door (ground through it) and the wall (a hole across it). Only the holes of the map count
/// (the sensed obstacles are not kept this way).
class OccupancyLod
{
public:
  /// `tile_size`: the RTIN tile (RtinErrorMap::tile_size(), >= dem.side - 1). `cell_valid`: per DEM
  /// cell, 1 = ground (empty: all ground). `max_leg_cells`: legs (in cells) of the largest NavCel
  /// the mesher can leave unsplit (cell_size.max / resolution); larger triangles are always split,
  /// so their topology is not looked at (0: all of them are).
  OccupancyLod(
    const DemGrid & dem, int tile_size, const std::vector<std::uint8_t> & cell_valid,
    int max_leg_cells = 0);

  OccupancyLod(const OccupancyLod &) = delete;
  OccupancyLod & operator=(const OccupancyLod &) = delete;

  /// The cells (i, j) of the DEM with a sensed obstacle now (any order, repeated ones allowed; out
  /// of the grid or on no ground: ignored). Returns true if any of them changed.
  bool set_obstacles(const std::vector<std::pair<int, int>> & cells);

  /// Something occupied under the node (x, y) of the RTIN tree (a hypotenuse midpoint)
  bool occupied(int x, int y) const
  {
    const auto m = masks_[idx(x, y)];
    return (m & kSensed) || ((m & kNoGround) && (m & kGround));
  }
  /// No ground at all under the node: nothing to split, nothing to show
  bool no_ground(int x, int y) const {return masks_[idx(x, y)] == kNoGround;}
  /// A triangle of the node (or of a node below it) would open a wall or close a door if it were
  /// not split (see TOPOLOGY above)
  bool breaks_topology(int x, int y) const {return (masks_[idx(x, y)] & kTopology) != 0;}
  /// Nodes flagged by their own triangles (log, tests)
  std::size_t topology_count() const;

  /// Kind of a cell: 0 ground, 1 sensed obstacle, 2 no ground (tests, debugging)
  int cell_kind(int i, int j) const;
  std::size_t obstacle_count() const {return obstacles_.size();}
  std::size_t memory_bytes() const;

private:
  static constexpr std::uint8_t kGround = 1, kSensed = 2, kNoGround = 4, kTopology = 8;

  std::size_t idx(int x, int y) const
  {
    return static_cast<std::size_t>(y) * static_cast<std::size_t>(grid_) +
           static_cast<std::size_t>(x);
  }
  std::uint8_t cell_mask(int i, int j) const;
  bool ground(int i, int j) const;   // real cell with ground
  void build_widths();               // width_
  /// The triangle (grid node coordinates; c: the right angle) breaks the topology
  bool breaks(int ax, int ay, int bx, int by, int cx, int cy) const;
  /// Flags the nodes of the RTIN tree below a triangle; returns whether it or one below is flagged
  bool flag_topology(int ax, int ay, int bx, int by, int cx, int cy);
  /// The masks of every node above the cells [cx0, cx1] x [cy0, cy1]
  void rebuild(int cx0, int cy0, int cx1, int cy1);

  int cells_{0};       // real cells per axis (dem.side - 1)
  int tile_{0};        // RTIN tile size
  int grid_{0};        // nodes per axis (tile + 1)
  std::vector<std::uint8_t> ground_;        // per real cell: 1 = ground
  std::vector<std::uint8_t> sensed_;        // per real cell: 1 = sensed obstacle
  std::vector<std::uint32_t> obstacles_;    // sorted cell indices with a sensed obstacle
  std::vector<std::uint8_t> masks_;         // per node
  std::vector<bool> topology_;              // per node: flagged by its own triangles
  std::vector<std::uint8_t> width_;         // per real cell: passage width / wall thickness (cells)
  int max_leg_{0};                          // cells; 0: no limit
};

}  // namespace bigmap_maps_manager

#endif  // EASYNAV_BIGMAP_MAPS_MANAGER__OCCUPANCY_LOD_HPP_
