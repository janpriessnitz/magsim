#include "LatticeGenerator.h"

#include <cmath>


HcpCobaltGenerator::HcpCobaltGenerator() {

}

real Co_anis = 2;

SpinLattice HcpCobaltGenerator::Generate() const {
  auto positions = GeneratePositions();

  PointLookup point_lookup(positions);

  SpinLattice res;
  res.positions_ = positions;
  res.anisotropy_ = Co_anis;
  res.exchange_ = GenerateExchange(point_lookup);

  res.spins_ = GenerateSpins(positions);
  res.Heffs_.resize(res.positions_.size());
  return res;
}

int nx = 10;
int ny = 5;
int nz = 5;
vec3d base1 = {1, 0, 0};
vec3d base2 = {0.5, 0.8660254037844386, 0};
vec3d base3 = {0, 0, 1.632993161855452};
std::vector<vec3d> spin_pos_list = {
  {0, 0, 0},
  (1/3.0)*(base1 + base2)
};

std::vector<vec3d> HcpCobaltGenerator::GeneratePositions() const {
  std::vector<vec3d> positions;

  for (int curx = 0; curx < nx; ++curx) {
    for (int cury = 0; cury < ny; ++cury) {
      for (int curz = 0; curz < nz; ++curz) {
        vec3d unit_cell_pos = curx*base1 + cury*base2 + curz*base3;
        for (const vec3d & spin_pos : spin_pos_list) {
          vec3d pos = unit_cell_pos + spin_pos;
          positions.emplace_back(pos);
        }
      }
    }
  }
  return positions;
}

real Co_J1 = 1;

std::vector<std::tuple<vec3d, real>> Js = {
  {base1, Co_J1},
  {-1*base1, Co_J1},
  {base2, Co_J1},
  {-1*base2, Co_J1},
  {base1-base2, Co_J1},
  {base2-base1, Co_J1},
};

  // #         [1, 1, 1, 0, 0, J1],
  // #         [1, 1, -1, 0, 0, J1],
  // #         [1, 1, 0, 1, 0, J1],
  // #         [1, 1, 0, -1, 0, J1],
  // #         [1, 1, 1, -1, 0, J1],
  // #         [1, 1, -1, 1, 0, J1],


std::vector<std::vector<std::tuple<size_t, real>>> HcpCobaltGenerator::GenerateExchange(const PointLookup & point_lookup) const {
  std::vector<std::vector<std::tuple<size_t, real>>> exch_list;
  for (size_t ind = 0; ind < point_lookup.points_.size(); ++ind) {
    std::vector<std::tuple<size_t, real>> one_exch;
    vec3d pos = point_lookup.points_[ind];
    for (const auto & [int_vec, int_energy] : Js) {
      vec3d partner_pos = pos + int_vec;
      auto partner_ind = point_lookup.GetExact(partner_pos);
      if (partner_ind) {
        one_exch.push_back({*partner_ind, int_energy});
      }
    }
    exch_list.push_back(one_exch);
  }
  return exch_list;
}

// Half of spins are {0, 0, -1}, half {0, 0, 1}
// Aim is to create a domain wall perpendicular to x-direction
std::vector<vec3d> HcpCobaltGenerator::GenerateSpins(std::vector<vec3d> positions) const {
  std::vector<vec3d> spins;
  spins.resize(positions.size());
  for (size_t ind = 0; ind < positions.size(); ++ind) {
    vec3d pos = positions[ind];
    if (std::get<0>(pos) < nx/2) {
      spins[ind] = {0, 0, -1};
    } else {
      spins[ind] = {0, 0, 1};
    }
  }
  return spins;
}


PointLookup::PointLookup(std::vector<vec3d> point_list)
  : points_(point_list)
{
  real minx, maxx, miny, maxy, minz, maxz;
  vec3d first_point = point_list[0];
  minx = maxx = std::get<0>(first_point);
  miny = maxy = std::get<1>(first_point);
  minz = maxz = std::get<2>(first_point);

  for (size_t ind = 0; ind < point_list.size(); ++ind) {
    vec3d pos = point_list[ind];

    real x, y, z;
    std::tie(x, y, z) = pos;

    minx = (x < minx) ? x : minx;
    miny = (y < miny) ? y : miny;
    minz = (z < minz) ? z : minz;

    maxx = (x > maxx) ? x : maxx;
    maxy = (y > maxy) ? y : maxy;
    maxz = (z > maxz) ? z : maxz;
  }

  real total_volume = (maxx - minx) * (maxy - miny) * (maxz - minz);
  // create a grid of cubes with ~ 1 point per cube
  real cube_volume = total_volume/point_list.size();
  cube_side_ = cbrt(cube_volume);
  origin_ = {minx - tol, miny - tol, minz - tol};
  n_x_ = ceil((maxx - minx + 2*tol)/cube_side_);
  n_y_ = ceil((maxy - miny + 2*tol)/cube_side_);
  n_z_ = ceil((maxz - minz + 2*tol)/cube_side_);

  grid_.resize(n_x_*n_y_*n_z_);

  for (size_t ind = 0; ind < point_list.size(); ++ind) {
    vec3d pos = point_list[ind];

    size_t cube_ind = GetCubeIndex(pos);
    grid_[cube_ind].push_back({ind, pos});
  }
}

size_t PointLookup::GetCubeIndex(int x_coord, int y_coord, int z_coord) const {
  size_t x, y, z;
  x = (x_coord + n_x_) % n_x_;
  y = (y_coord + n_y_) % n_y_;
  z = (z_coord + n_z_) % n_z_;
  return x + n_x_*y + n_x_*n_y_*z;
}

size_t PointLookup::GetCubeIndex(const vec3d & pos) const {
  int x_coord, y_coord, z_coord;
  std::tie(x_coord, y_coord, z_coord) = GetCubeCoords(pos);
  return GetCubeIndex(x_coord, y_coord, z_coord);
}

// TODO: periodic/vacuum boundary conditions
std::tuple<int, int, int> PointLookup::GetCubeCoords(const vec3d & pos) const {
  vec3d abspos = pos - origin_;
  real x, y, z;
  std::tie(x, y, z) = abspos;

  int x_coord, y_coord, z_coord;
  x_coord = x/cube_side_;
  y_coord = y/cube_side_;
  z_coord = z/cube_side_;
  return std::make_tuple(x_coord, y_coord, z_coord);
}


// TODO
// std::tuple<size_t, vec3d> PointLookup::GetClosest(const vec3d & pos) const {
//   return {};
// }

// search 3x3x3 cubes around the point
// return any point found within tolerance
// undefined behaviour if more points within tolerance
std::optional<size_t> PointLookup::GetExact(const vec3d & pos) const {
  int x_coord, y_coord, z_coord;
  std::tie(x_coord, y_coord, z_coord) = GetCubeCoords(pos);

  std::vector<std::tuple<size_t, vec3d>> close_points;

  // shortcut: try the center cube first
  size_t ind = GetCubeIndex(x_coord, y_coord, z_coord);
  for (const auto & point : grid_[ind]) {
    if (mag(pos - std::get<1>(point)) < tol) {
      return std::get<0>(point);
    }
  }

  for (int cur_x = x_coord - 1; cur_x <= x_coord + 1; ++cur_x) {
    for (int cur_y = y_coord - 1; cur_y <= y_coord + 1; ++cur_y) {
      for (int cur_z = z_coord - 1; cur_z <= z_coord + 1; ++cur_z) {
        size_t ind = GetCubeIndex(cur_x, cur_y, cur_z);
        for (const auto & point : grid_[ind]) {
          if (mag(pos - std::get<1>(point)) < tol) {
            return std::get<0>(point);
          }
        }
      }
    }
  }

  return {};
}

