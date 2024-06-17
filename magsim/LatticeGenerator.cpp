#include "LatticeGenerator.h"

#include "Constants.h"
#include "TupleReader.h"
#include "Config.h"

#include <cmath>

LatticeGenerator::LatticeGenerator(const Config & config) {
  nx_ = config.Get<int>("nx");
  ny_ = config.Get<int>("ny");
  nz_ = config.Get<int>("nz");

  periodic_x_ = config.Get<int>("periodic_x");
  periodic_y_ = config.Get<int>("periodic_y");
  periodic_z_ = config.Get<int>("periodic_z");

  cell_ = config.Get<mat3d>("cell");
  cell_inv_ = inverse(cell_);

  spin_pos_list_ = config.Get<std::vector<vec3d>>("positions");
}


std::vector<vec3d> LatticeGenerator::GeneratePositions() const {
  std::vector<vec3d> positions;

  for (int curz = 0; curz < nz_; ++curz) {
    for (int cury = 0; cury < ny_; ++cury) {
      for (int curx = 0; curx < nx_; ++curx) {
        vec3d unit_cell_pos = vec3d({curx, cury, curz})*cell_;
        for (const vec3d & direct_spin_pos : spin_pos_list_) {
          vec3d cart_spin_pos = direct_spin_pos*cell_;
          vec3d pos = unit_cell_pos + cart_spin_pos;
          positions.emplace_back(pos);
        }
      }
    }
  }
  return positions;
}

std::pair<std::optional<size_t>, int> LatticeGenerator::GetPoint(const PointLookup & lookup, const vec3d & pos) const {
  vec3d base_pos = pos + vec3d{tol/2, tol/2, tol/2};
  base_pos = base_pos*cell_inv_;
  auto [bx, by, bz] = base_pos;
  int phase = 1;
  if (periodic_x_ != 0) {
    real bx_p = fmod(bx+nx_, nx_);
    if (abs(bx - bx_p) > tol) {
      phase *= periodic_x_;
    }
    bx = bx_p;
  }
  if (periodic_y_ != 0) {
    real by_p = fmod(by+ny_, ny_);
    if (abs(by - by_p) > tol) {
      phase *= periodic_y_;
    }
    by = by_p;
  }
  if (periodic_z_ != 0) {
    real bz_p = fmod(bz+nz_, nz_);
    if (abs(bz - bz_p) > tol) {
      phase *= periodic_z_;
    }
    bz = bz_p;
  }

  vec3d new_pos = vec3d{bx, by, bz}*cell_;
  auto partner_ind = lookup.GetExact(new_pos);
  return std::make_pair(partner_ind, phase);
}

std::vector<vec3d> LatticeGenerator::ApplySymmetry(const vec3d & vec, const std::vector<mat3d> & syms) const {
  std::vector<vec3d> res;
  for (const auto & sym : syms) {
    res.push_back(vec*sym);
  }
  return res;
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

    auto [x, y, z] = pos;

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

    auto cube_ind = GetCubeIndex(pos);
    if (!cube_ind) {
      auto [x, y, z] = pos;
      fprintf(stderr, "fail in PointLookup init!!! Cannot insert point %lf %lf %lf\n", x, y, z);
    }
    grid_[*cube_ind].push_back({ind, pos});
  }
}

std::optional<size_t> PointLookup::GetCubeIndex(int x_coord, int y_coord, int z_coord) const {
  if (x_coord < 0 || x_coord >= n_x_) {
    return {};
  }
  if (y_coord < 0 || y_coord >= n_y_) {
    return {};
  }
  if (z_coord < 0 || z_coord >= n_z_) {
    return {};
  }

  return x_coord + n_x_*y_coord + n_x_*n_y_*z_coord;
}

std::optional<size_t> PointLookup::GetCubeIndex(const vec3d & pos) const {
  auto [x_coord, y_coord, z_coord] = GetCubeCoords(pos);
  return GetCubeIndex(x_coord, y_coord, z_coord);
}

std::tuple<int, int, int> PointLookup::GetCubeCoords(const vec3d & pos) const {
  vec3d abspos = pos - origin_;
  auto [x, y, z] = abspos;

  int x_coord, y_coord, z_coord;
  x_coord = x/cube_side_;
  y_coord = y/cube_side_;
  z_coord = z/cube_side_;
  return {x_coord, y_coord, z_coord};
}


// TODO
// std::tuple<size_t, vec3d> PointLookup::GetClosest(const vec3d & pos) const {
//   return {};
// }

// search 3x3x3 cubes around the point
// return any point found within tolerance
// undefined behaviour if more points within tolerance
std::optional<size_t> PointLookup::GetExact(const vec3d & pos) const {
  auto [x_coord, y_coord, z_coord] = GetCubeCoords(pos);

  std::vector<std::tuple<size_t, vec3d>> close_points;

  // shortcut: try the center cube first
  auto ind = GetCubeIndex(x_coord, y_coord, z_coord);
  if (ind) {
    for (const auto & point : grid_[*ind]) {
      if (mag(pos - std::get<1>(point)) < tol) {
        return std::get<0>(point);
      }
    }
  }

  for (int cur_x = x_coord - 1; cur_x <= x_coord + 1; ++cur_x) {
    for (int cur_y = y_coord - 1; cur_y <= y_coord + 1; ++cur_y) {
      for (int cur_z = z_coord - 1; cur_z <= z_coord + 1; ++cur_z) {
        auto ind = GetCubeIndex(cur_x, cur_y, cur_z);
        if (ind) {
          for (const auto & point : grid_[*ind]) {
            if (mag(pos - std::get<1>(point)) < tol) {
              return std::get<0>(point);
            }
          }
        }
      }
    }
  }

  return {};
}

