#include "LatticeGenerator.h"
#include "TupleReader.h"
#include "Constants.h"

#include <cmath>


HcpCobaltGenerator::HcpCobaltGenerator() {
}

real Co_anis = -5.83e-24;  // J
int nx = 500;
int ny = 30;
int nz = 30;
vec3d base1 = {1, 0, 0};
vec3d base2 = {0.5, 0.8660254037844386, 0};
vec3d base3 = {0, 0, 1.632993161855452};
mat3d base_mat = {base1, base2, base3};

std::vector<vec3d> spin_pos_list = {
  {0, 0, 0},
  (1/3.0)*(base1 + base2) + (1/2.0)*base3
};

SpinLattice HcpCobaltGenerator::Generate() const {
  printf("generating positions\n");
  auto positions = GeneratePositions();

  PointLookup point_lookup(positions);

  SpinLattice res;
  res.positions_ = positions;
  res.anisotropy_ = Co_anis;

  printf("loading symmetries\n");
  auto syms = LoadSymmetries("sym.mat");
  printf("loading exchange\n");
  auto exchange_ints = LoadExchange("exchange.in");

  printf("generating exchange\n");
  res.exchange_ = GenerateExchange(point_lookup, exchange_ints, syms);

  printf("generating spins\n");
  res.spins_ = GenerateSpins(positions);
  res.Heffs_.resize(res.positions_.size());
  return res;
}



std::vector<vec3d> HcpCobaltGenerator::GeneratePositions() const {
  std::vector<vec3d> positions;

  for (int curz = 0; curz < nz; ++curz) {
    for (int cury = 0; cury < ny; ++cury) {
      for (int curx = 0; curx < nx; ++curx) {
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


std::vector<std::vector<std::tuple<size_t, real>>> HcpCobaltGenerator::GenerateExchange(
  const PointLookup & point_lookup,
  const std::vector<std::tuple<vec3d, real>> & ints,
  const std::vector<mat3d> & syms) const
{
  std::vector<std::vector<std::tuple<size_t, real>>> exch_list;
  for (size_t ind = 0; ind < point_lookup.points_.size(); ++ind) {
    // std::vector<std::tuple<size_t, real>> one_exch;
    std::map<size_t, real> one_exch_map;
    vec3d pos = point_lookup.points_[ind];
    for (const auto & [int_vec, int_energy] : ints) {
      auto sym_vecs = ApplySymmetry(int_vec, syms);
      for (const auto & sym_vec : sym_vecs) {
        vec3d partner_pos = pos + sym_vec;
        auto partner_ind = GetPoint(point_lookup, partner_pos);
        if (partner_ind) {
          one_exch_map[*partner_ind] = int_energy;
        }
      }
    }
    std::vector<std::tuple<size_t, real>> one_exch;
    one_exch.insert(one_exch.end(), one_exch_map.begin(), one_exch_map.end());
    exch_list.push_back(one_exch);
  }
  return exch_list;
}

// Half of spins are {0, 0, -1}, half {0, 0, 1}
// Aim is to create a domain wall perpendicular to x-direction
std::vector<vec3d> HcpCobaltGenerator::GenerateSpins(const std::vector<vec3d> & positions) const {
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

// TODO: periodic boundary conditions
std::optional<size_t> HcpCobaltGenerator::GetPoint(const PointLookup & lookup, const vec3d & pos) const {
  auto partner_ind = lookup.GetExact(pos);
  return partner_ind;
}

// P6/mmc (no. 194)
std::vector<mat3d> hcp_syms = {
  {{1, 2, 3},
   {4, 5, 6},
   {7, 8, 9}}
};

std::vector<mat3d> HcpCobaltGenerator::LoadSymmetries(const std::string & fname) const {
  TupleReader reader = TupleReader(fname);
  std::vector<mat3d> res;
  int n_syms = reader.GetInt(0, 0);
  for (int i_sym = 0; i_sym < n_syms; ++i_sym) {
    int sr = i_sym*3 + 1;
    res.push_back({
      {reader.GetDouble(sr, 0), reader.GetDouble(sr, 1), reader.GetDouble(sr, 2)},
      {reader.GetDouble(sr+1, 0), reader.GetDouble(sr+1, 1), reader.GetDouble(sr+1, 2)},
      {reader.GetDouble(sr+2, 0), reader.GetDouble(sr+2, 1), reader.GetDouble(sr+2, 2)},
    });
  }
  return res;
}

std::vector<std::tuple<vec3d, real>> HcpCobaltGenerator::LoadExchange(const std::string & fname) const {
  TupleReader reader = TupleReader(fname);
  std::vector<std::tuple<vec3d, real>> res;
  for (int i = 0; i < reader.NumRows(); ++i) {
    vec3d vec = {reader.GetDouble(i, 0), reader.GetDouble(i, 1), reader.GetDouble(i, 2)};
    // maptype = 2 in UppASD
    vec = vec*base_mat;
    real energy = reader.GetDouble(i, 3)*constants::Ry*1e-3;
    energy *= 2;  // consistent with UppASD, TODO: Remove
    res.push_back({vec, energy});
  }
  return res;
}


std::vector<vec3d> HcpCobaltGenerator::ApplySymmetry(const vec3d & vec, const std::vector<mat3d> & syms) const {
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

