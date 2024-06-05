#include "LatticeGenerator.h"

#include "Constants.h"
#include "TupleReader.h"
#include "Config.h"

#include <cmath>

LatticeGenerator::LatticeGenerator(const Config & config) {
  nx_ = config.Get("nx");
  ny_ = config.Get("ny");
  nz_ = config.Get("nz");

  periodic_x_ = config.Get("periodic_x");
  periodic_y_ = config.Get("periodic_y");
  periodic_z_ = config.Get("periodic_z");

  cell_ = config.Get("cell");
  cell_inv_ = inverse(cell_);
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

// HcpCobaltGenerator::HcpCobaltGenerator(const Config & config) : LatticeGenerator(config) {
//   area_dims_ = nx_*base1_ + ny_*base2_ + nz_*base3_;
//   Co_anis_ = config.data_["anisotropy"];
//   symmetry_fname_ = config.data_["symmetry_file"];
//   exchange_fname_ = config.data_["exchange_file"];

//   domain_wall_direction_ = ((std::string)config.data_["domain_wall_direction"])[0];
//   middle_space_ = config.data_["middle_space"];
//   middle_offset_ = config.data_["middle_offset"];
// }

// SpinLattice HcpCobaltGenerator::Generate() const {
//   printf("generating positions\n");
//   auto positions = GeneratePositions();

//   PointLookup point_lookup(positions);

//   SpinLattice res;
//   res.positions_ = positions;
//   // res.anisotropy_ = Co_anis_;
//   std::fill(res.anisotropy_.begin(), res.anisotropy_.end(), Co_anis_);

//   printf("loading symmetries\n");
//   auto syms = ConfigReader::ReadSymmetries(symmetry_fname_);
//   printf("loading exchange\n");
//   auto exchange_ints = ConfigReader::ReadExchange(exchange_fname_);

//   printf("generating exchange\n");
//   res.exchange_ = GenerateExchange(point_lookup, exchange_ints, syms);

//   printf("generating spins\n");
//   res.spins_ = GenerateSpins(positions);
//   res.avg_spins_ = res.spins_;
//   res.n_avgs_ = 1;
//   return res;
// }



// std::vector<vec3d> HcpCobaltGenerator::GeneratePositions() const {
//   std::vector<vec3d> positions;

//   for (int curz = 0; curz < nz_; ++curz) {
//     for (int cury = 0; cury < ny_; ++cury) {
//       for (int curx = 0; curx < nx_; ++curx) {
//         vec3d unit_cell_pos = curx*base1_ + cury*base2_ + curz*base3_;
//         for (const vec3d & spin_pos : spin_pos_list) {
//           vec3d pos = unit_cell_pos + spin_pos;
//           if (domain_wall_direction_ == 'z') {
//             double center_z = std::get<2>(area_dims_)/2 + middle_offset_;
//             if (abs(center_z - std::get<2>(pos)) < middle_space_)
//               continue;
//           } else if (domain_wall_direction_ == 'x') {
//             double center_x = std::get<0>(area_dims_)/2 + middle_offset_;
//             if (abs(center_x - std::get<0>(pos)) < middle_space_)
//               continue;
//           } else {
//             fprintf(stderr, "GeneratePositions: unknown domain wall direction: %c\n", domain_wall_direction_);
//             exit(1);
//           }
//           positions.emplace_back(pos);

//         }
//       }
//     }
//   }
//   return positions;
// }


// std::vector<std::vector<std::tuple<size_t, real>>> HcpCobaltGenerator::GenerateExchange(
//   const PointLookup & point_lookup,
//   const std::vector<std::tuple<vec3d, real>> & ints,
//   const std::vector<mat3d> & syms) const
// {
//   auto start = std::chrono::high_resolution_clock::now();

//   std::vector<std::vector<std::tuple<size_t, real>>> exch_list;
//   exch_list.resize(point_lookup.points_.size());

//   #pragma omp parallel for
//   for (size_t ind = 0; ind < point_lookup.points_.size(); ++ind) {
//     // std::vector<std::tuple<size_t, real>> one_exch;
//     std::map<size_t, real> one_exch_map;
//     vec3d pos = point_lookup.points_[ind];
//     for (const auto & [int_vec, int_energy] : ints) {
//       auto sym_vecs = ApplySymmetry(int_vec, syms);
//       for (const auto & sym_vec : sym_vecs) {
//         vec3d partner_pos = pos + sym_vec;
//         auto [partner_ind, phase] = GetPoint(point_lookup, partner_pos);
//         if (partner_ind) {
//           one_exch_map[*partner_ind] = int_energy;
//         }
//       }
//     }
//     std::vector<std::tuple<size_t, real>> one_exch;
//     one_exch.insert(one_exch.end(), one_exch_map.begin(), one_exch_map.end());
//     exch_list[ind] = one_exch;
//   }
//   auto stop = std::chrono::high_resolution_clock::now();
//   global_timer.AddTime(stop - start, Timer::Section::GenExchange);
//   return exch_list;
// }

// // Half of spins are {0, 0, -1}, half {0, 0, 1}
// // Aim is to create a domain wall perpendicular to x- or z- direction
// std::vector<vec3d> HcpCobaltGenerator::GenerateSpins(const std::vector<vec3d> & positions) const {
//   std::vector<vec3d> spins;
//   spins.resize(positions.size());
//   for (size_t ind = 0; ind < positions.size(); ++ind) {
//     vec3d pos = positions[ind];
//     if (domain_wall_direction_ == 'z') {
//       if (std::get<2>(pos) < std::get<2>(area_dims_)/2) {
//         spins[ind] = {0, 0, -1};
//       } else {
//         spins[ind] = {0, 0, 1};
//       }
//     } else if (domain_wall_direction_ == 'x') {
//       if (std::get<0>(pos) < std::get<0>(area_dims_)/2) {
//         spins[ind] = {0, 0, -1};
//       } else {
//         spins[ind] = {0, 0, 1};
//       }
//     } else if (domain_wall_direction_ == '0') {
//       spins[ind] = {0, 0, 1};
//     } else {
//       fprintf(stderr, "unknown domain wall direction: %c\n", domain_wall_direction_);
//       exit(1);
//     }

//     // spins[ind] = {0, 0, 1};
//   }
//   return spins;
// }

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

