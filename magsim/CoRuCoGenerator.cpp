
#include "CoRuCoGenerator.h"

#include "ConfigReader.h"
#include "Constants.h"
#include "TupleReader.h"

#include <cmath>


CoRuCoGenerator::CoRuCoGenerator(const MapReader & config) {
  nx_ = config.GetInt("nx");
  ny_ = config.GetInt("ny");
  nz_ = config.GetInt("nz");
  area_dims_ = nx_*base1_ + ny_*base2_ + nz_*base3_;
  Co_anis_ = config.GetDouble("anisotropy");
  symmetry_fname_ = config.GetString("symmetry_file");
  exchange_fname_ = config.GetString("exchange_file");

  periodic_x_ = config.GetInt("periodic_x") != 0;
  periodic_y_ = config.GetInt("periodic_y") != 0;
  periodic_z_ = config.GetInt("periodic_z") != 0;
  interface_exchange_energy_ = config.GetDouble("interface_J")*constants::Ry;
  spin_direction_ = config.GetChar("spin_direction");
}

SpinLattice CoRuCoGenerator::Generate() const {
  printf("generating positions\n");
  auto positions = GeneratePositions();

  PointLookup point_lookup(positions);

  SpinLattice res;
  res.positions_ = positions;
  res.anisotropy_ = Co_anis_;

  printf("loading symmetries\n");
  auto syms = ConfigReader::ReadSymmetries(symmetry_fname_);
  printf("loading exchange\n");
  auto exchange_ints = ConfigReader::ReadExchange(exchange_fname_);

  printf("generating exchange\n");
  res.exchange_ = GenerateExchange(point_lookup, exchange_ints, syms);

  printf("generating spins\n");
  res.spins_ = GenerateSpins(positions);
  res.avg_spins_ = res.spins_;
  res.n_avgs_ = 1;
  return res;
}


std::vector<vec3d> CoRuCoGenerator::GeneratePositions() const {
  std::vector<vec3d> positions;
  double center_z = std::get<2>(area_dims_)/2;
  double middle_space = 10;
  for (int curz = 0; curz < nz_; ++curz) {
    for (int cury = 0; cury < ny_; ++cury) {
      for (int curx = 0; curx < nx_; ++curx) {
        vec3d unit_cell_pos = curx*base1_ + cury*base2_ + curz*base3_;
        for (const vec3d & spin_pos : spin_pos_list) {
          vec3d pos = unit_cell_pos + spin_pos;
          // make a hole on the interface
          if (abs(center_z - std::get<2>(pos)) < middle_space/2) {
              continue;
          }
          positions.emplace_back(pos);

        }
      }
    }
  }
  return positions;
}

std::vector<std::vector<std::tuple<size_t, real>>> CoRuCoGenerator::GenerateExchange(
  const PointLookup & point_lookup,
  const std::vector<std::tuple<vec3d, real>> & ints,
  const std::vector<mat3d> & syms) const
{
  auto start = std::chrono::high_resolution_clock::now();

  std::vector<std::vector<std::tuple<size_t, real>>> exch_list;
  exch_list.resize(point_lookup.points_.size());

  double middle_space = 10;
  double interface_A_z = nz_*std::get<2>(base3_)/2 - middle_space/2;

  #pragma omp parallel for
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
    exch_list[ind] = one_exch;
  }
  auto stop = std::chrono::high_resolution_clock::now();
  global_timer.AddTime(stop - start, Timer::Section::GenExchange);

  for (size_t ind = 0; ind < point_lookup.points_.size(); ++ind) {
    vec3d pos = point_lookup.points_[ind];
    double tol = 1;
    if (std::get<2>(pos) + tol > interface_A_z && std::get<2>(pos) - tol < interface_A_z) {
      vec3d partner_pos = pos + 7*base3_;
      auto partner_ind = GetPoint(point_lookup, partner_pos);
      if (partner_ind) {
        exch_list[ind].push_back({*partner_ind, interface_exchange_energy_});
        exch_list[*partner_ind].push_back({ind, interface_exchange_energy_});
      } else {
        printf("no partner across interface found!!!\n");
        exit(0);
      }
    }
  }

  return exch_list;
}

// Half of spins are {0, 0, -1}, half {0, 0, 1}
// Aim is to create a domain wall perpendicular to x- or z- direction
std::vector<vec3d> CoRuCoGenerator::GenerateSpins(const std::vector<vec3d> & positions) const {
  std::vector<vec3d> spins;
  spins.resize(positions.size());
  for (size_t ind = 0; ind < positions.size(); ++ind) {
    vec3d pos = positions[ind];
    if (spin_direction_ == 'z') {
      if (std::get<2>(pos) < std::get<2>(area_dims_)/2) {
        spins[ind] = {0, 0, -1};
      } else {
        spins[ind] = {0, 0, 1};
      }
    }
    else {
      spins[ind] = {0, 0, 1};
    }
  }
  return spins;
}

std::optional<size_t> CoRuCoGenerator::GetPoint(const PointLookup & lookup, const vec3d & pos) const {
  vec3d base_pos = pos + vec3d{tol/2, tol/2, tol/2};
  base_pos = base_pos*base_mat_inv_;
  auto [bx, by, bz] = base_pos;
  if (periodic_x_) {
    bx = fmod(bx+nx_, nx_);
  }
  if (periodic_y_) {
    by = fmod(by+ny_, ny_);
  }
  if (periodic_z_) {
    bz = fmod(bz+nz_, nz_);
  }

  vec3d new_pos = vec3d{bx, by, bz}*base_mat_;
  auto partner_ind = lookup.GetExact(new_pos);
  return partner_ind;
}

std::vector<vec3d> CoRuCoGenerator::ApplySymmetry(const vec3d & vec, const std::vector<mat3d> & syms) const {
  std::vector<vec3d> res;
  for (const auto & sym : syms) {
    res.push_back(vec*sym);
  }
  return res;
}