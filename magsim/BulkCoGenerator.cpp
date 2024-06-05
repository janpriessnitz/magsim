
#include "BulkCoGenerator.h"

#include "ConfigReader.h"
#include "Constants.h"
#include "TupleReader.h"

#include <cmath>


BulkCoGenerator::BulkCoGenerator(const Config & config) : LatticeGenerator(config) {
  area_dims_ = nx_*base1_ + ny_*base2_ + nz_*base3_;
  Co_anis_ = config.data_["anisotropy"];

  syms_ = config.GetSymmetries();
  exchange_ints_ = config.GetExchange();
}

SpinLattice BulkCoGenerator::Generate() const {
  printf("generating positions\n");
  auto positions = GeneratePositions();

  PointLookup point_lookup(positions);

  SpinLattice res;
  res.positions_ = positions;
  res.anisotropy_.resize(res.positions_.size());
  std::fill(res.anisotropy_.begin(), res.anisotropy_.end(), Co_anis_);

  printf("generating exchange\n");
  res.exchange_ = GenerateExchange(point_lookup, exchange_ints_, syms_);

  printf("generating spins\n");
  res.spins_ = GenerateSpins(positions);
  res.avg_spins_ = res.spins_;
  res.n_avgs_ = 1;
  return res;
}


std::vector<vec3d> BulkCoGenerator::GeneratePositions() const {
  std::vector<vec3d> positions;
  double center_z = std::get<2>(area_dims_)/2;
  for (int curz = 0; curz < nz_; ++curz) {
    for (int cury = 0; cury < ny_; ++cury) {
      for (int curx = 0; curx < nx_; ++curx) {
        vec3d unit_cell_pos = curx*base1_ + cury*base2_ + curz*base3_;
        for (const vec3d & spin_pos : spin_pos_list) {
          vec3d pos = unit_cell_pos + spin_pos;
          positions.emplace_back(pos);
        }
      }
    }
  }
  return positions;
}

std::vector<std::vector<std::tuple<size_t, real>>> BulkCoGenerator::GenerateExchange(
  const PointLookup & point_lookup,
  const std::vector<std::tuple<vec3d, real>> & ints,
  const std::vector<mat3d> & syms) const
{
  auto start = std::chrono::high_resolution_clock::now();

  std::vector<std::vector<std::tuple<size_t, real>>> exch_list;
  exch_list.resize(point_lookup.points_.size());

  #pragma omp parallel for
  for (size_t ind = 0; ind < point_lookup.points_.size(); ++ind) {
    // std::vector<std::tuple<size_t, real>> one_exch;
    std::map<size_t, real> one_exch_map;
    vec3d pos = point_lookup.points_[ind];
    for (const auto & [int_vec, int_energy] : ints) {
      auto sym_vecs = ApplySymmetry(int_vec, syms);
      for (const auto & sym_vec : sym_vecs) {
        vec3d partner_pos = pos + sym_vec;
        auto [partner_ind, phase] = GetPoint(point_lookup, partner_pos);
        if (partner_ind) {
          one_exch_map[*partner_ind] = int_energy*phase;
        }
      }
    }

    std::vector<std::tuple<size_t, real>> one_exch;
    one_exch.insert(one_exch.end(), one_exch_map.begin(), one_exch_map.end());
    exch_list[ind] = one_exch;
  }

  auto stop = std::chrono::high_resolution_clock::now();
  global_timer.AddTime(stop - start, Timer::Section::GenExchange);

  return exch_list;
}

// Half of spins are {0, 0, -1}, half {0, 0, 1}
// Aim is to create a domain wall perpendicular to z- direction
std::vector<vec3d> BulkCoGenerator::GenerateSpins(const std::vector<vec3d> & positions) const {
  std::vector<vec3d> spins;
  spins.resize(positions.size());
  for (size_t ind = 0; ind < positions.size(); ++ind) {
    vec3d pos = positions[ind];
    if (std::get<2>(pos) < std::get<2>(area_dims_)/2) {
      spins[ind] = {0, 0, -1};
    } else {
      spins[ind] = {0, 0, 1};
    }
  }
  return spins;
}
