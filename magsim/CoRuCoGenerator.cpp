
#include "CoRuCoGenerator.h"

#include "Constants.h"
#include "TupleReader.h"

#include <cmath>


CoRuCoGenerator::CoRuCoGenerator(const Config & config) : LatticeGenerator(config) {
  area_dims_ = vec3d({nx_, ny_, nz_})*cell_;
  Co_anis_ = config.Get<real>("anisotropy");

  syms_ = config.GetSymmetries();
  exchange_ints_ = config.GetExchange();

  interface_exchange_energy_ = (config.Get<real>("interface_J"))*constants::Ry;
  spin_direction_ = config.Get<std::string>("spin_direction")[0];
}

SpinLattice CoRuCoGenerator::Generate() const {
  printf("generating positions\n");
  auto positions = GeneratePositions();
  printf("aftergenerating positions %lu\n", positions.size());

  PointLookup point_lookup(positions);

  printf("after point lookup\n");

  SpinLattice res;
  res.positions_ = positions;
  res.anisotropy_.resize(res.positions_.size());
  // res.anisotropy_ = Co_anis_;
  std::fill(res.anisotropy_.begin(), res.anisotropy_.end(), Co_anis_);

  printf("generating exchange\n");
  res.exchange_ = GenerateExchange(point_lookup, exchange_ints_, syms_);

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
        vec3d unit_cell_pos = vec3d({curx, cury, curz})*cell_;
        for (const vec3d & direct_spin_pos : spin_pos_list_) {
          vec3d cart_spin_pos = direct_spin_pos*cell_;
          vec3d pos = unit_cell_pos + cart_spin_pos;
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
  double interface_A_z = nz_*std::get<2>(std::get<2>(cell_))/2 - middle_space/2;

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
          // TODO: deal with duplicates
          // especially for small supercells with periodic boundary conditions, multiple exchange interactions can map to a single pair of atoms -> we need to add them, not replace by the last read interaction
          // if (one_exch_map[*partner_ind] != 0 && ind == 0) {
          //   printf("conflict %lu;%lu;%s;%s;%lg,%lg\n", ind, *partner_ind, to_string(pos).c_str(), to_string(partner_pos).c_str(), int_energy, one_exch_map[*partner_ind]);
          // }
          // one_exch_map[*partner_ind] += int_energy*phase;

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

  for (size_t ind = 0; ind < point_lookup.points_.size(); ++ind) {
    vec3d pos = point_lookup.points_[ind];
    double tol = 1;
    if (std::get<2>(pos) + tol > interface_A_z && std::get<2>(pos) - tol < interface_A_z) {
      vec3d partner_pos = pos + 7*std::get<2>(cell_);
      auto [partner_ind, phase] = GetPoint(point_lookup, partner_pos);
      if (partner_ind) {
        exch_list[ind].push_back({*partner_ind, interface_exchange_energy_*phase});
        exch_list[*partner_ind].push_back({ind, interface_exchange_energy_*phase});
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
