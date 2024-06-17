
#include "Config.h"
#include "SimulationFactory.h"
#include "SpinDynamics.h"
#include "Metropolis.h"
#include "LatticeGenerator.h"
#include "BulkCoGenerator.h"

#include <cstdio>
#include <string>
#include <chrono>
#include <filesystem>

class SimpleInterfaceGenerator : public LatticeGenerator {
public:
  SimpleInterfaceGenerator(const Config & config) : LatticeGenerator(config) {
    pos_list_ = config.Get<std::vector<vec3d>>("positions");
    anis_ = config.Get<real>("anisotropy")*constants::Ry;
    J_ = config.Get<real>("exchange")*constants::Ry;
    J_interface_ = config.Get<real>("exchange_interface")*constants::Ry;
    spin_direction_ = config.Get<std::string>("spin_direction")[0];
  }

  std::tuple<std::vector<vec3d>, std::vector<vec3d>> GeneratePositionsAndSpins() const {
    std::vector<vec3d> positions;
    std::vector<vec3d> spins;
    for (int curz = 0; curz < nz_; ++curz) {
      for (int cury = 0; cury < ny_; ++cury) {
        for (int curx = 0; curx < nx_; ++curx) {
          vec3d unit_cell_pos = curx*std::get<0>(cell_) + cury*std::get<1>(cell_) + curz*std::get<2>(cell_);
          for (const vec3d & spin_pos : pos_list_) {
            vec3d pos = unit_cell_pos + spin_pos;
            vec3d spin;
            if (curz <= nz_/2 && spin_direction_ == 'z') {
              spin = {0, 0, -1};
            } else {
              spin = {0, 0, 1};
            }
            positions.emplace_back(pos);
            spins.emplace_back(spin);
          }
        }
      }
    }
    return {positions, spins};
  }

  std::vector<std::vector<std::tuple<size_t, real>>> GenerateExchange(
  const PointLookup & point_lookup) const
  {
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<vec3d> int_vecs = {
      std::get<0>(cell_),
      -1*std::get<0>(cell_),
      std::get<1>(cell_),
      -1*std::get<1>(cell_),
      std::get<2>(cell_),
      -1*std::get<2>(cell_),
    };

    std::vector<std::vector<std::tuple<size_t, real>>> exch_list;
    exch_list.resize(point_lookup.points_.size());

    // #pragma omp parallel for
    for (size_t ind = 0; ind < point_lookup.points_.size(); ++ind) {

      std::map<size_t, real> one_exch_map;
      vec3d pos = point_lookup.points_[ind];

      for (const auto & int_vec : int_vecs) {
        vec3d partner_pos = pos + int_vec;
        real int_energy;
        if ((std::get<2>(pos) <= nz_/2 && std::get<2>(partner_pos) > nz_/2) || (std::get<2>(pos) > nz_/2 && std::get<2>(partner_pos) <= nz_/2)) {
          int_energy = J_interface_;
        } else {
          int_energy = J_;
        }

        auto [partner_ind, phase] = GetPoint(point_lookup, partner_pos);
        if (partner_ind) {
          one_exch_map[*partner_ind] = int_energy*phase;
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

  SpinLattice Generate() const {
    SpinLattice res;

    auto [positions, spins] = GeneratePositionsAndSpins();
    res.positions_ = positions;
    res.spins_ = spins;
    res.anisotropy_.resize(res.positions_.size());
    std::fill(res.anisotropy_.begin(), res.anisotropy_.end(), anis_);

    PointLookup point_lookup(res.positions_);
    res.exchange_ = GenerateExchange(point_lookup);

    res.avg_spins_ = res.spins_;
    res.n_avgs_ = 1;

    return res;
  }

  std::vector<vec3d> pos_list_;
  real anis_;
  real J_;
  real J_interface_;
  char spin_direction_;
};

int main(int argc, char **argv) {
  Config c(argc, argv);
  SimpleInterfaceGenerator gen(c);
  SpinLattice lat = gen.Generate();

  auto sim = ConstructSimulation(c, &lat);

  lat.DumpPositions(c.out_dir_ + "/positions.out");
  lat.DumpXYZ(c.out_dir_ + "/positions.xyz");

  if (c.Get<int>("dump_exchange")) {
    lat.DumpExchange(c.out_dir_ + "/exchange.out");
  }

  if (!c.restart_file_.empty()) {
    printf("loading restart file\n");
    lat.LoadLattice(c.restart_file_);
  }

  printf("starting sim\n");
  sim->Run();
  return 0;
}
