
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
    pos_list_ = config.Get("positions");
    anis_ = config.Get("anisotropy").get<real>()*constants::Ry;
    J_ = config.Get("exchange").get<real>()*constants::Ry;
    J_interface_ = config.Get("exchange_interface").get<real>()*constants::Ry;
    spin_direction_ = config.Get("spin_direction").get<std::string>()[0];
  }

  std::tuple<std::vector<vec3d>, std::vector<vec3d>> GeneratePositions() const {
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
        if (std::get<2>(pos) <= nz_/2 && std::get<2>(partner_pos) > nz_/2) {
          int_energy = J_interface_;
        } else {
          int_energy = J_;
        }

        // printf("%s, %s\n", to_string(pos).c_str(), to_string(partner_pos).c_str());
        auto [partner_ind, phase] = GetPoint(point_lookup, partner_pos);
        if (partner_ind) {
          // printf("yes %lu %lu\n", ind, *partner_ind);
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

    auto [positions, spins] = GeneratePositions();
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
  std::string out_dir = "output/";
  if (argc > 2) {
    out_dir = argv[2];
  }
  std::string restart_fname = "";
  if (argc > 3) {
    restart_fname = argv[3];
  }
  if (!std::filesystem::create_directory(out_dir)) {
    fprintf(stderr, "failed to create output directory %s\n", out_dir.c_str());
  }

  Config c(argv[1]);
  SimpleInterfaceGenerator gen(c);
  printf("generating spin lattice\n");
  SpinLattice lat = gen.Generate();

  auto sim = ConstructSimulation(c, &lat);

  printf("dumping positions xyz\n");
  lat.DumpPositions(out_dir + "/positions.out");
  lat.DumpXYZ(out_dir + "/positions.xyz");

  if (c.data_["dump_exchange"].get<int>()) {
    lat.DumpExchange(out_dir + "/exchange.out");
  }

  if (restart_fname.length()) {
    printf("loading restart file\n");
    lat.LoadLattice(restart_fname);
  }

  printf("starting sim\n");

  int64_t num_step = c.Get("num_step");
  int64_t num_substep = c.Get("num_substep");

  auto start = std::chrono::high_resolution_clock::now();
  for (int j = 0; j < num_step; ++j) {
    for (int i = 0; i < num_substep; ++i) {
      sim->DoStep();
    }
    sim->lattice_->PrintEnergy();
    auto avgm = sim->lattice_->AvgM();
    printf("%s %lf\n", to_string(avgm).c_str(), mag(avgm));
    bool dump_avgs = true;
    sim->lattice_->DumpLattice(out_dir + "/lattice.out" + std::to_string(j), dump_avgs);
    sim->lattice_->DumpProfile(out_dir + "/profile.out" + std::to_string(j), c.Get("domain_wall_direction").get<std::string>()[0], dump_avgs);
    sim->lattice_->ResetAverages();
  }
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
  printf("took %lu ms\n", duration.count());
  global_timer.PrintStatistics();
  return 0;
}
