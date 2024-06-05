
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
    anis_ = config.Get("anisotropy");
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
            if (curz <= nz_/2) {
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

  SpinLattice Generate() const {
    SpinLattice res;

    auto [positions, spins] = GeneratePositions();
    res.positions_ = positions;
    res.spins_ = spins;
    res.anisotropy_.resize(res.positions_.size());
    std::fill(res.anisotropy_.begin(), res.anisotropy_.end(), anis_);

    PointLookup point_lookup(res.positions_);
    // res.exchange_ = GenerateExchange(point_lookup, exchange_ints_, syms_);
    
    res.avg_spins_ = res.spins_;
    res.n_avgs_ = 1;

    return res;
  }

  std::vector<vec3d> pos_list_;
  real anis_;
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
