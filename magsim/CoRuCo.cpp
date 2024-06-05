
#include "ConfigReader.h"
#include "SimulationFactory.h"
#include "SpinDynamics.h"
#include "Metropolis.h"
#include "LatticeGenerator.h"
#include "CoRuCoGenerator.h"

#include <cstdio>
#include <string>
#include <chrono>
#include <filesystem>

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

  MapReader reader(argv[1]);
  Config c(argv[1]);
  CoRuCoGenerator gen(c);
  printf("generating spin lattice\n");
  SpinLattice lat = gen.Generate();

  auto sim = ConstructSimulation(reader, &lat);

  printf("dumping positions xyz\n");
  lat.DumpPositions(out_dir + "/positions.out");
  lat.DumpXYZ(out_dir + "/positions.xyz");

  if (reader.GetInt("dump_exchange")) {
    lat.DumpExchange(out_dir + "/exchange.out");
  }

  if (restart_fname.length()) {
    printf("loading restart file\n");
    lat.LoadLattice(restart_fname);
  }

  printf("starting sim\n");

  int64_t num_step = reader.GetInt("num_step");
  int64_t num_substep = reader.GetInt("num_substep");

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
    sim->lattice_->DumpProfile(out_dir + "/profile.out" + std::to_string(j), reader.GetChar("domain_wall_direction"), dump_avgs);
    sim->lattice_->ResetAverages();
  }
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
  printf("took %lu ms\n", duration.count());
  global_timer.PrintStatistics();
  return 0;
}
