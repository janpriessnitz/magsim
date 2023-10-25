
#include "ConfigReader.h"
#include "SpinDynamics.h"
#include "LatticeGenerator.h"
#include "CoRuCoGenerator.h"

#include <cstdio>
#include <string>
#include <chrono>
#include <filesystem>

int main(int argc, char **argv) {
  Timer timer;

  std::string out_dir = "output/";
  if (argc > 2) {
    out_dir = argv[2];
  }
  if (!std::filesystem::create_directory(out_dir)) {
    fprintf(stderr, "failed to create output directory %s\n", out_dir.c_str());
  }

  MapReader reader(argv[1]);
  CoRuCoGenerator gen(reader, timer);
  printf("generating spin lattice\n");
  SpinLattice lat = gen.Generate();

  SpinDynamics* dyn = new SpinDynamics(&lat, timer);
  dyn->alpha_ = reader.GetFloat("damping");
  dyn->temperature_ = reader.GetFloat("temperature");
  dyn->timestep_ = reader.GetFloat("timestep");

  printf("dumping positions xyz\n");
  lat.DumpPositions(out_dir + "positions.out");
  lat.DumpXYZ(out_dir + "positions.xyz");

  if (reader.GetInt("dump_exchange")) {
    lat.DumpExchange(out_dir + "exchange.out");
  }

  printf("starting sim\n");

  int64_t num_step = reader.GetInt("num_step");
  int64_t num_substep = reader.GetInt("num_substep");

  auto start = std::chrono::high_resolution_clock::now();
  for (int j = 0; j < num_step; ++j) {
    for (int i = 0; i < num_substep; ++i) {
      dyn->DoStep();
    }
    dyn->lattice_->PrintEnergy();
    auto avgm = dyn->lattice_->AvgM();
    printf("%s %lf\n", to_string(avgm).c_str(), mag(avgm));
    dyn->lattice_->DumpLattice(out_dir + "lattice.out" + std::to_string(j));
    dyn->lattice_->DumpProfile(out_dir + "profile.out" + std::to_string(j), reader.GetChar("domain_wall_direction"));
  }
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
  printf("took %lu ms\n", duration.count());
  timer.PrintStatistics();
  return 0;
}
