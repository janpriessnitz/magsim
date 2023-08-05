
#include "ConfigReader.h"
#include "SpinDynamics.h"
#include "LatticeGenerator.h"

#include <cstdio>
#include <string>
#include <chrono>

int main(int argc, char **argv) {
  MapReader reader(argv[1]);
  HcpCobaltGenerator gen(reader);
  printf("generating spin lattice\n");
  SpinLattice lat = gen.Generate();
  printf("spin lattice generated\n");

  SpinDynamics* dyn = new SpinDynamics(&lat);
  dyn->alpha_ = reader.GetFloat("damping");
  dyn->temperature_ = reader.GetFloat("temperature");
  dyn->timestep_ = reader.GetFloat("timestep");

  printf("dumping positions\n");
  lat.DumpPositions("positions.out");

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
    dyn->lattice_->DumpLattice("lattice.out" + std::to_string(j));
  }
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
  printf("took %lu ms\n", duration.count());
  return 0;
}
