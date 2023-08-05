
#include "ConfigReader.h"
#include "SpinDynamics.h"
#include "LatticeGenerator.h"

#include <cstdio>
#include <string>
#include <chrono>

int main(int argc, char **argv) {
  HcpCobaltGenerator gen;
  printf("generating spin lattice\n");
  SpinLattice lat = gen.Generate();
  printf("spin lattice generated\n");

  SpinDynamics* dyn = new SpinDynamics(&lat);
  dyn->alpha_ = 0.1;
  dyn->temperature_ = 1;
  dyn->timestep_ = 1e-16;

  printf("dumping exchange\n");
  lat.DumpExchange("exchange.out");
  printf("dumping positions\n");
  lat.DumpPositions("positions.out");

  printf("starting sim\n");
  auto start = std::chrono::high_resolution_clock::now();
  for (int j = 0; j < 20; ++j) {
    for (int i = 0; i < 2000; ++i) {
      dyn->DoStep();
    }
    dyn->lattice_->PrintEnergy();
    auto avgm = dyn->lattice_->AvgM();
    printf("%s %lf\n", to_string(avgm).c_str(), mag(avgm));
  }
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
  printf("took %lu ms\n", duration.count());
  dyn->lattice_->DumpLattice("lattice.out");
  return 0;
}
