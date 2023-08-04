
#include "ConfigReader.h"
#include "SpinDynamics.h"
#include "LatticeGenerator.h"

#include <cstdio>
#include <string>
#include <chrono>

int main(int argc, char **argv) {
  // ConfigReader c(argv[1]);
  // Config conf = c.ReadConfig();

  // HcpCobaltGenerator gen;
  // SpinLattice lat = gen.Generate();

  SpinLattice lat;
  lat.anisotropy_ = -1e-24;
  lat.spins_ = {{1, 0, 1}};
  lat.Heffs_.resize(1);
  lat.exchange_.resize(1);

  SpinDynamics* dyn = new SpinDynamics(&lat);
  dyn->alpha_ = 0.1;
  dyn->temperature_ = 10;
  // dyn->timestep_ = 1e-15;

  lat.DumpExchange("exchange.out");
  lat.DumpPositions("positions.out");

  printf("starting sim\n");
  auto start = std::chrono::high_resolution_clock::now();
  for (int j = 0; j < 100; ++j) {
    for (int i = 0; i < 10000; ++i) {
      dyn->DoStep();
    }
    // dyn->lattice_->PrintEnergy();
    auto avgm = dyn->lattice_->AvgM();
    printf("%lg %s %lf\n", (j+1)*1000*dyn->timestep_, to_string(avgm).c_str(), mag(avgm));
  }
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
  printf("took %lu ms\n", duration.count());
  dyn->lattice_->DumpLattice("lattice.out");
  return 0;
}
