
#include "ConfigReader.h"
#include "SpinDynamics.h"
#include "LatticeGenerator.h"

#include <cstdio>
#include <string>
#include <chrono>

int main(int argc, char **argv) {
  // ConfigReader c(argv[1]);
  // Config conf = c.ReadConfig();

  HcpCobaltGenerator gen;
  SpinLattice lat = gen.Generate();
  SpinDynamics* dyn = new SpinDynamics(&lat);

  lat.DumpExchange("exchange.out");
  lat.DumpPositions("positions.out");

  printf("starting sim\n");
  auto start = std::chrono::high_resolution_clock::now();
  for (int j = 0; j < 100; ++j) {
    for (int i = 0; i < 1000; ++i) {
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
