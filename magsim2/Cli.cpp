
#include "ConfigReader.h"
#include "SpinDynamics.h"
#include "LatticeGenerator.h"

#include <cstdio>
#include <string>

int main(int argc, char **argv) {
  // ConfigReader c(argv[1]);
  // Config conf = c.ReadConfig();
  printf("main\n");

  HcpCobaltGenerator gen;
  SpinLattice lat = gen.Generate();
  SpinDynamics* dyn = new SpinDynamics(&lat);

  lat.DumpExchange("exchange.out");
  lat.DumpPositions("positions.out");

  for (int j = 0; j < 100; ++j) {
    for (int i = 0; i < 10000; ++i) {
      dyn->DoStep();
    }
    real x, y, z;
    std::tie(x, y, z) = dyn->lattice_->AvgM();
    printf("%lf %lf %lf\n", x, y, z);
  }
  dyn->lattice_->DumpLattice("lattice.out");
  dyn->lattice_->DumpHeffs("heffs.out");
  return 0;
}