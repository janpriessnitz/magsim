
#include "ConfigReader.h"
#include "SpinDynamics.h"
#include "LatticeGenerator.h"

#include <cstdio>
#include <string>

int main(int argc, char **argv) {
  // ConfigReader c(argv[1]);
  // Config conf = c.ReadConfig();

  HcpCobaltGenerator gen;
  SpinLattice lat = gen.Generate();
  SpinDynamics* dyn = new SpinDynamics(&lat);

  lat.DumpExchange("exchange.out");
  lat.DumpPositions("positions.out");

  for (int j = 0; j < 2000; ++j) {
    for (int i = 0; i < 10; ++i) {
      dyn->DoStep();
    }
    dyn->lattice_->PrintEnergy();
    auto [x, y, z] = dyn->lattice_->AvgM();
    printf("%lf %lf %lf\n", x, y, z);
  }
  dyn->lattice_->DumpLattice("lattice.out");
  dyn->lattice_->DumpHeffs("heffs.out");
  return 0;
}
