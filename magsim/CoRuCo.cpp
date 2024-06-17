
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
  Config c(argc, argv);
  CoRuCoGenerator gen(c);
  printf("generating spin lattice\n");
  SpinLattice lat = gen.Generate();

  auto sim = ConstructSimulation(c, &lat);

  printf("dumping positions xyz\n");
  lat.DumpPositions(c.out_dir_ + "/positions.out");
  lat.DumpXYZ(c.out_dir_ + "/positions.xyz");

  if (c.Get<int>("dump_exchange")) {
    lat.DumpExchange(c.out_dir_ + "/exchange.out");
  }

  if (!c.restart_file_.empty()) {
    printf("loading restart file\n");
    lat.LoadLattice(c.restart_file_);
  }

  printf("starting sim\n");

  sim->Run();
  return 0;
}
