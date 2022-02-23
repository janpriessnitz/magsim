
#include "ConfigReader.h"
#include "Metropolis.h"

#include <cstdio>
#include <string>

int main(int argc, char **argv) {
  ConfigReader c(argv[1]);
  Config conf = c.ReadConfig();
  Metropolis sim(conf);
  for (int j = 0; j < conf.T_steps; ++j) {
    sim.equilibrize();
    printf("T %le E %le\n", sim.T_, sim.lattice_.getEnergy());
    sim.lattice_.dump("lattice.dump" + std::to_string(j), sim.T_);

    sim.T_ *= conf.T_step_ratio;
  }
  return 0;
}