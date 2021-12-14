
#include "Config.h"
#include "Metropolis.h"

#include <cstdio>
#include <string>

int main() {
  Config conf;
  Metropolis sim(conf);
  printf("E0 %lf\n", sim.lattice_.getEnergy());
  for (int j = 0; j < 200; ++j) {
    for (int i = 0; i < 5000000; ++i) {
      sim.do_step();
    }
    printf("E100000 %lf\n", sim.lattice_.getEnergy());
    sim.T_ *= 0.97;
    sim.lattice_.dump("lattice.dump" + std::to_string(j));

  }
  return 0;
}