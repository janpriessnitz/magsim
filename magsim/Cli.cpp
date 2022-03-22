
#include "ConfigReader.h"
#include "Metropolis.h"
#include "HeatBath.h"

#include <cstdio>
#include <string>

int main(int argc, char **argv) {
  ConfigReader c(argv[1]);
  Config conf = c.ReadConfig();
  Simulation* sim;
  switch(conf.method) {
    case 'H':
      sim = new HeatBath(conf);
      break;
    case 'M':
      sim = new Metropolis(conf);
      break;
    default:
      fprintf(stderr, "unrecognized method in config file: %c", conf.method);
      return 1;
  }
  uint64_t steps_total = 0;
  for (int tstep = 0; tstep < conf.annealing_sched.size(); ++tstep) {
    uint64_t steps_T = 0;

    real T;
    uint64_t steps, reporting_period;
    std::tie(T, steps, reporting_period) = conf.annealing_sched[tstep];
    sim->setT(T);
    for (int i = 0; i < steps; ++i) {
      sim->doStep();
      ++steps_total;
      ++steps_T;
      if (reporting_period != 0 && steps_T % reporting_period == 0) {
        printf("T %le E %le\n", sim->getT(), sim->getLattice()->getEnergy());
        sim->getLattice()->dump(sim->getT(), steps_total);
      }
    }
  }
  return 0;
}