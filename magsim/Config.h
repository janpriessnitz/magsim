#ifndef CONFIG_H
#define CONFIG_H

#include <cstdint>
#include <vector>
#include <tuple>

#include "Util.h"

typedef std::tuple<real, uint64_t, uint64_t> annealing_step_t;
typedef std::vector<annealing_step_t> annealing_sched_t;

class Config {
public:
  char method = 'H';   // 'H' for heatbath, 'M' for metropolis
  uint64_t lattice_w = 48;
  uint64_t lattice_h = 48;

  real J = 1;
  real D = 1.4;
  vec3d B = {0, 0, 0.5};
  annealing_sched_t annealing_sched;  // temperature, number of steps, reporting period
  real deltaSpin = 0.1;
  std::string lattice_dump_file;
};

#endif // UTIL_H