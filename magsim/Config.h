#ifndef CONFIG_H
#define CONFIG_H

#include <cstdint>

#include "Util.h"

class Config {
public:
  uint64_t lattice_w = 48;
  uint64_t lattice_h = 48;

  real J = 1;
  real D = 1.4;
  vec3d B = {0, 0, 0.5};
  real T_init = 3;
  real T_step_ratio = 0.95;
  uint64_t T_steps = 100;
  real deltaSpin = 0.1;

  uint64_t metropolis_reporting_macrostep = 100000;
  uint64_t metropolis_equilibrium_macrosteps = 100;


  // NOT USED CURRENTLY - adaptive Metropolis equilibrium
  uint64_t equilibrium_E_sample_period = 100000;
  uint64_t equilibrium_E_sample_points = 10;
  real equilibrium_E_spread_thresh = 30;


};

#endif // UTIL_H