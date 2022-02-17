#ifndef CONFIG_H
#define CONFIG_H

#include <cstdint>

#include "Util.h"

class Config {
public:
  uint64_t lattice_w = 48;
  uint64_t lattice_h = 48;

  real J = 1;
  real D = 1.5;
  vec3d B = {0, 0, 0.2};
  real T_init = 3;
  real deltaSpin = 0.04;

  uint64_t equilibrium_E_sample_period = 5000;
  uint64_t equilibrium_E_sample_points = 100;
  real equilibrium_E_spread_thresh = 30;


};

#endif // UTIL_H