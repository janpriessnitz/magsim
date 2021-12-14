#ifndef CONFIG_H
#define CONFIG_H

#include <cstdint>

#include "Util.h"

class Config {
public:
  uint64_t lattice_w = 48;
  uint64_t lattice_h = 48;

  real J = 1;
  vec3d D = {0, 0, 1.5};
  vec3d B = {0, 0, 0.2};
  real T_init = 3;
  real deltaSpin = 0.04;

};

#endif // UTIL_H