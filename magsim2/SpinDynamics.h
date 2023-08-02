#ifndef SPINDYNAMICS_H
#define SPINDYNAMICS_H

#include "SpinLattice.h"


class SpinDynamics {
public:
  static constexpr real kGyromagneticRatio = 1;

  explicit SpinDynamics(SpinLattice *lattice);
  ~SpinDynamics();

  void DoStep();

  vec3d TemperatureField();

  SpinLattice* lattice_;

  real temperature_;
  real alpha_ = 0.1;
  // real timestep_ = 1e-16;
  real timestep_ = 1e-03;
};

#endif