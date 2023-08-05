#ifndef SPINDYNAMICS_H
#define SPINDYNAMICS_H

#include "SpinLattice.h"


class SpinDynamics {
public:
  static constexpr real kGyromagneticRatio = 1;

  explicit SpinDynamics(SpinLattice *lattice);
  ~SpinDynamics();

  void DoStep();

  SpinLattice* lattice_;

  real temperature_ = 300;  // K
  real alpha_ = 0.1;  // damping
  real timestep_ = 1e-15;

  std::vector<std::mt19937> rng_engs_;  // TODO: seed

};

#endif
