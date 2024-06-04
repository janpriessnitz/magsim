#ifndef METROPOLIS_H
#define METROPOLIS_H

#include "Simulation.h"

class Metropolis : public Simulation {
public:
  real kInitialRandomSpinMag = 0.02;
  real kRandomSpinMagReductionSpeed = 0.8; // value between 0, 1; governs how fast is the simulation adapting to the acceptance ratio - close to 1 means very slowly
  Metropolis(SpinLattice *lattice);

  virtual void DoStep();

  vec3d ProposeSpin(const vec3d & old_spin);
  bool ShouldAccept(const real & energy_diff);

  std::vector<vec3d> Heffs_;
  // magnitude of random spin added to the old spin when proposing a new spin in Metropolis
  real random_spin_mag_ = kInitialRandomSpinMag;
  real last_acceptance_ratio_ = 0.5;
};




#endif