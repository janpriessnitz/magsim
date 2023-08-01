#ifndef SPINDYNAMICS_H
#define SPINDYNAMICS_H

#include "Simulation.h"
#include "Config.h"
#include "SpinLattice.h"


class SpinDynamics : public Simulation {
public:
  SpinDynamics(const Config &conf);
  ~SpinDynamics();

  void setT(real T);
  real getT();
  SpinLattice* getLattice();
  real doStep();

  SpinLattice lattice_;
  real T_;

  const Config &conf_;
};

#endif