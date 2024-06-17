#ifndef HEATBATH_H
#define HEATBATH_H

#include "Simulation.h"

class Heatbath : public Simulation {
public:
  Heatbath(const Config & config, SpinLattice *lattice);

  virtual void DoStep();

  std::vector<vec3d> Heffs_;

};

#endif