#ifndef HEATBATH_H
#define HEATBATH_H

#include "Simulation.h"

class Heatbath : public Simulation {
public:
  Heatbath(SpinLattice *lattice);

  virtual void DoStep();

  std::vector<vec3d> Heffs_;

};

#endif