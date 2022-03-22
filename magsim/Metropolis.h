#ifndef METROPOLIS_H
#define METROPOLIS_H

#include "Simulation.h"
#include "Config.h"
#include "SpinLattice.h"

class Metropolis : public Simulation {
public:
  Metropolis(const Config &conf);
  ~Metropolis();

  real doStep();
  SpinLattice* getLattice();
  void setT(real T);
  real getT();


  real get_probability(real deltaE);
  bool should_accept(real deltaE);

  SpinLattice lattice_;
  real T_;

  const Config &conf_;
  real deltaSpin_;
  uint64_t equilibrium_E_sample_period_;
  uint64_t equilibrium_E_sample_points_;
  real equilibrium_E_spread_thresh_;

  FILE *E_history_fp_;
};


#endif  // METROPOLIS_H