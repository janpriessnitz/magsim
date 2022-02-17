#ifndef METROPOLIS_H
#define METROPOLIS_H

#include "Config.h"
#include "SpinLattice.h"

class Metropolis {
public:
  Metropolis(const Config &conf);
  ~Metropolis();

  real do_step();
  void equilibrize();

  bool is_in_equilibrium(const std::vector<real>& Es);

  real get_probability(real deltaE);
  bool should_accept(real deltaE);

  SpinLattice lattice_;
  real T_;
  real deltaSpin_;
  uint64_t equilibrium_E_sample_period_;
  uint64_t equilibrium_E_sample_points_;
  real equilibrium_E_spread_thresh_;

  FILE *E_history_fp_;
};


#endif  // METROPOLIS_H