#ifndef METROPOLIS_H
#define METROPOLIS_H

#include "Config.h"
#include "SpinLattice.h"

class Metropolis {
public:
  Metropolis(const Config &conf);

  void do_step();

  real get_probability(real deltaE);
  bool should_accept(real deltaE);

  SpinLattice lattice_;
  real T_;
  real deltaSpin_;
};


#endif  // METROPOLIS_H