
#include "Metropolis.h"

#include <cmath>

Metropolis::Metropolis(const Config& conf)
  : lattice_(conf)
  , T_(conf.T_init)
  , deltaSpin_(conf.deltaSpin)
{
}

void Metropolis::do_step() {
  int64_t rx = rnd_int(0, lattice_.w_);
  int64_t ry = rnd_int(0, lattice_.h_);
  vec3d oldSpin = lattice_.get(rx, ry);
  vec3d randSpin = deltaSpin_*rnd_vec();
  vec3d newSpin = oldSpin+randSpin;
  newSpin = (1/mag(newSpin))*newSpin;
  real deltaE = lattice_.getEnergyDelta(rx, ry, newSpin);
  if (should_accept(deltaE)) {
    lattice_.set(rx, ry, newSpin);
  }
}

real Metropolis::get_probability(real deltaE) {
  return (deltaE <= 0) ? 1 : exp(-deltaE/T_);
}

bool Metropolis::should_accept(real deltaE) {
  return rnd_uni(0, 1) < get_probability(deltaE);
}
