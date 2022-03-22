
#include "Metropolis.h"

#include <cmath>

Metropolis::Metropolis(const Config& conf)
  : lattice_(conf)
  , conf_(conf)
{
  E_history_fp_ = fopen("Ehistory.out", "w");
}

Metropolis::~Metropolis() {
  fclose(E_history_fp_);
}

void Metropolis::setT(real T) {
  T_ = T;
}

real Metropolis::getT() {
  return T_;
}

SpinLattice *Metropolis::getLattice() {
  return &lattice_;
}

real Metropolis::doStep() {
  int64_t rx = rnd_int(0, lattice_.w_);
  int64_t ry = rnd_int(0, lattice_.h_);
  vec3d oldSpin = lattice_.get(rx, ry);
  vec3d randSpin = conf_.deltaSpin*rnd_vec();
  vec3d newSpin = oldSpin+randSpin;
  newSpin = (1/mag(newSpin))*newSpin;
  real deltaE = lattice_.getEnergyDelta(rx, ry, newSpin);
  if (should_accept(deltaE)) {
    lattice_.set(rx, ry, newSpin);
    return deltaE;
  }
  return 0;
}

real Metropolis::get_probability(real deltaE) {
  return (deltaE <= 0) ? 1 : exp(-deltaE/T_);
}

bool Metropolis::should_accept(real deltaE) {
  return rnd_uni(0, 1) < get_probability(deltaE);
}
