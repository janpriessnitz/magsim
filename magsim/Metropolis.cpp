
#include "Metropolis.h"

#include <cmath>

Metropolis::Metropolis(const Config& conf)
  : lattice_(conf)
  , T_(conf.T_init)
  , deltaSpin_(conf.deltaSpin)
  , equilibrium_E_sample_period_(conf.equilibrium_E_sample_period)
  , equilibrium_E_sample_points_(conf.equilibrium_E_sample_points)
  , equilibrium_E_spread_thresh_(conf.equilibrium_E_spread_thresh)
{
  E_history_fp_ = fopen("Ehistory.out", "w");
}

Metropolis::~Metropolis() {
  fclose(E_history_fp_);
}

real Metropolis::do_step() {
  int64_t rx = rnd_int(0, lattice_.w_);
  int64_t ry = rnd_int(0, lattice_.h_);
  vec3d oldSpin = lattice_.get(rx, ry);
  vec3d randSpin = deltaSpin_*rnd_vec();
  vec3d newSpin = oldSpin+randSpin;
  newSpin = (1/mag(newSpin))*newSpin;
  real deltaE = lattice_.getEnergyDelta(rx, ry, newSpin);
  if (should_accept(deltaE)) {
    lattice_.set(rx, ry, newSpin);
    return deltaE;
  }
  return 0;
}

void Metropolis::equilibrize() {
  std::vector<real> E_history(equilibrium_E_sample_points_);
  E_history[equilibrium_E_sample_points_ - 1] = 1e9;  // hack
  for(int i = 0; ; i = (i + 1) % equilibrium_E_sample_points_) {
    for(int j = 0; j < equilibrium_E_sample_period_; ++j) {
      do_step();
    }
    real curE = lattice_.getEnergy();
    fprintf(E_history_fp_, "%lf %lf\n", T_, curE);
    E_history[i] = curE;
    printf("E %lf, T %lf\n", lattice_.getEnergy(), T_);
    if (is_in_equilibrium(E_history)) return;
  }
}

bool Metropolis::is_in_equilibrium(const std::vector<real>& Es) {
  real Emax = -1e9, Emin = 1e9;  // hacks
  for (int i = 0; i < equilibrium_E_sample_points_; ++i) {
    if (Es[i] > Emax) Emax = Es[i];
    if (Es[i] < Emin) Emin = Es[i];
  }
  real Espread = Emax - Emin;
  return Espread < equilibrium_E_spread_thresh_;
}


real Metropolis::get_probability(real deltaE) {
  return (deltaE <= 0) ? 1 : exp(-deltaE/T_);
}

bool Metropolis::should_accept(real deltaE) {
  return rnd_uni(0, 1) < get_probability(deltaE);
}
