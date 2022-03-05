#ifndef HEATBATH_H
#define HEATBATH_H

#include "Config.h"
#include "SpinLattice.h"


class HeatBath {
public:
  HeatBath(const Config &conf);
  ~HeatBath();

  real do_step();
  void equilibrize();

  SpinLattice lattice_;
  real T_;

  const Config &conf_;
  // uint64_t equilibrium_E_sample_period_;
  // uint64_t equilibrium_E_sample_points_;
  // real equilibrium_E_spread_thresh_;

  // FILE *E_history_fp_;
};

#endif