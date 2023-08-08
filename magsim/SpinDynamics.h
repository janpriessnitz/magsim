#ifndef SPINDYNAMICS_H
#define SPINDYNAMICS_H

#include "SpinLattice.h"


class SpinDynamics {
public:
  static constexpr real kGyromagneticRatio = 1;

  explicit SpinDynamics(SpinLattice *lattice, Timer & timer);
  ~SpinDynamics();

  void DoStep();
  void DoStep_Heun();
  void DoStep_Stupid();


  void GetTemperatureField();

  inline vec3d GetSpinUpdate(vec3d spin, vec3d Heff) const;

  SpinLattice* lattice_;

  real temperature_ = 300;  // K
  real alpha_ = 0.1;  // damping
  real timestep_ = 1e-15;

  std::vector<std::mt19937> rng_engs_;  // TODO: seed

  std::vector<vec3d> Heffs_;
  std::vector<vec3d> Heffs_prime_;
  std::vector<vec3d> spins_prime_;
  std::vector<vec3d> temp_field_;

  Timer & timer_;
};

#endif
