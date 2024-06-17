#ifndef SPINDYNAMICS_H
#define SPINDYNAMICS_H

#include "Simulation.h"
#include "SpinLattice.h"


class SpinDynamics : public Simulation {
public:
  static constexpr real kGyromagneticRatio = 1;

  explicit SpinDynamics(const Config & config, SpinLattice *lattice);

  virtual void DoStep();

  void DoStep_Heun();
  void DoStep_Stupid();

  void GetTemperatureField();

  inline vec3d GetSpinUpdate(vec3d spin, vec3d Heff) const;

  bool CheckTimestep();

  real alpha_ = 0.1;  // damping
  real timestep_ = 1e-15;

  bool check_timestep_;

  std::vector<vec3d> Heffs_;
  std::vector<vec3d> Heffs_prime_;
  std::vector<vec3d> spins_prime_;
  std::vector<vec3d> temp_field_;
};

#endif
