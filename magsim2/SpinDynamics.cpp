#include "SpinDynamics.h"

SpinDynamics::SpinDynamics(SpinLattice *lattice)
  : lattice_(lattice)
{
}

SpinDynamics::~SpinDynamics() {

}

void SpinDynamics::DoStep() {
  auto Heffs = lattice_->ComputeHeffs();

  #pragma omp parallel for simd
  for (size_t i = 0; i < lattice_->spins_.size(); ++i) {
    vec3d Heff = Heffs[i];
    vec3d spin = lattice_->spins_[i];
    // TODO: try generating temperature field earlier in batch
    Heff = Heff + TemperatureField();

    // TODO: better ODE solver
    vec3d derM = kGyromagneticRatio/(1 + alpha_*alpha_)*(vec_prod(spin, Heff)) - kGyromagneticRatio*alpha_/(1 + alpha_*alpha_)*(vec_prod(spin, vec_prod(spin, Heff)));

    spin = spin + timestep_*derM;
    lattice_->spins_[i] = spin;
  }
}

vec3d SpinDynamics::TemperatureField() {
  return {0, 0, 0};  // TODO
}