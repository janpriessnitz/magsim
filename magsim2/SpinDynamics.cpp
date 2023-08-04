#include "SpinDynamics.h"
#include "Constants.h"

SpinDynamics::SpinDynamics(SpinLattice *lattice)
  : lattice_(lattice)
{
}

SpinDynamics::~SpinDynamics() {

}

int n_step = 0;

void SpinDynamics::DoStep() {
  ++n_step;
  // printf("tempfield: %lg\n", mag(TemperatureField()));
  auto Heffs = lattice_->ComputeHeffs();

  real temp_field_mag = sqrt(2*alpha_*constants::boltzmann*temperature_/(timestep_*constants::gyromagnetic_ratio*constants::mu_B));
  auto norm_dist = std::normal_distribution<real>(0, 1);
  std::mt19937 rng_eng(n_step);  // TODO: seed
  #pragma omp parallel for private(norm_dist, rng_eng)
  for (size_t i = 0; i < lattice_->spins_.size(); ++i) {
    vec3d Heff = - (1/constants::mu_B)*Heffs[i];
    vec3d spin = lattice_->spins_[i];
    // TODO: try generating temperature field earlier in batch
    Heff = Heff + temp_field_mag*vec3d{norm_dist(rng_eng), norm_dist(rng_eng), norm_dist(rng_eng)};
    // printf("%lg %lg\n", mag(Heff), mag(TemperatureField()));

    // TODO: better ODE solver
    vec3d derM = - constants::gyromagnetic_ratio/(1 + alpha_*alpha_)*(vec_prod(spin, Heff))
               - constants::gyromagnetic_ratio*alpha_/(1 + alpha_*alpha_)*(vec_prod(spin, vec_prod(spin, Heff)));


    spin = spin + timestep_*derM;
    spin = (1/mag(spin))*spin;
    lattice_->spins_[i] = spin;
  }
}

vec3d SpinDynamics::TemperatureField() {
  // TODO: add magnetic moment magnitude
  real K = 2*alpha_*constants::boltzmann*temperature_/(timestep_*constants::gyromagnetic_ratio*constants::mu_B);
  return sqrt(K)*rnd_gauss_vec();
}
