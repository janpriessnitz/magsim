#include "SpinDynamics.h"
#include "Constants.h"

#include <omp.h>

int seed = 123456;

SpinDynamics::SpinDynamics(SpinLattice *lattice)
  : lattice_(lattice)
{
  int max_threads = omp_get_max_threads();
  rng_engs_.resize(max_threads);
  for (int i = 0; i < max_threads; ++i) {
    rng_engs_[i] = std::mt19937(seed + i);
  }
}

SpinDynamics::~SpinDynamics() {

}

int n_step = 0;

void SpinDynamics::DoStep() {
  ++n_step;
  auto Heffs = lattice_->ComputeHeffs();

  real temp_field_mag = sqrt(2*alpha_*constants::boltzmann*temperature_/(timestep_*constants::gyromagnetic_ratio*constants::mu_B));
  auto norm_dist = std::normal_distribution<real>(0, 1);
  #pragma omp parallel for private(norm_dist)
  for (size_t i = 0; i < lattice_->spins_.size(); ++i) {
    vec3d Heff = - (1/constants::mu_B)*Heffs[i];
    vec3d temp_field = temp_field_mag*vec3d{norm_dist(rng_engs_[omp_get_thread_num()]), norm_dist(rng_engs_[omp_get_thread_num()]), norm_dist(rng_engs_[omp_get_thread_num()])};
    // printf("%s %s\n", to_string(Heff).c_str(), to_string(temp_field).c_str());
    Heff = Heff + temp_field;
    // printf("%lg %lg\n", mag(Heff), mag(TemperatureField()));

    // TODO: better ODE solver
    vec3d spin = lattice_->spins_[i];
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
