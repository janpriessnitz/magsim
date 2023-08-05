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


// Heun
void SpinDynamics::DoStep() {
  DoStep_Heun();
  // DoStep_Stupid();
}

void SpinDynamics::DoStep_Heun() {
  size_t n_spins = lattice_->spins_.size();

  Heffs_.resize(n_spins);
  Heffs_prime_.resize(n_spins);
  spins_prime_.resize(n_spins);

  lattice_->ComputeHeffs(Heffs_);
  GetTemperatureField();


  #pragma omp parallel for
  for (size_t i = 0; i < n_spins; ++i) {
    vec3d Heff = (-1/constants::mu_B)*Heffs_[i];
    Heff = Heff + temp_field_[i];

    vec3d oldspin = lattice_->spins_[i];

    spins_prime_[i] = oldspin + timestep_*GetSpinUpdate(oldspin, Heff);
    spins_prime_[i] = (1/mag(spins_prime_[i]))*spins_prime_[i];
  }

  lattice_->ComputeHeffs(spins_prime_, Heffs_prime_);

  #pragma omp parallel for
  for (size_t i = 0; i < n_spins; ++i) {
    vec3d Heff = (-1/constants::mu_B)*Heffs_[i];
    Heff = Heff + temp_field_[i];

    vec3d Heff_prime = (-1/constants::mu_B)*Heffs_prime_[i];
    Heff_prime = Heff_prime + temp_field_[i];

    vec3d oldspin = lattice_->spins_[i];

    vec3d newspin = oldspin + (timestep_/2)*(GetSpinUpdate(oldspin, Heff) + GetSpinUpdate(spins_prime_[i], Heff_prime));
    lattice_->spins_[i] = normalize(newspin);
  }
}

void SpinDynamics::DoStep_Stupid() {
  size_t n_spins = lattice_->spins_.size();

  Heffs_.resize(n_spins);
  lattice_->ComputeHeffs(Heffs_);
  GetTemperatureField();

  #pragma omp parallel
  {
    #pragma omp for simd
    for (size_t i = 0; i < lattice_->spins_.size(); ++i) {
      vec3d Heff = - (1/constants::mu_B)*Heffs_[i];
      Heff = Heff + temp_field_[i];

      vec3d spin = lattice_->spins_[i];

      spin = spin + timestep_*GetSpinUpdate(spin, Heff);
      spin = (1/mag(spin))*spin;
      lattice_->spins_[i] = spin;
    }
  }
}

inline vec3d SpinDynamics::GetSpinUpdate(vec3d spin, vec3d Heff) const {
  vec3d derM = - constants::gyromagnetic_ratio/(1 + alpha_*alpha_)*(vec_prod(spin, Heff))
               - constants::gyromagnetic_ratio*alpha_/(1 + alpha_*alpha_)*(vec_prod(spin, vec_prod(spin, Heff)));
  return derM;
}

void SpinDynamics::GetTemperatureField() {
  temp_field_.resize(lattice_->spins_.size());
  real temp_field_mag = sqrt(2*alpha_*constants::boltzmann*temperature_/(timestep_*constants::gyromagnetic_ratio*constants::mu_B));
  auto norm_dist = std::normal_distribution<real>(0, 1);
  #pragma omp parallel for private(norm_dist)
  for (size_t i = 0; i < lattice_->spins_.size(); ++i) {
    temp_field_[i] = temp_field_mag*vec3d{norm_dist(rng_engs_[omp_get_thread_num()]),
                                         norm_dist(rng_engs_[omp_get_thread_num()]),
                                         norm_dist(rng_engs_[omp_get_thread_num()])};
  }
}
