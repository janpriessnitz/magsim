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
  auto Heffs = lattice_->ComputeHeffs();
  auto temp_field = GetTemperatureField();

  // result of Heun method's first part (predictor)
  std::vector<vec3d> spins_prime;
  spins_prime.resize(lattice_->spins_.size());

  #pragma omp parallel for
  for (size_t i = 0; i < lattice_->spins_.size(); ++i) {
    vec3d Heff = - (1/constants::mu_B)*Heffs[i];
    Heff = Heff + temp_field[i];

    vec3d oldspin = lattice_->spins_[i];

    spins_prime[i] = oldspin + timestep_*GetSpinUpdate(oldspin, Heff);
    spins_prime[i] = (1/mag(spins_prime[i]))*spins_prime[i];
  }

  auto Heffs_prime = lattice_->ComputeHeffs(spins_prime);

  #pragma omp parallel for
  for (size_t i = 0; i < lattice_->spins_.size(); ++i) {
    vec3d Heff = - (1/constants::mu_B)*Heffs[i];
    Heff = Heff + temp_field[i];

    vec3d Heff_prime = - (1/constants::mu_B)*Heffs_prime[i];
    Heff_prime = Heff_prime + temp_field[i];

    vec3d oldspin = lattice_->spins_[i];

    vec3d newspin = oldspin + (timestep_/2)*(GetSpinUpdate(oldspin, Heff) + GetSpinUpdate(spins_prime[i], Heff_prime));
    lattice_->spins_[i] = normalize(newspin);
  }
}

void SpinDynamics::DoStep_Stupid() {
  auto Heffs = lattice_->ComputeHeffs();
  auto temp_field = GetTemperatureField();

  #pragma omp parallel for
  for (size_t i = 0; i < lattice_->spins_.size(); ++i) {
    vec3d Heff = - (1/constants::mu_B)*Heffs[i];
    Heff = Heff + temp_field[i];

    vec3d spin = lattice_->spins_[i];

    spin = spin + timestep_*GetSpinUpdate(spin, Heff);
    spin = (1/mag(spin))*spin;
    lattice_->spins_[i] = spin;
  }
}

inline vec3d SpinDynamics::GetSpinUpdate(vec3d spin, vec3d Heff) const {
  vec3d derM = - constants::gyromagnetic_ratio/(1 + alpha_*alpha_)*(vec_prod(spin, Heff))
               - constants::gyromagnetic_ratio*alpha_/(1 + alpha_*alpha_)*(vec_prod(spin, vec_prod(spin, Heff)));
  return derM;
}

std::vector<vec3d> SpinDynamics::GetTemperatureField() {
  std::vector<vec3d> temp_vecs;
  temp_vecs.resize(lattice_->spins_.size());
  real temp_field_mag = sqrt(2*alpha_*constants::boltzmann*temperature_/(timestep_*constants::gyromagnetic_ratio*constants::mu_B));
  auto norm_dist = std::normal_distribution<real>(0, 1);
  #pragma omp parallel for private(norm_dist)
  for (size_t i = 0; i < lattice_->spins_.size(); ++i) {
    temp_vecs[i] = temp_field_mag*vec3d{norm_dist(rng_engs_[omp_get_thread_num()]),
                                         norm_dist(rng_engs_[omp_get_thread_num()]),
                                         norm_dist(rng_engs_[omp_get_thread_num()])};
  }
  return temp_vecs;
}
