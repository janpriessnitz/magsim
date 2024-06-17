#include "SpinDynamics.h"
#include "Constants.h"

#include <omp.h>
#include <chrono>

SpinDynamics::SpinDynamics(const Config & config, SpinLattice *lattice)
  : Simulation(config, lattice)
{
  alpha_ = config.Get<real>("damping");
  timestep_ = config.Get<real>("timestep");
  check_timestep_ = config.Get<int>("check_timestep", 0) != 0;

  size_t n_spins = lattice_->NumSpins();
  Heffs_.resize(n_spins);
  temp_field_.resize(n_spins);
}

// Heun
void SpinDynamics::DoStep() {
  DoStep_Heun();
  // DoStep_Stupid();
  lattice_->SampleAverages();
}

void SpinDynamics::DoStep_Heun() {
  size_t n_spins = lattice_->NumSpins();

  Heffs_prime_.resize(n_spins);
  spins_prime_.resize(n_spins);

  lattice_->ComputeHeffs(Heffs_);
  GetTemperatureField();

  auto start = std::chrono::high_resolution_clock::now();

  #pragma omp parallel for
  for (size_t i = 0; i < n_spins; ++i) {
    vec3d Heff = (-1/constants::mu_B)*Heffs_[i];
    Heff = Heff + temp_field_[i];

    vec3d oldspin = lattice_->spins_[i];

    spins_prime_[i] = oldspin + timestep_*GetSpinUpdate(oldspin, Heff);
    spins_prime_[i] = (1/mag(spins_prime_[i]))*spins_prime_[i];
  }
  auto stop = std::chrono::high_resolution_clock::now();
  global_timer.AddTime(stop - start, Timer::Section::Integrator);

  lattice_->ComputeHeffs(spins_prime_, Heffs_prime_);

  start = std::chrono::high_resolution_clock::now();
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

  stop = std::chrono::high_resolution_clock::now();
  global_timer.AddTime(stop - start, Timer::Section::Integrator);
}

void SpinDynamics::DoStep_Stupid() {
  size_t n_spins = lattice_->NumSpins();

  Heffs_.resize(n_spins);
  lattice_->ComputeHeffs(Heffs_);
  GetTemperatureField();

  #pragma omp parallel
  {
    #pragma omp for simd
    for (size_t i = 0; i < lattice_->NumSpins(); ++i) {
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
  auto start = std::chrono::high_resolution_clock::now();

  real temp_field_mag = sqrt(2*alpha_*constants::boltzmann*temperature_/(timestep_*constants::gyromagnetic_ratio*constants::mu_B));
  #pragma omp parallel for
  for (size_t i = 0; i < temp_field_.size(); ++i) {
    auto rnd_vec = rnd_gauss_vec(rng_engs_[omp_get_thread_num()]);
    temp_field_[i] = temp_field_mag*rnd_vec;
  }
  auto stop = std::chrono::high_resolution_clock::now();
  global_timer.AddTime(stop - start, Timer::Section::Temperature);
}

bool SpinDynamics::CheckTimestep() {
  size_t n_spins = lattice_->NumSpins();

  Heffs_prime_.resize(n_spins);
  spins_prime_.resize(n_spins);

  lattice_->ComputeHeffs(Heffs_);
  GetTemperatureField();

  #pragma omp parallel for
  for (size_t i = 0; i < n_spins; ++i) {
    vec3d Heff = (-1/constants::mu_B)*Heffs_[i];
    Heff = Heff + temp_field_[i];

    vec3d oldspin = lattice_->spins_[i];
    vec3d spin_update = timestep_*GetSpinUpdate(oldspin, Heff);
    real spin_update_mag = mag(spin_update);
  }
  // TODO
  return true;
}