#include "Metropolis.h"

#include "Constants.h"

Metropolis::Metropolis(SpinLattice *lattice)
  : Simulation(lattice)
{
  Heffs_.resize(lattice_->NumSpins());
}

void Metropolis::DoStep() {
  size_t n_spins = lattice_->NumSpins();
  lattice_->ComputeHeffs(Heffs_);

  auto start = std::chrono::high_resolution_clock::now();

  size_t accepted = 0;
  #pragma omp parallel for reduction(+:accepted)
  for (size_t i = 0; i < n_spins; ++i) {
    vec3d Heff = Heffs_[i];
    const vec3d & old_spin = lattice_->spins_[i];
    real old_energy = scal_prod(Heff, old_spin);
    vec3d new_spin = ProposeSpin(old_spin);
    real new_energy = scal_prod(Heff, new_spin);
    real energy_diff = new_energy - old_energy;
    if (ShouldAccept(energy_diff)) {
      ++accepted;
      lattice_->spins_[i] = new_spin;
    }
  }

  last_acceptance_ratio_ = accepted/(real)n_spins;
  auto stop = std::chrono::high_resolution_clock::now();
  global_timer.AddTime(stop - start, Timer::Section::Integrator);

  if (last_acceptance_ratio_  < 0.5) {
    random_spin_mag_ *= kRandomSpinMagReductionSpeed;
  } else if (last_acceptance_ratio_ > 0.5 && random_spin_mag_ < 0.69) {
    random_spin_mag_ /= kRandomSpinMagReductionSpeed;
  }
  lattice_->SampleAverages();
}

vec3d Metropolis::ProposeSpin(const vec3d & old_spin) {
  vec3d random_spin = random_spin_mag_*rnd_gauss_vec(rng_engs_[omp_get_thread_num()]);
  vec3d new_spin = old_spin + random_spin;
  return normalize(new_spin);
}

bool Metropolis::ShouldAccept(const real & energy_diff) {
  if (energy_diff <= 0) {
    return true;
  }

  real acceptance_probability = exp(-energy_diff/(constants::boltzmann*temperature_));
  return rnd_uni(0, 1, rng_engs_[omp_get_thread_num()]) <= acceptance_probability;
}