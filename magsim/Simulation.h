#ifndef SIMULATION_H
#define SIMULATION_H

#include <omp.h>
#include <chrono>

#include "SpinLattice.h"


class Simulation {
public:
  explicit Simulation(SpinLattice *lattice)
  : lattice_(lattice)
  {
    std::chrono::nanoseconds seed_ns = std::chrono::system_clock::now().time_since_epoch();
    uint32_t seed = seed_ns.count();

    int max_threads = omp_get_max_threads();
    rng_engs_.resize(max_threads);
    for (size_t i = 0; i < max_threads; ++i) {
      rng_engs_[i].seed(seed+i);
    }
  }

  virtual void DoStep() = 0;

  SpinLattice* lattice_;
  std::vector<std::mt19937> rng_engs_;
};

#endif
