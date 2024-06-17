#ifndef SIMULATION_H
#define SIMULATION_H

#include <omp.h>
#include <chrono>

#include "SpinLattice.h"


class Simulation {
public:
  explicit Simulation(const Config & config, SpinLattice *lattice)
  : lattice_(lattice)
  {
    temperature_ = config.Get<real>("temperature");
    out_dir_ = config.out_dir_;
    profile_direction_ = config.Get<std::string>("profile_direction")[0];
    std::chrono::nanoseconds seed_ns = std::chrono::system_clock::now().time_since_epoch();
    uint32_t seed = seed_ns.count();

    int max_threads = omp_get_max_threads();
    rng_engs_.resize(max_threads);
    for (size_t i = 0; i < max_threads; ++i) {
      rng_engs_[i].seed(seed+i);
    }

    num_step_ = config.Get<size_t>("num_step");
    num_substep_ = config.Get<size_t>("num_substep");
  }

  virtual void DoStep() = 0;
  virtual void Run() {
    auto start = std::chrono::high_resolution_clock::now();
    for (int j = 0; j < num_step_; ++j) {
      for (int i = 0; i < num_substep_; ++i) {
        DoStep();
      }
      lattice_->PrintEnergy();
      auto avgm = lattice_->AvgM();
      printf("%s %lf\n", to_string(avgm).c_str(), mag(avgm));
      bool dump_avgs = true;
      lattice_->DumpLattice(out_dir_ + "/lattice.out" + std::to_string(j), dump_avgs);
      lattice_->DumpProfile(out_dir_ + "/profile.out" + std::to_string(j), profile_direction_, dump_avgs);
      lattice_->ResetAverages();
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    printf("took %lu ms\n", duration.count());
    global_timer.PrintStatistics();
  }

  real temperature_;
  SpinLattice* lattice_;
  std::vector<std::mt19937_64> rng_engs_;
  size_t num_step_;
  size_t num_substep_;
  char profile_direction_;
  std::string out_dir_;
};

#endif
