#ifndef SPIN_LATTICE_H
#define SPIN_LATTICE_H

#include <vector>
#include <string>


#include "Config.h"
#include "Util.h"

// - periodic
// - Heisenberg
// - 2D squared
class SpinLattice {
public:
  SpinLattice(const Config &conf);
  ~SpinLattice();

  vec3d get(int64_t x, int64_t y);
  void set(int64_t x, int64_t y, vec3d spin);
  uint64_t index(int64_t x, int64_t y);

  real getEnergy();
  real getEnergyDelta(int64_t x, int64_t y, vec3d newspin);
  vec3d getLocalField(int64_t x, int64_t y);

  void dump(std::string filename, real T);

  void init_random();

  uint64_t w_;
  uint64_t h_;
  real J_;
  real D_;
  vec3d B_;
  std::vector<vec3d> lattice_;
};

#endif // SIMULATION_H