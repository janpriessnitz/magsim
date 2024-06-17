#ifndef LATTICEGENERATOR_H
#define LATTICEGENERATOR_H

#include <vector>
#include <unordered_map>
#include <optional>

#include "SpinLattice.h"
#include "Util.h"
#include "MapReader.h"
#include "Config.h"

static constexpr real tol = 0.1;

class PointLookup {
public:

  PointLookup(std::vector<vec3d> point_list);

  // std::tuple<size_t, vec3d> GetClosest(const vec3d & pos) const;
  std::optional<size_t> GetExact(const vec3d & pos) const;

  std::optional<size_t> GetCubeIndex(int x, int y, int z) const;
  std::optional<size_t> GetCubeIndex(const vec3d & pos) const;

  std::tuple<int, int, int> GetCubeCoords(const vec3d & pos) const;

  std::vector<vec3d> points_;

  vec3d origin_;
  real cube_side_;
  size_t n_x_, n_y_, n_z_;
  std::vector<std::vector<std::tuple<size_t, vec3d>>> grid_;
};

class LatticeGenerator {
public:
  LatticeGenerator(const Config & config);

  virtual SpinLattice Generate() const = 0;
  virtual std::vector<vec3d> GeneratePositions() const;
  std::pair<std::optional<size_t>, int> GetPoint(const PointLookup & lookup, const vec3d & pos) const;
  std::vector<vec3d> ApplySymmetry(const vec3d & vec, const std::vector<mat3d> & syms) const;

protected:
  int nx_, ny_, nz_;
  int periodic_x_, periodic_y_, periodic_z_;
  mat3d cell_;
  mat3d cell_inv_;
  std::vector<vec3d> spin_pos_list_;

};

#endif // LATTICEGENERATOR_H
