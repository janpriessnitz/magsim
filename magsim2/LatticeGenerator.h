#ifndef LATTICEGENERATOR_H
#define LATTICEGENERATOR_H

#include <vector>
#include <unordered_map>
#include <optional>

#include "SpinLattice.h"
#include "Util.h"


class PointLookup {
public:
  static constexpr real tol = 0.0001;

  PointLookup(std::vector<vec3d> point_list);

  // std::tuple<size_t, vec3d> GetClosest(const vec3d & pos) const;
  std::optional<size_t> GetExact(const vec3d & pos) const;

  size_t GetCubeIndex(int x, int y, int z) const;
  size_t GetCubeIndex(const vec3d & pos) const;

  std::tuple<int, int, int> GetCubeCoords(const vec3d & pos) const;

  std::vector<vec3d> points_;

  vec3d origin_;
  real cube_side_;
  size_t n_x_, n_y_, n_z_;
  std::vector<std::vector<std::tuple<size_t, vec3d>>> grid_;
};

class LatticeGenerator {
public:
  LatticeGenerator();
};

class HcpCobaltGenerator {
public:
  HcpCobaltGenerator();

  SpinLattice Generate() const;

  std::vector<vec3d> GeneratePositions() const;
  std::vector<std::vector<std::tuple<size_t, real>>> GenerateExchange(const PointLookup & point_lookup) const;
  std::vector<vec3d> GenerateSpins(std::vector<vec3d> positions) const;

};




#endif // LATTICEGENERATOR_H