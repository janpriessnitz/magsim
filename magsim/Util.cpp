
#include "Util.h"

#include <ctime>

vec3d operator+(const vec3d& a, const vec3d& b) {
  return {std::get<0>(a) + std::get<0>(b), std::get<1>(a) + std::get<1>(b), std::get<2>(a) + std::get<2>(b)};
}

vec3d operator-(const vec3d& a, const vec3d& b) {
  return {std::get<0>(a) - std::get<0>(b), std::get<1>(a) - std::get<1>(b), std::get<2>(a) - std::get<2>(b)};
}

vec3d operator*(const real& a, const vec3d& b) {
  return {a*std::get<0>(b), a*std::get<1>(b), a*std::get<2>(b)};
}

extern real scal_prod(const vec3d& a, const vec3d& b) {
  return std::get<0>(a)*std::get<0>(b) + std::get<1>(a)*std::get<1>(b) + std::get<2>(a)*std::get<2>(b);
}

extern vec3d vec_prod(const vec3d& a, const vec3d& b) {
  return {std::get<1>(a)*std::get<2>(b) - std::get<1>(b)*std::get<2>(a),
          std::get<2>(a)*std::get<0>(b) - std::get<2>(b)*std::get<0>(a),
          std::get<0>(a)*std::get<1>(b) - std::get<0>(b)*std::get<1>(a)};
}

real mag(const vec3d &vel) {
  real a, b, c;
  std::tie(a, b, c) = vel;
  return sqrt(a * a + b * b + c * c);
}

real rnd_uni(real min, real max) {
  return std::uniform_real_distribution<>(min, max)(rnd_eng);
}

real rnd_norm(real mean, real stddev) {
  return std::normal_distribution<real>(mean, stddev)(rnd_eng);
}

int64_t rnd_int(int64_t min, int64_t max) {
  return std::uniform_int_distribution<>(min, max-1)(rnd_eng);
}

vec3d rnd_vec() {
    real x, y, z;
    x = rnd_norm(0, 1);
    y = rnd_norm(0, 1);
    z = rnd_norm(0, 1);
    vec3d spin = {x, y, z};
    spin = (1/mag(spin))*spin;
    return spin;
}

std::default_random_engine rnd_eng(time(nullptr));
