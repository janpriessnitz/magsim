
#include "Util.h"

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
  auto [a, b, c] = vel;
  return sqrt(a * a + b * b + c * c);
}

vec3d operator*(const vec3d& v, const mat3d& A) {
  auto At = transpose(A);
  auto [A1, A2, A3] = At;
  return {scal_prod(v, A1), scal_prod(v, A2), scal_prod(v, A3)};
}

extern mat3d transpose(const mat3d& A) {
  auto [A1, A2, A3] = A;
  auto [A11, A12, A13] = A1;
  auto [A21, A22, A23] = A2;
  auto [A31, A32, A33] = A3;
  return {{A11, A21, A31}, {A12, A22, A32}, {A13, A23, A33}};
}

extern std::string to_string(const vec3d& v) {
  char buf[200];
  auto [x, y, z] = v;
  sprintf(buf, "%lg %lg %lg", x, y, z);
  return std::string(buf);
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

// TODO: not divide by sqrt(3) ?? to get unit variance
vec3d rnd_gauss_vec() {
  real x, y, z;
  x = rnd_norm(0, 1);
  y = rnd_norm(0, 1);
  z = rnd_norm(0, 1);
  vec3d vec = {x, y, z};
  return vec;
}

vec3d rnd_unit_vec() {
  vec3d vec = rnd_gauss_vec();
  return (1/mag(vec))*vec;
}

std::default_random_engine rnd_eng(time(nullptr));
