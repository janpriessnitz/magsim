
#include "Util.h"

#include <ctime>

vec3d operator+(const vec3d& a, const vec3d& b) {
  return {std::get<0>(a) + std::get<0>(b), std::get<1>(a) + std::get<1>(b), std::get<2>(a) + std::get<2>(b)};
}

vec3d& operator+=(vec3d& a, const vec3d& b) {
  a = a + b;
  return a;
}

vec3d operator-(const vec3d& a, const vec3d& b) {
  return {std::get<0>(a) - std::get<0>(b), std::get<1>(a) - std::get<1>(b), std::get<2>(a) - std::get<2>(b)};
}

vec3d& operator-=(vec3d& a, const vec3d& b) {
  a = a - b;
  return a;
}

vec3d operator*(const real& a, const vec3d& b) {
  return {a*std::get<0>(b), a*std::get<1>(b), a*std::get<2>(b)};
}

vec3d operator*(const vec3d& a, const real& b) {
  return b*a;
}

vec3d& operator*=(vec3d& a, const real& b) {
  a = b*a;
  return a;
}

vec3d operator/(const vec3d& a, const real& b) {
  return {std::get<0>(a)/b, std::get<1>(a)/b, std::get<2>(a)/b};
}

vec3d& operator/=(vec3d& a, const real& b) {
  a = a/b;
  return a;
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

extern vec3d normalize(const vec3d &vec) {
  return (1/mag(vec))*vec;
}

extern mat3d transpose(const mat3d& A) {
  auto [A1, A2, A3] = A;
  auto [A11, A12, A13] = A1;
  auto [A21, A22, A23] = A2;
  auto [A31, A32, A33] = A3;
  return {{A11, A21, A31}, {A12, A22, A32}, {A13, A23, A33}};
}

extern mat3d inverse(const mat3d& A) {
  auto [A1, A2, A3] = A;
  auto [a, b, c] = A1;
  auto [d, e, f] = A2;
  auto [g, h, i] = A3;
  real det = a*(e*i - f*h) - b*(d*i -f*g) + c*(d*h - e*g);
  return {(1/det)*vec3d{e*i - f*h, c*h - b*i, b*f - c*e},
          (1/det)*vec3d{f*g - d*i, a*i - c*g, c*d - a*f},
          (1/det)*vec3d{d*h - e*g, b*g - a*h, a*e - b*d}};
}


extern std::string to_string(const vec3d& v) {
  char buf[200];
  auto [x, y, z] = v;
  sprintf(buf, "%lg %lg %lg", x, y, z);
  return std::string(buf);
}

std::default_random_engine global_rng_eng(time(nullptr));
