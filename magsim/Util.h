#ifndef UTIL_H
#define UTIL_H

#include <tuple>
#include <cmath>
#include <random>

typedef double real;

typedef std::tuple<real, real, real> vec3d;

extern vec3d operator+(const vec3d& a, const vec3d& b);
extern vec3d& operator+=(vec3d& a, const vec3d& b);

extern vec3d operator-(const vec3d& a, const vec3d& b);
extern vec3d& operator-=(vec3d& a, const vec3d& b);

extern vec3d operator*(const real& a, const vec3d& b);
extern vec3d operator*(const vec3d& a, const real& b);
extern vec3d& operator*=(real& a, const vec3d& b);

extern vec3d operator/(const vec3d& a, const real& b);
extern vec3d& operator/=(vec3d& a, const real& b);

extern real scal_prod(const vec3d& a, const vec3d& b);

extern vec3d vec_prod(const vec3d& a, const vec3d& b);

extern real mag(const vec3d &vel);

extern vec3d normalize(const vec3d &vec);

extern std::string to_string(const vec3d& v);

typedef std::tuple<vec3d, vec3d, vec3d> mat3d;

extern vec3d operator*(const vec3d& v, const mat3d& A);

extern mat3d transpose(const mat3d& A);

extern mat3d inverse(const mat3d& A);

extern std::default_random_engine global_rng_eng;

template<typename T = std::default_random_engine>
real rnd_uni(real min, real max, T & rng_eng = global_rng_eng) {
  return std::uniform_real_distribution<>(min, max)(rng_eng);
}

template<typename T = std::default_random_engine>
real rnd_norm(real mean, real stddev, T & rng_eng = global_rng_eng) {
  return std::normal_distribution<real>(mean, stddev)(rng_eng);
}

template<typename T = std::default_random_engine>
int64_t rnd_int(int64_t min, int64_t max, T & rng_eng = global_rng_eng) {
  return std::uniform_int_distribution<>(min, max-1)(rng_eng);
}

// TODO: not divide by sqrt(3) ?? to get unit variance
// returns a Gaussian with unit variance
template<typename T = std::default_random_engine>
vec3d rnd_gauss_vec(T & rng_eng = global_rng_eng) {
  real x, y, z;
  x = rnd_norm(0, 1, rng_eng);
  y = rnd_norm(0, 1, rng_eng);
  z = rnd_norm(0, 1, rng_eng);
  vec3d vec = {x, y, z};
  return vec;
}

// returns an isotropically distributed unit vector
template<typename T = std::default_random_engine>
vec3d rnd_unit_vec(T & rng_eng = global_rng_eng) {
  vec3d vec = rnd_gauss_vec(rng_eng);
  return (1/mag(vec))*vec;
}


#endif // UTIL_H
