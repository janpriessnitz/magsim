#ifndef UTIL_H
#define UTIL_H

#include <tuple>
#include <cmath>
#include <random>

typedef double real;

typedef std::tuple<real, real, real> vec3d;

extern vec3d operator+(const vec3d& a, const vec3d& b);

extern vec3d operator-(const vec3d& a, const vec3d& b);

extern vec3d operator*(const real& a, const vec3d& b);

extern real scal_prod(const vec3d& a, const vec3d& b);

extern vec3d vec_prod(const vec3d& a, const vec3d& b);

extern real mag(const vec3d &vel);

extern vec3d normalize(const vec3d &vec);

extern std::string to_string(const vec3d& v);

typedef std::tuple<vec3d, vec3d, vec3d> mat3d;

extern vec3d operator*(const vec3d& v, const mat3d& A);

extern mat3d transpose(const mat3d& A);

extern mat3d inverse(const mat3d& A);

extern real rnd_uni(real min, real max);
extern real rnd_norm(real mean, real stddev);

// from interval [min, max)
extern int64_t rnd_int(int64_t min, int64_t max);

// returns a Gaussian with unit variance
extern vec3d rnd_gauss_vec();

// returns an isotropically distributed unit vector
extern vec3d rnd_unit_vec();

extern std::default_random_engine rnd_eng;

#endif // UTIL_H
