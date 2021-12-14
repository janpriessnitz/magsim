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

extern real rnd_uni(real min, real max);
extern real rnd_norm(real mean, real stddev);

// from interval [min, max)
extern int64_t rnd_int(int64_t min, int64_t max);


// returns an isotropically distributed unit vector
extern vec3d rnd_vec();

extern std::default_random_engine rnd_eng;

#endif // UTIL_H