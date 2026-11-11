/* Pseudo-random number generation: xoshiro256** (Blackman & Vigna, 2018). */
#ifndef RNG_H
#define RNG_H

#include <stdint.h>

typedef struct {
    uint64_t s[4];
} Rng;

void rng_seed(Rng *r, uint64_t seed);
double rng_uniform(Rng *r);              /* uniform in [0, 1) */
uint64_t rng_below(Rng *r, uint64_t n);  /* uniform integer in [0, n) */
int rng_poisson(Rng *r, double lambda);
double rng_normal(Rng *r);               /* standard normal */
double rng_gamma(Rng *r, double shape);  /* Gamma(shape, scale 1), shape > 0 */

#endif
