#include "rng.h"

#include <math.h>

static uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

static uint64_t next(Rng *r)
{
    uint64_t *s = r->s;
    uint64_t result = rotl(s[1] * 5, 7) * 9;
    uint64_t t = s[1] << 17;
    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];
    s[2] ^= t;
    s[3] = rotl(s[3], 45);
    return result;
}

/* The state is filled with splitmix64, as recommended by the xoshiro authors. */
void rng_seed(Rng *r, uint64_t seed)
{
    for (int i = 0; i < 4; i++) {
        uint64_t z = (seed += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        r->s[i] = z ^ (z >> 31);
    }
}

double rng_uniform(Rng *r) { return (double)(next(r) >> 11) * 0x1.0p-53; }

uint64_t rng_below(Rng *r, uint64_t n) { return (uint64_t)(rng_uniform(r) * (double)n); }

/* Knuth's multiplication method: count uniforms until their product drops
 * below e^-lambda. Large rates are split using Poisson(a + b) = Poisson(a) + Poisson(b). */
int rng_poisson(Rng *r, double lambda)
{
    int n = 0;
    while (lambda > 30.0) {
        n += rng_poisson(r, 30.0);
        lambda -= 30.0;
    }
    double limit = exp(-lambda), p = rng_uniform(r);
    while (p > limit) {
        n++;
        p *= rng_uniform(r);
    }
    return n;
}

/* Box-Muller transform. */
double rng_normal(Rng *r)
{
    double u = 1 - rng_uniform(r), v = rng_uniform(r);
    return sqrt(-2 * log(u)) * cos(2 * 3.14159265358979323846 * v);
}

/* Marsaglia & Tsang (2000) squeeze method. For shape < 1 it samples Gamma(shape + 1)
 * and scales by U^(1/shape). */
double rng_gamma(Rng *r, double shape)
{
    if (shape < 1)
        return rng_gamma(r, shape + 1) * pow(1 - rng_uniform(r), 1 / shape);
    double d = shape - 1.0 / 3, c = 1 / sqrt(9 * d);
    for (;;) {
        double z = rng_normal(r), v = 1 + c * z;
        if (v <= 0)
            continue;
        v = v * v * v;
        double u = 1 - rng_uniform(r);
        if (log(u) < 0.5 * z * z + d - d * v + d * log(v))
            return d * v;
    }
}
