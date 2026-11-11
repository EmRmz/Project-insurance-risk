/* Classification of drivers into 4 risk levels by expected annual claims. */
#ifndef CATEGORIES_H
#define CATEGORIES_H

#include <stddef.h>

#define N_LEVELS 4

extern const char *const LEVELS[N_LEVELS];

/* Quantile q of x (linear interpolation between order statistics). x is not modified. */
double quantile(const double *x, size_t n, double q);

/* Cut points on lambda at the quartiles, so each level holds a quarter of the population. */
int fit_thresholds(const double *lambda, size_t n, double *thresholds);

/* Level index 0..3 (0 = Preferential) for an expected annual claim rate. */
int categorize(double lambda, const double *thresholds);

#endif
