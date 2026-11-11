/* Markov chain over the 4 risk levels. Matrices are row-major N_LEVELS x N_LEVELS. */
#ifndef MARKOV_H
#define MARKOV_H

#include <stddef.h>

#include "categories.h"

/* Estimate P[i][j] = P(level j next year | level i) from the first simulated year
 * (states as filled by simulate_portfolio), so P describes one-year transitions
 * out of the current portfolio. */
void transition_matrix(const int *states, size_t n, double *p);

/* Share of the portfolio in each level in year t. */
void distribution(const int *states, size_t n, int t, double *pi);

/* Exact projection pi_t = pi_0 P^t for t = 0..n_years; out holds (n_years + 1) * N_LEVELS. */
void project(const double *pi0, const double *p, int n_years, double *out);

/* Solve pi = pi P with sum(pi) = 1. Returns 0 if the system is singular. */
int stationary_distribution(const double *p, double *pi);

#endif
