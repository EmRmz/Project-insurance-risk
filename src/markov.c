#include "markov.h"

#include <string.h>

#include "linalg.h"

#define K N_LEVELS

void transition_matrix(const int *states, size_t n, double *p)
{
    double counts[K * K] = {0};
    for (size_t i = 0; i < n; i++)
        counts[states[i] * K + states[n + i]] += 1;
    for (int a = 0; a < K; a++) {
        double row = 0;
        for (int b = 0; b < K; b++)
            row += counts[a * K + b];
        for (int b = 0; b < K; b++)
            p[a * K + b] = row > 0 ? counts[a * K + b] / row : (a == b);
    }
}

void distribution(const int *states, size_t n, int t, double *pi)
{
    memset(pi, 0, K * sizeof *pi);
    for (size_t i = 0; i < n; i++)
        pi[states[(size_t)t * n + i]] += 1.0 / (double)n;
}

void project(const double *pi0, const double *p, int n_years, double *out)
{
    memcpy(out, pi0, K * sizeof *out);
    for (int t = 1; t <= n_years; t++) {
        const double *prev = out + (t - 1) * K;
        double *cur = out + t * K;
        for (int j = 0; j < K; j++) {
            cur[j] = 0;
            for (int i = 0; i < K; i++)
                cur[j] += prev[i] * p[i * K + j];
        }
    }
}

/* pi (P - I) = 0 is the system (P' - I) pi' = 0; its last equation is redundant,
 * so it is replaced by the normalization sum(pi) = 1. */
int stationary_distribution(const double *p, double *pi)
{
    double a[K * K];
    for (int i = 0; i < K; i++)
        for (int j = 0; j < K; j++)
            a[i * K + j] = p[j * K + i] - (i == j);
    for (int j = 0; j < K; j++) {
        a[(K - 1) * K + j] = 1;
        pi[j] = 0;
    }
    pi[K - 1] = 1;
    return gauss_solve(a, pi, K);
}
