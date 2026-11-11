#include "categories.h"

#include <stdlib.h>
#include <string.h>

const char *const LEVELS[N_LEVELS] = {"Preferential", "Standard", "Moderate Risk", "High Risk"};

static int compare(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

static double sorted_quantile(const double *s, size_t n, double q)
{
    double pos = q * (double)(n - 1);
    size_t lo = (size_t)pos;
    if (lo + 1 >= n)
        return s[n - 1];
    return s[lo] + (pos - (double)lo) * (s[lo + 1] - s[lo]);
}

double quantile(const double *x, size_t n, double q)
{
    double *s = malloc(n * sizeof *s);
    if (!s)
        return 0;
    memcpy(s, x, n * sizeof *s);
    qsort(s, n, sizeof *s, compare);
    double v = sorted_quantile(s, n, q);
    free(s);
    return v;
}

int fit_thresholds(const double *lambda, size_t n, double *thresholds)
{
    double *s = malloc(n * sizeof *s);
    if (!s)
        return 0;
    memcpy(s, lambda, n * sizeof *s);
    qsort(s, n, sizeof *s, compare);
    for (int t = 0; t < N_LEVELS - 1; t++)
        thresholds[t] = sorted_quantile(s, n, (t + 1) / (double)N_LEVELS);
    free(s);
    return 1;
}

int categorize(double lambda, const double *thresholds)
{
    int level = 0;
    while (level < N_LEVELS - 1 && lambda >= thresholds[level])
        level++;
    return level;
}
