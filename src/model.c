#include "model.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "categories.h"
#include "linalg.h"

/* Knots of the age spline. Term b is the number of decades of the driver's life spent
 * between knots b and b + 1, so its coefficient is the slope of log lambda per decade
 * in that segment and the intercept is the log rate at age 18. Vehicle age 1-10 and
 * power class 4 are the reference levels of the categorical features. */
static const int AGE_KNOTS[N_AGE_SEGMENTS + 1] = {18, 20, 25, 30, 40, 50, 60, 70, AGE_MAX};

const char *const FEATURE_NAMES[N_FEATURES] = {
    "intercept", "age_slope_18_20", "age_slope_20_25", "age_slope_25_30", "age_slope_30_40",
    "age_slope_40_50", "age_slope_50_60", "age_slope_60_70", "age_slope_70_90", "veh_age_new", "veh_age_11plus",
    "power_5", "power_6", "power_7", "power_8", "power_9", "bonus_malus", "log_density",
    "diesel", "exposure_under_3m", "exposure_3_to_9m",
};

void features(const Driver *d, double *x)
{
    int k = 0;
    int power = d->veh_power < 9 ? d->veh_power : 9;
    double age = d->age < AGE_MAX ? d->age : AGE_MAX;
    x[k++] = 1.0;
    for (int b = 0; b < N_AGE_SEGMENTS; b++) {
        double hi = age < AGE_KNOTS[b + 1] ? age : AGE_KNOTS[b + 1];
        x[k++] = hi > AGE_KNOTS[b] ? (hi - AGE_KNOTS[b]) / 10 : 0;
    }
    x[k++] = d->veh_age == 0;
    x[k++] = d->veh_age > 10;
    for (int p = 5; p <= 9; p++)
        x[k++] = power == p;
    /* The model saturates at 150, as in the training data. */
    x[k++] = (d->bonus_malus < 150 ? d->bonus_malus : 150) / 100.0;
    x[k++] = log(d->density);
    x[k++] = d->diesel;
    x[k++] = d->exposure < 0.25;
    x[k++] = d->exposure >= 0.25 && d->exposure < 0.75;
}

static double linear(const double *coef, const double *x)
{
    double eta = 0;
    for (int j = 0; j < N_FEATURES; j++)
        eta += coef[j] * x[j];
    return eta;
}

double model_expected_claims(const Model *m, const Driver *d)
{
    double x[N_FEATURES];
    features(d, x);
    return d->exposure * exp(linear(m->coef, x));
}

double model_lambda(const Model *m, const Driver *d)
{
    double x[N_FEATURES];
    features(d, x);
    x[FEATURE_EXPOSURE] = x[FEATURE_EXPOSURE + 1] = 0;  /* a full year */
    return exp(linear(m->coef, x));
}

/* Accumulate the Hessian X'WX and the gradient X'(y - mu) of the Poisson log-likelihood. */
static void accumulate(const Model *m, const double *beta, const Driver *d, size_t n,
                       const unsigned char *use, double *h, double *g)
{
    double x[N_FEATURES];
    memset(h, 0, sizeof(double) * N_FEATURES * N_FEATURES);
    memset(g, 0, sizeof(double) * N_FEATURES);
    for (size_t i = 0; i < n; i++) {
        if (use && !use[i])
            continue;
        features(&d[i], x);
        for (int j = 0; j < N_FEATURES; j++)
            if (m->excluded[j])
                x[j] = 0;
        double mu = d[i].exposure * exp(linear(beta, x));
        for (int a = 0; a < N_FEATURES; a++) {
            if (x[a] == 0)
                continue;
            g[a] += (d[i].claims - mu) * x[a];
            for (int b = 0; b <= a; b++)
                h[a * N_FEATURES + b] += mu * x[a] * x[b];
        }
    }
    for (int a = 0; a < N_FEATURES; a++)
        for (int b = a + 1; b < N_FEATURES; b++)
            h[a * N_FEATURES + b] = h[b * N_FEATURES + a];
    /* An excluded feature gets a unit diagonal and zero gradient, so its step is always 0. */
    for (int a = 0; a < N_FEATURES; a++)
        if (m->excluded[a])
            h[a * N_FEATURES + a] = 1;
}

int model_fit(Model *m, const Driver *d, size_t n, const unsigned char *use)
{
    double h[N_FEATURES * N_FEATURES], g[N_FEATURES];
    double claims = 0, exposure = 0;
    for (size_t i = 0; i < n; i++) {
        if (use && !use[i])
            continue;
        claims += d[i].claims;
        exposure += d[i].exposure;
    }
    memset(m->coef, 0, sizeof m->coef);
    m->coef[0] = log(claims / exposure);

    for (int iter = 0; iter < 50; iter++) {
        accumulate(m, m->coef, d, n, use, h, g);
        if (!cholesky(h, N_FEATURES))
            return 0;
        cholesky_solve(h, N_FEATURES, g);   /* g becomes the Newton step */
        double largest = 0;
        for (int j = 0; j < N_FEATURES; j++) {
            m->coef[j] += g[j];
            largest = fmax(largest, fabs(g[j]));
        }
        if (largest < 1e-8)
            break;
    }

    /* Standard errors: square roots of the diagonal of (X'WX)^-1. */
    accumulate(m, m->coef, d, n, use, h, g);
    if (!cholesky(h, N_FEATURES))
        return 0;
    for (int j = 0; j < N_FEATURES; j++) {
        if (m->excluded[j]) {
            m->std_err[j] = 0;
            continue;
        }
        double e[N_FEATURES] = {0};
        e[j] = 1;
        cholesky_solve(h, N_FEATURES, e);
        m->std_err[j] = sqrt(e[j]);
    }
    return 1;
}

double poisson_deviance(const Model *m, double rate, const Driver *d, size_t n,
                        const unsigned char *use)
{
    double sum = 0;
    size_t count = 0;
    for (size_t i = 0; i < n; i++) {
        if (use && !use[i])
            continue;
        double mu = rate < 0 ? model_expected_claims(m, &d[i]) : d[i].exposure * rate;
        double y = d[i].claims;
        sum += (y > 0 ? y * log(y / mu) : 0) - (y - mu);
        count++;
    }
    return 2 * sum / count;
}

double model_dispersion(const Model *m, const Driver *d, size_t n, const unsigned char *use)
{
    double num = 0, den = 0;
    for (size_t i = 0; i < n; i++) {
        if (use && !use[i])
            continue;
        double mu = model_expected_claims(m, &d[i]), y = d[i].claims;
        num += (y - mu) * (y - mu) - y;
        den += mu * mu;
    }
    return num > 0 ? num / den : 0;
}

int model_save(const char *path, const Model *m, const double *thresholds)
{
    FILE *f = fopen(path, "w");
    if (!f)
        return 0;
    fprintf(f, "# Poisson GLM: feature coefficient\n");
    for (int j = 0; j < N_FEATURES; j++)
        fprintf(f, "%s %.17g\n", FEATURE_NAMES[j], m->coef[j]);
    fprintf(f, "dispersion %.17g\n", m->dispersion);
    fprintf(f, "# Risk level thresholds on lambda\nthresholds");
    for (int t = 0; t < N_LEVELS - 1; t++)
        fprintf(f, " %.17g", thresholds[t]);
    fprintf(f, "\n");
    return fclose(f) == 0;
}

int model_load(const char *path, Model *m, double *thresholds)
{
    FILE *f = fopen(path, "r");
    if (!f)
        return 0;
    char line[256], name[64];
    int found = 0;
    memset(m, 0, sizeof *m);
    while (fgets(line, sizeof line, f)) {
        double v[N_LEVELS - 1];
        if (line[0] == '#')
            continue;
        if (sscanf(line, "thresholds %lf %lf %lf", &v[0], &v[1], &v[2]) == N_LEVELS - 1) {
            memcpy(thresholds, v, sizeof v);
            found++;
            continue;
        }
        if (sscanf(line, "%63s %lf", name, &v[0]) != 2)
            continue;
        if (strcmp(name, "dispersion") == 0)
            m->dispersion = v[0];
        for (int j = 0; j < N_FEATURES; j++)
            if (strcmp(name, FEATURE_NAMES[j]) == 0) {
                m->coef[j] = v[0];
                found++;
            }
    }
    fclose(f);
    return found == N_FEATURES + 1;
}
