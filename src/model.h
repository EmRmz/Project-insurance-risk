/* Risk function: Poisson GLM for annual claim frequency.
 *
 * For a driver with features x observed during an exposure period e (years),
 * the number of claims N is modeled as
 *
 *     N ~ Poisson(e * lambda(x)),    log lambda(x) = beta_0 + x' beta
 *
 * so lambda(x) is the expected number of claims per year. This is the standard
 * frequency model in motor insurance pricing; the feature engineering follows
 * the GLM of Noll, Salzmann & Wuthrich (2018), "Case Study: French Motor
 * Third-Party Liability Claims" (SSRN 3164764), simplified to features a driver
 * can easily report about themselves, with two changes:
 *
 *  - Driver age enters as a piecewise-linear spline instead of age bands, so the
 *    risk changes smoothly with age instead of jumping at 31, 41, 51 and 71.
 *  - Two exposure controls (policy in force less than 3 months, or 3 to 9 months)
 *    absorb the higher claim frequency of short policies in freMTPL2freq. They are
 *    used only during the fit: lambda is always predicted for a full year.
 */
#ifndef MODEL_H
#define MODEL_H

#include <stddef.h>

#include "data.h"

#define N_FEATURES 21
#define N_AGE_SEGMENTS 8
#define FEATURE_AGE 1           /* first of the N_AGE_SEGMENTS age terms */
#define FEATURE_BONUS_MALUS 16
#define FEATURE_EXPOSURE 19     /* first of the 2 exposure controls */
#define AGE_MAX 90              /* the age curve is flat above this age */

extern const char *const FEATURE_NAMES[N_FEATURES];

typedef struct {
    double coef[N_FEATURES];
    double std_err[N_FEATURES];
    unsigned char excluded[N_FEATURES];  /* set before model_fit: coefficient fixed at 0 */
    double dispersion;  /* phi in Var(N) = mu + phi mu^2: spread of true risk between drivers */
} Model;

/* Fill x[N_FEATURES] with the design-matrix row of a driver (intercept included). */
void features(const Driver *d, double *x);

/* Fit by IRLS (Newton's method) on the drivers where use[i] != 0 (all if use is NULL).
 * Features marked in m->excluded are left out of the model. */
int model_fit(Model *m, const Driver *d, size_t n, const unsigned char *use);

/* Expected number of claims in one full year. */
double model_lambda(const Model *m, const Driver *d);

/* Fitted mean: expected claims during the driver's own exposure period, exposure
 * controls included. This is what the model compares against observed claims. */
double model_expected_claims(const Model *m, const Driver *d);

/* Moment estimate of the dispersion phi on the drivers where use[i] != 0:
 * phi = sum((y - mu)^2 - y) / sum(mu^2), floored at 0 (0 means plain Poisson). */
double model_dispersion(const Model *m, const Driver *d, size_t n, const unsigned char *use);

/* Mean Poisson deviance of the drivers where use[i] != 0; rate < 0 uses the model,
 * otherwise a constant annual rate (the null model). Lower is better. */
double poisson_deviance(const Model *m, double rate, const Driver *d, size_t n,
                        const unsigned char *use);

int model_save(const char *path, const Model *m, const double *thresholds);
int model_load(const char *path, Model *m, double *thresholds);

#endif
