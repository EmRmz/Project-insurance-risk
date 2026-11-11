/* Monte Carlo simulation of claims and of portfolio transitions. */
#ifndef MONTECARLO_H
#define MONTECARLO_H

#include <stddef.h>

#include "data.h"
#include "model.h"
#include "rng.h"

/* French bonus-malus rule (coefficient x100): -5% per claim-free year, +25% per claim. */
#define BONUS_FACTOR 0.95
#define MALUS_FACTOR 1.25
#define BM_MIN 50.0
#define BM_MAX 350.0
#define BM_NEW_DRIVER 100.0

typedef struct {
    double mean_claims;   /* average simulated claims per year */
    double prob_claim;    /* share of simulated years with at least one claim */
    double prob_claim_se; /* standard error of prob_claim, shrinks as 1/sqrt(n_years) */
} RiskScore;

/* Simulate n_years of driving for a driver with annual claim rate lambda. */
RiskScore risk_score(double lambda, int n_years, Rng *rng);

double update_bonus_malus(double bm, int claims, double bonus, double malus);

/* Multiplier of a driver's true claim rate: Gamma with mean 1 and variance phi
 * (always 1 when phi = 0). It stands for what the features do not see, such as how
 * carefully the person drives. */
double draw_frailty(double phi, Rng *rng);

/* Simulate a portfolio year by year and record each driver's risk level.
 *
 * Each driver's true claim rate is the full model's lambda with their bonus-malus of
 * year 0 (the best estimate of their risk given their real claim history), times a
 * frailty drawn with m->dispersion for the risk the model does not see. Claims do not
 * change it: they only change what the insurer sees, not how the person drives.
 * Each year the driver has Poisson(true rate) claims, their bonus-malus is updated with
 * the insurer's rule, they get one year older, and the insurer re-categorizes them with
 * the full model and the current bonus-malus. Vehicle, density and fuel are held fixed.
 *
 * states must hold (n_years + 1) * n ints; states[t * n + i] is driver i's level in year t. */
void simulate_portfolio(const Driver *drivers, size_t n, const Model *m, const double *thresholds,
                        int n_years, double bonus, double malus, Rng *rng, int *states);

/* Follow new drivers from start_age to end_age.
 *
 * Every driver starts at start_age with the new-driver bonus-malus (100). Their true
 * claim rate is the model `truth`, which has no bonus-malus, times a frailty drawn with
 * truth->dispersion. Each year they have Poisson(true rate) claims, their bonus-malus is
 * updated, they get one year older, and the insurer re-categorizes them with the full
 * model `m`. Vehicle, density and fuel are held fixed.
 *
 * share must hold (end_age - start_age + 1) * N_LEVELS doubles;
 * share[(age - start_age) * N_LEVELS + l] is the share of drivers in level l at that age.
 * claims, if not NULL, holds end_age - start_age + 1 doubles: mean claims at each age. */
void simulate_lifetime(const Driver *drivers, size_t n, const Model *m, const Model *truth,
                       const double *thresholds, int start_age, int end_age, Rng *rng,
                       double *share, double *claims);

#endif
