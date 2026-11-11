#include "montecarlo.h"

#include <math.h>

#include "categories.h"

RiskScore risk_score(double lambda, int n_years, Rng *rng)
{
    long total = 0, years_with_claim = 0;
    for (int y = 0; y < n_years; y++) {
        int claims = rng_poisson(rng, lambda);
        total += claims;
        years_with_claim += claims > 0;
    }
    RiskScore r;
    r.mean_claims = (double)total / n_years;
    r.prob_claim = (double)years_with_claim / n_years;
    r.prob_claim_se = sqrt(r.prob_claim * (1 - r.prob_claim) / n_years);
    return r;
}

double update_bonus_malus(double bm, int claims, double bonus, double malus)
{
    bm = claims == 0 ? bm * bonus : bm * pow(malus, claims);
    return bm < BM_MIN ? BM_MIN : bm > BM_MAX ? BM_MAX : bm;
}

double draw_frailty(double phi, Rng *rng)
{
    return phi > 0 ? phi * rng_gamma(rng, 1 / phi) : 1;
}

void simulate_portfolio(const Driver *drivers, size_t n, const Model *m, const double *thresholds,
                        int n_years, double bonus, double malus, Rng *rng, int *states)
{
    for (size_t i = 0; i < n; i++) {
        Driver seen = drivers[i], risk = drivers[i];
        double frailty = draw_frailty(m->dispersion, rng);
        for (int t = 0; t <= n_years; t++, seen.age++, risk.age++) {
            states[(size_t)t * n + i] = categorize(model_lambda(m, &seen), thresholds);
            if (t < n_years)
                seen.bonus_malus = update_bonus_malus(
                    seen.bonus_malus, rng_poisson(rng, frailty * model_lambda(m, &risk)), bonus,
                    malus);
        }
    }
}

void simulate_lifetime(const Driver *drivers, size_t n, const Model *m, const Model *truth,
                       const double *thresholds, int start_age, int end_age, Rng *rng,
                       double *share, double *claims)
{
    int n_ages = end_age - start_age + 1;
    for (int k = 0; k < n_ages * N_LEVELS; k++)
        share[k] = 0;
    for (int k = 0; claims && k < n_ages; k++)
        claims[k] = 0;
    for (size_t i = 0; i < n; i++) {
        Driver d = drivers[i];
        double frailty = draw_frailty(truth->dispersion, rng);
        d.age = start_age;
        d.bonus_malus = BM_NEW_DRIVER;
        for (int k = 0; k < n_ages; k++, d.age++) {
            share[k * N_LEVELS + categorize(model_lambda(m, &d), thresholds)] += 1.0 / (double)n;
            int c = rng_poisson(rng, frailty * model_lambda(truth, &d));
            if (claims)
                claims[k] += c / (double)n;
            d.bonus_malus = update_bonus_malus(d.bonus_malus, c, BONUS_FACTOR, MALUS_FACTOR);
        }
    }
}
