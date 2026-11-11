/* Unit tests. Build and run with `make test`. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../src/categories.h"
#include "../src/data.h"
#include "../src/markov.h"
#include "../src/model.h"
#include "../src/montecarlo.h"
#include "../src/rng.h"

static int failures = 0;

#define CHECK(cond)                                                    \
    do {                                                               \
        if (!(cond)) {                                                 \
            fprintf(stderr, "  FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
            failures++;                                                \
        }                                                              \
    } while (0)

static Driver random_driver(Rng *r)
{
    Driver d;
    d.age = 18 + (int)rng_below(r, 72);
    d.veh_age = (int)rng_below(r, 20);
    d.veh_power = 4 + (int)rng_below(r, 6);
    d.bonus_malus = 50 + 100 * rng_uniform(r);
    d.density = 10 + 20000 * rng_uniform(r);
    d.diesel = rng_uniform(r) < 0.5;
    d.exposure = 0.2 + 0.8 * rng_uniform(r);
    return d;
}

static void test_poisson_sampler_moments(void)
{
    Rng r;
    rng_seed(&r, 1);
    const double rates[] = {0.1, 2.0, 45.0};
    for (int k = 0; k < 3; k++) {
        int n = 200000;
        double sum = 0, sq = 0;
        for (int i = 0; i < n; i++) {
            int x = rng_poisson(&r, rates[k]);
            sum += x, sq += (double)x * x;
        }
        double mean = sum / n, var = sq / n - mean * mean;
        CHECK(fabs(mean - rates[k]) < 4 * sqrt(rates[k] / n));
        CHECK(fabs(var / rates[k] - 1) < 0.03);
    }
}

static void test_glm_recovers_known_coefficients(void)
{
    Rng r;
    rng_seed(&r, 2);
    size_t n = 200000;
    Driver *d = malloc(n * sizeof *d);
    Model truth = {0}, fit = {0};
    for (int j = 0; j < N_FEATURES; j++)
        truth.coef[j] = 0.2 * (rng_uniform(&r) - 0.5);
    truth.coef[0] = -2.5;
    /* Claims are simulated at the full-year rate, so the exposure controls are truly 0. */
    truth.coef[FEATURE_EXPOSURE] = truth.coef[FEATURE_EXPOSURE + 1] = 0;
    for (size_t i = 0; i < n; i++) {
        d[i] = random_driver(&r);
        d[i].claims = rng_poisson(&r, d[i].exposure * model_lambda(&truth, &d[i]));
    }
    CHECK(model_fit(&fit, d, n, NULL));
    for (int j = 0; j < N_FEATURES; j++)
        CHECK(fabs(fit.coef[j] - truth.coef[j]) < 4 * fit.std_err[j] + 1e-3);
    free(d);
}

static void test_gamma_sampler_moments(void)
{
    Rng r;
    rng_seed(&r, 6);
    const double phis[] = {0.3, 0.8, 2.0};  /* shape 1/phi on both sides of 1 */
    for (int k = 0; k < 3; k++) {
        int n = 200000;
        double sum = 0, sq = 0;
        for (int i = 0; i < n; i++) {
            double x = draw_frailty(phis[k], &r);
            CHECK(x >= 0);
            sum += x, sq += x * x;
        }
        double mean = sum / n, var = sq / n - mean * mean;
        CHECK(fabs(mean - 1) < 4 * sqrt(phis[k] / n));
        CHECK(fabs(var / phis[k] - 1) < 0.05);
    }
    CHECK(draw_frailty(0, &r) == 1);
}

/* With every age slope equal to 1, log lambda grows by exactly 0.1 per year of age up
 * to AGE_MAX and is flat after it: the spline segments join without jumps. */
static void test_age_spline_is_continuous(void)
{
    Model m = {0};
    for (int j = FEATURE_AGE; j < FEATURE_AGE + N_AGE_SEGMENTS; j++)
        m.coef[j] = 1;
    Driver d = {.exposure = 1, .age = MIN_DRIVER_AGE, .veh_age = 5, .veh_power = 4,
                .bonus_malus = 0, .density = 1, .diesel = 0};
    for (; d.age <= AGE_MAX + 5; d.age++) {
        int a = d.age < AGE_MAX ? d.age : AGE_MAX;
        CHECK(fabs(log(model_lambda(&m, &d)) - (a - MIN_DRIVER_AGE) / 10.0) < 1e-12);
    }
}

static void test_lambda_is_for_a_full_year(void)
{
    Model m = {0};
    m.coef[0] = -2;
    m.coef[FEATURE_EXPOSURE] = 1;
    m.coef[FEATURE_EXPOSURE + 1] = 0.5;
    Driver d = {.exposure = 0.1, .age = 40, .veh_age = 5, .veh_power = 4, .bonus_malus = 0,
                .density = 1, .diesel = 0};
    CHECK(fabs(model_lambda(&m, &d) - exp(-2)) < 1e-12);
    CHECK(fabs(model_expected_claims(&m, &d) - 0.1 * exp(-1)) < 1e-12);
}

static void test_data_keeps_only_adult_drivers(void)
{
    const char *path = "test_data.arff";
    FILE *f = fopen(path, "w");
    CHECK(f != NULL);
    if (!f)
        return;
    fprintf(f, "@relation test\n@data\n"
               "1,0,0.5,'A',5,3,17,100,'B1','Regular',100,'R11'\n"
               "2,1,1.0,'A',5,3,18,100,'B1','Diesel',100,'R11'\n"
               "3,0,0.8,'A',5,3,45,50,'B1','Regular',100,'R11'\n");
    fclose(f);
    size_t n = 0;
    Driver *d = data_load(path, &n);
    remove(path);
    CHECK(d != NULL && n == 2);
    for (size_t i = 0; d && i < n; i++)
        CHECK(d[i].age >= MIN_DRIVER_AGE);
    CHECK(d && d[0].diesel == 1);
    free(d);
}

static void test_lifetime_shares_sum_to_one(void)
{
    Rng r;
    rng_seed(&r, 7);
    Model m = {0}, truth = {0};
    m.coef[0] = -2.5;
    m.coef[FEATURE_BONUS_MALUS] = 2;
    truth.coef[0] = -2.3;
    truth.dispersion = 0.8;
    double th[N_LEVELS - 1] = {0.08, 0.1, 0.13}, share[11 * N_LEVELS], claims[11];
    Driver d[200];
    for (int i = 0; i < 200; i++)
        d[i] = random_driver(&r);
    simulate_lifetime(d, 200, &m, &truth, th, 18, 28, &r, share, claims);
    for (int k = 0; k < 11; k++) {
        double sum = 0;
        for (int l = 0; l < N_LEVELS; l++)
            sum += share[k * N_LEVELS + l];
        CHECK(fabs(sum - 1) < 1e-9);
        CHECK(claims[k] >= 0);
    }
    /* Everyone starts as a new driver, with the same bonus-malus and level. */
    CHECK(fabs(share[categorize(model_lambda(&m, &(Driver){.age = 18, .veh_age = 5,
        .veh_power = 4, .bonus_malus = BM_NEW_DRIVER, .density = 1}), th)] - 1) < 1e-9);
}

static void test_categorize_splits_into_quartiles(void)
{
    Rng r;
    rng_seed(&r, 3);
    size_t n = 10000;
    double lam[10000], th[N_LEVELS - 1];
    int counts[N_LEVELS] = {0};
    for (size_t i = 0; i < n; i++)
        lam[i] = -0.1 * log(1 - rng_uniform(&r));
    fit_thresholds(lam, n, th);
    for (size_t i = 0; i < n; i++)
        counts[categorize(lam[i], th)]++;
    for (int l = 0; l < N_LEVELS; l++)
        CHECK(abs(counts[l] - 2500) <= 100);
}

static void test_monte_carlo_score_converges_to_exact_probability(void)
{
    Rng r;
    rng_seed(&r, 4);
    const double rates[] = {0.05, 0.1, 0.3};
    for (int k = 0; k < 3; k++) {
        RiskScore s = risk_score(rates[k], 200000, &r);
        CHECK(fabs(s.prob_claim - (1 - exp(-rates[k]))) < 4 * s.prob_claim_se);
    }
}

static void test_bonus_malus_rule(void)
{
    CHECK(fabs(update_bonus_malus(100, 0, BONUS_FACTOR, MALUS_FACTOR) - 95) < 1e-9);
    CHECK(fabs(update_bonus_malus(100, 1, BONUS_FACTOR, MALUS_FACTOR) - 125) < 1e-9);
    CHECK(fabs(update_bonus_malus(52, 0, BONUS_FACTOR, MALUS_FACTOR) - 50) < 1e-9);
    CHECK(fabs(update_bonus_malus(300, 1, BONUS_FACTOR, MALUS_FACTOR) - 350) < 1e-9);
}

static void test_transition_matrix_is_stochastic(void)
{
    Rng r;
    rng_seed(&r, 5);
    size_t n = 5000;
    int *states = malloc(2 * n * sizeof *states);
    double P[N_LEVELS * N_LEVELS];
    for (size_t i = 0; i < 2 * n; i++)
        states[i] = (int)rng_below(&r, N_LEVELS);
    transition_matrix(states, n, P);
    for (int a = 0; a < N_LEVELS; a++) {
        double row = 0;
        for (int b = 0; b < N_LEVELS; b++) {
            CHECK(P[a * N_LEVELS + b] >= 0);
            row += P[a * N_LEVELS + b];
        }
        CHECK(fabs(row - 1) < 1e-12);
    }
    free(states);
}

static void test_stationary_distribution_is_fixed_point(void)
{
    double P[] = {0.9, 0.1, 0, 0, 0.1, 0.8, 0.1, 0, 0, 0.1, 0.8, 0.1, 0, 0, 0.2, 0.8};
    double pi[N_LEVELS], proj[6 * N_LEVELS], sum = 0;
    CHECK(stationary_distribution(P, pi));
    project(pi, P, 5, proj);
    for (int j = 0; j < N_LEVELS; j++) {
        sum += pi[j];
        CHECK(fabs(proj[1 * N_LEVELS + j] - pi[j]) < 1e-12);  /* pi P = pi */
        CHECK(fabs(proj[5 * N_LEVELS + j] - pi[j]) < 1e-12);  /* pi P^5 = pi */
    }
    CHECK(fabs(sum - 1) < 1e-12);
}

int main(void)
{
    struct { const char *name; void (*fn)(void); } tests[] = {
        {"poisson_sampler_moments", test_poisson_sampler_moments},
        {"gamma_sampler_moments", test_gamma_sampler_moments},
        {"age_spline_is_continuous", test_age_spline_is_continuous},
        {"lambda_is_for_a_full_year", test_lambda_is_for_a_full_year},
        {"data_keeps_only_adult_drivers", test_data_keeps_only_adult_drivers},
        {"lifetime_shares_sum_to_one", test_lifetime_shares_sum_to_one},
        {"glm_recovers_known_coefficients", test_glm_recovers_known_coefficients},
        {"categorize_splits_into_quartiles", test_categorize_splits_into_quartiles},
        {"monte_carlo_score_converges_to_exact_probability", test_monte_carlo_score_converges_to_exact_probability},
        {"bonus_malus_rule", test_bonus_malus_rule},
        {"transition_matrix_is_stochastic", test_transition_matrix_is_stochastic},
        {"stationary_distribution_is_fixed_point", test_stationary_distribution_is_fixed_point},
    };
    int n = sizeof tests / sizeof tests[0];
    for (int i = 0; i < n; i++) {
        int before = failures;
        tests[i].fn();
        printf("%s %s\n", failures == before ? "PASS" : "FAIL", tests[i].name);
    }
    printf("%d tests, %s\n", n, failures ? "FAILURES" : "all passed");
    return failures ? 1 : 0;
}
