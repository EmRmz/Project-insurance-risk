/* First approach: fit the risk function, categorize drivers and simulate transitions.
 *
 * Usage:  ./run_analysis [data/raw/freMTPL2freq.arff]
 * Outputs go to results/ (CSV tables and model.txt). The figures are drawn from
 * these tables by plots/plot_results.py.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "categories.h"
#include "data.h"
#include "markov.h"
#include "model.h"
#include "montecarlo.h"
#include "rng.h"

#define SEED 2026
#define N_PORTFOLIO 100000  /* drivers in the transition simulation */
#define N_YEARS 10          /* projection horizon */
#define N_MC_YEARS 5000     /* simulated years per driver for the risk score */
#define TEST_SHARE 0.10
#define N_SCENARIOS 3
#define N_LIFETIME 100000   /* new drivers followed from START_AGE to END_AGE */
#define START_AGE 18
#define END_AGE 80
#define N_UPGRADE_REPS 10   /* one-year simulations of every real driver for the upgrade table */

#define OUT "results/"

static const char *const SCENARIO_NAMES[N_SCENARIOS] = {
    "Current (-5% / +25%)", "Stricter malus (+40%)", "Faster bonus (-10%)"};
static const double SCENARIO_BONUS[N_SCENARIOS] = {0.95, 0.95, 0.90};
static const double SCENARIO_MALUS[N_SCENARIOS] = {1.25, 1.40, 1.25};

static FILE *open_out(const char *path)
{
    FILE *f = fopen(path, "w");
    if (!f) {
        fprintf(stderr, "Cannot write %s (does the results/ directory exist?)\n", path);
        exit(1);
    }
    return f;
}

static void *xmalloc(size_t size)
{
    void *p = malloc(size);
    if (!p) {
        fprintf(stderr, "Out of memory\n");
        exit(1);
    }
    return p;
}

/* ---- Data for the figures (plotted by plots/plot_results.py) ------------ */

/* Observed and fitted claim frequency at each driver age (ages above AGE_MAX pooled). */
static void write_age_fit(const Driver *d, size_t n, const Model *m)
{
    enum { AGES = AGE_MAX - MIN_DRIVER_AGE + 1 };
    double claims[AGES] = {0}, expected[AGES] = {0}, exposure[AGES] = {0};
    for (size_t i = 0; i < n; i++) {
        int a = (d[i].age < AGE_MAX ? d[i].age : AGE_MAX) - MIN_DRIVER_AGE;
        claims[a] += d[i].claims;
        expected[a] += model_expected_claims(m, &d[i]);
        exposure[a] += d[i].exposure;
    }
    FILE *f = open_out(OUT "age_fit.csv");
    fprintf(f, "age,exposure,observed_frequency,model_frequency\n");
    for (int a = 0; a < AGES; a++)
        if (exposure[a] > 0)
            fprintf(f, "%d,%.2f,%.5f,%.5f\n", a + MIN_DRIVER_AGE, exposure[a],
                    claims[a] / exposure[a], expected[a] / exposure[a]);
    fclose(f);
}

/* Multiplicative effect of age on lambda relative to a 45-year-old, all else equal,
 * for the full model and for the true-risk model (without bonus-malus). */
static void write_age_effect(const Model *m, const Model *truth)
{
    Driver ref = {.exposure = 1, .age = 45, .veh_age = 5, .veh_power = 6, .bonus_malus = 100,
                  .density = 1000, .diesel = 0};
    double base_m = model_lambda(m, &ref), base_t = model_lambda(truth, &ref);
    FILE *f = open_out(OUT "age_effect.csv");
    fprintf(f, "age,full_model,true_risk_model\n");
    for (ref.age = MIN_DRIVER_AGE; ref.age <= AGE_MAX; ref.age++)
        fprintf(f, "%d,%.5f,%.5f\n", ref.age, model_lambda(m, &ref) / base_m,
                model_lambda(truth, &ref) / base_t);
    fclose(f);
}

/* Histogram of lambda on log-spaced bins, with the count of each risk level per bin. */
static void write_lambda_histogram(const double *lambda, size_t n, const double *thresholds)
{
    enum { BINS = 80 };
    double lo = lambda[0], hi = quantile(lambda, n, 0.999);
    for (size_t i = 0; i < n; i++)
        lo = fmin(lo, lambda[i]);
    long count[BINS][N_LEVELS] = {{0}};
    for (size_t i = 0; i < n; i++) {
        if (lambda[i] > hi)
            continue;
        int b = (int)(BINS * log(lambda[i] / lo) / log(hi / lo));
        count[b < BINS ? b : BINS - 1][categorize(lambda[i], thresholds)]++;
    }
    FILE *f = open_out(OUT "lambda_histogram.csv");
    fprintf(f, "bin_low,bin_high,%s,%s,%s,%s\n", LEVELS[0], LEVELS[1], LEVELS[2], LEVELS[3]);
    for (int b = 0; b < BINS; b++) {
        fprintf(f, "%.6f,%.6f", lo * pow(hi / lo, (double)b / BINS), lo * pow(hi / lo, (b + 1.0) / BINS));
        for (int l = 0; l < N_LEVELS; l++)
            fprintf(f, ",%ld", count[b][l]);
        fprintf(f, "\n");
    }
    fclose(f);
}

/* Running Monte Carlo estimate of P(claim in a year) for one driver per level,
 * recorded at about 200 log-spaced numbers of simulated years. */
static void write_mc_convergence(const double *examples, Rng *rng)
{
    enum { POINTS = 200 };
    static double est[N_MC_YEARS + 1][N_LEVELS];
    for (int l = 0; l < N_LEVELS; l++) {
        long with_claim = 0;
        for (int yr = 1; yr <= N_MC_YEARS; yr++) {
            with_claim += rng_poisson(rng, examples[l]) > 0;
            est[yr][l] = (double)with_claim / yr;
        }
    }
    FILE *f = open_out(OUT "mc_convergence.csv");
    fprintf(f, "years,%s,%s,%s,%s\n", LEVELS[0], LEVELS[1], LEVELS[2], LEVELS[3]);
    int last = 0;
    for (int k = 0; k < POINTS; k++) {
        int yr = (int)round(pow(N_MC_YEARS, (double)k / (POINTS - 1)));
        if (yr <= last)
            continue;
        last = yr;
        fprintf(f, "%d", yr);
        for (int l = 0; l < N_LEVELS; l++)
            fprintf(f, ",%.5f", est[yr][l]);
        fprintf(f, "\n");
    }
    fclose(f);
}

/* Probability that a driver moves up at least one level next year, by current age and
 * level. Every real driver is simulated for one year N_UPGRADE_REPS times, with the same
 * rules as the portfolio (true risk fixed, frailty, current bonus-malus rule).
 * prob_up_without_claims is the share who move up even with a claim-free year, only
 * because they get one year older. */
static void write_upgrade_by_age(const Driver *d, size_t n, const Model *m,
                                 const double *thresholds, Rng *rng)
{
    enum { AGES = AGE_MAX - MIN_DRIVER_AGE + 1 };
    static double total[AGES][N_LEVELS], up[AGES][N_LEVELS], aging[AGES][N_LEVELS];
    int *states = xmalloc(2 * n * sizeof *states);
    for (int r = 0; r < N_UPGRADE_REPS; r++) {
        simulate_portfolio(d, n, m, thresholds, 1, BONUS_FACTOR, MALUS_FACTOR, rng, states);
        for (size_t i = 0; i < n; i++) {
            int a = (d[i].age < AGE_MAX ? d[i].age : AGE_MAX) - MIN_DRIVER_AGE;
            total[a][states[i]] += 1;
            up[a][states[i]] += states[n + i] > states[i];
        }
    }
    free(states);
    for (size_t i = 0; i < n; i++) {
        Driver next = d[i];
        next.age++;
        next.bonus_malus = update_bonus_malus(next.bonus_malus, 0, BONUS_FACTOR, MALUS_FACTOR);
        int a = (d[i].age < AGE_MAX ? d[i].age : AGE_MAX) - MIN_DRIVER_AGE;
        int now = categorize(model_lambda(m, &d[i]), thresholds);
        aging[a][now] += categorize(model_lambda(m, &next), thresholds) > now;
    }
    FILE *f = open_out(OUT "upgrade_by_age.csv");
    fprintf(f, "age,level,drivers,prob_up,prob_up_without_claims\n");
    for (int a = 0; a < AGES; a++)
        for (int l = 0; l < N_LEVELS - 1; l++)  /* High Risk cannot move up */
            if (total[a][l] > 0)
                fprintf(f, "%d,%s,%.0f,%.5f,%.5f\n", a + MIN_DRIVER_AGE, LEVELS[l],
                        total[a][l] / N_UPGRADE_REPS, up[a][l] / total[a][l],
                        aging[a][l] * N_UPGRADE_REPS / total[a][l]);
    fclose(f);
}

/* ---- Pipeline ----------------------------------------------------------- */

int main(int argc, char **argv)
{
    Rng rng;
    rng_seed(&rng, SEED);

    /* 1-2. Risk function fitted on a real driver population */
    size_t n;
    Driver *d = data_load(argc > 1 ? argv[1] : DATA_PATH, &n);
    if (!d)
        return 1;
    double claims = 0, exposure = 0;
    for (size_t i = 0; i < n; i++)
        claims += d[i].claims, exposure += d[i].exposure;
    printf("Loaded %zu policies, observed frequency %.4f claims/year\n", n, claims / exposure);

    unsigned char *train = xmalloc(n), *test = xmalloc(n);
    double train_claims = 0, train_exposure = 0;
    for (size_t i = 0; i < n; i++) {
        test[i] = rng_uniform(&rng) < TEST_SHARE;
        train[i] = !test[i];
        if (train[i])
            train_claims += d[i].claims, train_exposure += d[i].exposure;
    }
    Model m = {0};
    if (!model_fit(&m, d, n, train)) {
        fprintf(stderr, "Model fit failed\n");
        return 1;
    }
    double dev_null = poisson_deviance(&m, train_claims / train_exposure, d, n, test);
    double dev_model = poisson_deviance(&m, -1, d, n, test);

    FILE *f = open_out(OUT "model_coefficients.csv");
    fprintf(f, "feature,coef,std_err,relative_effect\n");
    printf("\n%-16s %8s %8s %8s\n", "feature", "coef", "std_err", "exp(coef)");
    for (int j = 0; j < N_FEATURES; j++) {
        fprintf(f, "%s,%.4f,%.4f,%.4f\n", FEATURE_NAMES[j], m.coef[j], m.std_err[j], exp(m.coef[j]));
        printf("%-16s %8.3f %8.3f %8.3f\n", FEATURE_NAMES[j], m.coef[j], m.std_err[j], exp(m.coef[j]));
    }
    fclose(f);

    /* Risk the features do not explain, used to draw each simulated driver's frailty */
    m.dispersion = model_dispersion(&m, d, n, train);

    /* 3-4. Risk of every driver and its level */
    double *lambda = xmalloc(n * sizeof *lambda), thresholds[N_LEVELS - 1];
    for (size_t i = 0; i < n; i++)
        lambda[i] = model_lambda(&m, &d[i]);
    fit_thresholds(lambda, n, thresholds);
    if (!model_save(OUT "model.txt", &m, thresholds)) {
        fprintf(stderr, "Cannot write " OUT "model.txt\n");
        return 1;
    }
    write_age_fit(d, n, &m);
    write_lambda_histogram(lambda, n, thresholds);

    /* Monte Carlo risk score on one median driver per level */
    double examples[N_LEVELS];
    f = open_out(OUT "risk_levels.csv");
    fprintf(f, "level,drivers,mean_lambda,min_lambda,max_lambda,example_lambda,"
               "mc_prob_claim,mc_prob_claim_se,exact_prob_claim\n");
    printf("\n%-14s %8s %8s %8s %8s | example: %7s %8s %8s\n", "level", "drivers", "mean_lam",
           "min_lam", "max_lam", "lambda", "MC_prob", "exact");
    for (int l = 0; l < N_LEVELS; l++) {
        size_t count = 0;
        double sum = 0, lo = INFINITY, hi = 0;
        double *in_level = xmalloc(n * sizeof *in_level);
        for (size_t i = 0; i < n; i++)
            if (categorize(lambda[i], thresholds) == l) {
                in_level[count++] = lambda[i];
                sum += lambda[i];
                lo = fmin(lo, lambda[i]);
                hi = fmax(hi, lambda[i]);
            }
        examples[l] = quantile(in_level, count, 0.5);
        free(in_level);
        RiskScore r = risk_score(examples[l], N_MC_YEARS, &rng);
        double exact = 1 - exp(-examples[l]);
        fprintf(f, "%s,%zu,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n", LEVELS[l], count, sum / count, lo,
                hi, examples[l], r.prob_claim, r.prob_claim_se, exact);
        printf("%-14s %8zu %8.4f %8.4f %8.4f | %16.4f %8.4f %8.4f\n", LEVELS[l], count, sum / count,
               lo, hi, examples[l], r.prob_claim, exact);
    }
    fclose(f);
    write_mc_convergence(examples, &rng);

    /* Transitions: Monte Carlo portfolio vs Markov chain.
     * The portfolio is a random sample without replacement (partial Fisher-Yates shuffle). */
    size_t *idx = xmalloc(n * sizeof *idx);
    for (size_t i = 0; i < n; i++)
        idx[i] = i;
    Driver *portfolio = xmalloc(N_PORTFOLIO * sizeof *portfolio);
    for (size_t i = 0; i < N_PORTFOLIO; i++) {
        size_t j = i + (size_t)rng_below(&rng, n - i), t = idx[i];
        idx[i] = idx[j], idx[j] = t;
        portfolio[i] = d[idx[i]];
    }
    int *states = xmalloc((size_t)(N_YEARS + 1) * N_PORTFOLIO * sizeof *states);
    double P[N_LEVELS * N_LEVELS], sim[(N_YEARS + 1) * N_LEVELS], proj[(N_YEARS + 1) * N_LEVELS];
    double stationary[N_LEVELS];

    Rng scenario_rng;
    rng_seed(&scenario_rng, SEED);
    simulate_portfolio(portfolio, N_PORTFOLIO, &m, thresholds, N_YEARS, BONUS_FACTOR, MALUS_FACTOR,
                       &scenario_rng, states);
    transition_matrix(states, N_PORTFOLIO, P);
    for (int t = 0; t <= N_YEARS; t++)
        distribution(states, N_PORTFOLIO, t, sim + t * N_LEVELS);
    project(sim, P, N_YEARS, proj);
    if (!stationary_distribution(P, stationary)) {
        fprintf(stderr, "Transition matrix has no unique stationary distribution\n");
        return 1;
    }

    f = open_out(OUT "transition_matrix.csv");
    fprintf(f, "from\\to,%s,%s,%s,%s\n", LEVELS[0], LEVELS[1], LEVELS[2], LEVELS[3]);
    printf("\nTransition matrix P (rows: from, columns: to)\n");
    for (int a = 0; a < N_LEVELS; a++) {
        fprintf(f, "%s", LEVELS[a]);
        printf("%-14s", LEVELS[a]);
        for (int b = 0; b < N_LEVELS; b++) {
            fprintf(f, ",%.4f", P[a * N_LEVELS + b]);
            printf(" %7.3f", P[a * N_LEVELS + b]);
        }
        fprintf(f, "\n");
        printf("\n");
    }
    fclose(f);

    double max_gap = 0;
    f = open_out(OUT "portfolio_projection.csv");
    fprintf(f, "year");
    for (int l = 0; l < N_LEVELS; l++)
        fprintf(f, ",sim_%s", LEVELS[l]);
    for (int l = 0; l < N_LEVELS; l++)
        fprintf(f, ",markov_%s", LEVELS[l]);
    fprintf(f, "\n");
    for (int t = 0; t <= N_YEARS; t++) {
        fprintf(f, "%d", t);
        for (int l = 0; l < N_LEVELS; l++)
            fprintf(f, ",%.4f", sim[t * N_LEVELS + l]);
        for (int l = 0; l < N_LEVELS; l++)
            fprintf(f, ",%.4f", proj[t * N_LEVELS + l]);
        fprintf(f, "\n");
    }
    fclose(f);
    for (int l = 0; l < N_LEVELS; l++)
        max_gap = fmax(max_gap, fabs(sim[N_YEARS * N_LEVELS + l] - proj[N_YEARS * N_LEVELS + l]));

    /* Drivers followed from age 18. Their vehicle, density and fuel are sampled from the
     * real 18-20-year-old policyholders in the data. */
    size_t n_young = 0;
    for (size_t i = 0; i < n; i++)
        if (d[i].age <= 20)
            idx[n_young++] = i;
    Driver *young = xmalloc(N_LIFETIME * sizeof *young);
    for (size_t i = 0; i < N_LIFETIME; i++)
        young[i] = d[idx[rng_below(&rng, n_young)]];
    /* The driver's true risk comes from a model without bonus-malus: the bonus-malus is what
     * the insurer observes, not a cause of claims. Feeding it back into the claim rate makes
     * drivers with early claims spiral to lambda ~ 1 claim/year and never recover. */
    Model truth = {0};
    truth.excluded[FEATURE_BONUS_MALUS] = 1;
    if (!model_fit(&truth, d, n, train)) {
        fprintf(stderr, "Model fit (without bonus-malus) failed\n");
        return 1;
    }
    /* Drivers differ in risk beyond their features; the spread is estimated from the data. */
    truth.dispersion = model_dispersion(&truth, d, n, train);
    if (!model_save(OUT "model_truth.txt", &truth, thresholds)) {
        fprintf(stderr, "Cannot write " OUT "model_truth.txt\n");
        return 1;
    }
    write_age_effect(&m, &truth);
    double life[(END_AGE - START_AGE + 1) * N_LEVELS], life_claims[END_AGE - START_AGE + 1];
    simulate_lifetime(young, N_LIFETIME, &m, &truth, thresholds, START_AGE, END_AGE, &rng, life,
                      life_claims);
    f = open_out(OUT "levels_by_age.csv");
    fprintf(f, "age,%s,%s,%s,%s,mean_claims\n", LEVELS[0], LEVELS[1], LEVELS[2], LEVELS[3]);
    for (int a = START_AGE; a <= END_AGE; a++) {
        fprintf(f, "%d", a);
        for (int l = 0; l < N_LEVELS; l++)
            fprintf(f, ",%.4f", life[(a - START_AGE) * N_LEVELS + l]);
        fprintf(f, ",%.4f\n", life_claims[a - START_AGE]);
    }
    fclose(f);
    printf("\nNew drivers from age %d (%zu real 18-20-year-old profiles, %d simulated drivers, "
           "dispersion %.3f)\n", START_AGE, n_young, N_LIFETIME, truth.dispersion);
    for (int a = START_AGE; a <= END_AGE; a += 5) {
        printf("age %2d:", a);
        for (int l = 0; l < N_LEVELS; l++)
            printf("  %s %5.1f%%", LEVELS[l], 100 * life[(a - START_AGE) * N_LEVELS + l]);
        printf("\n");
    }
    free(young);

    /* Policy scenarios: same portfolio and same random numbers, different rule */
    double final[N_SCENARIOS][N_LEVELS];
    for (int s = 0; s < N_SCENARIOS; s++) {
        rng_seed(&scenario_rng, SEED);
        simulate_portfolio(portfolio, N_PORTFOLIO, &m, thresholds, N_YEARS, SCENARIO_BONUS[s],
                           SCENARIO_MALUS[s], &scenario_rng, states);
        distribution(states, N_PORTFOLIO, N_YEARS, final[s]);
    }
    f = open_out(OUT "policy_scenarios.csv");
    fprintf(f, "level,%s,%s,%s\n", SCENARIO_NAMES[0], SCENARIO_NAMES[1], SCENARIO_NAMES[2]);
    printf("\nShare of portfolio in year %d by policy\n%-14s", N_YEARS, "");
    for (int s = 0; s < N_SCENARIOS; s++)
        printf(" %22s", SCENARIO_NAMES[s]);
    printf("\n");
    for (int l = 0; l < N_LEVELS; l++) {
        fprintf(f, "%s", LEVELS[l]);
        printf("%-14s", LEVELS[l]);
        for (int s = 0; s < N_SCENARIOS; s++) {
            fprintf(f, ",%.4f", final[s][l]);
            printf(" %22.3f", final[s][l]);
        }
        fprintf(f, "\n");
        printf("\n");
    }
    fclose(f);

    /* Own generator, so this table does not shift the random numbers of the other results */
    Rng upgrade_rng;
    rng_seed(&upgrade_rng, SEED + 1);
    write_upgrade_by_age(d, n, &m, thresholds, &upgrade_rng);

    f = open_out(OUT "metrics.txt");
    size_t n_test = 0;
    for (size_t i = 0; i < n; i++)
        n_test += test[i];
    fprintf(f, "policies %zu\ntrain_policies %zu\ntest_policies %zu\n", n, n - n_test, n_test);
    fprintf(f, "observed_frequency %.4f\n", claims / exposure);
    fprintf(f, "test_deviance_null %.4f\ntest_deviance_model %.4f\ndeviance_reduction_pct %.2f\n",
            dev_null, dev_model, 100 * (1 - dev_model / dev_null));
    fprintf(f, "thresholds %.4f %.4f %.4f\n", thresholds[0], thresholds[1], thresholds[2]);
    fprintf(f, "stationary_distribution %.4f %.4f %.4f %.4f\n", stationary[0], stationary[1],
            stationary[2], stationary[3]);
    fprintf(f, "max_gap_sim_vs_markov_year10 %.4f\n", max_gap);
    fprintf(f, "dispersion_full_model %.4f\ndispersion_true_risk_model %.4f\n", m.dispersion,
            truth.dispersion);
    fclose(f);

    printf("\nTest deviance: null %.4f, model %.4f (%.1f%% reduction)\n", dev_null, dev_model,
           100 * (1 - dev_model / dev_null));
    printf("Thresholds: %.4f %.4f %.4f\n", thresholds[0], thresholds[1], thresholds[2]);
    printf("Stationary distribution: %.3f %.3f %.3f %.3f\n", stationary[0], stationary[1],
           stationary[2], stationary[3]);
    printf("Max gap simulation vs Markov in year %d: %.4f\n", N_YEARS, max_gap);
    printf("Dispersion: full model %.3f, true-risk model %.3f\n", m.dispersion, truth.dispersion);

    free(d), free(train), free(test), free(lambda), free(idx), free(portfolio), free(states);
    return 0;
}
