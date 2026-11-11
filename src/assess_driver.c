/* Assess the risk level of a single driver with the fitted model, and project it.
 *
 * The projection follows the driver for --horizon years in many simulated futures:
 * each year they have claims at their current risk, their bonus-malus is updated with
 * the French rule, they get one year older, and the insurer re-categorizes them.
 *
 * Example:
 *   ./assess_driver --age 23 --veh-age 3 --veh-power 6 --bonus-malus 90 --density 1200 --gas Regular
 *
 * Run ./run_analysis first to create results/model.txt.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "categories.h"
#include "markov.h"
#include "model.h"
#include "montecarlo.h"
#include "rng.h"

#define MODEL_PATH "results/model.txt"
#define MAX_HORIZON 50

static void usage(const char *prog)
{
    fprintf(stderr,
            "Usage: %s --age N --veh-age N [options]\n"
            "  --age N            driver age (years)\n"
            "  --veh-age N        vehicle age (years)\n"
            "  --veh-power N      vehicle power class, 4-15 (default 6)\n"
            "  --bonus-malus X    bonus-malus coefficient x100: 50 = best, 100 = new driver (default 100)\n"
            "  --density X        population density where you live, inhabitants/km2 (default 1000)\n"
            "  --gas Regular|Diesel  fuel type (default Regular)\n"
            "  --years N          Monte Carlo simulated years and futures (default 10000)\n"
            "  --horizon N        years to project, 0-%d (default 10)\n"
            "  --seed N           random seed (default: current time)\n",
            prog, MAX_HORIZON);
    exit(2);
}

int main(int argc, char **argv)
{
    Driver d = {.claims = 0, .exposure = 1, .age = -1, .veh_age = -1, .veh_power = 6,
                .bonus_malus = 100, .density = 1000, .diesel = 0};
    int years = 10000, horizon = 10;
    unsigned long seed = (unsigned long)time(NULL);

    for (int i = 1; i < argc; i++) {
        const char *opt = argv[i], *val = i + 1 < argc ? argv[i + 1] : NULL;
        if (!val)
            usage(argv[0]);
        i++;
        if (strcmp(opt, "--age") == 0)
            d.age = atoi(val);
        else if (strcmp(opt, "--veh-age") == 0)
            d.veh_age = atoi(val);
        else if (strcmp(opt, "--veh-power") == 0)
            d.veh_power = atoi(val);
        else if (strcmp(opt, "--bonus-malus") == 0)
            d.bonus_malus = atof(val);
        else if (strcmp(opt, "--density") == 0)
            d.density = atof(val);
        else if (strcmp(opt, "--gas") == 0)
            d.diesel = strcmp(val, "Diesel") == 0;
        else if (strcmp(opt, "--years") == 0)
            years = atoi(val);
        else if (strcmp(opt, "--horizon") == 0)
            horizon = atoi(val);
        else if (strcmp(opt, "--seed") == 0)
            seed = strtoul(val, NULL, 10);
        else
            usage(argv[0]);
    }
    if (d.age < MIN_DRIVER_AGE || d.veh_age < 0 || d.veh_power < 4 || d.density <= 0 ||
        d.bonus_malus < BM_MIN || d.bonus_malus > BM_MAX || years <= 0 || horizon < 0 ||
        horizon > MAX_HORIZON)
        usage(argv[0]);

    Model m;
    double thresholds[N_LEVELS - 1];
    if (!model_load(MODEL_PATH, &m, thresholds)) {
        fprintf(stderr, "Cannot read %s; run ./run_analysis first.\n", MODEL_PATH);
        return 1;
    }

    Rng rng;
    rng_seed(&rng, seed);
    double lambda = model_lambda(&m, &d);
    RiskScore r = risk_score(lambda, years, &rng);

    printf("Expected claims per year (lambda):  %.4f\n", lambda);
    printf("P(at least one claim), Monte Carlo: %.4f +/- %.4f (%d simulated years)\n",
           r.prob_claim, 1.96 * r.prob_claim_se, years);
    printf("Exact P(at least one claim):        %.4f\n", 1 - exp(-lambda));
    printf("Risk level:                         %s\n", LEVELS[categorize(lambda, thresholds)]);
    printf("Thresholds (lambda):                ");
    for (int t = 0; t < N_LEVELS - 1; t++)
        printf("%s < %.4f | ", LEVELS[t], thresholds[t]);
    printf("%s above\n", LEVELS[N_LEVELS - 1]);
    if (horizon == 0)
        return 0;

    /* Many futures of the same driver are a portfolio of identical copies. */
    Driver *copies = malloc((size_t)years * sizeof *copies);
    int *states = malloc((size_t)(horizon + 1) * (size_t)years * sizeof *states);
    if (!copies || !states) {
        fprintf(stderr, "Out of memory\n");
        return 1;
    }
    for (int i = 0; i < years; i++)
        copies[i] = d;
    simulate_portfolio(copies, (size_t)years, &m, thresholds, horizon, BONUS_FACTOR, MALUS_FACTOR,
                       &rng, states);
    printf("\nProjection (%d simulated futures, French bonus-malus rule)\n", years);
    printf("%5s %4s", "year", "age");
    for (int l = 0; l < N_LEVELS; l++)
        printf(" %14s", LEVELS[l]);
    printf("\n");
    for (int t = 1; t <= horizon; t++) {
        if (t != 1 && t != horizon && t % 5 != 0)
            continue;
        double share[N_LEVELS];
        distribution(states, (size_t)years, t, share);
        printf("%5d %4d", t, d.age + t);
        for (int l = 0; l < N_LEVELS; l++)
            printf(" %13.1f%%", 100 * share[l]);
        printf("\n");
    }
    free(copies);
    free(states);
    return 0;
}
