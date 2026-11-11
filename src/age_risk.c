#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "age_risk.h"

#define AGE_GROUPS 16
#define MAX_AGE 100
#define MAX_YEARS 1000000

typedef struct {
    long policies;
    double claims;
    double exposure;
} AgeStats;

typedef struct {
    uint64_t state;
} Random;

static const char *const GROUP_NAMES[AGE_GROUPS] = {
    "18-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54",
    "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85-89", "90+"
};
static const double GROUP_CENTERS[AGE_GROUPS] = {
    18.5, 22, 27, 32, 37, 42, 47, 52, 57, 62, 67, 72, 77, 82, 87, 95
};

int age_group_index(int age)
{
    if (age < 18 || age > MAX_AGE)
        return -1;
    if (age < 20)
        return 0;
    if (age >= 90)
        return AGE_GROUPS - 1;
    return 1 + (age - 20) / 5;
}

static uint64_t random_next(Random *random)
{
    uint64_t value = random->state;
    value ^= value >> 12;
    value ^= value << 25;
    value ^= value >> 27;
    random->state = value;
    return value * UINT64_C(2685821657736338717);
}

static double random_uniform(Random *random)
{
    return (double)(random_next(random) >> 11) * 0x1.0p-53;
}

static int poisson_sample(Random *random, double rate)
{
    int claims = 0;
    while (rate > 30) {
        claims += poisson_sample(random, 30);
        rate -= 30;
    }

    double limit = exp(-rate);
    double product = random_uniform(random);
    while (product > limit) {
        claims++;
        product *= random_uniform(random);
    }
    return claims;
}

static void confidence_interval(int successes, int trials, double *lower, double *upper)
{
    const double z = 1.96;
    double probability = (double)successes / trials;
    double denominator = 1 + z * z / trials;
    double center = (probability + z * z / (2 * trials)) / denominator;
    double margin = z * sqrt(probability * (1 - probability) / trials +
                             z * z / (4.0 * trials * trials)) / denominator;
    *lower = fmax(0, center - margin);
    *upper = fmin(1, center + margin);
}

static int parse_number(const char *text, double *number)
{
    char *end;
    errno = 0;
    double value = strtod(text, &end);
    if (errno || end == text || *end != '\0' || !isfinite(value))
        return 0;
    *number = value;
    return 1;
}

static int read_data(const char *path, AgeStats stats[AGE_GROUPS], FILE *errors)
{
    FILE *file = fopen(path, "r");
    if (!file) {
        fprintf(errors, "No se pudo abrir el dataset: %s\n", path);
        return 0;
    }

    char line[1024];
    long valid_policies = 0;
    while (fgets(line, sizeof line, file)) {
        if (line[0] == '@' || line[0] == '\n' || line[0] == '\r')
            continue;

        char *fields[12];
        char *field = line;
        int field_count = 0;
        while (field_count < 12) {
            fields[field_count++] = field;
            char *comma = strchr(field, ',');
            if (!comma)
                break;
            *comma = '\0';
            field = comma + 1;
        }
        if (field_count != 12)
            continue;

        char *line_end = strpbrk(fields[11], "\r\n");
        if (line_end)
            *line_end = '\0';

        double claims, exposure, age_value;
        if (!parse_number(fields[1], &claims) ||
            !parse_number(fields[2], &exposure) ||
            !parse_number(fields[6], &age_value) ||
            claims < 0 || exposure <= 0 ||
            age_value < 18 || age_value > MAX_AGE || floor(age_value) != age_value)
            continue;

        int index = age_group_index((int)age_value);
        stats[index].policies++;
        stats[index].claims += claims;
        stats[index].exposure += exposure;
        valid_policies++;
    }

    int read_ok = !ferror(file);
    fclose(file);
    if (!read_ok) {
        fprintf(errors, "Error al leer el dataset: %s\n", path);
        return 0;
    }
    if (valid_policies == 0) {
        fprintf(errors, "El dataset no contiene pólizas adultas válidas.\n");
        return 0;
    }
    return 1;
}

int run_analysis(const char *path, int selected_age, int years, uint64_t seed,
                 FILE *output, FILE *errors)
{
    if (!path || !output || !errors || age_group_index(selected_age) < 0 ||
        years < 1 || years > MAX_YEARS) {
        if (errors)
            fprintf(errors, "Edad o cantidad de años simulados fuera de rango.\n");
        return 0;
    }

    AgeStats stats[AGE_GROUPS] = {{0}};
    if (!read_data(path, stats, errors))
        return 0;

    int selected_group = age_group_index(selected_age);
    if (stats[selected_group].policies == 0) {
        fprintf(errors, "No hay pólizas válidas para el grupo de edad %s.\n",
                GROUP_NAMES[selected_group]);
        return 0;
    }

    Random random = {seed ? seed : UINT64_C(88172645463325252)};
    fprintf(output, "age_group,age_center,policies,claims,exposure,rate,exact_probability,"
                    "mc_probability,ci_lower,ci_upper,error_plus,error_minus,selected\n");
    for (int i = 0; i < AGE_GROUPS; i++) {
        if (stats[i].policies == 0 || stats[i].exposure <= 0)
            continue;

        double rate = stats[i].claims / stats[i].exposure;
        int years_with_claim = 0;
        for (int year = 0; year < years; year++)
            years_with_claim += poisson_sample(&random, rate) > 0;

        double probability = (double)years_with_claim / years;
        double lower, upper;
        confidence_interval(years_with_claim, years, &lower, &upper);
        fprintf(output, "%s,%.1f,%ld,%.0f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%d\n",
                GROUP_NAMES[i], GROUP_CENTERS[i], stats[i].policies, stats[i].claims,
                stats[i].exposure, rate, 1 - exp(-rate), probability, lower, upper,
                upper - probability, probability - lower, i == selected_group);
    }
    if (ferror(output)) {
        fprintf(errors, "Error al escribir los resultados.\n");
        return 0;
    }
    return 1;
}

#ifndef AGE_RISK_NO_MAIN
static void usage(const char *program)
{
    fprintf(stderr,
            "Uso: %s --data ARCHIVO --age EDAD --years ANIOS [--seed SEMILLA]\n"
            "  EDAD: 18-100; ANIOS: 1-%d\n",
            program, MAX_YEARS);
}

static int parse_integer(const char *text, long min, long max, long *value)
{
    char *end;
    errno = 0;
    long parsed = strtol(text, &end, 10);
    if (errno || end == text || *end != '\0' || parsed < min || parsed > max)
        return 0;
    *value = parsed;
    return 1;
}

int main(int argc, char **argv)
{
    const char *path = NULL;
    long age = -1, years = 10000, seed = (long)time(NULL);
    for (int i = 1; i < argc; i++) {
        if (i + 1 >= argc) {
            usage(argv[0]);
            return 2;
        }
        const char *option = argv[i++];
        const char *value = argv[i];
        if (strcmp(option, "--data") == 0)
            path = value;
        else if (strcmp(option, "--age") == 0) {
            if (!parse_integer(value, 18, MAX_AGE, &age)) {
                fprintf(stderr, "La edad debe estar entre 18 y %d.\n", MAX_AGE);
                return 2;
            }
        } else if (strcmp(option, "--years") == 0) {
            if (!parse_integer(value, 1, MAX_YEARS, &years)) {
                fprintf(stderr, "Los años simulados deben estar entre 1 y %d.\n", MAX_YEARS);
                return 2;
            }
        } else if (strcmp(option, "--seed") == 0) {
            if (!parse_integer(value, 0, 2147483647, &seed)) {
                fprintf(stderr, "La semilla debe ser un entero no negativo.\n");
                return 2;
            }
        } else {
            usage(argv[0]);
            return 2;
        }
    }
    if (!path || age < 0) {
        usage(argv[0]);
        return 2;
    }
    return run_analysis(path, (int)age, (int)years, (uint64_t)seed, stdout, stderr) ? 0 : 1;
}
#endif
