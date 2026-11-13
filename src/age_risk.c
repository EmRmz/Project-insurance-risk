#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "age_risk.h"

#define AGE_GROUPS 16     /* cantidad de grupos de edad */
#define MAX_AGE 100       /* edad máxima aceptada */
#define MAX_YEARS 1000000 /* máximo de años que se pueden simular */
#define Z_95 1.96         /* valor z para un intervalo de confianza del 95% */

/* Totales de un grupo de edad, acumulados al leer el dataset. */
typedef struct {
    long policies;   /* número de pólizas del grupo */
    double claims;   /* suma de siniestros */
    double exposure; /* suma de años asegurados */
} AgeStats;

static const char *const GROUP_NAMES[AGE_GROUPS] = {
    "18-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54",
    "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85-89", "90+"
};

/* Edad que representa a cada grupo en el eje x de la gráfica. */
static const double GROUP_CENTERS[AGE_GROUPS] = {
    18.5, 22, 27, 32, 37, 42, 47, 52, 57, 62, 67, 72, 77, 82, 87, 95
};

/* Devuelve el número de grupo (0 a 15) de una edad, o -1 si está fuera de rango. */
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

/* Número pseudoaleatorio uniforme entre 0 y 1. */
static double random_uniform(void)
{
    return rand() / ((double)RAND_MAX + 1.0);
}

/*
 * Simula cuántos siniestros tiene un conductor en un año.
 * Método de Knuth: se multiplican números aleatorios hasta que el producto
 * baja de e^(-tasa). La cantidad de multiplicaciones extra sigue una Poisson.
 */
static int poisson_sample(double rate)
{
    double limit = exp(-rate);
    double product = random_uniform();
    int claims = 0;
    while (product > limit) {
        claims++;
        product *= random_uniform();
    }
    return claims;
}

/* Margen de error del 95%: z * raíz(p(1 - p) / N). */
static double margin_of_error(double probability, int trials)
{
    return Z_95 * sqrt(probability * (1 - probability) / trials);
}

/*
 * Lee el archivo .arff y suma pólizas, siniestros y exposición por grupo de edad.
 * Columnas que se usan: 1 = ClaimNb (siniestros), 2 = Exposure, 6 = DrivAge (edad).
 */
static int read_data(const char *path, AgeStats *stats, FILE *errors)
{
    FILE *file = fopen(path, "r");
    if (file == NULL) {
        fprintf(errors, "No se pudo abrir el dataset: %s\n", path);
        return 0;
    }

    char line[1024];
    long valid_policies = 0;
    while (fgets(line, sizeof line, file) != NULL) {
        /* Saltar encabezados (@...), comentarios (%...) y líneas vacías. */
        if (line[0] == '@' || line[0] == '%' || line[0] == '\n' || line[0] == '\r')
            continue;

        double claims = -1, exposure = 0;
        int age = 0;
        char *field = strtok(line, ",");
        for (int column = 0; field != NULL && column <= 6; column++) {
            if (column == 1)
                claims = atof(field);
            else if (column == 2)
                exposure = atof(field);
            else if (column == 6)
                age = atoi(field);
            field = strtok(NULL, ",");
        }

        int group = age_group_index(age);
        if (group < 0 || claims < 0 || exposure <= 0)
            continue; /* registro incompleto o fuera de rango */

        stats[group].policies++;
        stats[group].claims += claims;
        stats[group].exposure += exposure;
        valid_policies++;
    }
    fclose(file);

    if (valid_policies == 0) {
        fprintf(errors, "El dataset no contiene pólizas adultas válidas.\n");
        return 0;
    }
    return 1;
}

/*
 * Calcula la tasa de cada grupo, ejecuta Monte Carlo y escribe una tabla CSV.
 * years = cuántos años se simulan por grupo (N de Monte Carlo).
 */
int run_analysis(const char *path, int selected_age, int years, unsigned int seed,
                 FILE *output, FILE *errors)
{
    int selected_group = age_group_index(selected_age);
    if (selected_group < 0 || years < 1 || years > MAX_YEARS) {
        fprintf(errors, "Edad o cantidad de años simulados fuera de rango.\n");
        return 0;
    }

    AgeStats stats[AGE_GROUPS] = {{0}};
    if (!read_data(path, stats, errors))
        return 0;

    if (stats[selected_group].policies == 0) {
        fprintf(errors, "No hay pólizas válidas para el grupo de edad %s.\n",
                GROUP_NAMES[selected_group]);
        return 0;
    }

    srand(seed);
    fprintf(output, "age_group,age_center,policies,claims,exposure,rate,"
                    "exact_probability,mc_probability,margin,selected\n");

    for (int i = 0; i < AGE_GROUPS; i++) {
        if (stats[i].policies == 0)
            continue;

        /* Modelo: tasa = siniestros / exposición, probabilidad exacta = 1 - e^(-tasa). */
        double rate = stats[i].claims / stats[i].exposure;
        double exact_probability = 1 - exp(-rate);

        /* Monte Carlo: contar en cuántos años simulados hubo al menos un siniestro. */
        int years_with_claim = 0;
        for (int year = 0; year < years; year++) {
            if (poisson_sample(rate) > 0)
                years_with_claim++;
        }
        double mc_probability = (double)years_with_claim / years;
        double margin = margin_of_error(mc_probability, years);

        fprintf(output, "%s,%.1f,%ld,%.0f,%.6f,%.6f,%.6f,%.6f,%.6f,%d\n",
                GROUP_NAMES[i], GROUP_CENTERS[i], stats[i].policies, stats[i].claims,
                stats[i].exposure, rate, exact_probability, mc_probability, margin,
                i == selected_group);
    }
    return 1;
}

#ifndef AGE_RISK_NO_MAIN
static void usage(const char *program)
{
    fprintf(stderr,
            "Uso: %s --data ARCHIVO --age EDAD --years ANIOS [--seed SEMILLA]\n"
            "  EDAD: 18-%d; ANIOS: 1-%d\n",
            program, MAX_AGE, MAX_YEARS);
}

int main(int argc, char **argv)
{
    const char *path = NULL;
    int age = -1;
    int years = 10000;
    unsigned int seed = (unsigned int)time(NULL);

    /* Leer los argumentos en pares: --opcion valor */
    for (int i = 1; i + 1 < argc; i += 2) {
        if (strcmp(argv[i], "--data") == 0)
            path = argv[i + 1];
        else if (strcmp(argv[i], "--age") == 0)
            age = atoi(argv[i + 1]);
        else if (strcmp(argv[i], "--years") == 0)
            years = atoi(argv[i + 1]);
        else if (strcmp(argv[i], "--seed") == 0)
            seed = (unsigned int)atoi(argv[i + 1]);
        else {
            usage(argv[0]);
            return 2;
        }
    }
    if (path == NULL || argc % 2 == 0) {
        usage(argv[0]);
        return 2;
    }
    return run_analysis(path, age, years, seed, stdout, stderr) ? 0 : 1;
}
#endif
