#include "data.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static double min_d(double a, double b) { return a < b ? a : b; }
static int min_i(int a, int b) { return a < b ? a : b; }

/* Cap extreme values, following Noll, Salzmann & Wuthrich (2018):
 * more than 4 claims is treated as a data error, exposure is capped at one year,
 * and vehicle age, power and bonus-malus are capped where data becomes sparse. */
static void clean(Driver *d)
{
    d->claims = min_i(d->claims, 4);
    d->exposure = min_d(d->exposure, 1.0);
    d->veh_age = min_i(d->veh_age, 20);
    d->veh_power = min_i(d->veh_power, 9);
    d->bonus_malus = min_d(d->bonus_malus, 150.0);
}

/* One data row of the ARFF file:
 * IDpol,ClaimNb,Exposure,'Area',VehPower,VehAge,DrivAge,BonusMalus,'VehBrand','VehGas',Density,'Region' */
static int parse_row(char *line, Driver *d)
{
    char *field[12];
    int n = 0;
    for (char *tok = strtok(line, ",\r\n"); tok && n < 12; tok = strtok(NULL, ",\r\n"))
        field[n++] = tok;
    if (n != 12)
        return 0;
    d->claims = atoi(field[1]);
    d->exposure = atof(field[2]);
    d->veh_power = atoi(field[4]);
    d->veh_age = atoi(field[5]);
    d->age = atoi(field[6]);
    d->bonus_malus = atof(field[7]);
    d->diesel = strstr(field[9], "Diesel") != NULL;
    d->density = atof(field[10]);
    return 1;
}

Driver *data_load(const char *path, size_t *n)
{
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "Cannot open %s.\nDownload it first:\n  curl -L -o %s %s\n",
                path, path, DATA_URL);
        return NULL;
    }
    size_t cap = 1 << 20, count = 0;
    Driver *drivers = malloc(cap * sizeof *drivers);
    char line[512];
    int in_data = 0;
    while (drivers && fgets(line, sizeof line, f)) {
        if (!in_data) {
            in_data = strncmp(line, "@data", 5) == 0 || strncmp(line, "@DATA", 5) == 0;
            continue;
        }
        Driver d;
        if (!parse_row(line, &d))
            continue;
        clean(&d);
        if (d.exposure <= 0 || d.age < MIN_DRIVER_AGE)
            continue;
        if (count == cap) {
            Driver *grown = realloc(drivers, (cap *= 2) * sizeof *drivers);
            if (!grown) {
                free(drivers);
                drivers = NULL;
                break;
            }
            drivers = grown;
        }
        drivers[count++] = d;
    }
    fclose(f);
    if (!drivers || count == 0) {
        fprintf(stderr, "No driver records read from %s\n", path);
        free(drivers);
        return NULL;
    }
    *n = count;
    return drivers;
}
