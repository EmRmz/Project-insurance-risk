/* Loading and cleaning of the freMTPL2freq driver population.
 *
 * freMTPL2freq contains 678,013 French motor third-party liability policies
 * (CASdatasets, distributed by OpenML under CC0). Each row is a policy with the
 * number of claims observed during its exposure period (in years).
 */
#ifndef DATA_H
#define DATA_H

#include <stddef.h>

#define DATA_URL "https://openml.org/data/v1/download/20649148/freMTPL2freq.arff"
#define DATA_PATH "data/raw/freMTPL2freq.arff"
#define MIN_DRIVER_AGE 18  /* only adult drivers are kept */

typedef struct {
    int claims;          /* number of claims in the exposure period */
    double exposure;     /* years the policy was in force */
    int age;             /* driver age */
    int veh_age;         /* vehicle age */
    int veh_power;       /* vehicle power class */
    double bonus_malus;  /* claim-history coefficient x100: 50 best, 100 new driver */
    double density;      /* inhabitants per km2 where the driver lives */
    int diesel;          /* 1 = diesel, 0 = regular */
} Driver;

/* Load and clean the dataset, keeping only drivers aged MIN_DRIVER_AGE or more.
 * Returns a malloc'd array (n in *n), or NULL on error. */
Driver *data_load(const char *path, size_t *n);

#endif
