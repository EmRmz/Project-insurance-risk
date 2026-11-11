#ifndef AGE_RISK_H
#define AGE_RISK_H

#include <stdint.h>
#include <stdio.h>

int age_group_index(int age);
int run_analysis(const char *path, int selected_age, int years, uint64_t seed,
                 FILE *output, FILE *errors);

#endif
