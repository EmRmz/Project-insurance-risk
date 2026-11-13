#ifndef AGE_RISK_H
#define AGE_RISK_H

#include <stdio.h>

int age_group_index(int age);
int run_analysis(const char *path, int selected_age, int years, unsigned int seed,
                 FILE *output, FILE *errors);

#endif
