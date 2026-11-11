#define AGE_RISK_NO_MAIN
#include "../src/age_risk.c"

static int failures;

#define CHECK(condition)                                                        \
    do {                                                                         \
        if (!(condition)) {                                                      \
            fprintf(stderr, "FAIL line %d: %s\n", __LINE__, #condition);        \
            failures++;                                                         \
        }                                                                        \
    } while (0)

static void test_age_groups(void)
{
    CHECK(age_group_index(17) == -1);
    CHECK(age_group_index(18) == 0);
    CHECK(age_group_index(19) == 0);
    CHECK(age_group_index(20) == 1);
    CHECK(age_group_index(24) == 1);
    CHECK(age_group_index(25) == 2);
    CHECK(age_group_index(89) == 14);
    CHECK(age_group_index(90) == 15);
    CHECK(age_group_index(100) == 15);
    CHECK(age_group_index(101) == -1);
}

static void test_poisson_simulation(void)
{
    const int trials = 200000;
    const double rate = 0.1;
    Random random = {12345};
    int years_with_claim = 0;
    long total_claims = 0;
    for (int i = 0; i < trials; i++) {
        int claims = poisson_sample(&random, rate);
        total_claims += claims;
        years_with_claim += claims > 0;
    }
    double observed = (double)years_with_claim / trials;
    CHECK(fabs((double)total_claims / trials - rate) < 0.004);
    CHECK(fabs(observed - (1 - exp(-rate))) < 0.005);
}

static void test_dataset_rates(void)
{
    const char *path = "test_age_risk_fixture.arff";
    FILE *fixture = fopen(path, "w");
    CHECK(fixture != NULL);
    if (!fixture)
        return;
    fputs("@relation test\n@data\n"
          "1,1,0.5,'D',5,0,20,50,'B12','Regular',1217,'R82'\n"
          "2,0,0.5,'D',5,0,22,50,'B12','Regular',1217,'R82'\n"
          "3,2,1.0,'D',5,0,25,50,'B12','Regular',1217,'R82'\n"
          "4,4,1.0,'D',5,0,17,50,'B12','Regular',1217,'R82'\n"
          "5,4,0.0,'D',5,0,23,50,'B12','Regular',1217,'R82'\n", fixture);
    fclose(fixture);

    FILE *output = tmpfile();
    FILE *errors = tmpfile();
    CHECK(output != NULL && errors != NULL);
    if (!output || !errors) {
        if (output)
            fclose(output);
        if (errors)
            fclose(errors);
        remove(path);
        return;
    }

    CHECK(run_analysis(path, 22, 100000, 42, output, errors));
    rewind(output);
    char result[4096] = {0};
    size_t size = fread(result, 1, sizeof result - 1, output);
    result[size] = '\0';
    char *selected_row = strstr(result, "20-24,22.0,2,1,1.000000,1.000000,0.632121,");
    CHECK(selected_row != NULL);
    CHECK(selected_row && strstr(selected_row, ",1\n") != NULL);
    CHECK(strstr(result, "25-29,27.0,1,2,1.000000,2.000000,0.864665,") != NULL);

    fclose(output);
    fclose(errors);
    remove(path);
}

int main(void)
{
    test_age_groups();
    test_poisson_simulation();
    test_dataset_rates();
    printf("%s (%d fallas)\n", failures ? "Pruebas fallidas" : "Todas las pruebas pasaron",
           failures);
    return failures ? 1 : 0;
}
