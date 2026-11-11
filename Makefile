CC ?= cc
CFLAGS ?= -std=c11 -O2 -Wall -Wextra -pedantic
LDLIBS = -lm

age_risk: src/age_risk.c src/age_risk.h
	$(CC) $(CFLAGS) -o $@ src/age_risk.c $(LDLIBS)

test_age_risk: tests/test_age_risk.c src/age_risk.c src/age_risk.h
	$(CC) $(CFLAGS) -o $@ tests/test_age_risk.c $(LDLIBS)

test: test_age_risk
	./test_age_risk

run: age_risk
	python -m streamlit run app.py

clean:
	rm -f age_risk age_risk.exe test_age_risk test_age_risk.exe

.PHONY: test run clean
