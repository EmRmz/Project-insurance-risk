CC ?= cc
CFLAGS ?= -std=c11 -O2 -Wall -Wextra -pedantic
LDLIBS = -lm
# The project's virtual environment if it exists, otherwise the system Python
PYTHON ?= $(if $(wildcard .venv/bin/python),.venv/bin/python,python3)

LIB = src/rng.c src/data.c src/linalg.c src/model.c src/categories.c src/montecarlo.c src/markov.c
HDR = $(wildcard src/*.h)
DATA = data/raw/freMTPL2freq.arff
DATA_URL = https://openml.org/data/v1/download/20649148/freMTPL2freq.arff

all: run_analysis assess_driver

run_analysis: src/run_analysis.c $(LIB) $(HDR)
	$(CC) $(CFLAGS) -o $@ src/run_analysis.c $(LIB) $(LDLIBS)

assess_driver: src/assess_driver.c $(LIB) $(HDR)
	$(CC) $(CFLAGS) -o $@ src/assess_driver.c $(LIB) $(LDLIBS)

test_risk: tests/test_risk.c $(LIB) $(HDR)
	$(CC) $(CFLAGS) -o $@ tests/test_risk.c $(LIB) $(LDLIBS)

$(DATA):
	mkdir -p data/raw
	curl -L -o $@ $(DATA_URL)

run: run_analysis $(DATA)
	mkdir -p results
	./run_analysis $(DATA)

# Figures are drawn in Python from the CSV tables written by run_analysis.
# plot: risk level by age for a driver who starts at 18, in a window.
plot:
	$(PYTHON) plots/plot_results.py

# The 8 detailed figures of the reports, plus the age figure, saved as PNG in results/figures/.
figures:
	$(PYTHON) plots/report_figures.py
	$(PYTHON) plots/plot_results.py --save

test: test_risk
	./test_risk

clean:
	rm -f run_analysis assess_driver test_risk

.PHONY: all run plot figures test clean
