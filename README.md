# Project-insurance-risk

## Project Description
Insurance companies need to predict the probability that a customer will change risk categories (e.g., moving from a safe-driver bonus tier to higher rates due to accumulated claims). 

This model uses **Monte Carlo simulation** to estimate driver risk and to simulate the annual transitions of policyholders between different risk states, combining historical data, probability, and linear algebra (Markov chains).

## Objective
Establish risk categories based on claim history (e.g., *Preferential, Standard, Moderate Risk, High Risk*).

## Methodology
1. **Risk function:** Study a concrete risk function from the literature (paper) and use it as the basis of the model.
2. **Driver populations:** Find real driver population data to feed the risk function.
3. **Risk calculation:** Compute the risk score for each driver in the population using Monte Carlo simulation (see below).
4. **Categorization:** Classify drivers into 4 risk levels (*Preferential, Standard, Moderate Risk, High Risk*).
5. **Individual assessment:** Given your own data, use the model to determine which risk category you belong to.

## Monte Carlo Model
Monte Carlo simulation is used in two stages of the project:

### 1. Risk estimation
For each driver, the risk function from the paper (e.g., an expected claim rate λ that depends on age, driving experience, vehicle type, etc.) is used to simulate thousands of years of driving:

- In each simulated year, the number of claims is drawn at random from the driver's claim distribution (e.g., Poisson(λ)).
- Repeating this *N* times produces a distribution of claims per driver, not just a single number.
- The driver's risk score is taken from this distribution (e.g., expected claims, or the probability of having one or more claims per year), and the score determines the driver's risk category.

### 2. Transition simulation
Once drivers are categorized, the 4 risk levels become the states of a Markov chain with transition matrix **P**, where `P[i][j]` is the probability of moving from category *i* to category *j* in one year.

- Simulate a portfolio of *N* drivers year by year: for each driver, draw a random number and use the row of **P** for their current category to choose next year's category.
- Project the portfolio over **5 and 10 years**.
- Compare the simulated results with the exact result `π_t = π_0 · P^t` and with the stationary distribution (`π = π · P`) to validate the simulation.
- Evaluate scenarios by changing **P** (e.g., stricter or more lenient insurer policies) and comparing the resulting portfolio distributions.
