# Finite Sample Properties of FGLS and White's Heteroskedasticity Procedures: A Monte Carlo Study

**Simulation-based evidence on estimation, testing and robust variance estimation**

## Overview

This project studies the finite-sample properties of econometric estimation and inference procedures in the presence of heteroskedasticity using Monte Carlo simulations in R.

The analysis focuses on three complementary procedures:

- Feasible Generalized Least Squares (FGLS)
- White's test for heteroskedasticity
- White's heteroskedasticity-robust variance estimator

The simulations examine how sample size, heteroskedasticity, and the distribution of the disturbances affect coefficient estimation, hypothesis testing, and variance estimation.

---

## Research Objective

The main objective is to evaluate the finite-sample behavior of commonly used econometric procedures under different forms of heteroskedasticity.

The study considers both estimation and inference. In particular, it examines the sampling properties of FGLS when the heteroskedastic covariance structure is known, the ability of White's test to detect heteroskedasticity, and the finite-sample accuracy of White's robust variance estimator.

---

## Monte Carlo Experiments

All simulations are based on repeated samples generated from controlled data-generating processes.

Each experimental configuration is replicated **5,000 times**, using **seed 2026** to ensure reproducibility.

### 1. Feasible Generalized Least Squares (FGLS)

The first experiment evaluates the finite-sample properties of FGLS under a known heteroskedastic covariance structure.

The model is:

$$y_i = \beta_0 + \beta_1 x_i + u_i$$

with:

$$\beta_0=-3, \qquad \beta_1=0.8$$

and:

$$x_i \sim U(1,50)$$

The disturbance variance takes five values:

$$Var(u_i|X) \in \{4,9,16,25,36\}.$$

The sample sizes considered are:

$$n \in \{5,10,30,100,200,500\}.$$

The experiment evaluates:

- Empirical size
- Empirical power
- Monte Carlo means and medians
- Sampling variability of the coefficient estimators
- Equivalence between FGLS and OLS after a Cholesky transformation

---

### 2. White's Test

The second experiment evaluates the finite-sample size and power of White's test for heteroskedasticity.

The regression model is:

$$y_i =\beta_0+\beta_1x_{1i}+\beta_2x_{2i}+\sqrt{\nu_i}\varepsilon_i$$

with:

$$\beta_0=\beta_1=\beta_2=1.$$

Three error designs are considered.

#### Design 0 — Homoskedastic Normal Errors

$$\nu_i=1,\qquad \varepsilon_i\sim N(0,1)$$

This design is used to evaluate the empirical size of White's test.

#### Design 1 — Heteroskedastic Normal Errors

$$\nu_i=\exp(0.25x_{1i}+0.25x_{2i}),\qquad\varepsilon_i\sim N(0,1)$$

This design is used to evaluate empirical power under heteroskedasticity with normally distributed disturbances.

#### Design 2 — Heteroskedastic Non-Normal Errors

$$\nu_i=\exp(0.25x_{1i}+0.25x_{2i}),\qquad\varepsilon_i\sim t_5$$

This design combines heteroskedasticity with non-normal disturbances.

The sample sizes are:

$$n \in \{20,60,100,200,400,600\}.$$

The analysis evaluates empirical rejection frequencies at the 1%, 5%, and 10% significance levels.

---

### 3. White's Robust Variance Estimator

The third experiment studies the finite-sample accuracy of White's heteroskedasticity-robust covariance estimator under Designs 1 and 2.

The standard HC0 estimator is calculated using the OLS residuals.

The estimated variances are compared with the Monte Carlo covariance matrix of the OLS coefficient estimates.

Two approaches are considered:

1. White's covariance estimator using estimated OLS residuals
2. An artificial benchmark using the true simulated disturbances

Relative bias is calculated for the variance estimates of:

$$\hat{\beta}_0,\qquad\hat{\beta}_1,\qquad\hat{\beta}_2.$$

The aggregate measure is:

$$RB_{\mathrm{total}}=|RB_0|+|RB_1|+|RB_2|.$$

This comparison allows the finite-sample effect of replacing the unobserved disturbances with estimated residuals to be examined.

---

## Main Findings

### FGLS

FGLS performs well when the covariance structure of the disturbances is correctly specified. The Monte Carlo estimates remain close to the true parameter values, while sampling variability decreases substantially with sample size.

The empirical size remains close to the nominal significance levels, while statistical power increases rapidly as the sample size grows.

The FGLS estimates are also numerically equivalent to those obtained by applying OLS to the appropriately Cholesky-transformed model.

### White's Test

White's test exhibits reasonable empirical size under the homoskedastic design.

Its power depends strongly on sample size. Under heteroskedasticity with normally distributed disturbances, power increases substantially as the sample size grows.

The presence of non-normal disturbances in Design 2 reduces the finite-sample power of the test relative to Design 1.

This result highlights the importance of considering both sample size and the distribution of the disturbances when interpreting heteroskedasticity tests.

### White's Robust Variance

White's robust variance estimator exhibits noticeable finite-sample bias, particularly in small samples.

The comparison using true disturbances shows that replacing the unobserved errors with estimated residuals can generate substantially larger discrepancies in the smallest samples.

However, the difference between the two approaches is not uniform across all sample sizes and designs. The results therefore illustrate the importance of distinguishing finite-sample behavior from asymptotic properties.

---

## Repository Structure

```
econometrics_analysis/
│
├── README.md
├── econometrics_analysis.R
│
├── report/
│   └── Econometrics_Analysis.pdf
│
├── figures/
│   ├── fgls_test_size.png
│   ├── fgls_power.png
│   ├── fgls_sampling_variability.png
│   ├── fgls_vs_cholesky.png
│   ├── white_test_size.png
│   ├── white_test_power.png
│   ├── total_bias_residuals.png
│   ├── individual_bias_residuals.png
│   ├── total_bias_true_errors.png
│   └── individual_bias_true_errors.png
│
└── outputs/
    ├── fgls_results.csv
    ├── fgls_cholesky_comparison.csv
    ├── white_test_results.csv
    ├── white_variance_residuals.csv
    └── white_variance_true_errors.csv
```
---

## Reproducibility

The complete analysis is implemented in:

```
econometrics_analysis.R
```
The script performs the Monte Carlo simulations, computes the statistical results, generates the figures, and saves the output tables.

The simulations use:

- R
- 5,000 Monte Carlo replications per configuration
- Seed: 2026

To reproduce the analysis, clone or download the repository and run the R script from the project directory.

---
## Software and Packages

The analysis was conducted in **R**.

The main packages used are:

- `ggplot2` — data visualization
- `dplyr` — data manipulation
- `tidyr` — data reshaping
- `sandwich` — heteroskedasticity-robust covariance estimation

Base R functionality is also used for simulation, estimation, and statistical calculations.

---

## Report

The complete methodological discussion, simulation results, tables, and figures are presented in the report:

`Econometrics_Analysis.pdf`

The report is organized into:

1. Introduction
2. Theoretical Framework
3. Simulation Design
4. Finite-Sample Properties of FGLS
5. Finite-Sample Properties of White's Test
6. Finite-Sample Properties of White's Robust Variance
7. Estimated Errors versus True Errors
8. Comparative Discussion
9. Conclusion

---

## Author
Agustina Gimenez

2026
