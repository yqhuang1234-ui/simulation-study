# simulation-study

**Author:** Yanqi Huang  
**Course:** BIOS 6618 Advanced Biostatistical Methods — University of Colorado (Fall 2025)

## Overview

This repository contains a Monte Carlo simulation study investigating the effects of **heteroscedasticity** on Ordinary Least Squares (OLS) regression inference. By varying the degree of non-constant error variance, the study evaluates how OLS estimation, standard error accuracy, confidence interval coverage, and hypothesis testing behave as the homoscedasticity assumption is violated.

## Research Question

How does heteroscedasticity affect OLS inference on the slope parameter β₁? Specifically:
- Is the OLS estimator biased under heteroscedasticity?
- How accurate are the model-reported standard errors compared to the empirical standard error?
- Do 95% confidence intervals maintain their nominal coverage probability?

## Simulation Design

| Parameter | Value |
|---|---|
| Sample size (n) | 100 |
| Monte Carlo repetitions | 1,000 |
| Heteroscedasticity levels (c) | 0, 0.3, 1, 6 |
| True intercept (β₀) | 14 |
| True slope (β₁) | 0.4 |
| Predictor (X) | Uniform[1, 6] |
| Random seed | 1234 |

**Variance model:** Error variance is modeled as an exponential function of X, centered at μ = 3.5 (the midpoint of the X range):

```
Var(ε | X) = exp(c · (X − μ)) / normalization_factor
```

The normalization ensures the average variance equals 1 across all values of c, enabling fair comparisons across heteroscedasticity levels. When c = 0, errors are homoscedastic (constant variance). As c increases, variance becomes increasingly non-constant.

## Repository Structure

```
simulation-study/
├── code/
│   ├── 01_simulate-data-fit-ols-heteroscedasticity.r    # Data simulation and OLS fitting
│   ├── 02_analyze-results-ols-heteroscedasticity-simulation.r  # Results analysis and plots
│   └── simulation_shell.Rmd                             # Full reproducible R Markdown report
├── data/                    # RDS files with simulated datasets and model fits
├── doc/                     # Simulation study worksheet (Word)
├── results/
│   └── 2025-10-24/          # PNG outputs: summary table and multi-panel plots
├── PDF Reports/             # Final knitted PDF reports
└── README.md
```

## Code Overview

| File | Purpose |
|---|---|
| `01_simulate-data-fit-ols-heteroscedasticity.r` | Defines helper functions for variance generation, data simulation, and OLS fitting; runs the full simulation loop across all c values and repetitions; saves RDS output files |
| `02_analyze-results-ols-heteroscedasticity-simulation.r` | Loads simulation results; computes bias, empirical SE, model SE, relative SE error, and CI coverage (with Monte Carlo SEs); generates publication-quality tables and plots |
| `simulation_shell.Rmd` | R Markdown document that integrates narrative, code, and results; knits to a PDF report using `bookdown::pdf_document2` |

## Requirements

- **R** (>= 3.5.0)
- Required packages (installed automatically if missing):

```r
install.packages(c("dplyr", "ggplot2", "kableExtra", "webshot2", "magick",
                   "knitr", "tidyverse", "tidyr", "ggh4x", "bookdown"))
```

## Usage

### Option A — Run scripts individually

```bash
# Step 1: Simulate data and fit OLS models
Rscript code/01_simulate-data-fit-ols-heteroscedasticity.r

# Step 2: Analyze results and generate plots
Rscript code/02_analyze-results-ols-heteroscedasticity-simulation.r
```

### Option B — Render the full report

Open R or RStudio and run:

```r
rmarkdown::render(
  "code/simulation_shell.Rmd",
  output_format = "bookdown::pdf_document2"
)
```

## Outputs

| Location | Contents |
|---|---|
| `data/` | RDS files with simulated datasets and OLS fit results for each c value |
| `results/2025-10-24/` | PNG: summary table with MC-SE, multi-panel plot (bias, SE, coverage, relative error) |
| `PDF Reports/` | Final knitted PDF reports — [view report](PDF%20Reports/10-24-25_yanqi-huang_simulation-report_only.pdf) |

## Key Findings

- **Unbiasedness:** OLS estimates of β₁ remain unbiased across all heteroscedasticity levels.
- **Standard error underestimation:** As c increases, the model-reported SE increasingly underestimates the true empirical SE (relative error becomes negative).
- **Coverage degradation:** The 95% CI coverage probability drops below the nominal 95% level as heteroscedasticity strengthens, most notably at c = 6.
- **Implication:** Classical OLS inference is unreliable under heteroscedasticity; robust standard error corrections (e.g., HC2, HC3) are needed to restore valid inference.
