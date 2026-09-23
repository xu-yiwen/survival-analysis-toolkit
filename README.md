# Survival Analysis Toolkit

An interactive R Shiny application providing an end-to-end survival analysis workflow for user-uploaded datasets. Built as the final project for **HBDS 5008 Biostatistics II, Spring 2026** (Weill Cornell Medicine, Instructor: Yushu Shi, PhD).

## Features

The app covers all required analyses from the project specification, plus several extensions:

- Kaplan-Meier estimation, with optional stratification and risk tables
- Log-rank test between groups
- Multivariable Cox proportional hazards regression with a protective check (>= 10 events per covariate)
- Forest plot of hazard ratios
- Proportional hazards assumption diagnostics (cox.zph)
- Optimal cutpoint analysis for continuous variables (surv_cutpoint)
- Covariate-adjusted survival curves (ggadjustedcurves)
- Competing risks analysis with Gray's test (cmprsk)
- Multi-state model state occupation probabilities
- Automated plain-language interpretation under every output
- PNG/CSV download for all figures and tables

## Getting Started

The app requires R with the following packages:

```r
install.packages(c("shiny", "readr", "survival", "survminer", "broom", "dplyr", "ggplot2", "cmprsk"))
```

Run locally from the project directory:

```r
shiny::runApp(".")
```

Then open the app in a browser, upload a CSV file, and map the columns:

- **Time variable**: numeric follow-up time
- **Event indicator**: 0 = censored, 1 = event (the app also auto-converts common encodings such as 1/2 or yes/no/dead)
- **Covariates**: any numeric or categorical predictors

## Data

Two sample CSV files, `veteran.csv` and `survival_test_data.csv`, are included in the repository so the app can be tried immediately. The app works with any CSV in the format described above.
