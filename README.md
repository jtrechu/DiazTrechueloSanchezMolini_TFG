# CSCEN Experiments Repository

This repository contains the implementation code for the experiments in my thesis on th Cost Sensitive Constrained Elastic Net (CSCEN). The code is organized by experiments with supporting utility files.

## Experiment 1: Prostate Cancer Data Analysis

**Main Implementation File:**
- `ProstateCSCEN.R` - Primary implementation for prostate cancer dataset analysis

**Supporting Files:**
- `ECM_star_function.R` - Computes MSE among sensitive groups
- `tau_min_function.R` - Calculates lower bound for τ hyperparameter
- `tau_max_function.R` - Calculates upper bound for τ hyperparameter

**Visualization:**
- `BetaHeatMaps.R` - Generates coefficient heatmaps
- `AsymptoticBetas.R` - Produces asymptotic beta plots

## Experiment 2: Collinearity Analysis

**Implementation Files:**
- `CorrelationCSCEN.R` - Restricted version implementation
- `CorrelationUnconstrained.R` - Unconstrained version implementation

**Supporting Files:**
- `tau_min_function_2.R`, `tau_max_function_2.R` - Hyperparameter bound calculations
- `ECM_star_function.R` - MSE computation among groups
- `CorrelationDataSet.R` - Synthetic data generation

**Note:** Visualization is handled within the implementation files themselves.

## Experiment 3: Lasso Saturation Analysis

**Implementation File:**
- `HighDimensionalitySelection.R` - Examines lasso saturation effects

## Experiment 4: Tabular Results Experiment

**Implementation Files:**
- `TableExperiment.R` - Main analysis implementation
- `TableData.R` - Dataset generation for this experiment
