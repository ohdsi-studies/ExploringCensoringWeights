# Baseline Simulation Studies

Low-dimensional (baseline) survival analysis simulations with informative censoring.

## Overview

This project implements Monte Carlo simulations to evaluate the performance of different survival analysis methods under various censoring scenarios. The simulations compare unadjusted, IPTW (Inverse Probability of Treatment Weighting), IPCW (Inverse Probability of Censoring Weighting), and combined (IPTW×IPCW) estimators.

## Scenarios

The simulation includes three main scenarios:

1. **Sim 1: Non-differential Censoring** (`phi_Z = 0`)
   - Censoring mechanism does not depend on treatment assignment
   
2. **Sim 2: Differential Censoring** (`phi_Z = 1.5`)
   - Censoring mechanism depends on treatment assignment
   
3. **Sim 3: Proportional Hazards Violation**
   - Non-proportional hazards with time-varying treatment effects
   - Early treatment effect (`beta_Z_early`) before `t_PH`
   - Late treatment effect (`beta_Z_late`) after `t_PH`

Each scenario is run with three levels of informative censoring (IC):
- **None** (`ic_scale = 0`)
- **Weak** (`ic_scale = 1.0`)
- **Strong** (`ic_scale = 1.5`)

## Data Structure

The simulation uses 4 covariates:
- **X1 (Age)**: Continuous, Normal(0, 1)
- **X2 (Sex)**: Binary, Bernoulli(0.5)
- **X3 (Cens1)**: Binary, Bernoulli(0.5)
- **X4 (Cens2)**: Continuous, Normal(0, 1)

## Main Files

### Configuration
- **`sim_small_config.R`**: Simulation parameters, sample sizes, coefficients, and scenario settings

### Core Simulation
- **`run_small.R`**: Main simulation runner (sequential)
  - Generates data for each scenario and IC level
  - Computes unadjusted, IPTW, IPCW, and combined estimators
  - Saves results as `.rdata` files in `results_sim_small/`

- **`run_small_parallel.R`**: Parallel simulation runner
  - Same logic as `run_small.R` but uses parallel processing for replications (4–8× speedup)
  - Sources `run_small_optimized.R` for the optimized single-replication and inline RMST SE logic
  - Set `use_parallel = FALSE` in the main call to disable parallel execution

- **`run_small_optimized.R`**: Optimized simulation machinery
  - Single-replication function for parallel execution
  - Inline RMST SE computation (no second pass)
  - Can be used standalone or via `run_small_parallel.R`

### Utility Functions
- **`sim_small_utils.R`**: 
  - Covariate generation
  - True parameter calculation (counterfactual truth)
  - Weibull and piecewise exponential event time generation
  - IPTW computation
  - True DGP construction (`build_true_etas`)

- **`helper_D_censVars.R`**: 
  - Survival curve calculations (unadjusted, IPTW, IPCW, combined)
  - Cox model fitting
  - Data transformation to long format for IPCW
  - IPCW weight calculation

- **`stratified_survival_curves.R`**: 
  - Stratified survival curve estimation by treatment arm
  - Functions for each method (unadjusted, IPTW, IPCW, combined)

### Plotting Scripts
- **`sim_small_plot_HR.R`**: Hazard ratio plots across scenarios and IC levels
- **`sim_small_plot_KM.R`**: Kaplan-Meier survival curves with truth overlay
- **`sim_small_plot_RMST.R`**: Restricted mean survival time (RMST) difference plots

### Analysis
- **`MCSE_calculations.R`**: Monte Carlo standard error calculations for bias estimates

## Usage

### Running Simulations

1. **Set parameters** in `sim_small_config.R`:
   - `n_per_rep`: Sample size per replication (default: 2500)
   - `n_reps`: Number of replications (default: 200)
   - Adjust other parameters as needed

2. **Run the simulation** (choose one):
   ```r
   # Sequential (simplest)
   source("run_small.R")

   # Parallel (faster; uses run_small_optimized.R internally)
   source("run_small_parallel.R")
   ```
   
   Both produce the same result files:
   - `results_sim_small/res_sim1_{none,weak,strong}.rdata`
   - `results_sim_small/res_sim2_{none,weak,strong}.rdata`
   - `results_sim_small/res_sim3_{none,weak,strong}.rdata`

### Generating Plots

After running simulations, generate plots:

```r
source("sim_small_plot_HR.R")      # Hazard ratio plots
source("sim_small_plot_KM.R")      # Kaplan-Meier curves
source("sim_small_plot_RMST.R")    # RMST differences
```

Plots are saved in `results_sim_small/` as PNG and PDF files.

## Output Structure

Each result file (`.rdata`) contains:
- `tt`: Time grid
- `ic_label`: IC level ("none", "weak", "strong")
- `phiZ`: Treatment effect on censoring
- `surv_strat`: List of survival curve matrices (n_reps × n_timepoints)
  - `real_0`, `real_1`: True counterfactual survival curves
  - `unadj_0`, `unadj_1`: Unadjusted estimates
  - `iptw_0`, `iptw_1`: IPTW estimates
  - `ipcw_0`, `ipcw_1`: IPCW estimates
  - `combo_0`, `combo_1`: Combined (IPTW×IPCW) estimates
- `beta`: List of log-hazard ratio estimates
  - `true`: True causal log-HR
  - `unadj`, `iptw`, `ipcw`, `combo`: Estimates from each method
- `rmst_true_diff`: True RMST difference
- `cens_prop`: Censoring proportion per replication

## Dependencies

Required R packages:
- `survival`: Survival analysis
- `glmnet`: Lasso for propensity score estimation
- `dplyr`, `tidyr`: Data manipulation
- `ggplot2`: Plotting
- `purrr`: Functional programming
- `stringr`: String manipulation
- `MASS`: Statistical functions
- `parallel`: Parallel execution (optional; used by `run_small_parallel.R`, part of base R)

## Notes

- The simulation uses a fixed seed (20251124) for reproducibility
- Default sample size is 2500 per replication with 200 replications
- RMST is calculated at `tau_primary = 4.5` (configurable in `sim_small_config.R`)
- Results are saved in `results_sim_small/` (directory created automatically)
- For faster runs, use `run_small_parallel.R`; it requires the `parallel` package (built-in in R)

---

**Last Updated:** Jan 2026
