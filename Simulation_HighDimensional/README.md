# Simulation 4: High-Dimensional Informative Censoring Study

## 📘 Overview

This project implements Monte Carlo simulations to evaluate the performance of different survival analysis methods under various censoring scenarios in **high dimensional** settings. The simulations compare unadjusted, IPTW (Inverse Probability of Treatment Weighting), IPCW (Inverse Probability of Censoring Weighting), and combined (IPTW×IPCW) estimators.


---

## 📦 Required R Packages

The simulation requires the following R packages:
- `survival` - Survival analysis and Cox models
- `glmnet` - Penalized regression for high-dimensional propensity score and IPCW models
- `MASS` - Multivariate normal generation
- `Matrix` - Sparse matrix operations
- `dplyr`, `tidyr`, `purrr` - Data manipulation
- `ggplot2` - Plotting
- `stringr` - String manipulation

Install missing packages with:
```r
install.packages(c("survival", "glmnet", "MASS", "Matrix", "dplyr", "tidyr", "purrr", "ggplot2", "stringr"))
```

---

## 🎯 Objectives
- Incorporate **high-dimensional covariates** (continuous + binary; correlated AR(1) structure).  
- Introduce **informative and differential censoring** of varying strengths.  
- Evaluate **confounding and unmeasured confounding** mechanisms.  
- Implement **penalized models** (via *glmnet*) for IPTW and IPCW estimation.  
- Compare performance of:
  1. Unadjusted (naïve) estimates  
  2. IPTW  
  3. IPCW  
  4. Combined IPTW + IPCW  

Outcomes are summarized by:
- Marginal HR bias and Monte Carlo error  
- Restricted Mean Survival Time (RMST) differences  
- Kaplan–Meier (KM) curves by treatment arm  
- Effective sample size (ESS) diagnostics

---

## 📁 File Structure

| File | Description |
|------|--------------|
| **`sim4_config.R`** | Defines all global simulation parameters (sample size, replications, number and type of covariates, target censoring, confounding indices, etc.). Includes tunable knobs for binary prevalence, confounding strength, and censoring strength. |
| **`sim4_utils_hd.R`** | Utility functions for high-dimensional data generation, treatment calibration, censor intercept solving, and Weibull event-time sampling. |
| **`helper_D_censVars.R`** | Core IPCW machinery for transforming data into long format, computing stabilized weights, truncation, and weight diagnostics. |
| **`stratified_survival_curves.R`** | Functions for estimating stratified survival curves (unadjusted, IPTW, IPCW, and combined), and for computing RMST differences by treatment arm. |
| **`sim4_run.R`** | Main simulation driver. Iterates over scenarios, runs Monte Carlo replications, computes survival estimates, log-HRs, and diagnostics, and saves `.rdata` results. |
| **`sim4_scenarios.R`** | Sources `sim4_run.R` and runs all scenarios (4A–4E) across all informative-censoring levels. Can override sample sizes and replications. |
| **`sim4_plot_HR.R`** | Produces a single HR comparison plot across scenarios and informative-censoring levels. |
| **`sim4_plot_KM.R`** | Generates Kaplan–Meier survival curves stratified by treatment arm for each scenario, comparing truth vs. estimated curves across methods. |
| **`sim4_plot_RMST.R`** | Produces a single RMST difference plot with confidence intervals and red reference lines for the true RMST difference. |
| **`MCSE_calculations_sim4.R`** | Computes Monte Carlo bias and Monte Carlo standard error (MCSE) for HR and RMST differences. Produces formatted tables with bias (MCSE), coverage, and mean SE for reporting. |
| **`bootstrap_se_calculations.R`** | Optional bootstrap SE for HR and RMST (re-estimates weights in each bootstrap sample). Used when `use_bootstrap_se = TRUE` in config. |

---

## 🧪 Simulation Scenarios

| Scenario | Description | Key Features |
|-----------|--------------|---------------|
| **4A – Baseline** | Standard high-dimensional confounding; no informative censoring. | Reference case |
| **4B – Poor Overlap** | Stronger correlation among treatment covariates; overlap violation. | Tests IPTW stability |
| **4C – Differential Informative Censoring** | Treatment affects both event and censoring processes (φZ set via `phi_Z_levels`), creating differential informative censoring. | Multiple φZ values; IC strength varied |
| **4D – Misspecified Models** | Nonlinear terms and interactions added to true DGP but omitted in analysis models. | Robustness check |
| **4E – Unmeasured Confounding** | Some true confounders omitted from analysis (`idx_unmeasured`). | Tests estimator bias under hidden confounding |

Each scenario is further run under **three informative-censoring strengths** (configurable in `sim4_config.R`):  
- `none` (0.0)  
- `weak` (0.5)  
- `strong` (1.0)  

Scenario **4C** additionally varies **differential censoring** via `phi_Z_levels` in `sim4_config.R` (e.g. `c(-0.5, -0.8)`); filenames use suffixes like `phiZneg0_5`, `phiZneg0_8`.

---

## ⚙️ Key Parameters

| Parameter | Description |
|------------|-------------|
| `n_per_rep`, `n_reps` | Sample size per replication and number of replications. Default in `sim4_config.R` is 5000 and 50; `sim4_scenarios.R` overrides to 2500 and 200. |
| `p_total = 100` | Total covariates (3 continuous, 97 binary). |
| `p_continuous = 3`, `p_binary = 97` | Number of continuous and binary covariates. |
| `ar1_rho = 0.4` | AR(1) correlation among continuous covariates. |
| `ar1_rho_bin = 0.1` | AR(1) correlation among binary covariates (via probit thresholding). |
| `bin_prevalence = 0.30` | Desired prevalence for binary covariates. |
| `target_censoring = 0.45` | Target proportion censored. |
| `target_treated_prop = 0.45` | Target proportion treated. |
| `conf_strength = 1.0` | Controls degree of treatment–outcome confounding. |
| `ic_strength_levels` | Informative censoring dial (e.g. `none=0.0`, `weak=0.5`, `strong=1.0`; scales shared event–censor covariates). |
| `phi_Z_levels` | Differential censoring coefficients for Scenario 4C (e.g. `c(-0.5, -0.8)`); set in `sim4_config.R`. |
| `trunc_lo = 0.01`, `trunc_hi = 0.99` | Weight truncation percentiles for stabilized weights. |
| `tau_primary = 10` | RMST horizon (time point for restricted mean survival time calculation). |
| `use_bootstrap_se` | If `TRUE`, use bootstrap SE for HR/RMST instead of analytical SE; requires `bootstrap_se_calculations.R`. |

---

## 📊 Outputs
Each scenario produces `.rdata` files under `results_sim4/`. Each file name includes an informative-censoring suffix: `_ic-none`, `_ic-weak`, or `_ic-strong`.

| Output | Example Filename | Description |
|---------|------------------|-------------|
| Baseline | `Sim4A_baseline_ic-none.rdata` | Reference simulation without informative censoring. |
| Poor overlap | `Sim4B_poorOverlap_ic-weak.rdata` | Low overlap scenario. |
| Differential censoring | `Sim4C_diffCens_phiZneg0_5_ic-weak.rdata` | 4C with a given φZ (e.g. −0.5); suffix format `phiZ` or `phiZneg` for negative values. |
| Misspecified | `Sim4D_misspecified_ic-weak.rdata` | Model misspecification scenario. |
| Unmeasured confounding | `Sim4E_unmeasured_ic-strong.rdata` | Hidden confounders scenario. |

Checkpoint files (e.g. `*_checkpoint.rdata`) may appear during runs and are removed after successful completion.

### Generated Figures
- **`Sim4_combined_HR_plot.png/pdf`** — Hazard ratio comparison across methods, scenarios, and IC levels (with truth reference lines)  
- **`Sim4_KM_*.png/pdf`** — Kaplan–Meier survival curves for each scenario, showing truth (dashed) vs. estimated (solid) curves by treatment arm, method, and IC level  
- **`Sim4_RMSTdiff.png/pdf`** — RMST difference by method and IC level with confidence intervals and truth reference lines  
- **Diagnostics in result objects:** Effective sample size (ESS) and censoring proportions per replication  

---

## ▶️ Running the Simulation

1. **Run all scenarios (4A–4E):**  
   The easiest way is to use the convenience script:
   ```r
   source("sim4_scenarios.R")
   ```
   This runs all scenarios across all informative-censoring levels (and all `phi_Z_levels` for 4C), saves `.rdata` results to `results_sim4/`, then sources the plot scripts and `MCSE_calculations_sim4.R` so that figures and MCSE tables are generated automatically.
   
   Runs support **checkpointing**: if a scenario is interrupted, re-running will resume from the last saved checkpoint (see `checkpoint_file` and `checkpoint_freq` in `run_scenario()`).
   
   For a single scenario, source `sim4_run.R` and call `run_scenario()` directly (e.g. `run_scenario("4A_baseline", ic_key = "weak")`).

2. **Generate plots (if not run via `sim4_scenarios.R`):**
   ```r
   source("sim4_plot_HR.R")
   source("sim4_plot_KM.R")
   source("sim4_plot_RMST.R")
   source("MCSE_calculations_sim4.R")   # MCSE tables
   ```
   
   Figures and tables are saved to `results_sim4/`.

---

## 🧩 Extending / Customizing

You can easily adjust settings in `sim4_config.R` to run new experiments:
- Modify sample size, number of replications, or covariate count.
- Change binary prevalence (`bin_prevalence`), AR(1) correlation (`ar1_rho`, `ar1_rho_bin`), or censoring targets.
- Tune `conf_strength` or `ic_strength_levels` for sensitivity analysis.
- Adjust `phi_Z_levels` for differential censoring scenarios.
- Override sample sizes in `sim4_scenarios.R` (e.g., set `n_per_rep <- 2000` and `n_reps <- 100` for final runs).
- Add new scenarios by:
  - Adding a new scenario key to `sim4_files` in `sim4_config.R`
  - Extending `build_true_etas()` or adding scenario-specific tweaks (e.g., `inflate_overlap()`, `add_misspec()`)
  - Adding a new block in `sim4_scenarios.R` to call `run_scenario()` with your scenario key
- Set `use_bootstrap_se <- TRUE` in `sim4_config.R` to use bootstrap SE for HR/RMST (and optionally adjust `n_bootstrap_iterations`).


## 📚 Citation

If you use or adapt this simulation framework, please cite the following manuscript:

> **Chen HY**, Anand TV, Zhang L, Hripcsak G.  
> *When Does IPCW Help? Simulation and Real-World Evidence on Censoring Adjustment in Survival Analysis.*  
> medRxiv 2025. doi: [https://doi.org/10.1101/2025.08.22.25334253](https://doi.org/10.1101/2025.08.22.25334253)

This Simulation 4 framework extends the original study by introducing **high-dimensional confounding and informative censoring mechanisms**, with new scenarios for model misspecification and unmeasured confounding.

---

*Last updated: January 2026*  