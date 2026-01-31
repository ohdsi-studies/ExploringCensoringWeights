############################################################
## run_small_parallel.R — Optimized version with parallel execution
## - Uses parallel processing for replications (4-8x speedup)
## - Inline RMST SE computation (eliminates second pass)
## - Same output structure as run_small.R
##
## Usage: Same as run_small.R, but with parallel execution enabled
## To disable parallel: set use_parallel = FALSE in the function call
############################################################

rm(list = ls())

## ---- Load small config ----
source("sim_small_config.R")

## ---- Load core machinery ----
source("helper_D_censVars.R")
source("stratified_survival_curves.R")
source("sim_small_utils.R")

## ---- Load optimized simulation function ----
source("run_small_optimized.R")

library(survival)
library(glmnet)
library(dplyr)
library(tidyr)
library(parallel)


############################################################
## Gold-standard "truth" under PH violation (small sim)
## Uses rpiecewise_exp_loglin() with early / late tx effects
############################################################
get_true_counterfactual_PH <- function(
    Xdf           = Xdf_truth,
    tt,
    scenario      = "4A",
    k_event,
    lambda_event,
    beta_Z_early,
    beta_Z_late,
    t_PH
) {
  message("=== computing marginal beta (PH-violation truth) ===")
  
  if (abs(k_event - 1) > 1e-8) {
    stop("get_true_counterfactual_PH assumes k_event = 1 (exponential baseline).")
  }
  
  n <- nrow(Xdf)
  
  ## 1. Build baseline LP (no treatment) via build_true_etas
  Z0 <- rep(0, n)
  etas0 <- build_true_etas(
    Xdf     = Xdf,
    Z       = Z0,
    scenario = scenario,
    phi_Z   = 0,
    ic_scale = 0
  )
  
  if (!("eta_base" %in% names(etas0))) {
    stop("build_true_etas must return 'eta_base' (baseline event LP).")
  }
  
  eta_base <- etas0$eta_base  # X * beta_event, no treatment
  
  ## 2. Define early / late LPs for each treatment arm
  # Control arm (D = 0): no treatment effect at any time
  eta0_early <- eta_base
  eta0_late  <- eta_base
  
  # Treated arm (D = 1): beta_Z_early before t_PH, beta_Z_late after t_PH
  eta1_early <- eta_base + beta_Z_early
  eta1_late  <- eta_base + beta_Z_late
  
  ## 3. Simulate counterfactual event times under non-PH
  T0 <- rpiecewise_exp_loglin(
    n      = n,
    lambda = lambda_event,
    eta1   = eta0_early,
    eta2   = eta0_late,
    t_star = t_PH
  )
  
  T1 <- rpiecewise_exp_loglin(
    n      = n,
    lambda = lambda_event,
    eta1   = eta1_early,
    eta2   = eta1_late,
    t_star = t_PH
  )
  
  ## 4. Build a counterfactual trial dataset
  dat_cf <- rbind(
    data.frame(D = 0, time = T0, status = 1),
    data.frame(D = 1, time = T1, status = 1)
  )
  
  ## 5. Cox model for "causal" HR under this DGP
  cox_cf   <- coxph(Surv(time, status) ~ D, data = dat_cf)
  beta_true <- coef(cox_cf)[["D"]]
  
  ## 6. Causal RMST difference and survival curves
  fit_cf <- survfit(Surv(time, status) ~ D, data = dat_cf)
  sfit   <- summary(fit_cf, times = tt, extend = TRUE)
  
  surv_real <- sfit$surv
  surv_0    <- sfit$surv[sfit$strata == "D=0"]
  surv_1    <- sfit$surv[sfit$strata == "D=1"]
  
  dt       <- diff(c(0, tt))
  rmst_0   <- sum(surv_0 * dt)
  rmst_1   <- sum(surv_1 * dt)
  rmst_diff <- rmst_1 - rmst_0
  
  list(
    beta_true      = beta_true,
    rmst_true_diff = rmst_diff,
    surv_real      = surv_real,
    surv_0         = surv_0,
    surv_1         = surv_1
  )
}


############################################################
## Compute truth ONCE
############################################################
n_truth <- 1e6

Xdf_truth <- gen_covariates(
  n        = n_truth,
  p_total  = p_total,
  p_cont   = p_continuous
)

# truth for sim 1 and 2
truth_small <- get_true_counterfactual(
  Xdf = Xdf_truth,
  tt  = tt,
  scenario = "4A",
  k_event      = k_event,
  lambda_event = lambda_event
)

beta_true_fixed <- truth_small$beta_true
surv_true_0_fixed <- truth_small$surv_0
surv_true_1_fixed <- truth_small$surv_1
rmst_true_diff    <- truth_small$rmst_true_diff

## Truth for PH-violation scenario (Sim3)
truth_small_PH <- get_true_counterfactual_PH(
  Xdf          = Xdf_truth,
  tt           = tt,
  scenario     = "4A",           # irrelevant for sim_small, but required
  k_event      = k_event,
  lambda_event = lambda_event,
  beta_Z_early = beta_Z_early,   # you defined these above
  beta_Z_late  = beta_Z_late,
  t_PH         = t_PH
)

beta_true_PH        <- truth_small_PH$beta_true
surv_true_0_PH      <- truth_small_PH$surv_0
surv_true_1_PH      <- truth_small_PH$surv_1
rmst_true_diff_PH   <- truth_small_PH$rmst_true_diff


############################################################
## Combined estimator
############################################################
Wcombo_fun <- function(W_IPTW, W_IPCW) {
  W <- W_IPTW * W_IPCW
  # W / mean(W) # don't normalize
  return(W)
}


############################################################
## Configuration for parallel execution
############################################################
# Set to TRUE to use parallel execution (recommended)
# Set to FALSE to run sequentially (for debugging or if memory is limited)
USE_PARALLEL <- TRUE

# Number of cores to use (NULL = auto-detect, uses all cores - 1)
# You can also specify a number, e.g., N_CORES <- 4
N_CORES <- NULL

# Display parallel configuration
if (USE_PARALLEL) {
  if (is.null(N_CORES)) {
    n_cores_actual <- max(1, detectCores() - 1)
  } else {
    n_cores_actual <- N_CORES
  }
  message(sprintf("Parallel execution enabled: using %d cores", n_cores_actual))
} else {
  message("Parallel execution disabled: running sequentially")
}


############################################################
## Run all 3 scenarios (sim1, sim2, sim3) across 3 IC levels
############################################################

ic_levels <- c("none" = 0, "weak" = 0.5, "strong" = 1.2) # prev: 1.5, 3

for (ic_label in names(ic_levels)) {
  # for (scen in c("sim2")) {
  for (scen in c("sim1", "sim2", "sim3")) {
    
    # ------------------------------------------------------
    # Set phiZ depending on scenario
    # sim1 → nondiff censoring
    # sim2 → diff censoring
    # sim3 → PH violation scenario 
    # ------------------------------------------------------
    phiZ_val <- switch(
      scen,
      "sim1" = phi_Z_levels["sim1"],
      "sim2" = phi_Z_levels["sim2"],
      "sim3" = phi_Z_levels["sim2"],  # same as sim2 unless you want otherwise
      stop("Invalid scenario")
    )
    
    # ------------------------------------------------------
    # PH violation applies ONLY to scenario 3
    # ------------------------------------------------------
    ph_flag <- (scen == "sim3")
    
    # ------------------------------------------------------
    # Run the scenario using optimized parallel function
    # ------------------------------------------------------
    message(sprintf("\n=== Running %s (%s IC) ===", scen, ic_label))
    
    res <- run_small_scenario_optimized(
      phiZ        = phiZ_val,
      ic_scale    = ic_levels[ic_label],
      scen_name   = paste0("Scenario_", scen),
      ic_label    = ic_label,
      ph_violation = ph_flag,
      use_parallel = USE_PARALLEL,
      n_cores     = N_CORES
    )
    
    # ------------------------------------------------------
    # Save results
    # ------------------------------------------------------
    outfile <- file.path(out_dir, paste0("res_", scen, "_", ic_label, ".rdata"))
    save(res, file = outfile)
    message(sprintf("Saved results to %s", outfile))
  }
}

message("\n=== All results saved in ", out_dir, " ===")
