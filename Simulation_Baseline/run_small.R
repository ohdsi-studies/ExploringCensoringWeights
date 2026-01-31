############################################################
## run_small.R — Sim 1 & 2 with 4 covariates, including:
## - 3 IC levels (none, weak, strong)
## - Combined estimator IPTW*IPCW
## - Truth stored in each result object
## - Supports Sim4-style RMST/HR plots
############################################################

rm(list = ls())

## ---- Load small config ----
source("sim_small_config.R")

## ---- Load core machinery ----
source("helper_D_censVars.R")
source("stratified_survival_curves.R")
source("sim_small_utils.R")

library(survival)
library(glmnet)
library(dplyr)
library(tidyr)


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
## Runner with IC levels
############################################################
run_small_scenario <- function(phiZ, ic_scale, scen_name, ic_label,
                               ph_violation = FALSE) {
  
  message(sprintf("Running %s (%s IC) phiZ=%.2f",
                  scen_name, ic_label, phiZ))
  
  surv_unadj_0 <- matrix(NA,n_reps,length(tt))
  surv_unadj_1 <- matrix(NA,n_reps,length(tt))
  surv_IPTW_0  <- matrix(NA,n_reps,length(tt))
  surv_IPTW_1  <- matrix(NA,n_reps,length(tt))
  surv_IPCW_0  <- matrix(NA,n_reps,length(tt))
  surv_IPCW_1  <- matrix(NA,n_reps,length(tt))
  surv_combo_0 <- matrix(NA,n_reps,length(tt))
  surv_combo_1 <- matrix(NA,n_reps,length(tt))
  
  # Storage for survival curve SEs
  se_surv_unadj_0 <- matrix(NA,n_reps,length(tt))
  se_surv_unadj_1 <- matrix(NA,n_reps,length(tt))
  se_surv_IPTW_0 <- matrix(NA,n_reps,length(tt))
  se_surv_IPTW_1 <- matrix(NA,n_reps,length(tt))
  se_surv_IPCW_0 <- matrix(NA,n_reps,length(tt))
  se_surv_IPCW_1 <- matrix(NA,n_reps,length(tt))
  se_surv_combo_0 <- matrix(NA,n_reps,length(tt))
  se_surv_combo_1 <- matrix(NA,n_reps,length(tt))
  
  beta_unadj <- beta_iptw <- beta_ipcw <- beta_combo <- numeric(n_reps)
  se_unadj <- se_iptw <- se_ipcw <- se_combo <- numeric(n_reps)
  se_rmst_unadj <- se_rmst_iptw <- se_rmst_ipcw <- se_rmst_combo <- numeric(n_reps)
  
  cens_prop <- numeric(n_reps)
  
  pb <- txtProgressBar(0,n_reps,style=3)
  
  for (r in 1:n_reps) {
    
    Xdf <- gen_covariates(n_per_rep, p_total, p_continuous)
    
    # Treatment
    eta_tmp <- build_true_etas(Xdf, rep(0,n_per_rep),
                               phi_Z=phiZ, ic_scale=ic_scale)
    a0 <- calibrate_logit_intercept(eta_tmp$etaZ_wo_a0, target_treated_prop)
    pZ <- plogis(a0 + eta_tmp$etaZ_wo_a0)
    Z  <- rbinom(n_per_rep,1,pZ)
    
    etas    <- build_true_etas(Xdf, Z, phi_Z = phiZ, ic_scale = ic_scale)
    etaT    <- etas$etaT              # PH version (for PH runs)
    etaC    <- etas$etaC_wo_g0
    eta_base <- etas$eta_base         # X * beta (no Z)
    
    if (ph_violation) {
      # Piecewise treatment effect: early vs late
      etaT_early <- eta_base + beta_Z_early * Z
      etaT_late  <- eta_base + beta_Z_late  * Z
      etaT_for_calib <- etaT_early         # use early effect for censor calibration
    } else {
      # Standard PH case
      etaT_early     <- etaT
      etaT_late      <- etaT
      etaT_for_calib <- etaT
    }
    
    g0 <- calibrate_censor_intercept(
      etaC, k_censor, lambda_censor,
      etaT_for_calib, k_event, lambda_event,
      target_censoring
    )
    
    if (ph_violation && abs(k_event - 1) < 1e-8) {
      # Piecewise exponential with different treatment effect after t_PH
      Ttime <- rpiecewise_exp_loglin(
        n       = n_per_rep,
        lambda  = lambda_event,
        eta1    = etaT_early,
        eta2    = etaT_late,
        t_star  = t_PH
      )
    } else {
      # Default PH Weibull generator (your original)
      Ttime <- rweibull_loglin(n_per_rep, k_event, lambda_event, etaT)
    }
    
    Ctime <- rweibull_loglin(n_per_rep, k_censor, lambda_censor, g0 + etaC)
    
    ti <- pmin(Ttime,Ctime)
    di <- as.numeric(Ttime<=Ctime)
    
    cens_prop[r] <- mean(di==0)
    
    dat <- data.frame(ID=1:n_per_rep, D=Z, ti, di, xi=Ttime, ci=Ctime, Xdf)
    
    ## ----- Unadjusted -----
    s_un <- calc.surv.unadj.stratified(dat$ti,dat$di,tt,dat,tau=tau_primary)
    surv_unadj_0[r,] <- s_un[s_un$D==0,"survival"]
    surv_unadj_1[r,] <- s_un[s_un$D==1,"survival"]
    se_surv_unadj_0[r,] <- s_un[s_un$D==0,"se"]
    se_surv_unadj_1[r,] <- s_un[s_un$D==1,"se"]
    # Extract rmean SEs for RMST difference SE calculation
    rmean_se_unadj_0 <- s_un[s_un$D==0,"rmean_se"][1]
    rmean_se_unadj_1 <- s_un[s_un$D==1,"rmean_se"][1]
    if (is.finite(rmean_se_unadj_0) && is.finite(rmean_se_unadj_1)) {
      se_rmst_unadj[r] <- sqrt(rmean_se_unadj_0^2 + rmean_se_unadj_1^2)
    } else {
      se_rmst_unadj[r] <- NA_real_
    }
    beta_result <- calc.beta(dat$ti,dat$di,dat)
    beta_unadj[r] <- beta_result["coef"]
    se_unadj[r] <- beta_result["se"]
    
    ## ----- IPTW -----
    X_obs_ps <- Xdf[,c("X1","X2"),drop=FALSE]
    iptw <- compute_iptw_variants(dat$D,X_obs_ps,trunc_lo,trunc_hi)
    dat$W_IPTW <- iptw$stab_trunc
    
    s_i <- calc.surv.IPW.stratified(dat$ti,dat$di,tt,dat$W_IPTW,dat,tau=tau_primary)
    surv_IPTW_0[r,] <- s_i[s_i$D==0,"survival"]
    surv_IPTW_1[r,] <- s_i[s_i$D==1,"survival"]
    se_surv_IPTW_0[r,] <- s_i[s_i$D==0,"se"]
    se_surv_IPTW_1[r,] <- s_i[s_i$D==1,"se"]
    # Extract rmean SEs for RMST difference SE calculation
    rmean_se_iptw_0 <- s_i[s_i$D==0,"rmean_se"][1]
    rmean_se_iptw_1 <- s_i[s_i$D==1,"rmean_se"][1]
    if (is.finite(rmean_se_iptw_0) && is.finite(rmean_se_iptw_1)) {
      se_rmst_iptw[r] <- sqrt(rmean_se_iptw_0^2 + rmean_se_iptw_1^2)
    } else {
      se_rmst_iptw[r] <- NA_real_
    }
    beta_result <- calc.beta.IPW(dat$ti,dat$di,dat$W_IPTW,dat)
    beta_iptw[r] <- beta_result["coef"]
    se_iptw[r] <- beta_result["se"]
    
    ## ----- IPCW -----
    cut.times <- sort(unique(round(dat$ti, 2))) # slightly coarser cut times than all unique times
    dat.long <- transform.data(dat,cut.times)
    
    if (phiZ>0) {
      formC <- Surv(Tstart,ti,censored,type="counting") ~ D + X3 + X4
    } else {
      formC <- Surv(Tstart,ti,censored,type="counting") ~ X3 + X4
    }
    
    C0 <- coxph(Surv(Tstart,ti,censored,type="counting") ~ 1,
                data=dat.long, control=coxph.control(timefix=FALSE))
    CZ <- coxph(formC,data=dat.long, control=coxph.control(timefix=FALSE))
    
    dat.long <- calc.IPCW_fast(C0,CZ,dat.long,p_trunc=trunc_hi)
    dat.long$W_IPCW <- dat.long$WStab
    
    s_c <- calc.surv.IPCW.stratified(dat.long$Tstart,dat.long$ti,dat.long$di,
                                     tt,dat.long$W_IPCW,dat.long,tau=tau_primary)
    surv_IPCW_0[r,] <- s_c[s_c$D==0,"survival"]
    surv_IPCW_1[r,] <- s_c[s_c$D==1,"survival"]
    se_surv_IPCW_0[r,] <- s_c[s_c$D==0,"se"]
    se_surv_IPCW_1[r,] <- s_c[s_c$D==1,"se"]
    # Extract rmean SEs for RMST difference SE calculation
    rmean_se_ipcw_0 <- s_c[s_c$D==0,"rmean_se"][1]
    rmean_se_ipcw_1 <- s_c[s_c$D==1,"rmean_se"][1]
    if (is.finite(rmean_se_ipcw_0) && is.finite(rmean_se_ipcw_1)) {
      se_rmst_ipcw[r] <- sqrt(rmean_se_ipcw_0^2 + rmean_se_ipcw_1^2)
    } else {
      se_rmst_ipcw[r] <- NA_real_
    }
    
    beta_result <- calc.beta.IPCW(
      dat.long$Tstart,dat.long$ti,dat.long$di,
      dat.long$W_IPCW,dat.long
    )
    beta_ipcw[r] <- beta_result["coef"]
    se_ipcw[r] <- beta_result["se"]
    
    ## ----- Combined -----
    dat.long$W_combo <- Wcombo_fun(
      W_IPTW=dat$W_IPTW[dat.long$ID],
      W_IPCW=dat.long$W_IPCW
    )
    
    s_comb <- calc.surv.comb.stratified(
      dat.long$Tstart,dat.long$ti,dat.long$di,
      tt, dat.long$W_combo, dat.long, tau=tau_primary
    )
    
    surv_combo_0[r,] <- s_comb[s_comb$D==0,"survival"]
    surv_combo_1[r,] <- s_comb[s_comb$D==1,"survival"]
    se_surv_combo_0[r,] <- s_comb[s_comb$D==0,"se"]
    se_surv_combo_1[r,] <- s_comb[s_comb$D==1,"se"]
    # Extract rmean SEs for RMST difference SE calculation
    rmean_se_combo_0 <- s_comb[s_comb$D==0,"rmean_se"][1]
    rmean_se_combo_1 <- s_comb[s_comb$D==1,"rmean_se"][1]
    if (is.finite(rmean_se_combo_0) && is.finite(rmean_se_combo_1)) {
      se_rmst_combo[r] <- sqrt(rmean_se_combo_0^2 + rmean_se_combo_1^2)
    } else {
      se_rmst_combo[r] <- NA_real_
    }
    
    beta_result <- calc.beta.comb(
      dat.long$Tstart,dat.long$ti,dat.long$di,
      dat.long$W_combo,dat.long
    )
    beta_combo[r] <- beta_result["coef"]
    se_combo[r] <- beta_result["se"]
    
    setTxtProgressBar(pb,r)
  }
  close(pb)
  
  ## RMST SEs are now computed directly from survfit rmean.std.err during the main loop above
  ## No need for separate computation using compute_rmst_with_se
  
  ## Choose which "truth" to attach based on PH vs non-PH scenario
  if (ph_violation) {
    surv_true_0_use   <- surv_true_0_PH
    surv_true_1_use   <- surv_true_1_PH
    beta_true_use     <- beta_true_PH
    rmst_true_use     <- rmst_true_diff_PH
  } else {
    surv_true_0_use   <- surv_true_0_fixed
    surv_true_1_use   <- surv_true_1_fixed
    beta_true_use     <- beta_true_fixed
    rmst_true_use     <- rmst_true_diff
  }
  
  
  list(
    tt = tt,
    ic_label = ic_label,
    phiZ = phiZ,
    
    surv_strat = list(
      real_0 = matrix(rep(surv_true_0_use, n_reps), nrow = n_reps, byrow = TRUE),
      real_1 = matrix(rep(surv_true_1_use, n_reps), nrow = n_reps, byrow = TRUE),
      
      
      unadj_0 = surv_unadj_0,
      unadj_1 = surv_unadj_1,
      iptw_0  = surv_IPTW_0,
      iptw_1  = surv_IPTW_1,
      ipcw_0  = surv_IPCW_0,
      ipcw_1  = surv_IPCW_1,
      combo_0 = surv_combo_0,
      combo_1 = surv_combo_1
    ),
    
    beta = list(
      true  = rep(beta_true_use, n_reps),
      unadj = beta_unadj,
      iptw  = beta_iptw,
      ipcw  = beta_ipcw,
      combo = beta_combo
    ),
    
    se = list(
      unadj = se_unadj,
      iptw  = se_iptw,
      ipcw  = se_ipcw,
      combo = se_combo
    ),
    
    se_rmst = list(
      unadj = se_rmst_unadj,
      iptw  = se_rmst_iptw,
      ipcw  = se_rmst_ipcw,
      combo = se_rmst_combo
    ),

    # Pointwise SEs for the survival curves at each time in tt
    # (useful for survival curve CIs/bands later)
    se_surv_strat = list(
      unadj_0 = se_surv_unadj_0,
      unadj_1 = se_surv_unadj_1,
      iptw_0  = se_surv_IPTW_0,
      iptw_1  = se_surv_IPTW_1,
      ipcw_0  = se_surv_IPCW_0,
      ipcw_1  = se_surv_IPCW_1,
      combo_0 = se_surv_combo_0,
      combo_1 = se_surv_combo_1
    ),
    
    rmst_true_diff = rmst_true_use,
    cens_prop      = cens_prop
  )
}


############################################################
## Run all 3 scenarios (sim1, sim2, sim3) across 3 IC levels
############################################################

ic_levels <- c("none" = 0, "weak" = 0.5, "strong" = 1.2)

for (ic_label in names(ic_levels)) {
  for (scen in c("sim1", "sim2", "sim3")) {
  # for (scen in c("sim2")) { line for debugging
    
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
    # Run the scenario
    # ------------------------------------------------------
    res <- run_small_scenario(
      phiZ        = phiZ_val,
      ic_scale    = ic_levels[ic_label],
      scen_name   = paste0("Scenario_", scen),
      ic_label    = ic_label,
      ph_violation = ph_flag     # <--- NEW argument
    )
    
    # ------------------------------------------------------
    # Save results
    # ------------------------------------------------------
    outfile <- file.path(out_dir, paste0("res_", scen, "_", ic_label, ".rdata"))
    save(res, file = outfile)
  }
}

message("All results saved in ", out_dir)
