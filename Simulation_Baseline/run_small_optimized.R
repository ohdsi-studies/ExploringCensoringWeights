############################################################
## run_small_optimized.R — Optimized version with:
## - Parallel execution option
## - Inline RMST SE computation
## - Better memory management
############################################################

## ---- Load required packages ----
if (!require(parallel)) install.packages("parallel")
library(parallel)

## ---- Combined estimator function ----
Wcombo_fun <- function(W_IPTW, W_IPCW) {
  W <- W_IPTW * W_IPCW
  # W / mean(W) # don't normalize
  return(W)
}

## ---- Single replication function (for parallel execution) ----
run_single_replication <- function(r, phiZ, ic_scale, ph_violation, 
                                   n_per_rep, p_total, p_continuous,
                                   target_treated_prop, target_censoring,
                                   k_event, lambda_event, k_censor, lambda_censor,
                                   beta_Z_early, beta_Z_late, t_PH,
                                   trunc_lo, trunc_hi, tt, tau_primary) {
  
  # Set seed for reproducibility (using L'Ecuyer-CMRG for parallel)
  set.seed(20251124 + r)
  
  Xdf <- gen_covariates(n_per_rep, p_total, p_continuous)
  
  # Treatment
  eta_tmp <- build_true_etas(Xdf, rep(0, n_per_rep),
                             phi_Z = phiZ, ic_scale = ic_scale)
  a0 <- calibrate_logit_intercept(eta_tmp$etaZ_wo_a0, target_treated_prop)
  pZ <- plogis(a0 + eta_tmp$etaZ_wo_a0)
  Z  <- rbinom(n_per_rep, 1, pZ)
  
  etas    <- build_true_etas(Xdf, Z, phi_Z = phiZ, ic_scale = ic_scale)
  etaT    <- etas$etaT
  etaC    <- etas$etaC_wo_g0
  eta_base <- etas$eta_base
  
  if (ph_violation) {
    etaT_early <- eta_base + beta_Z_early * Z
    etaT_late  <- eta_base + beta_Z_late  * Z
    etaT_for_calib <- etaT_early
  } else {
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
    Ttime <- rpiecewise_exp_loglin(
      n       = n_per_rep,
      lambda  = lambda_event,
      eta1    = etaT_early,
      eta2    = etaT_late,
      t_star  = t_PH
    )
  } else {
    Ttime <- rweibull_loglin(n_per_rep, k_event, lambda_event, etaT)
  }
  
  Ctime <- rweibull_loglin(n_per_rep, k_censor, lambda_censor, g0 + etaC)
  
  ti <- pmin(Ttime, Ctime)
  di <- as.numeric(Ttime <= Ctime)
  
  cens_prop <- mean(di == 0)
  
  dat <- data.frame(ID = 1:n_per_rep, D = Z, ti, di, xi = Ttime, ci = Ctime, Xdf)
  
  # Initialize result structure
  result <- list(
    surv_unadj_0 = numeric(length(tt)),
    surv_unadj_1 = numeric(length(tt)),
    surv_IPTW_0  = numeric(length(tt)),
    surv_IPTW_1  = numeric(length(tt)),
    surv_IPCW_0  = numeric(length(tt)),
    surv_IPCW_1  = numeric(length(tt)),
    surv_combo_0 = numeric(length(tt)),
    surv_combo_1 = numeric(length(tt)),
    se_surv_unadj_0 = numeric(length(tt)),
    se_surv_unadj_1 = numeric(length(tt)),
    se_surv_IPTW_0  = numeric(length(tt)),
    se_surv_IPTW_1  = numeric(length(tt)),
    se_surv_IPCW_0  = numeric(length(tt)),
    se_surv_IPCW_1  = numeric(length(tt)),
    se_surv_combo_0 = numeric(length(tt)),
    se_surv_combo_1 = numeric(length(tt)),
    beta_unadj = NA_real_,
    beta_iptw   = NA_real_,
    beta_ipcw   = NA_real_,
    beta_combo  = NA_real_,
    se_unadj    = NA_real_,
    se_iptw     = NA_real_,
    se_ipcw     = NA_real_,
    se_combo    = NA_real_,
    se_rmst_unadj = NA_real_,
    se_rmst_iptw  = NA_real_,
    se_rmst_ipcw  = NA_real_,
    se_rmst_combo = NA_real_,
    cens_prop = cens_prop
  )
  
  ## ----- Unadjusted -----
  s_un <- calc.surv.unadj.stratified(dat$ti, dat$di, tt, dat, tau = tau_primary)
  result$surv_unadj_0 <- s_un[s_un$D == 0, "survival"]
  result$surv_unadj_1 <- s_un[s_un$D == 1, "survival"]
  result$se_surv_unadj_0 <- s_un[s_un$D == 0, "se"]
  result$se_surv_unadj_1 <- s_un[s_un$D == 1, "se"]
  
  beta_result <- calc.beta(dat$ti, dat$di, dat)
  result$beta_unadj <- beta_result["coef"]
  result$se_unadj <- beta_result["se"]
  
  # Compute RMST SE from survfit rmean.std.err
  rmean_se_unadj_0 <- s_un[s_un$D == 0, "rmean_se"][1]
  rmean_se_unadj_1 <- s_un[s_un$D == 1, "rmean_se"][1]
  if (is.finite(rmean_se_unadj_0) && is.finite(rmean_se_unadj_1)) {
    result$se_rmst_unadj <- sqrt(rmean_se_unadj_0^2 + rmean_se_unadj_1^2)
  }
  
  ## ----- IPTW -----
  X_obs_ps <- Xdf[, c("X1", "X2"), drop = FALSE]
  iptw <- compute_iptw_variants(dat$D, X_obs_ps, trunc_lo, trunc_hi)
  dat$W_IPTW <- iptw$stab_trunc
  
  s_i <- calc.surv.IPW.stratified(dat$ti, dat$di, tt, dat$W_IPTW, dat, tau = tau_primary)
  result$surv_IPTW_0 <- s_i[s_i$D == 0, "survival"]
  result$surv_IPTW_1 <- s_i[s_i$D == 1, "survival"]
  result$se_surv_IPTW_0 <- s_i[s_i$D == 0, "se"]
  result$se_surv_IPTW_1 <- s_i[s_i$D == 1, "se"]
  
  beta_result <- calc.beta.IPW(dat$ti, dat$di, dat$W_IPTW, dat)
  result$beta_iptw <- beta_result["coef"]
  result$se_iptw <- beta_result["se"]
  
  # Compute RMST SE from survfit rmean.std.err
  rmean_se_iptw_0 <- s_i[s_i$D == 0, "rmean_se"][1]
  rmean_se_iptw_1 <- s_i[s_i$D == 1, "rmean_se"][1]
  if (is.finite(rmean_se_iptw_0) && is.finite(rmean_se_iptw_1)) {
    result$se_rmst_iptw <- sqrt(rmean_se_iptw_0^2 + rmean_se_iptw_1^2)
  }
  
  ## ----- IPCW -----
  cut.times <- sort(unique(round(dat$ti, 2)))
  dat.long <- transform.data(dat, cut.times)
  
  if (phiZ > 0) {
    formC <- Surv(Tstart, ti, censored, type = "counting") ~ D + X3 + X4
  } else {
    formC <- Surv(Tstart, ti, censored, type = "counting") ~ X3 + X4
  }
  
  C0 <- coxph(Surv(Tstart, ti, censored, type = "counting") ~ 1,
              data = dat.long, control = coxph.control(timefix = FALSE))
  CZ <- coxph(formC, data = dat.long, control = coxph.control(timefix = FALSE))
  
  dat.long <- calc.IPCW_fast(C0, CZ, dat.long, p_trunc = trunc_hi)
  dat.long$W_IPCW <- dat.long$WStab
  
  s_c <- calc.surv.IPCW.stratified(dat.long$Tstart, dat.long$ti, dat.long$di,
                                   tt, dat.long$W_IPCW, dat.long, tau = tau_primary)
  result$surv_IPCW_0 <- s_c[s_c$D == 0, "survival"]
  result$surv_IPCW_1 <- s_c[s_c$D == 1, "survival"]
  result$se_surv_IPCW_0 <- s_c[s_c$D == 0, "se"]
  result$se_surv_IPCW_1 <- s_c[s_c$D == 1, "se"]
  
  beta_result <- calc.beta.IPCW(
    dat.long$Tstart, dat.long$ti, dat.long$di,
    dat.long$W_IPCW, dat.long
  )
  result$beta_ipcw <- beta_result["coef"]
  result$se_ipcw <- beta_result["se"]
  
  # Compute RMST SE from survfit rmean.std.err
  rmean_se_ipcw_0 <- s_c[s_c$D == 0, "rmean_se"][1]
  rmean_se_ipcw_1 <- s_c[s_c$D == 1, "rmean_se"][1]
  if (is.finite(rmean_se_ipcw_0) && is.finite(rmean_se_ipcw_1)) {
    result$se_rmst_ipcw <- sqrt(rmean_se_ipcw_0^2 + rmean_se_ipcw_1^2)
  }
  
  ## ----- Combined -----
  dat.long$W_combo <- Wcombo_fun(
    W_IPTW = dat$W_IPTW[dat.long$ID],
    W_IPCW = dat.long$W_IPCW
  )
  
  s_comb <- calc.surv.comb.stratified(
    dat.long$Tstart, dat.long$ti, dat.long$di,
    tt, dat.long$W_combo, dat.long, tau = tau_primary
  )
  
  result$surv_combo_0 <- s_comb[s_comb$D == 0, "survival"]
  result$surv_combo_1 <- s_comb[s_comb$D == 1, "survival"]
  result$se_surv_combo_0 <- s_comb[s_comb$D == 0, "se"]
  result$se_surv_combo_1 <- s_comb[s_comb$D == 1, "se"]
  
  beta_result <- calc.beta.comb(
    dat.long$Tstart, dat.long$ti, dat.long$di,
    dat.long$W_combo, dat.long
  )
  result$beta_combo <- beta_result["coef"]
  result$se_combo <- beta_result["se"]
  
  # Compute RMST SE from survfit rmean.std.err
  rmean_se_combo_0 <- s_comb[s_comb$D == 0, "rmean_se"][1]
  rmean_se_combo_1 <- s_comb[s_comb$D == 1, "rmean_se"][1]
  if (is.finite(rmean_se_combo_0) && is.finite(rmean_se_combo_1)) {
    result$se_rmst_combo <- sqrt(rmean_se_combo_0^2 + rmean_se_combo_1^2)
  }
  
  return(result)
}

############################################################
## Optimized runner with parallel option
############################################################
run_small_scenario_optimized <- function(phiZ, ic_scale, scen_name, ic_label,
                                         ph_violation = FALSE,
                                         use_parallel = TRUE,
                                         n_cores = NULL) {
  
  message(sprintf("Running %s (%s IC) phiZ=%.2f",
                  scen_name, ic_label, phiZ))
  
  # Detect cores if not specified
  if (use_parallel && is.null(n_cores)) {
    n_cores <- max(1, detectCores() - 1)
    message(sprintf("Using %d cores for parallel execution", n_cores))
  }
  
  # Set up parallel RNG for reproducibility
  if (use_parallel && n_cores > 1) {
    RNGkind("L'Ecuyer-CMRG")
  }
  
  # Prepare arguments for single replication function
  rep_args <- list(
    phiZ = phiZ,
    ic_scale = ic_scale,
    ph_violation = ph_violation,
    n_per_rep = n_per_rep,
    p_total = p_total,
    p_continuous = p_continuous,
    target_treated_prop = target_treated_prop,
    target_censoring = target_censoring,
    k_event = k_event,
    lambda_event = lambda_event,
    k_censor = k_censor,
    lambda_censor = lambda_censor,
    beta_Z_early = beta_Z_early,
    beta_Z_late = beta_Z_late,
    t_PH = t_PH,
    trunc_lo = trunc_lo,
    trunc_hi = trunc_hi,
    tt = tt,
    tau_primary = tau_primary
  )
  
  # Run replications
  if (use_parallel && n_cores > 1) {
    message("Running replications in parallel...")
    results_list <- mclapply(1:n_reps, function(r) {
      do.call(run_single_replication, c(list(r = r), rep_args))
    }, mc.cores = n_cores)
  } else {
    message("Running replications sequentially...")
    pb <- txtProgressBar(0, n_reps, style = 3)
    results_list <- vector("list", n_reps)
    for (r in 1:n_reps) {
      results_list[[r]] <- do.call(run_single_replication, c(list(r = r), rep_args))
      setTxtProgressBar(pb, r)
    }
    close(pb)
  }
  
  # Combine results
  message("Combining results...")
  surv_unadj_0 <- do.call(rbind, lapply(results_list, function(x) x$surv_unadj_0))
  surv_unadj_1 <- do.call(rbind, lapply(results_list, function(x) x$surv_unadj_1))
  surv_IPTW_0  <- do.call(rbind, lapply(results_list, function(x) x$surv_IPTW_0))
  surv_IPTW_1  <- do.call(rbind, lapply(results_list, function(x) x$surv_IPTW_1))
  surv_IPCW_0  <- do.call(rbind, lapply(results_list, function(x) x$surv_IPCW_0))
  surv_IPCW_1  <- do.call(rbind, lapply(results_list, function(x) x$surv_IPCW_1))
  surv_combo_0 <- do.call(rbind, lapply(results_list, function(x) x$surv_combo_0))
  surv_combo_1 <- do.call(rbind, lapply(results_list, function(x) x$surv_combo_1))
  
  se_surv_unadj_0 <- do.call(rbind, lapply(results_list, function(x) x$se_surv_unadj_0))
  se_surv_unadj_1 <- do.call(rbind, lapply(results_list, function(x) x$se_surv_unadj_1))
  se_surv_IPTW_0  <- do.call(rbind, lapply(results_list, function(x) x$se_surv_IPTW_0))
  se_surv_IPTW_1  <- do.call(rbind, lapply(results_list, function(x) x$se_surv_IPTW_1))
  se_surv_IPCW_0  <- do.call(rbind, lapply(results_list, function(x) x$se_surv_IPCW_0))
  se_surv_IPCW_1  <- do.call(rbind, lapply(results_list, function(x) x$se_surv_IPCW_1))
  se_surv_combo_0 <- do.call(rbind, lapply(results_list, function(x) x$se_surv_combo_0))
  se_surv_combo_1 <- do.call(rbind, lapply(results_list, function(x) x$se_surv_combo_1))
  
  beta_unadj <- sapply(results_list, function(x) x$beta_unadj)
  beta_iptw  <- sapply(results_list, function(x) x$beta_iptw)
  beta_ipcw  <- sapply(results_list, function(x) x$beta_ipcw)
  beta_combo <- sapply(results_list, function(x) x$beta_combo)
  
  se_unadj <- sapply(results_list, function(x) x$se_unadj)
  se_iptw  <- sapply(results_list, function(x) x$se_iptw)
  se_ipcw  <- sapply(results_list, function(x) x$se_ipcw)
  se_combo <- sapply(results_list, function(x) x$se_combo)
  
  se_rmst_unadj <- sapply(results_list, function(x) x$se_rmst_unadj)
  se_rmst_iptw  <- sapply(results_list, function(x) x$se_rmst_iptw)
  se_rmst_ipcw  <- sapply(results_list, function(x) x$se_rmst_ipcw)
  se_rmst_combo <- sapply(results_list, function(x) x$se_rmst_combo)
  
  cens_prop <- sapply(results_list, function(x) x$cens_prop)
  
  # Choose which "truth" to attach
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
