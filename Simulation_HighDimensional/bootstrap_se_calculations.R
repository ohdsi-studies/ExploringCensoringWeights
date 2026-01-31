## Bootstrap SE Calculation for HR and RMST
## Re-estimates IPTW and IPCW weights in each bootstrap iteration
##
## Note: This file should be sourced from sim4_run.R which already has:
##   - helper_D_censVars.R (transform.data, calc.IPCW_fast, etc.)
##   - stratified_survival_curves.R (calc.surv.*.stratified functions)
##   - sim4_utils_hd.R (compute_iptw_variants)
##   - All necessary libraries (survival, glmnet, etc.)

# Helper function to compute RMST from survival curve
compute_rmst <- function(tt, surv, tau) {
  idx <- which(tt <= tau)
  if (length(idx) == 0L) return(NA_real_)
  
  if (tau %in% tt) {
    time_trunc <- tt[idx]
    surv_trunc <- surv[idx]
  } else {
    surv_tau <- approx(tt, surv, xout = tau, method = "linear")$y
    time_trunc <- c(tt[idx], tau)
    surv_trunc <- c(surv[idx], surv_tau)
  }
  
  sum(diff(time_trunc) * (head(surv_trunc, -1) + tail(surv_trunc, -1)) / 2)
}

# Single bootstrap iteration
bootstrap_single_iteration <- function(dat, method, tt, tau, 
                                        X_obs_ps, X_obs_cens, scen_key,
                                        trunc_lo, trunc_hi) {
  
  n <- nrow(dat)
  
  # Resample with replacement
  boot_idx <- sample(n, replace = TRUE)
  dat_boot <- dat[boot_idx, ]
  
  # Check for edge cases
  if (length(unique(dat_boot$D)) < 2) {
    return(list(hr = NA_real_, rmst_diff = NA_real_))
  }
  
  if (sum(dat_boot$di) == 0) {
    return(list(hr = NA_real_, rmst_diff = NA_real_))
  }
  
  # Transform to long format (needed for IPCW and comb methods)
  dat_boot_long <- tryCatch({
    transform.data(dat_boot, cut.times = dat_boot$ti)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(dat_boot_long)) {
    return(list(hr = NA_real_, rmst_diff = NA_real_))
  }
  
  # Initialize weights
  iptw_weights <- NULL
  ipcw_weights <- NULL
  comb_weights <- NULL
  
  # Compute weights based on method
  if (method == "unadj") {
    # No weights needed
    
  } else if (method == "iptw") {
    # Re-estimate IPTW
    X_obs_ps_boot <- X_obs_ps[boot_idx, , drop = FALSE]
    iptw_boot <- tryCatch({
      compute_iptw_variants(
        Z = dat_boot$D,
        X_df = X_obs_ps_boot,
        trunc_lo = trunc_lo,
        trunc_hi = trunc_hi
      )
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(iptw_boot)) {
      return(list(hr = NA_real_, rmst_diff = NA_real_))
    }
    
    iptw_weights <- iptw_boot$stab_trunc
    
  } else if (method == "ipcw") {
    # Re-estimate IPCW
    # Step 1: Feature selection on wide format
    dat_boot$censored <- 1 - dat_boot$di
    x_wide <- as.matrix(dat_boot[, c("D", colnames(X_obs_cens)), drop = FALSE])
    y_wide <- with(dat_boot, survival::Surv(ti, censored))
    
    cvfit <- tryCatch({
      glmnet::cv.glmnet(
        x = x_wide, y = y_wide,
        family = "cox",
        alpha = 1,
        nfolds = 3
      )
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(cvfit)) {
      return(list(hr = NA_real_, rmst_diff = NA_real_))
    }
    
    b <- as.matrix(coef(cvfit, s = "lambda.1se"))
    keep <- rownames(b)[as.numeric(b[, 1]) != 0]
    
    form_C0 <- as.formula("Surv(Tstart, ti, censored, type='counting') ~ 1")
    
    if (length(keep) == 0) {
      rhs <- if (startsWith(scen_key, "4C")) "D" else "1"
    } else {
      rhs <- paste(keep, collapse = " + ")
    }
    
    form_CZ_sel <- as.formula(
      paste0("Surv(Tstart, ti, censored, type='counting') ~ ", rhs)
    )
    
    # Step 2: Fit C0 and CZ models on long format
    C0 <- tryCatch({
      coxph(form_C0, data = dat_boot_long, control = coxph.control(timefix = FALSE))
    }, error = function(e) {
      return(NULL)
    })
    
    CZ <- tryCatch({
      coxph(form_CZ_sel, data = dat_boot_long, control = coxph.control(timefix = FALSE))
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(C0) || is.null(CZ)) {
      return(list(hr = NA_real_, rmst_diff = NA_real_))
    }
    
    # Step 3: Compute IPCW weights
    dat_boot_long <- tryCatch({
      calc.IPCW_fast(C0, CZ, dat_boot_long, p_trunc = trunc_hi)
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(dat_boot_long)) {
      return(list(hr = NA_real_, rmst_diff = NA_real_))
    }
    
    ipcw_weights <- dat_boot_long$WStab
    
  } else if (method == "comb") {
    # Re-estimate both IPTW and IPCW, then combine
    
    # IPTW
    X_obs_ps_boot <- X_obs_ps[boot_idx, , drop = FALSE]
    iptw_boot <- tryCatch({
      compute_iptw_variants(
        Z = dat_boot$D,
        X_df = X_obs_ps_boot,
        trunc_lo = trunc_lo,
        trunc_hi = trunc_hi
      )
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(iptw_boot)) {
      return(list(hr = NA_real_, rmst_diff = NA_real_))
    }
    
    # IPCW (same as above)
    dat_boot$censored <- 1 - dat_boot$di
    x_wide <- as.matrix(dat_boot[, c("D", colnames(X_obs_cens)), drop = FALSE])
    y_wide <- with(dat_boot, survival::Surv(ti, censored))
    
    cvfit <- tryCatch({
      glmnet::cv.glmnet(
        x = x_wide, y = y_wide,
        family = "cox",
        alpha = 1,
        nfolds = 3
      )
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(cvfit)) {
      return(list(hr = NA_real_, rmst_diff = NA_real_))
    }
    
    b <- as.matrix(coef(cvfit, s = "lambda.1se"))
    keep <- rownames(b)[as.numeric(b[, 1]) != 0]
    
    form_C0 <- as.formula("Surv(Tstart, ti, censored, type='counting') ~ 1")
    
    if (length(keep) == 0) {
      rhs <- if (startsWith(scen_key, "4C")) "D" else "1"
    } else {
      rhs <- paste(keep, collapse = " + ")
    }
    
    form_CZ_sel <- as.formula(
      paste0("Surv(Tstart, ti, censored, type='counting') ~ ", rhs)
    )
    
    C0 <- tryCatch({
      coxph(form_C0, data = dat_boot_long, control = coxph.control(timefix = FALSE))
    }, error = function(e) {
      return(NULL)
    })
    
    CZ <- tryCatch({
      coxph(form_CZ_sel, data = dat_boot_long, control = coxph.control(timefix = FALSE))
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(C0) || is.null(CZ)) {
      return(list(hr = NA_real_, rmst_diff = NA_real_))
    }
    
    dat_boot_long <- tryCatch({
      calc.IPCW_fast(C0, CZ, dat_boot_long, p_trunc = trunc_hi)
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(dat_boot_long)) {
      return(list(hr = NA_real_, rmst_diff = NA_real_))
    }
    
    # Combine weights
    pZ <- iptw_boot$pZ
    # Match IPTW weights to long format by ID
    dat_boot_long$IPW_unstab_untrunc <- iptw_boot$unstab_untrunc[match(dat_boot_long$ID, dat_boot$ID)]
    # IPCW weights already computed in dat_boot_long$KZti
    dat_boot_long$IPCW_raw <- 1 / pmax(dat_boot_long$KZti, 1e-8)
    dat_boot_long$comb_raw <- dat_boot_long$IPW_unstab_untrunc * dat_boot_long$IPCW_raw
    dat_boot_long$IPTW_stab_factor <- ifelse(dat_boot_long$D == 1, pZ, 1 - pZ)
    dat_boot_long$comb_stab_untrunc <- dat_boot_long$comb_raw * dat_boot_long$IPTW_stab_factor * dat_boot_long$K0ti
    
    cap <- quantile(dat_boot_long$comb_stab_untrunc, trunc_hi, na.rm = TRUE)
    comb_weights <- pmin(dat_boot_long$comb_stab_untrunc, cap)
  }
  
  # Compute HR
  hr <- tryCatch({
    if (method == "unadj") {
      Cox <- coxph(Surv(ti, di) ~ D, data = dat_boot)
      coef(Cox)[["D"]]
    } else if (method == "iptw") {
      Cox <- coxph(Surv(ti, di) ~ as.factor(D), data = dat_boot, weights = iptw_weights)
      coef(Cox)[["as.factor(D)1"]]
    } else if (method == "ipcw") {
      Cox <- coxph(Surv(Tstart, ti, di, type = "counting") ~ D, 
                   data = dat_boot_long, weights = ipcw_weights,
                   control = coxph.control(timefix = FALSE))
      coef(Cox)[["D"]]
    } else if (method == "comb") {
      Cox <- coxph(Surv(Tstart, ti, di, type = "counting") ~ D,
                   data = dat_boot_long, weights = comb_weights,
                   control = coxph.control(timefix = FALSE))
      coef(Cox)[["D"]]
    } else {
      NA_real_
    }
  }, error = function(e) {
    NA_real_
  })
  
  # Compute RMST difference
  rmst_diff <- tryCatch({
    if (method == "unadj") {
      surv_strat_df <- calc.surv.unadj.stratified(
        times = dat_boot$ti,
        status = dat_boot$di,
        tt = tt,
        data = dat_boot
      )
      surv_0 <- surv_strat_df[surv_strat_df$D == 0, "survival"]
      surv_1 <- surv_strat_df[surv_strat_df$D == 1, "survival"]
      
    } else if (method == "iptw") {
      surv_strat_df <- calc.surv.IPW.stratified(
        times = dat_boot$ti,
        status = dat_boot$di,
        tt = tt,
        IPW.weights = iptw_weights,
        data = dat_boot
      )
      surv_0 <- surv_strat_df[surv_strat_df$D == 0, "survival"]
      surv_1 <- surv_strat_df[surv_strat_df$D == 1, "survival"]
      
    } else if (method == "ipcw") {
      surv_strat_df <- calc.surv.IPCW.stratified(
        Tstart = dat_boot_long$Tstart,
        Tstop = dat_boot_long$ti,
        status = dat_boot_long$di,
        tt = tt,
        IPCW.weights = ipcw_weights,
        data.long = dat_boot_long
      )
      surv_0 <- surv_strat_df[surv_strat_df$D == 0, "survival"]
      surv_1 <- surv_strat_df[surv_strat_df$D == 1, "survival"]
      
    } else if (method == "comb") {
      surv_strat_df <- calc.surv.comb.stratified(
        Tstart = dat_boot_long$Tstart,
        Tstop = dat_boot_long$ti,
        status = dat_boot_long$di,
        tt = tt,
        comb.weights = comb_weights,
        data.long = dat_boot_long
      )
      surv_0 <- surv_strat_df[surv_strat_df$D == 0, "survival"]
      surv_1 <- surv_strat_df[surv_strat_df$D == 1, "survival"]
      
    } else {
      return(list(hr = hr, rmst_diff = NA_real_))
    }
    
    rmst_0 <- compute_rmst(tt, surv_0, tau)
    rmst_1 <- compute_rmst(tt, surv_1, tau)
    rmst_1 - rmst_0
    
  }, error = function(e) {
    NA_real_
  })
  
  list(hr = hr, rmst_diff = rmst_diff)
}

# Main bootstrap SE calculation function
compute_bootstrap_se <- function(dat, method, tt, tau,
                                  X_obs_ps, X_obs_cens, scen_key,
                                  trunc_lo, trunc_hi,
                                  n_boot = 200, seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  hr_boot <- numeric(n_boot)
  rmst_boot <- numeric(n_boot)
  
  for (b in seq_len(n_boot)) {
    result <- bootstrap_single_iteration(
      dat = dat,
      method = method,
      tt = tt,
      tau = tau,
      X_obs_ps = X_obs_ps,
      X_obs_cens = X_obs_cens,
      scen_key = scen_key,
      trunc_lo = trunc_lo,
      trunc_hi = trunc_hi
    )
    
    hr_boot[b] <- result$hr
    rmst_boot[b] <- result$rmst_diff
  }
  
  # Compute SEs from bootstrap distribution
  se_hr <- sd(hr_boot, na.rm = TRUE)
  se_rmst <- sd(rmst_boot, na.rm = TRUE)
  
  # Check for too many failures
  n_failures_hr <- sum(!is.finite(hr_boot))
  n_failures_rmst <- sum(!is.finite(rmst_boot))
  
  if (n_failures_hr > n_boot * 0.1) {
    warning("More than 10% of bootstrap iterations failed for HR calculation (method: ", method, ")")
  }
  
  if (n_failures_rmst > n_boot * 0.1) {
    warning("More than 10% of bootstrap iterations failed for RMST calculation (method: ", method, ")")
  }
  
  # Return NA if all iterations failed
  if (n_failures_hr == n_boot) {
    se_hr <- NA_real_
  }
  
  if (n_failures_rmst == n_boot) {
    se_rmst <- NA_real_
  }
  
  list(
    se_hr = se_hr,
    se_rmst = se_rmst,
    hr_boot = hr_boot,
    rmst_boot = rmst_boot
  )
}
