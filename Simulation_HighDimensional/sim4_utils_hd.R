# ---- Utilities for HD simulation ----
library(MASS)
library(glmnet)
library(Matrix)
library(survival)
library(dplyr)

# Build HD covariate matrix with AR(1) correlation among continuous block and binary via probit

# this old version has no binary correlation
gen_hd_covariates_old <- function(n, p_total, p_cont, ar1_rho, p_bin, bin_prev) {
  stopifnot(p_total == p_cont + p_bin)
  X <- matrix(0, n, p_total)
  colnames(X) <- paste0("X", seq_len(p_total))
  
  if (p_cont > 0) {
    Sigma <- outer(1:p_cont, 1:p_cont, function(i,j) ar1_rho^abs(i-j))
    V <- mvrnorm(n, mu = rep(0, p_cont), Sigma = Sigma)
    X[, 1:p_cont] <- V
  }
  if (p_bin > 0) {
    # Probit thresholding for desired prevalence
    # threshold tau s.t. P(Z>tau) ~ bin_prev; Z~N(0,1) ⇒ tau = qnorm(1 - bin_prev)
    tau <- qnorm(1 - bin_prev)
    Zb  <- matrix(rnorm(n * p_bin), n, p_bin)
    X[, (p_cont+1):p_total] <- (Zb > tau) * 1
  }
  X <- as.data.frame(X)
  X
}

# Build HD covariate matrix with:
#  - AR(1) correlation among continuous block
#  - AR(1) correlation among binary block via probit thresholding
gen_hd_covariates <- function(n,
                              p_total,
                              p_cont,
                              ar1_rho,
                              p_bin,
                              bin_prev,
                              ar1_rho_bin = 0) {
  stopifnot(p_total == p_cont + p_bin)
  
  X <- matrix(0, n, p_total)
  colnames(X) <- paste0("X", seq_len(p_total))
  
  ## ----- Continuous block (multivariate normal with AR(1)) -----
  if (p_cont > 0) {
    Sigma_cont <- outer(1:p_cont, 1:p_cont,
                        function(i, j) ar1_rho^abs(i - j))
    V_cont <- mvrnorm(n, mu = rep(0, p_cont), Sigma = Sigma_cont)
    X[, 1:p_cont] <- V_cont
  }
  
  ## ----- Binary block (latent MVN with AR(1), probit threshold) -----
  if (p_bin > 0) {
    # threshold tau s.t. P(Z > tau) ~ bin_prev; Z ~ N(0,1)
    tau <- qnorm(1 - bin_prev)
    
    if (p_bin == 1 || ar1_rho_bin == 0) {
      # independent binaries (original behaviour)
      Zb <- matrix(rnorm(n * p_bin), n, p_bin)
    } else {
      # AR(1) correlation on latent scale
      Sigma_bin <- outer(1:p_bin, 1:p_bin,
                         function(i, j) ar1_rho_bin^abs(i - j))
      Zb <- mvrnorm(n, mu = rep(0, p_bin), Sigma = Sigma_bin)
    }
    
    X[, (p_cont + 1):p_total] <- (Zb > tau) * 1L
  }
  
  X <- as.data.frame(X)
  X
}


# Calibrate logistic intercept to hit target treated proportion
calibrate_logit_intercept <- function(eta_without_intercept, target) {
  f <- function(a0) mean(plogis(a0 + eta_without_intercept)) - target
  uniroot(f, interval = c(-10, 10))$root
}

# Compute IPTW variants (glmnet logistic)
compute_iptw_variants <- function(Z, X_df, trunc_lo = 0.01, trunc_hi = 0.99) {
  x_mat <- as.matrix(X_df)
  fit <- cv.glmnet(x = x_mat, y = Z, family = "binomial", alpha = 1, nfolds = 5)  # Reduced from 10 for speed
  ps  <- as.numeric(predict(fit, x_mat, s = "lambda.min", type = "response"))
  ps  <- pmin(pmax(ps, 1e-6), 1 - 1e-6)
  
  pZ  <- mean(Z)
  w_unstab <- ifelse(Z == 1, 1/ps, 1/(1-ps))
  w_stab   <- ifelse(Z == 1, pZ/ps, (1-pZ)/(1-ps))
  
  # Truncation
  lo_u <- quantile(w_unstab, trunc_lo, na.rm = TRUE)
  hi_u <- quantile(w_unstab, trunc_hi, na.rm = TRUE)
  lo_s <- quantile(w_stab,   trunc_lo, na.rm = TRUE)
  hi_s <- quantile(w_stab,   trunc_hi, na.rm = TRUE)
  
  list(
    ps = ps,
    pZ = pZ,
    unstab_untrunc = w_unstab,
    unstab_trunc   = pmin(pmax(w_unstab, lo_u), hi_u),
    stab_untrunc   = w_stab,
    stab_trunc     = pmin(pmax(w_stab,   lo_s), hi_s)
  )
}

# Effective Sample Size (ESS)
ess <- function(w) {
  s1 <- sum(w, na.rm = TRUE)
  s2 <- sum(w^2, na.rm = TRUE)
  if (s2 == 0) return(NA_real_)
  (s1^2)/s2
}

ess_by_arm <- function(w, Z) {
  tibble(
    arm = c("Z=0","Z=1"),
    ESS = c(ess(w[Z==0]), ess(w[Z==1])),
    n   = c(sum(Z==0), sum(Z==1))
  )
}

# Weibull event and censoring time generation under log-linear hazard
# hazard_T(t|X) = h0_T(t; k_event, lambda_event) * exp(eta_T)
# hazard_C(t|X) = h0_C(t; k_cens, lambda_cens)  * exp(eta_C)
# Inverse transform: T = [ -log(U) * lambda * exp(-eta) ]^(1/k)
rweibull_loglin <- function(n, k, lambda, eta) {
  U <- runif(n)
  (-log(U) * lambda * exp(-eta))^(1/k)
}

# Calibrate censoring intercept gamma0 so the censoring proportion ~= target_censoring
calibrate_censor_intercept <- function(etaC_without_intercept, k_censor, lambda_censor,
                                       etaT, k_event, lambda_event,
                                       target_censoring, n_cal = 10000) {  # Reduced from 20000 for speed
  # draw a fresh sample of etas for calibration
  idx <- sample(seq_along(etaC_without_intercept), size = min(n_cal, length(etaC_without_intercept)))
  eC <- etaC_without_intercept[idx]
  eT <- etaT[idx]
  f <- function(g0) {
    C <- rweibull_loglin(length(idx), k_censor, lambda_censor, g0 + eC)
    T <- rweibull_loglin(length(idx), k_event,  lambda_event,  eT)
    mean(T > C) - target_censoring  # P(censored) = P(C < T)
  }
  uniroot(f, interval = c(-8, 8))$root
}

# get marginal HRs and curves by simulating counterfactual event times
get_true_counterfactual <- function(Xdf = Xdf_truth, tt, scenario, k_event, lambda_event) {
  
  message("=== computing marginal beta ===")
  
  # Step 1. Simulate counterfactual event times
  Z0 <- rep(0, nrow(Xdf))
  Z1 <- rep(1, nrow(Xdf))
  
  etas0 <- build_true_etas(Xdf = Xdf, Z = Z0,
                           scenario = scenario,
                           phi_Z = 0,
                           ic_scale = 0)
  etas1 <- build_true_etas(Xdf = Xdf, Z = Z1,
                           scenario = scenario,
                           phi_Z = 0,
                           ic_scale = 0)
  
  T0 <- rweibull_loglin(nrow(Xdf), k_event, lambda_event, etas0$etaT)
  T1 <- rweibull_loglin(nrow(Xdf), k_event, lambda_event, etas1$etaT)
  
  # Step 2. Build a counterfactual trial dataset
  dat_cf <- rbind(
    data.frame(D = 0, time = T0, status = 1),
    data.frame(D = 1, time = T1, status = 1)
  )
  
  # Step 3. Cox model for causal HR
  cox_cf <- coxph(Surv(time, status) ~ D, data = dat_cf)
  beta_true <- coef(cox_cf)[["D"]]
  
  # Step 4. Compute causal RMST difference
  fit_cf <- survfit(Surv(time, status) ~ D, data = dat_cf)
  sfit <- summary(fit_cf, times = tt, extend = TRUE)
  
  surv_real <- sfit$surv
  surv_0 <- sfit$surv[sfit$strata == "D=0"]
  surv_1 <- sfit$surv[sfit$strata == "D=1"]
  dt     <- diff(c(0, tt))
  rmst_0 <- sum(surv_0 * dt)
  rmst_1 <- sum(surv_1 * dt)
  rmst_diff <- rmst_1 - rmst_0
  
  list(beta_true = beta_true,
       rmst_true_diff = rmst_diff,
       surv_real = surv_real,
       surv_0 = surv_0,
       surv_1 = surv_1)
}


# Make a model matrix (no intercept) for a given set of columns
make_mm <- function(df, cols) {
  if (length(cols) == 0) return(matrix(0, nrow(df), 0))
  as.matrix(df[ , cols, drop = FALSE])
}

# Compute weight diagnostics: mean and SD of IPTW, IPCW(t), and combined weights
# by replicate, arm, and time point
compute_weight_diagnostics <- function(dat, dat.long, tt, rep, iptw) {
  
  # Initialize result list
  diag_list <- list()
  
  # Get IPTW weights (stabilized truncated)
  iptw_weights <- iptw$stab_trunc
  
  # For each time point, compute diagnostics
  for (t_idx in seq_along(tt)) {
    t <- tt[t_idx]
    
    # IPTW diagnostics (time-constant, same for all time points)
    for (arm_val in c(0, 1)) {
      arm_idx <- dat$D == arm_val
      if (sum(arm_idx) > 0) {
        iptw_arm <- iptw_weights[arm_idx]
        iptw_mean <- mean(iptw_arm, na.rm = TRUE)
        iptw_sd <- sd(iptw_arm, na.rm = TRUE)
        
        diag_list[[length(diag_list) + 1]] <- data.frame(
          rep = rep,
          arm = arm_val,
          time = t,
          weight_type = "IPTW",
          mean = iptw_mean,
          sd = iptw_sd,
          stringsAsFactors = FALSE
        )
      }
    }
    
    # IPCW diagnostics (time-varying, using risk set)
    # Risk set at time t: rows where Tstart < t <= ti
    risk_set <- dat.long$Tstart < t & t <= dat.long$ti
    
    if (sum(risk_set) > 0) {
      for (arm_val in c(0, 1)) {
        arm_risk_set <- risk_set & dat.long$D == arm_val
        if (sum(arm_risk_set) > 0) {
          ipcw_arm <- dat.long$WStab[arm_risk_set]
          ipcw_mean <- mean(ipcw_arm, na.rm = TRUE)
          ipcw_sd <- sd(ipcw_arm, na.rm = TRUE)
          
          diag_list[[length(diag_list) + 1]] <- data.frame(
            rep = rep,
            arm = arm_val,
            time = t,
            weight_type = "IPCW",
            mean = ipcw_mean,
            sd = ipcw_sd,
            stringsAsFactors = FALSE
          )
        }
      }
    }
    
    # Combined weight diagnostics (time-varying, using risk set)
    if (sum(risk_set) > 0) {
      for (arm_val in c(0, 1)) {
        arm_risk_set <- risk_set & dat.long$D == arm_val
        if (sum(arm_risk_set) > 0) {
          comb_arm <- dat.long$comb_stab_trunc[arm_risk_set]
          comb_mean <- mean(comb_arm, na.rm = TRUE)
          comb_sd <- sd(comb_arm, na.rm = TRUE)
          
          diag_list[[length(diag_list) + 1]] <- data.frame(
            rep = rep,
            arm = arm_val,
            time = t,
            weight_type = "Combined",
            mean = comb_mean,
            sd = comb_sd,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  
  # Combine all diagnostics
  if (length(diag_list) > 0) {
    bind_rows(diag_list)
  } else {
    data.frame(
      rep = integer(0),
      arm = integer(0),
      time = numeric(0),
      weight_type = character(0),
      mean = numeric(0),
      sd = numeric(0)
    )
  }
}
