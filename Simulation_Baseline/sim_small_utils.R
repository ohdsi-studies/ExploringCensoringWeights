# ---- Utilities for HD simulation ----
library(MASS)
library(glmnet)
library(Matrix)
library(survival)
library(dplyr)

## ---------------------------
## Covariate generation (4X)
## ---------------------------
## X1 = Age    (continuous, Normal)
## X2 = Sex    (binary, Bernoulli)
## X3 = Cens1  (binary)
## X4 = Cens2  (continuous)

# This overrides the HD generator from sim4_utils_hd.R, but keeps the same signature.
gen_covariates <- function(n,
                           p_total,
                           p_cont) {
  stopifnot(p_total == 4)
  
  age   <- rnorm(n, mean = 0, sd = 1)   # Normalized Age
  sex   <- rbinom(n, size = 1, prob = 0.5)  # Sex
  cens1 <- rbinom(n, size = 1, prob = 0.5)  # Cens1 (binary)
  cens2 <- rnorm(n, mean = 0, sd = 1)       # Cens2 (continuous)
  
  X <- data.frame(
    X1 = age,
    X2 = sex,
    X3 = cens1,
    X4 = cens2
  )
  
  return(X)
}

# Calibrate logistic intercept to hit target treated proportion
calibrate_logit_intercept <- function(eta_without_intercept, target) {
  f <- function(a0) mean(plogis(a0 + eta_without_intercept)) - target
  uniroot(f, interval = c(-10, 10))$root
}

# Compute IPTW variants (glmnet logistic)
compute_iptw_variants <- function(Z, X_df, trunc_lo = 0.01, trunc_hi = 0.99) {
  x_mat <- as.matrix(X_df)
  fit <- cv.glmnet(x = x_mat, y = Z, family = "binomial", alpha = 1, nfolds = 10)
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

# Piecewise exponential with log-linear LP, change at t_star
# Assumes k_event = 1, exponential baseline.
# hazard1(t | X) = exp(eta1) / lambda    for t <= t_star
# hazard2(t | X) = exp(eta2) / lambda    for t >  t_star
rpiecewise_exp_loglin <- function(n, lambda, eta1, eta2, t_star) {
  U <- runif(n)
  rate1 <- exp(eta1) / lambda
  rate2 <- exp(eta2) / lambda
  
  H      <- -log(U)
  H_star <- rate1 * t_star
  
  T <- numeric(n)
  idx_early <- (H <= H_star)
  
  # Event happens before the change time
  T[idx_early] <- H[idx_early] / rate1[idx_early]
  
  # Event happens after the change time
  T[!idx_early] <- t_star + (H[!idx_early] - H_star[!idx_early]) / rate2[!idx_early]
  
  T
}

# Compute RMST with SE using delta method from survival curve SEs
compute_rmst_with_se <- function(tt, surv, surv_se, tau) {
  idx <- which(tt <= tau)
  if (length(idx) == 0L) return(list(rmst = NA_real_, se = NA_real_))
  
  # Handle NA/infinite values in surv_se - replace with 0 for SE calculation
  # (this is conservative - if SE is unknown, we can't compute variance)
  surv_se_clean <- surv_se
  surv_se_clean[!is.finite(surv_se_clean)] <- 0
  
  if (tau %in% tt) {
    time_trunc <- tt[idx]
    surv_trunc <- surv[idx]
    surv_se_trunc <- surv_se_clean[idx]
  } else {
    surv_tau <- approx(tt, surv, xout = tau, method = "linear")$y
    surv_se_tau <- approx(tt, surv_se_clean, xout = tau, method = "linear")$y
    time_trunc <- c(tt[idx], tau)
    surv_trunc <- c(surv[idx], surv_tau)
    surv_se_trunc <- c(surv_se_clean[idx], surv_se_tau)
  }
  
  # RMST (trapezoidal integration)
  dt <- diff(time_trunc)
  surv_avg <- (head(surv_trunc, -1) + tail(surv_trunc, -1)) / 2
  rmst <- sum(dt * surv_avg)
  
  # SE using delta method for trapezoidal integration
  # RMST = sum_i dt_i * (S_i + S_{i+1})/2
  # 
  # For variance calculation, we use the approximation:
  # Var(RMST) ≈ sum_i dt_i^2 * Var((S_i + S_{i+1})/2)
  #
  # Note: This assumes approximate independence between intervals (common approximation)
  # The variance of the average of two correlated estimates is complex, but we approximate:
  # Var((S_i + S_{i+1})/2) ≈ (Var(S_i) + Var(S_{i+1}))/4 when correlation is moderate
  # However, a simpler and more conservative approximation uses the average SE:
  # Var((S_i + S_{i+1})/2) ≈ ((se_i + se_{i+1})/2)^2
  #
  # This gives: Var(RMST) ≈ sum_i dt_i^2 * ((se_i + se_{i+1})/2)^2
  surv_se_avg <- (head(surv_se_trunc, -1) + tail(surv_se_trunc, -1)) / 2
  rmst_se <- sqrt(sum(dt^2 * surv_se_avg^2))
  
  # If all SEs were NA/infinite, return NA for SE
  if (all(!is.finite(surv_se[idx]))) {
    rmst_se <- NA_real_
  }
  
  list(rmst = rmst, se = rmst_se)
}


# Calibrate censoring intercept gamma0 so the censoring proportion ~= target_censoring
calibrate_censor_intercept <- function(etaC_without_intercept, k_censor, lambda_censor,
                                       etaT, k_event, lambda_event,
                                       target_censoring, n_cal = 20000) {
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

############################################################
## Small-case build_true_etas
############################################################
build_true_etas <- function(Xdf, Z, scenario=NULL, phi_Z=0, ic_scale=1) {
  n <- nrow(Xdf)
  
  ## ----- Treatment LP -----
  Xtreat <- make_mm(Xdf, paste0("X", idx_treat))
  alpha_treat <- alpha_vals[seq_along(idx_treat)]
  names(alpha_treat) <- as.character(idx_treat)
  etaZ_wo_a0 <- if (ncol(Xtreat)>0) as.numeric(Xtreat %*% alpha_treat) else rep(0,n)
  
  ## ----- Event LP -----
  Xevent <- make_mm(Xdf, paste0("X", idx_event))
  
  beta_event <- numeric(length(idx_event))
  names(beta_event) <- as.character(idx_event)
  
  conf_pos    <- which(idx_event %in% idx_confound)
  nonconf_pos <- setdiff(seq_along(idx_event), conf_pos)
  
  ## confounders get mapped strength * alpha
  for (k in conf_pos) {
    idx_j <- idx_event[k]
    beta_event[k] <- conf_strength * alpha_treat[as.character(idx_j)]
  }
  
  ## non-confounders get beta_vals (deterministic or sampled)
  if (length(nonconf_pos) > 0) {
    beta_event[nonconf_pos] <- beta_vals[seq_along(nonconf_pos)]
  }
  
  ## Construct baseline (no-treatment) LP
  if (ncol(Xevent) > 0) {
    eta_base <- as.numeric(Xevent %*% beta_event)
  } else {
    eta_base <- rep(0, nrow(Xdf))
  }
  
  ## Full event LP = baseline + treatment effect
  etaT <- beta_Z * Z + eta_base
  
  ## ----- Censoring LP -----
  shared_idx <- intersect(idx_event, idx_censo)
  Xc_shared  <- make_mm(Xdf, paste0("X", shared_idx))
  
  phi_shared <- phi_vals[seq_along(shared_idx)]
  names(phi_shared) <- as.character(shared_idx)
  
  lp_shared <- if (ncol(Xc_shared)>0) as.numeric(Xc_shared %*% phi_shared) else 0
  
  ## uniq (none here)
  lp_uniq <- 0
  
  etaC_wo_g0 <- ic_scale * lp_shared + lp_uniq + phi_Z * Z
  
  list(
    etaZ_wo_a0 = etaZ_wo_a0,   # treatment LP (no intercept)
    eta_base   = eta_base,     # baseline event LP (no treatment)
    etaT       = etaT,         # full event LP (baseline + beta_Z * Z)
    etaC_wo_g0 = etaC_wo_g0    # censoring LP (no intercept)
  )
  
}


# Make a model matrix (no intercept) for a given set of columns
make_mm <- function(df, cols) {
  cols <- cols[cols %in% colnames(df)]
  if (length(cols) == 0) return(matrix(0, nrow(df), 0))
  as.matrix(df[, cols, drop = FALSE])
}
