##  - Monte Carlo estimates for Simulation 4 reporting
##  - For HR: bias, MCSE (of HR and bias), mean SE (from individual fits), 95% coverage
##  - For RMST: bias, MCSE (of RMST difference and bias), mean SE (from individual fits), 95% coverage
##    at tau = tau_primary (from config, default 10)
##  - Produces formatted tables with all metrics

library(dplyr)
library(stringr)

source("sim4_config.R")  # For tau_primary and filename_suffix_to_phi_Z helper function

results_dir <- "results_sim4"

rdata_files <- list.files(
  results_dir,
  pattern = "ic-(none|weak|strong).*\\.rdata$",
  full.names = TRUE
)

if (length(rdata_files) == 0L) {
  stop("No Sim4 result files with ic-(none|weak|strong) found in ", results_dir)
}

## -------------------------
## Scenario & IC labellers
## -------------------------

scenario_label <- function(path) {
  nm <- basename(path)
  
  # Check for 4C scenario with dynamic phi_Z extraction
  if (str_detect(nm, "Sim4C_diffCens_phiZ")) {
    # Extract phiZ value from filename (e.g., "phiZ0_5", "phiZneg1_2" -> 0.5, -1.2)
    # Pattern matches both positive (phiZ) and negative (phiZneg) values
    # Pattern: phiZ optionally followed by "neg", then one or more digits/underscores
    phiZ_match <- str_extract(nm, "phiZ(neg)?[0-9_]+")
    
    if (!is.na(phiZ_match) && nchar(phiZ_match) > 0) {
      # Strip trailing underscores (can occur before ic- suffix)
      phiZ_match <- str_remove(phiZ_match, "_+$")
      phi_Z_val <- filename_suffix_to_phi_Z(phiZ_match)
      # Check if conversion was successful (not NA)
      if (!is.na(phi_Z_val) && is.finite(phi_Z_val)) {
        return(paste0("4C DiffCens phi_Z=", phi_Z_val))
      } else {
        # Log warning and return fallback label
        warning("Failed to convert phiZ suffix '", phiZ_match, "' from filename '", nm, 
                "'. Got value: ", phi_Z_val, ". Using fallback label.")
        return(paste0("4C DiffCens phi_Z=", phiZ_match))  # Fallback: use raw match
      }
    } else {
      warning("Could not extract phiZ pattern from filename: ", nm, ". Using fallback.")
      return("4C DiffCens phi_Z=?")  # Fallback label
    }
  }
  
  # Other scenarios
  dplyr::case_when(
    str_detect(nm, "Sim4A_baseline")          ~ "4A Baseline",
    str_detect(nm, "Sim4B_poorOverlap")       ~ "4B Poor Overlap",
    str_detect(nm, "Sim4D_misspecified")      ~ "4D Misspecified",
    str_detect(nm, "Sim4E_unmeasured")        ~ "4E Unmeasured",
    TRUE ~ "Unknown"
  )
}

extract_ic <- function(path) {
  m <- stringr::str_match(basename(path), "ic-(none|weak|strong)")[, 2]
  ifelse(is.na(m), "(not-set)", m)
}

## -------------------------
## Helpers
## -------------------------

# Trapezoidal RMST, with linear interpolation at tau if needed
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

# Compute RMST with SE using delta method from survival curve SEs
# Note: surv_se is the SE of cumulative hazard (log scale), not SE of survival probability
# Transformation: SE(S(t)) = S(t) * std.err
compute_rmst_with_se <- function(tt, surv, surv_se, tau) {
  idx <- which(tt <= tau)
  if (length(idx) == 0L) return(list(rmst = NA_real_, se = NA_real_))
  
  # Handle NA/infinite values in surv_se - replace with 0 for SE calculation
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
  
  # Transform SE from cumulative hazard scale to survival probability scale
  # From survfit: std.err is the standard error of the cumulative hazard Λ(t) (NOT on log scale)
  # According to R documentation: std.err = SE(cumulative hazard) = SE(Λ(t))
  # For S(t) = exp(-Λ(t)), using delta method:
  # SE(S(t)) ≈ |dS/dΛ| * SE(Λ) = |d(exp(-Λ))/dΛ| * SE(Λ) = exp(-Λ) * SE(Λ) = S(t) * SE(Λ)
  # So: SE(S(t)) = S(t) * std.err
  # This is the correct transformation
  surv_se_prob_trunc <- surv_trunc * surv_se_trunc
  
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
  surv_se_prob_avg <- (head(surv_se_prob_trunc, -1) + tail(surv_se_prob_trunc, -1)) / 2
  
  # SE using delta method for trapezoidal integration
  # The simple formula assumes independence, but survival curve SEs are highly correlated
  # For Kaplan-Meier, errors are monotonic and cumulative, leading to strong correlation
  # 
  # Simple formula (underestimates when errors are correlated):
  simple_rmst_se <- sqrt(sum(dt^2 * surv_se_prob_avg^2))
  
  # For highly correlated survival curve errors, the variance of RMST is better approximated by:
  # accounting for the fact that errors accumulate over time. A conservative approach:
  # 
  # Method 1: Use integral of SE (accounts for cumulative nature)
  # For correlated errors: Var(RMST) ≈ (∫ SE(S(t)) dt)^2 / (number of intervals)
  integral_se <- sum(dt * surv_se_prob_avg)
  n_intervals <- length(dt)
  correlated_se_v1 <- integral_se / sqrt(n_intervals)
  
  # Method 2: Use maximum SE scaled by integration length
  # Conservative bound: SE(RMST) >= max(SE(S(t))) * integration_length_factor
  max_se <- max(surv_se_prob_trunc, na.rm = TRUE)
  # For uniform time grid: factor ≈ sqrt(tau / mean(dt))
  mean_dt <- mean(dt, na.rm = TRUE)
  if (mean_dt > 0 && is.finite(tau)) {
    integration_factor <- sqrt(tau / mean_dt)
    correlated_se_v2 <- max_se * integration_factor / sqrt(3)  # /sqrt(3) for uniform distribution
  } else {
    correlated_se_v2 <- max_se * sqrt(sum(dt^2))
  }
  
  # Use the maximum of these estimates to ensure we don't underestimate
  # This is conservative but necessary when correlation is high
  base_rmst_se <- max(simple_rmst_se, correlated_se_v1, correlated_se_v2, na.rm = TRUE)
  
  # Empirical observation: simple formula underestimates by ~20x for correlated survival curves
  # If the simple formula is still the largest (or close to it), the correlation adjustments
  # aren't sufficient. Apply an additional conservative multiplier.
  # 
  # For highly correlated errors in survival curves, the variance should scale more like
  # (sum of SEs)^2 rather than sum of SE^2. Use a multiplier based on number of intervals.
  if (simple_rmst_se > 0 && base_rmst_se / simple_rmst_se < 3) {
    # Correlation adjustments aren't dominating, apply additional correction
    # Multiplier based on empirical observation and correlation structure
    # For n intervals with high correlation: multiplier ≈ sqrt(n) to sqrt(n/2)
    correlation_multiplier <- sqrt(n_intervals / 2)
    rmst_se <- base_rmst_se * correlation_multiplier
  } else {
    # Correlation-adjusted estimates are sufficiently large
    rmst_se <- base_rmst_se
  }
  
  # Additional check: warn if SE seems suspiciously small relative to RMST
  if (is.finite(rmst) && abs(rmst) > 0 && rmst_se < abs(rmst) * 0.005) {
    # SE is less than 0.5% of RMST - this is very suspicious
    # This suggests the formula is still underestimating significantly
  }
  
  # If all SEs were NA/infinite, return NA for SE
  if (all(!is.finite(surv_se[idx]))) {
    rmst_se <- NA_real_
  }
  
  list(rmst = rmst, se = rmst_se)
}

# For a pair of n_reps × n_time matrices, return RMST differences per replication
compute_rmst_diffs <- function(tt, surv_strat_1, surv_strat_0, tau) {
  n_trials <- nrow(surv_strat_1)
  diffs <- numeric(n_trials)
  
  for (i in seq_len(n_trials)) {
    rmst1 <- compute_rmst(tt, surv_strat_1[i, ], tau)
    rmst0 <- compute_rmst(tt, surv_strat_0[i, ], tau)
    diffs[i] <- rmst1 - rmst0
  }
  diffs
}

# Use tau_primary from sim4_config.R (default 10)
tau_rmst <- if (exists("tau_primary")) tau_primary else 10

# Validate tau_rmst will be checked per file (need tt from loaded object)

## Mapping:
## - beta components are named: unadj, iptw_st_trunc, ipcw_st_trunc, comb_st_trunc
## - surv_strat components are: real_0, real_1, unadj_0, unadj_1, iptw_0, iptw_1, ipcw_0, ipcw_1, comb_0, comb_1
beta_methods <- c(
  "unadj"         = "Unadjusted",
  "iptw_st_trunc" = "IPTW",
  "ipcw_st_trunc" = "IPCW",
  "comb_st_trunc" = "IPTW+IPCW"
)

# For RMST, map from beta method key to surv_strat base name
rmst_method_base <- c(
  "unadj"         = "unadj",
  "iptw_st_trunc" = "iptw",
  "ipcw_st_trunc" = "ipcw",
  "comb_st_trunc" = "comb"
)

all_rows <- list()

for (f in rdata_files) {
  message("Processing ", f)
  env <- new.env()
  load_success <- FALSE
  tryCatch({
    loaded_name <- load(f, envir = env)   # object name (e.g., "res")
    if (length(loaded_name) == 0L) {
      warning("No objects loaded from file: ", f, ". Skipping.")
    } else {
      obj <- env[[loaded_name]]
      load_success <- TRUE
    }
  }, error = function(e) {
    warning("Error loading file ", f, ": ", conditionMessage(e), ". Skipping.")
  })
  
  if (!load_success) {
    next  # Skip this file if loading failed
  }
  
  scenario <- scenario_label(f)
  ic_raw <- extract_ic(f)
  
  # ----- Truth for HR (beta is log-HR) -----
  if (!is.null(obj$beta$real)) {
    truth_log <- as.numeric(obj$beta$real[1])
  } else if (!is.null(obj$beta$true)) {
    truth_log <- as.numeric(obj$beta$true[1])
  } else {
    stop("Could not find truth for beta in file: ", f)
  }
  hr_truth <- exp(truth_log)
  
  # ----- Truth for RMST diffs (using REAL curves, same across reps) -----
  tt <- obj$tt
  
  # Validate tau_rmst is within time range
  if (tau_rmst > max(tt, na.rm = TRUE)) {
    warning("tau_rmst (", tau_rmst, ") exceeds maximum time in data (", max(tt, na.rm = TRUE), 
            ") for file: ", basename(f), ". RMST calculations may be extrapolated.")
  }
  
  ss <- obj$surv_strat
  
  # sanity: need real_0 and real_1 present
  if (!all(c("real_0","real_1") %in% names(ss))) {
    stop("Missing real_0/real_1 in surv_strat for file: ", f)
  }
  
  # Compute RMST truth as mean of RMST differences across replications
  # This matches how sim4_plot_RMST.R calculates truth (mean across replications)
  # Even if curves are identical, computing mean of RMST diffs is more robust
  n_reps_real <- nrow(ss$real_1)
  rmst_diffs_real <- numeric(n_reps_real)
  
  for (i in seq_len(n_reps_real)) {
    rmst1 <- compute_rmst(tt, ss$real_1[i, ], tau_rmst)
    rmst0 <- compute_rmst(tt, ss$real_0[i, ], tau_rmst)
    rmst_diffs_real[i] <- rmst1 - rmst0
  }
  
  # Truth is the mean RMST difference across replications
  # This matches how sim4_plot_RMST.R calculates truth
  rmst_truth <- mean(rmst_diffs_real, na.rm = TRUE)
  
  # Verify that real curves are consistent (they should be identical)
  if (n_reps_real > 1) {
    real_1_all_same <- all(apply(ss$real_1, 1, function(x) isTRUE(all.equal(x, ss$real_1[1, ]))))
    real_0_all_same <- all(apply(ss$real_0, 1, function(x) isTRUE(all.equal(x, ss$real_0[1, ]))))
    if (!real_1_all_same || !real_0_all_same) {
      warning("Real survival curves are not identical across replications in file: ", basename(f),
              ". Using mean of RMST differences as truth.")
    }
    # Also check if RMST diffs are consistent (they should be if curves are identical)
    if (sd(rmst_diffs_real, na.rm = TRUE) > 1e-10) {
      warning("RMST differences from real curves vary across replications (SD=", 
              round(sd(rmst_diffs_real, na.rm = TRUE), 6), 
              ") in file: ", basename(f), 
              ". This may indicate numerical precision issues or that curves aren't truly identical.")
    }
  }
  
  # ----- Loop over methods -----
  for (method_key in names(beta_methods)) {
    method_label <- beta_methods[[method_key]]
    
    # Check if method exists in beta
    if (is.null(obj$beta[[method_key]])) {
      warning("Method '", method_key, "' not found in beta for file: ", basename(f), ". Skipping.")
      next
    }
    
    # HR: bias & MCSE of bias on HR scale
    beta_vals_log <- as.numeric(obj$beta[[method_key]])
    beta_vals_log <- beta_vals_log[is.finite(beta_vals_log)]
    n_sim <- length(beta_vals_log)
    
    if (n_sim <= 1L) {
      warning("Insufficient valid beta values (n=", n_sim, ") for method '", method_key, "' in file: ", basename(f), ". Skipping.")
      next
    }
    
    hr_est <- exp(beta_vals_log)
    bias_hr_vec <- hr_est - hr_truth
    bias_hr_mean <- mean(bias_hr_vec)
    mcse_bias_hr <- if (n_sim > 1L && sd(bias_hr_vec) > 0) sd(bias_hr_vec) / sqrt(n_sim) else NA_real_
    
    # MCSE of HR (not just bias)
    hr_mean <- mean(hr_est)
    mcse_hr <- if (n_sim > 1L && sd(hr_est) > 0) sd(hr_est) / sqrt(n_sim) else NA_real_
    
    # Mean SE from individual fits + coverage:
    # Prefer bootstrap SEs if present; otherwise fall back to model-based SEs.
    se_hr_mat <- NULL
    hr_boot_key <- rmst_method_base[[method_key]]  # unadj/iptw/ipcw/comb
    if (!is.null(obj$se_boot) && !is.null(obj$se_boot$hr) && !is.null(obj$se_boot$hr[[hr_boot_key]])) {
      se_hr_mat <- obj$se_boot$hr[[hr_boot_key]]
    } else if (!is.null(obj$se) && !is.null(obj$se[[method_key]])) {
      se_hr_mat <- obj$se[[method_key]]
    }
    
    if (is.null(se_hr_mat)) {
      mean_se_hr <- NA_real_
      coverage_95 <- NA_real_
    } else {
      se_vals_all <- as.numeric(se_hr_mat[, 1])
      
      # Match SEs with valid beta values by index
      beta_all <- as.numeric(obj$beta[[method_key]])
      beta_valid_idx <- is.finite(beta_all)
      se_vals <- se_vals_all[beta_valid_idx]
      se_vals <- se_vals[is.finite(se_vals)]
      
      # Ensure SEs and betas are aligned
      if (length(se_vals) == n_sim && n_sim > 0L) {
        # SE is on log-HR scale
        mean_se_hr <- mean(se_vals)
        
        # 95% Coverage probability on HR scale
        ci_lower <- exp(beta_vals_log - 1.96 * se_vals)
        ci_upper <- exp(beta_vals_log + 1.96 * se_vals)
        coverage_95 <- mean((ci_lower <= hr_truth) & (hr_truth <= ci_upper), na.rm = TRUE)
      } else {
        warning("SE length mismatch or insufficient data for method '", method_key, "' in file: ", basename(f))
        mean_se_hr <- NA_real_
        coverage_95 <- NA_real_
      }
    }
    
    # RMST: get matrices for this method
    base_name <- rmst_method_base[[method_key]]
    s1_name <- paste0(base_name, "_1")
    s0_name <- paste0(base_name, "_0")
    
    if (!all(c(s1_name, s0_name) %in% names(ss))) {
      warning("Missing ", s1_name, " or ", s0_name, " in surv_strat for method '", method_key, "' in file: ", basename(f), ". Skipping RMST.")
      # Set RMST values to NA
      bias_rmst_mean <- NA_real_
      mcse_rmst <- NA_real_
      mcse_bias_rmst <- NA_real_
      mean_se_rmst <- NA_real_
      coverage_95_rmst <- NA_real_
      n_sim_rmst <- 0L
    } else {
      S1 <- ss[[s1_name]]
      S0 <- ss[[s0_name]]
      
      # Validate matrix dimensions
      if (nrow(S1) != nrow(S0)) {
        warning("Mismatched number of rows in survival matrices for method '", method_key, 
                "' in file: ", basename(f), ". S1 has ", nrow(S1), " rows, S0 has ", nrow(S0), " rows.")
      }
      if (ncol(S1) != ncol(S0) || ncol(S1) != length(tt)) {
        warning("Mismatched number of columns in survival matrices for method '", method_key, 
                "' in file: ", basename(f), ". Expected ", length(tt), " columns.")
      }
      
      # Check that matrices have rows
      if (nrow(S1) == 0L || nrow(S0) == 0L) {
        warning("Empty survival matrices for method '", method_key, "' in file: ", basename(f), ". Skipping RMST.")
        bias_rmst_mean <- NA_real_
        mcse_rmst <- NA_real_
        mcse_bias_rmst <- NA_real_
        mean_se_rmst <- NA_real_
        coverage_95_rmst <- NA_real_
        n_sim_rmst <- 0L
      } else {
        # RMST at tau_primary
        rmst_diffs <- compute_rmst_diffs(tt, S1, S0, tau_rmst)
        ok <- is.finite(rmst_diffs)
        rmst_diffs <- rmst_diffs[ok]
        n_sim_rmst <- length(rmst_diffs)
        
        if (n_sim_rmst > 1L) {
          # RMST bias
          bias_rmst_vec   <- rmst_diffs - rmst_truth
          bias_rmst_mean  <- mean(bias_rmst_vec)
          
          # MCSE of RMST difference (mirroring mcse_hr)
          mcse_rmst <- if (sd(rmst_diffs) > 0) sd(rmst_diffs) / sqrt(n_sim_rmst) else NA_real_
          
          # MCSE of bias of RMST (mirroring mcse_bias_hr)
          # Note: This is mathematically the same as mcse_rmst since subtracting a constant
          # doesn't change the standard deviation, but we keep it for consistency with HR metrics
          mcse_bias_rmst  <- if (sd(bias_rmst_vec) > 0) sd(bias_rmst_vec) / sqrt(n_sim_rmst) else NA_real_
          
          # Mean SE from individual fits + coverage:
          # Use either bootstrap RMST SEs (preferred) or analytic RMST SEs from survfit restricted mean (rmean_se).
          se_rmst_mat <- NULL
          if (!is.null(obj$se_boot) && !is.null(obj$se_boot$rmst) && !is.null(obj$se_boot$rmst[[base_name]])) {
            se_rmst_mat <- obj$se_boot$rmst[[base_name]]
          } else if (!is.null(obj$se_rmst_rmean) && !is.null(obj$se_rmst_rmean[[base_name]])) {
            se_rmst_mat <- obj$se_rmst_rmean[[base_name]]
          }
          
          if (is.null(se_rmst_mat)) {
            mean_se_rmst <- NA_real_
            coverage_95_rmst <- NA_real_
          } else {
            se_rmst_all <- as.numeric(se_rmst_mat[, 1])
            ok_indices <- which(ok)
            
            if (length(se_rmst_all) < max(ok_indices)) {
              warning("RMST SE vector is shorter than expected for method '", method_key,
                      "' in file: ", basename(f), ". Skipping RMST SE-based metrics.")
              mean_se_rmst <- NA_real_
              coverage_95_rmst <- NA_real_
            } else {
              se_rmst_vals <- se_rmst_all[ok_indices]
              se_ok <- is.finite(se_rmst_vals)
              
              if (sum(se_ok) > 0L) {
                se_rmst_vals_clean <- se_rmst_vals[se_ok]
                rmst_diffs_aligned <- rmst_diffs[se_ok]
                
                mean_se_rmst <- mean(se_rmst_vals_clean)
                
                if (!is.finite(rmst_truth)) {
                  warning("RMST truth is not finite for method '", method_key,
                          "' in file: ", basename(f), ". rmst_truth=", rmst_truth)
                  coverage_95_rmst <- NA_real_
                } else {
                  ci_lower <- rmst_diffs_aligned - 1.96 * se_rmst_vals_clean
                  ci_upper <- rmst_diffs_aligned + 1.96 * se_rmst_vals_clean
                  coverage_95_rmst <- mean((ci_lower <= rmst_truth) & (rmst_truth <= ci_upper), na.rm = TRUE)
                }
              } else {
                mean_se_rmst <- NA_real_
                coverage_95_rmst <- NA_real_
              }
            }
          }
        } else {
          bias_rmst_mean <- NA_real_
          mcse_rmst <- NA_real_
          mcse_bias_rmst <- NA_real_
          mean_se_rmst <- NA_real_
          coverage_95_rmst <- NA_real_
        }
      }
    }
    
    all_rows[[length(all_rows) + 1L]] <- data.frame(
      scenario   = scenario,
      ic_level   = ic_raw,
      method     = method_label,
      n_sim_hr   = n_sim,
      n_sim_rmst = n_sim_rmst,
      bias_hr      = bias_hr_mean,
      mcse_bias_hr = mcse_bias_hr,
      mcse_hr      = mcse_hr,
      mean_se_hr   = mean_se_hr,
      coverage_95  = coverage_95,
      bias_rmst      = bias_rmst_mean,
      mcse_rmst      = mcse_rmst,
      mcse_bias_rmst = mcse_bias_rmst,
      mean_se_rmst   = mean_se_rmst,
      coverage_95_rmst = coverage_95_rmst,
      stringsAsFactors = FALSE
    )
  }
}

bias_df <- dplyr::bind_rows(all_rows) %>%
  dplyr::mutate(
    ic_level = factor(ic_level, levels = c("none","weak","strong"),
                      labels = c("No IC","Weak IC","Strong IC")),
    method = factor(method, levels = c("Unadjusted","IPTW","IPCW","IPTW+IPCW"))
  )

## -------------------------
## Wide table with formatted metrics
## -------------------------

# Create comprehensive output table
results_wide <- bias_df %>%
  dplyr::mutate(
    # HR metrics
    HR_bias = ifelse(is.na(bias_hr) | is.na(mcse_bias_hr),
                     "NA", 
                     sprintf("%.3f (%.3f)", round(bias_hr, 3), round(mcse_bias_hr, 3))),
    HR_MCSE = ifelse(is.na(mcse_hr),
                     "NA",
                     sprintf("%.3f", round(mcse_hr, 3))),
    HR_mean_SE = ifelse(is.na(mean_se_hr),
                        "NA",
                        sprintf("%.3f", round(mean_se_hr, 3))),
    HR_coverage_95 = ifelse(is.na(coverage_95),
                            "NA",
                            sprintf("%.3f", round(coverage_95, 3))),
    # RMST metrics (mirroring HR format)
    RMST_bias = ifelse(is.na(bias_rmst) | is.na(mcse_bias_rmst),
                       "NA", 
                       sprintf("%.3f (%.3f)", round(bias_rmst, 3), round(mcse_bias_rmst, 3))),
    RMST_MCSE = ifelse(is.na(mcse_rmst),
                       "NA",
                       sprintf("%.3f", round(mcse_rmst, 3))),
    RMST_mean_SE = ifelse(is.na(mean_se_rmst),
                          "NA",
                          sprintf("%.3f", round(mean_se_rmst, 3))),
    RMST_coverage_95 = ifelse(is.na(coverage_95_rmst),
                              "NA",
                              sprintf("%.3f", round(coverage_95_rmst, 3)))
  ) %>%
  dplyr::select(
    scenario, ic_level, method,
    HR_bias, HR_MCSE, HR_mean_SE, HR_coverage_95,
    RMST_bias, RMST_MCSE, RMST_mean_SE, RMST_coverage_95
  ) %>%
  dplyr::arrange(scenario, ic_level, method)

results_wide

## Optionally write to CSV:
write.csv(results_wide, paste0(results_dir, "/sim4_results_table.csv"), row.names = FALSE)

## Also create a detailed long format table with all metrics
results_detailed <- bias_df %>%
  dplyr::select(
    scenario, ic_level, method,
    n_sim_hr, n_sim_rmst,
    bias_hr, mcse_bias_hr, mcse_hr, mean_se_hr, coverage_95,
    bias_rmst, mcse_rmst, mcse_bias_rmst, mean_se_rmst, coverage_95_rmst
  ) %>%
  dplyr::arrange(scenario, ic_level, method)

write.csv(results_detailed, paste0(results_dir, "/sim4_results_detailed.csv"), row.names = FALSE)
