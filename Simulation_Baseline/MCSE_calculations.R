##  - Monte Carlo estimates for simulation study reporting
##  - For HR: bias, MCSE (of HR and bias), mean SE (from individual fits), 95% coverage
##  - For RMST: bias, MCSE at tau = tau_primary (from config, default 4.5)
##  - Produces formatted tables with all metrics

library(dplyr)
library(stringr)

source("sim_small_config.R")

results_dir <- "results_sim_small"

rdata_files <- list.files(
  results_dir,
  pattern = "^res_sim[123]_(none|weak|strong)\\.rdata$",
  full.names = TRUE
)

if (length(rdata_files) == 0L) {
  stop("No sim_small result files found in ", results_dir)
}

## -------------------------
## Scenario & IC labellers
## -------------------------

scenario_label <- function(path) {
  nm <- basename(path)
  
  # Extract sim number from filename (e.g., "res_sim1_none.rdata" -> "sim1")
  scen_match <- str_match(nm, "^res_(sim[123])_")
  
  if (!is.na(scen_match[1, 1]) && !is.na(scen_match[1, 2])) {
    scen_code <- scen_match[1, 2]
    dplyr::case_when(
      scen_code == "sim1" ~ "Sim 1: Nondiff Censoring",
      scen_code == "sim2" ~ "Sim 2: Differential Censoring",
      scen_code == "sim3" ~ "Sim 3: PH Violation",
      TRUE ~ paste0("Sim ", scen_code)
    )
  } else {
    "Unknown"
  }
}

extract_ic <- function(path) {
  # Extract IC level from filename (e.g., "res_sim1_none.rdata" -> "none")
  m <- stringr::str_match(basename(path), "^res_sim[123]_(none|weak|strong)\\.rdata$")
  if (!is.na(m[1, 1]) && !is.na(m[1, 2])) {
    m[1, 2]
  } else {
    "(not-set)"
  }
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
# NOTE: This function is no longer used in the simulation code (run_small.R).
# RMST SEs are now computed directly from survfit's rmean.std.err component.
# This function is kept here for reference/compatibility but is not called.
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

# Use tau_primary from config (default 4.5)
tau_rmst <- if (exists("tau_primary")) tau_primary else 4.5

## Mapping:
## - beta components are named: unadj, iptw, ipcw, combo
## - surv_strat components are: real_0, real_1, unadj_0, unadj_1, iptw_0, iptw_1, ipcw_0, ipcw_1, combo_0, combo_1
beta_methods <- c(
  "unadj" = "Unadjusted",
  "iptw"  = "IPTW",
  "ipcw"  = "IPCW",
  "combo" = "IPTW+IPCW"
)

# For RMST, map from beta method key to surv_strat base name
rmst_method_base <- c(
  "unadj" = "unadj",
  "iptw"  = "iptw",
  "ipcw"  = "ipcw",
  "combo" = "combo"  # Note: surv_strat actually uses "combo" (not "comb")
)

all_rows <- list()

for (f in rdata_files) {
  message("Processing ", f)
  env <- new.env()
  loaded_name <- load(f, envir = env)   # object name (e.g., "res")
  obj <- env[[loaded_name]]
  
  scenario <- scenario_label(f)
  # Try to get IC from filename, fallback to object's ic_label
  ic_raw <- extract_ic(f)
  if (ic_raw == "(not-set)" && !is.null(obj$ic_label)) {
    ic_raw <- obj$ic_label
  }
  
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
  ss <- obj$surv_strat
  
  # sanity: need real_0 and real_1 present
  if (!all(c("real_0","real_1") %in% names(ss))) {
    stop("Missing real_0/real_1 in surv_strat for file: ", f)
  }
  
  # Use first row of REAL stratified curves as truth (all rows are identical copies)
  surv_real_1 <- ss$real_1[1, ]
  surv_real_0 <- ss$real_0[1, ]
  
  rmst_truth <- compute_rmst(tt, surv_real_1, tau_rmst) - compute_rmst(tt, surv_real_0, tau_rmst)
  
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
    
    # Mean SE from individual fits
    # Check if SE is available (may not be in old result files)
    if (is.null(obj$se) || is.null(obj$se[[method_key]])) {
      mean_se_hr <- NA_real_
      coverage_95 <- NA_real_
    } else {
      se_vals_all <- as.numeric(obj$se[[method_key]])
      beta_all <- as.numeric(obj$beta[[method_key]])
      
      # Keep only pairs where both beta and SE are finite (properly aligned)
      both_finite <- is.finite(beta_all) & is.finite(se_vals_all)
      beta_vals_log_aligned <- beta_all[both_finite]
      se_vals <- se_vals_all[both_finite]
      
      n_sim_aligned <- length(beta_vals_log_aligned)
      
      if (n_sim_aligned > 0L) {
        # SE is on log scale (SE of log-HR from Cox model)
        # Standard practice: report SE on log scale for HR
        mean_se_log <- mean(se_vals)
        mean_se_hr <- mean_se_log  # Note: despite name, this is SE on log scale (standard practice)
        
        # 95% Coverage probability
        # All calculations on log scale first, then convert to HR scale
        # beta_vals_log_aligned: log-HR estimates (log scale)
        # se_vals: SE of log-HR (log scale)
        # CI on log scale: beta ± 1.96 * se
        ci_lower_log <- beta_vals_log_aligned - 1.96 * se_vals
        ci_upper_log <- beta_vals_log_aligned + 1.96 * se_vals
        
        # Convert CI bounds to HR scale
        ci_lower <- exp(ci_lower_log)  # HR scale
        ci_upper <- exp(ci_upper_log)  # HR scale
        
        # hr_truth is on HR scale, compare with CI bounds (also on HR scale)
        coverage_95 <- mean((ci_lower <= hr_truth) & (hr_truth <= ci_upper), na.rm = TRUE)
      } else {
        warning("No valid beta-SE pairs for method '", method_key, "' in file: ", basename(f))
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
      mcse_rmst_bias <- NA_real_
      mean_se_rmst <- NA_real_
      coverage_95_rmst <- NA_real_
      n_sim_rmst <- 0L
    } else {
      S1 <- ss[[s1_name]]
      S0 <- ss[[s0_name]]
      
      # Check that matrices have rows
      if (nrow(S1) == 0L || nrow(S0) == 0L) {
        warning("Empty survival matrices for method '", method_key, "' in file: ", basename(f), ". Skipping RMST.")
        bias_rmst_mean <- NA_real_
        mcse_rmst_bias <- NA_real_
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
          bias_rmst_vec   <- rmst_diffs - rmst_truth
          bias_rmst_mean  <- mean(bias_rmst_vec)
          mcse_rmst_bias  <- if (sd(bias_rmst_vec) > 0) sd(bias_rmst_vec) / sqrt(n_sim_rmst) else NA_real_
          
          # Mean SE from individual fits
          if (is.null(obj$se_rmst) || is.null(obj$se_rmst[[method_key]])) {
            mean_se_rmst <- NA_real_
            coverage_95_rmst <- NA_real_
          } else {
            se_rmst_all <- as.numeric(obj$se_rmst[[method_key]])
            # Match SEs with valid RMST differences by index
            rmst_all <- compute_rmst_diffs(tt, S1, S0, tau_rmst)
            rmst_valid_idx <- is.finite(rmst_all)
            se_rmst_vals <- se_rmst_all[rmst_valid_idx]
            se_rmst_vals <- se_rmst_vals[is.finite(se_rmst_vals)]
            
            if (length(se_rmst_vals) == length(rmst_diffs) && length(se_rmst_vals) > 0L) {
              mean_se_rmst <- mean(se_rmst_vals)
              
              # 95% Coverage probability
              ci_lower <- rmst_diffs - 1.96 * se_rmst_vals
              ci_upper <- rmst_diffs + 1.96 * se_rmst_vals
              coverage_95_rmst <- mean((ci_lower <= rmst_truth) & (rmst_truth <= ci_upper), na.rm = TRUE)
            } else {
              mean_se_rmst <- NA_real_
              coverage_95_rmst <- NA_real_
            }
          }
        } else {
          bias_rmst_mean <- NA_real_
          mcse_rmst_bias <- NA_real_
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
      mcse_bias_rmst = mcse_rmst_bias,
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
    # RMST metrics
    RMST_MCSE = ifelse(is.na(mcse_bias_rmst),
                       "NA",
                       sprintf("%.3f", round(mcse_bias_rmst, 3))),
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
    RMST_MCSE, RMST_mean_SE, RMST_coverage_95
  ) %>%
  dplyr::arrange(scenario, ic_level, method)

results_wide

## Optionally write to CSV:
write.csv(results_wide, "sim_small_results_table.csv", row.names = FALSE)

## Also create a detailed long format table with all metrics
results_detailed <- bias_df %>%
  dplyr::select(
    scenario, ic_level, method,
    n_sim_hr, n_sim_rmst,
    bias_hr, mcse_bias_hr, mcse_hr, mean_se_hr, coverage_95,
    bias_rmst, mcse_bias_rmst, mean_se_rmst, coverage_95_rmst
  ) %>%
  dplyr::arrange(scenario, ic_level, method)

write.csv(results_detailed, "sim_small_results_detailed.csv", row.names = FALSE)
