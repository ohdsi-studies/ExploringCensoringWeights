# ---- Simulation 4 main ----
rm(list = ls())

# Your IPCW machinery (apply the small fixes noted):
source("helper_D_censVars.R")
source("stratified_survival_curves.R")

# Sim 4 config and utils
source("sim4_config.R")
source("sim4_utils_hd.R")

sim4_config_env <- environment()

library(survival)
library(ggplot2)
library(dplyr)
library(tidyr)

# ---- Scenario builders ----

# Build true linear predictors for treatment/event/censoring
# Xdf : full HD covariate matrix (X1, ..., Xp)
# Z   : realized treatment (0/1)
# scenario : "4A","4B","4C","4D","4E_unmeasured" etc. (we only check if == "4C")
# phi_Z   : coefficient on Z in censoring model (differential censoring)
# ic_scale: scales SHARED part between event & censoring (informative censoring)
build_true_etas <- function(Xdf, Z, scenario, phi_Z = 0.0, ic_scale = 1.0) {
  
  n <- nrow(Xdf)
  
  ## -------------------------------
  ## 1. Treatment linear predictor
  ## -------------------------------
  Xtreat <- make_mm(Xdf, paste0("X", idx_treat))
  
  if (ncol(Xtreat) > 0) {
    # one alpha per idx_treat, keep names so we can align with event model
    # alpha_treat <- sample(alpha_vals, length(idx_treat), replace = TRUE)
    # names(alpha_treat) <- as.character(idx_treat)
    alpha_treat <- alpha_treat_map[as.character(idx_treat)]
    etaZ_wo_a0 <- as.numeric(Xtreat %*% alpha_treat)
  } else {
    alpha_treat <- numeric(0)
    etaZ_wo_a0 <- rep(0, n)
  }
  
  ## -------------------------------
  ## 2. Event linear predictor
  ## -------------------------------
  Xevent <- make_mm(Xdf, paste0("X", idx_event))
  
  if (ncol(Xevent) > 0) {
    beta_event <- numeric(length(idx_event))
    names(beta_event) <- as.character(idx_event)
    
    for (j in seq_along(idx_event)) {
      idx_j <- idx_event[j]
      idx_ch <- as.character(idx_j)
      
      if (idx_j %in% idx_confound && idx_ch %in% names(alpha_treat_map)) {
        # shared confounder: event coef is conf_strength × treatment coef
        beta_event[j] <- conf_strength * alpha_treat_map[idx_ch]
      } else {
        # pure outcome predictor: draw from beta_vals
        # beta_event[j] <- sample(beta_vals, 1)
        beta_event[j] <- beta_event_map[idx_ch]
      }
    }
    
    etaT <- beta_Z * Z + as.numeric(Xevent %*% beta_event)
  } else {
    etaT <- beta_Z * Z
  }
  
  ## -------------------------------
  ## 3. Censoring linear predictor
  ## -------------------------------
  shared_idx <- intersect(idx_event, idx_censo)
  uniq_idx   <- setdiff(idx_censo, shared_idx)
  
  Xc_shared <- make_mm(Xdf, paste0("X", shared_idx))
  Xc_uniq   <- make_mm(Xdf, paste0("X", uniq_idx))
  
  if (ncol(Xc_shared) > 0) {
    # p_shared <- sample(phi_vals, ncol(Xc_shared), replace = TRUE)
    p_shared  <- phi_censo_map[as.character(shared_idx)]
    lp_shared <- as.numeric(Xc_shared %*% p_shared)
  } else {
    lp_shared <- 0
  }
  
  if (ncol(Xc_uniq) > 0) {
    # p_uniq <- sample(phi_vals, ncol(Xc_uniq), replace = TRUE)
    p_uniq  <- phi_censo_map[as.character(uniq_idx)]
    lp_uniq <- as.numeric(Xc_uniq %*% p_uniq)
  } else {
    lp_uniq <- 0
  }
  
  # Informative censoring dial: only scales the SHARED part
  etaC_wo_g0 <- ic_scale * lp_shared + lp_uniq +
    if (scenario == "4C") phi_Z * Z else 0.0
  
  list(
    etaZ_wo_a0 = etaZ_wo_a0,  # treatment LP without intercept
    etaT       = etaT,        # event LP (includes beta_Z * Z)
    etaC_wo_g0 = etaC_wo_g0   # censoring LP without intercept
  )
}

# Poor overlap tweak: inflate a few treatment coefficients
inflate_overlap <- function(etaZ_wo_a0, Xdf, inflate_idx = c(1,3,5), inflate_factor = 1.8) {
  # crude inflation by adding a multiple of selected features
  adj <- rowSums(as.matrix(Xdf[ , paste0("X", inflate_idx), drop=FALSE])) * (inflate_factor - 1)
  etaZ_wo_a0 + adj
}

# Misspecification tweak: add nonlinearities/interactions to true DGP only
add_misspec <- function(etaT, etaC_wo_g0, Xdf) {
  # add a couple of unmodeled interactions/quadratics
  if (ncol(Xdf) >= 22) {
    x2 <- Xdf$X2; x10 <- Xdf$X10; x22 <- if ("X22" %in% names(Xdf)) Xdf$X22 else 0
    etaT <- etaT + 0.25*(x2^2) + 0.20*(x10 * x22)
    etaC_wo_g0 <- etaC_wo_g0 + 0.20*(x2 * x10)
  }
  list(etaT = etaT, etaC_wo_g0 = etaC_wo_g0)
}

# ---- Runner for one scenario ----
run_scenario <- function(scen_key, phi_Z = 0.0, ic_key = "weak", checkpoint_file = NULL, checkpoint_freq = 5,
                          use_bootstrap_se = NULL) {
  
  # Get bootstrap setting from config if not provided
  if (is.null(use_bootstrap_se)) {
    use_bootstrap_se <- if (exists("use_bootstrap_se", envir = sim4_config_env)) {
      get("use_bootstrap_se", envir = sim4_config_env)
    } else {
      FALSE
    }
  }
  
  # Source bootstrap file if needed
  if (use_bootstrap_se) {
    source("bootstrap_se_calculations.R")
  }
  
  ic_scale <- sim4_config_env$ic_strength_levels[[ic_key]]
  
  # Checkpoint management: try to load existing checkpoint
  start_rep <- 1
  checkpoint_data <- NULL
  if (!is.null(checkpoint_file) && file.exists(checkpoint_file)) {
    message(sprintf("  Found checkpoint file: %s", checkpoint_file))
    message("  Loading checkpoint and resuming from saved state...")
    checkpoint_data <- tryCatch({
      # Load checkpoint into a new environment
      checkpoint_env <- new.env()
      load(checkpoint_file, envir = checkpoint_env)
      
      # Extract checkpoint_start_rep
      if (exists("checkpoint_start_rep", envir = checkpoint_env)) {
        start_rep <- get("checkpoint_start_rep", envir = checkpoint_env) + 1  # Start from next rep
        message(sprintf("  Resuming from replication %d (checkpoint had %d completed)", 
                       start_rep, get("checkpoint_start_rep", envir = checkpoint_env)))
      } else {
        message("  Warning: Checkpoint file exists but missing start_rep. Starting from beginning.")
        start_rep <- 1
      }
      
      # Return loaded environment as list
      as.list(checkpoint_env)
    }, error = function(e) {
      message(sprintf("  Error loading checkpoint: %s. Starting fresh.", conditionMessage(e)))
      start_rep <- 1
      NULL
    })
  } 
  
  # find true marginal beta
  # Retry logic for exp overflow errors in get_true_counterfactual
  max_truth_retries <- 3
  truth_retry_count <- 0
  truth_success <- FALSE
  
  while (!truth_success && truth_retry_count <= max_truth_retries) {
    if (truth_retry_count > 0) {
      message(sprintf("  Retry attempt %d/%d for computing marginal beta...", 
                      truth_retry_count, max_truth_retries))
    }
    
    tryCatch({
      n_truth <- 1e6  
      
      Xdf_truth <- gen_hd_covariates(
        n        = n_truth,
        p_total  = p_total,
        p_cont   = p_continuous,
        ar1_rho  = ar1_rho,
        p_bin    = p_binary,
        bin_prev = bin_prevalence,
        ar1_rho_bin = ar1_rho_bin
      )
      
      truth <- get_true_counterfactual(
        Xdf          = Xdf_truth,
        tt           = tt,
        scenario     = scen_key,
        k_event      = k_event,
        lambda_event = lambda_event
      )
      
      beta_true_fixed <- truth$beta_true
      truth_success <- TRUE
      
    }, error = function(e) {
      error_msg <- conditionMessage(e)
      if (grepl("exp overflow", error_msg, ignore.case = TRUE)) {
        truth_retry_count <<- truth_retry_count + 1
        if (truth_retry_count > max_truth_retries) {
          stop(sprintf(
            "Failed to compute marginal beta after %d retries due to exp overflow. ",
            max_truth_retries),
            "Error: ", error_msg)
        } else {
          message("  Exp overflow error in get_true_counterfactual. Will retry...")
        }
      } else {
        # For other errors, don't retry, just propagate
        stop(e)
      }
    })
  }
  
  ## ----------------------------------------------------------
  ##  Overall survival (marginal) – what you already had
  ## ----------------------------------------------------------
  surv.real      <- matrix(NA_real_, nrow = n_reps, ncol = length(tt))
  surv.unadj     <- matrix(NA_real_, nrow = n_reps, ncol = length(tt))
  surv.iptw_st_u <- matrix(NA_real_, nrow = n_reps, ncol = length(tt))
  surv.iptw_st_t <- matrix(NA_real_, nrow = n_reps, ncol = length(tt))
  surv.ipcw_st_u <- matrix(NA_real_, nrow = n_reps, ncol = length(tt))
  surv.ipcw_st_t <- matrix(NA_real_, nrow = n_reps, ncol = length(tt))
  surv.combo_st_u<- matrix(NA_real_, nrow = n_reps, ncol = length(tt))
  surv.combo_st_t<- matrix(NA_real_, nrow = n_reps, ncol = length(tt))
  
  # Load from checkpoint if available
  if (!is.null(checkpoint_data)) {
    if ("surv.real" %in% names(checkpoint_data)) surv.real <- checkpoint_data$surv.real
    if ("surv.unadj" %in% names(checkpoint_data)) surv.unadj <- checkpoint_data$surv.unadj
    if ("surv.iptw_st_u" %in% names(checkpoint_data)) surv.iptw_st_u <- checkpoint_data$surv.iptw_st_u
    if ("surv.iptw_st_t" %in% names(checkpoint_data)) surv.iptw_st_t <- checkpoint_data$surv.iptw_st_t
    if ("surv.ipcw_st_u" %in% names(checkpoint_data)) surv.ipcw_st_u <- checkpoint_data$surv.ipcw_st_u
    if ("surv.ipcw_st_t" %in% names(checkpoint_data)) surv.ipcw_st_t <- checkpoint_data$surv.ipcw_st_t
    if ("surv.combo_st_u" %in% names(checkpoint_data)) surv.combo_st_u <- checkpoint_data$surv.combo_st_u
    if ("surv.combo_st_t" %in% names(checkpoint_data)) surv.combo_st_t <- checkpoint_data$surv.combo_st_t
  }
  
  ## ----------------------------------------------------------
  ##  progress bar
  ## ----------------------------------------------------------
  if (start_rep > 1) {
    message(
      sprintf(
        "Running scenario %s (IC = %s, phi_Z = %.2f), %d reps (resuming from rep %d)...",
        scen_key, ic_key, phi_Z, n_reps, start_rep
      )
    )
  } else {
    message(
      sprintf(
        "Running scenario %s (IC = %s, phi_Z = %.2f), %d reps...",
        scen_key, ic_key, phi_Z, n_reps
      )
    )
  }
  pb <- txtProgressBar(min = 0, max = n_reps, style = 3, initial = start_rep - 1)
  on.exit({
    close(pb)
    cat("\n")  # newline after the bar
  }, add = TRUE)
  
  ## ----------------------------------------------------------
  ##  NEW: stratified survival by treatment arm
  ##       (rows = repetition, cols = time points)
  ## ----------------------------------------------------------
  surv_real_strat_1  <- matrix(NA_real_, nrow = n_reps, ncol = length(tt))
  surv_real_strat_0  <- matrix(NA_real_, nrow = n_reps, ncol = length(tt))
  surv_unadj_strat_1 <- matrix(NA_real_, nrow = n_reps, ncol = length(tt))
  surv_unadj_strat_0 <- matrix(NA_real_, nrow = n_reps, ncol = length(tt))
  surv_IPW_strat_1   <- matrix(NA_real_, nrow = n_reps, ncol = length(tt))
  surv_IPW_strat_0   <- matrix(NA_real_, nrow = n_reps, ncol = length(tt))
  surv_IPCW_strat_1  <- matrix(NA_real_, nrow = n_reps, ncol = length(tt))
  surv_IPCW_strat_0  <- matrix(NA_real_, nrow = n_reps, ncol = length(tt))
  surv_comb_strat_1  <- matrix(NA_real_, nrow = n_reps, ncol = length(tt))
  surv_comb_strat_0  <- matrix(NA_real_, nrow = n_reps, ncol = length(tt))
  
  # Load from checkpoint if available
  if (!is.null(checkpoint_data)) {
    if ("surv_real_strat_1" %in% names(checkpoint_data)) surv_real_strat_1 <- checkpoint_data$surv_real_strat_1
    if ("surv_real_strat_0" %in% names(checkpoint_data)) surv_real_strat_0 <- checkpoint_data$surv_real_strat_0
    if ("surv_unadj_strat_1" %in% names(checkpoint_data)) surv_unadj_strat_1 <- checkpoint_data$surv_unadj_strat_1
    if ("surv_unadj_strat_0" %in% names(checkpoint_data)) surv_unadj_strat_0 <- checkpoint_data$surv_unadj_strat_0
    if ("surv_IPW_strat_1" %in% names(checkpoint_data)) surv_IPW_strat_1 <- checkpoint_data$surv_IPW_strat_1
    if ("surv_IPW_strat_0" %in% names(checkpoint_data)) surv_IPW_strat_0 <- checkpoint_data$surv_IPW_strat_0
    if ("surv_IPCW_strat_1" %in% names(checkpoint_data)) surv_IPCW_strat_1 <- checkpoint_data$surv_IPCW_strat_1
    if ("surv_IPCW_strat_0" %in% names(checkpoint_data)) surv_IPCW_strat_0 <- checkpoint_data$surv_IPCW_strat_0
    if ("surv_comb_strat_1" %in% names(checkpoint_data)) surv_comb_strat_1 <- checkpoint_data$surv_comb_strat_1
    if ("surv_comb_strat_0" %in% names(checkpoint_data)) surv_comb_strat_0 <- checkpoint_data$surv_comb_strat_0
  }
  
  ## ----------------------------------------------------------
  ##  NEW: stratified survival SEs by treatment arm
  ##       (rows = repetition, cols = time points)
  ## ----------------------------------------------------------
  surv_se_unadj_strat_1 <- matrix(NA_real_, nrow = n_reps, ncol = length(tt))
  surv_se_unadj_strat_0 <- matrix(NA_real_, nrow = n_reps, ncol = length(tt))
  surv_se_IPW_strat_1   <- matrix(NA_real_, nrow = n_reps, ncol = length(tt))
  surv_se_IPW_strat_0   <- matrix(NA_real_, nrow = n_reps, ncol = length(tt))
  surv_se_IPCW_strat_1  <- matrix(NA_real_, nrow = n_reps, ncol = length(tt))
  surv_se_IPCW_strat_0  <- matrix(NA_real_, nrow = n_reps, ncol = length(tt))
  surv_se_comb_strat_1  <- matrix(NA_real_, nrow = n_reps, ncol = length(tt))
  surv_se_comb_strat_0  <- matrix(NA_real_, nrow = n_reps, ncol = length(tt))
  
  # Analytic RMST SEs from survfit (restricted mean), per replication
  # Stored as SE of RMST difference (Tx - Ctrl), approximated as sqrt(se1^2 + se0^2)
  se_rmst_rmean_unadj <- matrix(NA_real_, nrow = n_reps, ncol = 1)
  se_rmst_rmean_iptw  <- matrix(NA_real_, nrow = n_reps, ncol = 1)
  se_rmst_rmean_ipcw  <- matrix(NA_real_, nrow = n_reps, ncol = 1)
  se_rmst_rmean_comb  <- matrix(NA_real_, nrow = n_reps, ncol = 1)
  
  # Load from checkpoint if available
  if (!is.null(checkpoint_data)) {
    if ("surv_se_unadj_strat_1" %in% names(checkpoint_data)) surv_se_unadj_strat_1 <- checkpoint_data$surv_se_unadj_strat_1
    if ("surv_se_unadj_strat_0" %in% names(checkpoint_data)) surv_se_unadj_strat_0 <- checkpoint_data$surv_se_unadj_strat_0
    if ("surv_se_IPW_strat_1" %in% names(checkpoint_data)) surv_se_IPW_strat_1 <- checkpoint_data$surv_se_IPW_strat_1
    if ("surv_se_IPW_strat_0" %in% names(checkpoint_data)) surv_se_IPW_strat_0 <- checkpoint_data$surv_se_IPW_strat_0
    if ("surv_se_IPCW_strat_1" %in% names(checkpoint_data)) surv_se_IPCW_strat_1 <- checkpoint_data$surv_se_IPCW_strat_1
    if ("surv_se_IPCW_strat_0" %in% names(checkpoint_data)) surv_se_IPCW_strat_0 <- checkpoint_data$surv_se_IPCW_strat_0
    if ("surv_se_comb_strat_1" %in% names(checkpoint_data)) surv_se_comb_strat_1 <- checkpoint_data$surv_se_comb_strat_1
    if ("surv_se_comb_strat_0" %in% names(checkpoint_data)) surv_se_comb_strat_0 <- checkpoint_data$surv_se_comb_strat_0
    if ("se_rmst_rmean_unadj" %in% names(checkpoint_data)) se_rmst_rmean_unadj <- checkpoint_data$se_rmst_rmean_unadj
    if ("se_rmst_rmean_iptw" %in% names(checkpoint_data)) se_rmst_rmean_iptw <- checkpoint_data$se_rmst_rmean_iptw
    if ("se_rmst_rmean_ipcw" %in% names(checkpoint_data)) se_rmst_rmean_ipcw <- checkpoint_data$se_rmst_rmean_ipcw
    if ("se_rmst_rmean_comb" %in% names(checkpoint_data)) se_rmst_rmean_comb <- checkpoint_data$se_rmst_rmean_comb
  }
  
  ## ----------------------------------------------------------
  ##  log-HR (cox) by method – your existing structure
  ## ----------------------------------------------------------
  beta.real         <- matrix(NA_real_, nrow = n_reps, ncol = 1)
  beta.unadj        <- matrix(NA_real_, nrow = n_reps, ncol = 1)
  beta.iptw_st_u    <- matrix(NA_real_, nrow = n_reps, ncol = 1)
  beta.iptw_st_t    <- matrix(NA_real_, nrow = n_reps, ncol = 1)
  beta.ipcw_st_u    <- matrix(NA_real_, nrow = n_reps, ncol = 1)
  beta.ipcw_st_t    <- matrix(NA_real_, nrow = n_reps, ncol = 1)
  beta.combo_st_u   <- matrix(NA_real_, nrow = n_reps, ncol = 1)
  beta.combo_st_t   <- matrix(NA_real_, nrow = n_reps, ncol = 1)
  
  ## ----------------------------------------------------------
  ##  Bootstrap SEs for HR and RMST (if use_bootstrap_se == TRUE)
  ## ----------------------------------------------------------
  if (use_bootstrap_se) {
    n_boot <- if (exists("n_bootstrap_iterations", envir = sim4_config_env)) {
      get("n_bootstrap_iterations", envir = sim4_config_env)
    } else {
      200
    }
    
    # Bootstrap SEs for HR (log-HR)
    se_hr_boot_unadj     <- matrix(NA_real_, nrow = n_reps, ncol = 1)
    se_hr_boot_iptw      <- matrix(NA_real_, nrow = n_reps, ncol = 1)
    se_hr_boot_ipcw      <- matrix(NA_real_, nrow = n_reps, ncol = 1)
    se_hr_boot_comb      <- matrix(NA_real_, nrow = n_reps, ncol = 1)
    
    # Bootstrap SEs for RMST difference
    se_rmst_boot_unadj  <- matrix(NA_real_, nrow = n_reps, ncol = 1)
    se_rmst_boot_iptw   <- matrix(NA_real_, nrow = n_reps, ncol = 1)
    se_rmst_boot_ipcw   <- matrix(NA_real_, nrow = n_reps, ncol = 1)
    se_rmst_boot_comb   <- matrix(NA_real_, nrow = n_reps, ncol = 1)
    
    # Load from checkpoint if available
    if (!is.null(checkpoint_data)) {
      if ("se_hr_boot_unadj" %in% names(checkpoint_data)) se_hr_boot_unadj <- checkpoint_data$se_hr_boot_unadj
      if ("se_hr_boot_iptw" %in% names(checkpoint_data)) se_hr_boot_iptw <- checkpoint_data$se_hr_boot_iptw
      if ("se_hr_boot_ipcw" %in% names(checkpoint_data)) se_hr_boot_ipcw <- checkpoint_data$se_hr_boot_ipcw
      if ("se_hr_boot_comb" %in% names(checkpoint_data)) se_hr_boot_comb <- checkpoint_data$se_hr_boot_comb
      if ("se_rmst_boot_unadj" %in% names(checkpoint_data)) se_rmst_boot_unadj <- checkpoint_data$se_rmst_boot_unadj
      if ("se_rmst_boot_iptw" %in% names(checkpoint_data)) se_rmst_boot_iptw <- checkpoint_data$se_rmst_boot_iptw
      if ("se_rmst_boot_ipcw" %in% names(checkpoint_data)) se_rmst_boot_ipcw <- checkpoint_data$se_rmst_boot_ipcw
      if ("se_rmst_boot_comb" %in% names(checkpoint_data)) se_rmst_boot_comb <- checkpoint_data$se_rmst_boot_comb
    }
  }
  
  # Load from checkpoint if available
  if (!is.null(checkpoint_data)) {
    if ("beta.real" %in% names(checkpoint_data)) beta.real <- checkpoint_data$beta.real
    if ("beta.unadj" %in% names(checkpoint_data)) beta.unadj <- checkpoint_data$beta.unadj
    if ("beta.iptw_st_u" %in% names(checkpoint_data)) beta.iptw_st_u <- checkpoint_data$beta.iptw_st_u
    if ("beta.iptw_st_t" %in% names(checkpoint_data)) beta.iptw_st_t <- checkpoint_data$beta.iptw_st_t
    if ("beta.ipcw_st_u" %in% names(checkpoint_data)) beta.ipcw_st_u <- checkpoint_data$beta.ipcw_st_u
    if ("beta.ipcw_st_t" %in% names(checkpoint_data)) beta.ipcw_st_t <- checkpoint_data$beta.ipcw_st_t
    if ("beta.combo_st_u" %in% names(checkpoint_data)) beta.combo_st_u <- checkpoint_data$beta.combo_st_u
    if ("beta.combo_st_t" %in% names(checkpoint_data)) beta.combo_st_t <- checkpoint_data$beta.combo_st_t
  }
  
  ## ----------------------------------------------------------
  ##  Standard errors (SE) by method – matching beta structure
  ## ----------------------------------------------------------
  se.unadj        <- matrix(NA_real_, nrow = n_reps, ncol = 1)
  se.iptw_st_u    <- matrix(NA_real_, nrow = n_reps, ncol = 1)
  se.iptw_st_t    <- matrix(NA_real_, nrow = n_reps, ncol = 1)
  se.ipcw_st_u    <- matrix(NA_real_, nrow = n_reps, ncol = 1)
  se.ipcw_st_t    <- matrix(NA_real_, nrow = n_reps, ncol = 1)
  se.combo_st_u   <- matrix(NA_real_, nrow = n_reps, ncol = 1)
  se.combo_st_t   <- matrix(NA_real_, nrow = n_reps, ncol = 1)
  
  # Load from checkpoint if available
  if (!is.null(checkpoint_data)) {
    if ("se.unadj" %in% names(checkpoint_data)) se.unadj <- checkpoint_data$se.unadj
    if ("se.iptw_st_u" %in% names(checkpoint_data)) se.iptw_st_u <- checkpoint_data$se.iptw_st_u
    if ("se.iptw_st_t" %in% names(checkpoint_data)) se.iptw_st_t <- checkpoint_data$se.iptw_st_t
    if ("se.ipcw_st_u" %in% names(checkpoint_data)) se.ipcw_st_u <- checkpoint_data$se.ipcw_st_u
    if ("se.ipcw_st_t" %in% names(checkpoint_data)) se.ipcw_st_t <- checkpoint_data$se.ipcw_st_t
    if ("se.combo_st_u" %in% names(checkpoint_data)) se.combo_st_u <- checkpoint_data$se.combo_st_u
    if ("se.combo_st_t" %in% names(checkpoint_data)) se.combo_st_t <- checkpoint_data$se.combo_st_t
  }
  
  ## ----------------------------------------------------------
  ##  Diagnostics
  ## ----------------------------------------------------------
  ess_tab      <- list()
  weight_diag_list <- list()
  cens_prop    <- numeric(n_reps)
  
  # Load from checkpoint if available
  if (!is.null(checkpoint_data)) {
    if ("ess_tab" %in% names(checkpoint_data)) ess_tab <- checkpoint_data$ess_tab
    if ("weight_diag_list" %in% names(checkpoint_data)) weight_diag_list <- checkpoint_data$weight_diag_list
    if ("cens_prop" %in% names(checkpoint_data)) cens_prop <- checkpoint_data$cens_prop
  }
  
  # Maximum number of retries per replication if exp overflow occurs
  max_retries <- 3
  
  # Helper function to save checkpoint
  save_checkpoint <- function(rep_num) {
    if (!is.null(checkpoint_file)) {
      tryCatch({
        # Save all result matrices and current replication number
        checkpoint_start_rep <- rep_num
        
        # Build list of object names to save
        obj_names <- c(
          "checkpoint_start_rep",
          "surv.real", "surv.unadj", "surv.iptw_st_u", "surv.iptw_st_t",
          "surv.ipcw_st_u", "surv.ipcw_st_t", "surv.combo_st_u", "surv.combo_st_t",
          "surv_real_strat_1", "surv_real_strat_0",
          "surv_unadj_strat_1", "surv_unadj_strat_0",
          "surv_IPW_strat_1", "surv_IPW_strat_0",
          "surv_IPCW_strat_1", "surv_IPCW_strat_0",
          "surv_comb_strat_1", "surv_comb_strat_0",
          "surv_se_unadj_strat_1", "surv_se_unadj_strat_0",
          "surv_se_IPW_strat_1", "surv_se_IPW_strat_0",
          "surv_se_IPCW_strat_1", "surv_se_IPCW_strat_0",
          "surv_se_comb_strat_1", "surv_se_comb_strat_0",
          "se_rmst_rmean_unadj", "se_rmst_rmean_iptw", "se_rmst_rmean_ipcw", "se_rmst_rmean_comb",
          "beta.real", "beta.unadj", "beta.iptw_st_u", "beta.iptw_st_t",
          "beta.ipcw_st_u", "beta.ipcw_st_t", "beta.combo_st_u", "beta.combo_st_t",
          "se.unadj", "se.iptw_st_u", "se.iptw_st_t",
          "se.ipcw_st_u", "se.ipcw_st_t", "se.combo_st_u", "se.combo_st_t",
          "ess_tab", "weight_diag_list", "cens_prop"
        )
        
        # Conditionally add bootstrap SE object names if they exist
        if (use_bootstrap_se) {
          obj_names <- c(obj_names, c(
            "se_hr_boot_unadj", "se_hr_boot_iptw", "se_hr_boot_ipcw", "se_hr_boot_comb",
            "se_rmst_boot_unadj", "se_rmst_boot_iptw", "se_rmst_boot_ipcw", "se_rmst_boot_comb"
          ))
        }
        
        # Save using list parameter
        save(list = obj_names, file = checkpoint_file)
      
        message(sprintf("  Checkpoint saved after replication %d", rep_num))
      }, error = function(e) {
        warning(sprintf("  Failed to save checkpoint: %s", conditionMessage(e)))
      })
    }
  }
  
  for (r in start_rep:n_reps) {
    retry_count <- 0
    success <- FALSE
    
    while (!success && retry_count <= max_retries) {
      if (retry_count > 0) {
        message(sprintf("  Retry attempt %d/%d for replication %d...", retry_count, max_retries, r))
      }
      
      tryCatch({
        ## 1) Generate HD covariates
        Xdf <- gen_hd_covariates(
      n        = n_per_rep,
      p_total  = p_total,
      p_cont   = p_continuous,
      ar1_rho  = ar1_rho,
      p_bin    = p_binary,
      bin_prev = bin_prevalence,
      ar1_rho_bin = ar1_rho_bin
    )
    
    X_full <- Xdf
    
    ## 1) Which covariates are "allowed" to be observed?
    if (scen_key == "4E_unmeasured") {
      idx_allowed <- idx_observed          # only observed subset
    } else {
      idx_allowed <- 1:p_total             # all covariates observed
    }
    
    ## Identify informative censoring variables: shared event–censor
    idx_inf_cens <- intersect(idx_event, idx_censo)
    
    ## Final PS covariate set: allowed minus informative censor vars
    idx_ps <- setdiff(idx_allowed, idx_inf_cens)
    idx_censoring = setdiff(idx_allowed, idx_confound) # without confounders
    
    ## Build X_obs that will go into the PS model
    X_obs_ps <- X_full[ , paste0("X", idx_ps), drop = FALSE]
    
    ## X_obs that will go into the *censoring* (IPCW) model
    X_obs_cens <- X_full[ , paste0("X", idx_censoring), drop = FALSE]
    
    
    ## 2) Build true etas (Z unknown yet)
    eta_tmp <- build_true_etas(
      Xdf = Xdf,
      Z = rep(0, n_per_rep),
      scenario = if (startsWith(scen_key, "4C")) "4C" else "4A",
      phi_Z = phi_Z,
      ic_scale = ic_scale
    )
    etaZ_wo_a0 <- eta_tmp$etaZ_wo_a0
    
    ## 3) Calibrate treatment intercept a0
    a0 <- calibrate_logit_intercept(etaZ_wo_a0, target_treated_prop)
    pZ <- plogis(a0 + etaZ_wo_a0)
    Z  <- rbinom(n_per_rep, 1, pZ)
    
    ## 4) Rebuild etas with realized Z
    etas <- build_true_etas(
      Xdf = Xdf,
      Z = Z,
      scenario = if (startsWith(scen_key, "4C")) "4C" else "4A",
      phi_Z = phi_Z,
      ic_scale = ic_scale
    )
    etaT      <- etas$etaT
    etaC_wo_g0 <- etas$etaC_wo_g0
    
    ## Scenario tweaks
    if (scen_key == "4B_poor_overlap") {
      etaZ_wo_a0 <- inflate_overlap(etaZ_wo_a0, Xdf, inflate_idx = c(1, 3, 5), inflate_factor = 2.0)
      a0  <- calibrate_logit_intercept(etaZ_wo_a0, target_treated_prop)
      pZ  <- plogis(a0 + etaZ_wo_a0)
      Z   <- rbinom(n_per_rep, 1, pZ)
      etas <- build_true_etas(Xdf = Xdf, Z = Z, scenario = "4A", phi_Z = 0.0, ic_scale = ic_scale)
      etaT      <- etas$etaT
      etaC_wo_g0<- etas$etaC_wo_g0
    }
    if (scen_key == "4D_misspecified") {
      mm        <- add_misspec(etaT, etaC_wo_g0, Xdf)
      etaT      <- mm$etaT
      etaC_wo_g0<- mm$etaC_wo_g0
    }
    
    ## 5) Calibrate censoring intercept g0
    g0 <- calibrate_censor_intercept(
      etaC_without_intercept = etaC_wo_g0,
      k_censor      = k_censor,
      lambda_censor = lambda_censor,
      etaT          = etaT,
      k_event       = k_event,
      lambda_event  = lambda_event,
      target_censoring = target_censoring
    )
    
    ## 6) Generate event and censoring times
    Ttime <- rweibull_loglin(n_per_rep, k = k_event,  lambda = lambda_event,  eta = etaT)
    Ctime <- rweibull_loglin(n_per_rep, k = k_censor, lambda = lambda_censor, eta = g0 + etaC_wo_g0)
    time  <- pmin(Ttime, Ctime)
    event <- as.numeric(Ttime <= Ctime)
    cens_prop[r] <- mean(event == 0)
    
    ## 7) Build analysis data.frames
    # xi = event; ci = censoring; ti = observed
    dat <- cbind(
      ID = 1:n_per_rep, D = Z, ti = time, di = event
    ) %>% as.data.frame() 
    dat <- cbind(dat, Xdf)
    dat$xi <- Ttime # event time
    dat$ci <- Ctime # censoring time
    dat$one <- 1

    
    beta.real[r, ] <- beta_true_fixed # same through all iterations
    # surv.real[r, ] <- (truth$surv_0 + truth$surv_1) / 2
    
    
    # optional: store stratified curves (by D) 
    surv_real_strat_0[r, ] <- truth$surv_0
    surv_real_strat_1[r, ] <- truth$surv_1
    
    ## ---- unadjusted (overall) ----
    surv.unadj[r, ] <- calc.surv(times = dat$ti, status = dat$di, tt = tt, data = dat)
    beta.unadj[r, ] <- calc.beta(times = dat$ti, status = dat$di, data = dat)
    se.unadj[r, ] <- calc.beta.se(times = dat$ti, status = dat$di, data = dat)
    
    ## ---- NEW: Unadjusted stratified curves (by D) ----
    surv_unadj_strat_df <- calc.surv.unadj.stratified(
      times  = dat$ti,
      status = dat$di,
      tt     = tt,
      data   = dat,
      tau    = tau_primary
    )
    
    surv_unadj_strat_0[r, ] <- surv_unadj_strat_df[surv_unadj_strat_df$D == 0, "survival"]
    surv_unadj_strat_1[r, ] <- surv_unadj_strat_df[surv_unadj_strat_df$D == 1, "survival"]
    surv_se_unadj_strat_0[r, ] <- surv_unadj_strat_df[surv_unadj_strat_df$D == 0, "std_err"]
    surv_se_unadj_strat_1[r, ] <- surv_unadj_strat_df[surv_unadj_strat_df$D == 1, "std_err"]
    se_rmst_rmean_unadj[r, 1] <- {
      se0 <- unique(surv_unadj_strat_df[surv_unadj_strat_df$D == 0, "rmean_se"])
      se1 <- unique(surv_unadj_strat_df[surv_unadj_strat_df$D == 1, "rmean_se"])
      if (length(se0) == 1 && length(se1) == 1 && is.finite(se0) && is.finite(se1)) sqrt(se0^2 + se1^2) else NA_real_
    }
    
    ## ---- Bootstrap SE for unadj (if enabled) ----
    if (use_bootstrap_se) {
      boot_result_unadj <- tryCatch({
        compute_bootstrap_se(
          dat = dat,
          method = "unadj",
          tt = tt,
          tau = tau_primary,
          X_obs_ps = X_obs_ps,
          X_obs_cens = X_obs_cens,
          scen_key = scen_key,
          trunc_lo = trunc_lo,
          trunc_hi = trunc_hi,
          n_boot = n_boot,
          seed = NULL  # Let each replication have different bootstrap samples
        )
      }, error = function(e) {
        warning("Bootstrap SE calculation failed for unadj (rep ", r, "): ", conditionMessage(e))
        list(se_hr = NA_real_, se_rmst = NA_real_)
      })
      
      se_hr_boot_unadj[r, 1] <- boot_result_unadj$se_hr
      se_rmst_boot_unadj[r, 1] <- boot_result_unadj$se_rmst
    }
    
    ## ---- IPTW (glmnet) variants ----
    iptw <- compute_iptw_variants(
      Z     = dat$D,
      X_df  = X_obs_ps, # X_obs,      # <-- only observed covariates enter the PS model
      trunc_lo = trunc_lo,
      trunc_hi = trunc_hi
    )
    
    ## We'll use stabilized *truncated* IPTW for stratified curves (matches your Sim 1–3 plots)
    dat$IPW.Stab <- iptw$stab_trunc
    
    # Overall IPTW curves (untrunc/trunc) – what you had
    surv.iptw_st_u[r, ] <- calc.surv.IPW(times = dat$ti, status = dat$di, tt = tt,
                                         IPW.weights = iptw$stab_untrunc, data = dat)
    beta.iptw_st_u[r, ] <- calc.beta.IPW(times = dat$ti, status = dat$di,
                                         IPW.weights = iptw$stab_untrunc, data = dat)
    se.iptw_st_u[r, ] <- calc.beta.se.IPW(times = dat$ti, status = dat$di,
                                          IPW.weights = iptw$stab_untrunc, data = dat)
    
    surv.iptw_st_t[r, ] <- calc.surv.IPW(times = dat$ti, status = dat$di, tt = tt,
                                         IPW.weights = iptw$stab_trunc, data = dat)
    beta.iptw_st_t[r, ] <- calc.beta.IPW(times = dat$ti, status = dat$di,
                                         IPW.weights = iptw$stab_trunc, data = dat)
    se.iptw_st_t[r, ] <- calc.beta.se.IPW(times = dat$ti, status = dat$di,
                                          IPW.weights = iptw$stab_trunc, data = dat)
    
    ## ---- Bootstrap SE for IPTW (if enabled) ----
    if (use_bootstrap_se) {
      boot_result_iptw <- tryCatch({
        compute_bootstrap_se(
          dat = dat,
          method = "iptw",
          tt = tt,
          tau = tau_primary,
          X_obs_ps = X_obs_ps,
          X_obs_cens = X_obs_cens,
          scen_key = scen_key,
          trunc_lo = trunc_lo,
          trunc_hi = trunc_hi,
          n_boot = n_boot,
          seed = NULL
        )
      }, error = function(e) {
        warning("Bootstrap SE calculation failed for IPTW (rep ", r, "): ", conditionMessage(e))
        list(se_hr = NA_real_, se_rmst = NA_real_)
      })
      
      se_hr_boot_iptw[r, 1] <- boot_result_iptw$se_hr
      se_rmst_boot_iptw[r, 1] <- boot_result_iptw$se_rmst
    }
    
    ## ---- IPTW stratified curves ----
    surv_ipw_strat_df <- calc.surv.IPW.stratified(
      times = dat$ti,
      status = dat$di,
      tt = tt,
      IPW.weights = dat$IPW.Stab,
      data = dat,
      tau  = tau_primary
    )
    surv_IPW_strat_0[r, ] <- surv_ipw_strat_df[surv_ipw_strat_df$D == 0, "survival"]
    surv_IPW_strat_1[r, ] <- surv_ipw_strat_df[surv_ipw_strat_df$D == 1, "survival"]
    surv_se_IPW_strat_0[r, ] <- surv_ipw_strat_df[surv_ipw_strat_df$D == 0, "std_err"]
    surv_se_IPW_strat_1[r, ] <- surv_ipw_strat_df[surv_ipw_strat_df$D == 1, "std_err"]
    se_rmst_rmean_iptw[r, 1] <- {
      se0 <- unique(surv_ipw_strat_df[surv_ipw_strat_df$D == 0, "rmean_se"])
      se1 <- unique(surv_ipw_strat_df[surv_ipw_strat_df$D == 1, "rmean_se"])
      if (length(se0) == 1 && length(se1) == 1 && is.finite(se0) && is.finite(se1)) sqrt(se0^2 + se1^2) else NA_real_
    }
    
    ## ---- IPCW using your machinery ----
    #cut.times.new <- sort(unique(dat$ti)) # slower
    cut.times.new <- sort(unique(round(dat$ti, 2)))  # More aggressive rounding for speed
    dat.long <- transform.data(dat, cut.times = cut.times.new)
    
    # observed covars but without the conf
    # without_conf = setdiff(idx_observed, idx_confound)
    # X_without_conf <- X_obs[ , paste0("X", without_conf), drop = FALSE]
    
    # X_terms_obs <- paste(colnames(X_obs_cens), collapse = " + ") # replace iwth X_obs if you want all X
    # form_C0 <- as.formula(Surv(Tstart, ti, censored, type = "counting") ~ 1)
    # if (startsWith(scen_key, "4C")) {
    #   # differential censoring: include D + observed X
    #   form_CZ <- as.formula(
    #     paste("Surv(Tstart, ti, censored, type='counting') ~", paste(c("D", X_terms_obs), collapse = " + "))
    #   )
    # } else {
    #   # non-differential censoring: only observed X
    #   form_CZ <- as.formula(
    #     paste("Surv(Tstart, ti, censored, type='counting') ~", X_terms_obs)
    #   )
    # }
    
    
    # C0 <- coxph(form_C0, data = dat.long, control = coxph.control(timefix = FALSE))
    # CZ <- coxph(form_CZ, data = dat.long, control = coxph.control(timefix = FALSE)) # NO regularization
    # 
    
    # --- lasso Cox on WIDE for feature selection ---
    dat$censored <- 1 - dat$di
    x_wide <- as.matrix(dat[, c("D", colnames(X_obs_cens)), drop = FALSE])
    y_wide <- with(dat, survival::Surv(ti, censored)) # observed time and censored status
    
    cvfit <- glmnet::cv.glmnet(
      x = x_wide, y = y_wide,
      family = "cox",
      alpha = 1,
      nfolds = 3  # Reduced from 5 for speed
    )
    
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
    
    C0 <- coxph(form_C0, data = dat.long, control = coxph.control(timefix = FALSE))
    CZ <- coxph(form_CZ_sel, data = dat.long, control = coxph.control(timefix = FALSE))

    
    ## Untruncated & truncated IPCW
    dat.long_u <- calc.IPCW_fast(C0, CZ, dat.long, p_trunc = 1.0)
    dat.long_t <- calc.IPCW_fast(C0, CZ, dat.long, p_trunc = trunc_hi)
    
    ## ---- Copy IPTW to long format & build combo weights ----
    # Get raw unstabilized IPTW
    dat.long_u$IPW_unstab_untrunc <- iptw$unstab_untrunc[match(dat.long_u$ID, dat$ID)]
    dat.long_t$IPW_unstab_untrunc <- iptw$unstab_untrunc[match(dat.long_t$ID, dat$ID)]
    
    # Get raw unstabilized IPCW (before truncation)
    # Note: KZti is computed in calc.IPCW_fast, so we can compute raw IPCW
    dat.long_u$IPCW_raw <- 1 / pmax(dat.long_u$KZti, 1e-8)
    dat.long_t$IPCW_raw <- 1 / pmax(dat.long_t$KZti, 1e-8)
    
    # Multiply raw IPTW * raw IPCW
    dat.long_u$comb_raw <- dat.long_u$IPW_unstab_untrunc * dat.long_u$IPCW_raw
    dat.long_t$comb_raw <- dat.long_t$IPW_unstab_untrunc * dat.long_t$IPCW_raw
    
    # Stabilize: multiply by IPTW stabilization factor (pZ or 1-pZ) and IPCW stabilization factor (K0ti)
    pZ <- iptw$pZ
    dat.long_u$IPTW_stab_factor <- ifelse(dat.long_u$D == 1, pZ, 1 - pZ)
    dat.long_t$IPTW_stab_factor <- ifelse(dat.long_t$D == 1, pZ, 1 - pZ)
    
    dat.long_u$comb_stab_untrunc <- dat.long_u$comb_raw * dat.long_u$IPTW_stab_factor * dat.long_u$K0ti
    dat.long_t$comb_stab_untrunc <- dat.long_t$comb_raw * dat.long_t$IPTW_stab_factor * dat.long_t$K0ti
    
    # Truncate
    cap_u <- quantile(dat.long_u$comb_stab_untrunc, trunc_hi, na.rm = TRUE)
    cap_t <- quantile(dat.long_t$comb_stab_untrunc, trunc_hi, na.rm = TRUE)
    dat.long_u$comb_stab_trunc <- pmin(dat.long_u$comb_stab_untrunc, cap_u)
    dat.long_t$comb_stab_trunc <- pmin(dat.long_t$comb_stab_untrunc, cap_t)
    
    # Also keep the old stabilized IPTW for other uses if needed
    dat.long_u$IPW_stab_untrunc <- iptw$stab_untrunc[match(dat.long_u$ID, dat$ID)]
    dat.long_u$IPW_stab_trunc   <- iptw$stab_trunc[match(dat.long_u$ID, dat$ID)]
    dat.long_t$IPW_stab_untrunc <- iptw$stab_untrunc[match(dat.long_t$ID, dat$ID)]
    dat.long_t$IPW_stab_trunc   <- iptw$stab_trunc[match(dat.long_t$ID, dat$ID)]
    
    ## ---- IPCW overall (stabilized, untrunc/trunc) ----
    surv.ipcw_st_u[r, ] <- calc.surv.IPCW(
      Tstart = dat.long_u$Tstart, Tstop = dat.long_u$ti,
      status = dat.long_u$di, tt = tt,
      IPCW.weights = dat.long_u$WStab, data.long = dat.long_u
    )
    beta.ipcw_st_u[r, ] <- calc.beta.IPCW(
      Tstart = dat.long_u$Tstart, Tstop = dat.long_u$ti,
      status = dat.long_u$di, IPCW.weights = dat.long_u$WStab,
      data.long = dat.long_u
    )
    se.ipcw_st_u[r, ] <- calc.beta.se.IPCW(
      Tstart = dat.long_u$Tstart, Tstop = dat.long_u$ti,
      status = dat.long_u$di, IPCW.weights = dat.long_u$WStab,
      data.long = dat.long_u
    )
    
    surv.ipcw_st_t[r, ] <- calc.surv.IPCW(
      Tstart = dat.long_t$Tstart, Tstop = dat.long_t$ti,
      status = dat.long_t$di, tt = tt,
      IPCW.weights = dat.long_t$WStab, data.long = dat.long_t
    )
    beta.ipcw_st_t[r, ] <- calc.beta.IPCW(
      Tstart = dat.long_t$Tstart, Tstop = dat.long_t$ti,
      status = dat.long_t$di, IPCW.weights = dat.long_t$WStab,
      data.long = dat.long_t
    )
    se.ipcw_st_t[r, ] <- calc.beta.se.IPCW(
      Tstart = dat.long_t$Tstart, Tstop = dat.long_t$ti,
      status = dat.long_t$di, IPCW.weights = dat.long_t$WStab,
      data.long = dat.long_t
    )
    
    ## ---- Bootstrap SE for IPCW (if enabled) ----
    if (use_bootstrap_se) {
      boot_result_ipcw <- tryCatch({
        compute_bootstrap_se(
          dat = dat,
          method = "ipcw",
          tt = tt,
          tau = tau_primary,
          X_obs_ps = X_obs_ps,
          X_obs_cens = X_obs_cens,
          scen_key = scen_key,
          trunc_lo = trunc_lo,
          trunc_hi = trunc_hi,
          n_boot = n_boot,
          seed = NULL
        )
      }, error = function(e) {
        warning("Bootstrap SE calculation failed for IPCW (rep ", r, "): ", conditionMessage(e))
        list(se_hr = NA_real_, se_rmst = NA_real_)
      })
      
      se_hr_boot_ipcw[r, 1] <- boot_result_ipcw$se_hr
      se_rmst_boot_ipcw[r, 1] <- boot_result_ipcw$se_rmst
    }
    
    ## ---- COMBO overall (stabilized, untrunc/trunc) ----
    surv.combo_st_u[r, ] <- calc.surv.comb(
      Tstart = dat.long_u$Tstart, Tstop = dat.long_u$ti,
      status = dat.long_u$di, tt = tt,
      comb.weights = dat.long_u$comb_stab_untrunc, data.long = dat.long_u
    )
    beta.combo_st_u[r, ] <- calc.beta.comb(
      Tstart = dat.long_u$Tstart, Tstop = dat.long_u$ti,
      status = dat.long_u$di, comb.weights = dat.long_u$comb_stab_untrunc,
      data.long = dat.long_u
    )
    se.combo_st_u[r, ] <- calc.beta.se.comb(
      Tstart = dat.long_u$Tstart, Tstop = dat.long_u$ti,
      status = dat.long_u$di, comb.weights = dat.long_u$comb_stab_untrunc,
      data.long = dat.long_u
    )
    
    surv.combo_st_t[r, ] <- calc.surv.comb(
      Tstart = dat.long_t$Tstart, Tstop = dat.long_t$ti,
      status = dat.long_t$di, tt = tt,
      comb.weights = dat.long_t$comb_stab_trunc, data.long = dat.long_t
    )
    beta.combo_st_t[r, ] <- calc.beta.comb(
      Tstart = dat.long_t$Tstart, Tstop = dat.long_t$ti,
      status = dat.long_t$di, comb.weights = dat.long_t$comb_stab_trunc,
      data.long = dat.long_t
    )
    se.combo_st_t[r, ] <- calc.beta.se.comb(
      Tstart = dat.long_t$Tstart, Tstop = dat.long_t$ti,
      status = dat.long_t$di, comb.weights = dat.long_t$comb_stab_trunc,
      data.long = dat.long_t
    )
    
    ## ---- IPCW & COMBO stratified (use truncated IPCW / COMBO) ----
    # These functions in stratified_survival_curves.R use WStab / comb.Stab inside, 
    # so we define comb.Stab explicitly for the truncated object:
    dat.long_t$comb.Stab <- dat.long_t$comb_stab_trunc
    
    surv_ipcw_strat_df <- calc.surv.IPCW.stratified(
      Tstart = dat.long_t$Tstart,
      Tstop  = dat.long_t$ti,
      status = dat.long_t$di,
      tt     = tt,
      IPCW.weights = dat.long_t$WStab,
      data.long    = dat.long_t,
      tau          = tau_primary
    )
    surv_IPCW_strat_0[r, ] <- surv_ipcw_strat_df[surv_ipcw_strat_df$D == 0, "survival"]
    surv_IPCW_strat_1[r, ] <- surv_ipcw_strat_df[surv_ipcw_strat_df$D == 1, "survival"]
    surv_se_IPCW_strat_0[r, ] <- surv_ipcw_strat_df[surv_ipcw_strat_df$D == 0, "std_err"]
    surv_se_IPCW_strat_1[r, ] <- surv_ipcw_strat_df[surv_ipcw_strat_df$D == 1, "std_err"]
    se_rmst_rmean_ipcw[r, 1] <- {
      se0 <- unique(surv_ipcw_strat_df[surv_ipcw_strat_df$D == 0, "rmean_se"])
      se1 <- unique(surv_ipcw_strat_df[surv_ipcw_strat_df$D == 1, "rmean_se"])
      if (length(se0) == 1 && length(se1) == 1 && is.finite(se0) && is.finite(se1)) sqrt(se0^2 + se1^2) else NA_real_
    }
    
    surv_comb_strat_df <- calc.surv.comb.stratified(
      Tstart = dat.long_t$Tstart,
      Tstop  = dat.long_t$ti,
      status = dat.long_t$di,
      tt     = tt,
      comb.weights = dat.long_t$comb.Stab,
      data.long    = dat.long_t,
      tau          = tau_primary
    )
    surv_comb_strat_0[r, ] <- surv_comb_strat_df[surv_comb_strat_df$D == 0, "survival"]
    surv_comb_strat_1[r, ] <- surv_comb_strat_df[surv_comb_strat_df$D == 1, "survival"]
    surv_se_comb_strat_0[r, ] <- surv_comb_strat_df[surv_comb_strat_df$D == 0, "std_err"]
    surv_se_comb_strat_1[r, ] <- surv_comb_strat_df[surv_comb_strat_df$D == 1, "std_err"]
    se_rmst_rmean_comb[r, 1] <- {
      se0 <- unique(surv_comb_strat_df[surv_comb_strat_df$D == 0, "rmean_se"])
      se1 <- unique(surv_comb_strat_df[surv_comb_strat_df$D == 1, "rmean_se"])
      if (length(se0) == 1 && length(se1) == 1 && is.finite(se0) && is.finite(se1)) sqrt(se0^2 + se1^2) else NA_real_
    }
    
    ## ---- Bootstrap SE for comb (if enabled) ----
    if (use_bootstrap_se) {
      boot_result_comb <- tryCatch({
        compute_bootstrap_se(
          dat = dat,
          method = "comb",
          tt = tt,
          tau = tau_primary,
          X_obs_ps = X_obs_ps,
          X_obs_cens = X_obs_cens,
          scen_key = scen_key,
          trunc_lo = trunc_lo,
          trunc_hi = trunc_hi,
          n_boot = n_boot,
          seed = NULL
        )
      }, error = function(e) {
        warning("Bootstrap SE calculation failed for comb (rep ", r, "): ", conditionMessage(e))
        list(se_hr = NA_real_, se_rmst = NA_real_)
      })
      
      se_hr_boot_comb[r, 1] <- boot_result_comb$se_hr
      se_rmst_boot_comb[r, 1] <- boot_result_comb$se_rmst
    }
    
    ## ---- ESS diagnostics by arm ----
    ess_tab[[r]] <- bind_rows(
      ess_by_arm(iptw$stab_untrunc, Z) %>% mutate(method = "IPTW_stab_untrunc"),
      ess_by_arm(iptw$stab_trunc,   Z) %>% mutate(method = "IPTW_stab_trunc"),
      ess_by_arm(dat.long_u$WStab[match(dat$ID, dat.long_u$ID)], Z) %>% mutate(method = "IPCW_stab_untrunc"),
      ess_by_arm(dat.long_t$WStab[match(dat$ID, dat.long_t$ID)], Z) %>% mutate(method = "IPCW_stab_trunc"),
      ess_by_arm(tapply(dat.long_u$comb_stab_untrunc, dat.long_u$ID, mean)[dat$ID], Z) %>% mutate(method = "COMB_stab_untrunc"),
      ess_by_arm(tapply(dat.long_t$comb_stab_trunc,   dat.long_t$ID, mean)[dat$ID], Z) %>% mutate(method = "COMB_stab_trunc")
    ) %>% mutate(rep = r)
    
    ## ---- Weight diagnostics ----
    weight_diag_list[[r]] <- compute_weight_diagnostics(
      dat = dat,
      dat.long = dat.long_t,
      tt = tt,
      rep = r,
      iptw = iptw
    )
    
        ## --- update progress bar at the end of the iteration ---
        setTxtProgressBar(pb, r)
        flush.console()   # <-- forces the update to render
        
        # Save checkpoint if enabled and frequency matches
        if (!is.null(checkpoint_file) && (r %% checkpoint_freq == 0 || r == n_reps)) {
          save_checkpoint(r)
        }
        
        # If we get here, the replication was successful
        success <- TRUE
        
      }, error = function(e) {
        # Check if this is the exp overflow error
        error_msg <- conditionMessage(e)
        if (grepl("exp overflow", error_msg, ignore.case = TRUE)) {
          retry_count <<- retry_count + 1
          if (retry_count > max_retries) {
            warning(sprintf(
              "Replication %d failed after %d retries due to exp overflow. ",
              r, max_retries),
              "Skipping this replication. Error: ", error_msg)
            # Mark as success to exit the retry loop, but this rep will have NA values
            success <<- TRUE
          } else {
            # Will retry on next iteration of while loop
            message(sprintf("  Exp overflow error in replication %d. Will retry...", r))
          }
        } else {
          # For other errors, don't retry, just propagate
          stop(e)
        }
      })
    } # end while retry loop
    
  } # end reps
  

  
  ## ----------------------------------------------------------
  ##  Return result object
  ## ----------------------------------------------------------
  list(
    tt = tt,
    surv = list(
      real            = surv.real,
      unadj           = surv.unadj,
      iptw_st_untrunc = surv.iptw_st_u,
      iptw_st_trunc   = surv.iptw_st_t,
      ipcw_st_untrunc = surv.ipcw_st_u,
      ipcw_st_trunc   = surv.ipcw_st_t,
      comb_st_untrunc = surv.combo_st_u,
      comb_st_trunc   = surv.combo_st_t
    ),
    ## NEW: stratified curves by arm for plotting KM / RMST
    surv_strat = list(
      real_1  = surv_real_strat_1,
      real_0  = surv_real_strat_0,
      unadj_1 = surv_unadj_strat_1,
      unadj_0 = surv_unadj_strat_0,
      iptw_1  = surv_IPW_strat_1,
      iptw_0  = surv_IPW_strat_0,
      ipcw_1  = surv_IPCW_strat_1,
      ipcw_0  = surv_IPCW_strat_0,
      comb_1  = surv_comb_strat_1,
      comb_0  = surv_comb_strat_0
    ),
    ## NEW: stratified survival curve SEs by arm for RMST SE calculations
    surv_se_strat = list(
      unadj_1 = surv_se_unadj_strat_1,
      unadj_0 = surv_se_unadj_strat_0,
      iptw_1  = surv_se_IPW_strat_1,
      iptw_0  = surv_se_IPW_strat_0,
      ipcw_1  = surv_se_IPCW_strat_1,
      ipcw_0  = surv_se_IPCW_strat_0,
      comb_1  = surv_se_comb_strat_1,
      comb_0  = surv_se_comb_strat_0
    ),
    beta = list(
      real            = beta.real,
      unadj           = beta.unadj,
      iptw_st_untrunc = beta.iptw_st_u,
      iptw_st_trunc   = beta.iptw_st_t,
      ipcw_st_untrunc = beta.ipcw_st_u,
      ipcw_st_trunc   = beta.ipcw_st_t,
      comb_st_untrunc = beta.combo_st_u,
      comb_st_trunc   = beta.combo_st_t
    ),
    se = list(
      unadj           = se.unadj,
      iptw_st_untrunc = se.iptw_st_u,
      iptw_st_trunc   = se.iptw_st_t,
      ipcw_st_untrunc = se.ipcw_st_u,
      ipcw_st_trunc   = se.ipcw_st_t,
      comb_st_untrunc = se.combo_st_u,
      comb_st_trunc   = se.combo_st_t
    ),
    ## Bootstrap SEs (if use_bootstrap_se == TRUE)
    se_boot = if (use_bootstrap_se) {
      list(
        hr = list(
          unadj = se_hr_boot_unadj,
          iptw  = se_hr_boot_iptw,
          ipcw  = se_hr_boot_ipcw,
          comb  = se_hr_boot_comb
        ),
        rmst = list(
          unadj = se_rmst_boot_unadj,
          iptw  = se_rmst_boot_iptw,
          ipcw  = se_rmst_boot_ipcw,
          comb  = se_rmst_boot_comb
        )
      )
    } else {
      NULL
    },
    ## Analytic RMST SEs derived from survfit restricted mean (always stored; may be NA)
    se_rmst_rmean = list(
      unadj = se_rmst_rmean_unadj,
      iptw  = se_rmst_rmean_iptw,
      ipcw  = se_rmst_rmean_ipcw,
      comb  = se_rmst_rmean_comb
    ),
    ess       = bind_rows(ess_tab),
    weight_diag = bind_rows(weight_diag_list),
    cens_prop = cens_prop
  )

}

