# ---- Simulation 4: Config ----

# Reproducibility
set.seed(20251103)

# Sample sizes and replications
n_per_rep  <- 5000
n_reps     <- 50  # use 500–1000 for final MCSE runs

# Time grid and RMST horizon
tt <- c(seq(0, 9.95, by = 0.05), seq(10, 25, by = 0.1))
tau_primary <- 10 #180   # RMST horizon
# tau_sens    <- 365

# High-dimensional covariates
p_total         <- 100
p_continuous    <- 3          # number continuous
p_binary        <- p_total - p_continuous
ar1_rho         <- 0.4         # correlation among continuous block (AR1)
ar1_rho_bin     <- 0.1         # correlation among binary block
bin_prevalence  <- 0.30        # <-- param so you can change (0–1); Desired marginal prevalence for each binary covariate.

# Treatment prevalence target
target_treated_prop <- 0.45

# Event (Weibull) parameters
k_event     <- 1.0     # shape
lambda_event <- 7      # scale (we’ll adjust overall level via coefficients)

# Censoring target + knobs (Weibull censoring)
target_censoring <- 0.45  # <-- parameterized
k_censor        <- 1.2    # shape (fixed)
lambda_censor   <- 8    # scale (fixed)
# We’ll calibrate the censoring intercept gamma0 to hit target_censoring for each scenario

# Index sets for sparse truth (non-nested vs. Sim 1–3)
# NOTE: these are indices into the HD feature matrix X1..Xp created below
idx_treat <- c(5,7,9,12,15,18,61,64, 85, 86)             # 8 features, 4 confounders
idx_event <- c(1,2,3,4,5,8,10,12,15,17,21,28,33,40,61,66, 80, 81, 82, 83) # 16 features 
idx_censo <- c(1,2,3,4,6,8,10,11,16,22,29,35,62,71, 80, 81, 82, 83)      # 14 features, 10 censoring conf 

# ---- Deterministic coefficient maps (instead of sampling) ----

# helper to recycle a base vector over a set of indices
recycle_map <- function(indices, base_vals) {
  out <- base_vals[(seq_along(indices) - 1L) %% length(base_vals) + 1L]
  names(out) <- as.character(indices)
  out
}

# Coefficient magnitudes (true DGP; sparse)
beta_Z    <- 0.5  # ~ HR is exp(0.50) conditionally
beta_vals <- c(0.35, 0.80, -0.20, 1.2, -0.3, -0.1)  # pool for event covariates
alpha_vals <- c(-0.20, -0.10, -0.40, -1, -0.8)  # pool for treatment covariates
phi_vals   <- c(-0.20, 0.5, 0.30, 0.1, 0.6, 0.8, 1.2) # pool for censor covariates


# this would round robin the coefficients
# alpha_treat_map <- recycle_map(idx_treat, alpha_vals)  # treatment coefs α_j
# beta_event_map  <- recycle_map(idx_event, beta_vals)   # baseline event coefs β_j
# phi_censo_map   <- recycle_map(idx_censo, phi_vals)    # censoring coefs φ_j

# sample the coefficients instead
#set.seed(1001)
alpha_treat_map <- sample(alpha_vals, length(idx_treat), replace = TRUE)
names(alpha_treat_map) <- as.character(idx_treat)
beta_event_map <- sample(beta_vals, length(idx_event), replace = TRUE)
names(beta_event_map) <- as.character(idx_event)
phi_censo_map <- sample(phi_vals, length(idx_censo), replace = TRUE)
names(phi_censo_map) <- as.character(idx_censo)

# Strength levels for informative censoring (scales the SHARED part between event & censor)
ic_strength_levels <- c(none = 0.0, weak = 0.5, strong = 1.0)

# -- deterministic (instead of sampling from coefficients)
set_alpha_from_vector <- function(coefs, idx) {
  stopifnot(length(coefs) == length(idx))
  out <- coefs
  names(out) <- as.character(idx)
  out
}

# alpha_treat_map <- set_alpha_from_vector(
#   coefs = c( -0.2,  0.1,  0.4, -0.2,  0.1,  0.4, -0.2,  0.1 ),
#   idx   = idx_treat
# )
# beta_event_map <- set_alpha_from_vector(
#   coefs = c(0.35,  0.80, -0.20,  2.00,  1.40,
#             -1.50, -0.50,  0.35,  0.80, -0.20,  2.00,
#              1.40, -1.50, -0.50,  0.35,  0.80),
#   idx   = idx_event
# )
# phi_censo_map <- set_alpha_from_vector( # 1, 2, 3, 4, 8, 10 are the shared variables
#   coefs = c(-0.6,  0.5,  0.8,  1.2, -3.0, -0.6,
#             0.5,  0.8,  1.2, -3.0, -0.6,  0.5,  0.8,  1.2 ),
#   idx   = idx_censo
# )



# ---- Confounding structure ----
# True confounders are those that affect BOTH treatment and event
idx_confound <- intersect(idx_treat, idx_event)
idx_inf_cens <- intersect(idx_censo, idx_event)

# Safety fallback in case someone edits idx_treat/idx_event to be disjoint
if (length(idx_confound) == 0L) {
  idx_confound <- 1:5
}

# ---- Confounding strength dial ----
# 0   = no confounding
# 1.0 = baseline level
# >1  = stronger treatment–outcome association via shared covariates
conf_strength <- 1

# ---- Observed vs. unmeasured covariates (for Scenario 5) ----
# For Sim 4 (baseline), we treat ALL p_total covariates as observed.
# In Scenario 5, we will:
#   - let only idx_observed enter the fitted PS / IPCW models
#   - allow some of idx_unmeasured to act as *true* confounders in the DGP.

# By default, say first 80 are observed, last 20 are unmeasured.
# You can easily change p_observed if you want a different split.
p_observed <- 80

idx_observed   <- 1:p_observed
idx_unmeasured <- setdiff(1:p_total, idx_observed)

# Example: in Scenario 4E we will override idx_unmeasured to a subset of true confounders:
# e.g., drop the first 3 shared confounders (if available)
idx_unmeasured_default <- head(idx_confound, 3L)
if (length(idx_unmeasured_default) > 0L) {
  idx_unmeasured <- idx_unmeasured_default
  idx_observed        <- setdiff(1:p_total, idx_unmeasured)
}

# Differential censoring lever (only used in scenario 4C)
# Change these values to set phi_Z for scenario 4C simulations
# File names and plot labels will automatically update to match these values
phi_Z_levels <- c(-0.5, -0.8)  # changing differential informative censoring strength

# Truncation percentiles
trunc_lo <- 0.01
trunc_hi <- 0.99

# Bootstrap SE calculation
use_bootstrap_se <- FALSE  # Toggle to use bootstrap SEs instead of analytical SEs
n_bootstrap_iterations <- 25  # Number of bootstrap iterations

# Output files
out_dir <- "results_sim4"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Optional: suffix helpers for filenames related to informative censoring
ic_suffix <- function(key) paste0("_ic-", key)

# Helper functions for phi_Z value <-> filename conversion
# Converts phi_Z value (e.g., 0.5, -1.2, 3, 1.0) to filename suffix 
# Examples: "phiZ0_5", "phiZneg1_2", "phiZ3", "phiZ1_0"
phi_Z_to_filename_suffix <- function(phi_Z) {
  phi_Z_num <- as.numeric(phi_Z)
  
  # Handle negative values: replace "-" with "neg" prefix
  is_negative <- phi_Z_num < 0
  if (is_negative) {
    phi_Z_str <- as.character(abs(phi_Z_num))
    prefix <- "phiZneg"
  } else {
    phi_Z_str <- as.character(phi_Z_num)
    prefix <- "phiZ"
  }
  
  # Replace decimal point with underscore for filename compatibility
  # For integers (no decimal), keep as-is (e.g., "3" stays "3")
  # For decimals, convert "." to "_" (e.g., "0.5" becomes "0_5")
  if (grepl("\\.", phi_Z_str)) {
    phi_Z_str <- gsub("\\.", "_", phi_Z_str)
    # Clean up trailing zeros: "1_00" -> "1_0", "1_50" stays "1_50"
    phi_Z_str <- gsub("_0+$", "_0", phi_Z_str)
  }
  
  paste0(prefix, phi_Z_str)
}

# Converts filename suffix back to phi_Z value
# Handles both formats: "phiZ3" (integer) and "phiZ3_0" (with decimal)
# Examples: "phiZ0_5" -> 0.5, "phiZneg1_2" -> -1.2, "phiZ3" -> 3, "phiZ1_0" -> 1.0
filename_suffix_to_phi_Z <- function(suffix) {
  if (is.na(suffix) || length(suffix) == 0 || nchar(suffix) == 0) {
    return(NA_real_)
  }
  
  # Remove "phiZ" or "phiZneg" prefix and determine sign
  is_negative <- grepl("^phiZneg", suffix)
  if (is_negative) {
    phi_Z_str <- gsub("^phiZneg", "", suffix)
    sign_mult <- -1
  } else {
    phi_Z_str <- gsub("^phiZ", "", suffix)
    sign_mult <- 1
  }
  
  # Check if we successfully extracted a numeric part
  if (is.na(phi_Z_str) || length(phi_Z_str) == 0 || nchar(phi_Z_str) == 0) {
    return(NA_real_)
  }
  
  # Replace underscore with decimal point
  # If no underscore (integer format like "phiZ3"), it stays as "3"
  # If underscore exists (like "phiZ0_5"), becomes "0.5"
  if (grepl("_", phi_Z_str)) {
    phi_Z_str <- gsub("_", ".", phi_Z_str)
  }
  
  # Convert to numeric and apply sign
  phi_Z_num <- as.numeric(phi_Z_str)
  if (is.na(phi_Z_num)) {
    return(NA_real_)
  }
  
  phi_Z_num * sign_mult
}

# Generate file path for scenario 4C with given phi_Z value
get_4C_file_path <- function(phi_Z) {
  suffix <- phi_Z_to_filename_suffix(phi_Z)
  file.path(out_dir, paste0("Sim4C_diffCens_", suffix, ".rdata"))
}

# Scenario names → file names (4C files are generated dynamically)
sim4_files <- list(
  "4A_baseline"      = file.path(out_dir, "Sim4A_baseline.rdata"),
  "4B_poor_overlap"  = file.path(out_dir, "Sim4B_poorOverlap.rdata"),
  "4D_misspecified"  = file.path(out_dir, "Sim4D_misspecified.rdata"),
  "4E_unmeasured"    = file.path(out_dir, "Sim4E_unmeasured.rdata")
)

combined_fig <- file.path(out_dir, "Sim4_combined_HR_plot.png")
combined_pdf <- file.path(out_dir, "Sim4_combined_HR_plot.pdf")


print_sim4_treatment_coefs <- function(phi_Z = NULL) {
  
  cat("\n==============================\n")
  cat("Simulation 4: True Coefficients\n")
  cat("==============================\n\n")
  
  ## ---------------------------
  ## Treatment model (alpha)
  ## ---------------------------
  cat("α (Treatment assignment coefficients)\n")
  print(alpha_treat_map)
  cat("\n")
  
  ## ---------------------------
  ## Event model (beta)
  ## ---------------------------
  cat("β (Event model coefficients)\n")
  cat(sprintf("  β_Z (treatment effect on hazard): %.3f  [HR = %.3f]\n",
              beta_Z, exp(beta_Z)))
  cat("\nCovariate effects:\n")
  print(beta_event_map)
  cat("\n")
  
  ## ---------------------------
  ## Censoring model (phi)
  ## ---------------------------
  cat("φ (Censoring model coefficients)\n")
  cat("Covariate effects:\n")
  print(phi_censo_map)
  
  if (!is.null(phi_Z)) {
    cat(sprintf("\n  φ_Z (treatment effect on censoring): %.3f\n", phi_Z))
  }
  
  cat("\n==============================\n\n")
  invisible(NULL)
}

