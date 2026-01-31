############################################################
## sim_small_config.R — Minimal Sim 1 & 2 (4 covariates) ##
############################################################

## Reproducibility
set.seed(20251124)

## Sample sizes and replications
n_per_rep <- 2500          # or smaller for quick debugging (e.g., 500)
n_reps    <- 200           # lower for debug, higher for final

# Time grid and RMST horizon
tt <- c(seq(0, 4.95, by = 0.05), seq(5, 25, by = 0.5))
tau_primary <- 4.5


## ---------------------------
## Dimensionality knobs
## ---------------------------
p_total      <- 4
p_continuous <- 2   # (Age, Cens2)
p_binary     <- 2   # (Sex, Cens1)


## ---------------------------
## Treatment prevalence
## ---------------------------
target_treated_prop <- 0.45   # ~45% treated

## ---------------------------
## Event (Weibull) parameters
## ---------------------------
k_event      <- 1.0    # shape
lambda_event <- 7.3    # scale

## ---------------------------
## Censoring target + Weibull knobs
## ---------------------------
target_censoring <- 0.45
k_censor        <- 1.2
lambda_censor   <- 8.0

## ---------------------------
## Coefficient pools (true DGP)
## ---------------------------

## Conditional treatment effect of Z on event time
beta_Z <- 0.7   # HR ~ exp(0.70) conditionally

## Event covariate pool (for age, sex, cens1, cens2)
beta_vals  <- c(0.5, 0.25, -2, -1.5) # prev: 1.5, 0.25, 2, 0.1

## Treatment covariate pool (for age, sex)
alpha_vals <- c(0.75, 0.25)

## Censor covariate pool (for cens1, cens2)
phi_vals   <- c(-0.8, -1.2) # prev: 0.8, 0.2

############################################################
## PH-violation knobs for small sim (Sim3)
############################################################

# Early treatment effect: use your usual log-HR beta_Z
beta_Z_early <- beta_Z          # e.g., log(0.7) or whatever you set beta_Z to

# Late treatment effect: here set to 0 (no effect after t_PH),
# you can change this if you want "reversal" or attenuation.
beta_Z_late  <- 1.2      # (HR = 3.32)        # log-HR after t_PH

# Change-point in time (consistent with tau = 500)
t_PH <- 1.2                      # mid-point of follow-up


## ---------------------------
## Index sets (Sim 1 & 2)
## ---------------------------
## Interpretation:
##   X1 = age
##   X2 = sex
##   X3 = cens1
##   X4 = cens2

idx_treat <- c(1, 2)         # age, sex
idx_event <- c(1, 2, 3, 4)   # age, sex, cens1, cens2
idx_censo <- c(3, 4)         # cens1, cens2

## True confounders: affect BOTH treatment and event
idx_confound <- intersect(idx_treat, idx_event)  # here: c(1, 2)

## Informative censoring variables: shared between event and censor DGP
idx_inf_cens <- intersect(idx_censo, idx_event)  # here: c(3, 4)

## Sanity fallback (shouldn't trigger here)
if (length(idx_confound) == 0L) {
  idx_confound <- 1:2
}

## ---------------------------
## Confounding strength
## ---------------------------
## 0   = no confounding
## 1.0 = baseline
## >1  = stronger link via shared covariates
conf_strength <- 0.35

## ---------------------------
## Observed vs unmeasured (none here)
## ---------------------------
idx_observed   <- 1:p_total  # All 4 are observed
idx_unmeasured <- integer(0)

## ---------------------------
## Differential censoring lever (Sim 1 vs Sim 2)
## ---------------------------
## Sim 1: no differential censoring → phi_Z = 0
## Sim 2: differential censoring → phi_Z > 0

phi_Z_levels <- c(
  sim1 = 0.0,   # censoring DGP does NOT depend on Z
  sim2 = 2 #1.5    # censoring DGP includes Z with coefficient 0.8 (set as you like)
)

## ---------------------------
## Truncation percentiles for weights
## ---------------------------
trunc_lo <- 0.01
trunc_hi <- 0.99

## ---------------------------
## Output directory / filenames
## ---------------------------
out_dir <- "results_sim_small"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

small_files <- list(
  "Sim1_nondiff" = file.path(out_dir, "Sim1_nondiff.rdata"),
  "Sim2_diff"    = file.path(out_dir, "Sim2_diffCens.rdata")
)

