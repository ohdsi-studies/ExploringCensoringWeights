# sim4_plot_RMST.R

# detach("package:MASS", unload = TRUE)

library(dplyr)
library(ggplot2)
library(stringr)
library(purrr)
library(tidyr)

source("sim4_config.R")

results_dir <- "results_sim4"
rdata_files <- list.files(results_dir, pattern = "\\.rdata$", full.names = TRUE)
rdata_files <- rdata_files[
  str_detect(basename(rdata_files), "ic-(none|weak|strong)") &
  !str_detect(basename(rdata_files), "_checkpoint")
]

# -------------------------
# Scenario & IC labellers
# -------------------------

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
        return(paste0("4C DiffCens phiZ=", phi_Z_val))
      } else {
        # Log warning and return fallback label
        warning("Failed to convert phiZ suffix '", phiZ_match, "' from filename '", nm, 
                "'. Got value: ", phi_Z_val, ". Using fallback label.")
        return(paste0("4C DiffCens phiZ=", phiZ_match))  # Fallback: use raw match
      }
    } else {
      warning("Could not extract phiZ pattern from filename: ", nm, ". Using fallback.")
      return("4C DiffCens \u03D5Z=?")  # Fallback label
    }
  }
  
  # Other scenarios
  case_when(
    str_detect(nm, "Sim4A_baseline")          ~ "4A Baseline",
    str_detect(nm, "Sim4B_poorOverlap")       ~ "4B Poor Overlap",
    str_detect(nm, "Sim4D_misspecified")      ~ "4D Misspecified",
    str_detect(nm, "Sim4E_unmeasured")        ~ "4E Unmeasured",
    TRUE ~ "Unknown"
  )
}

extract_ic <- function(path) {
  m <- str_match(basename(path), "ic-(none|weak|strong)")[, 2]
  ifelse(is.na(m), "(not-set)", m)
}

# -------------------------
# RMST helpers
# -------------------------

compute_rmst <- function(tt, surv, tau) {
  idx <- which(tt <= tau)
  if (length(idx) == 0) return(NA_real_)
  
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

compute_rmst_diff_se <- function(tt, surv_strat_1, surv_strat_0, tau) {
  n_trials <- nrow(surv_strat_1)
  diffs <- numeric(n_trials)
  
  for (i in 1:n_trials) {
    rmst1 <- compute_rmst(tt, surv_strat_1[i, ], tau)
    rmst0 <- compute_rmst(tt, surv_strat_0[i, ], tau)
    diffs[i] <- rmst1 - rmst0
  }
  
  list(
    rmst_diff_mean = mean(diffs, na.rm = TRUE),
    rmst_diff_se   = sd(diffs,   na.rm = TRUE)
  )
}

selected_tau <- tau_primary

all_results <- list()

for (f in rdata_files) {
  message("Processing ", f)
  
  # Each .rdata contains a single object, saved as `obj`
  loaded_name <- load(f)          # e.g. "obj"
  obj <- get(loaded_name)         # retrieve it
  
  # pull time grid and stratified curves out of obj
  tt  <- obj$tt
  ss  <- obj$surv_strat
  
  # sanity check: required elements
  needed <- c("real_1","real_0",
              "unadj_1","unadj_0",
              "iptw_1","iptw_0",
              "ipcw_1","ipcw_0",
              "comb_1","comb_0")
  if (!all(needed %in% names(ss))) {
    stop("Missing stratified survival matrices in file: ", f)
  }
  
  scenario  <- scenario_label(f)
  ic_level  <- extract_ic(f)
  
  # Compute RMST differences for each method
  res_real <- compute_rmst_diff_se(tt, ss$real_1,  ss$real_0,  tau = selected_tau)
  res_unad <- compute_rmst_diff_se(tt, ss$unadj_1, ss$unadj_0, tau = selected_tau)
  res_iptw <- compute_rmst_diff_se(tt, ss$iptw_1,  ss$iptw_0,  tau = selected_tau)
  res_ipcw <- compute_rmst_diff_se(tt, ss$ipcw_1,  ss$ipcw_0,  tau = selected_tau)
  res_comb <- compute_rmst_diff_se(tt, ss$comb_1,  ss$comb_0,  tau = selected_tau)
  
  df <- data.frame(
    method = c("Real", "Unadjusted", "IPTW", "IPCW", "IPTW+IPCW"),
    rmst_diff = c(
      res_real$rmst_diff_mean,
      res_unad$rmst_diff_mean,
      res_iptw$rmst_diff_mean,
      res_ipcw$rmst_diff_mean,
      res_comb$rmst_diff_mean
    ),
    se = c(
      res_real$rmst_diff_se,
      res_unad$rmst_diff_se,
      res_iptw$rmst_diff_se,
      res_ipcw$rmst_diff_se,
      res_comb$rmst_diff_se
    ),
    scenario = scenario,
    ic_level = ic_level,
    stringsAsFactors = FALSE
  )
  
  df$upper_ci <- df$rmst_diff + 1.96 * df$se
  df$lower_ci <- df$rmst_diff - 1.96 * df$se
  
  all_results[[f]] <- df
}

results_df <- bind_rows(all_results)

# factor ordering
results_df$method <- factor(
  results_df$method,
  levels = c("Real", "Unadjusted", "IPTW", "IPCW", "IPTW+IPCW")
)

results_df$ic_level <- factor(
  results_df$ic_level,
  levels = c("none", "weak", "strong"),
  labels = c("No IC", "Weak IC", "Strong IC")
)

# Extract "truth" (Real) per scenario × IC combo
truth_df <- results_df %>%
  filter(method == "Real") %>%
  dplyr::select(scenario, ic_level, real_mean = rmst_diff) %>%
  distinct()

rmst_pl <- ggplot(results_df %>% filter(method != "Real"),
                  aes(x = method, y = rmst_diff)) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
  geom_hline(data = truth_df,
             aes(yintercept = real_mean),
             color = "red", linewidth = 0.6, alpha = 0.5) +
  facet_grid(ic_level ~ scenario) +
  ylab(paste0("RMST Difference (Tx - Ctrl), Tau = ", selected_tau)) +
  xlab("Method") +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    strip.text   = element_text(size = 10),
    legend.position = "none"
  )

ggsave(file.path(results_dir, "Sim4_RMSTdiff.png"), rmst_pl,
       width = 10, height = 8, dpi = 300)
ggsave(file.path(results_dir, "Sim4_RMSTdiff.pdf"), rmst_pl,
       width = 10, height = 8)
