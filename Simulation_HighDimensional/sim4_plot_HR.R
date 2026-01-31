# sim4_plot.R
# Hazard Ratio plots for Simulation 4

library(dplyr)
library(ggplot2)
library(stringr)
library(purrr)
library(tidyr)

source("sim4_config.R")  # For filename_suffix_to_phi_Z helper function

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
      return("4C DiffCens phi_Z=?")  # Fallback label
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
# HR helpers
# -------------------------

# Given an n_reps x 1 matrix/vector of log-HRs, return HR summary
compute_hr_stats <- function(beta_mat) {
  log_vals <- as.numeric(beta_mat)
  log_vals <- log_vals[is.finite(log_vals)]
  if (length(log_vals) == 0) {
    return(list(
      hr_mean  = NA_real_,
      hr_selog = NA_real_,
      hr_lcl   = NA_real_,
      hr_ucl   = NA_real_
    ))
  }
  
  m_log  <- mean(log_vals)
  se_log <- sd(log_vals)          # Monte Carlo SD on log-scale
  hr_mean <- exp(m_log)
  hr_lcl  <- exp(m_log - 1.96 * se_log)
  hr_ucl  <- exp(m_log + 1.96 * se_log)
  
  list(
    hr_mean  = hr_mean,
    hr_selog = se_log,
    hr_lcl   = hr_lcl,
    hr_ucl   = hr_ucl
  )
}

all_results <- list()

for (f in rdata_files) {
  message("Processing ", f)
  
  loaded_name <- load(f)   # e.g. "obj"
  obj <- get(loaded_name)
  
  # sanity check: beta list elements we need
  needed_beta <- c(
    "real", "unadj",
    "iptw_st_trunc",
    "ipcw_st_trunc",
    "comb_st_trunc"
  )
  if (!all(needed_beta %in% names(obj$beta))) {
    stop("Missing required beta components in file: ", f)
  }
  
  scenario <- scenario_label(f)
  ic_level <- extract_ic(f)
  
  # compute HR stats for each method
  res_real <- compute_hr_stats(obj$beta$real)
  res_unad <- compute_hr_stats(obj$beta$unadj)
  res_iptw <- compute_hr_stats(obj$beta$iptw_st_trunc)
  res_ipcw <- compute_hr_stats(obj$beta$ipcw_st_trunc)
  res_comb <- compute_hr_stats(obj$beta$comb_st_trunc)
  
  df <- data.frame(
    method = c("Real", "Unadjusted", "IPTW", "IPCW", "IPTW+IPCW"),
    hr     = c(
      res_real$hr_mean,
      res_unad$hr_mean,
      res_iptw$hr_mean,
      res_ipcw$hr_mean,
      res_comb$hr_mean
    ),
    se_log = c(
      res_real$hr_selog,
      res_unad$hr_selog,
      res_iptw$hr_selog,
      res_ipcw$hr_selog,
      res_comb$hr_selog
    ),
    lower_ci = c(
      res_real$hr_lcl,
      res_unad$hr_lcl,
      res_iptw$hr_lcl,
      res_ipcw$hr_lcl,
      res_comb$hr_lcl
    ),
    upper_ci = c(
      res_real$hr_ucl,
      res_unad$hr_ucl,
      res_iptw$hr_ucl,
      res_ipcw$hr_ucl,
      res_comb$hr_ucl
    ),
    scenario = scenario,
    ic_level = ic_level,
    stringsAsFactors = FALSE
  )
  
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

# Extract "truth" (Real HR) per scenario × IC combo
truth_df <- results_df %>%
  filter(method == "Real") %>%
  dplyr::select(scenario, ic_level, hr_true = hr) %>%
  distinct()

# -------------------------
# Plot
# -------------------------

hr_pl <- ggplot(results_df %>% filter(method != "Real"),
                aes(x = method, y = hr)) +
  #geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.6) +
  geom_hline(data = truth_df,
             aes(yintercept = hr_true),
             color = "red", linewidth = 0.6, alpha = 0.5) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci),
                width = 0.2) +
  facet_grid(ic_level ~ scenario) +
  ylab("Hazard Ratio (Tx vs Ctrl)") +
  xlab("Method") +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x    = element_text(angle = 45, hjust = 1),
    strip.text     = element_text(size = 8),
    legend.position = "none"
  )

# Save combined plot (matches sim4_config paths/names)
ggsave(file.path(results_dir, "Sim4_combined_HR_plot.png"), hr_pl,
       width = 10, height = 8, dpi = 300)
ggsave(file.path(results_dir, "Sim4_combined_HR_plot.pdf"), hr_pl,
       width = 10, height = 8)

