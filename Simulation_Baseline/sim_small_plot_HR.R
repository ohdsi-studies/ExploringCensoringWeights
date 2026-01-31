# sim_small_plot_HR.R
# Hazard Ratio plots for sim_small (Sim1, Sim2, Sim3)

library(dplyr)
library(ggplot2)
library(stringr)
library(purrr)
library(tidyr)

results_dir <- "results_sim_small"

# Grab all .rdata files in that folder
rdata_files <- list.files(results_dir, pattern = "\\.rdata$", full.names = TRUE)

# ------------------------------------------------
# Scenario & IC label helpers for sim_small files
# Filenames look like: res_sim1_none.rdata, res_sim2_weak.rdata, ...
# ------------------------------------------------

scenario_label <- function(path) {
  nm <- basename(path)
  scen <- str_match(nm, "res_(sim[0-9]+)_")[, 2]
  
  case_when(
    scen == "sim1" ~ "Sim1: Nondiff",
    scen == "sim2" ~ "Sim2: Diff censoring",
    scen == "sim3" ~ "Sim3: Non-PH",
    TRUE           ~ "Unknown"
  )
}

extract_ic <- function(path) {
  # matches "_none.rdata", "_weak.rdata", "_strong.rdata"
  m <- str_match(basename(path), "_(none|weak|strong)\\.rdata$")[, 2]
  ifelse(is.na(m), "(not-set)", m)
}

# ------------------------------------------------
# HR helpers
# ------------------------------------------------

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
  
  loaded_name <- load(f)   # should be "res"
  obj <- get(loaded_name)
  
  # For sim_small, beta list elements are:
  #   "true", "unadj", "iptw", "ipcw", "combo"
  needed_beta <- c("true", "unadj", "iptw", "ipcw", "combo")
  if (!all(needed_beta %in% names(obj$beta))) {
    stop("Missing required beta components in file: ", f)
  }
  
  scenario <- scenario_label(f)
  ic_level <- extract_ic(f)
  
  # compute HR stats for each method
  res_true <- compute_hr_stats(obj$beta$true)
  res_unad <- compute_hr_stats(obj$beta$unadj)
  res_iptw <- compute_hr_stats(obj$beta$iptw)
  res_ipcw <- compute_hr_stats(obj$beta$ipcw)
  res_comb <- compute_hr_stats(obj$beta$combo)
  
  df <- data.frame(
    method = c("True", "Unadjusted", "IPTW", "IPCW", "IPTW+IPCW"),
    hr     = c(
      res_true$hr_mean,
      res_unad$hr_mean,
      res_iptw$hr_mean,
      res_ipcw$hr_mean,
      res_comb$hr_mean
    ),
    se_log = c(
      res_true$hr_selog,
      res_unad$hr_selog,
      res_iptw$hr_selog,
      res_ipcw$hr_selog,
      res_comb$hr_selog
    ),
    lower_ci = c(
      res_true$hr_lcl,
      res_unad$hr_lcl,
      res_iptw$hr_lcl,
      res_ipcw$hr_lcl,
      res_comb$hr_lcl
    ),
    upper_ci = c(
      res_true$hr_ucl,
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
  levels = c("True", "Unadjusted", "IPTW", "IPCW", "IPTW+IPCW")
)

results_df$ic_level <- factor(
  results_df$ic_level,
  levels = c("none", "weak", "strong"),
  labels = c("No IC", "Weak IC", "Strong IC")
)

# Extract "truth" (True HR) per scenario × IC combo
truth_df <- results_df %>%
  filter(method == "True") %>%
  dplyr::select(scenario, ic_level, hr_true = hr) %>%
  distinct()

# ------------------------------------------------
# Plot
# ------------------------------------------------

hr_pl <- ggplot(results_df %>% filter(method != "True"),
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
    strip.text     = element_text(size = 14),
    legend.position = "none"
  )

# Save combined plot
ggsave(file.path(results_dir, "sim_small_combined_HR_plot.png"), hr_pl,
       width = 10, height = 8, dpi = 300)
ggsave(file.path(results_dir, "sim_small_combined_HR_plot.pdf"), hr_pl,
       width = 10, height = 8)
