# sim_small_plot_RMST.R
# RMST (Tx - Ctrl) plots for sim_small (Sim1, Sim2, Sim3)

library(dplyr)
library(ggplot2)
library(stringr)
library(purrr)
library(tidyr)

source("sim_small_config.R")

results_dir <- "results_sim_small"
rdata_files <- list.files(results_dir, pattern = "\\.rdata$", full.names = TRUE)

# -------------------------
# Scenario & IC labellers (sim_small)
# -------------------------

scenario_label <- function(path) {
  nm <- basename(path)
  scen <- str_match(nm, "res_(sim[0-9]+)_")[, 2]
  
  dplyr::case_when(
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
  
  # trapezoidal integration
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

# choose tau (you can change this)
selected_tau <- tau_primary # change this in config file

all_results <- list()

for (f in rdata_files) {
  message("Processing ", f)
  
  loaded_name <- load(f)          # should be "res"
  obj <- get(loaded_name)         # retrieve it
  
  # pull time grid and stratified curves out of obj
  tt  <- obj$tt
  ss  <- obj$surv_strat
  
  # sanity check: required elements
  needed <- c("real_1","real_0",
              "unadj_1","unadj_0",
              "iptw_1","iptw_0",
              "ipcw_1","ipcw_0",
              "combo_1","combo_0")
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
  res_comb <- compute_rmst_diff_se(tt, ss$combo_1,  ss$combo_0,  tau = selected_tau)
  
  df <- data.frame(
    method = c("True", "Unadjusted", "IPTW", "IPCW", "IPTW+IPCW"),
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
  levels = c("True", "Unadjusted", "IPTW", "IPCW", "IPTW+IPCW")
)

results_df$ic_level <- factor(
  results_df$ic_level,
  levels = c("none", "weak", "strong"),
  labels = c("No IC", "Weak IC", "Strong IC")
)

# Extract "truth" (True RMST diff) per scenario × IC combo
truth_df <- results_df %>%
  filter(method == "True") %>%
  dplyr::select(scenario, ic_level, true_mean = rmst_diff) %>%
  distinct()

rmst_pl <- ggplot(results_df %>% filter(method != "True"),
                  aes(x = method, y = rmst_diff)) +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
  geom_hline(data = truth_df,
             aes(yintercept = true_mean),
             color = "red", linewidth = 0.6, alpha = 0.5) +
  facet_grid(ic_level ~ scenario) +
  ylab(paste0("RMST Difference (Tx - Ctrl), Tau = ", selected_tau)) +
  xlab("Method") +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    strip.text   = element_text(size = 14),
    legend.position = "none"
  )

ggsave(file.path(results_dir, "sim_small_RMSTdiff.png"), rmst_pl,
       width = 10, height = 8, dpi = 300)
ggsave(file.path(results_dir, "sim_small_RMSTdiff.pdf"), rmst_pl,
       width = 10, height = 8)
