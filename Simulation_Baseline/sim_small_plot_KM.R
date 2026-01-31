############################################################
## sim_small_plot_KM.R
##
## For each sim_small scenario (sim1, sim2, sim3),
## generate ONE PNG/PDF:
##   - rows: IC level (none, weak, strong)
##   - cols: estimator (Unadjusted, IPTW, IPCW, IPTW+IPCW)
##   - each panel: 4 KM curves
##        * blue   = Untreated (D=0)
##        * red    = Treated   (D=1)
##        * dashed = REAL (truth from counterfactual)
##        * solid  = ESTIMATED (method-specific)
##
## Uses the 'surv_strat' matrices saved by run_small.R:
##   real_0, real_1, unadj_0, unadj_1, iptw_0, ...
############################################################

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(stringr)

results_dir <- "results_sim_small"

# Files should look like: res_sim1_none.rdata, res_sim2_weak.rdata, ...
rdata_files <- list.files(
  results_dir,
  pattern = "^res_sim[123]_.*\\.rdata$",
  full.names = TRUE
)

if (length(rdata_files) == 0L) {
  stop("No sim_small result files res_sim[123]_*.rdata found in ", results_dir)
}

# -------------------------
# Scenario & IC labellers
# -------------------------

scenario_label <- function(path) {
  nm <- basename(path)
  scen_code <- str_match(nm, "^res_(sim[123])_")[, 2]
  case_when(
    scen_code == "sim1" ~ "Sim 1: Nondiff censoring",
    scen_code == "sim2" ~ "Sim 2: Diff censoring",
    scen_code == "sim3" ~ "Sim 3: PH violation",
    TRUE                ~ "Unknown"
  )
}

extract_ic <- function(path) {
  m <- stringr::str_match(basename(path), "_(none|weak|strong)\\.rdata$")[, 2]
  ifelse(is.na(m), "(not-set)", m)
}

# Method mapping: internal name → pretty label
method_map <- c(
  unadj = "Unadjusted",
  iptw  = "IPTW",
  ipcw  = "IPCW",
  combo = "IPTW+IPCW"
)
method_levels <- unname(method_map)

# Helper to take column means safely
col_mean_safe <- function(mat) {
  if (is.null(mat) || length(mat) == 0L) return(NULL)
  colMeans(mat, na.rm = TRUE)
}

# ---------------------------------------------------------
# 1. Collect truth (REAL) and estimated survival curves
# ---------------------------------------------------------

truth_list  <- list()
est_list    <- list()

for (f in rdata_files) {
  message("Processing ", f)
  env <- new.env()
  loaded_name <- load(f, envir = env)  # usually "res"
  obj <- env[[loaded_name]]
  
  tt <- obj$tt
  ss <- obj$surv_strat
  
  scen   <- scenario_label(f)
  ic_raw <- extract_ic(f)
  
  # REAL curves (truth; averaged over reps)
  s_real_0 <- col_mean_safe(ss$real_0)
  s_real_1 <- col_mean_safe(ss$real_1)
  
  if (is.null(s_real_0) || is.null(s_real_1)) {
    warning("Missing real_0 / real_1 in file ", f, ", skipping.")
    next
  }
  
  truth_list[[length(truth_list) + 1L]] <- bind_rows(
    data.frame(
      scenario = scen,
      ic_level = ic_raw,
      time     = tt,
      surv     = s_real_0,
      arm      = "Untreated",
      type     = "Truth",
      stringsAsFactors = FALSE
    ),
    data.frame(
      scenario = scen,
      ic_level = ic_raw,
      time     = tt,
      surv     = s_real_1,
      arm      = "Treated",
      type     = "Truth",
      stringsAsFactors = FALSE
    )
  )
  
  # Estimated curves per method, averaged over reps
  for (m in names(method_map)) {
    s0_name <- paste0(m, "_0")
    s1_name <- paste0(m, "_1")
    
    s_est_0 <- col_mean_safe(ss[[s0_name]])
    s_est_1 <- col_mean_safe(ss[[s1_name]])
    
    if (is.null(s_est_0) || is.null(s_est_1)) {
      next
    }
    
    est_list[[length(est_list) + 1L]] <- bind_rows(
      data.frame(
        scenario = scen,
        ic_level = ic_raw,
        method   = method_map[[m]],
        time     = tt,
        surv     = s_est_0,
        arm      = "Untreated",
        type     = "Estimate",
        stringsAsFactors = FALSE
      ),
      data.frame(
        scenario = scen,
        ic_level = ic_raw,
        method   = method_map[[m]],
        time     = tt,
        surv     = s_est_1,
        arm      = "Treated",
        type     = "Estimate",
        stringsAsFactors = FALSE
      )
    )
  }
}

truth_df <- bind_rows(truth_list)
est_df   <- bind_rows(est_list)


if (nrow(truth_df) == 0L || nrow(est_df) == 0L) {
  stop("No truth or estimate curves were constructed – check result files.")
}

# ---------------------------------------------------------
# Factor ordering / labels
# ---------------------------------------------------------

truth_df$ic_level <- factor(
  truth_df$ic_level,
  levels = c("none", "weak", "strong"),
  labels = c("No IC", "Weak IC", "Strong IC")
)

est_df$ic_level <- factor(
  est_df$ic_level,
  levels = c("none", "weak", "strong"),
  labels = c("No IC", "Weak IC", "Strong IC")
)

est_df$method <- factor(
  est_df$method,
  levels = method_levels
)

truth_df$arm <- factor(truth_df$arm, levels = c("Untreated", "Treated"))
est_df$arm   <- factor(est_df$arm,   levels = c("Untreated", "Treated"))

truth_df$type <- factor(truth_df$type, levels = c("Truth", "Estimate"))
est_df$type   <- factor(est_df$type,   levels = c("Truth", "Estimate"))

# We want TRUTH curves to appear in every method-column panel.
truth_df_rep <- truth_df %>%
  tidyr::crossing(
    method = factor(method_levels, levels = method_levels)
  )

# ---------------------------------------------------------
# Plotting function: one PNG/PDF per scenario
# ---------------------------------------------------------

make_km_plot_for_scenario <- function(scen_label) {
  df_est   <- est_df %>% filter(scenario == scen_label)
  df_truth <- truth_df_rep %>% filter(scenario == scen_label)
  
  if (nrow(df_est) == 0L || nrow(df_truth) == 0L) {
    warning("No data for scenario '", scen_label, "', skipping.")
    return(invisible(NULL))
  }
  
  pl <- ggplot() +
    # REAL curves (dashed)
    geom_step(
      data = df_truth,
      aes(x = time, y = surv,
          color = arm,
          linetype = type),
      linewidth = 0.7,
      alpha = 0.9
    ) +
    # ESTIMATED curves (solid)
    geom_step(
      data = df_est,
      aes(x = time, y = surv,
          color = arm,
          linetype = type),
      linewidth = 0.7
    ) +
    facet_grid(ic_level ~ method) +
    scale_color_manual(
      values = c("Untreated" = "blue", "Treated" = "red"),
      name   = "Treatment Arm"
    ) +
    scale_linetype_manual(
      values = c("Truth" = "dashed", "Estimate" = "solid"),
      name   = ""
    ) +
    coord_cartesian(ylim = c(0, 1), xlim = c(0, 5)) +  # sim_small horizon
    labs(
      x = "Time",
      y = "Survival Probability",
      title = scen_label
    ) +
    theme_minimal(base_size = 14) +
    theme(
      strip.text      = element_text(size = 12),
      legend.position = "bottom"
    )
  
  scen_safe <- gsub("[^A-Za-z0-9]+", "_", scen_label)
  png_file  <- file.path(results_dir, paste0("sim_small_KM_", scen_safe, ".png"))
  pdf_file  <- file.path(results_dir, paste0("sim_small_KM_", scen_safe, ".pdf"))
  
  ggsave(png_file, pl, width = 10, height = 7, dpi = 300)
  ggsave(pdf_file, pl, width = 10, height = 7)
  
  message("Saved KM figure for scenario '", scen_label, "' to: ",
          "\n  ", png_file,
          "\n  ", pdf_file)
}

# ---------------------------------------------------------
# Generate ONE PNG/PDF per scenario
# ---------------------------------------------------------
all_scenarios <- sort(unique(est_df$scenario))

for (sc in all_scenarios) {
  make_km_plot_for_scenario(sc)
}

message("Done. KM plots saved in ", results_dir)
