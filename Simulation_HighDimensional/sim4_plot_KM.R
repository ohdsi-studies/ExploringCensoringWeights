############################################################
## sim4_plot_KM.R
##
## For each Sim 4 scenario (4A, 4B, 4C φZ=0.5, 4C φZ=1.0,
## 4D, 4E), generate ONE PNG/PDF:
##   - rows: IC level (none, weak, strong)
##   - cols: estimator (Unadjusted, IPTW, IPCW, IPTW+IPCW)
##   - each panel: 4 KM curves
##        * blue  = Untreated (D=0)
##        * red   = Treated   (D=1)
##        * dashed = REAL (truth from counterfactual)
##        * solid  = ESTIMATED (method-specific)
##
## Uses the 'surv_strat' matrices saved by sim4_run.R:
##   real_0, real_1, unadj_0, unadj_1, iptw_0, ...
############################################################

library(dplyr)
library(tidyr)
library(purrr)

source("sim4_config.R")  # For filename_suffix_to_phi_Z helper function
library(ggplot2)
library(stringr)

results_dir <- "results_sim4"

# Only files that have IC label in name (exclude checkpoint files)
rdata_files <- list.files(
  results_dir,
  pattern = "\\.rdata$",
  full.names = TRUE
)
rdata_files <- rdata_files[
  str_detect(basename(rdata_files), "ic-(none|weak|strong)") &
  !str_detect(basename(rdata_files), "_checkpoint")
]

if (length(rdata_files) == 0L) {
  stop("No Sim4 result files with ic-(none|weak|strong) found in ", results_dir)
}

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
        return(paste0("4C DiffCens =", phiZ_match))  # Fallback: use raw match
      }
    } else {
      warning("Could not extract phiZ pattern from filename: ", nm, ". Using fallback.")
      return("4C DiffCens φZ=?")  # Fallback label
    }
  }
  
  # Other scenarios
  dplyr::case_when(
    str_detect(nm, "Sim4A_baseline")          ~ "4A Baseline",
    str_detect(nm, "Sim4B_poorOverlap")       ~ "4B Poor Overlap",
    str_detect(nm, "Sim4D_misspecified")      ~ "4D Misspecified",
    str_detect(nm, "Sim4E_unmeasured")        ~ "4E Unmeasured Conf.",
    TRUE ~ "Unknown"
  )
}

extract_ic <- function(path) {
  m <- stringr::str_match(basename(path), "ic-(none|weak|strong)")[, 2]
  ifelse(is.na(m), "(not-set)", m)
}

# Method mapping: internal name → pretty column label
method_map <- c(
  unadj = "Unadjusted",
  iptw  = "IPTW",
  ipcw  = "IPCW",
  comb  = "IPTW+IPCW"
)
method_levels <- unname(method_map)

# Helper to take column means safely
col_mean_safe <- function(mat) {
  if (is.null(mat) || length(mat) == 0L) return(NULL)
  colMeans(mat, na.rm = TRUE)
}

# ---------------------------------------------------------
# 1. Collect truth (REAL) survival curves by scenario × IC
# ---------------------------------------------------------

truth_list  <- list()
est_list    <- list()

for (f in rdata_files) {
  message("Processing ", f)
  env <- new.env()
  loaded_name <- load(f, envir = env)
  obj <- env[[loaded_name]]
  
  tt <- obj$tt
  ss <- obj$surv_strat
  
  scen <- scenario_label(f)
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
# So replicate across 'method' factor.
truth_df_rep <- truth_df %>%
  tidyr::crossing(
    method = factor(method_levels, levels = method_levels)
  )

# ---------------------------------------------------------
# Plotting function: one PNG/PDF per scenario label
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
    # scale_x_continuous(
    #   name = "Time",
    #   limits = c(0, 5),
    #   breaks = seq(0, 5, by = 1),
    #   labels = seq(0, 5, by = 1) * 100
    # ) +
    scale_color_manual(
      values = c("Untreated" = "blue", "Treated" = "red"),
      name   = "Treatment Arm"
    ) +
    scale_linetype_manual(
      values = c("Truth" = "dashed", "Estimate" = "solid"),
      name   = ""
    ) +
    coord_cartesian(ylim = c(0, 1), xlim = c(0, 10)) +
    labs(
      x = "Time",
      y = "Survival Probability",
      title = paste0("Sim 4 KM Curves: ", scen_label)
    ) +
    theme_minimal(base_size = 14) +
    theme(
      strip.text   = element_text(size = 12),
      legend.position = "bottom"
    )
  
  # Sanitize filename based on scenario label
  # Preserve decimal points and negative signs for phi_Z values
  # Replace spaces and special chars (except . and -) with underscores
  scen_safe <- scen_label
  # Replace spaces and equals signs with underscores
  scen_safe <- gsub("\\s+", "_", scen_safe)
  scen_safe <- gsub("=", "_", scen_safe)
  # Replace Unicode φ with "phi" for filename compatibility
  scen_safe <- gsub("\u03C6", "phi", scen_safe)
  scen_safe <- gsub("φ", "phi", scen_safe)
  # Remove any remaining non-alphanumeric except . and - and _
  scen_safe <- gsub("[^A-Za-z0-9._-]+", "_", scen_safe)
  # Clean up multiple consecutive underscores
  scen_safe <- gsub("_+", "_", scen_safe)
  scen_safe <- gsub("^_|_$", "", scen_safe)  # Remove leading/trailing underscores
  
  png_file  <- file.path(results_dir, paste0("Sim4_KM_", scen_safe, ".png"))
  pdf_file  <- file.path(results_dir, paste0("Sim4_KM_", scen_safe, ".pdf"))
  
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
