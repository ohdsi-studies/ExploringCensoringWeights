source("sim4_run.R")

# override sample sizes and replications
n_per_rep  <- 2500
n_reps     <- 200

# ---- Run all Sim 4 scenarios ----
ic_keys <- names(ic_strength_levels)  # c("none","weak","strong")

save_with_ic <- function(obj, base_path, ic_key) {
  # Remove .rdata extension if present, add IC suffix
  base_no_ext <- sub("\\.rdata$", "", base_path, ignore.case = TRUE)
  new_path <- paste0(base_no_ext, ic_suffix(ic_key), ".rdata")
  save(obj, file = new_path)
}

# Helper function to get checkpoint file path
get_checkpoint_path <- function(base_path, ic_key) {
  base_no_ext <- sub("\\.rdata$", "", base_path, ignore.case = TRUE)
  paste0(base_no_ext, ic_suffix(ic_key), "_checkpoint.rdata")
}

# Helper function to clean up checkpoint after successful completion
cleanup_checkpoint <- function(checkpoint_path) {
  if (file.exists(checkpoint_path)) {
    file.remove(checkpoint_path)
    message(sprintf("  Cleaned up checkpoint file: %s", checkpoint_path))
  }
}

# Baseline / non-differential censoring
print_sim4_treatment_coefs()

# Differential censoring scenario (e.g., Sim 4C)
# print_sim4_treatment_coefs(phi_Z = 1.5)

print(out_dir)

# 4A
use_bootstrap_se <- F  # Toggle to use bootstrap SEs instead of analytical SEs

message("=== Starting 4A Baseline ===")
for (k in ic_keys) {
  checkpoint_path <- get_checkpoint_path(sim4_files[["4A_baseline"]], k)
  res_4A <- run_scenario("4A_baseline", ic_key = k, checkpoint_file = checkpoint_path)
  save_with_ic(res_4A, sim4_files[["4A_baseline"]], k)
  cleanup_checkpoint(checkpoint_path)
}

# 4B
message("=== Starting 4B Poor Overlap ===")
for (k in ic_keys) {
  checkpoint_path <- get_checkpoint_path(sim4_files[["4B_poor_overlap"]], k)
  res_4B <- run_scenario("4B_poor_overlap", ic_key = k, checkpoint_file = checkpoint_path)
  save_with_ic(res_4B, sim4_files[["4B_poor_overlap"]], k)
  cleanup_checkpoint(checkpoint_path)
}

# 4C: differential informative censoring
# phi_Z values are defined in sim4_config.R as phi_Z_levels
message(sprintf("=== Starting 4C Informative Censoring (phi_Z = %s) ===", 
                paste(phi_Z_levels, collapse = ", ")))

for (k in ic_keys) {
  for (i in (1:seq_along(phi_Z_levels))){s 
    phi_Z_val <- phi_Z_levels[i]
    scen_key <- paste0("4C_diff_censor_", i - 1)  # 0-indexed for backward compatibility
    file_path <- get_4C_file_path(phi_Z_val)
    
    message(sprintf("  Running 4C with phi_Z = %.2f (IC = %s)", phi_Z_val, k))
    checkpoint_path <- get_checkpoint_path(file_path, k)
    res_4C <- run_scenario(scen_key, phi_Z = phi_Z_val, ic_key = k, checkpoint_file = checkpoint_path)
    save_with_ic(res_4C, file_path, k)
    cleanup_checkpoint(checkpoint_path)
  }
}

# 4D
message("=== Starting 4D Model Misspecification ===")
for (k in ic_keys) {
  checkpoint_path <- get_checkpoint_path(sim4_files[["4D_misspecified"]], k)
  res_4D <- run_scenario("4D_misspecified", ic_key = k, checkpoint_file = checkpoint_path)
  save_with_ic(res_4D, sim4_files[["4D_misspecified"]], k)
  cleanup_checkpoint(checkpoint_path)
}

# 4E
message("=== Starting 4E Unmeasured Confounders ===")
# (same IC grid as others, but analysis omits some true confounders)
for (k in ic_keys) {
  checkpoint_path <- get_checkpoint_path(sim4_files[["4E_unmeasured"]], k)
  res_4E <- run_scenario("4E_unmeasured", ic_key = k, checkpoint_file = checkpoint_path)
  save_with_ic(res_4E, sim4_files[["4E_unmeasured"]], k)
  cleanup_checkpoint(checkpoint_path)
}

source("sim4_plot_HR.R")
source("sim4_plot_KM.R")
source("sim4_plot_RMST.R")
source("MCSE_calculations_sim4.R")
