# stratified survival curves

library(survival)

#####################################
## Stratified REAL Survival Curves ##
#####################################

calc.surv.real.stratified <- function(times, status, tt, data) {
  
  # Initialize an empty list to store survival curves for D = 0 and D = 1
  surv_results <- list()
  
  for (d in unique(data$D)) {
    # Subset data for each treatment group
    data.subset <- data[data$D == d, ]
    
    # Fit the unadjusted survival model
    surv <- survfit(Surv(xi, 1) ~ 1, data = data.subset)
    
    # Estimate survival probabilities at specified time points
    ssf <- summary(surv, times = tt, extend = TRUE)
    
    # Store results with treatment group label
    surv_results[[as.character(d)]] <- data.frame(
      time = tt,
      survival = ssf$surv,
      D = d  # Add treatment group
    )
  }
  
  # Combine results into a single dataframe
  surv_df <- do.call(rbind, surv_results)
  
  return(surv_df)
}


#####################################
## Stratified unadjusted Survival Curves ##
#####################################

calc.surv.unadj.stratified <- function(times, status, tt, data, tau = max(tt, na.rm = TRUE)) {
  
  # Initialize an empty list to store survival curves for D = 0 and D = 1
  surv_results <- list()
  
  for (d in unique(data$D)) {
    # Subset data for each treatment group
    data.subset <- data[data$D == d, ]
    
    # Fit the unadjusted survival model
    surv <- survfit(Surv(ti, di) ~ 1, data = data.subset)
    
    # Estimate survival probabilities at specified time points
    ssf <- summary(surv, times = tt, extend = TRUE)

    # RMST and SE from survfit (restricted mean)
    rmst_tab <- survival:::survmean(surv, rmean = tau, scale = 1)$matrix
    rmst_val <- unname(rmst_tab["rmean"])
    rmst_se  <- unname(rmst_tab["se(rmean)"])
    
    # Store results with treatment group label
    surv_results[[as.character(d)]] <- data.frame(
      time = tt,
      survival = ssf$surv,
      std_err = if (!is.null(ssf$std.err)) ssf$std.err else rep(NA_real_, length(tt)),
      rmean    = rmst_val,
      rmean_se = rmst_se,
      D = d  # Add treatment group
    )
  }
  
  # Combine results into a single dataframe
  surv_df <- do.call(rbind, surv_results)
  
  return(surv_df)
}
  
# surv_unadj_strat <- calc.surv.unadj.stratified(
#   times = data$ti, 
#   status = data$di, 
#   tt = tt, 
#   data = data
# )
# 
# 
# ggplot(surv_unadj_strat, aes(x = time, y = survival, color = as.factor(D))) +
#   geom_step() +
#   labs(x = "Time", y = "Survival Probability", color = "Treatment (D)") +
#   theme_minimal() +
#   scale_color_manual(values = c("blue", "red"))


#####################################
## Stratified IPTW Survival Curves ##
#####################################

calc.surv.IPW.stratified <- function(times, status, tt, IPW.weights, data, tau = max(tt, na.rm = TRUE)) {
  
  # times/status are not actually needed here; we use data$ti/di
  surv_results <- list()
  
  for (d in unique(data$D)) {
    # Subset data and weights for each treatment group
    idx <- data$D == d
    data.subset <- data[idx, ]
    w.subset <- IPW.weights[idx]
    
    # Fit weighted KM
    surv <- survfit(
      Surv(ti, di) ~ 1,
      weights = w.subset,
      data    = data.subset
    )
    
    ssf <- summary(surv, times = tt, extend = TRUE)

    # RMST and SE from survfit (restricted mean)
    rmst_tab <- survival:::survmean(surv, rmean = tau, scale = 1)$matrix
    rmst_val <- unname(rmst_tab["rmean"])
    rmst_se  <- unname(rmst_tab["se(rmean)"])
    
    surv_results[[as.character(d)]] <- data.frame(
      time     = tt,
      survival = ssf$surv,
      std_err  = if (!is.null(ssf$std.err)) ssf$std.err else rep(NA_real_, length(tt)),
      rmean    = rmst_val,
      rmean_se = rmst_se,
      D        = d
    )
  }
  
  do.call(rbind, surv_results)
}

#####################################
## Stratified IPCW Survival Curves ##
#####################################
calc.surv.IPCW.stratified <- function(Tstart, Tstop, status, tt, IPCW.weights, data.long, tau = max(tt, na.rm = TRUE)) {
  
  surv_results <- list()
  
  for (d in unique(data.long$D)) {
    idx <- data.long$D == d
    data.subset <- data.long[idx, ]
    w.subset    <- IPCW.weights[idx]
    
    surv.IPCW <- survfit(
      Surv(Tstart, ti, di) ~ 1,
      data    = data.subset,
      weights = w.subset,
      timefix = FALSE
    )
    
    ssurv.IPCW <- summary(surv.IPCW, times = tt, extend = TRUE)

    # RMST and SE from survfit (restricted mean)
    rmst_tab <- survival:::survmean(surv.IPCW, rmean = tau, scale = 1)$matrix
    rmst_val <- unname(rmst_tab["rmean"])
    rmst_se  <- unname(rmst_tab["se(rmean)"])
    
    surv_results[[as.character(d)]] <- data.frame(
      time     = tt,
      survival = ssurv.IPCW$surv,
      std_err  = if (!is.null(ssurv.IPCW$std.err)) ssurv.IPCW$std.err else rep(NA_real_, length(tt)),
      rmean    = rmst_val,
      rmean_se = rmst_se,
      D        = d
    )
  }
  
  do.call(rbind, surv_results)
}


#####################################
## Stratified COMB Survival Curves ##
#####################################
calc.surv.comb.stratified <- function(Tstart, Tstop, status, tt, comb.weights, data.long, tau = max(tt, na.rm = TRUE)) {
  
  surv_results <- list()
  
  for (d in unique(data.long$D)) {
    idx <- data.long$D == d
    data.subset <- data.long[idx, ]
    w.subset    <- comb.weights[idx]
    
    surv.comb <- survfit(
      Surv(Tstart, ti, di) ~ 1,
      data    = data.subset,
      weights = w.subset,
      timefix = FALSE
    )
    
    ssurv.comb <- summary(surv.comb, times = tt, extend = TRUE)

    # RMST and SE from survfit (restricted mean)
    rmst_tab <- survival:::survmean(surv.comb, rmean = tau, scale = 1)$matrix
    rmst_val <- unname(rmst_tab["rmean"])
    rmst_se  <- unname(rmst_tab["se(rmean)"])
    
    surv_results[[as.character(d)]] <- data.frame(
      time     = tt,
      survival = ssurv.comb$surv,
      std_err  = if (!is.null(ssurv.comb$std.err)) ssurv.comb$std.err else rep(NA_real_, length(tt)),
      rmean    = rmst_val,
      rmean_se = rmst_se,
      D        = d
    )
  }
  
  do.call(rbind, surv_results)
}


