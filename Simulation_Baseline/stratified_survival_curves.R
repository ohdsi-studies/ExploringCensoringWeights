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
      # se = ssf$std.err,  # Add SE column
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

calc.surv.unadj.stratified <- function(times, status, tt, data, tau = NULL) {
  
  # Initialize an empty list to store survival curves for D = 0 and D = 1
  surv_results <- list()
  if (is.null(tau)) tau <- max(tt)
  
  for (d in unique(data$D)) {
    # Subset data for each treatment group
    data.subset <- data[data$D == d, ]
    
    # Fit the unadjusted survival model
    surv <- survfit(Surv(ti, di) ~ 1, data = data.subset)
    
    # Estimate survival probabilities at specified time points
    ssf <- summary(surv, times = tt, extend = TRUE)
    
    # std.err from survfit is on the log(-log(S)) scale; convert to survival scale
    se_surv <- ssf$surv * ssf$std.err
    
    # RMST and SE from survfit (restricted mean), at tau
    rmst_tab <- survival:::survmean(surv, rmean = tau, scale = 1)$matrix
    rmst_val <- unname(rmst_tab["rmean"])
    rmst_se  <- unname(rmst_tab["se(rmean)"])
    
    # Store results with treatment group label
    surv_results[[as.character(d)]] <- data.frame(
      time = tt,
      survival = ssf$surv,
      se = se_surv,  # SE on survival scale
      D = d,  # Add treatment group
      rmean = rmst_val,
      rmean_se = rmst_se
    )
  }
  
  # Combine results into a single dataframe
  surv_df <- do.call(rbind, surv_results)
  
  return(surv_df)
}
  

#####################################
## Stratified IPTW Survival Curves ##
#####################################

calc.surv.IPW.stratified <- function(times, status, tt, IPW.weights, data, tau = NULL) {
  
  # times/status are not actually needed here; we use data$ti/di
  surv_results <- list()
  if (is.null(tau)) tau <- max(tt)
  
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
    
    se_surv <- ssf$surv * ssf$std.err
    
    # RMST and SE from survfit (restricted mean)
    rmst_tab <- survival:::survmean(surv, rmean = tau, scale = 1)$matrix
    rmst_val <- unname(rmst_tab["rmean"])
    rmst_se  <- unname(rmst_tab["se(rmean)"])
    
    surv_results[[as.character(d)]] <- data.frame(
      time     = tt,
      survival = ssf$surv,
      se       = se_surv,  # SE on survival scale
      D        = d,
      rmean    = rmst_val,
      rmean_se = rmst_se
    )
  }
  
  do.call(rbind, surv_results)
}

#####################################
## Stratified IPCW Survival Curves ##
#####################################
calc.surv.IPCW.stratified <- function(Tstart, Tstop, status, tt, IPCW.weights, data.long, tau = NULL) {
  
  surv_results <- list()
  if (is.null(tau)) tau <- max(tt)
  
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
    
    se_surv <- ssurv.IPCW$surv * ssurv.IPCW$std.err
    
    # RMST and SE from survfit (restricted mean)
    rmst_tab <- survival:::survmean(surv.IPCW, rmean = tau, scale = 1)$matrix
    rmst_val <- unname(rmst_tab["rmean"])
    rmst_se  <- unname(rmst_tab["se(rmean)"])
    
    surv_results[[as.character(d)]] <- data.frame(
      time     = tt,
      survival = ssurv.IPCW$surv,
      se       = se_surv,  # SE on survival scale
      D        = d,
      rmean    = rmst_val,
      rmean_se = rmst_se
    )
  }
  
  do.call(rbind, surv_results)
}


#####################################
## Stratified COMB Survival Curves ##
#####################################
calc.surv.comb.stratified <- function(Tstart, Tstop, status, tt, comb.weights, data.long, tau = NULL) {
  
  surv_results <- list()
  if (is.null(tau)) tau <- max(tt)
  
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
    
    se_surv <- ssurv.comb$surv * ssurv.comb$std.err
    
    # RMST and SE from survfit (restricted mean)
    rmst_tab <- survival:::survmean(surv.comb, rmean = tau, scale = 1)$matrix
    rmst_val <- unname(rmst_tab["rmean"])
    rmst_se  <- unname(rmst_tab["se(rmean)"])
    
    surv_results[[as.character(d)]] <- data.frame(
      time     = tt,
      survival = ssurv.comb$surv,
      se       = se_surv,  # SE on survival scale
      D        = d,
      rmean    = rmst_val,
      rmean_se = rmst_se
    )
  }
  
  do.call(rbind, surv_results)
}


