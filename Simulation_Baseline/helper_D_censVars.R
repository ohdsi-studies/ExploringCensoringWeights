## replicating Willems et al
# step 2 to replicating code -> add treamtnet variable D 
# on top of age and gender

# function to estimate survival probabilities at time points tt
# can be used to estimate the real survival curve
calc.surv <- function(times, status, tt, data)
{
  # Fit model
  #cox <- coxpSurv(times, status) ~ ZAge + Gender, data = data)
  surv <- survfit(Surv(times, status) ~ 1, data = data)
  # Estimate survival probabilities at time points tt
  ssf <- summary(surv, times = tt, extend = TRUE)
  # Return results:
  return(ssf$surv)
}

# fx to fit Cox model
calc.beta <- function(times, status, data)
{
  # Fit model:
  Cox <- coxph(Surv(times, status) ~ D, data = data)
  # Return coefficient and SE:
  coef_summary <- summary(Cox)$coef
  return(c(coef = coef_summary[1, 1], se = coef_summary[1, 3]))
}

calc.surv.IPW <- function(times, status, tt, IPW.weights, data)
{
  # Fit model
  #cox <- coxpSurv(times, status) ~ ZAge + Gender, data = data)
  surv <- survfit(Surv(times, status) ~ 1, weights = IPW.weights, 
                  data = data)
  # Estimate survival probabilities at time points tt
  ssf <- summary(surv, times = tt, extend = TRUE)
  # Return results:
  return(ssf$surv)
}

# fx to fit IPW weights - Cox model
calc.beta.IPW <- function(times, status, IPW.weights, data)
{
  # Fit model:
  Cox <- coxph(Surv(times, status) ~ as.factor(D),
               data = data, weights = IPW.weights)
  # Return coefficient and SE:
  coef_summary <- summary(Cox)$coef
  return(c(coef = coef_summary[1, 1], se = coef_summary[1, 3]))
}

# transform data into long format
transform.data <- function(data, cut.times = data$ti)
{
  # Define Tstart and the indicator "censored":
  data$Tstart <- 0
  data$censored <- 1 - data$di
  
  # Times at which to split the intervals:
  # cut.times <- data$ti
  
  # Split data with event = di (event):
  data.long <- survSplit(data = data,
                         cut = cut.times,
                         end = "ti",
                         start = "Tstart",
                         event = "di")
  data.long <- data.long[order(data.long$ID,
                               data.long$ti),]
  
  # Split data with event = censored (censoring):
  data.long.cens <- survSplit(data,
                              cut=cut.times,
                              end="ti",
                              start="Tstart",
                              event="censored")
  data.long.cens <- data.long.cens[order(data.long.cens$ID,
                                         data.long.cens$ti),]
  
  # Add "censored" indicator to long data format:
  data.long$censored <- data.long.cens$censored
  data.long$ID <- as.numeric(data.long$ID)
  # Return long data format:
  return(data.long)
}


# calculate IPCW stabilized and unstabilized weights
calc.IPCW_fast <- function(C0, CZ, data.long, p_trunc = 1) {
  
  # ⁠K0(t)⁠ – null model survival
  baseC0  <- basehaz(C0, centered = FALSE)
  H0_fun  <- stats::approxfun(baseC0$time, baseC0$hazard, rule = 2)
  data.long$K0ti <- exp(-H0_fun(data.long$Tstart))
  
  # ⁠KZ(t)⁠ – covariate model survival
  baseCZ  <- basehaz(CZ, centered = FALSE)
  HZ_fun  <- stats::approxfun(baseCZ$time, baseCZ$hazard, rule = 2)
  linpred <- predict(CZ, type = "lp", newdata = data.long)
  data.long$KZti <- exp(-HZ_fun(data.long$Tstart) * exp(linpred))
  
  # Weights
  data.long$WUnStab <- 1 / pmax(data.long$KZti, 1e-8)
  data.long$WStab   <- data.long$K0ti / pmax(data.long$KZti, 1e-8)
  
  # Truncate
  cap <- quantile(data.long$WStab, p_trunc, na.rm = TRUE)
  data.long$WStab   <- pmin(data.long$WStab, cap)
  data.long$WUnStab <- pmin(data.long$WUnStab, cap * mean(data.long$WUnStab / data.long$WStab, na.rm = TRUE))
  
  return(data.long)
}

# function to estimate survival curve with weighted subjects
calc.surv.IPCW <- function(Tstart, Tstop, status, tt, IPCW.weights,
                           data.long)
{
  # Fit the Product-Limit estimator with weighted subjects
  surv.IPCW <- survfit(Surv(Tstart, Tstop, status) ~ 1, data = data.long,
                       weights = IPCW.weights, timefix=FALSE)
  # Estimate survival probabilities at time points tt
  ssurv.IPCW <- summary(surv.IPCW, times=tt, extend = TRUE)
  # Return resulting survival probabilities:
  return(ssurv.IPCW$surv)
}

# fit cox model for time to event with weighted subjects
calc.beta.IPCW <- function(Tstart, Tstop, status, IPCW.weights, data.long)
{
  # Fit model:
  Cox <- coxph(Surv(Tstart, Tstop, status,
                    type = "counting") ~ D,
               data = data.long,
               weights = IPCW.weights, 
               control = coxph.control(timefix = FALSE))
  # Return coefficient and SE:
  coef_summary <- summary(Cox)$coef
  return(c(coef = coef_summary[1, 1], se = coef_summary[1, 3]))
}


# function to calculate IPW + IPCW (comb.weights)
calc.surv.comb <- function(Tstart, Tstop, status, tt, comb.weights,
                           data.long)
{
  # Fit the Product-Limit estimator with weighted subjects
  surv.comb <- survfit(Surv(Tstart, Tstop, status) ~ 1, data = data.long,
                       weights = comb.weights, timefix=FALSE)
  # Estimate survival probabilities at time points tt
  ssurv.comb <- summary(surv.comb, times=tt, extend = TRUE)
  # Return resulting survival probabilities:
  return(ssurv.comb$surv)
}

# fx to fit IPW + IPCW weights - Cox model
calc.beta.comb <- function(Tstart, Tstop, status, comb.weights, data.long)
{
  # Fit model:
  Cox <- coxph(Surv(Tstart, Tstop, status,
                    type = "counting") ~ D,
               data = data.long,
               weights = comb.weights,
               control = coxph.control(timefix = FALSE))
  # Return coefficient and SE:
  coef_summary <- summary(Cox)$coef
  return(c(coef = coef_summary[1, 1], se = coef_summary[1, 3]))
}
