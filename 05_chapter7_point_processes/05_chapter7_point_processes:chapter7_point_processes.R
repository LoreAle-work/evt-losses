# ============================================================
# Chapter 7: Point Process Models for Extremes
# Point Process Characterization of SPY Daily Losses
# ============================================================
#
# This script applies the Chapter 7 point process approach to
# daily SPY losses.
#
# Earlier chapters used:
#
#   Chapter 3: GEV model for block maxima
#   Chapter 4: GPD model for threshold exceedances
#   Chapter 5: declustering for dependent extremes
#   Chapter 6: non-stationary GPD models with volatility
#
# Chapter 7 connects block maxima and threshold exceedances through
# the point process representation of extremes.
#
# Main idea:
#
#   Extreme observations are represented as points in time-value space:
#
#       (time, loss)
#
#   above a high threshold.
#
# The stationary point process model estimates GEV-compatible parameters:
#
#   mu, sigma, xi
#
# using threshold exceedances.
#
# This revised version fits four point process models:
#
#   1. Stationary point process
#      mu constant, sigma constant, xi constant
#
#   2. Volatility-dependent scale point process
#      mu constant
#      log(sigma_t) = beta_0 + beta_vol * z_t
#      xi constant
#
#   3. Volatility-dependent location point process
#      mu_t = mu_0 + mu_vol * z_t
#      sigma constant
#      xi constant
#
#   4. Volatility-dependent location + scale point process
#      mu_t = mu_0 + mu_vol * z_t
#      log(sigma_t) = beta_0 + beta_vol * z_t
#      xi constant
#
# Why add location?
#
#   The previous volatility scale-only model improved fit, but produced
#   a negative volatility effect on scale and decreasing conditional
#   return levels. That may happen because volatility shifts the overall
#   level of extremes, not only their spread. So now we let mu move too.
#
# This script includes:
#
#   1. Data download and loss construction
#   2. Volatility covariate construction
#   3. Threshold exceedance identification
#   4. Stationary point process model
#   5. Three volatility-dependent point process models
#   6. AIC, BIC, and likelihood ratio comparisons
#   7. Return level estimation
#   8. Comparison with Chapter 3 GEV and Chapter 4 GPD
#   9. Diagnostic plots
#
# ============================================================


# ============================================================
# 0. Package Setup
# ============================================================

# Run once if needed:
# install.packages(c("quantmod", "ismev", "xts", "zoo"))

library(quantmod)
library(ismev)
library(xts)
library(zoo)


# ============================================================
# 1. Project Folder Setup
# ============================================================


setwd("~/Desktop/UNI/Projects/EVT/evt-losses")

cat("Current working directory:\n")
print(getwd())

if (basename(getwd()) == "05_chapter7_point_processes") {
  chapter_dir <- "."
} else {
  chapter_dir <- "05_chapter7_point_processes"
}

fig_dir <- file.path(chapter_dir, "figures")
res_dir <- file.path(chapter_dir, "results")

if (chapter_dir != ".") {
  dir.create(chapter_dir, showWarnings = FALSE, recursive = TRUE)
}

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)


# ============================================================
# 2. Download SPY Prices
# ============================================================

# We use SPY adjusted closing prices from Yahoo Finance.
# Adjusted prices account for dividends and splits.
#
# For perfect reproducibility later, you can set a fixed end date:
#
# getSymbols("SPY", src = "yahoo", from = "1993-01-01",
#            to = "2026-06-25", auto.assign = TRUE)

getSymbols("SPY", src = "yahoo", from = "1993-01-01", auto.assign = TRUE)

prices <- Ad(SPY)

cat("\nPrice data:\n")
print(head(prices))
print(tail(prices))

png(file.path(fig_dir, "spy_adjusted_prices.png"), width = 1000, height = 600)
plot(prices, main = "SPY Adjusted Closing Prices")
dev.off()


# ============================================================
# 3. Compute Daily Log Returns
# ============================================================

# Daily log returns:
#
#   R_t = log(P_t) - log(P_{t-1})

returns <- diff(log(prices))
returns <- na.omit(returns)
colnames(returns) <- "return"

ret_num <- as.numeric(returns)
dates <- index(returns)

cat("\nReturn summary:\n")
print(summary(ret_num))

png(file.path(fig_dir, "spy_daily_log_returns.png"), width = 1000, height = 600)
plot(returns, main = "SPY Daily Log Returns")
dev.off()


# ============================================================
# 4. Convert Returns into Losses
# ============================================================

# Losses:
#
#   L_t = -R_t
#
# This turns large negative returns into large positive losses.
# EVT is then applied to the upper tail of the loss distribution.

losses <- -returns
colnames(losses) <- "loss"

losses_num <- as.numeric(losses)

cat("\nLoss summary:\n")
print(summary(losses_num))

cat("\nLargest daily losses (%):\n")
print(100 * sort(losses_num, decreasing = TRUE)[1:10])

png(file.path(fig_dir, "spy_daily_losses.png"), width = 1000, height = 600)
plot(losses, main = "SPY Daily Losses")
dev.off()


# ============================================================
# 5. Build Lagged Volatility Covariate
# ============================================================

# We construct the same volatility covariate used in Chapter 6:
#
#   lagged 21-day realized volatility
#
# The rolling standard deviation uses the previous 21 trading days.
# Then we lag the volatility by one day to avoid look-ahead bias.

rolling_sd_21 <- rollapply(
  ret_num,
  width = 21,
  FUN = sd,
  align = "right",
  fill = NA
)

rolling_sd_21_lag <- c(NA, rolling_sd_21[-length(rolling_sd_21)])

df <- data.frame(
  date = dates,
  return = ret_num,
  loss = -ret_num,
  rv21_lag = rolling_sd_21_lag
)

df <- df[complete.cases(df), ]

# Log-transform volatility because raw volatility is positive and skewed.
# Standardize so one unit equals one standard deviation.

df$log_vol <- log(df$rv21_lag)
df$z_vol <- as.numeric(scale(df$log_vol))

# Scaled time covariate, useful for diagnostics or later extensions.

df$t_scaled <- as.numeric(scale(seq_len(nrow(df))))

n_total <- nrow(df)
trading_days_per_year <- 252
n_years <- n_total / trading_days_per_year

cat("\nUsable observations after volatility construction:", n_total, "\n")
cat("Approximate number of trading years:", n_years, "\n")

write.csv(
  df,
  file.path(res_dir, "analysis_dataset_with_covariates.csv"),
  row.names = FALSE
)

png(file.path(fig_dir, "lagged_21day_volatility.png"), width = 1000, height = 600)

plot(
  df$date,
  df$rv21_lag,
  type = "l",
  col = "black",
  main = "Lagged 21-Day Realized Volatility",
  xlab = "Date",
  ylab = "Lagged 21-Day Volatility"
)

dev.off()


# ============================================================
# 6. Threshold Selection
# ============================================================

# We use the same 97.5% empirical loss threshold used in Chapters 4-6.
#
# This makes Chapter 7 directly comparable with the previous threshold
# models.

threshold_probability <- 0.975
u <- as.numeric(quantile(df$loss, threshold_probability, names = FALSE))

df$exceed <- df$loss > u
df$excess <- ifelse(df$exceed, df$loss - u, NA)

n_exceedances <- sum(df$exceed)
exceedance_fraction <- mean(df$exceed)

exceed_df <- df[df$exceed, ]

cat("\nThreshold summary:\n")
cat("Threshold probability:", threshold_probability, "\n")
cat("Threshold:", u, "\n")
cat("Threshold percent:", 100 * u, "\n")
cat("Number of exceedances:", n_exceedances, "\n")
cat("Exceedance fraction:", exceedance_fraction, "\n")

threshold_summary <- data.frame(
  threshold_probability = threshold_probability,
  threshold_value = u,
  threshold_percent = 100 * u,
  total_observations = n_total,
  approximate_years = n_years,
  number_exceedances = n_exceedances,
  exceedance_fraction = exceedance_fraction
)

write.csv(
  threshold_summary,
  file.path(res_dir, "threshold_summary.csv"),
  row.names = FALSE
)


# ============================================================
# 7. Plot Extreme Points
# ============================================================

# Chapter 7 views extremes as points:
#
#   (time, loss)
#
# above a high threshold.

png(file.path(fig_dir, "point_process_extreme_points.png"), width = 1100, height = 650)

plot(
  df$date,
  df$loss,
  type = "h",
  col = "gray70",
  lwd = 0.5,
  ylim = c(0, max(df$loss) * 1.05),
  main = "Point Process View of SPY Extreme Losses",
  xlab = "Date",
  ylab = "Daily Loss"
)

abline(h = u, col = "red", lwd = 3, lty = 2)

points(
  exceed_df$date,
  exceed_df$loss,
  pch = 19,
  col = "blue"
)

legend(
  "topright",
  legend = c(
    paste0("97.5% threshold = ", round(100 * u, 2), "%"),
    "Extreme points"
  ),
  col = c("red", "blue"),
  lty = c(2, NA),
  lwd = c(3, NA),
  pch = c(NA, 19),
  bty = "n"
)

dev.off()


png(file.path(fig_dir, "extreme_points_time_value_space.png"), width = 900, height = 600)

plot(
  as.numeric(exceed_df$date),
  exceed_df$loss,
  pch = 19,
  col = "blue",
  xlab = "Time",
  ylab = "Loss",
  main = "Extreme Points Above the Threshold",
  xaxt = "n"
)

axis.Date(
  side = 1,
  at = pretty(exceed_df$date),
  format = "%Y"
)

abline(h = u, col = "red", lwd = 2, lty = 2)

dev.off()


# ============================================================
# 8. Helper Functions for Point Process Likelihood
# ============================================================

# ------------------------------------------------------------
# Stationary point process likelihood
# ------------------------------------------------------------
#
# Parameters:
#
#   mu, sigma, xi
#
# with sigma > 0.
#
# Annual tail measure above x:
#
#   Lambda(x) = [1 + xi * ((x - mu) / sigma)]^(-1 / xi)
#
# for xi != 0.
#
# Integrated intensity above threshold u over the sample:
#
#   sum_t Lambda_t(u) / 252
#
# In the stationary case Lambda_t(u) is constant, so:
#
#   n * Lambda(u) / 252
#
# Each exceedance x contributes:
#
#   lambda(x) = (1 / sigma)
#               [1 + xi * ((x - mu) / sigma)]^(-1 / xi - 1)
#
# The same time-scale convention is used across all models, so AIC,
# BIC, and LR comparisons are meaningful.

pp_stationary_nll <- function(par, data, threshold, trading_days_per_year = 252) {
  
  mu <- par[1]
  log_sigma <- par[2]
  xi <- par[3]
  
  sigma <- exp(log_sigma)
  
  if (!is.finite(sigma) || sigma <= 0) {
    return(1e10)
  }
  
  x_exceed <- data$loss[data$loss > threshold]
  
  if (abs(xi) < 1e-6) {
    
    tail_intensity_all <- exp(-((threshold - mu) / sigma))
    
    integrated_intensity <- nrow(data) * tail_intensity_all / trading_days_per_year
    
    log_intensity <- -log(sigma) - ((x_exceed - mu) / sigma)
    
  } else {
    
    threshold_support <- 1 + xi * ((threshold - mu) / sigma)
    exceed_support <- 1 + xi * ((x_exceed - mu) / sigma)
    
    if (
      !is.finite(threshold_support) ||
      threshold_support <= 0 ||
      any(!is.finite(exceed_support)) ||
      any(exceed_support <= 0)
    ) {
      return(1e10)
    }
    
    tail_intensity_all <- threshold_support^(-1 / xi)
    
    integrated_intensity <- nrow(data) * tail_intensity_all / trading_days_per_year
    
    log_intensity <- -log(sigma) -
      (1 / xi + 1) * log(exceed_support)
  }
  
  nll <- integrated_intensity - sum(log_intensity)
  
  if (!is.finite(nll)) {
    return(1e10)
  }
  
  return(nll)
}


# ------------------------------------------------------------
# Volatility-dependent scale point process likelihood
# ------------------------------------------------------------
#
# Model:
#
#   mu constant
#   log(sigma_t) = beta_0 + beta_vol * z_t
#   xi constant

pp_vol_scale_nll <- function(par, data, threshold, trading_days_per_year = 252) {
  
  mu <- par[1]
  beta_0 <- par[2]
  beta_vol <- par[3]
  xi <- par[4]
  
  sigma_t <- exp(beta_0 + beta_vol * data$z_vol)
  
  if (any(!is.finite(sigma_t)) || any(sigma_t <= 0)) {
    return(1e10)
  }
  
  exceed_indicator <- data$loss > threshold
  
  x_exceed <- data$loss[exceed_indicator]
  sigma_exceed <- sigma_t[exceed_indicator]
  
  if (abs(xi) < 1e-6) {
    
    tail_intensity_all <- exp(-((threshold - mu) / sigma_t))
    
    integrated_intensity <- sum(tail_intensity_all) / trading_days_per_year
    
    log_intensity <- -log(sigma_exceed) -
      ((x_exceed - mu) / sigma_exceed)
    
  } else {
    
    threshold_support_all <- 1 + xi * ((threshold - mu) / sigma_t)
    exceed_support <- 1 + xi * ((x_exceed - mu) / sigma_exceed)
    
    if (
      any(!is.finite(threshold_support_all)) ||
      any(threshold_support_all <= 0) ||
      any(!is.finite(exceed_support)) ||
      any(exceed_support <= 0)
    ) {
      return(1e10)
    }
    
    tail_intensity_all <- threshold_support_all^(-1 / xi)
    
    integrated_intensity <- sum(tail_intensity_all) / trading_days_per_year
    
    log_intensity <- -log(sigma_exceed) -
      (1 / xi + 1) * log(exceed_support)
  }
  
  nll <- integrated_intensity - sum(log_intensity)
  
  if (!is.finite(nll)) {
    return(1e10)
  }
  
  return(nll)
}


# ------------------------------------------------------------
# Volatility-dependent location point process likelihood
# ------------------------------------------------------------
#
# Model:
#
#   mu_t = mu_0 + mu_vol * z_t
#   sigma constant
#   xi constant
#
# This tests whether volatility shifts the level/location of the
# extreme-loss process.

pp_vol_location_nll <- function(par, data, threshold, trading_days_per_year = 252) {
  
  mu_0 <- par[1]
  mu_vol <- par[2]
  log_sigma <- par[3]
  xi <- par[4]
  
  mu_t <- mu_0 + mu_vol * data$z_vol
  sigma <- exp(log_sigma)
  
  if (!is.finite(sigma) || sigma <= 0) {
    return(1e10)
  }
  
  exceed_indicator <- data$loss > threshold
  
  x_exceed <- data$loss[exceed_indicator]
  mu_exceed <- mu_t[exceed_indicator]
  
  if (abs(xi) < 1e-6) {
    
    tail_intensity_all <- exp(-((threshold - mu_t) / sigma))
    
    integrated_intensity <- sum(tail_intensity_all) / trading_days_per_year
    
    log_intensity <- -log(sigma) -
      ((x_exceed - mu_exceed) / sigma)
    
  } else {
    
    threshold_support_all <- 1 + xi * ((threshold - mu_t) / sigma)
    exceed_support <- 1 + xi * ((x_exceed - mu_exceed) / sigma)
    
    if (
      any(!is.finite(threshold_support_all)) ||
      any(threshold_support_all <= 0) ||
      any(!is.finite(exceed_support)) ||
      any(exceed_support <= 0)
    ) {
      return(1e10)
    }
    
    tail_intensity_all <- threshold_support_all^(-1 / xi)
    
    integrated_intensity <- sum(tail_intensity_all) / trading_days_per_year
    
    log_intensity <- -log(sigma) -
      (1 / xi + 1) * log(exceed_support)
  }
  
  nll <- integrated_intensity - sum(log_intensity)
  
  if (!is.finite(nll)) {
    return(1e10)
  }
  
  return(nll)
}


# ------------------------------------------------------------
# Volatility-dependent location + scale point process likelihood
# ------------------------------------------------------------
#
# Model:
#
#   mu_t = mu_0 + mu_vol * z_t
#   log(sigma_t) = beta_0 + beta_vol * z_t
#   xi constant
#
# This is the most flexible volatility point process model here.
# It can separate shifts in the level of extremes from changes in
# the spread of extremes.

pp_vol_location_scale_nll <- function(par, data, threshold, trading_days_per_year = 252) {
  
  mu_0 <- par[1]
  mu_vol <- par[2]
  beta_0 <- par[3]
  beta_vol <- par[4]
  xi <- par[5]
  
  mu_t <- mu_0 + mu_vol * data$z_vol
  sigma_t <- exp(beta_0 + beta_vol * data$z_vol)
  
  if (any(!is.finite(sigma_t)) || any(sigma_t <= 0)) {
    return(1e10)
  }
  
  exceed_indicator <- data$loss > threshold
  
  x_exceed <- data$loss[exceed_indicator]
  mu_exceed <- mu_t[exceed_indicator]
  sigma_exceed <- sigma_t[exceed_indicator]
  
  if (abs(xi) < 1e-6) {
    
    tail_intensity_all <- exp(-((threshold - mu_t) / sigma_t))
    
    integrated_intensity <- sum(tail_intensity_all) / trading_days_per_year
    
    log_intensity <- -log(sigma_exceed) -
      ((x_exceed - mu_exceed) / sigma_exceed)
    
  } else {
    
    threshold_support_all <- 1 + xi * ((threshold - mu_t) / sigma_t)
    exceed_support <- 1 + xi * ((x_exceed - mu_exceed) / sigma_exceed)
    
    if (
      any(!is.finite(threshold_support_all)) ||
      any(threshold_support_all <= 0) ||
      any(!is.finite(exceed_support)) ||
      any(exceed_support <= 0)
    ) {
      return(1e10)
    }
    
    tail_intensity_all <- threshold_support_all^(-1 / xi)
    
    integrated_intensity <- sum(tail_intensity_all) / trading_days_per_year
    
    log_intensity <- -log(sigma_exceed) -
      (1 / xi + 1) * log(exceed_support)
  }
  
  nll <- integrated_intensity - sum(log_intensity)
  
  if (!is.finite(nll)) {
    return(1e10)
  }
  
  return(nll)
}


# ------------------------------------------------------------
# Standard error helper
# ------------------------------------------------------------

get_standard_errors <- function(hessian_matrix) {
  
  vcov_matrix <- tryCatch(
    solve(hessian_matrix),
    error = function(e) NULL
  )
  
  if (is.null(vcov_matrix)) {
    return(rep(NA, nrow(hessian_matrix)))
  }
  
  diagonal_values <- diag(vcov_matrix)
  diagonal_values[diagonal_values < 0] <- NA
  
  standard_errors <- sqrt(diagonal_values)
  
  return(standard_errors)
}


# ============================================================
# 9. Initial Values from GEV and GPD Fits
# ============================================================

# To help numerical optimization, we use earlier EVT models as
# starting values.
#
# GEV monthly maxima provides starting values for mu, sigma, xi.
# GPD threshold modeling provides a check on xi.

monthly_max_losses <- apply.monthly(
  xts(df$loss, order.by = df$date),
  max
)

monthly_max_losses_num <- as.numeric(monthly_max_losses)

capture.output(
  gev_monthly <- gev.fit(monthly_max_losses_num)
)

capture.output(
  gpd_threshold <- gpd.fit(df$loss, threshold = u, npy = trading_days_per_year)
)

gev_mu_start <- gev_monthly$mle[1]
gev_sigma_start <- gev_monthly$mle[2]
gev_xi_start <- gev_monthly$mle[3]

gpd_sigma_start <- gpd_threshold$mle[1]
gpd_xi_start <- gpd_threshold$mle[2]

cat("\nInitial values:\n")
cat("GEV monthly mu:", gev_mu_start, "\n")
cat("GEV monthly sigma:", gev_sigma_start, "\n")
cat("GEV monthly xi:", gev_xi_start, "\n")
cat("GPD sigma:", gpd_sigma_start, "\n")
cat("GPD xi:", gpd_xi_start, "\n")

initial_values <- data.frame(
  source = c(
    "GEV monthly maxima",
    "GEV monthly maxima",
    "GEV monthly maxima",
    "GPD threshold",
    "GPD threshold"
  ),
  parameter = c("mu", "sigma", "xi", "sigma", "xi"),
  estimate = c(
    gev_mu_start,
    gev_sigma_start,
    gev_xi_start,
    gpd_sigma_start,
    gpd_xi_start
  )
)

write.csv(
  initial_values,
  file.path(res_dir, "initial_values_from_gev_gpd.csv"),
  row.names = FALSE
)


# ============================================================
# 10. Fit Point Process Models
# ============================================================

# ------------------------------------------------------------
# 10.1 Stationary point process
# ------------------------------------------------------------

start_stationary_pp <- c(
  gev_mu_start,
  log(gev_sigma_start),
  gev_xi_start
)

fit_stationary_pp <- optim(
  par = start_stationary_pp,
  fn = pp_stationary_nll,
  data = df,
  threshold = u,
  trading_days_per_year = trading_days_per_year,
  method = "BFGS",
  hessian = TRUE,
  control = list(maxit = 10000)
)

stationary_pp_se <- get_standard_errors(fit_stationary_pp$hessian)

stationary_pp_estimates <- data.frame(
  model = "stationary point process",
  parameter = c("mu", "log_sigma", "xi"),
  estimate = fit_stationary_pp$par,
  standard_error = stationary_pp_se
)

stationary_pp_estimates$natural_scale <- NA
stationary_pp_estimates$natural_scale[
  stationary_pp_estimates$parameter == "log_sigma"
] <- exp(stationary_pp_estimates$estimate[
  stationary_pp_estimates$parameter == "log_sigma"
])

cat("\nStationary point process fit:\n")
print(stationary_pp_estimates)

write.csv(
  stationary_pp_estimates,
  file.path(res_dir, "stationary_point_process_parameter_estimates.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 10.2 Volatility-dependent scale point process
# ------------------------------------------------------------

start_vol_scale_pp <- c(
  fit_stationary_pp$par[1],
  fit_stationary_pp$par[2],
  0,
  fit_stationary_pp$par[3]
)

fit_vol_scale_pp <- optim(
  par = start_vol_scale_pp,
  fn = pp_vol_scale_nll,
  data = df,
  threshold = u,
  trading_days_per_year = trading_days_per_year,
  method = "BFGS",
  hessian = TRUE,
  control = list(maxit = 10000)
)

vol_scale_pp_se <- get_standard_errors(fit_vol_scale_pp$hessian)

vol_scale_pp_estimates <- data.frame(
  model = "volatility-dependent point process scale",
  parameter = c("mu", "beta_0", "beta_vol", "xi"),
  estimate = fit_vol_scale_pp$par,
  standard_error = vol_scale_pp_se
)

cat("\nVolatility-dependent scale point process fit:\n")
print(vol_scale_pp_estimates)

write.csv(
  vol_scale_pp_estimates,
  file.path(res_dir, "volatility_scale_point_process_parameter_estimates.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 10.3 Volatility-dependent location point process
# ------------------------------------------------------------

start_vol_location_pp <- c(
  fit_stationary_pp$par[1],
  0,
  fit_stationary_pp$par[2],
  fit_stationary_pp$par[3]
)

fit_vol_location_pp <- optim(
  par = start_vol_location_pp,
  fn = pp_vol_location_nll,
  data = df,
  threshold = u,
  trading_days_per_year = trading_days_per_year,
  method = "BFGS",
  hessian = TRUE,
  control = list(maxit = 10000)
)

vol_location_pp_se <- get_standard_errors(fit_vol_location_pp$hessian)

vol_location_pp_estimates <- data.frame(
  model = "volatility-dependent point process location",
  parameter = c("mu_0", "mu_vol", "log_sigma", "xi"),
  estimate = fit_vol_location_pp$par,
  standard_error = vol_location_pp_se
)

vol_location_pp_estimates$natural_scale <- NA
vol_location_pp_estimates$natural_scale[
  vol_location_pp_estimates$parameter == "log_sigma"
] <- exp(vol_location_pp_estimates$estimate[
  vol_location_pp_estimates$parameter == "log_sigma"
])

cat("\nVolatility-dependent location point process fit:\n")
print(vol_location_pp_estimates)

write.csv(
  vol_location_pp_estimates,
  file.path(res_dir, "volatility_location_point_process_parameter_estimates.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 10.4 Volatility-dependent location + scale point process
# ------------------------------------------------------------

start_vol_location_scale_pp <- c(
  fit_stationary_pp$par[1],
  0,
  fit_stationary_pp$par[2],
  0,
  fit_stationary_pp$par[3]
)

fit_vol_location_scale_pp <- optim(
  par = start_vol_location_scale_pp,
  fn = pp_vol_location_scale_nll,
  data = df,
  threshold = u,
  trading_days_per_year = trading_days_per_year,
  method = "BFGS",
  hessian = TRUE,
  control = list(maxit = 10000)
)

vol_location_scale_pp_se <- get_standard_errors(fit_vol_location_scale_pp$hessian)

vol_location_scale_pp_estimates <- data.frame(
  model = "volatility-dependent point process location and scale",
  parameter = c("mu_0", "mu_vol", "beta_0", "beta_vol", "xi"),
  estimate = fit_vol_location_scale_pp$par,
  standard_error = vol_location_scale_pp_se
)

cat("\nVolatility-dependent location + scale point process fit:\n")
print(vol_location_scale_pp_estimates)

write.csv(
  vol_location_scale_pp_estimates,
  file.path(res_dir, "volatility_location_scale_point_process_parameter_estimates.csv"),
  row.names = FALSE
)


# ============================================================
# 11. Point Process Model Comparison
# ============================================================

nll_stationary_pp <- fit_stationary_pp$value
nll_vol_scale_pp <- fit_vol_scale_pp$value
nll_vol_location_pp <- fit_vol_location_pp$value
nll_vol_location_scale_pp <- fit_vol_location_scale_pp$value

k_stationary_pp <- length(fit_stationary_pp$par)
k_vol_scale_pp <- length(fit_vol_scale_pp$par)
k_vol_location_pp <- length(fit_vol_location_pp$par)
k_vol_location_scale_pp <- length(fit_vol_location_scale_pp$par)

pp_model_comparison <- data.frame(
  model = c(
    "stationary point process",
    "volatility-dependent scale point process",
    "volatility-dependent location point process",
    "volatility-dependent location + scale point process"
  ),
  n_parameters = c(
    k_stationary_pp,
    k_vol_scale_pp,
    k_vol_location_pp,
    k_vol_location_scale_pp
  ),
  n_observations = rep(n_exceedances, 4),
  nll = c(
    nll_stationary_pp,
    nll_vol_scale_pp,
    nll_vol_location_pp,
    nll_vol_location_scale_pp
  ),
  AIC = c(
    2 * k_stationary_pp + 2 * nll_stationary_pp,
    2 * k_vol_scale_pp + 2 * nll_vol_scale_pp,
    2 * k_vol_location_pp + 2 * nll_vol_location_pp,
    2 * k_vol_location_scale_pp + 2 * nll_vol_location_scale_pp
  ),
  BIC = c(
    log(n_exceedances) * k_stationary_pp + 2 * nll_stationary_pp,
    log(n_exceedances) * k_vol_scale_pp + 2 * nll_vol_scale_pp,
    log(n_exceedances) * k_vol_location_pp + 2 * nll_vol_location_pp,
    log(n_exceedances) * k_vol_location_scale_pp + 2 * nll_vol_location_scale_pp
  ),
  xi = c(
    fit_stationary_pp$par[3],
    fit_vol_scale_pp$par[4],
    fit_vol_location_pp$par[4],
    fit_vol_location_scale_pp$par[5]
  ),
  convergence = c(
    fit_stationary_pp$convergence,
    fit_vol_scale_pp$convergence,
    fit_vol_location_pp$convergence,
    fit_vol_location_scale_pp$convergence
  )
)

cat("\nPoint process model comparison:\n")
print(pp_model_comparison)

write.csv(
  pp_model_comparison,
  file.path(res_dir, "point_process_model_comparison.csv"),
  row.names = FALSE
)


# ============================================================
# 12. Likelihood Ratio Tests for Nested PP Models
# ============================================================

pp_lr_tests <- data.frame(
  comparison = c(
    "stationary vs volatility scale",
    "stationary vs volatility location",
    "volatility scale vs volatility location + scale",
    "volatility location vs volatility location + scale"
  ),
  LR_statistic = c(
    2 * (nll_stationary_pp - nll_vol_scale_pp),
    2 * (nll_stationary_pp - nll_vol_location_pp),
    2 * (nll_vol_scale_pp - nll_vol_location_scale_pp),
    2 * (nll_vol_location_pp - nll_vol_location_scale_pp)
  ),
  df = c(1, 1, 1, 1)
)

pp_lr_tests$p_value <- pchisq(
  pp_lr_tests$LR_statistic,
  df = pp_lr_tests$df,
  lower.tail = FALSE
)

cat("\nPoint process likelihood ratio tests:\n")
print(pp_lr_tests)

write.csv(
  pp_lr_tests,
  file.path(res_dir, "point_process_likelihood_ratio_tests.csv"),
  row.names = FALSE
)


# ============================================================
# 13. Correct Return-Level Functions
# ============================================================

# ------------------------------------------------------------
# Point process return levels
# ------------------------------------------------------------
#
# In the point process model, the annual tail measure is:
#
#   Lambda(z) = [1 + xi * ((z - mu) / sigma)]^(-1 / xi)
#
# A T-year return level z_T is defined by:
#
#   Lambda(z_T) = 1 / T
#
# Therefore:
#
#   z_T = mu + sigma / xi * (T^xi - 1)
#
# for xi != 0.
#
# For xi = 0:
#
#   z_T = mu + sigma * log(T)

pp_return_level <- function(T_years, mu, sigma, xi) {
  
  if (abs(xi) < 1e-6) {
    z_T <- mu + sigma * log(T_years)
  } else {
    z_T <- mu + (sigma / xi) * (T_years^xi - 1)
  }
  
  return(z_T)
}


# ------------------------------------------------------------
# Monthly GEV return levels
# ------------------------------------------------------------
#
# Chapter 3 used monthly block maxima.
#
# For a T-year return level:
#
#   m = 12 * T
#
# monthly blocks are involved.

gev_monthly_return_level <- function(T_years, mu, sigma, xi, blocks_per_year = 12) {
  
  m_blocks <- T_years * blocks_per_year
  p <- 1 - 1 / m_blocks
  
  if (abs(xi) < 1e-6) {
    z_T <- mu - sigma * log(-log(p))
  } else {
    z_T <- mu + (sigma / xi) * ((-log(p))^(-xi) - 1)
  }
  
  return(z_T)
}


# ------------------------------------------------------------
# GPD return levels
# ------------------------------------------------------------
#
# For the Chapter 4 threshold model:
#
#   x_m = u + sigma / xi * { (m * zeta)^xi - 1 }
#
# where:
#
#   m     = number of trading days in the return period
#   zeta  = P(L > u)

gpd_return_level <- function(m_days, threshold, sigma, xi, zeta) {
  
  if (m_days * zeta <= 1) {
    return(NA)
  }
  
  if (abs(xi) < 1e-6) {
    x_m <- threshold + sigma * log(m_days * zeta)
  } else {
    x_m <- threshold + (sigma / xi) * ((m_days * zeta)^xi - 1)
  }
  
  return(x_m)
}


# ------------------------------------------------------------
# Conditional return levels for different PP models
# ------------------------------------------------------------

pp_scale_conditional_return_level <- function(T_years, z_vol, par) {
  
  mu <- par[1]
  beta_0 <- par[2]
  beta_vol <- par[3]
  xi <- par[4]
  
  sigma_z <- exp(beta_0 + beta_vol * z_vol)
  
  return(pp_return_level(T_years, mu, sigma_z, xi))
}


pp_location_conditional_return_level <- function(T_years, z_vol, par) {
  
  mu_0 <- par[1]
  mu_vol <- par[2]
  log_sigma <- par[3]
  xi <- par[4]
  
  mu_z <- mu_0 + mu_vol * z_vol
  sigma <- exp(log_sigma)
  
  return(pp_return_level(T_years, mu_z, sigma, xi))
}


pp_location_scale_conditional_return_level <- function(T_years, z_vol, par) {
  
  mu_0 <- par[1]
  mu_vol <- par[2]
  beta_0 <- par[3]
  beta_vol <- par[4]
  xi <- par[5]
  
  mu_z <- mu_0 + mu_vol * z_vol
  sigma_z <- exp(beta_0 + beta_vol * z_vol)
  
  return(pp_return_level(T_years, mu_z, sigma_z, xi))
}


# ============================================================
# 14. Stationary Point Process Return Levels
# ============================================================

return_periods_years <- c(1, 5, 10, 20, 50, 100)

stationary_mu <- fit_stationary_pp$par[1]
stationary_sigma <- exp(fit_stationary_pp$par[2])
stationary_xi <- fit_stationary_pp$par[3]

stationary_pp_return_levels <- data.frame(
  model = "stationary point process",
  return_period_years = return_periods_years,
  return_level_decimal = sapply(
    return_periods_years,
    pp_return_level,
    mu = stationary_mu,
    sigma = stationary_sigma,
    xi = stationary_xi
  )
)

stationary_pp_return_levels$return_level_percent <-
  100 * stationary_pp_return_levels$return_level_decimal

cat("\nStationary point process return levels:\n")
print(stationary_pp_return_levels)

write.csv(
  stationary_pp_return_levels,
  file.path(res_dir, "stationary_point_process_return_levels.csv"),
  row.names = FALSE
)


# ============================================================
# 15. Volatility Scenarios
# ============================================================

# Same volatility states as Chapter 6:
#
#   25th percentile = low volatility
#   50th percentile = median volatility
#   75th percentile = high volatility
#   90th percentile = crisis volatility
#   95th percentile = extreme crisis volatility

vol_scenario_probs <- c(0.25, 0.50, 0.75, 0.90, 0.95)

vol_scenarios <- data.frame(
  volatility_state = c(
    "low volatility",
    "median volatility",
    "high volatility",
    "crisis volatility",
    "extreme crisis volatility"
  ),
  volatility_quantile = vol_scenario_probs,
  z_vol = as.numeric(quantile(df$z_vol, vol_scenario_probs, names = FALSE))
)

write.csv(
  vol_scenarios,
  file.path(res_dir, "volatility_scenarios.csv"),
  row.names = FALSE
)


# ============================================================
# 16. Conditional Return Levels: Volatility Scale Model
# ============================================================

vol_scale_scenarios <- vol_scenarios

vol_scale_scenarios$mu_hat <- fit_vol_scale_pp$par[1]
vol_scale_scenarios$sigma_hat <- exp(
  fit_vol_scale_pp$par[2] +
    fit_vol_scale_pp$par[3] * vol_scale_scenarios$z_vol
)
vol_scale_scenarios$xi_hat <- fit_vol_scale_pp$par[4]

vol_scale_conditional_return_levels <- do.call(
  rbind,
  lapply(seq_len(nrow(vol_scale_scenarios)), function(i) {
    
    levels <- sapply(
      return_periods_years,
      pp_scale_conditional_return_level,
      z_vol = vol_scale_scenarios$z_vol[i],
      par = fit_vol_scale_pp$par
    )
    
    data.frame(
      model = "volatility-dependent scale point process",
      volatility_state = vol_scale_scenarios$volatility_state[i],
      volatility_quantile = vol_scale_scenarios$volatility_quantile[i],
      z_vol = vol_scale_scenarios$z_vol[i],
      mu_hat = vol_scale_scenarios$mu_hat[i],
      sigma_hat = vol_scale_scenarios$sigma_hat[i],
      xi_hat = vol_scale_scenarios$xi_hat[i],
      return_period_years = return_periods_years,
      return_level_decimal = as.numeric(levels),
      return_level_percent = 100 * as.numeric(levels)
    )
  })
)

write.csv(
  vol_scale_scenarios,
  file.path(res_dir, "volatility_scale_scenarios.csv"),
  row.names = FALSE
)

write.csv(
  vol_scale_conditional_return_levels,
  file.path(res_dir, "volatility_scale_conditional_return_levels.csv"),
  row.names = FALSE
)


# ============================================================
# 17. Conditional Return Levels: Volatility Location Model
# ============================================================

vol_location_scenarios <- vol_scenarios

vol_location_scenarios$mu_hat <- 
  fit_vol_location_pp$par[1] +
  fit_vol_location_pp$par[2] * vol_location_scenarios$z_vol

vol_location_scenarios$sigma_hat <- exp(fit_vol_location_pp$par[3])
vol_location_scenarios$xi_hat <- fit_vol_location_pp$par[4]

vol_location_conditional_return_levels <- do.call(
  rbind,
  lapply(seq_len(nrow(vol_location_scenarios)), function(i) {
    
    levels <- sapply(
      return_periods_years,
      pp_location_conditional_return_level,
      z_vol = vol_location_scenarios$z_vol[i],
      par = fit_vol_location_pp$par
    )
    
    data.frame(
      model = "volatility-dependent location point process",
      volatility_state = vol_location_scenarios$volatility_state[i],
      volatility_quantile = vol_location_scenarios$volatility_quantile[i],
      z_vol = vol_location_scenarios$z_vol[i],
      mu_hat = vol_location_scenarios$mu_hat[i],
      sigma_hat = vol_location_scenarios$sigma_hat[i],
      xi_hat = vol_location_scenarios$xi_hat[i],
      return_period_years = return_periods_years,
      return_level_decimal = as.numeric(levels),
      return_level_percent = 100 * as.numeric(levels)
    )
  })
)

write.csv(
  vol_location_scenarios,
  file.path(res_dir, "volatility_location_scenarios.csv"),
  row.names = FALSE
)

write.csv(
  vol_location_conditional_return_levels,
  file.path(res_dir, "volatility_location_conditional_return_levels.csv"),
  row.names = FALSE
)


# ============================================================
# 18. Conditional Return Levels: Volatility Location + Scale
# ============================================================

vol_location_scale_scenarios <- vol_scenarios

vol_location_scale_scenarios$mu_hat <- 
  fit_vol_location_scale_pp$par[1] +
  fit_vol_location_scale_pp$par[2] * vol_location_scale_scenarios$z_vol

vol_location_scale_scenarios$sigma_hat <- exp(
  fit_vol_location_scale_pp$par[3] +
    fit_vol_location_scale_pp$par[4] * vol_location_scale_scenarios$z_vol
)

vol_location_scale_scenarios$xi_hat <- fit_vol_location_scale_pp$par[5]

vol_location_scale_conditional_return_levels <- do.call(
  rbind,
  lapply(seq_len(nrow(vol_location_scale_scenarios)), function(i) {
    
    levels <- sapply(
      return_periods_years,
      pp_location_scale_conditional_return_level,
      z_vol = vol_location_scale_scenarios$z_vol[i],
      par = fit_vol_location_scale_pp$par
    )
    
    data.frame(
      model = "volatility-dependent location + scale point process",
      volatility_state = vol_location_scale_scenarios$volatility_state[i],
      volatility_quantile = vol_location_scale_scenarios$volatility_quantile[i],
      z_vol = vol_location_scale_scenarios$z_vol[i],
      mu_hat = vol_location_scale_scenarios$mu_hat[i],
      sigma_hat = vol_location_scale_scenarios$sigma_hat[i],
      xi_hat = vol_location_scale_scenarios$xi_hat[i],
      return_period_years = return_periods_years,
      return_level_decimal = as.numeric(levels),
      return_level_percent = 100 * as.numeric(levels)
    )
  })
)

write.csv(
  vol_location_scale_scenarios,
  file.path(res_dir, "volatility_location_scale_scenarios.csv"),
  row.names = FALSE
)

write.csv(
  vol_location_scale_conditional_return_levels,
  file.path(res_dir, "volatility_location_scale_conditional_return_levels.csv"),
  row.names = FALSE
)


# ============================================================
# 19. Combined Conditional Return Level Table
# ============================================================

combined_pp_conditional_return_levels <- rbind(
  vol_scale_conditional_return_levels,
  vol_location_conditional_return_levels,
  vol_location_scale_conditional_return_levels
)

write.csv(
  combined_pp_conditional_return_levels,
  file.path(res_dir, "combined_volatility_pp_conditional_return_levels.csv"),
  row.names = FALSE
)


# ============================================================
# 20. Return Level Plots
# ============================================================

# ------------------------------------------------------------
# 20.1 Stationary PP return levels
# ------------------------------------------------------------

png(file.path(fig_dir, "stationary_point_process_return_levels.png"), width = 900, height = 600)

plot(
  stationary_pp_return_levels$return_period_years,
  stationary_pp_return_levels$return_level_percent,
  type = "b",
  pch = 19,
  col = "black",
  xlab = "Return Period in Years",
  ylab = "Return Level (%)",
  main = "Stationary Point Process Return Levels"
)

dev.off()


# ------------------------------------------------------------
# 20.2 Scale-only volatility PP return levels
# ------------------------------------------------------------

plot_colors <- c("darkgreen", "blue", "orange", "red", "darkred")
plot_pch <- c(19, 17, 15, 18, 8)

png(file.path(fig_dir, "volatility_scale_point_process_return_levels.png"), width = 1000, height = 650)

plot(
  return_periods_years,
  vol_scale_conditional_return_levels$return_level_percent[
    vol_scale_conditional_return_levels$volatility_state == "low volatility"
  ],
  type = "b",
  pch = plot_pch[1],
  col = plot_colors[1],
  ylim = range(vol_scale_conditional_return_levels$return_level_percent, na.rm = TRUE),
  xlab = "Return Period in Years",
  ylab = "Return Level (%)",
  main = "Scale-Only Volatility Point Process Return Levels"
)

for (i in 2:nrow(vol_scenarios)) {
  
  state <- vol_scenarios$volatility_state[i]
  
  lines(
    return_periods_years,
    vol_scale_conditional_return_levels$return_level_percent[
      vol_scale_conditional_return_levels$volatility_state == state
    ],
    type = "b",
    pch = plot_pch[i],
    col = plot_colors[i]
  )
}

legend(
  "topleft",
  legend = vol_scenarios$volatility_state,
  col = plot_colors,
  pch = plot_pch,
  lty = 1,
  bty = "n"
)

dev.off()


# ------------------------------------------------------------
# 20.3 Location-only volatility PP return levels
# ------------------------------------------------------------

png(file.path(fig_dir, "volatility_location_point_process_return_levels.png"), width = 1000, height = 650)

plot(
  return_periods_years,
  vol_location_conditional_return_levels$return_level_percent[
    vol_location_conditional_return_levels$volatility_state == "low volatility"
  ],
  type = "b",
  pch = plot_pch[1],
  col = plot_colors[1],
  ylim = range(vol_location_conditional_return_levels$return_level_percent, na.rm = TRUE),
  xlab = "Return Period in Years",
  ylab = "Return Level (%)",
  main = "Location-Only Volatility Point Process Return Levels"
)

for (i in 2:nrow(vol_scenarios)) {
  
  state <- vol_scenarios$volatility_state[i]
  
  lines(
    return_periods_years,
    vol_location_conditional_return_levels$return_level_percent[
      vol_location_conditional_return_levels$volatility_state == state
    ],
    type = "b",
    pch = plot_pch[i],
    col = plot_colors[i]
  )
}

legend(
  "topleft",
  legend = vol_scenarios$volatility_state,
  col = plot_colors,
  pch = plot_pch,
  lty = 1,
  bty = "n"
)

dev.off()


# ------------------------------------------------------------
# 20.4 Location + scale volatility PP return levels
# ------------------------------------------------------------

png(file.path(fig_dir, "volatility_location_scale_point_process_return_levels.png"), width = 1000, height = 650)

plot(
  return_periods_years,
  vol_location_scale_conditional_return_levels$return_level_percent[
    vol_location_scale_conditional_return_levels$volatility_state == "low volatility"
  ],
  type = "b",
  pch = plot_pch[1],
  col = plot_colors[1],
  ylim = range(vol_location_scale_conditional_return_levels$return_level_percent, na.rm = TRUE),
  xlab = "Return Period in Years",
  ylab = "Return Level (%)",
  main = "Location + Scale Volatility Point Process Return Levels"
)

for (i in 2:nrow(vol_scenarios)) {
  
  state <- vol_scenarios$volatility_state[i]
  
  lines(
    return_periods_years,
    vol_location_scale_conditional_return_levels$return_level_percent[
      vol_location_scale_conditional_return_levels$volatility_state == state
    ],
    type = "b",
    pch = plot_pch[i],
    col = plot_colors[i]
  )
}

legend(
  "topleft",
  legend = vol_scenarios$volatility_state,
  col = plot_colors,
  pch = plot_pch,
  lty = 1,
  bty = "n"
)

dev.off()


# ============================================================
# 21. Parameter vs Volatility Plots
# ============================================================

z_grid <- seq(min(df$z_vol), max(df$z_vol), length.out = 200)

# Scale-only sigma(z)

sigma_grid_scale <- exp(
  fit_vol_scale_pp$par[2] +
    fit_vol_scale_pp$par[3] * z_grid
)

png(file.path(fig_dir, "scale_only_pp_scale_vs_volatility.png"), width = 900, height = 600)

plot(
  z_grid,
  sigma_grid_scale,
  type = "l",
  lwd = 2,
  col = "blue",
  xlab = "Standardized Lagged Log Volatility",
  ylab = "Fitted Point Process Scale",
  main = "Scale-Only PP: Scale vs Volatility"
)

rug(df$z_vol[df$exceed], col = "red")
rug(df$z_vol[!df$exceed], col = "gray70")

legend(
  "topleft",
  legend = c("Fitted scale", "Exceedance volatility values"),
  col = c("blue", "red"),
  lwd = c(2, NA),
  pch = c(NA, "|"),
  bty = "n"
)

dev.off()


# Location-only mu(z)

mu_grid_location <- fit_vol_location_pp$par[1] +
  fit_vol_location_pp$par[2] * z_grid

png(file.path(fig_dir, "location_only_pp_location_vs_volatility.png"), width = 900, height = 600)

plot(
  z_grid,
  mu_grid_location,
  type = "l",
  lwd = 2,
  col = "blue",
  xlab = "Standardized Lagged Log Volatility",
  ylab = "Fitted Point Process Location",
  main = "Location-Only PP: Location vs Volatility"
)

rug(df$z_vol[df$exceed], col = "red")
rug(df$z_vol[!df$exceed], col = "gray70")

legend(
  "topleft",
  legend = c("Fitted location", "Exceedance volatility values"),
  col = c("blue", "red"),
  lwd = c(2, NA),
  pch = c(NA, "|"),
  bty = "n"
)

dev.off()


# Location + scale mu(z) and sigma(z)

mu_grid_location_scale <- fit_vol_location_scale_pp$par[1] +
  fit_vol_location_scale_pp$par[2] * z_grid

sigma_grid_location_scale <- exp(
  fit_vol_location_scale_pp$par[3] +
    fit_vol_location_scale_pp$par[4] * z_grid
)

png(file.path(fig_dir, "location_scale_pp_location_vs_volatility.png"), width = 900, height = 600)

plot(
  z_grid,
  mu_grid_location_scale,
  type = "l",
  lwd = 2,
  col = "blue",
  xlab = "Standardized Lagged Log Volatility",
  ylab = "Fitted Point Process Location",
  main = "Location + Scale PP: Location vs Volatility"
)

rug(df$z_vol[df$exceed], col = "red")
rug(df$z_vol[!df$exceed], col = "gray70")

dev.off()


png(file.path(fig_dir, "location_scale_pp_scale_vs_volatility.png"), width = 900, height = 600)

plot(
  z_grid,
  sigma_grid_location_scale,
  type = "l",
  lwd = 2,
  col = "blue",
  xlab = "Standardized Lagged Log Volatility",
  ylab = "Fitted Point Process Scale",
  main = "Location + Scale PP: Scale vs Volatility"
)

rug(df$z_vol[df$exceed], col = "red")
rug(df$z_vol[!df$exceed], col = "gray70")

dev.off()


# ============================================================
# 22. Stationary Point Process Diagnostics
# ============================================================

# Diagnostic logic:
#
# Under the fitted point process model, the conditional distribution
# of exceedances above u is:
#
#   P(X <= x | X > u) = 1 - Lambda(x) / Lambda(u)
#
# where Lambda is the annual tail measure.
#
# If the model fits well, the fitted CDF values should look Uniform(0,1),
# and:
#
#   -log(1 - fitted CDF)
#
# should look Exponential(1).

pp_tail_measure <- function(x, mu, sigma, xi) {
  
  if (abs(xi) < 1e-6) {
    tail <- exp(-((x - mu) / sigma))
  } else {
    support <- 1 + xi * ((x - mu) / sigma)
    
    tail <- rep(NA, length(x))
    tail[support > 0] <- support[support > 0]^(-1 / xi)
  }
  
  return(tail)
}

x_exceed <- exceed_df$loss

stationary_tail_u <- pp_tail_measure(
  u,
  stationary_mu,
  stationary_sigma,
  stationary_xi
)

stationary_tail_x <- pp_tail_measure(
  x_exceed,
  stationary_mu,
  stationary_sigma,
  stationary_xi
)

stationary_fitted_cdf <- 1 - stationary_tail_x / stationary_tail_u
stationary_fitted_cdf <- pmin(pmax(stationary_fitted_cdf, 1e-10), 1 - 1e-10)

stationary_exp_residuals <- -log(1 - stationary_fitted_cdf)

diagnostic_data_stationary <- data.frame(
  loss = x_exceed,
  fitted_cdf = stationary_fitted_cdf,
  exponential_residual = stationary_exp_residuals
)

write.csv(
  diagnostic_data_stationary,
  file.path(res_dir, "stationary_point_process_diagnostics.csv"),
  row.names = FALSE
)

png(file.path(fig_dir, "diagnostics_stationary_point_process.png"), width = 1000, height = 800)

par(mfrow = c(2, 2))

n_diag <- length(stationary_fitted_cdf)
pp <- ppoints(n_diag)

plot(
  pp,
  sort(stationary_fitted_cdf),
  pch = 19,
  main = "Stationary PP Probability Plot",
  xlab = "Theoretical probability",
  ylab = "Empirical fitted probability"
)
abline(0, 1, col = "red", lwd = 2)

plot(
  qexp(pp),
  sort(stationary_exp_residuals),
  pch = 19,
  main = "Stationary PP Exponential QQ Plot",
  xlab = "Theoretical exponential quantiles",
  ylab = "Model residual quantiles"
)
abline(0, 1, col = "red", lwd = 2)

hist(
  stationary_exp_residuals,
  breaks = 30,
  probability = TRUE,
  main = "Stationary PP Residual Density",
  xlab = "Exponential residuals"
)
curve(dexp(x), add = TRUE, col = "red", lwd = 2)

plot(
  stationary_exp_residuals,
  type = "h",
  main = "Stationary PP Residuals by Exceedance Order",
  xlab = "Exceedance index",
  ylab = "Exponential residual"
)

par(mfrow = c(1, 1))

dev.off()


# ============================================================
# 23. Shape Comparison with Earlier Chapters
# ============================================================

chapter_shape_comparison <- data.frame(
  model = c(
    "Chapter 3 GEV monthly maxima",
    "Chapter 4 stationary GPD",
    "Chapter 7 stationary point process",
    "Chapter 7 volatility scale point process",
    "Chapter 7 volatility location point process",
    "Chapter 7 volatility location + scale point process"
  ),
  xi = c(
    gev_monthly$mle[3],
    gpd_threshold$mle[2],
    fit_stationary_pp$par[3],
    fit_vol_scale_pp$par[4],
    fit_vol_location_pp$par[4],
    fit_vol_location_scale_pp$par[5]
  ),
  se_xi = c(
    gev_monthly$se[3],
    gpd_threshold$se[2],
    stationary_pp_se[3],
    vol_scale_pp_se[4],
    vol_location_pp_se[4],
    vol_location_scale_pp_se[5]
  )
)

cat("\nShape comparison:\n")
print(chapter_shape_comparison)

write.csv(
  chapter_shape_comparison,
  file.path(res_dir, "chapter_shape_comparison.csv"),
  row.names = FALSE
)

png(file.path(fig_dir, "chapter_shape_comparison.png"), width = 1100, height = 650)

barplot(
  chapter_shape_comparison$xi,
  names.arg = c("Ch3", "Ch4", "Ch7 stat", "Scale", "Loc", "Loc+Scale"),
  ylab = expression(hat(xi)),
  main = "Shape Parameter Comparison Across EVT Models",
  ylim = c(
    min(0, min(chapter_shape_comparison$xi, na.rm = TRUE)),
    max(chapter_shape_comparison$xi + chapter_shape_comparison$se_xi, na.rm = TRUE) * 1.4
  )
)

abline(h = 0, lty = 2, col = "red")

dev.off()


# ============================================================
# 24. Return Level Comparison: GEV, GPD, and Stationary PP
# ============================================================

# Compare:
#
#   Chapter 3 monthly GEV
#   Chapter 4 stationary GPD
#   Chapter 7 stationary point process
#
# The Chapter 3 comparison correctly uses 12 monthly blocks per year.
# The Chapter 7 comparison uses the point process return-level formula.

comparison_return_periods_years <- c(1, 5, 10)
comparison_return_periods_days <- trading_days_per_year * comparison_return_periods_years

gev_return_levels <- data.frame(
  model = "Chapter 3 GEV monthly maxima",
  return_period_years = comparison_return_periods_years,
  return_level_decimal = sapply(
    comparison_return_periods_years,
    gev_monthly_return_level,
    mu = gev_monthly$mle[1],
    sigma = gev_monthly$mle[2],
    xi = gev_monthly$mle[3],
    blocks_per_year = 12
  )
)

gpd_return_levels <- data.frame(
  model = "Chapter 4 stationary GPD",
  return_period_years = comparison_return_periods_years,
  return_level_decimal = sapply(
    comparison_return_periods_days,
    gpd_return_level,
    threshold = u,
    sigma = gpd_threshold$mle[1],
    xi = gpd_threshold$mle[2],
    zeta = exceedance_fraction
  )
)

pp_return_levels_comparison <- stationary_pp_return_levels[
  stationary_pp_return_levels$return_period_years %in% comparison_return_periods_years,
  c("model", "return_period_years", "return_level_decimal")
]

return_level_comparison <- rbind(
  gev_return_levels,
  gpd_return_levels,
  pp_return_levels_comparison
)

return_level_comparison$return_level_percent <-
  100 * return_level_comparison$return_level_decimal

cat("\nReturn level comparison:\n")
print(return_level_comparison)

write.csv(
  return_level_comparison,
  file.path(res_dir, "return_level_comparison_gev_gpd_pp.csv"),
  row.names = FALSE
)

png(file.path(fig_dir, "return_level_comparison_gev_gpd_pp.png"), width = 1000, height = 650)

plot(
  comparison_return_periods_years,
  return_level_comparison$return_level_percent[
    return_level_comparison$model == "Chapter 3 GEV monthly maxima"
  ],
  type = "b",
  pch = 19,
  col = "black",
  ylim = range(return_level_comparison$return_level_percent, na.rm = TRUE),
  xlab = "Return Period in Years",
  ylab = "Return Level (%)",
  main = "Return Level Comparison: GEV, GPD, and Point Process"
)

lines(
  comparison_return_periods_years,
  return_level_comparison$return_level_percent[
    return_level_comparison$model == "Chapter 4 stationary GPD"
  ],
  type = "b",
  pch = 17,
  col = "blue"
)

lines(
  comparison_return_periods_years,
  return_level_comparison$return_level_percent[
    return_level_comparison$model == "stationary point process"
  ],
  type = "b",
  pch = 15,
  col = "red"
)

legend(
  "topleft",
  legend = c("Chapter 3 GEV", "Chapter 4 GPD", "Chapter 7 PP"),
  col = c("black", "blue", "red"),
  pch = c(19, 17, 15),
  lty = 1,
  bty = "n"
)

dev.off()


# ============================================================
# 25. Save Complete Chapter 7 Summary
# ============================================================

sink(file.path(res_dir, "chapter7_model_summary.txt"))

cat("Chapter 7: Point Process Models for Extremes\n")
cat("============================================\n\n")

cat("Data:\n")
cat("Asset: SPY\n")
cat("Returns: daily log returns\n")
cat("Losses: negative daily log returns\n")
cat("Covariate: lagged 21-day realized volatility\n")
cat("Usable observations:", n_total, "\n")
cat("Approximate trading years:", n_years, "\n\n")

cat("Threshold summary:\n")
print(threshold_summary)
cat("\n")

cat("Initial values from GEV and GPD fits:\n")
print(initial_values)
cat("\n")

cat("Stationary point process parameter estimates:\n")
print(stationary_pp_estimates)
cat("\n")

cat("Volatility scale point process parameter estimates:\n")
print(vol_scale_pp_estimates)
cat("\n")

cat("Volatility location point process parameter estimates:\n")
print(vol_location_pp_estimates)
cat("\n")

cat("Volatility location + scale point process parameter estimates:\n")
print(vol_location_scale_pp_estimates)
cat("\n")

cat("Point process model comparison:\n")
print(pp_model_comparison)
cat("\n")

cat("Point process likelihood ratio tests:\n")
print(pp_lr_tests)
cat("\n")

cat("Stationary point process return levels:\n")
print(stationary_pp_return_levels)
cat("\n")

cat("Volatility scenarios:\n")
print(vol_scenarios)
cat("\n")

cat("Volatility scale scenarios:\n")
print(vol_scale_scenarios)
cat("\n")

cat("Volatility scale conditional return levels:\n")
print(vol_scale_conditional_return_levels)
cat("\n")

cat("Volatility location scenarios:\n")
print(vol_location_scenarios)
cat("\n")

cat("Volatility location conditional return levels:\n")
print(vol_location_conditional_return_levels)
cat("\n")

cat("Volatility location + scale scenarios:\n")
print(vol_location_scale_scenarios)
cat("\n")

cat("Volatility location + scale conditional return levels:\n")
print(vol_location_scale_conditional_return_levels)
cat("\n")

cat("Shape comparison:\n")
print(chapter_shape_comparison)
cat("\n")

cat("Return level comparison: GEV, GPD, and point process:\n")
print(return_level_comparison)
cat("\n")

cat("Interpretation guide:\n")
cat("The stationary point process model uses exceedances above a high threshold while estimating GEV-compatible parameters.\n")
cat("It provides a unified framework connecting Chapter 3 block maxima and Chapter 4 threshold exceedances.\n")
cat("The point process return levels use Lambda(z_T) = 1 / T.\n")
cat("The Chapter 3 GEV return levels are adjusted for monthly blocks using 12 blocks per year.\n")
cat("Positive xi estimates support heavy-tailed extreme losses.\n")
cat("The volatility-dependent point process models test whether volatility affects scale, location, or both.\n")
cat("If allowing location to vary improves the model and gives positive mu_vol, volatility shifts the extreme-loss process upward.\n")
cat("The best volatility specification should be chosen using AIC, BIC, LR tests, convergence, and interpretability.\n")

sink()


sink(file.path(res_dir, "session_info.txt"))
print(sessionInfo())
sink()


# ============================================================
# 26. Final Console Summary
# ============================================================

cat("\n================ CHAPTER 7 SUMMARY ================\n")
cat("Point process EVT analysis completed.\n\n")

cat("Threshold:", round(100 * u, 3), "%\n")
cat("Number of exceedances:", n_exceedances, "\n")
cat("Approximate trading years:", round(n_years, 2), "\n\n")

cat("Stationary point process estimates:\n")
print(stationary_pp_estimates)

cat("\nVolatility scale PP estimates:\n")
print(vol_scale_pp_estimates)

cat("\nVolatility location PP estimates:\n")
print(vol_location_pp_estimates)

cat("\nVolatility location + scale PP estimates:\n")
print(vol_location_scale_pp_estimates)

cat("\nPoint process model comparison:\n")
print(pp_model_comparison)

cat("\nLikelihood ratio tests:\n")
print(pp_lr_tests)

cat("\nStationary point process return levels:\n")
print(stationary_pp_return_levels)

cat("\nVolatility scale conditional return levels:\n")
print(vol_scale_conditional_return_levels)

cat("\nVolatility location conditional return levels:\n")
print(vol_location_conditional_return_levels)

cat("\nVolatility location + scale conditional return levels:\n")
print(vol_location_scale_conditional_return_levels)

cat("\nShape comparison:\n")
print(chapter_shape_comparison)

cat("\nReturn level comparison:\n")
print(return_level_comparison)

cat("\nFigures saved in:\n")
cat(fig_dir, "\n")

cat("\nResults saved in:\n")
cat(res_dir, "\n")

cat("\nSaved figures:\n")
print(list.files(fig_dir))

cat("\nSaved result files:\n")
print(list.files(res_dir))

cat("===================================================\n")