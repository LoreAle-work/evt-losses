# ============================================================
# Chapter 6: Extremes of Non-Stationary Sequences
# Non-stationary GPD Models for SPY Daily Losses
# ============================================================
#
# This script extends the Chapter 4 and Chapter 5 threshold models
# by allowing the scale parameter of the GPD to vary with market volatility.
#
# Main idea:
#   Chapter 4: threshold exceedances with constant GPD parameters
#   Chapter 5: threshold exceedances with clustering / declustering
#   Chapter 6: threshold exceedances with covariate-dependent parameters
#
# This revised version also adds a Chapter 5-style robustness check:
#
#   all exceedances -> clusters -> cluster maxima -> non-stationary GPD
#
# Main models:
#   1. Stationary GPD on all exceedances
#   2. Volatility-dependent GPD on all exceedances
#   3. Volatility + time GPD on all exceedances
#   4. Stationary GPD on declustered cluster maxima
#   5. Volatility-dependent GPD on declustered cluster maxima
#
# ============================================================


# ============================================================
# 0. Package Setup
# ============================================================

# Run this only once if packages are not installed:
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

if (basename(getwd()) == "04_chapter6_nonstationary_sequences") {
  chapter_dir <- "."
} else {
  chapter_dir <- "04_chapter6_nonstationary_sequences"
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

getSymbols("SPY", src = "yahoo", from = "1993-01-01", auto.assign = TRUE)

prices <- Ad(SPY)

head(prices)
tail(prices)

png(file.path(fig_dir, "spy_adjusted_prices.png"), width = 1000, height = 600)
plot(prices, main = "SPY Adjusted Closing Prices")
dev.off()


# ============================================================
# 3. Compute Daily Log Returns
# ============================================================

# Daily log returns:
#
# R_t = log(P_t) - log(P_{t-1})

returns <- diff(log(prices))
returns <- na.omit(returns)
colnames(returns) <- "return"

ret_num <- as.numeric(returns)
dates <- index(returns)

head(returns)
summary(returns)

png(file.path(fig_dir, "spy_daily_log_returns.png"), width = 1000, height = 600)
plot(returns, main = "SPY Daily Log Returns")
dev.off()


# ============================================================
# 4. Convert Returns into Losses
# ============================================================

# Losses:
#
# L_t = -R_t

losses <- -returns
colnames(losses) <- "loss"

losses_num <- as.numeric(losses)
loss_dates <- index(losses)

head(losses)
summary(losses)

png(file.path(fig_dir, "spy_daily_losses.png"), width = 1000, height = 600)
plot(losses, main = "SPY Daily Losses")
dev.off()

cat("\nLargest daily losses (%):\n")
print(100 * sort(losses_num, decreasing = TRUE)[1:10])


# ============================================================
# 5. Build Volatility Covariate
# ============================================================

# We use lagged 21-day realized volatility.
#
# rolling_sd_21 estimates volatility using the previous 21 trading days.
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

df$log_vol <- log(df$rv21_lag)
df$z_vol <- as.numeric(scale(df$log_vol))
df$t_scaled <- as.numeric(scale(seq_len(nrow(df))))

head(df)
summary(df)

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
# 6. Threshold Choice
# ============================================================

# Use the 97.5% empirical quantile of daily losses.

u <- as.numeric(quantile(df$loss, 0.975, names = FALSE))

df$exceed <- df$loss > u
df$excess <- ifelse(df$exceed, df$loss - u, NA)

n_total <- nrow(df)
n_exceedances <- sum(df$exceed)
exceedance_rate <- mean(df$exceed)

cat("\nThreshold summary:\n")
cat("97.5% threshold =", round(100 * u, 3), "%\n")
cat("Total observations:", n_total, "\n")
cat("Number of exceedances:", n_exceedances, "\n")
cat("Exceedance rate:", exceedance_rate, "\n")

threshold_summary <- data.frame(
  threshold_probability = 0.975,
  threshold_value = u,
  threshold_percent = 100 * u,
  total_observations = n_total,
  number_exceedances = n_exceedances,
  exceedance_fraction = exceedance_rate
)

write.csv(
  threshold_summary,
  file.path(res_dir, "threshold_summary.csv"),
  row.names = FALSE
)

exceed_df <- df[df$exceed, ]

y <- exceed_df$excess
z_exceed <- exceed_df$z_vol
t_exceed <- exceed_df$t_scaled


# ============================================================
# 7. Plot Losses, Threshold, and Exceedances
# ============================================================

png(file.path(fig_dir, "daily_losses_with_threshold_and_exceedances.png"), width = 1100, height = 650)

plot(
  df$date,
  df$loss,
  type = "h",
  lwd = 0.5,
  col = "gray65",
  ylim = c(0, max(df$loss) * 1.05),
  main = "SPY Daily Losses with 97.5% Threshold",
  xlab = "Date",
  ylab = "Daily Loss"
)

abline(h = u, col = "red", lwd = 3, lty = 2)

points(
  exceed_df$date,
  exceed_df$loss,
  pch = 20,
  col = "red"
)

legend(
  "topright",
  legend = c(
    paste0("97.5% threshold = ", round(100 * u, 2), "%"),
    "Exceedances"
  ),
  col = c("red", "red"),
  lty = c(2, NA),
  lwd = c(3, NA),
  pch = c(NA, 20),
  bty = "n"
)

dev.off()


png(file.path(fig_dir, "volatility_with_exceedance_days.png"), width = 1100, height = 650)

plot(
  df$date,
  df$rv21_lag,
  type = "l",
  col = "black",
  main = "Lagged 21-Day Volatility with Exceedance Days",
  xlab = "Date",
  ylab = "Lagged 21-Day Volatility"
)

points(
  exceed_df$date,
  exceed_df$rv21_lag,
  pch = 20,
  col = "red"
)

legend(
  "topright",
  legend = c("Lagged volatility", "Exceedance days"),
  col = c("black", "red"),
  lty = c(1, NA),
  pch = c(NA, 20),
  bty = "n"
)

dev.off()


# ============================================================
# 8. Logistic Model for Exceedance Probability
# ============================================================

# This models:
#
# P(L_t > u | z_t)
#
# where z_t is standardized lagged log-volatility.

logit_const <- glm(exceed ~ 1, data = df, family = binomial)
logit_vol <- glm(exceed ~ z_vol, data = df, family = binomial)
logit_vol_time <- glm(exceed ~ z_vol + t_scaled, data = df, family = binomial)

logit_model_comparison <- data.frame(
  model = c(
    "constant exceedance probability",
    "volatility-dependent exceedance probability",
    "volatility + time exceedance probability"
  ),
  n_parameters = c(1, 2, 3),
  AIC = c(AIC(logit_const), AIC(logit_vol), AIC(logit_vol_time)),
  BIC = c(BIC(logit_const), BIC(logit_vol), BIC(logit_vol_time))
)

print(logit_model_comparison)

write.csv(
  logit_model_comparison,
  file.path(res_dir, "logistic_exceedance_model_comparison.csv"),
  row.names = FALSE
)

logit_vol_coefficients <- as.data.frame(coef(summary(logit_vol)))
logit_vol_coefficients$term <- row.names(logit_vol_coefficients)
row.names(logit_vol_coefficients) <- NULL

write.csv(
  logit_vol_coefficients,
  file.path(res_dir, "logistic_volatility_coefficients.csv"),
  row.names = FALSE
)

z_grid <- seq(min(df$z_vol), max(df$z_vol), length.out = 200)

logit_pred <- predict(
  logit_vol,
  newdata = data.frame(z_vol = z_grid),
  type = "response"
)

png(file.path(fig_dir, "exceedance_probability_vs_volatility.png"), width = 900, height = 600)

plot(
  z_grid,
  logit_pred,
  type = "l",
  lwd = 2,
  col = "blue",
  xlab = "Standardized Lagged Log Volatility",
  ylab = "Predicted Exceedance Probability",
  main = "Threshold Exceedance Probability vs Volatility"
)

rug(df$z_vol[df$exceed], col = "red")
rug(df$z_vol[!df$exceed], col = "gray70")

legend(
  "topleft",
  legend = c("Predicted probability", "Exceedance days"),
  col = c("blue", "red"),
  lwd = c(2, NA),
  pch = c(NA, "|"),
  bty = "n"
)

dev.off()


# ============================================================
# 9. Custom GPD Likelihood Functions
# ============================================================

gpd_negative_loglik <- function(par, y, X) {
  
  k <- ncol(X)
  
  beta <- par[1:k]
  xi <- par[k + 1]
  
  eta <- as.vector(X %*% beta)
  sigma <- exp(eta)
  
  if (any(!is.finite(sigma)) || any(sigma <= 0)) {
    return(1e10)
  }
  
  if (abs(xi) < 1e-6) {
    
    nll <- sum(log(sigma) + y / sigma)
    
  } else {
    
    support <- 1 + xi * y / sigma
    
    if (any(!is.finite(support)) || any(support <= 0)) {
      return(1e10)
    }
    
    nll <- sum(log(sigma) + (1 + 1 / xi) * log(support))
  }
  
  if (!is.finite(nll)) {
    return(1e10)
  }
  
  return(nll)
}


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


fit_gpd_model <- function(y, X, start_par, model_name) {
  
  fit <- optim(
    par = start_par,
    fn = gpd_negative_loglik,
    y = y,
    X = X,
    method = "BFGS",
    hessian = TRUE,
    control = list(maxit = 10000)
  )
  
  k <- ncol(X)
  beta_hat <- fit$par[1:k]
  xi_hat <- fit$par[k + 1]
  
  se <- get_standard_errors(fit$hessian)
  
  output <- list(
    model_name = model_name,
    par = fit$par,
    beta = beta_hat,
    xi = xi_hat,
    nll = fit$value,
    convergence = fit$convergence,
    hessian = fit$hessian,
    se = se,
    se_beta = se[1:k],
    se_xi = se[k + 1],
    X = X,
    y = y
  )
  
  return(output)
}


# ============================================================
# 10. Initial Values
# ============================================================

capture.output(
  initial_gpd <- gpd.fit(df$loss, threshold = u, npy = 252)
)

initial_sigma <- initial_gpd$mle[1]
initial_xi <- initial_gpd$mle[2]

cat("\nInitial GPD estimates from ismev:\n")
cat("sigma:", initial_sigma, "\n")
cat("xi:", initial_xi, "\n")


# ============================================================
# 11. Fit Stationary and Non-Stationary GPD Models
#     Using All Threshold Exceedances
# ============================================================

X_stationary <- matrix(1, nrow = length(y), ncol = 1)
colnames(X_stationary) <- "intercept"

X_vol <- cbind(1, z_exceed)
colnames(X_vol) <- c("intercept", "z_vol")

X_vol_time <- cbind(1, z_exceed, t_exceed)
colnames(X_vol_time) <- c("intercept", "z_vol", "t_scaled")

start_stationary <- c(log(initial_sigma), initial_xi)
start_vol <- c(log(initial_sigma), 0, initial_xi)
start_vol_time <- c(log(initial_sigma), 0, 0, initial_xi)

fit_stationary <- fit_gpd_model(
  y = y,
  X = X_stationary,
  start_par = start_stationary,
  model_name = "stationary GPD, all exceedances"
)

fit_vol <- fit_gpd_model(
  y = y,
  X = X_vol,
  start_par = start_vol,
  model_name = "volatility-dependent GPD, all exceedances"
)

fit_vol_time <- fit_gpd_model(
  y = y,
  X = X_vol_time,
  start_par = start_vol_time,
  model_name = "volatility + time GPD, all exceedances"
)

cat("\nConvergence codes, all exceedances:\n")
cat("Stationary:", fit_stationary$convergence, "\n")
cat("Volatility:", fit_vol$convergence, "\n")
cat("Volatility + time:", fit_vol_time$convergence, "\n")


# ============================================================
# 12. Model Comparison: All Exceedances
# ============================================================

n_excess <- length(y)

model_comparison <- data.frame(
  model = c(
    fit_stationary$model_name,
    fit_vol$model_name,
    fit_vol_time$model_name
  ),
  n_parameters = c(
    length(fit_stationary$par),
    length(fit_vol$par),
    length(fit_vol_time$par)
  ),
  n_observations = c(n_excess, n_excess, n_excess),
  nll = c(
    fit_stationary$nll,
    fit_vol$nll,
    fit_vol_time$nll
  ),
  AIC = c(
    2 * length(fit_stationary$par) + 2 * fit_stationary$nll,
    2 * length(fit_vol$par) + 2 * fit_vol$nll,
    2 * length(fit_vol_time$par) + 2 * fit_vol_time$nll
  ),
  BIC = c(
    log(n_excess) * length(fit_stationary$par) + 2 * fit_stationary$nll,
    log(n_excess) * length(fit_vol$par) + 2 * fit_vol$nll,
    log(n_excess) * length(fit_vol_time$par) + 2 * fit_vol_time$nll
  ),
  xi = c(
    fit_stationary$xi,
    fit_vol$xi,
    fit_vol_time$xi
  ),
  se_xi = c(
    fit_stationary$se_xi,
    fit_vol$se_xi,
    fit_vol_time$se_xi
  ),
  convergence = c(
    fit_stationary$convergence,
    fit_vol$convergence,
    fit_vol_time$convergence
  )
)

print(model_comparison)

write.csv(
  model_comparison,
  file.path(res_dir, "gpd_nonstationary_model_comparison.csv"),
  row.names = FALSE
)

lr_stationary_vs_vol <- 2 * (fit_stationary$nll - fit_vol$nll)
p_stationary_vs_vol <- pchisq(lr_stationary_vs_vol, df = 1, lower.tail = FALSE)

lr_vol_vs_vol_time <- 2 * (fit_vol$nll - fit_vol_time$nll)
p_vol_vs_vol_time <- pchisq(lr_vol_vs_vol_time, df = 1, lower.tail = FALSE)

lr_tests <- data.frame(
  comparison = c(
    "stationary vs volatility-dependent scale, all exceedances",
    "volatility-dependent scale vs volatility + time scale, all exceedances"
  ),
  LR_statistic = c(lr_stationary_vs_vol, lr_vol_vs_vol_time),
  df = c(1, 1),
  p_value = c(p_stationary_vs_vol, p_vol_vs_vol_time)
)

print(lr_tests)

write.csv(
  lr_tests,
  file.path(res_dir, "likelihood_ratio_tests.csv"),
  row.names = FALSE
)


# ============================================================
# 13. Parameter Tables: All Exceedances
# ============================================================

make_parameter_table <- function(fit, parameter_names) {
  
  data.frame(
    model = fit$model_name,
    parameter = parameter_names,
    estimate = fit$par,
    standard_error = fit$se
  )
}

parameter_table <- rbind(
  make_parameter_table(fit_stationary, c("beta_0", "xi")),
  make_parameter_table(fit_vol, c("beta_0", "beta_vol", "xi")),
  make_parameter_table(fit_vol_time, c("beta_0", "beta_vol", "beta_time", "xi"))
)

print(parameter_table)

write.csv(
  parameter_table,
  file.path(res_dir, "gpd_nonstationary_parameter_estimates.csv"),
  row.names = FALSE
)


# ============================================================
# 14. Fitted Scale as a Function of Volatility
# ============================================================

z_grid <- seq(min(df$z_vol), max(df$z_vol), length.out = 200)

beta_vol_model <- fit_vol$beta

sigma_grid <- exp(beta_vol_model[1] + beta_vol_model[2] * z_grid)

png(file.path(fig_dir, "gpd_scale_vs_volatility.png"), width = 900, height = 600)

plot(
  z_grid,
  sigma_grid,
  type = "l",
  lwd = 2,
  col = "blue",
  xlab = "Standardized Lagged Log Volatility",
  ylab = "Fitted GPD Scale",
  main = "Non-Stationary GPD Scale vs Volatility"
)

rug(z_exceed, col = "red")

legend(
  "topleft",
  legend = c("Fitted scale", "Exceedance volatility values"),
  col = c("blue", "red"),
  lwd = c(2, NA),
  pch = c(NA, "|"),
  bty = "n"
)

dev.off()

sigma_exceed_vol <- exp(as.vector(X_vol %*% fit_vol$beta))

exceed_df$fitted_sigma_vol_model <- sigma_exceed_vol

png(file.path(fig_dir, "fitted_scale_at_exceedances_over_time.png"), width = 1100, height = 650)

plot(
  exceed_df$date,
  exceed_df$fitted_sigma_vol_model,
  type = "b",
  pch = 19,
  col = "blue",
  main = "Fitted GPD Scale at Exceedance Dates",
  xlab = "Date",
  ylab = "Fitted Scale"
)

dev.off()

write.csv(
  exceed_df,
  file.path(res_dir, "exceedance_data_with_fitted_scale.csv"),
  row.names = FALSE
)


# ============================================================
# 15. Custom Diagnostic Functions
# ============================================================

pgpd_custom <- function(y, sigma, xi) {
  
  if (abs(xi) < 1e-6) {
    cdf <- 1 - exp(-y / sigma)
  } else {
    support <- 1 + xi * y / sigma
    cdf <- 1 - support^(-1 / xi)
  }
  
  cdf <- pmin(pmax(cdf, 1e-10), 1 - 1e-10)
  
  return(cdf)
}


save_custom_gpd_diagnostics <- function(fit, file_name, plot_title) {
  
  sigma_hat <- exp(as.vector(fit$X %*% fit$beta))
  xi_hat <- fit$xi
  
  fitted_cdf <- pgpd_custom(fit$y, sigma_hat, xi_hat)
  exp_residuals <- -log(1 - fitted_cdf)
  
  n <- length(fitted_cdf)
  pp <- ppoints(n)
  
  png(file_name, width = 1000, height = 800)
  
  par(mfrow = c(2, 2))
  
  plot(
    pp,
    sort(fitted_cdf),
    pch = 19,
    main = paste(plot_title, "Probability Plot"),
    xlab = "Theoretical probability",
    ylab = "Empirical fitted probability"
  )
  abline(0, 1, col = "red", lwd = 2)
  
  plot(
    qexp(pp),
    sort(exp_residuals),
    pch = 19,
    main = paste(plot_title, "Exponential QQ Plot"),
    xlab = "Theoretical exponential quantiles",
    ylab = "Model residual quantiles"
  )
  abline(0, 1, col = "red", lwd = 2)
  
  hist(
    exp_residuals,
    breaks = 30,
    probability = TRUE,
    main = paste(plot_title, "Residual Density"),
    xlab = "Exponential residuals"
  )
  curve(dexp(x), add = TRUE, col = "red", lwd = 2)
  
  plot(
    exp_residuals,
    type = "h",
    main = paste(plot_title, "Residuals by Observation"),
    xlab = "Exceedance index",
    ylab = "Exponential residual"
  )
  
  par(mfrow = c(1, 1))
  
  dev.off()
  
  return(exp_residuals)
}


resid_stationary <- save_custom_gpd_diagnostics(
  fit_stationary,
  file.path(fig_dir, "diagnostics_stationary_gpd.png"),
  "Stationary GPD"
)

resid_vol <- save_custom_gpd_diagnostics(
  fit_vol,
  file.path(fig_dir, "diagnostics_volatility_gpd.png"),
  "Volatility GPD"
)

resid_vol_time <- save_custom_gpd_diagnostics(
  fit_vol_time,
  file.path(fig_dir, "diagnostics_volatility_time_gpd.png"),
  "Volatility + Time GPD"
)


# ============================================================
# 16. Conditional Return Levels with Crisis Volatility Scenarios
#     Main Model: All Exceedances
# ============================================================

gpd_return_level <- function(m, u, sigma, xi, zeta) {
  
  if (m * zeta <= 1) {
    return(NA)
  }
  
  if (abs(xi) < 1e-6) {
    x_m <- u + sigma * log(m * zeta)
  } else {
    x_m <- u + (sigma / xi) * ((m * zeta)^xi - 1)
  }
  
  return(x_m)
}

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

vol_scenarios$zeta_hat <- predict(
  logit_vol,
  newdata = vol_scenarios,
  type = "response"
)

vol_scenarios$sigma_hat <- exp(
  fit_vol$beta[1] + fit_vol$beta[2] * vol_scenarios$z_vol
)

vol_scenarios$xi_hat <- fit_vol$xi

return_periods_days <- c(252, 1260, 2520)
return_labels <- c("1 year", "5 years", "10 years")

conditional_return_levels <- do.call(
  rbind,
  lapply(seq_len(nrow(vol_scenarios)), function(i) {
    
    levels <- sapply(
      return_periods_days,
      gpd_return_level,
      u = u,
      sigma = vol_scenarios$sigma_hat[i],
      xi = vol_scenarios$xi_hat[i],
      zeta = vol_scenarios$zeta_hat[i]
    )
    
    data.frame(
      model = "all exceedances volatility-dependent GPD",
      volatility_state = vol_scenarios$volatility_state[i],
      volatility_quantile = vol_scenarios$volatility_quantile[i],
      z_vol = vol_scenarios$z_vol[i],
      zeta_hat = vol_scenarios$zeta_hat[i],
      sigma_hat = vol_scenarios$sigma_hat[i],
      xi_hat = vol_scenarios$xi_hat[i],
      return_period = return_labels,
      m_trading_days = return_periods_days,
      return_level_decimal = as.numeric(levels),
      return_level_percent = 100 * as.numeric(levels)
    )
  })
)

print(vol_scenarios)
print(conditional_return_levels)

write.csv(
  vol_scenarios,
  file.path(res_dir, "volatility_scenarios.csv"),
  row.names = FALSE
)

write.csv(
  conditional_return_levels,
  file.path(res_dir, "conditional_return_levels_by_volatility.csv"),
  row.names = FALSE
)

stationary_sigma <- exp(fit_stationary$beta[1])
stationary_xi <- fit_stationary$xi
stationary_zeta <- mean(df$exceed)

stationary_return_levels <- data.frame(
  model = "stationary GPD, all exceedances",
  return_period = return_labels,
  m_trading_days = return_periods_days,
  return_level_decimal = sapply(
    return_periods_days,
    gpd_return_level,
    u = u,
    sigma = stationary_sigma,
    xi = stationary_xi,
    zeta = stationary_zeta
  )
)

stationary_return_levels$return_level_percent <- 100 * stationary_return_levels$return_level_decimal

write.csv(
  stationary_return_levels,
  file.path(res_dir, "stationary_return_levels.csv"),
  row.names = FALSE
)

plot_years <- c(1, 5, 10)

plot_colors <- c("darkgreen", "blue", "orange", "red", "darkred")
plot_pch <- c(19, 17, 15, 18, 8)

png(file.path(fig_dir, "conditional_return_levels_by_volatility.png"), width = 1000, height = 650)

plot(
  plot_years,
  conditional_return_levels$return_level_percent[
    conditional_return_levels$volatility_state == "low volatility"
  ],
  type = "b",
  pch = plot_pch[1],
  col = plot_colors[1],
  ylim = range(conditional_return_levels$return_level_percent, na.rm = TRUE),
  xlab = "Return Period in Years",
  ylab = "Return Level (%)",
  main = "Conditional Return Levels by Volatility State"
)

for (i in 2:nrow(vol_scenarios)) {
  
  state <- vol_scenarios$volatility_state[i]
  
  lines(
    plot_years,
    conditional_return_levels$return_level_percent[
      conditional_return_levels$volatility_state == state
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


png(file.path(fig_dir, "stationary_vs_conditional_return_levels.png"), width = 1000, height = 650)

plot(
  plot_years,
  stationary_return_levels$return_level_percent,
  type = "b",
  pch = 19,
  col = "black",
  ylim = range(
    stationary_return_levels$return_level_percent,
    conditional_return_levels$return_level_percent,
    na.rm = TRUE
  ),
  xlab = "Return Period in Years",
  ylab = "Return Level (%)",
  main = "Stationary vs Conditional Non-Stationary Return Levels"
)

for (i in seq_len(nrow(vol_scenarios))) {
  
  state <- vol_scenarios$volatility_state[i]
  
  lines(
    plot_years,
    conditional_return_levels$return_level_percent[
      conditional_return_levels$volatility_state == state
    ],
    type = "b",
    pch = plot_pch[i],
    col = plot_colors[i]
  )
}

legend(
  "topleft",
  legend = c("stationary GPD", vol_scenarios$volatility_state),
  col = c("black", plot_colors),
  pch = c(19, plot_pch),
  lty = 1,
  bty = "n"
)

dev.off()


# ============================================================
# 17. Chapter 5 Runs Declustering Robustness Extension
# ============================================================

# This section imports the Chapter 5 logic into Chapter 6.
#
# Instead of fitting the non-stationary GPD to all exceedances,
# we first group nearby exceedances into clusters.
#
# Then we keep only the maximum loss in each cluster and fit:
#
#   cluster maxima excesses -> non-stationary GPD with volatility
#
# This checks whether the volatility effect remains after accounting
# for clustering of extremes.

decluster_runs_with_covariates <- function(data, threshold, run_length = 5) {
  
  exceedance_indices <- which(data$loss > threshold)
  
  if (length(exceedance_indices) == 0) {
    return(list(exceedance_data = data.frame(), cluster_summary = data.frame()))
  }
  
  cluster_id <- 1
  cluster_ids <- rep(NA, length(exceedance_indices))
  cluster_ids[1] <- cluster_id
  
  for (i in 2:length(exceedance_indices)) {
    
    gap_non_exceedances <- exceedance_indices[i] - exceedance_indices[i - 1] - 1
    
    if (gap_non_exceedances >= run_length) {
      cluster_id <- cluster_id + 1
    }
    
    cluster_ids[i] <- cluster_id
  }
  
  exceedance_data <- data[exceedance_indices, ]
  exceedance_data$index <- exceedance_indices
  exceedance_data$cluster_id <- cluster_ids
  exceedance_data$excess <- exceedance_data$loss - threshold
  
  cluster_list <- split(exceedance_data, exceedance_data$cluster_id)
  
  cluster_summary <- do.call(
    rbind,
    lapply(cluster_list, function(cluster_df) {
      
      max_row <- cluster_df[which.max(cluster_df$loss), ]
      
      data.frame(
        cluster_id = unique(cluster_df$cluster_id),
        start_index = min(cluster_df$index),
        end_index = max(cluster_df$index),
        start_date = min(cluster_df$date),
        end_date = max(cluster_df$date),
        cluster_size = nrow(cluster_df),
        cluster_duration_trading_days = max(cluster_df$index) - min(cluster_df$index) + 1,
        cluster_max_loss = max_row$loss,
        cluster_max_excess = max_row$loss - threshold,
        cluster_max_date = max_row$date,
        z_vol_at_cluster_max = max_row$z_vol,
        t_scaled_at_cluster_max = max_row$t_scaled,
        rv21_lag_at_cluster_max = max_row$rv21_lag
      )
    })
  )
  
  row.names(cluster_summary) <- NULL
  
  return(list(
    exceedance_data = exceedance_data,
    cluster_summary = cluster_summary
  ))
}

main_run_length <- 5

declustered <- decluster_runs_with_covariates(
  data = df,
  threshold = u,
  run_length = main_run_length
)

cluster_summary <- declustered$cluster_summary
cluster_exceedance_data <- declustered$exceedance_data

n_clusters <- nrow(cluster_summary)
theta_hat <- n_clusters / n_exceedances

cat("\nDeclustering robustness summary:\n")
cat("Run length:", main_run_length, "\n")
cat("Raw exceedances:", n_exceedances, "\n")
cat("Clusters:", n_clusters, "\n")
cat("Extremal index:", theta_hat, "\n")

declustering_ch6_summary <- data.frame(
  run_length = main_run_length,
  threshold = u,
  threshold_percent = 100 * u,
  raw_exceedances = n_exceedances,
  clusters = n_clusters,
  extremal_index = theta_hat,
  mean_cluster_size = mean(cluster_summary$cluster_size),
  max_cluster_size = max(cluster_summary$cluster_size),
  cluster_rate_per_day = n_clusters / n_total,
  cluster_rate_per_year = 252 * n_clusters / n_total
)

print(declustering_ch6_summary)

write.csv(
  declustering_ch6_summary,
  file.path(res_dir, "chapter6_declustering_summary.csv"),
  row.names = FALSE
)

write.csv(
  cluster_summary,
  file.path(res_dir, "chapter6_clusters_run_length_5.csv"),
  row.names = FALSE
)

# Plot declustered exceedances and cluster maxima.

png(file.path(fig_dir, "chapter6_declustered_cluster_maxima.png"), width = 1100, height = 650)

plot(
  df$date,
  df$loss,
  type = "h",
  lwd = 0.5,
  col = "gray70",
  ylim = c(0, max(df$loss) * 1.05),
  main = "Chapter 6 Declustered Cluster Maxima, Run Length = 5",
  xlab = "Date",
  ylab = "Daily Loss"
)

abline(h = u, col = "red", lwd = 3, lty = 2)

points(
  cluster_exceedance_data$date,
  cluster_exceedance_data$loss,
  pch = 20,
  col = "gray30"
)

points(
  cluster_summary$cluster_max_date,
  cluster_summary$cluster_max_loss,
  pch = 19,
  col = "blue"
)

legend(
  "topright",
  legend = c("97.5% threshold", "Raw exceedances", "Cluster maxima"),
  col = c("red", "gray30", "blue"),
  lty = c(2, NA, NA),
  lwd = c(3, NA, NA),
  pch = c(NA, 20, 19),
  bty = "n"
)

dev.off()


# ============================================================
# 18. Non-Stationary GPD on Declustered Cluster Maxima
# ============================================================

# Cluster maxima excesses:
#
# Y_j = C_j - u
#
# where C_j is the maximum loss in cluster j.

y_cluster <- cluster_summary$cluster_max_excess
z_cluster <- cluster_summary$z_vol_at_cluster_max
t_cluster <- cluster_summary$t_scaled_at_cluster_max

# Initial stationary GPD fit to cluster maxima.

capture.output(
  initial_gpd_cluster <- gpd.fit(cluster_summary$cluster_max_loss, threshold = u)
)

initial_cluster_sigma <- initial_gpd_cluster$mle[1]
initial_cluster_xi <- initial_gpd_cluster$mle[2]

X_cluster_stationary <- matrix(1, nrow = length(y_cluster), ncol = 1)
colnames(X_cluster_stationary) <- "intercept"

X_cluster_vol <- cbind(1, z_cluster)
colnames(X_cluster_vol) <- c("intercept", "z_vol")

start_cluster_stationary <- c(log(initial_cluster_sigma), initial_cluster_xi)
start_cluster_vol <- c(log(initial_cluster_sigma), 0, initial_cluster_xi)

fit_cluster_stationary <- fit_gpd_model(
  y = y_cluster,
  X = X_cluster_stationary,
  start_par = start_cluster_stationary,
  model_name = "stationary GPD, declustered cluster maxima"
)

fit_cluster_vol <- fit_gpd_model(
  y = y_cluster,
  X = X_cluster_vol,
  start_par = start_cluster_vol,
  model_name = "volatility-dependent GPD, declustered cluster maxima"
)

cat("\nConvergence codes, declustered cluster maxima:\n")
cat("Cluster stationary:", fit_cluster_stationary$convergence, "\n")
cat("Cluster volatility:", fit_cluster_vol$convergence, "\n")

n_cluster_excess <- length(y_cluster)

cluster_model_comparison <- data.frame(
  model = c(
    fit_cluster_stationary$model_name,
    fit_cluster_vol$model_name
  ),
  n_parameters = c(
    length(fit_cluster_stationary$par),
    length(fit_cluster_vol$par)
  ),
  n_observations = c(n_cluster_excess, n_cluster_excess),
  nll = c(
    fit_cluster_stationary$nll,
    fit_cluster_vol$nll
  ),
  AIC = c(
    2 * length(fit_cluster_stationary$par) + 2 * fit_cluster_stationary$nll,
    2 * length(fit_cluster_vol$par) + 2 * fit_cluster_vol$nll
  ),
  BIC = c(
    log(n_cluster_excess) * length(fit_cluster_stationary$par) + 2 * fit_cluster_stationary$nll,
    log(n_cluster_excess) * length(fit_cluster_vol$par) + 2 * fit_cluster_vol$nll
  ),
  xi = c(
    fit_cluster_stationary$xi,
    fit_cluster_vol$xi
  ),
  se_xi = c(
    fit_cluster_stationary$se_xi,
    fit_cluster_vol$se_xi
  ),
  convergence = c(
    fit_cluster_stationary$convergence,
    fit_cluster_vol$convergence
  )
)

print(cluster_model_comparison)

write.csv(
  cluster_model_comparison,
  file.path(res_dir, "cluster_maxima_nonstationary_model_comparison.csv"),
  row.names = FALSE
)

cluster_lr_stat <- 2 * (fit_cluster_stationary$nll - fit_cluster_vol$nll)
cluster_lr_p_value <- pchisq(cluster_lr_stat, df = 1, lower.tail = FALSE)

cluster_lr_test <- data.frame(
  comparison = "stationary vs volatility-dependent scale, declustered cluster maxima",
  LR_statistic = cluster_lr_stat,
  df = 1,
  p_value = cluster_lr_p_value
)

print(cluster_lr_test)

write.csv(
  cluster_lr_test,
  file.path(res_dir, "cluster_maxima_likelihood_ratio_test.csv"),
  row.names = FALSE
)

cluster_parameter_table <- rbind(
  make_parameter_table(fit_cluster_stationary, c("beta_0", "xi")),
  make_parameter_table(fit_cluster_vol, c("beta_0", "beta_vol", "xi"))
)

print(cluster_parameter_table)

write.csv(
  cluster_parameter_table,
  file.path(res_dir, "cluster_maxima_parameter_estimates.csv"),
  row.names = FALSE
)


# Diagnostics for declustered cluster maxima models.

resid_cluster_stationary <- save_custom_gpd_diagnostics(
  fit_cluster_stationary,
  file.path(fig_dir, "diagnostics_cluster_stationary_gpd.png"),
  "Cluster Stationary GPD"
)

resid_cluster_vol <- save_custom_gpd_diagnostics(
  fit_cluster_vol,
  file.path(fig_dir, "diagnostics_cluster_volatility_gpd.png"),
  "Cluster Volatility GPD"
)


# ============================================================
# 19. Conditional Return Levels:
#     Declustered Cluster Maxima Robustness
# ============================================================

# For the declustered model, the event rate is the cluster rate,
# not the raw exceedance rate.
#
# We approximate the cluster probability by:
#
#   cluster_rate(z) = theta_hat * P(raw exceedance | z)
#
# where theta_hat is the runs estimate of the extremal index.
#
# This is a simple and transparent robustness approximation.
# It keeps the volatility-dependent exceedance probability but adjusts
# it downward to reflect clustering.

cluster_vol_scenarios <- vol_scenarios

cluster_vol_scenarios$raw_zeta_hat <- cluster_vol_scenarios$zeta_hat
cluster_vol_scenarios$cluster_zeta_hat <- theta_hat * cluster_vol_scenarios$raw_zeta_hat

cluster_vol_scenarios$cluster_sigma_hat <- exp(
  fit_cluster_vol$beta[1] + fit_cluster_vol$beta[2] * cluster_vol_scenarios$z_vol
)

cluster_vol_scenarios$cluster_xi_hat <- fit_cluster_vol$xi

cluster_conditional_return_levels <- do.call(
  rbind,
  lapply(seq_len(nrow(cluster_vol_scenarios)), function(i) {
    
    levels <- sapply(
      return_periods_days,
      gpd_return_level,
      u = u,
      sigma = cluster_vol_scenarios$cluster_sigma_hat[i],
      xi = cluster_vol_scenarios$cluster_xi_hat[i],
      zeta = cluster_vol_scenarios$cluster_zeta_hat[i]
    )
    
    data.frame(
      model = "declustered cluster maxima volatility-dependent GPD",
      volatility_state = cluster_vol_scenarios$volatility_state[i],
      volatility_quantile = cluster_vol_scenarios$volatility_quantile[i],
      z_vol = cluster_vol_scenarios$z_vol[i],
      raw_zeta_hat = cluster_vol_scenarios$raw_zeta_hat[i],
      cluster_zeta_hat = cluster_vol_scenarios$cluster_zeta_hat[i],
      sigma_hat = cluster_vol_scenarios$cluster_sigma_hat[i],
      xi_hat = cluster_vol_scenarios$cluster_xi_hat[i],
      return_period = return_labels,
      m_trading_days = return_periods_days,
      return_level_decimal = as.numeric(levels),
      return_level_percent = 100 * as.numeric(levels)
    )
  })
)

print(cluster_vol_scenarios)
print(cluster_conditional_return_levels)

write.csv(
  cluster_vol_scenarios,
  file.path(res_dir, "cluster_volatility_scenarios.csv"),
  row.names = FALSE
)

write.csv(
  cluster_conditional_return_levels,
  file.path(res_dir, "cluster_conditional_return_levels_by_volatility.csv"),
  row.names = FALSE
)

# Compare all-exceedance vs declustered conditional return levels.

combined_conditional_return_levels <- rbind(
  conditional_return_levels[, c(
    "model",
    "volatility_state",
    "volatility_quantile",
    "z_vol",
    "return_period",
    "m_trading_days",
    "return_level_percent"
  )],
  cluster_conditional_return_levels[, c(
    "model",
    "volatility_state",
    "volatility_quantile",
    "z_vol",
    "return_period",
    "m_trading_days",
    "return_level_percent"
  )]
)

write.csv(
  combined_conditional_return_levels,
  file.path(res_dir, "combined_conditional_return_levels_all_vs_cluster.csv"),
  row.names = FALSE
)

# Plot 10-year conditional return levels across volatility states.

ten_year_all <- conditional_return_levels[
  conditional_return_levels$return_period == "10 years",
]

ten_year_cluster <- cluster_conditional_return_levels[
  cluster_conditional_return_levels$return_period == "10 years",
]

png(file.path(fig_dir, "ten_year_return_levels_all_vs_cluster.png"), width = 1000, height = 650)

plot(
  vol_scenarios$volatility_quantile,
  ten_year_all$return_level_percent,
  type = "b",
  pch = 19,
  col = "black",
  ylim = range(
    ten_year_all$return_level_percent,
    ten_year_cluster$return_level_percent,
    na.rm = TRUE
  ),
  xlab = "Volatility Quantile",
  ylab = "10-Year Return Level (%)",
  main = "10-Year Conditional Return Levels: All vs Declustered"
)

lines(
  vol_scenarios$volatility_quantile,
  ten_year_cluster$return_level_percent,
  type = "b",
  pch = 17,
  col = "blue"
)

legend(
  "topleft",
  legend = c("All exceedances", "Declustered cluster maxima"),
  col = c("black", "blue"),
  pch = c(19, 17),
  lty = 1,
  bty = "n"
)

dev.off()


# ============================================================
# 20. Main Model vs Declustered Robustness Summary
# ============================================================

main_vs_cluster_model_summary <- data.frame(
  model = c(
    "all exceedances stationary GPD",
    "all exceedances volatility-dependent GPD",
    "cluster maxima stationary GPD",
    "cluster maxima volatility-dependent GPD"
  ),
  observations_used = c(
    length(y),
    length(y),
    length(y_cluster),
    length(y_cluster)
  ),
  xi = c(
    fit_stationary$xi,
    fit_vol$xi,
    fit_cluster_stationary$xi,
    fit_cluster_vol$xi
  ),
  se_xi = c(
    fit_stationary$se_xi,
    fit_vol$se_xi,
    fit_cluster_stationary$se_xi,
    fit_cluster_vol$se_xi
  ),
  beta_vol = c(
    NA,
    fit_vol$beta[2],
    NA,
    fit_cluster_vol$beta[2]
  ),
  se_beta_vol = c(
    NA,
    fit_vol$se_beta[2],
    NA,
    fit_cluster_vol$se_beta[2]
  ),
  AIC = c(
    model_comparison$AIC[1],
    model_comparison$AIC[2],
    cluster_model_comparison$AIC[1],
    cluster_model_comparison$AIC[2]
  ),
  BIC = c(
    model_comparison$BIC[1],
    model_comparison$BIC[2],
    cluster_model_comparison$BIC[1],
    cluster_model_comparison$BIC[2]
  )
)

print(main_vs_cluster_model_summary)

write.csv(
  main_vs_cluster_model_summary,
  file.path(res_dir, "main_vs_declustered_nonstationary_summary.csv"),
  row.names = FALSE
)

png(file.path(fig_dir, "main_vs_declustered_shape_comparison.png"), width = 900, height = 600)

barplot(
  main_vs_cluster_model_summary$xi,
  names.arg = c("All stat.", "All vol.", "Clust stat.", "Clust vol."),
  ylab = expression(hat(xi)),
  main = "Shape Parameter: All Exceedances vs Declustered",
  ylim = c(
    0,
    max(
      main_vs_cluster_model_summary$xi + main_vs_cluster_model_summary$se_xi,
      na.rm = TRUE
    ) * 1.4
  )
)

abline(h = 0, lty = 2, col = "red")

dev.off()


png(file.path(fig_dir, "main_vs_declustered_beta_vol_comparison.png"), width = 900, height = 600)

barplot(
  c(fit_vol$beta[2], fit_cluster_vol$beta[2]),
  names.arg = c("All exceedances", "Cluster maxima"),
  ylab = expression(hat(beta)[vol]),
  main = "Volatility Coefficient: All vs Declustered"
)

abline(h = 0, lty = 2, col = "red")

dev.off()


# ============================================================
# 21. Chapter 3, 4, 5, 6 Shape Comparison
# ============================================================

monthly_max_losses <- apply.monthly(losses, max)
monthly_max_losses_num <- as.numeric(monthly_max_losses)

capture.output(
  gev_ch3 <- gev.fit(monthly_max_losses_num)
)

capture.output(
  gpd_ch5_cluster <- gpd.fit(cluster_summary$cluster_max_loss, threshold = u)
)

chapter_shape_comparison <- data.frame(
  model = c(
    "Chapter 3 GEV monthly maxima",
    "Chapter 4 stationary GPD",
    "Chapter 5 declustered GPD",
    "Chapter 6 volatility-dependent GPD, all exceedances",
    "Chapter 6 volatility-dependent GPD, declustered cluster maxima"
  ),
  xi = c(
    gev_ch3$mle[3],
    fit_stationary$xi,
    gpd_ch5_cluster$mle[2],
    fit_vol$xi,
    fit_cluster_vol$xi
  ),
  se_xi = c(
    gev_ch3$se[3],
    fit_stationary$se_xi,
    gpd_ch5_cluster$se[2],
    fit_vol$se_xi,
    fit_cluster_vol$se_xi
  )
)

print(chapter_shape_comparison)

write.csv(
  chapter_shape_comparison,
  file.path(res_dir, "chapter3_to_chapter6_shape_comparison.csv"),
  row.names = FALSE
)

png(file.path(fig_dir, "chapter3_to_chapter6_shape_comparison.png"), width = 1000, height = 650)

barplot(
  chapter_shape_comparison$xi,
  names.arg = c("Ch3", "Ch4", "Ch5", "Ch6 all", "Ch6 clust."),
  ylab = expression(hat(xi)),
  main = "Shape Parameter Comparison Across Chapters",
  ylim = c(0, max(chapter_shape_comparison$xi + chapter_shape_comparison$se_xi, na.rm = TRUE) * 1.4)
)

abline(h = 0, lty = 2, col = "red")

dev.off()


# ============================================================
# 22. Save Summary Files
# ============================================================

sink(file.path(res_dir, "chapter6_model_summary.txt"))

cat("Chapter 6: Extremes of Non-Stationary Sequences\n")
cat("===============================================\n\n")

cat("Data:\n")
cat("Asset: SPY\n")
cat("Returns: daily log returns\n")
cat("Losses: negative daily log returns\n")
cat("Covariate: lagged 21-day realized volatility\n")
cat("Number of usable daily observations:", n_total, "\n\n")

cat("Threshold summary:\n")
print(threshold_summary)
cat("\n")

cat("Logistic exceedance model comparison:\n")
print(logit_model_comparison)
cat("\n")

cat("Logistic volatility coefficients:\n")
print(logit_vol_coefficients)
cat("\n")

cat("GPD non-stationary model comparison, all exceedances:\n")
print(model_comparison)
cat("\n")

cat("Likelihood ratio tests, all exceedances:\n")
print(lr_tests)
cat("\n")

cat("GPD parameter estimates, all exceedances:\n")
print(parameter_table)
cat("\n")

cat("Volatility scenarios, all exceedances:\n")
print(vol_scenarios)
cat("\n")

cat("Conditional return levels, all exceedances:\n")
print(conditional_return_levels)
cat("\n")

cat("Stationary return levels, all exceedances:\n")
print(stationary_return_levels)
cat("\n")

cat("Declustering robustness summary:\n")
print(declustering_ch6_summary)
cat("\n")

cat("Cluster maxima non-stationary model comparison:\n")
print(cluster_model_comparison)
cat("\n")

cat("Cluster maxima likelihood ratio test:\n")
print(cluster_lr_test)
cat("\n")

cat("Cluster maxima parameter estimates:\n")
print(cluster_parameter_table)
cat("\n")

cat("Cluster conditional return levels:\n")
print(cluster_conditional_return_levels)
cat("\n")

cat("Main vs declustered non-stationary summary:\n")
print(main_vs_cluster_model_summary)
cat("\n")

cat("Chapter 3 to Chapter 6 shape comparison:\n")
print(chapter_shape_comparison)
cat("\n")

cat("Interpretation guide:\n")
cat("The main Chapter 6 model uses all threshold exceedances and allows GPD scale to depend on volatility.\n")
cat("The robustness extension first applies Chapter 5 runs declustering and then fits a non-stationary GPD to cluster maxima.\n")
cat("If the volatility coefficient remains positive after declustering, volatility affects the severity of independent extreme episodes.\n")
cat("If the volatility effect weakens after declustering, part of the all-exceedance effect may reflect clustered crisis periods.\n")

sink()

sink(file.path(res_dir, "session_info.txt"))
print(sessionInfo())
sink()


# ============================================================
# 23. Final Console Summary
# ============================================================

cat("\n================ CHAPTER 6 SUMMARY ================\n")
cat("Non-stationary EVT analysis completed.\n\n")

cat("Threshold:", round(100 * u, 3), "%\n")
cat("Number of raw exceedances:", n_exceedances, "\n")
cat("Number of declustered clusters:", n_clusters, "\n")
cat("Extremal index estimate:", round(theta_hat, 4), "\n\n")

cat("Logistic exceedance model comparison:\n")
print(logit_model_comparison)

cat("\nGPD model comparison, all exceedances:\n")
print(model_comparison)

cat("\nLikelihood ratio tests, all exceedances:\n")
print(lr_tests)

cat("\nVolatility scenarios, all exceedances:\n")
print(vol_scenarios)

cat("\nConditional return levels, all exceedances:\n")
print(conditional_return_levels)

cat("\nStationary return levels, all exceedances:\n")
print(stationary_return_levels)

cat("\nCluster maxima model comparison:\n")
print(cluster_model_comparison)

cat("\nCluster maxima likelihood ratio test:\n")
print(cluster_lr_test)

cat("\nCluster conditional return levels:\n")
print(cluster_conditional_return_levels)

cat("\nMain vs declustered non-stationary summary:\n")
print(main_vs_cluster_model_summary)

cat("\nChapter 3 to Chapter 6 shape comparison:\n")
print(chapter_shape_comparison)

cat("\nFigures saved in:\n")
cat(fig_dir, "\n")

cat("\nResults saved in:\n")
cat(res_dir, "\n")

cat("\nSaved figures:\n")
print(list.files(fig_dir))

cat("\nSaved result files:\n")
print(list.files(res_dir))

cat("===================================================\n")