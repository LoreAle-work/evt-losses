# ============================================================
# Chapter 4: Threshold Models for SPY Daily Losses
# Peaks Over Threshold + Generalized Pareto Distribution
# ============================================================
#
# This script applies the threshold exceedance approach from
# Chapter 4 of Coles to daily SPY losses.
#
# Main steps:
#   1. Download SPY adjusted prices
#   2. Compute daily log returns
#   3. Convert returns into losses
#   4. Choose high thresholds
#   5. Fit Generalized Pareto Distribution models
#   6. Check threshold stability
#   7. Estimate return levels
#   8. Compare Chapter 4 GPD results with Chapter 3 GEV results
#   9. Save figures and numerical results
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

# Run this script from the root of your GitHub repository.
# If needed, uncomment and modify the following line:
#
# setwd("~/Desktop/UNI/Projects/EVT/evt-losses")

cat("Current working directory:\n")
print(getwd())

# If the script is run from the repository root, outputs go inside:
# 02_chapter4_threshold_models/figures
# 02_chapter4_threshold_models/results
#
# If the script is run from inside 02_chapter4_threshold_models,
# outputs go inside ./figures and ./results.

if (basename(getwd()) == "02_chapter4_threshold_models") {
  chapter_dir <- "."
} else {
  chapter_dir <- "02_chapter4_threshold_models"
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

# getSymbols downloads financial data from Yahoo Finance.
# auto.assign = TRUE creates an object called SPY in the environment.

getSymbols("SPY", src = "yahoo", from = "1993-01-01", auto.assign = TRUE)

# Ad extracts the adjusted closing price.
# Adjusted prices account for dividends and stock splits.

prices <- Ad(SPY)

head(prices)
tail(prices)

# Save adjusted price plot.

png(file.path(fig_dir, "spy_adjusted_prices.png"), width = 1000, height = 600)
plot(prices, main = "SPY Adjusted Closing Prices")
dev.off()


# ============================================================
# 3. Compute Daily Log Returns
# ============================================================

# Daily log returns are:
#
# R_t = log(P_t) - log(P_{t-1})
#
# diff(log(prices)) computes this difference.

returns <- diff(log(prices))

# The first return is NA because there is no previous price.
# na.omit removes this missing value.

returns <- na.omit(returns)

# Rename the column.

colnames(returns) <- "return"

head(returns)
summary(returns)

# Save returns plot.

png(file.path(fig_dir, "spy_daily_log_returns.png"), width = 1000, height = 600)
plot(returns, main = "SPY Daily Log Returns")
dev.off()


# ============================================================
# 4. Convert Returns into Losses
# ============================================================

# We study extreme negative returns.
# Define losses as:
#
# L_t = -R_t
#
# Then large negative returns become large positive losses.

losses <- -returns
colnames(losses) <- "loss"

# Convert the xts object into a numeric vector.
# Many EVT functions want a plain numeric vector.

losses_num <- as.numeric(losses)
loss_dates <- index(losses)

head(losses)
summary(losses)

# Save loss plot.

png(file.path(fig_dir, "spy_daily_losses.png"), width = 1000, height = 600)
plot(losses, main = "SPY Daily Losses")
dev.off()

# Print largest daily losses in percent.

cat("\nLargest daily losses (%):\n")
print(100 * sort(losses_num, decreasing = TRUE)[1:10])


# ============================================================
# 5. Threshold Choice
# ============================================================

# In Chapter 4, we choose a high threshold u and model
# exceedances above u.
#
# We start with three common thresholds:
#   95%
#   97.5%
#   99%
#
# These are empirical quantiles of the daily loss distribution.

u_95 <- as.numeric(quantile(losses_num, 0.95, names = FALSE))
u_975 <- as.numeric(quantile(losses_num, 0.975, names = FALSE))
u_99 <- as.numeric(quantile(losses_num, 0.99, names = FALSE))

cat("\nThresholds:\n")
cat("95% threshold:   ", round(100 * u_95, 3), "%\n", sep = "")
cat("97.5% threshold: ", round(100 * u_975, 3), "%\n", sep = "")
cat("99% threshold:   ", round(100 * u_99, 3), "%\n", sep = "")

# Count exceedances above each threshold.

threshold_summary <- data.frame(
  threshold_probability = c(0.95, 0.975, 0.99),
  threshold_value = c(u_95, u_975, u_99),
  threshold_percent = 100 * c(u_95, u_975, u_99),
  number_exceedances = c(
    sum(losses_num > u_95),
    sum(losses_num > u_975),
    sum(losses_num > u_99)
  ),
  exceedance_fraction = c(
    mean(losses_num > u_95),
    mean(losses_num > u_975),
    mean(losses_num > u_99)
  )
)

print(threshold_summary)

write.csv(
  threshold_summary,
  file.path(res_dir, "threshold_summary.csv"),
  row.names = FALSE
)


# ============================================================
# 6. Plot Daily Losses with Thresholds
# ============================================================

# Full loss plot with thresholds.

png(file.path(fig_dir, "daily_losses_with_thresholds.png"), width = 1100, height = 650)

plot(
  loss_dates,
  losses_num,
  type = "l",
  lwd = 0.5,
  col = "gray35",
  main = "SPY Daily Losses with Thresholds",
  xlab = "Date",
  ylab = "Daily Loss"
)

abline(h = 0, col = "black", lwd = 1)
abline(h = u_95, col = "red", lwd = 3, lty = 1)
abline(h = u_975, col = "blue", lwd = 3, lty = 2)
abline(h = u_99, col = "darkgreen", lwd = 3, lty = 3)

legend(
  "topright",
  legend = c(
    paste0("95% threshold = ", round(100 * u_95, 2), "%"),
    paste0("97.5% threshold = ", round(100 * u_975, 2), "%"),
    paste0("99% threshold = ", round(100 * u_99, 2), "%")
  ),
  col = c("red", "blue", "darkgreen"),
  lwd = 3,
  lty = c(1, 2, 3),
  bty = "n"
)

dev.off()


# Positive-tail zoom plot.
# This makes the thresholds and exceedances more visible.

png(file.path(fig_dir, "daily_losses_thresholds_tail_zoom.png"), width = 1100, height = 650)

plot(
  loss_dates,
  losses_num,
  type = "h",
  lwd = 0.6,
  col = "gray40",
  ylim = c(0, max(losses_num) * 1.05),
  main = "SPY Positive Daily Losses with Thresholds",
  xlab = "Date",
  ylab = "Daily Loss"
)

abline(h = u_95, col = "red", lwd = 3, lty = 1)
abline(h = u_975, col = "blue", lwd = 3, lty = 2)
abline(h = u_99, col = "darkgreen", lwd = 3, lty = 3)

points(
  loss_dates[losses_num > u_95],
  losses_num[losses_num > u_95],
  pch = 20,
  col = "red"
)

legend(
  "topright",
  legend = c(
    paste0("95% threshold = ", round(100 * u_95, 2), "%"),
    paste0("97.5% threshold = ", round(100 * u_975, 2), "%"),
    paste0("99% threshold = ", round(100 * u_99, 2), "%"),
    "Exceedances above 95%"
  ),
  col = c("red", "blue", "darkgreen", "red"),
  lwd = c(3, 3, 3, NA),
  lty = c(1, 2, 3, NA),
  pch = c(NA, NA, NA, 20),
  bty = "n"
)

dev.off()


# ============================================================
# 7. Mean Residual Life Plot
# ============================================================

# The default mrl.plot can be visually unhelpful for financial losses
# because it includes too much of the non-extreme distribution.
# We save it anyway for reference.

png(file.path(fig_dir, "mean_residual_life_plot_default.png"), width = 900, height = 600)
mrl.plot(losses_num)
dev.off()


# Custom mean residual life plot for high positive thresholds.
# This focuses on thresholds from the 90th percentile upward.

mrl_probs <- seq(0.90, 0.995, by = 0.005)
mrl_thresholds <- as.numeric(quantile(losses_num, mrl_probs, names = FALSE))

mean_excess <- sapply(mrl_thresholds, function(u) {
  excesses <- losses_num[losses_num > u] - u
  mean(excesses)
})

num_exceedances_mrl <- sapply(mrl_thresholds, function(u) {
  sum(losses_num > u)
})

mrl_table <- data.frame(
  threshold_probability = mrl_probs,
  threshold_value = mrl_thresholds,
  threshold_percent = 100 * mrl_thresholds,
  mean_excess = mean_excess,
  number_exceedances = num_exceedances_mrl
)

write.csv(
  mrl_table,
  file.path(res_dir, "mean_residual_life_table.csv"),
  row.names = FALSE
)

png(file.path(fig_dir, "mean_residual_life_high_thresholds.png"), width = 900, height = 600)

plot(
  100 * mrl_thresholds,
  mean_excess,
  type = "b",
  pch = 19,
  xlab = "Threshold (%)",
  ylab = "Mean Excess",
  main = "Mean Residual Life Plot for High Thresholds"
)

abline(v = 100 * u_95, col = "red", lwd = 2, lty = 1)
abline(v = 100 * u_975, col = "blue", lwd = 2, lty = 2)
abline(v = 100 * u_99, col = "darkgreen", lwd = 2, lty = 3)

legend(
  "topright",
  legend = c("95%", "97.5%", "99%"),
  col = c("red", "blue", "darkgreen"),
  lwd = 2,
  lty = c(1, 2, 3),
  bty = "n"
)

dev.off()


# ============================================================
# 8. Helper Function for Safe GPD Fitting
# ============================================================

# gpd.fit estimates the GPD parameters by maximum likelihood.
#
# For ismev::gpd.fit:
#   mle[1] = sigma, scale parameter
#   mle[2] = xi, shape parameter
#
# This helper wraps gpd.fit in tryCatch so the script does not
# crash if a threshold produces a numerical issue.

fit_gpd_safely <- function(x, threshold, npy = 252) {
  
  fit <- tryCatch(
    {
      tmp <- NULL
      capture.output(
        tmp <- gpd.fit(x, threshold = threshold, npy = npy)
      )
      tmp
    },
    error = function(e) {
      return(NULL)
    }
  )
  
  return(fit)
}


# ============================================================
# 9. Parameter Stability Across Thresholds
# ============================================================

# Fit GPD models across a range of thresholds.
# The goal is to see whether xi remains stable.
#
# If xi changes wildly across nearby thresholds, the model is sensitive.

threshold_probs <- seq(0.90, 0.99, by = 0.005)
threshold_values <- as.numeric(quantile(losses_num, threshold_probs, names = FALSE))

stability_results <- data.frame(
  threshold_probability = threshold_probs,
  threshold_value = threshold_values,
  threshold_percent = 100 * threshold_values,
  exceedances = NA,
  sigma = NA,
  xi = NA,
  modified_scale = NA,
  se_sigma = NA,
  se_xi = NA
)

for (i in seq_along(threshold_values)) {
  
  u <- threshold_values[i]
  k <- sum(losses_num > u)
  
  # Skip thresholds with too few exceedances.
  # Too few exceedances means unstable estimation.
  
  if (k < 30) {
    next
  }
  
  fit <- fit_gpd_safely(losses_num, threshold = u, npy = 252)
  
  if (is.null(fit)) {
    next
  }
  
  sigma_hat <- fit$mle[1]
  xi_hat <- fit$mle[2]
  
  stability_results$exceedances[i] <- k
  stability_results$sigma[i] <- sigma_hat
  stability_results$xi[i] <- xi_hat
  stability_results$modified_scale[i] <- sigma_hat - xi_hat * u
  stability_results$se_sigma[i] <- fit$se[1]
  stability_results$se_xi[i] <- fit$se[2]
}

print(stability_results)

write.csv(
  stability_results,
  file.path(res_dir, "parameter_stability.csv"),
  row.names = FALSE
)


# ============================================================
# 10. Plot Parameter Stability
# ============================================================

valid <- !is.na(stability_results$xi)

# Shape parameter stability without confidence intervals.

png(file.path(fig_dir, "shape_parameter_stability.png"), width = 900, height = 600)

plot(
  stability_results$threshold_probability[valid],
  stability_results$xi[valid],
  type = "b",
  pch = 19,
  col = "black",
  xlab = "Threshold Quantile",
  ylab = expression(hat(xi)),
  main = "GPD Shape Parameter Stability"
)

abline(h = 0, lty = 2, col = "red")

dev.off()


# Shape parameter stability with approximate 95% confidence intervals.

xi_lower <- stability_results$xi[valid] - 1.96 * stability_results$se_xi[valid]
xi_upper <- stability_results$xi[valid] + 1.96 * stability_results$se_xi[valid]

png(file.path(fig_dir, "shape_parameter_stability_with_ci.png"), width = 900, height = 600)

plot(
  stability_results$threshold_probability[valid],
  stability_results$xi[valid],
  type = "b",
  pch = 19,
  ylim = range(c(xi_lower, xi_upper), na.rm = TRUE),
  xlab = "Threshold Quantile",
  ylab = expression(hat(xi)),
  main = "GPD Shape Parameter Stability with 95% CI"
)

arrows(
  stability_results$threshold_probability[valid],
  xi_lower,
  stability_results$threshold_probability[valid],
  xi_upper,
  angle = 90,
  code = 3,
  length = 0.04
)

abline(h = 0, lty = 2, col = "red")

dev.off()


# Modified scale parameter stability.

png(file.path(fig_dir, "modified_scale_stability.png"), width = 900, height = 600)

plot(
  stability_results$threshold_probability[valid],
  stability_results$modified_scale[valid],
  type = "b",
  pch = 19,
  col = "black",
  xlab = "Threshold Quantile",
  ylab = "Modified scale estimate",
  main = "Modified Scale Parameter Stability"
)

dev.off()


# Number of exceedances by threshold.

png(file.path(fig_dir, "number_exceedances_by_threshold.png"), width = 900, height = 600)

plot(
  stability_results$threshold_probability,
  stability_results$exceedances,
  type = "b",
  pch = 19,
  col = "black",
  xlab = "Threshold Quantile",
  ylab = "Number of Exceedances",
  main = "Number of Exceedances by Threshold"
)

dev.off()


# ============================================================
# 11. Fit Main GPD Models
# ============================================================

# Fit GPD models at 95%, 97.5%, and 99% thresholds.

gpd_fit_95 <- fit_gpd_safely(losses_num, threshold = u_95, npy = 252)
gpd_fit_975 <- fit_gpd_safely(losses_num, threshold = u_975, npy = 252)
gpd_fit_99 <- fit_gpd_safely(losses_num, threshold = u_99, npy = 252)

if (is.null(gpd_fit_95) || is.null(gpd_fit_975) || is.null(gpd_fit_99)) {
  stop("At least one main GPD fit failed.")
}

# Build a table of parameter estimates.

gpd_estimates <- data.frame(
  threshold = c("95%", "97.5%", "99%"),
  threshold_value = c(u_95, u_975, u_99),
  threshold_percent = 100 * c(u_95, u_975, u_99),
  exceedances = c(
    sum(losses_num > u_95),
    sum(losses_num > u_975),
    sum(losses_num > u_99)
  ),
  sigma = c(
    gpd_fit_95$mle[1],
    gpd_fit_975$mle[1],
    gpd_fit_99$mle[1]
  ),
  xi = c(
    gpd_fit_95$mle[2],
    gpd_fit_975$mle[2],
    gpd_fit_99$mle[2]
  ),
  se_sigma = c(
    gpd_fit_95$se[1],
    gpd_fit_975$se[1],
    gpd_fit_99$se[1]
  ),
  se_xi = c(
    gpd_fit_95$se[2],
    gpd_fit_975$se[2],
    gpd_fit_99$se[2]
  )
)

print(gpd_estimates)

write.csv(
  gpd_estimates,
  file.path(res_dir, "gpd_estimates_by_threshold.csv"),
  row.names = FALSE
)


# ============================================================
# 12. GPD Diagnostic Plots
# ============================================================

# gpd.diag creates:
#   1. Probability plot
#   2. Quantile plot
#   3. Return level plot
#   4. Density plot

png(file.path(fig_dir, "gpd_diagnostics_95.png"), width = 1000, height = 800)
gpd.diag(gpd_fit_95)
dev.off()

png(file.path(fig_dir, "gpd_diagnostics_975.png"), width = 1000, height = 800)
gpd.diag(gpd_fit_975)
dev.off()

png(file.path(fig_dir, "gpd_diagnostics_99.png"), width = 1000, height = 800)
gpd.diag(gpd_fit_99)
dev.off()


# ============================================================
# 13. Manual GPD Return Level Function
# ============================================================

# For daily data, return periods are measured in trading days.
#
# 252 trading days approx 1 year
# 1260 trading days approx 5 years
# 2520 trading days approx 10 years
#
# Let zeta_u = P(X > u), estimated by exceedances / total observations.
#
# For xi != 0:
#
# x_m = u + (sigma / xi) * ((m * zeta_u)^xi - 1)

gpd_return_level <- function(m, u, sigma, xi, zeta_u) {
  
  if (abs(xi) < 1e-6) {
    x_m <- u + sigma * log(m * zeta_u)
  } else {
    x_m <- u + (sigma / xi) * ((m * zeta_u)^xi - 1)
  }
  
  return(x_m)
}


# ============================================================
# 14. Compute GPD Return Levels
# ============================================================

return_periods_daily <- c(252, 1260, 2520)
return_labels <- c("1 year", "5 years", "10 years")

compute_return_levels_for_fit <- function(fit, threshold, losses_num, label) {
  
  sigma_hat <- fit$mle[1]
  xi_hat <- fit$mle[2]
  zeta_hat <- mean(losses_num > threshold)
  
  return_levels <- sapply(
    return_periods_daily,
    gpd_return_level,
    u = threshold,
    sigma = sigma_hat,
    xi = xi_hat,
    zeta_u = zeta_hat
  )
  
  output <- data.frame(
    threshold = label,
    return_period = return_labels,
    m_trading_days = return_periods_daily,
    return_level_decimal = as.numeric(return_levels),
    return_level_percent = 100 * as.numeric(return_levels)
  )
  
  return(output)
}

rl_95 <- compute_return_levels_for_fit(gpd_fit_95, u_95, losses_num, "95%")
rl_975 <- compute_return_levels_for_fit(gpd_fit_975, u_975, losses_num, "97.5%")
rl_99 <- compute_return_levels_for_fit(gpd_fit_99, u_99, losses_num, "99%")

return_levels_gpd <- rbind(rl_95, rl_975, rl_99)

print(return_levels_gpd)

write.csv(
  return_levels_gpd,
  file.path(res_dir, "gpd_return_levels.csv"),
  row.names = FALSE
)


# ============================================================
# 15. Plot Return Levels Across Thresholds
# ============================================================

png(file.path(fig_dir, "gpd_return_levels_comparison.png"), width = 900, height = 600)

plot(
  return_periods_daily,
  rl_95$return_level_percent,
  type = "b",
  pch = 19,
  col = "red",
  ylim = range(return_levels_gpd$return_level_percent),
  xlab = "Return Period in Trading Days",
  ylab = "Return Level (%)",
  main = "GPD Return Levels Across Thresholds"
)

lines(return_periods_daily, rl_975$return_level_percent, type = "b", pch = 17, col = "blue")
lines(return_periods_daily, rl_99$return_level_percent, type = "b", pch = 15, col = "darkgreen")

legend(
  "topleft",
  legend = c("95%", "97.5%", "99%"),
  col = c("red", "blue", "darkgreen"),
  pch = c(19, 17, 15),
  lty = 1,
  bty = "n"
)

dev.off()


# ============================================================
# 16. Histograms of Exceedances
# ============================================================

# These plots show the excesses:
#
# Y = X - u | X > u

excesses_95 <- losses_num[losses_num > u_95] - u_95
excesses_975 <- losses_num[losses_num > u_975] - u_975
excesses_99 <- losses_num[losses_num > u_99] - u_99

png(file.path(fig_dir, "excesses_histogram_95.png"), width = 900, height = 600)
hist(
  excesses_95,
  breaks = 40,
  main = "Threshold Excesses above 95% Threshold",
  xlab = "Excess over Threshold"
)
dev.off()

png(file.path(fig_dir, "excesses_histogram_975.png"), width = 900, height = 600)
hist(
  excesses_975,
  breaks = 35,
  main = "Threshold Excesses above 97.5% Threshold",
  xlab = "Excess over Threshold"
)
dev.off()

png(file.path(fig_dir, "excesses_histogram_99.png"), width = 900, height = 600)
hist(
  excesses_99,
  breaks = 25,
  main = "Threshold Excesses above 99% Threshold",
  xlab = "Excess over Threshold"
)
dev.off()


# ============================================================
# 17. Chapter 3 Comparison: Monthly Block Maxima GEV
# ============================================================

# To compare Chapter 4 with Chapter 3 using the same data,
# we refit the Chapter 3 monthly block maxima GEV model here.
#
# Chapter 3:
#   monthly maxima -> GEV
#
# Chapter 4:
#   daily exceedances over threshold -> GPD

monthly_max_losses <- apply.monthly(losses, max)
colnames(monthly_max_losses) <- "monthly_max_loss"
monthly_max_losses_num <- as.numeric(monthly_max_losses)

capture.output(
  gev_ch3_fit <- gev.fit(monthly_max_losses_num)
)

gev_ch3_mu <- gev_ch3_fit$mle[1]
gev_ch3_sigma <- gev_ch3_fit$mle[2]
gev_ch3_xi <- gev_ch3_fit$mle[3]

gev_return_level <- function(m, mu, sigma, xi) {
  
  p <- 1 / m
  
  if (abs(xi) < 1e-6) {
    z <- mu - sigma * log(-log(1 - p))
  } else {
    z <- mu - (sigma / xi) * (1 - (-log(1 - p))^(-xi))
  }
  
  return(z)
}

# Monthly return periods:
# 12 months = 1 year
# 60 months = 5 years
# 120 months = 10 years

gev_ch3_return_levels <- data.frame(
  method = "Chapter 3 GEV monthly block maxima",
  return_period = return_labels,
  m_months = c(12, 60, 120),
  return_level_decimal = c(
    gev_return_level(12, gev_ch3_mu, gev_ch3_sigma, gev_ch3_xi),
    gev_return_level(60, gev_ch3_mu, gev_ch3_sigma, gev_ch3_xi),
    gev_return_level(120, gev_ch3_mu, gev_ch3_sigma, gev_ch3_xi)
  )
)

gev_ch3_return_levels$return_level_percent <- 100 * gev_ch3_return_levels$return_level_decimal

# Main Chapter 4 model: 97.5% threshold.
# This is chosen as a balance between enough exceedances and tail focus.

gpd_ch4_main_return_levels <- rl_975
gpd_ch4_main_return_levels$method <- "Chapter 4 GPD threshold model, 97.5% threshold"

# Build comparison table.

chapter3_chapter4_comparison <- data.frame(
  return_period = return_labels,
  chapter3_gev_percent = gev_ch3_return_levels$return_level_percent,
  chapter4_gpd_975_percent = gpd_ch4_main_return_levels$return_level_percent
)

print(chapter3_chapter4_comparison)

write.csv(
  chapter3_chapter4_comparison,
  file.path(res_dir, "chapter3_chapter4_return_level_comparison.csv"),
  row.names = FALSE
)

shape_comparison <- data.frame(
  model = c(
    "Chapter 3 GEV monthly block maxima",
    "Chapter 4 GPD 95% threshold",
    "Chapter 4 GPD 97.5% threshold",
    "Chapter 4 GPD 99% threshold"
  ),
  xi = c(
    gev_ch3_xi,
    gpd_fit_95$mle[2],
    gpd_fit_975$mle[2],
    gpd_fit_99$mle[2]
  ),
  se_xi = c(
    gev_ch3_fit$se[3],
    gpd_fit_95$se[2],
    gpd_fit_975$se[2],
    gpd_fit_99$se[2]
  )
)

print(shape_comparison)

write.csv(
  shape_comparison,
  file.path(res_dir, "chapter3_chapter4_shape_comparison.csv"),
  row.names = FALSE
)


# Plot Chapter 3 vs Chapter 4 return levels.

png(file.path(fig_dir, "chapter3_vs_chapter4_return_levels.png"), width = 900, height = 600)

plot(
  c(1, 5, 10),
  chapter3_chapter4_comparison$chapter3_gev_percent,
  type = "b",
  pch = 19,
  col = "black",
  ylim = range(
    chapter3_chapter4_comparison$chapter3_gev_percent,
    chapter3_chapter4_comparison$chapter4_gpd_975_percent
  ),
  xlab = "Return Period in Years",
  ylab = "Return Level (%)",
  main = "Chapter 3 GEV vs Chapter 4 GPD Return Levels"
)

lines(
  c(1, 5, 10),
  chapter3_chapter4_comparison$chapter4_gpd_975_percent,
  type = "b",
  pch = 17,
  col = "blue"
)

legend(
  "topleft",
  legend = c("Chapter 3 GEV", "Chapter 4 GPD, 97.5%"),
  col = c("black", "blue"),
  pch = c(19, 17),
  lty = 1,
  bty = "n"
)

dev.off()


# Plot shape parameter comparison.

png(file.path(fig_dir, "chapter3_vs_chapter4_shape_comparison.png"), width = 900, height = 600)

barplot(
  shape_comparison$xi,
  names.arg = c("GEV", "GPD 95", "GPD 97.5", "GPD 99"),
  ylab = expression(hat(xi)),
  main = "Shape Parameter Comparison",
  ylim = c(0, max(shape_comparison$xi + shape_comparison$se_xi, na.rm = TRUE) * 1.3)
)

abline(h = 0, lty = 2, col = "red")

dev.off()


# ============================================================
# 18. Save Plain Text Model Summary
# ============================================================

sink(file.path(res_dir, "chapter4_model_summary.txt"))

cat("Chapter 4: Threshold Models for SPY Daily Losses\n")
cat("================================================\n\n")

cat("Data:\n")
cat("Asset: SPY\n")
cat("Returns: daily log returns\n")
cat("Losses: negative daily log returns\n")
cat("Number of daily observations:", length(losses_num), "\n\n")

cat("Largest daily losses (%):\n")
print(100 * sort(losses_num, decreasing = TRUE)[1:10])
cat("\n")

cat("Threshold summary:\n")
print(threshold_summary)
cat("\n")

cat("GPD estimates by threshold:\n")
print(gpd_estimates)
cat("\n")

cat("GPD return levels:\n")
print(return_levels_gpd)
cat("\n")

cat("Parameter stability results:\n")
print(stability_results)
cat("\n")

cat("Chapter 3 GEV comparison:\n")
cat("GEV xi:", gev_ch3_xi, "\n")
cat("GEV SE xi:", gev_ch3_fit$se[3], "\n\n")

cat("Shape comparison:\n")
print(shape_comparison)
cat("\n")

cat("Return level comparison:\n")
print(chapter3_chapter4_comparison)
cat("\n")

cat("Main interpretation:\n")
cat("The threshold model estimates the tail behavior of daily SPY losses\n")
cat("using exceedances above high thresholds. The key parameter is xi.\n")
cat("Positive xi estimates suggest heavy-tailed behavior in extreme losses.\n")
cat("The Chapter 4 GPD results are broadly consistent with the Chapter 3 GEV results.\n")

sink()


# ============================================================
# 19. Save Session Info
# ============================================================

sink(file.path(res_dir, "session_info.txt"))
print(sessionInfo())
sink()


# ============================================================
# 20. Final Console Summary
# ============================================================

cat("\n================ CHAPTER 4 SUMMARY ================\n")
cat("Threshold models completed.\n\n")

cat("Figures saved in:\n")
cat(fig_dir, "\n\n")

cat("Results saved in:\n")
cat(res_dir, "\n\n")

cat("Threshold summary:\n")
print(threshold_summary)

cat("\nGPD shape estimates:\n")
print(gpd_estimates[, c("threshold", "xi", "se_xi")])

cat("\nGPD return levels:\n")
print(return_levels_gpd)

cat("\nChapter 3 vs Chapter 4 shape comparison:\n")
print(shape_comparison)

cat("\nChapter 3 vs Chapter 4 return level comparison:\n")
print(chapter3_chapter4_comparison)

cat("\nSaved figures:\n")
print(list.files(fig_dir))

cat("\nSaved result files:\n")
print(list.files(res_dir))

cat("===================================================\n")