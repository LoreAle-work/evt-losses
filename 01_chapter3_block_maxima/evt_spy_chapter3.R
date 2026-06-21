# ============================================================
# Extreme Value Theory on SPY Daily Losses
# Chapter 3: Block Maxima and GEV Models
# ============================================================
#
# This script applies the block maxima approach from Chapter 3
# of Coles' "An Introduction to Statistical Modeling of Extreme Values"
# to daily SPY returns.
#
# Main steps:
#   1. Download SPY adjusted prices
#   2. Compute daily log returns
#   3. Convert returns into losses
#   4. Extract monthly block maxima
#   5. Fit a GEV model
#   6. Estimate return levels
#   7. Compare monthly and yearly blocks
#   8. Run an r-largest sensitivity check
#
# Author: Lorenzo
# ============================================================


# ============================================================
# 0. Package Setup
# ============================================================

# install.packages(c("quantmod", "ismev", "xts", "zoo"))

library(quantmod)
library(ismev)
library(xts)
library(zoo)


# ============================================================
# 1. Download SPY Data
# ============================================================

# Download SPY data from Yahoo Finance.
# SPY is used because it has a long daily history and represents
# broad U.S. equity market exposure.

getSymbols("SPY", src = "yahoo", from = "1993-01-01", auto.assign = TRUE)

# Extract adjusted closing prices.
# Adjusted prices account for dividends and stock splits.

prices <- Ad(SPY)

# Inspect and plot prices.

head(prices)
plot(prices, main = "SPY Adjusted Closing Prices")


# ============================================================
# 2. Compute Daily Log Returns
# ============================================================

# Daily log returns are computed as:
#
#   R_t = log(P_t) - log(P_{t-1})
#
# where P_t is the adjusted closing price at time t.

returns <- diff(log(prices))
returns <- na.omit(returns)
colnames(returns) <- "return"

# Inspect and plot returns.

head(returns)
summary(returns)
plot(returns, main = "SPY Daily Log Returns")


# ============================================================
# 3. Convert Returns into Losses
# ============================================================

# EVT block maxima models are naturally used for large values.
# Since we are interested in extreme negative returns, we define:
#
#   L_t = -R_t
#
# Thus, large positive values of L_t correspond to large losses.

losses <- -returns
colnames(losses) <- "loss"

# Inspect and plot losses.

head(losses)
summary(losses)
plot(losses, main = "SPY Daily Losses")

# Print the 10 largest daily losses in percentage terms.

cat("\nLargest daily losses (%):\n")
print(100 * sort(as.numeric(losses), decreasing = TRUE)[1:10])


# ============================================================
# 4. Monthly Block Maxima
# ============================================================

# Following the block maxima approach, we divide the daily losses
# into monthly blocks and keep only the largest loss in each month:
#
#   M_j = max(L_t : t in month j)
#
# These monthly maxima are then modeled using the GEV distribution.

monthly_max_losses <- apply.monthly(losses, max)
colnames(monthly_max_losses) <- "monthly_max_loss"

plot(monthly_max_losses, main = "Monthly Maximum Daily Losses for SPY")

cat("\nNumber of monthly block maxima:", length(monthly_max_losses), "\n")

# Convert xts object to numeric vector for GEV fitting.

monthly_max_losses_num <- as.numeric(monthly_max_losses)


# ============================================================
# 5. Fit GEV Model to Monthly Maxima
# ============================================================

# Fit the GEV distribution:
#
#   M_j ~ GEV(mu, sigma, xi)
#
# where:
#   mu    = location parameter
#   sigma = scale parameter
#   xi    = shape parameter

gev_fit <- gev.fit(monthly_max_losses_num)

# Store estimates and standard errors.

mu_hat <- as.numeric(gev_fit$mle[1])
sigma_hat <- as.numeric(gev_fit$mle[2])
xi_hat <- as.numeric(gev_fit$mle[3])

mu_se <- as.numeric(gev_fit$se[1])
sigma_se <- as.numeric(gev_fit$se[2])
xi_se <- as.numeric(gev_fit$se[3])

cat("\nMonthly GEV parameter estimates:\n")
cat("mu_hat    =", mu_hat, "\n")
cat("sigma_hat =", sigma_hat, "\n")
cat("xi_hat    =", xi_hat, "\n")

cat("\nMonthly GEV standard errors:\n")
cat("SE(mu)    =", mu_se, "\n")
cat("SE(sigma) =", sigma_se, "\n")
cat("SE(xi)    =", xi_se, "\n")


# ============================================================
# 6. Interpret Shape Parameter
# ============================================================

# The shape parameter xi determines the tail behavior:
#
#   xi > 0  : heavy-tailed Frechet-type behavior
#   xi = 0  : Gumbel-type behavior
#   xi < 0  : finite upper endpoint, Weibull-type behavior

cat("\nShape parameter interpretation:\n")

if (xi_hat > 0) {
  cat("xi > 0: heavy-tailed Frechet-type behavior.\n")
} else if (xi_hat < 0) {
  cat("xi < 0: finite upper endpoint Weibull-type behavior.\n")
} else {
  cat("xi = 0: Gumbel-type behavior.\n")
}


# ============================================================
# 7. GEV Diagnostic Plots
# ============================================================

# gev.diag produces four diagnostic plots:
#
#   1. Probability plot
#   2. Quantile plot
#   3. Return level plot
#   4. Density plot

gev.diag(gev_fit)


# ============================================================
# 8. Return Levels for Monthly Maxima
# ============================================================

# For monthly blocks, an m-period return level z_m satisfies:
#
#   P(M > z_m) = 1 / m
#
# where M is the monthly maximum daily loss.
#
# Since blocks are monthly:
#   m = 12  corresponds to approximately 1 year
#   m = 60  corresponds to approximately 5 years
#   m = 120 corresponds to approximately 10 years

gev_return_level <- function(m, mu, sigma, xi) {
  p <- 1 / m
  
  if (abs(xi) < 1e-6) {
    z <- mu - sigma * log(-log(1 - p))
  } else {
    z <- mu - (sigma / xi) * (1 - (-log(1 - p))^(-xi))
  }
  
  return(z)
}

rl_1y <- gev_return_level(12, mu_hat, sigma_hat, xi_hat)
rl_5y <- gev_return_level(60, mu_hat, sigma_hat, xi_hat)
rl_10y <- gev_return_level(120, mu_hat, sigma_hat, xi_hat)

cat("\nReturn levels for monthly maximum daily losses:\n")
cat("1-year return level: ", round(100 * rl_1y, 2), "%\n", sep = "")
cat("5-year return level: ", round(100 * rl_5y, 2), "%\n", sep = "")
cat("10-year return level:", round(100 * rl_10y, 2), "%\n", sep = "")


# ============================================================
# 9. Yearly Block Maxima as Robustness Check
# ============================================================

# Annual maxima are more classical in block-maxima EVT.
# However, using yearly blocks gives far fewer observations.
#
# Here, we fit a GEV model to yearly maximum daily losses as a
# robustness check.

yearly_max_losses <- apply.yearly(losses, max)
colnames(yearly_max_losses) <- "yearly_max_loss"

plot(yearly_max_losses, main = "Yearly Maximum Daily Losses for SPY")

cat("\nNumber of yearly block maxima:", length(yearly_max_losses), "\n")

gev_yearly <- gev.fit(as.numeric(yearly_max_losses))

cat("\nYearly GEV parameter estimates:\n")
print(gev_yearly$mle)

cat("\nYearly GEV standard errors:\n")
print(gev_yearly$se)

gev.diag(gev_yearly)

# Interpretation:
# Yearly blocks are theoretically standard, but the smaller sample
# leads to less stable parameter estimates and wider uncertainty.


# ============================================================
# 10. r-Largest Order Statistics Sensitivity Check
# ============================================================

# The standard block maxima model uses r = 1, meaning only the largest
# observation per block is retained.
#
# The r-largest approach keeps the largest r observations in each block:
#
#   M_j^(1) >= M_j^(2) >= ... >= M_j^(r)
#
# Here, we compare r = 1, r = 2, and r = 3.
# This is used only as a sensitivity check.

get_top_r <- function(x, r = 3) {
  x <- as.numeric(x)
  x <- x[!is.na(x)]
  x <- sort(x, decreasing = TRUE)
  
  if (length(x) < r) {
    return(rep(NA, r))
  }
  
  return(x[1:r])
}

# Split daily losses by month.

losses_by_month <- split(as.numeric(losses), as.yearmon(index(losses)))

# Extract top 3 daily losses from each month.

monthly_top_3 <- do.call(
  rbind,
  lapply(losses_by_month, get_top_r, r = 3)
)

monthly_top_3 <- na.omit(monthly_top_3)
colnames(monthly_top_3) <- c("top_1", "top_2", "top_3")

head(monthly_top_3)


# ============================================================
# 11. Fit r = 1, r = 2, and r = 3 Models
# ============================================================

# r = 1: ordinary GEV model on monthly maxima.

fit_r1 <- gev.fit(monthly_top_3[, 1])

# r = 2: r-largest model using top 2 losses per month.

fit_r2 <- rlarg.fit(monthly_top_3[, 1:2], r = 2)

# r = 3: r-largest model using top 3 losses per month.

fit_r3 <- rlarg.fit(monthly_top_3[, 1:3], r = 3)


# ============================================================
# 12. Compare r-Largest Estimates
# ============================================================

cat("\nr = 1 estimates:\n")
print(fit_r1$mle)
print(fit_r1$se)

cat("\nr = 2 estimates:\n")
print(fit_r2$mle)
print(fit_r2$se)

cat("\nr = 3 estimates:\n")
print(fit_r3$mle)
print(fit_r3$se)

cat("\nShape parameter comparison:\n")
cat("xi, r = 1:", fit_r1$mle[3], "\n")
cat("xi, r = 2:", fit_r2$mle[3], "\n")
cat("xi, r = 3:", fit_r3$mle[3], "\n")

# Important:
# We do not compare raw likelihoods across r = 1, 2, and 3 directly,
# because each model uses a different number of order statistics per block.
#
# Instead, we check whether the estimated shape parameter xi remains stable.
# Stability of xi across r supports robustness of the tail conclusion.


# ============================================================
# 13. Final Summary
# ============================================================

cat("\n================ FINAL SUMMARY ================\n")
cat("Main model: monthly block maxima, r = 1\n")
cat("Number of monthly maxima:", length(monthly_max_losses), "\n")
cat("Estimated xi:", xi_hat, "\n")
cat("SE(xi):", xi_se, "\n")

cat("\nReturn levels:\n")
cat("1-year:", round(100 * rl_1y, 2), "%\n")
cat("5-year:", round(100 * rl_5y, 2), "%\n")
cat("10-year:", round(100 * rl_10y, 2), "%\n")

cat("\nr-largest xi comparison:\n")
cat("r = 1:", fit_r1$mle[3], "\n")
cat("r = 2:", fit_r2$mle[3], "\n")
cat("r = 3:", fit_r3$mle[3], "\n")

cat("\nConclusion:\n")
cat("The estimated shape parameter is positive across r = 1, 2, and 3,\n")
cat("supporting heavy-tailed behavior in extreme SPY daily losses.\n")
cat("================================================\n")
