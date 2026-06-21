# ============================================================
# Chapter 5: Extremes of Dependent Sequences
# SPY Daily Losses, Clustering, Declustering, Extremal Index
# ============================================================
#
# This script extends the Chapter 4 threshold model by allowing
# for dependence among extremes.
#
# Main idea:
#   Chapter 4 treated threshold exceedances as approximately independent.
#   Chapter 5 checks whether exceedances cluster over time.
#
# Main steps:
#   1. Download SPY adjusted prices
#   2. Compute daily log returns
#   3. Convert returns into losses
#   4. Choose a high threshold
#   5. Identify exceedances
#   6. Study clustering of exceedances
#   7. Estimate the extremal index
#   8. Apply runs declustering
#   9. Fit GPD to cluster maxima
#   10. Compare Chapter 4 iid threshold model with Chapter 5 declustered model
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
# If needed, uncomment and modify:
#
setwd("~/Desktop/UNI/Projects/EVT/evt-losses")

cat("Current working directory:\n")
print(getwd())

# If the script is run from the repository root, outputs go inside:
# 03_chapter5_dependent_sequences/figures
# 03_chapter5_dependent_sequences/results
#
# If the script is run from inside 03_chapter5_dependent_sequences,
# outputs go inside ./figures and ./results.

if (basename(getwd()) == "03_chapter5_dependent_sequences") {
  chapter_dir <- "."
} else {
  chapter_dir <- "03_chapter5_dependent_sequences"
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

# getSymbols downloads SPY data from Yahoo Finance.
# auto.assign = TRUE creates an object called SPY.

getSymbols("SPY", src = "yahoo", from = "1993-01-01", auto.assign = TRUE)

# Extract adjusted closing prices.

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

head(returns)
summary(returns)

png(file.path(fig_dir, "spy_daily_log_returns.png"), width = 1000, height = 600)
plot(returns, main = "SPY Daily Log Returns")
dev.off()


# ============================================================
# 4. Convert Returns into Losses
# ============================================================

# We study extreme negative returns.
# Define losses:
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
# 5. Choose Threshold
# ============================================================

# Chapter 4 considered 95%, 97.5%, and 99% thresholds.
# For Chapter 5, we use 97.5% as the main threshold.
#
# Why 97.5%?
#   - More tail-focused than 95%
#   - More observations than 99%
#   - Good compromise between bias and variance

u <- as.numeric(quantile(losses_num, 0.975, names = FALSE))

cat("\nMain threshold:\n")
cat("97.5% threshold =", round(100 * u, 3), "%\n")

exceedance_indicator <- losses_num > u
exceedance_indices <- which(exceedance_indicator)
exceedance_dates <- loss_dates[exceedance_indices]
exceedance_losses <- losses_num[exceedance_indices]

n_total <- length(losses_num)
n_exceedances <- length(exceedance_indices)
exceedance_rate <- n_exceedances / n_total

cat("Number of daily observations:", n_total, "\n")
cat("Number of exceedances:", n_exceedances, "\n")
cat("Exceedance fraction:", exceedance_rate, "\n")

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


# ============================================================
# 6. Plot Exceedances Over Time
# ============================================================

png(file.path(fig_dir, "daily_losses_with_exceedances.png"), width = 1100, height = 650)

plot(
  loss_dates,
  losses_num,
  type = "h",
  lwd = 0.5,
  col = "gray50",
  ylim = c(0, max(losses_num) * 1.05),
  main = "SPY Daily Losses and 97.5% Threshold Exceedances",
  xlab = "Date",
  ylab = "Daily Loss"
)

abline(h = u, col = "red", lwd = 3, lty = 2)

points(
  exceedance_dates,
  exceedance_losses,
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


# ============================================================
# 7. Dependence Diagnostics
# ============================================================

# Chapter 5 asks whether extremes are dependent.
#
# We look at:
#   - ACF of returns
#   - ACF of losses
#   - ACF of absolute returns
#   - ACF of exceedance indicators
#
# Absolute returns and exceedance indicators are especially useful
# because financial volatility tends to cluster.

png(file.path(fig_dir, "acf_returns.png"), width = 900, height = 600)
acf(as.numeric(returns), main = "ACF of SPY Daily Log Returns", na.action = na.pass)
dev.off()

png(file.path(fig_dir, "acf_losses.png"), width = 900, height = 600)
acf(losses_num, main = "ACF of SPY Daily Losses", na.action = na.pass)
dev.off()

png(file.path(fig_dir, "acf_absolute_returns.png"), width = 900, height = 600)
acf(abs(as.numeric(returns)), main = "ACF of Absolute SPY Returns", na.action = na.pass)
dev.off()

png(file.path(fig_dir, "acf_exceedance_indicator.png"), width = 900, height = 600)
acf(as.numeric(exceedance_indicator), main = "ACF of Threshold Exceedance Indicator", na.action = na.pass)
dev.off()


# ============================================================
# 8. Runs Declustering Function
# ============================================================

# Runs declustering groups nearby exceedances into clusters.
#
# A new cluster starts only after at least 'run_length' consecutive
# non-exceedance days.
#
# Example:
#   run_length = 5
#
# means that exceedances are considered part of the same cluster
# unless they are separated by at least 5 non-exceedance trading days.
#
# This is useful because financial losses often cluster during crises.

decluster_runs <- function(x, dates, threshold, run_length = 5) {
  
  exceedance_indices <- which(x > threshold)
  
  if (length(exceedance_indices) == 0) {
    return(data.frame())
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
  
  exceedance_data <- data.frame(
    cluster_id = cluster_ids,
    index = exceedance_indices,
    date = dates[exceedance_indices],
    loss = x[exceedance_indices],
    excess = x[exceedance_indices] - threshold
  )
  
  cluster_list <- split(exceedance_data, exceedance_data$cluster_id)
  
  cluster_summary <- do.call(
    rbind,
    lapply(cluster_list, function(df) {
      
      max_row <- df[which.max(df$loss), ]
      
      data.frame(
        cluster_id = unique(df$cluster_id),
        start_index = min(df$index),
        end_index = max(df$index),
        start_date = min(df$date),
        end_date = max(df$date),
        cluster_size = nrow(df),
        cluster_duration_trading_days = max(df$index) - min(df$index) + 1,
        cluster_max_loss = max(df$loss),
        cluster_max_excess = max(df$loss) - threshold,
        cluster_max_date = max_row$date
      )
    })
  )
  
  row.names(cluster_summary) <- NULL
  
  return(list(
    exceedance_data = exceedance_data,
    cluster_summary = cluster_summary
  ))
}


# ============================================================
# 9. Apply Declustering for Several Run Lengths
# ============================================================

# We compare different run lengths because cluster definitions are partly subjective.
# A longer run length usually creates fewer, larger clusters.

run_lengths <- c(1, 3, 5, 10)

declustering_results <- data.frame(
  run_length = run_lengths,
  threshold = u,
  threshold_percent = 100 * u,
  number_exceedances = n_exceedances,
  number_clusters = NA,
  mean_cluster_size = NA,
  max_cluster_size = NA,
  extremal_index = NA,
  cluster_rate_per_day = NA,
  cluster_rate_per_year = NA
)

cluster_objects <- list()

for (j in seq_along(run_lengths)) {
  
  r <- run_lengths[j]
  
  dec <- decluster_runs(
    x = losses_num,
    dates = loss_dates,
    threshold = u,
    run_length = r
  )
  
  cluster_objects[[as.character(r)]] <- dec
  
  clusters <- dec$cluster_summary
  
  number_clusters <- nrow(clusters)
  mean_cluster_size <- mean(clusters$cluster_size)
  max_cluster_size <- max(clusters$cluster_size)
  
  # Simple runs estimate of extremal index:
  #
  # theta_hat = number of clusters / number of exceedances
  
  theta_hat <- number_clusters / n_exceedances
  
  declustering_results$number_clusters[j] <- number_clusters
  declustering_results$mean_cluster_size[j] <- mean_cluster_size
  declustering_results$max_cluster_size[j] <- max_cluster_size
  declustering_results$extremal_index[j] <- theta_hat
  declustering_results$cluster_rate_per_day[j] <- number_clusters / n_total
  declustering_results$cluster_rate_per_year[j] <- 252 * number_clusters / n_total
}

print(declustering_results)

write.csv(
  declustering_results,
  file.path(res_dir, "declustering_results.csv"),
  row.names = FALSE
)


# ============================================================
# 10. Save Cluster Tables
# ============================================================

# Save cluster summaries for each run length.

for (r in run_lengths) {
  
  clusters <- cluster_objects[[as.character(r)]]$cluster_summary
  
  write.csv(
    clusters,
    file.path(res_dir, paste0("clusters_run_length_", r, ".csv")),
    row.names = FALSE
  )
}


# ============================================================
# 11. Main Declustering Choice
# ============================================================

# We use run_length = 5 as the main specification.
#
# Interpretation:
#   A cluster ends only after at least 5 consecutive non-exceedance
#   trading days.
#
# This roughly corresponds to one trading week.

main_run_length <- 5
main_clusters <- cluster_objects[[as.character(main_run_length)]]$cluster_summary
main_exceedance_data <- cluster_objects[[as.character(main_run_length)]]$exceedance_data

main_theta <- nrow(main_clusters) / n_exceedances

cat("\nMain runs declustering choice:\n")
cat("Run length:", main_run_length, "\n")
cat("Number of exceedances:", n_exceedances, "\n")
cat("Number of clusters:", nrow(main_clusters), "\n")
cat("Estimated extremal index:", main_theta, "\n")


# ============================================================
# 12. Plot Clusters Over Time
# ============================================================

png(file.path(fig_dir, "declustered_exceedances_run_length_5.png"), width = 1100, height = 650)

plot(
  loss_dates,
  losses_num,
  type = "h",
  lwd = 0.5,
  col = "gray70",
  ylim = c(0, max(losses_num) * 1.05),
  main = "Declustered SPY Loss Exceedances, Run Length = 5",
  xlab = "Date",
  ylab = "Daily Loss"
)

abline(h = u, col = "red", lwd = 3, lty = 2)

points(
  main_exceedance_data$date,
  main_exceedance_data$loss,
  pch = 20,
  col = "gray30"
)

points(
  main_clusters$cluster_max_date,
  main_clusters$cluster_max_loss,
  pch = 19,
  col = "blue"
)

legend(
  "topright",
  legend = c("97.5% threshold", "Exceedances", "Cluster maxima"),
  col = c("red", "gray30", "blue"),
  lty = c(2, NA, NA),
  lwd = c(3, NA, NA),
  pch = c(NA, 20, 19),
  bty = "n"
)

dev.off()


# ============================================================
# 13. Plot Cluster Size Distribution
# ============================================================

png(file.path(fig_dir, "cluster_size_distribution_run_length_5.png"), width = 900, height = 600)

hist(
  main_clusters$cluster_size,
  breaks = seq(0.5, max(main_clusters$cluster_size) + 0.5, by = 1),
  main = "Cluster Size Distribution, Run Length = 5",
  xlab = "Number of Exceedances in Cluster",
  ylab = "Frequency"
)

dev.off()


# ============================================================
# 14. Plot Extremal Index by Run Length
# ============================================================

png(file.path(fig_dir, "extremal_index_by_run_length.png"), width = 900, height = 600)

plot(
  declustering_results$run_length,
  declustering_results$extremal_index,
  type = "b",
  pch = 19,
  xlab = "Run Length",
  ylab = expression(hat(theta)),
  main = "Runs Estimate of Extremal Index"
)

abline(h = 1, lty = 2, col = "red")

dev.off()


# ============================================================
# 15. Fit Chapter 4 iid GPD Model
# ============================================================

# This is the Chapter 4 threshold model using all exceedances.
# It treats exceedances as approximately independent.

capture.output(
  gpd_iid <- gpd.fit(losses_num, threshold = u, npy = 252)
)

cat("\nChapter 4 iid GPD fit:\n")
print(gpd_iid$mle)
print(gpd_iid$se)


# ============================================================
# 16. Fit GPD to Declustered Cluster Maxima
# ============================================================

# For dependent extremes, we model the cluster maxima rather than
# every exceedance.
#
# Each cluster contributes only its largest loss.

cluster_maxima <- main_clusters$cluster_max_loss

capture.output(
  gpd_cluster <- gpd.fit(cluster_maxima, threshold = u)
)

cat("\nChapter 5 declustered GPD fit:\n")
print(gpd_cluster$mle)
print(gpd_cluster$se)

gpd_fit_comparison <- data.frame(
  model = c("Chapter 4 iid exceedances", "Chapter 5 declustered cluster maxima"),
  observations_used = c(n_exceedances, length(cluster_maxima)),
  sigma = c(gpd_iid$mle[1], gpd_cluster$mle[1]),
  xi = c(gpd_iid$mle[2], gpd_cluster$mle[2]),
  se_sigma = c(gpd_iid$se[1], gpd_cluster$se[1]),
  se_xi = c(gpd_iid$se[2], gpd_cluster$se[2])
)

print(gpd_fit_comparison)

write.csv(
  gpd_fit_comparison,
  file.path(res_dir, "gpd_iid_vs_declustered_comparison.csv"),
  row.names = FALSE
)


# ============================================================
# 17. Diagnostics for Declustered GPD
# ============================================================

png(file.path(fig_dir, "gpd_diagnostics_declustered_run_length_5.png"), width = 1000, height = 800)
gpd.diag(gpd_cluster)
dev.off()


# ============================================================
# 18. Return Level Functions
# ============================================================

# Chapter 4 iid return level:
#
# x_m = u + sigma / xi * ((m * zeta_u)^xi - 1)
#
# where zeta_u is the raw exceedance rate.

gpd_return_level <- function(m, u, sigma, xi, rate) {
  
  if (abs(xi) < 1e-6) {
    x_m <- u + sigma * log(m * rate)
  } else {
    x_m <- u + (sigma / xi) * ((m * rate)^xi - 1)
  }
  
  return(x_m)
}


# ============================================================
# 19. Compute Return Levels: iid vs Declustered
# ============================================================

return_periods_days <- c(252, 1260, 2520)
return_labels <- c("1 year", "5 years", "10 years")

# Chapter 4 iid rate:
# raw exceedance probability per trading day.

raw_exceedance_rate <- n_exceedances / n_total

# Chapter 5 declustered rate:
# cluster occurrence probability per trading day.

cluster_rate <- nrow(main_clusters) / n_total

iid_return_levels <- sapply(
  return_periods_days,
  gpd_return_level,
  u = u,
  sigma = gpd_iid$mle[1],
  xi = gpd_iid$mle[2],
  rate = raw_exceedance_rate
)

declustered_return_levels <- sapply(
  return_periods_days,
  gpd_return_level,
  u = u,
  sigma = gpd_cluster$mle[1],
  xi = gpd_cluster$mle[2],
  rate = cluster_rate
)

return_level_comparison <- data.frame(
  return_period = return_labels,
  m_trading_days = return_periods_days,
  chapter4_iid_percent = 100 * iid_return_levels,
  chapter5_declustered_percent = 100 * declustered_return_levels
)

print(return_level_comparison)

write.csv(
  return_level_comparison,
  file.path(res_dir, "return_level_iid_vs_declustered.csv"),
  row.names = FALSE
)


# ============================================================
# 20. Plot Return Level Comparison
# ============================================================

png(file.path(fig_dir, "return_levels_iid_vs_declustered.png"), width = 900, height = 600)

plot(
  c(1, 5, 10),
  return_level_comparison$chapter4_iid_percent,
  type = "b",
  pch = 19,
  col = "black",
  ylim = range(
    return_level_comparison$chapter4_iid_percent,
    return_level_comparison$chapter5_declustered_percent
  ),
  xlab = "Return Period in Years",
  ylab = "Return Level (%)",
  main = "Return Levels: iid vs Declustered Threshold Models"
)

lines(
  c(1, 5, 10),
  return_level_comparison$chapter5_declustered_percent,
  type = "b",
  pch = 17,
  col = "blue"
)

legend(
  "topleft",
  legend = c("Chapter 4 iid exceedances", "Chapter 5 declustered clusters"),
  col = c("black", "blue"),
  pch = c(19, 17),
  lty = 1,
  bty = "n"
)

dev.off()


# ============================================================
# 21. Chapter 3, 4, 5 Shape Comparison
# ============================================================

# Refit Chapter 3 monthly GEV model for comparison.

monthly_max_losses <- apply.monthly(losses, max)
monthly_max_losses_num <- as.numeric(monthly_max_losses)

capture.output(
  gev_ch3 <- gev.fit(monthly_max_losses_num)
)

shape_comparison <- data.frame(
  model = c(
    "Chapter 3 GEV monthly maxima",
    "Chapter 4 iid GPD, 97.5% threshold",
    "Chapter 5 declustered GPD, run length 5"
  ),
  xi = c(
    gev_ch3$mle[3],
    gpd_iid$mle[2],
    gpd_cluster$mle[2]
  ),
  se_xi = c(
    gev_ch3$se[3],
    gpd_iid$se[2],
    gpd_cluster$se[2]
  )
)

print(shape_comparison)

write.csv(
  shape_comparison,
  file.path(res_dir, "chapter3_chapter4_chapter5_shape_comparison.csv"),
  row.names = FALSE
)

png(file.path(fig_dir, "chapter3_chapter4_chapter5_shape_comparison.png"), width = 900, height = 600)

barplot(
  shape_comparison$xi,
  names.arg = c("Ch3 GEV", "Ch4 iid GPD", "Ch5 declust."),
  ylab = expression(hat(xi)),
  main = "Shape Parameter Comparison Across Chapters",
  ylim = c(0, max(shape_comparison$xi + shape_comparison$se_xi, na.rm = TRUE) * 1.4)
)

abline(h = 0, lty = 2, col = "red")

dev.off()


# ============================================================
# 22. Save Model Summary
# ============================================================

sink(file.path(res_dir, "chapter5_model_summary.txt"))

cat("Chapter 5: Extremes of Dependent Sequences\n")
cat("==========================================\n\n")

cat("Data:\n")
cat("Asset: SPY\n")
cat("Returns: daily log returns\n")
cat("Losses: negative daily log returns\n")
cat("Number of daily observations:", n_total, "\n\n")

cat("Threshold:\n")
print(threshold_summary)
cat("\n")

cat("Declustering results:\n")
print(declustering_results)
cat("\n")

cat("Main run length:", main_run_length, "\n")
cat("Number of exceedances:", n_exceedances, "\n")
cat("Number of clusters:", nrow(main_clusters), "\n")
cat("Extremal index estimate:", main_theta, "\n\n")

cat("GPD comparison:\n")
print(gpd_fit_comparison)
cat("\n")

cat("Return level comparison:\n")
print(return_level_comparison)
cat("\n")

cat("Shape comparison across Chapters 3, 4, and 5:\n")
print(shape_comparison)
cat("\n")

cat("Interpretation:\n")
cat("Chapter 5 accounts for dependence by grouping nearby threshold exceedances into clusters.\n")
cat("The extremal index estimate is below 1 if exceedances cluster.\n")
cat("The declustered GPD model uses cluster maxima rather than all exceedances.\n")

sink()


# ============================================================
# 23. Save Session Info
# ============================================================

sink(file.path(res_dir, "session_info.txt"))
print(sessionInfo())
sink()


# ============================================================
# 24. Final Console Summary
# ============================================================

cat("\n================ CHAPTER 5 SUMMARY ================\n")
cat("Dependent extremes analysis completed.\n\n")

cat("Threshold:", round(100 * u, 3), "%\n")
cat("Number of exceedances:", n_exceedances, "\n")
cat("Main run length:", main_run_length, "\n")
cat("Number of clusters:", nrow(main_clusters), "\n")
cat("Estimated extremal index:", round(main_theta, 4), "\n\n")

cat("Declustering results:\n")
print(declustering_results)

cat("\nGPD iid vs declustered comparison:\n")
print(gpd_fit_comparison)

cat("\nReturn level comparison:\n")
print(return_level_comparison)

cat("\nShape comparison:\n")
print(shape_comparison)

cat("\nFigures saved in:\n")
cat(fig_dir, "\n")

cat("\nResults saved in:\n")
cat(res_dir, "\n")

cat("\nSaved figures:\n")
print(list.files(fig_dir))

cat("\nSaved result files:\n")
print(list.files(res_dir))

cat("===================================================\n")