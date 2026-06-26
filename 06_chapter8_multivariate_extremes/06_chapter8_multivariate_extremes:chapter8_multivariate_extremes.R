# ============================================================
# Chapter 8: Multivariate Extremes
# Joint Extreme Losses in SPY, QQQ, TLT, and VIX
# ============================================================

library(quantmod)
library(ismev)
library(xts)
library(zoo)

setwd("~/Desktop/UNI/Projects/EVT/evt-losses")
cat("Current working directory:\n")
print(getwd())

if (basename(getwd()) == "06_chapter8_multivariate_extremes") {
  chapter_dir <- "."
} else {
  chapter_dir <- "06_chapter8_multivariate_extremes"
}

fig_dir <- file.path(chapter_dir, "figures")
res_dir <- file.path(chapter_dir, "results")

if (chapter_dir != ".") {
  dir.create(chapter_dir, showWarnings = FALSE, recursive = TRUE)
}

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)


# ============================================================
# Helper Functions
# ============================================================

fit_gpd_safely <- function(x, threshold, npy = 252) {
  fit <- tryCatch(
    {
      capture.output(
        fitted_model <- gpd.fit(x, threshold = threshold, npy = npy)
      )
      fitted_model
    },
    error = function(e) {
      message("GPD fit failed: ", e$message)
      NULL
    }
  )
  
  return(fit)
}


gpd_return_level <- function(m_observations, threshold, sigma, xi, zeta) {
  if (m_observations * zeta <= 1) {
    return(NA_real_)
  }
  
  if (abs(xi) < 1e-6) {
    return(threshold + sigma * log(m_observations * zeta))
  }
  
  return(threshold + (sigma / xi) * ((m_observations * zeta)^xi - 1))
}


empirical_frechet_transform <- function(x) {
  n <- length(x)
  u <- rank(x, ties.method = "average") / (n + 1)
  z <- -1 / log(u)
  
  return(data.frame(u = u, z = z))
}


estimate_chi_diagnostics <- function(u_x, u_y, u_grid) {
  output <- data.frame(
    u = u_grid,
    joint_survival = NA_real_,
    single_survival = NA_real_,
    chi = NA_real_,
    chi_bar = NA_real_,
    n_joint = NA_integer_
  )
  
  for (i in seq_along(u_grid)) {
    level <- u_grid[i]
    
    joint <- mean(u_x > level & u_y > level)
    single <- 1 - level
    n_joint <- sum(u_x > level & u_y > level)
    
    chi_value <- joint / single
    
    if (joint > 0 && joint < 1 && single > 0 && single < 1) {
      chi_bar_value <- (2 * log(single) / log(joint)) - 1
    } else {
      chi_bar_value <- NA_real_
    }
    
    output$joint_survival[i] <- joint
    output$single_survival[i] <- single
    output$chi[i] <- chi_value
    output$chi_bar[i] <- chi_bar_value
    output$n_joint[i] <- n_joint
  }
  
  return(output)
}


estimate_implied_logistic_from_chi <- function(chi_values) {
  chi_hat <- mean(chi_values, na.rm = TRUE)
  chi_hat <- max(min(chi_hat, 1), 0)
  
  implied_alpha_hat <- log(2 - chi_hat) / log(2)
  
  return(
    data.frame(
      chi_hat = chi_hat,
      implied_alpha_hat = implied_alpha_hat
    )
  )
}


bootstrap_chi_hat_iid <- function(u_x, u_y, u_grid,
                                  high_chi_cutoff = 0.95,
                                  n_bootstrap = 500,
                                  seed = 123) {
  set.seed(seed)
  
  n <- length(u_x)
  boot_chi <- rep(NA_real_, n_bootstrap)
  high_grid <- u_grid[u_grid >= high_chi_cutoff]
  
  for (b in seq_len(n_bootstrap)) {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    
    diag_b <- estimate_chi_diagnostics(
      u_x = u_x[idx],
      u_y = u_y[idx],
      u_grid = high_grid
    )
    
    boot_chi[b] <- mean(diag_b$chi, na.rm = TRUE)
  }
  
  boot_chi <- boot_chi[is.finite(boot_chi)]
  
  if (length(boot_chi) == 0) {
    return(
      data.frame(
        bootstrap_type = "iid",
        block_length = NA_integer_,
        chi_boot_mean = NA_real_,
        chi_boot_sd = NA_real_,
        chi_ci_lower = NA_real_,
        chi_ci_upper = NA_real_,
        n_bootstrap_success = 0
      )
    )
  }
  
  return(
    data.frame(
      bootstrap_type = "iid",
      block_length = NA_integer_,
      chi_boot_mean = mean(boot_chi),
      chi_boot_sd = sd(boot_chi),
      chi_ci_lower = as.numeric(quantile(boot_chi, 0.025, names = FALSE)),
      chi_ci_upper = as.numeric(quantile(boot_chi, 0.975, names = FALSE)),
      n_bootstrap_success = length(boot_chi)
    )
  )
}


bootstrap_chi_hat_block <- function(u_x, u_y, u_grid,
                                    high_chi_cutoff = 0.95,
                                    n_bootstrap = 500,
                                    block_length = 5,
                                    seed = 123) {
  set.seed(seed)
  
  n <- length(u_x)
  high_grid <- u_grid[u_grid >= high_chi_cutoff]
  boot_chi <- rep(NA_real_, n_bootstrap)
  
  if (n <= block_length) {
    return(
      data.frame(
        bootstrap_type = "block",
        block_length = block_length,
        chi_boot_mean = NA_real_,
        chi_boot_sd = NA_real_,
        chi_ci_lower = NA_real_,
        chi_ci_upper = NA_real_,
        n_bootstrap_success = 0
      )
    )
  }
  
  possible_starts <- seq_len(n - block_length + 1)
  
  for (b in seq_len(n_bootstrap)) {
    sampled_indices <- integer(0)
    
    while (length(sampled_indices) < n) {
      start <- sample(possible_starts, size = 1)
      block <- start:(start + block_length - 1)
      sampled_indices <- c(sampled_indices, block)
    }
    
    sampled_indices <- sampled_indices[seq_len(n)]
    
    diag_b <- estimate_chi_diagnostics(
      u_x = u_x[sampled_indices],
      u_y = u_y[sampled_indices],
      u_grid = high_grid
    )
    
    boot_chi[b] <- mean(diag_b$chi, na.rm = TRUE)
  }
  
  boot_chi <- boot_chi[is.finite(boot_chi)]
  
  if (length(boot_chi) == 0) {
    return(
      data.frame(
        bootstrap_type = "block",
        block_length = block_length,
        chi_boot_mean = NA_real_,
        chi_boot_sd = NA_real_,
        chi_ci_lower = NA_real_,
        chi_ci_upper = NA_real_,
        n_bootstrap_success = 0
      )
    )
  }
  
  return(
    data.frame(
      bootstrap_type = "block",
      block_length = block_length,
      chi_boot_mean = mean(boot_chi),
      chi_boot_sd = sd(boot_chi),
      chi_ci_lower = as.numeric(quantile(boot_chi, 0.025, names = FALSE)),
      chi_ci_upper = as.numeric(quantile(boot_chi, 0.975, names = FALSE)),
      n_bootstrap_success = length(boot_chi)
    )
  )
}


summarize_gpd_fit <- function(asset_name, x, threshold_probability = 0.975,
                              npy = 252) {
  threshold <- as.numeric(quantile(x, threshold_probability, names = FALSE))
  exceedances <- x[x > threshold]
  exceedance_fraction <- mean(x > threshold)
  
  fit <- fit_gpd_safely(x, threshold = threshold, npy = npy)
  
  if (is.null(fit)) {
    return(
      list(
        fit = NULL,
        summary = data.frame(
          asset = asset_name,
          threshold_probability = threshold_probability,
          threshold = threshold,
          threshold_percent = 100 * threshold,
          n_observations = length(x),
          n_exceedances = length(exceedances),
          exceedance_fraction = exceedance_fraction,
          sigma = NA_real_,
          xi = NA_real_,
          se_sigma = NA_real_,
          se_xi = NA_real_
        )
      )
    )
  }
  
  summary <- data.frame(
    asset = asset_name,
    threshold_probability = threshold_probability,
    threshold = threshold,
    threshold_percent = 100 * threshold,
    n_observations = length(x),
    n_exceedances = length(exceedances),
    exceedance_fraction = exceedance_fraction,
    sigma = fit$mle[1],
    xi = fit$mle[2],
    se_sigma = fit$se[1],
    se_xi = fit$se[2]
  )
  
  return(list(fit = fit, summary = summary))
}


analyze_structure_variable <- function(z, structure_name,
                                       threshold_probability = 0.975,
                                       return_periods_years = c(1, 5, 10),
                                       trading_days_per_year = 252) {
  threshold <- as.numeric(quantile(z, threshold_probability, names = FALSE))
  zeta <- mean(z > threshold)
  
  fit <- fit_gpd_safely(z, threshold = threshold, npy = trading_days_per_year)
  
  if (is.null(fit)) {
    parameter_summary <- data.frame(
      structure_variable = structure_name,
      threshold_probability = threshold_probability,
      threshold = threshold,
      threshold_percent = 100 * threshold,
      n_observations = length(z),
      n_exceedances = sum(z > threshold),
      exceedance_fraction = zeta,
      sigma = NA_real_,
      xi = NA_real_,
      se_sigma = NA_real_,
      se_xi = NA_real_
    )
    
    return_levels <- data.frame(
      structure_variable = structure_name,
      return_period_years = return_periods_years,
      return_level_decimal = NA_real_,
      return_level_percent = NA_real_
    )
    
    return(
      list(
        fit = NULL,
        parameter_summary = parameter_summary,
        return_levels = return_levels
      )
    )
  }
  
  sigma <- fit$mle[1]
  xi <- fit$mle[2]
  
  return_levels_decimal <- sapply(
    return_periods_years,
    function(T_years) {
      gpd_return_level(
        m_observations = trading_days_per_year * T_years,
        threshold = threshold,
        sigma = sigma,
        xi = xi,
        zeta = zeta
      )
    }
  )
  
  parameter_summary <- data.frame(
    structure_variable = structure_name,
    threshold_probability = threshold_probability,
    threshold = threshold,
    threshold_percent = 100 * threshold,
    n_observations = length(z),
    n_exceedances = sum(z > threshold),
    exceedance_fraction = zeta,
    sigma = sigma,
    xi = xi,
    se_sigma = fit$se[1],
    se_xi = fit$se[2]
  )
  
  return_levels <- data.frame(
    structure_variable = structure_name,
    return_period_years = return_periods_years,
    return_level_decimal = as.numeric(return_levels_decimal),
    return_level_percent = 100 * as.numeric(return_levels_decimal)
  )
  
  return(
    list(
      fit = fit,
      parameter_summary = parameter_summary,
      return_levels = return_levels
    )
  )
}


compute_angular_sensitivity <- function(pair_df, pair_name,
                                        probabilities = c(0.90, 0.95, 0.975)) {
  output <- do.call(
    rbind,
    lapply(probabilities, function(p) {
      r_threshold <- as.numeric(quantile(pair_df$R, p, names = FALSE))
      angular_df <- pair_df[pair_df$R > r_threshold, ]
      
      data.frame(
        pair = pair_name,
        angular_probability = p,
        radial_threshold = r_threshold,
        n_radial_extremes = nrow(angular_df),
        mean_W = mean(angular_df$W),
        median_W = median(angular_df$W),
        sd_W = sd(angular_df$W),
        angular_balance_mean = mean(angular_df$W) - 0.5,
        angular_balance_median = median(angular_df$W) - 0.5,
        share_W_between_0_4_and_0_6 = mean(angular_df$W >= 0.4 & angular_df$W <= 0.6),
        share_W_below_0_25 = mean(angular_df$W < 0.25),
        share_W_above_0_75 = mean(angular_df$W > 0.75),
        share_endpoint_heavy = mean(angular_df$W < 0.25 | angular_df$W > 0.75),
        endpoint_asymmetry = mean(angular_df$W > 0.75) - mean(angular_df$W < 0.25)
      )
    })
  )
  
  return(output)
}


# ============================================================
# Pair Analysis Function
# ============================================================

analyze_pair <- function(data, asset_x, asset_y, pair_name,
                         fig_dir, res_dir,
                         threshold_probability = 0.975,
                         angular_probability = 0.95,
                         angular_sensitivity_probabilities = c(0.90, 0.95, 0.975),
                         chi_grid = seq(0.80, 0.975, by = 0.005),
                         high_chi_cutoff = 0.95,
                         n_bootstrap = 500,
                         block_length = 5,
                         trading_days_per_year = 252) {
  cat("\n====================================================\n")
  cat("Analyzing pair:", asset_x, "and", asset_y, "\n")
  cat("Pair name:", pair_name, "\n")
  cat("====================================================\n")
  
  x_col <- paste0(asset_x, "_extreme")
  y_col <- paste0(asset_y, "_extreme")
  
  pair_df <- data[, c("date", x_col, y_col)]
  colnames(pair_df) <- c("date", "x", "y")
  pair_df <- pair_df[complete.cases(pair_df), ]
  
  x <- pair_df$x
  y <- pair_df$y
  n <- nrow(pair_df)
  
  pair_summary <- data.frame(
    pair = pair_name,
    asset_x = asset_x,
    asset_y = asset_y,
    n_observations = n,
    start_date = min(pair_df$date),
    end_date = max(pair_df$date),
    mean_x = mean(x),
    mean_y = mean(y),
    sd_x = sd(x),
    sd_y = sd(y),
    pearson_correlation = cor(x, y),
    threshold_probability = threshold_probability
  )
  
  write.csv(
    pair_summary,
    file.path(res_dir, paste0(pair_name, "_pair_summary.csv")),
    row.names = FALSE
  )
  
  png(file.path(fig_dir, paste0(pair_name, "_extreme_direction_series.png")),
      width = 1100, height = 650)
  
  plot(
    pair_df$date,
    x,
    type = "l",
    col = "black",
    ylim = range(c(x, y), na.rm = TRUE),
    xlab = "Date",
    ylab = "Extreme-direction value",
    main = paste0(asset_x, " and ", asset_y, " Extreme-Direction Series")
  )
  
  lines(pair_df$date, y, col = "blue")
  
  legend(
    "topright",
    legend = c(asset_x, asset_y),
    col = c("black", "blue"),
    lty = 1,
    bty = "n"
  )
  
  dev.off()
  
  gpd_x <- summarize_gpd_fit(
    asset_name = asset_x,
    x = x,
    threshold_probability = threshold_probability,
    npy = trading_days_per_year
  )
  
  gpd_y <- summarize_gpd_fit(
    asset_name = asset_y,
    x = y,
    threshold_probability = threshold_probability,
    npy = trading_days_per_year
  )
  
  marginal_gpd_summary <- rbind(gpd_x$summary, gpd_y$summary)
  
  write.csv(
    marginal_gpd_summary,
    file.path(res_dir, paste0(pair_name, "_marginal_gpd_fits.csv")),
    row.names = FALSE
  )
  
  q_levels <- c(0.90, 0.95, 0.975, 0.99)
  
  joint_exceedance_table <- do.call(
    rbind,
    lapply(q_levels, function(q) {
      qx <- as.numeric(quantile(x, q, names = FALSE))
      qy <- as.numeric(quantile(y, q, names = FALSE))
      
      x_extreme <- x > qx
      y_extreme <- y > qy
      
      px <- mean(x_extreme)
      py <- mean(y_extreme)
      pxy <- mean(x_extreme & y_extreme)
      independence_prediction <- px * py
      
      ratio <- ifelse(independence_prediction > 0,
                      pxy / independence_prediction,
                      NA_real_)
      
      p_y_given_x <- ifelse(px > 0, pxy / px, NA_real_)
      p_x_given_y <- ifelse(py > 0, pxy / py, NA_real_)
      
      data.frame(
        pair = pair_name,
        quantile_level = q,
        x_threshold = qx,
        y_threshold = qy,
        x_threshold_percent = 100 * qx,
        y_threshold_percent = 100 * qy,
        p_x_exceed = px,
        p_y_exceed = py,
        p_joint_exceed = pxy,
        independence_prediction = independence_prediction,
        joint_to_independence_ratio = ratio,
        p_y_extreme_given_x_extreme = p_y_given_x,
        p_x_extreme_given_y_extreme = p_x_given_y,
        conditional_coexceedance = p_y_given_x,
        n_joint_exceedances = sum(x_extreme & y_extreme)
      )
    })
  )
  
  conditional_coexceedance_table <- joint_exceedance_table[
    , c(
      "pair",
      "quantile_level",
      "conditional_coexceedance",
      "p_y_extreme_given_x_extreme",
      "p_x_extreme_given_y_extreme",
      "n_joint_exceedances"
    )
  ]
  
  write.csv(
    joint_exceedance_table,
    file.path(res_dir, paste0(pair_name, "_joint_exceedance_table.csv")),
    row.names = FALSE
  )
  
  write.csv(
    conditional_coexceedance_table,
    file.path(res_dir, paste0(pair_name, "_conditional_coexceedance_table.csv")),
    row.names = FALSE
  )
  
  u_x <- marginal_gpd_summary$threshold[marginal_gpd_summary$asset == asset_x]
  u_y <- marginal_gpd_summary$threshold[marginal_gpd_summary$asset == asset_y]
  
  png(file.path(fig_dir, paste0(pair_name, "_scatter_thresholds.png")),
      width = 800, height = 800)
  
  plot(
    x,
    y,
    pch = 19,
    col = rgb(0, 0, 0, 0.25),
    xlab = paste0(asset_x, " extreme-direction value"),
    ylab = paste0(asset_y, " extreme-direction value"),
    main = paste0(pair_name, ": Scatter with 97.5% Thresholds")
  )
  
  abline(v = u_x, col = "red", lwd = 2, lty = 2)
  abline(h = u_y, col = "blue", lwd = 2, lty = 2)
  
  legend(
    "topleft",
    legend = c(
      paste0(asset_x, " 97.5% threshold"),
      paste0(asset_y, " 97.5% threshold")
    ),
    col = c("red", "blue"),
    lwd = 2,
    lty = 2,
    bty = "n"
  )
  
  dev.off()
  
  frechet_x <- empirical_frechet_transform(x)
  frechet_y <- empirical_frechet_transform(y)
  
  pair_df$u_x <- frechet_x$u
  pair_df$u_y <- frechet_y$u
  pair_df$z_x <- frechet_x$z
  pair_df$z_y <- frechet_y$z
  
  write.csv(
    pair_df,
    file.path(res_dir, paste0(pair_name, "_frechet_transformed_data.csv")),
    row.names = FALSE
  )
  
  chi_diagnostics <- estimate_chi_diagnostics(
    u_x = pair_df$u_x,
    u_y = pair_df$u_y,
    u_grid = chi_grid
  )
  
  chi_diagnostics$pair <- pair_name
  
  chi_diagnostics <- chi_diagnostics[
    , c("pair", "u", "joint_survival", "single_survival",
        "chi", "chi_bar", "n_joint")
  ]
  
  write.csv(
    chi_diagnostics,
    file.path(res_dir, paste0(pair_name, "_extremal_dependence_diagnostics.csv")),
    row.names = FALSE
  )
  
  high_chi_values <- chi_diagnostics$chi[chi_diagnostics$u >= high_chi_cutoff]
  implied_logistic_summary <- estimate_implied_logistic_from_chi(high_chi_values)
  
  implied_logistic_summary$pair <- pair_name
  implied_logistic_summary$high_chi_cutoff <- high_chi_cutoff
  implied_logistic_summary$n_high_threshold_points <- sum(chi_diagnostics$u >= high_chi_cutoff)
  implied_logistic_summary$model_note <- "implied_from_empirical_chi_not_full_likelihood"
  
  implied_logistic_summary <- implied_logistic_summary[
    , c("pair", "high_chi_cutoff", "n_high_threshold_points",
        "chi_hat", "implied_alpha_hat", "model_note")
  ]
  
  chi_bootstrap_iid <- bootstrap_chi_hat_iid(
    u_x = pair_df$u_x,
    u_y = pair_df$u_y,
    u_grid = chi_grid,
    high_chi_cutoff = high_chi_cutoff,
    n_bootstrap = n_bootstrap,
    seed = 123
  )
  
  chi_bootstrap_block <- bootstrap_chi_hat_block(
    u_x = pair_df$u_x,
    u_y = pair_df$u_y,
    u_grid = chi_grid,
    high_chi_cutoff = high_chi_cutoff,
    n_bootstrap = n_bootstrap,
    block_length = block_length,
    seed = 123
  )
  
  chi_bootstrap <- rbind(chi_bootstrap_iid, chi_bootstrap_block)
  chi_bootstrap$pair <- pair_name
  chi_bootstrap$high_chi_cutoff <- high_chi_cutoff
  chi_bootstrap$n_bootstrap_requested <- n_bootstrap
  
  chi_bootstrap <- chi_bootstrap[
    , c(
      "pair",
      "bootstrap_type",
      "block_length",
      "high_chi_cutoff",
      "n_bootstrap_requested",
      "n_bootstrap_success",
      "chi_boot_mean",
      "chi_boot_sd",
      "chi_ci_lower",
      "chi_ci_upper"
    )
  ]
  
  write.csv(
    implied_logistic_summary,
    file.path(res_dir, paste0(pair_name, "_implied_logistic_dependence_summary.csv")),
    row.names = FALSE
  )
  
  write.csv(
    chi_bootstrap,
    file.path(res_dir, paste0(pair_name, "_chi_bootstrap_intervals.csv")),
    row.names = FALSE
  )
  
  block_ci <- chi_bootstrap[chi_bootstrap$bootstrap_type == "block", ]
  
  png(file.path(fig_dir, paste0(pair_name, "_chi_chibar_diagnostics.png")),
      width = 1000, height = 700)
  
  par(mfrow = c(2, 1))
  
  plot(
    chi_diagnostics$u,
    chi_diagnostics$chi,
    type = "b",
    pch = 19,
    col = "black",
    ylim = c(0, 1),
    xlab = "Quantile level u",
    ylab = "chi(u)",
    main = paste0(pair_name, ": Empirical Extremal Dependence chi(u)")
  )
  
  abline(h = implied_logistic_summary$chi_hat, col = "red", lwd = 2, lty = 2)
  abline(v = high_chi_cutoff, col = "gray50", lwd = 2, lty = 3)
  
  if (nrow(block_ci) == 1 &&
      is.finite(block_ci$chi_ci_lower) &&
      is.finite(block_ci$chi_ci_upper)) {
    abline(h = block_ci$chi_ci_lower, col = "red", lwd = 1, lty = 3)
    abline(h = block_ci$chi_ci_upper, col = "red", lwd = 1, lty = 3)
  }
  
  legend(
    "topright",
    legend = c("empirical chi(u)", "high-threshold average", "block bootstrap CI"),
    col = c("black", "red", "red"),
    lty = c(1, 2, 3),
    pch = c(19, NA, NA),
    bty = "n"
  )
  
  plot(
    chi_diagnostics$u,
    chi_diagnostics$chi_bar,
    type = "b",
    pch = 19,
    col = "blue",
    ylim = c(-1, 1),
    xlab = "Quantile level u",
    ylab = "chi_bar(u)",
    main = paste0(pair_name, ": Empirical Extremal Dependence chi_bar(u)")
  )
  
  abline(h = 0, col = "gray50", lwd = 2, lty = 2)
  abline(h = 1, col = "red", lwd = 2, lty = 2)
  abline(v = high_chi_cutoff, col = "gray50", lwd = 2, lty = 3)
  
  par(mfrow = c(1, 1))
  
  dev.off()
  
  pair_df$R <- pair_df$z_x + pair_df$z_y
  pair_df$W <- pair_df$z_x / (pair_df$z_x + pair_df$z_y)
  
  r_threshold <- as.numeric(
    quantile(pair_df$R, angular_probability, names = FALSE)
  )
  
  angular_df <- pair_df[pair_df$R > r_threshold, ]
  
  angular_summary <- data.frame(
    pair = pair_name,
    angular_probability = angular_probability,
    radial_threshold = r_threshold,
    n_radial_extremes = nrow(angular_df),
    mean_W = mean(angular_df$W),
    median_W = median(angular_df$W),
    sd_W = sd(angular_df$W),
    angular_balance_mean = mean(angular_df$W) - 0.5,
    angular_balance_median = median(angular_df$W) - 0.5,
    share_W_between_0_4_and_0_6 = mean(angular_df$W >= 0.4 & angular_df$W <= 0.6),
    share_W_below_0_25 = mean(angular_df$W < 0.25),
    share_W_above_0_75 = mean(angular_df$W > 0.75),
    share_endpoint_heavy = mean(angular_df$W < 0.25 | angular_df$W > 0.75),
    endpoint_asymmetry = mean(angular_df$W > 0.75) - mean(angular_df$W < 0.25)
  )
  
  angular_sensitivity <- compute_angular_sensitivity(
    pair_df = pair_df,
    pair_name = pair_name,
    probabilities = angular_sensitivity_probabilities
  )
  
  write.csv(
    angular_summary,
    file.path(res_dir, paste0(pair_name, "_angular_summary.csv")),
    row.names = FALSE
  )
  
  write.csv(
    angular_sensitivity,
    file.path(res_dir, paste0(pair_name, "_angular_sensitivity.csv")),
    row.names = FALSE
  )
  
  write.csv(
    angular_df,
    file.path(res_dir, paste0(pair_name, "_angular_extreme_points.csv")),
    row.names = FALSE
  )
  
  z_max_frechet <- pmax(pair_df$z_x, pair_df$z_y)
  z_min_frechet <- pmin(pair_df$z_x, pair_df$z_y)
  
  frechet_structure_summary <- data.frame(
    pair = pair_name,
    max_frechet_q90 = as.numeric(quantile(z_max_frechet, 0.90, names = FALSE)),
    max_frechet_q95 = as.numeric(quantile(z_max_frechet, 0.95, names = FALSE)),
    max_frechet_q975 = as.numeric(quantile(z_max_frechet, 0.975, names = FALSE)),
    max_frechet_q99 = as.numeric(quantile(z_max_frechet, 0.99, names = FALSE)),
    min_frechet_q90 = as.numeric(quantile(z_min_frechet, 0.90, names = FALSE)),
    min_frechet_q95 = as.numeric(quantile(z_min_frechet, 0.95, names = FALSE)),
    min_frechet_q975 = as.numeric(quantile(z_min_frechet, 0.975, names = FALSE)),
    min_frechet_q99 = as.numeric(quantile(z_min_frechet, 0.99, names = FALSE))
  )
  
  write.csv(
    frechet_structure_summary,
    file.path(res_dir, paste0(pair_name, "_frechet_structure_summary.csv")),
    row.names = FALSE
  )
  
  png(file.path(fig_dir, paste0(pair_name, "_frechet_extreme_points.png")),
      width = 850, height = 800)
  
  plot(
    pair_df$z_x,
    pair_df$z_y,
    pch = 19,
    col = rgb(0, 0, 0, 0.20),
    xlab = paste0(asset_x, " standard Frechet value"),
    ylab = paste0(asset_y, " standard Frechet value"),
    main = paste0(pair_name, ": Standard Frechet Extreme Points")
  )
  
  points(
    angular_df$z_x,
    angular_df$z_y,
    pch = 19,
    col = "red"
  )
  
  legend(
    "topleft",
    legend = c("all observations", "radial extremes"),
    col = c("black", "red"),
    pch = 19,
    bty = "n"
  )
  
  dev.off()
  
  png(file.path(fig_dir, paste0(pair_name, "_angular_density.png")),
      width = 900, height = 650)
  
  hist(
    angular_df$W,
    breaks = 20,
    probability = TRUE,
    col = "gray80",
    border = "white",
    xlab = "W = Z_X / (Z_X + Z_Y)",
    main = paste0(pair_name, ": Angular Distribution of Radial Extremes")
  )
  
  lines(density(angular_df$W), col = "blue", lwd = 2)
  abline(v = 0.5, col = "red", lwd = 2, lty = 2)
  
  legend(
    "topright",
    legend = c("density", "equal contribution W = 0.5"),
    col = c("blue", "red"),
    lwd = 2,
    lty = c(1, 2),
    bty = "n"
  )
  
  dev.off()
  
  png(file.path(fig_dir, paste0(pair_name, "_angular_sensitivity.png")),
      width = 1000, height = 650)
  
  plot(
    angular_sensitivity$angular_probability,
    angular_sensitivity$share_W_between_0_4_and_0_6,
    type = "b",
    pch = 19,
    col = "black",
    ylim = c(0, 1),
    xlab = "Radial threshold probability",
    ylab = "Share of angular extremes",
    main = paste0(pair_name, ": Angular Sensitivity")
  )
  
  lines(
    angular_sensitivity$angular_probability,
    angular_sensitivity$share_endpoint_heavy,
    type = "b",
    pch = 17,
    col = "blue"
  )
  
  legend(
    "topright",
    legend = c("W between 0.4 and 0.6", "W below 0.25 or above 0.75"),
    col = c("black", "blue"),
    pch = c(19, 17),
    lty = 1,
    bty = "n"
  )
  
  dev.off()
  
  z_max <- pmax(x, y)
  z_min <- pmin(x, y)
  
  structure_max <- analyze_structure_variable(
    z = z_max,
    structure_name = paste0("max_", asset_x, "_", asset_y),
    threshold_probability = threshold_probability,
    return_periods_years = c(1, 5, 10),
    trading_days_per_year = trading_days_per_year
  )
  
  structure_min <- analyze_structure_variable(
    z = z_min,
    structure_name = paste0("min_", asset_x, "_", asset_y),
    threshold_probability = threshold_probability,
    return_periods_years = c(1, 5, 10),
    trading_days_per_year = trading_days_per_year
  )
  
  structure_parameter_summary <- rbind(
    structure_max$parameter_summary,
    structure_min$parameter_summary
  )
  
  structure_return_levels <- rbind(
    structure_max$return_levels,
    structure_min$return_levels
  )
  
  write.csv(
    structure_parameter_summary,
    file.path(res_dir, paste0(pair_name, "_structure_variable_gpd_fits.csv")),
    row.names = FALSE
  )
  
  write.csv(
    structure_return_levels,
    file.path(res_dir, paste0(pair_name, "_structure_variable_return_levels.csv")),
    row.names = FALSE
  )
  
  png(file.path(fig_dir, paste0(pair_name, "_structure_return_levels.png")),
      width = 900, height = 650)
  
  max_name <- paste0("max_", asset_x, "_", asset_y)
  min_name <- paste0("min_", asset_x, "_", asset_y)
  
  plot(
    structure_return_levels$return_period_years[
      structure_return_levels$structure_variable == max_name
    ],
    structure_return_levels$return_level_percent[
      structure_return_levels$structure_variable == max_name
    ],
    type = "b",
    pch = 19,
    col = "black",
    ylim = range(structure_return_levels$return_level_percent, na.rm = TRUE),
    xlab = "Return period in years",
    ylab = "Return level (%)",
    main = paste0(pair_name, ": Raw Structure Variable Return Levels")
  )
  
  lines(
    structure_return_levels$return_period_years[
      structure_return_levels$structure_variable == min_name
    ],
    structure_return_levels$return_level_percent[
      structure_return_levels$structure_variable == min_name
    ],
    type = "b",
    pch = 17,
    col = "blue"
  )
  
  legend(
    "topleft",
    legend = c("max extreme", "joint stress min"),
    col = c("black", "blue"),
    pch = c(19, 17),
    lty = 1,
    bty = "n"
  )
  
  dev.off()
  
  joint_975 <- joint_exceedance_table[
    joint_exceedance_table$quantile_level == 0.975,
  ]
  
  block_boot <- chi_bootstrap[chi_bootstrap$bootstrap_type == "block", ]
  
  dependence_dashboard <- data.frame(
    pair = pair_name,
    asset_x = asset_x,
    asset_y = asset_y,
    pearson_correlation = cor(x, y),
    chi_hat = implied_logistic_summary$chi_hat,
    chi_ci_lower_block = block_boot$chi_ci_lower,
    chi_ci_upper_block = block_boot$chi_ci_upper,
    implied_alpha_hat = implied_logistic_summary$implied_alpha_hat,
    joint_ratio_975 = joint_975$joint_to_independence_ratio,
    conditional_coexceedance_975 = joint_975$conditional_coexceedance,
    n_joint_975 = joint_975$n_joint_exceedances,
    angular_share_joint_0_4_0_6 = angular_summary$share_W_between_0_4_and_0_6,
    angular_share_endpoint_heavy = angular_summary$share_endpoint_heavy,
    angular_balance_mean = angular_summary$angular_balance_mean,
    angular_balance_median = angular_summary$angular_balance_median,
    endpoint_asymmetry = angular_summary$endpoint_asymmetry
  )
  
  write.csv(
    dependence_dashboard,
    file.path(res_dir, paste0(pair_name, "_dependence_dashboard.csv")),
    row.names = FALSE
  )
  
  sink(file.path(res_dir, paste0(pair_name, "_summary.txt")))
  
  cat("Chapter 8 Pair Summary\n")
  cat("======================\n\n")
  cat("Pair:", pair_name, "\n")
  cat("Asset X:", asset_x, "\n")
  cat("Asset Y:", asset_y, "\n")
  cat("Observations:", n, "\n")
  cat("Start date:", as.character(min(pair_df$date)), "\n")
  cat("End date:", as.character(max(pair_df$date)), "\n")
  cat("Pearson correlation:", cor(x, y), "\n\n")
  
  cat("Marginal GPD fits:\n")
  print(marginal_gpd_summary)
  cat("\n")
  
  cat("Joint exceedance table:\n")
  print(joint_exceedance_table)
  cat("\n")
  
  cat("Conditional co-exceedance table:\n")
  print(conditional_coexceedance_table)
  cat("\n")
  
  cat("Empirical implied logistic dependence summary:\n")
  print(implied_logistic_summary)
  cat("\n")
  
  cat("Chi bootstrap intervals:\n")
  print(chi_bootstrap)
  cat("\n")
  
  cat("Angular summary:\n")
  print(angular_summary)
  cat("\n")
  
  cat("Angular sensitivity:\n")
  print(angular_sensitivity)
  cat("\n")
  
  cat("Frechet structure summary:\n")
  print(frechet_structure_summary)
  cat("\n")
  
  cat("Raw structure variable GPD fits:\n")
  print(structure_parameter_summary)
  cat("\n")
  
  cat("Raw structure variable return levels:\n")
  print(structure_return_levels)
  cat("\n")
  
  cat("Dependence dashboard:\n")
  print(dependence_dashboard)
  cat("\n")
  
  cat("Methodological notes:\n")
  cat("1. Conditional co-exceedance probabilities are computed at equal marginal quantile levels.\n")
  cat("2. Therefore, they are interpreted as co-exceedance strength, not directional asymmetry.\n")
  cat("3. Directional structure is assessed using angular summaries based on W.\n")
  cat("4. The implied logistic alpha is derived from empirical chi and is not a full likelihood estimate.\n")
  cat("5. Block bootstrap intervals are included because financial extremes cluster.\n")
  cat("6. Frechet structure summaries are better for cross-pair comparison than raw SPY-VIX structure levels.\n")
  
  sink()
  
  return(
    list(
      pair_summary = pair_summary,
      marginal_gpd_summary = marginal_gpd_summary,
      joint_exceedance_table = joint_exceedance_table,
      conditional_coexceedance_table = conditional_coexceedance_table,
      chi_diagnostics = chi_diagnostics,
      implied_logistic_summary = implied_logistic_summary,
      chi_bootstrap = chi_bootstrap,
      angular_summary = angular_summary,
      angular_sensitivity = angular_sensitivity,
      frechet_structure_summary = frechet_structure_summary,
      structure_parameter_summary = structure_parameter_summary,
      structure_return_levels = structure_return_levels,
      dependence_dashboard = dependence_dashboard
    )
  )
}


# ============================================================
# Download Price Data
# ============================================================

symbols_yahoo <- c("SPY", "QQQ", "TLT", "^VIX")
data_env <- new.env()

getSymbols(
  Symbols = symbols_yahoo,
  src = "yahoo",
  from = "1993-01-01",
  env = data_env,
  auto.assign = TRUE
)

spy_prices <- Ad(get("SPY", envir = data_env))
qqq_prices <- Ad(get("QQQ", envir = data_env))
tlt_prices <- Ad(get("TLT", envir = data_env))

if (exists("VIX", envir = data_env)) {
  vix_object <- get("VIX", envir = data_env)
} else if (exists("^VIX", envir = data_env)) {
  vix_object <- get("^VIX", envir = data_env)
} else {
  stop("Could not find VIX object after Yahoo download.")
}

vix_prices <- Ad(vix_object)

adjusted_prices <- merge(spy_prices, qqq_prices, join = "inner")
adjusted_prices <- merge(adjusted_prices, tlt_prices, join = "inner")
adjusted_prices <- merge(adjusted_prices, vix_prices, join = "inner")

colnames(adjusted_prices) <- c("SPY", "QQQ", "TLT", "VIX")
adjusted_prices <- na.omit(adjusted_prices)

cat("\nAdjusted price/index data:\n")
print(head(adjusted_prices))
print(tail(adjusted_prices))

adjusted_prices_df <- data.frame(
  date = index(adjusted_prices),
  coredata(adjusted_prices)
)

write.csv(
  adjusted_prices_df,
  file.path(res_dir, "adjusted_prices_and_vix.csv"),
  row.names = FALSE
)

png(file.path(fig_dir, "adjusted_prices_and_vix.png"), width = 1100, height = 650)

plot(
  adjusted_prices_df$date,
  adjusted_prices_df$SPY,
  type = "l",
  col = "black",
  ylim = range(adjusted_prices_df[, c("SPY", "QQQ", "TLT", "VIX")], na.rm = TRUE),
  xlab = "Date",
  ylab = "Level",
  main = "Adjusted Prices and VIX Index"
)

lines(adjusted_prices_df$date, adjusted_prices_df$QQQ, col = "blue")
lines(adjusted_prices_df$date, adjusted_prices_df$TLT, col = "red")
lines(adjusted_prices_df$date, adjusted_prices_df$VIX, col = "darkgreen")

legend(
  "topleft",
  legend = c("SPY", "QQQ", "TLT", "VIX"),
  col = c("black", "blue", "red", "darkgreen"),
  lty = 1,
  bty = "n"
)

dev.off()

# Indexed plot. This is better for visual comparison because all series
# start at 100. The raw-level plot is kept for transparency, but this one
# is better for the README.

indexed_prices_df <- adjusted_prices_df
indexed_prices_df$SPY <- 100 * indexed_prices_df$SPY / indexed_prices_df$SPY[1]
indexed_prices_df$QQQ <- 100 * indexed_prices_df$QQQ / indexed_prices_df$QQQ[1]
indexed_prices_df$TLT <- 100 * indexed_prices_df$TLT / indexed_prices_df$TLT[1]
indexed_prices_df$VIX <- 100 * indexed_prices_df$VIX / indexed_prices_df$VIX[1]

write.csv(
  indexed_prices_df,
  file.path(res_dir, "indexed_prices_and_vix.csv"),
  row.names = FALSE
)

png(file.path(fig_dir, "indexed_prices_and_vix.png"), width = 1100, height = 650)

plot(
  indexed_prices_df$date,
  indexed_prices_df$SPY,
  type = "l",
  col = "black",
  ylim = range(indexed_prices_df[, c("SPY", "QQQ", "TLT", "VIX")], na.rm = TRUE),
  xlab = "Date",
  ylab = "Index, first observation = 100",
  main = "Indexed Prices and VIX"
)

lines(indexed_prices_df$date, indexed_prices_df$QQQ, col = "blue")
lines(indexed_prices_df$date, indexed_prices_df$TLT, col = "red")
lines(indexed_prices_df$date, indexed_prices_df$VIX, col = "darkgreen")

legend(
  "topleft",
  legend = c("SPY", "QQQ", "TLT", "VIX"),
  col = c("black", "blue", "red", "darkgreen"),
  lty = 1,
  bty = "n"
)

dev.off()


# ============================================================
# Compute Daily Log Returns and Extreme-Direction Variables
# ============================================================

returns <- diff(log(adjusted_prices))
returns <- na.omit(returns)

spy_extreme <- -returns$SPY
qqq_extreme <- -returns$QQQ
tlt_extreme <- -returns$TLT
vix_extreme <- returns$VIX

colnames(spy_extreme) <- "SPY_extreme"
colnames(qqq_extreme) <- "QQQ_extreme"
colnames(tlt_extreme) <- "TLT_extreme"
colnames(vix_extreme) <- "VIX_extreme"

extreme_data_xts <- merge(spy_extreme, qqq_extreme, join = "inner")
extreme_data_xts <- merge(extreme_data_xts, tlt_extreme, join = "inner")
extreme_data_xts <- merge(extreme_data_xts, vix_extreme, join = "inner")
extreme_data_xts <- na.omit(extreme_data_xts)

extreme_data <- data.frame(
  date = index(extreme_data_xts),
  coredata(extreme_data_xts)
)

cat("\nExtreme-direction data summary:\n")
print(summary(extreme_data))

write.csv(
  extreme_data,
  file.path(res_dir, "daily_extreme_direction_data.csv"),
  row.names = FALSE
)

png(file.path(fig_dir, "daily_extreme_direction_all_series.png"),
    width = 1100, height = 650)

plot(
  extreme_data$date,
  extreme_data$SPY_extreme,
  type = "l",
  col = "black",
  ylim = range(extreme_data[, c("SPY_extreme", "QQQ_extreme", "TLT_extreme", "VIX_extreme")],
               na.rm = TRUE),
  xlab = "Date",
  ylab = "Extreme-direction daily log return",
  main = "Extreme-Direction Series: SPY, QQQ, TLT, and VIX"
)

lines(extreme_data$date, extreme_data$QQQ_extreme, col = "blue")
lines(extreme_data$date, extreme_data$TLT_extreme, col = "red")
lines(extreme_data$date, extreme_data$VIX_extreme, col = "darkgreen")

legend(
  "topright",
  legend = c("SPY loss", "QQQ loss", "TLT loss", "VIX increase"),
  col = c("black", "blue", "red", "darkgreen"),
  lty = 1,
  bty = "n"
)

dev.off()


# ============================================================
# Run Pair Analyses
# ============================================================

result_spy_qqq <- analyze_pair(
  data = extreme_data,
  asset_x = "SPY",
  asset_y = "QQQ",
  pair_name = "spy_qqq",
  fig_dir = fig_dir,
  res_dir = res_dir
)

result_spy_tlt <- analyze_pair(
  data = extreme_data,
  asset_x = "SPY",
  asset_y = "TLT",
  pair_name = "spy_tlt",
  fig_dir = fig_dir,
  res_dir = res_dir
)

result_spy_vix <- analyze_pair(
  data = extreme_data,
  asset_x = "SPY",
  asset_y = "VIX",
  pair_name = "spy_vix",
  fig_dir = fig_dir,
  res_dir = res_dir
)


# ============================================================
# Combined Cross-Pair Results
# ============================================================

combined_pair_summary <- rbind(
  result_spy_qqq$pair_summary,
  result_spy_tlt$pair_summary,
  result_spy_vix$pair_summary
)

combined_marginal_gpd <- rbind(
  result_spy_qqq$marginal_gpd_summary,
  result_spy_tlt$marginal_gpd_summary,
  result_spy_vix$marginal_gpd_summary
)

combined_joint_exceedance <- rbind(
  result_spy_qqq$joint_exceedance_table,
  result_spy_tlt$joint_exceedance_table,
  result_spy_vix$joint_exceedance_table
)

combined_conditional_coexceedance <- rbind(
  result_spy_qqq$conditional_coexceedance_table,
  result_spy_tlt$conditional_coexceedance_table,
  result_spy_vix$conditional_coexceedance_table
)

combined_implied_logistic_summary <- rbind(
  result_spy_qqq$implied_logistic_summary,
  result_spy_tlt$implied_logistic_summary,
  result_spy_vix$implied_logistic_summary
)

combined_chi_bootstrap <- rbind(
  result_spy_qqq$chi_bootstrap,
  result_spy_tlt$chi_bootstrap,
  result_spy_vix$chi_bootstrap
)

combined_angular_summary <- rbind(
  result_spy_qqq$angular_summary,
  result_spy_tlt$angular_summary,
  result_spy_vix$angular_summary
)

combined_angular_sensitivity <- rbind(
  result_spy_qqq$angular_sensitivity,
  result_spy_tlt$angular_sensitivity,
  result_spy_vix$angular_sensitivity
)

combined_frechet_structure_summary <- rbind(
  result_spy_qqq$frechet_structure_summary,
  result_spy_tlt$frechet_structure_summary,
  result_spy_vix$frechet_structure_summary
)

combined_structure_return_levels <- rbind(
  result_spy_qqq$structure_return_levels,
  result_spy_tlt$structure_return_levels,
  result_spy_vix$structure_return_levels
)

combined_dependence_dashboard <- rbind(
  result_spy_qqq$dependence_dashboard,
  result_spy_tlt$dependence_dashboard,
  result_spy_vix$dependence_dashboard
)

combined_tail_dependence_ranking <- combined_dependence_dashboard[
  order(-combined_dependence_dashboard$chi_hat),
]

combined_tail_dependence_ranking$rank_by_chi <- seq_len(
  nrow(combined_tail_dependence_ranking)
)

write.csv(combined_pair_summary, file.path(res_dir, "combined_pair_summary.csv"), row.names = FALSE)
write.csv(combined_marginal_gpd, file.path(res_dir, "combined_marginal_gpd_fits.csv"), row.names = FALSE)
write.csv(combined_joint_exceedance, file.path(res_dir, "combined_joint_exceedance_table.csv"), row.names = FALSE)
write.csv(combined_conditional_coexceedance, file.path(res_dir, "combined_conditional_coexceedance_table.csv"), row.names = FALSE)
write.csv(combined_implied_logistic_summary, file.path(res_dir, "combined_implied_logistic_dependence_summary.csv"), row.names = FALSE)
write.csv(combined_chi_bootstrap, file.path(res_dir, "combined_chi_bootstrap_intervals.csv"), row.names = FALSE)
write.csv(combined_angular_summary, file.path(res_dir, "combined_angular_summary.csv"), row.names = FALSE)
write.csv(combined_angular_sensitivity, file.path(res_dir, "combined_angular_sensitivity.csv"), row.names = FALSE)
write.csv(combined_frechet_structure_summary, file.path(res_dir, "combined_frechet_structure_summary.csv"), row.names = FALSE)
write.csv(combined_structure_return_levels, file.path(res_dir, "combined_structure_variable_return_levels.csv"), row.names = FALSE)
write.csv(combined_dependence_dashboard, file.path(res_dir, "combined_dependence_dashboard.csv"), row.names = FALSE)
write.csv(combined_tail_dependence_ranking, file.path(res_dir, "combined_tail_dependence_ranking.csv"), row.names = FALSE)


# ============================================================
# Cross-Pair Figures
# ============================================================

png(file.path(fig_dir, "comparison_chi_all_pairs.png"), width = 1000, height = 650)

plot(
  result_spy_qqq$chi_diagnostics$u,
  result_spy_qqq$chi_diagnostics$chi,
  type = "b",
  pch = 19,
  col = "black",
  ylim = c(0, 1),
  xlab = "Quantile level u",
  ylab = "chi(u)",
  main = "Empirical Extremal Dependence chi(u)"
)

lines(result_spy_tlt$chi_diagnostics$u, result_spy_tlt$chi_diagnostics$chi,
      type = "b", pch = 17, col = "blue")

lines(result_spy_vix$chi_diagnostics$u, result_spy_vix$chi_diagnostics$chi,
      type = "b", pch = 15, col = "red")

legend(
  "topright",
  legend = c("SPY-QQQ", "SPY-TLT", "SPY-VIX"),
  col = c("black", "blue", "red"),
  pch = c(19, 17, 15),
  lty = 1,
  bty = "n"
)

dev.off()


png(file.path(fig_dir, "comparison_implied_logistic_alpha.png"), width = 850, height = 600)

barplot(
  combined_implied_logistic_summary$implied_alpha_hat,
  names.arg = combined_implied_logistic_summary$pair,
  ylab = "Implied logistic alpha",
  main = "Implied Logistic Dependence Parameter",
  ylim = c(0, 1)
)

abline(h = 1, col = "red", lty = 2, lwd = 2)

dev.off()


png(file.path(fig_dir, "comparison_chi_hat_block_bootstrap.png"), width = 900, height = 650)

bar_centers <- barplot(
  combined_dependence_dashboard$chi_hat,
  names.arg = combined_dependence_dashboard$pair,
  ylab = "Empirical chi_hat",
  main = "Empirical Extremal Dependence Coefficient with Block Bootstrap CI",
  ylim = c(0, 1)
)

arrows(
  x0 = bar_centers,
  y0 = combined_dependence_dashboard$chi_ci_lower_block,
  x1 = bar_centers,
  y1 = combined_dependence_dashboard$chi_ci_upper_block,
  angle = 90,
  code = 3,
  length = 0.06
)

abline(h = 0, col = "red", lty = 2, lwd = 2)

dev.off()


png(file.path(fig_dir, "comparison_joint_exceedance_ratios.png"), width = 1000, height = 650)

ratio_subset <- combined_joint_exceedance

plot(
  ratio_subset$quantile_level[ratio_subset$pair == "spy_qqq"],
  ratio_subset$joint_to_independence_ratio[ratio_subset$pair == "spy_qqq"],
  type = "b",
  pch = 19,
  col = "black",
  ylim = range(ratio_subset$joint_to_independence_ratio, na.rm = TRUE),
  xlab = "Quantile level",
  ylab = "Joint / independence probability ratio",
  main = "Joint Exceedance Amplification vs Independence"
)

lines(
  ratio_subset$quantile_level[ratio_subset$pair == "spy_tlt"],
  ratio_subset$joint_to_independence_ratio[ratio_subset$pair == "spy_tlt"],
  type = "b",
  pch = 17,
  col = "blue"
)

lines(
  ratio_subset$quantile_level[ratio_subset$pair == "spy_vix"],
  ratio_subset$joint_to_independence_ratio[ratio_subset$pair == "spy_vix"],
  type = "b",
  pch = 15,
  col = "red"
)

abline(h = 1, col = "gray50", lty = 2, lwd = 2)

legend(
  "topleft",
  legend = c("SPY-QQQ", "SPY-TLT", "SPY-VIX"),
  col = c("black", "blue", "red"),
  pch = c(19, 17, 15),
  lty = 1,
  bty = "n"
)

dev.off()


coex_975 <- combined_conditional_coexceedance[
  combined_conditional_coexceedance$quantile_level == 0.975,
]

png(file.path(fig_dir, "comparison_conditional_coexceedance_975.png"),
    width = 950, height = 650)

bar_positions <- barplot(
  coex_975$conditional_coexceedance,
  names.arg = coex_975$pair,
  ylim = c(0, 1),
  ylab = "Conditional co-exceedance probability",
  main = "Conditional Co-Exceedance at 97.5%"
)

text(
  x = bar_positions,
  y = coex_975$conditional_coexceedance + 0.04,
  labels = paste0("n = ", coex_975$n_joint_exceedances)
)

dev.off()


png(file.path(fig_dir, "comparison_angular_sensitivity_joint_share.png"),
    width = 1000, height = 650)

plot(
  combined_angular_sensitivity$angular_probability[
    combined_angular_sensitivity$pair == "spy_qqq"
  ],
  combined_angular_sensitivity$share_W_between_0_4_and_0_6[
    combined_angular_sensitivity$pair == "spy_qqq"
  ],
  type = "b",
  pch = 19,
  col = "black",
  ylim = c(0, 1),
  xlab = "Radial threshold probability",
  ylab = "Share with W between 0.4 and 0.6",
  main = "Angular Joint-Contribution Sensitivity"
)

lines(
  combined_angular_sensitivity$angular_probability[
    combined_angular_sensitivity$pair == "spy_tlt"
  ],
  combined_angular_sensitivity$share_W_between_0_4_and_0_6[
    combined_angular_sensitivity$pair == "spy_tlt"
  ],
  type = "b",
  pch = 17,
  col = "blue"
)

lines(
  combined_angular_sensitivity$angular_probability[
    combined_angular_sensitivity$pair == "spy_vix"
  ],
  combined_angular_sensitivity$share_W_between_0_4_and_0_6[
    combined_angular_sensitivity$pair == "spy_vix"
  ],
  type = "b",
  pch = 15,
  col = "red"
)

legend(
  "topright",
  legend = c("SPY-QQQ", "SPY-TLT", "SPY-VIX"),
  col = c("black", "blue", "red"),
  pch = c(19, 17, 15),
  lty = 1,
  bty = "n"
)

dev.off()


png(file.path(fig_dir, "comparison_endpoint_asymmetry.png"), width = 900, height = 650)

barplot(
  combined_angular_summary$endpoint_asymmetry,
  names.arg = combined_angular_summary$pair,
  ylab = "Endpoint asymmetry: share(W > 0.75) - share(W < 0.25)",
  main = "Angular Endpoint Asymmetry"
)

abline(h = 0, col = "red", lty = 2, lwd = 2)

dev.off()


# ============================================================
# Chapter Summary
# ============================================================

sink(file.path(res_dir, "chapter8_model_summary.txt"))

cat("Chapter 8: Multivariate Extremes\n")
cat("================================\n\n")

cat("Objective:\n")
cat("This chapter studies joint extreme events using multivariate EVT.\n")
cat("The main pair is SPY-QQQ, representing equity-equity joint crash risk.\n")
cat("The contrast pair is SPY-TLT, representing equity-bond joint stress behavior.\n")
cat("The stress-volatility pair is SPY-VIX, representing equity losses and volatility spikes.\n\n")

cat("Data:\n")
cat("Series: SPY, QQQ, TLT, VIX\n")
cat("Start date:", as.character(min(extreme_data$date)), "\n")
cat("End date:", as.character(max(extreme_data$date)), "\n")
cat("Observations:", nrow(extreme_data), "\n\n")

cat("Variable construction:\n")
cat("SPY_extreme = - SPY daily log return\n")
cat("QQQ_extreme = - QQQ daily log return\n")
cat("TLT_extreme = - TLT daily log return\n")
cat("VIX_extreme = VIX daily log return\n")
cat("All variables are oriented so that large positive values are stress events.\n\n")

cat("Combined pair summary:\n")
print(combined_pair_summary)
cat("\n")

cat("Combined marginal GPD fits:\n")
print(combined_marginal_gpd)
cat("\n")

cat("Combined joint exceedance table:\n")
print(combined_joint_exceedance)
cat("\n")

cat("Combined conditional co-exceedance table:\n")
print(combined_conditional_coexceedance)
cat("\n")

cat("Combined implied logistic dependence summary:\n")
print(combined_implied_logistic_summary)
cat("\n")

cat("Combined chi bootstrap intervals:\n")
print(combined_chi_bootstrap)
cat("\n")

cat("Combined angular summary:\n")
print(combined_angular_summary)
cat("\n")

cat("Combined angular sensitivity:\n")
print(combined_angular_sensitivity)
cat("\n")

cat("Combined Frechet structure summary:\n")
print(combined_frechet_structure_summary)
cat("\n")

cat("Combined dependence dashboard:\n")
print(combined_dependence_dashboard)
cat("\n")

cat("Tail dependence ranking:\n")
print(combined_tail_dependence_ranking)
cat("\n")

cat("Combined raw structure variable return levels:\n")
print(combined_structure_return_levels)
cat("\n")

cat("Methodological notes:\n")
cat("1. Conditional co-exceedance probabilities are computed at equal marginal quantile levels.\n")
cat("2. They are interpreted as co-exceedance strength, not directional asymmetry.\n")
cat("3. Directional structure is assessed through angular summaries based on W.\n")
cat("4. The implied logistic alpha is derived from empirical chi and is not a full likelihood estimate.\n")
cat("5. Both iid and block bootstrap intervals are reported; block bootstrap is preferred because extremes cluster.\n")
cat("6. Raw structure-variable return levels are interpreted within each pair's units.\n")
cat("7. Frechet-scale structure summaries are better for cross-pair comparison, especially for SPY-VIX.\n")
cat("8. Results are empirical evidence over the observed tail range, not definitive proof of asymptotic dependence.\n")
cat("9. The 99% threshold results are treated as robustness checks because joint counts can be small.\n")
cat("10. The main reported threshold is 97.5%, balancing tail focus and sample size.\n")

sink()


sink(file.path(res_dir, "session_info.txt"))
print(sessionInfo())
sink()


# ============================================================
# Final Console Summary
# ============================================================

cat("\n================ CHAPTER 8 SUMMARY ================\n")
cat("Multivariate EVT analysis completed.\n\n")

cat("Pairs analyzed:\n")
cat("1. SPY-QQQ: equity-equity joint extremes\n")
cat("2. SPY-TLT: equity-bond joint stress\n")
cat("3. SPY-VIX: equity crash vs volatility spike\n\n")

cat("Combined pair summary:\n")
print(combined_pair_summary)

cat("\nEmpirical implied logistic dependence summary:\n")
print(combined_implied_logistic_summary)

cat("\nChi bootstrap intervals:\n")
print(combined_chi_bootstrap)

cat("\nAngular summary:\n")
print(combined_angular_summary)

cat("\nDependence dashboard:\n")
print(combined_dependence_dashboard)

cat("\nTail dependence ranking:\n")
print(combined_tail_dependence_ranking)

cat("\nFrechet structure summary:\n")
print(combined_frechet_structure_summary)

cat("\nRaw structure variable return levels:\n")
print(combined_structure_return_levels)

cat("\nFigures saved in:\n")
cat(fig_dir, "\n")

cat("\nResults saved in:\n")
cat(res_dir, "\n")

cat("\nSaved figures:\n")
print(list.files(fig_dir))

cat("\nSaved result files:\n")
print(list.files(res_dir))

cat("===================================================\n")