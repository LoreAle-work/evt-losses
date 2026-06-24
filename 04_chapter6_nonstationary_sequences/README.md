# Chapter 6: Extremes of Non-Stationary Sequences

This chapter applies the ideas from Chapter 6 of Coles' *An Introduction to Statistical Modeling of Extreme Values* to daily SPY losses.

The goal is to move beyond stationary EVT models and allow extreme loss behavior to depend on market conditions.

In Chapters 3, 4, and 5, the tail model parameters were treated as constant. Chapter 6 relaxes that assumption.

Financial markets are not stationary. Volatility changes, crisis regimes appear, liquidity conditions shift, and the probability of extreme losses varies over time. A stationary model forces one fixed tail model across both calm and crisis periods, which is convenient but often unrealistic.

## Objective

The objective is to study whether extreme SPY losses change with market volatility.

The main questions are:

- Does the probability of exceeding a high loss threshold depend on volatility?
- Does the severity of threshold exceedances depend on volatility?
- Does a non-stationary GPD model fit better than a stationary GPD model?
- Do return levels change between low-volatility and crisis-volatility states?
- Does the volatility effect survive after accounting for clustering?
- Is the heavy-tail conclusion from Chapters 3, 4, and 5 still present after accounting for volatility?

## Data

The analysis uses daily adjusted closing prices for SPY.

Let P<sub>t</sub> denote the adjusted closing price on day t. Daily log returns are computed as:

$$
R_t = \log(P_t) - \log(P_{t-1})
$$

Losses are defined as:

$$
L_t = -R_t
$$

This transformation turns large negative returns into large positive losses.

After constructing the lagged volatility covariate, the usable dataset contains **8384 daily observations**.

## Covariate: Lagged Realized Volatility

Chapter 6 allows model parameters to depend on covariates.

In this project, the covariate is **lagged 21-day realized volatility**.

For each day t, volatility is estimated using the standard deviation of daily returns over the previous 21 trading days.

The covariate is lagged by one day to avoid look-ahead bias. This means that today's extreme loss is modeled using volatility information that would have been available before today.

The volatility covariate is then log-transformed and standardized:

```text
raw volatility -> log volatility -> standardized log volatility
```

The standardized volatility covariate is denoted z<sub>t</sub>.

## Threshold

The threshold is the 97.5% empirical quantile of daily losses.

An exceedance occurs when:

$$
L_t > u
$$

where u is the threshold.

The excess above the threshold is:

$$
Y_t = L_t - u \mid L_t > u
$$

The threshold summary is:

| Quantity | Value |
|---|---:|
| Threshold probability | 97.5% |
| Threshold value | 0.0242 |
| Threshold in percent | 2.42% |
| Usable observations | 8384 |
| Number of exceedances | 210 |
| Exceedance fraction | 2.50% |

This is the same threshold modeling framework used in Chapter 4, but Chapter 6 allows the model to change with volatility.

## Stationary GPD Model

The stationary threshold model assumes that exceedance sizes follow a GPD with constant parameters:

$$
Y_t \sim GPD(\sigma, \xi)
$$

where:

- sigma is the scale parameter
- xi is the shape parameter

In this model, the tail behavior is constant over time.

This is essentially the Chapter 4 model.

## Non-Stationary GPD Model

The main Chapter 6 model allows the GPD scale parameter to depend on volatility.

The model is:

$$
Y_t \sim GPD(\sigma_t, \xi)
$$

with:

$$
\log(\sigma_t) = \beta_0 + \beta_1 z_t
$$

where:

- sigma<sub>t</sub> is the time-varying scale parameter
- z<sub>t</sub> is standardized lagged log-volatility
- beta<sub>0</sub> is the intercept
- beta<sub>1</sub> measures the effect of volatility
- xi is the shape parameter

The logarithmic link ensures that:

$$
\sigma_t > 0
$$

If beta<sub>1</sub> is positive, high-volatility periods are associated with larger exceedance scale.

## Volatility + Time Model

A second non-stationary model also includes a time trend:

$$
\log(\sigma_t) = \beta_0 + \beta_1 z_t + \beta_2 t
$$

This tests whether there is additional time variation in exceedance severity beyond volatility.

The three GPD models compared are:

| Model | Scale structure | Shape |
|---|---|---|
| Stationary GPD | constant sigma | constant xi |
| Volatility GPD | sigma depends on volatility | constant xi |
| Volatility + Time GPD | sigma depends on volatility and time | constant xi |

## Exceedance Probability Model

The GPD models the size of exceedances conditional on exceeding the threshold.

For return levels, the model also needs the probability of exceeding the threshold.

This chapter models threshold exceedance probability using logistic regression:

$$
P(L_t > u \mid z_t)
$$

The main logistic model is:

$$
\log \left( \frac{P(L_t > u)}{1 - P(L_t > u)} \right)
=
\alpha_0 + \alpha_1 z_t
$$

If alpha<sub>1</sub> is positive, high-volatility periods have a higher probability of producing threshold exceedances.

This separates two effects:

| Component | What it models |
|---|---|
| Logistic exceedance model | Probability of crossing the threshold |
| GPD excess model | Size of the loss after crossing the threshold |

Together, these describe how extreme risk changes with volatility.

## Exceedance Probability Results

The logistic model comparison is:

| Model | Parameters | AIC | BIC |
|---|---:|---:|---:|
| Constant exceedance probability | 1 | 1965.22 | 1972.26 |
| Volatility-dependent exceedance probability | 2 | 1699.60 | 1713.66 |
| Volatility + time exceedance probability | 3 | 1701.37 | 1722.47 |

The volatility-dependent model strongly improves fit relative to the constant exceedance probability model.

The estimated volatility coefficient is:

$$
\hat{\alpha}_1 = 1.0096
$$

with standard error:

$$
0.0613
$$

The corresponding odds multiplier is:

$$
e^{1.0096} \approx 2.74
$$

Thus, a one-standard-deviation increase in lagged log-volatility multiplies the odds of a threshold exceedance by about **2.74**.

Adding a time trend does not improve the model. Both AIC and BIC prefer the volatility-only logistic model over the volatility + time logistic model.

## Likelihood Estimation

The GPD parameters are estimated by maximum likelihood.

For exceedances y<sub>i</sub>, the GPD negative log-likelihood is minimized numerically.

For xi not equal to zero, the contribution of each exceedance is based on:

$$
\log(\sigma_i)
+
\left(1 + \frac{1}{\xi}\right)
\log\left(
1 + \xi \frac{y_i}{\sigma_i}
\right)
$$

where sigma<sub>i</sub> may depend on covariates.

The support condition is:

$$
1 + \xi \frac{y_i}{\sigma_i} > 0
$$

The script estimates the models using `optim` in R.

## GPD Model Comparison: All Exceedances

The stationary GPD model is compared with two non-stationary alternatives using all threshold exceedances.

| Model | Parameters | Observations | AIC | BIC | xi |
|---|---:|---:|---:|---:|---:|
| Stationary GPD, all exceedances | 2 | 210 | -1462.23 | -1455.54 | 0.2840 |
| Volatility-dependent GPD, all exceedances | 3 | 210 | -1492.76 | -1482.72 | 0.1208 |
| Volatility + time GPD, all exceedances | 4 | 210 | -1493.62 | -1480.23 | 0.1299 |

The volatility-dependent scale model substantially improves over the stationary model.

The likelihood ratio test comparing the stationary GPD with the volatility-dependent GPD gives:

$$
LR = 32.53
$$

with p-value:

$$
p \approx 1.17 \times 10^{-8}
$$

This provides strong evidence that exceedance severity changes with volatility.

The comparison between the volatility-only model and the volatility + time model gives:

$$
p \approx 0.091
$$

so the additional time trend is not significant at the 5% level.

The volatility + time model has a slightly lower AIC, but BIC favors the simpler volatility-only model. Since the likelihood ratio test does not support the time term at the 5% level, the main Chapter 6 specification is the volatility-dependent scale GPD:

$$
\log(\sigma_t) = \beta_0 + \beta_1 z_t
$$

## Volatility Effect on GPD Scale

The estimated volatility coefficient in the GPD scale model is:

$$
\hat{\beta}_1 = 0.3771
$$

with standard error:

$$
0.0628
$$

Since the coefficient is positive, the fitted GPD scale parameter increases with volatility.

The scale multiplier for a one-standard-deviation increase in lagged log-volatility is:

$$
e^{0.3771} \approx 1.46
$$

So a one-standard-deviation increase in volatility raises the GPD scale parameter by about **46%**.

This means that high-volatility periods are associated not only with more frequent threshold exceedances, but also with more severe exceedances.

## Conditional Return Levels

In stationary models, return levels are fixed.

In non-stationary models, return levels can depend on covariates.

This chapter estimates return levels under five volatility states:

| Volatility state | Definition |
|---|---|
| Low volatility | 25th percentile of volatility |
| Median volatility | 50th percentile of volatility |
| High volatility | 75th percentile of volatility |
| Crisis volatility | 90th percentile of volatility |
| Extreme crisis volatility | 95th percentile of volatility |

The 90th and 95th percentile scenarios are included because the 75th percentile is high, but not necessarily crisis-level high. Finance apparently needs several layers of “bad,” because one level of disaster was too simple.

For a return period m, the conditional return level is:

$$
x_m(z) =
u +
\frac{\sigma(z)}{\xi}
\left[
(m \zeta(z))^\xi - 1
\right]
$$

where:

- sigma(z) is the volatility-dependent GPD scale
- zeta(z) is the volatility-dependent threshold exceedance probability
- xi is the shape parameter
- m is the return period in trading days

The return periods are:

| Horizon | Trading days |
|---|---:|
| 1 year | 252 |
| 5 years | 1260 |
| 10 years | 2520 |

## Volatility Scenarios

The fitted volatility scenarios are:

| Volatility state | Quantile | z-vol | Exceedance probability | GPD scale |
|---|---:|---:|---:|---:|
| Low volatility | 25% | -0.7183 | 0.0072 | 0.00445 |
| Median volatility | 50% | -0.0719 | 0.0137 | 0.00568 |
| High volatility | 75% | 0.6575 | 0.0283 | 0.00748 |
| Crisis volatility | 90% | 1.3019 | 0.0528 | 0.00954 |
| Extreme crisis volatility | 95% | 1.6520 | 0.0736 | 0.01089 |

Both the exceedance probability and the fitted GPD scale increase strongly with volatility.

## Conditional Return Level Results: All Exceedances

Using the volatility-dependent exceedance probability and volatility-dependent GPD scale model, the estimated return levels are:

| Volatility state | 1-year | 5-year | 10-year |
|---|---:|---:|---:|
| Low volatility | 2.70% | 3.55% | 3.97% |
| Median volatility | 3.18% | 4.35% | 4.93% |
| High volatility | 4.08% | 5.76% | 6.60% |
| Crisis volatility | 5.32% | 7.64% | 8.78% |
| Extreme crisis volatility | 6.23% | 8.98% | 10.34% |

The return levels increase sharply with volatility.

The 10-year conditional return level rises from about **3.97%** in the low-volatility state to about **10.34%** in the extreme crisis-volatility state.

This shows that extreme loss risk is not constant over time. It depends strongly on market conditions.

## Stationary vs Conditional Return Levels

The stationary GPD model estimates fixed return levels:

| Model | 1-year | 5-year | 10-year |
|---|---:|---:|---:|
| Stationary GPD, all exceedances | 4.46% | 7.37% | 9.09% |

The stationary 10-year return level is about **9.09%**.

The 75th percentile high-volatility 10-year conditional return level is about **6.60%**, which may look surprising at first. The reason is that the 75th percentile of volatility is high, but not crisis-level high.

Once 90th and 95th percentile volatility states are included, the crisis-state return levels rise closer to or above the stationary estimate:

| State | 10-year return level |
|---|---:|
| High volatility, 75th percentile | 6.60% |
| Crisis volatility, 90th percentile | 8.78% |
| Extreme crisis volatility, 95th percentile | 10.34% |
| Stationary GPD | 9.09% |

This clarifies the interpretation:

> The stationary model averages across calm and crisis periods, while the non-stationary model gives conditional risk estimates depending on volatility.

## Robustness: Combining Chapter 5 Declustering with Chapter 6 Non-Stationarity

Chapter 5 showed that threshold exceedances cluster over time.

To check whether the Chapter 6 volatility effect is only driven by repeated exceedances inside crisis periods, this chapter adds a robustness model using Chapter 5-style runs declustering.

Runs declustering is applied with run length 5. Nearby exceedances are grouped into clusters, and only the maximum loss from each cluster is retained.

The declustering summary is:

| Quantity | Value |
|---|---:|
| Run length | 5 |
| Raw exceedances | 210 |
| Clusters | 129 |
| Extremal index | 0.6143 |
| Mean cluster size | 1.6279 |
| Max cluster size | 14 |
| Cluster rate per year | 3.8774 |

The extremal index is clearly below 1, confirming that extreme losses cluster even in the Chapter 6 sample.

The declustered robustness model fits a non-stationary GPD to the cluster maxima.

Let C<sub>j</sub> be the maximum loss in cluster j. The cluster excess is:

$$
Y_j = C_j - u
$$

The volatility-dependent cluster-maxima model is:

$$
Y_j \sim GPD(\sigma_j, \xi)
$$

with:

$$
\log(\sigma_j) = \beta_0 + \beta_1 z_j
$$

where z<sub>j</sub> is lagged volatility on the date of the cluster maximum.

## Declustered Cluster-Maxima Model Results

The comparison between stationary and volatility-dependent GPD models for declustered cluster maxima is:

| Model | Parameters | Observations | AIC | BIC | xi |
|---|---:|---:|---:|---:|---:|
| Stationary GPD, cluster maxima | 2 | 129 | -928.04 | -922.32 | 0.3478 |
| Volatility-dependent GPD, cluster maxima | 3 | 129 | -942.53 | -933.95 | 0.1630 |

The volatility-dependent cluster-maxima model improves both AIC and BIC.

The likelihood ratio test gives:

$$
LR = 16.49
$$

with p-value:

$$
p \approx 4.90 \times 10^{-5}
$$

The estimated volatility coefficient is:

$$
\hat{\beta}_{vol} = 0.4041
$$

with standard error:

$$
0.0900
$$

This coefficient is positive and statistically meaningful.

Therefore, volatility remains important even after accounting for clustering. This suggests that the Chapter 6 volatility effect is not only caused by repeated exceedances inside crisis episodes. Volatility also affects the severity of independent extreme-loss episodes.

## Main Model vs Declustered Robustness

The main all-exceedance model and the declustered cluster-maxima robustness model are compared below.

| Model | Observations | xi | SE(xi) | beta_vol | SE(beta_vol) | AIC | BIC |
|---|---:|---:|---:|---:|---:|---:|---:|
| All exceedances stationary GPD | 210 | 0.2840 | 0.0942 | NA | NA | -1462.23 | -1455.54 |
| All exceedances volatility GPD | 210 | 0.1208 | 0.0849 | 0.3771 | 0.0628 | -1492.76 | -1482.72 |
| Cluster maxima stationary GPD | 129 | 0.3478 | 0.1187 | NA | NA | -928.04 | -922.32 |
| Cluster maxima volatility GPD | 129 | 0.1630 | 0.1153 | 0.4041 | 0.0900 | -942.53 | -933.95 |

In both the all-exceedance model and the declustered cluster-maxima model, adding volatility improves the model.

The volatility coefficient is also similar across the two non-stationary models:

| Model | beta_vol |
|---|---:|
| All exceedances volatility GPD | 0.3771 |
| Cluster maxima volatility GPD | 0.4041 |

This supports the robustness of the volatility effect.

## Return Level Robustness

The estimated 10-year return levels from the main model and the declustered robustness model are:

| Volatility state | All exceedances | Declustered cluster maxima |
|---|---:|---:|
| Low volatility | 3.97% | 3.69% |
| Median volatility | 4.93% | 4.63% |
| High volatility | 6.60% | 6.32% |
| Crisis volatility | 8.78% | 8.66% |
| Extreme crisis volatility | 10.34% | 10.39% |

The crisis-volatility return levels are very similar after declustering.

This suggests that the main Chapter 6 conditional return level results are not driven only by clustered repeated exceedances.

## Diagnostics

The script produces diagnostic plots based on transformed residuals.

For a fitted model, the fitted GPD CDF value is:

$$
H_i(y_i)
$$

If the model is appropriate, these fitted probabilities should behave approximately like Uniform(0, 1) random variables.

Equivalently:

$$
-\log(1 - H_i(y_i))
$$

should behave approximately like Exponential(1) random variables.

Diagnostics are saved for:

- stationary GPD
- volatility-dependent GPD
- volatility + time GPD
- stationary GPD on declustered cluster maxima
- volatility-dependent GPD on declustered cluster maxima

## Comparison with Previous Chapters

The project progression so far is:

| Chapter | Method | Main model | Main idea |
|---|---|---|---|
| Chapter 3 | Block maxima | GEV | Model monthly maximum daily losses |
| Chapter 4 | Threshold exceedances | GPD | Model all losses above a high threshold |
| Chapter 5 | Dependent extremes | Declustered GPD | Account for clustering of extremes |
| Chapter 6 | Non-stationary extremes | Covariate-dependent GPD | Allow tail behavior to change with volatility |

The estimated shape parameters are:

| Chapter | Model | xi | SE(xi) |
|---|---|---:|---:|
| Chapter 3 | GEV monthly maxima | 0.1999 | 0.0433 |
| Chapter 4 | Stationary GPD | 0.2840 | 0.0942 |
| Chapter 5 | Declustered GPD | 0.3478 | 0.1176 |
| Chapter 6 | Volatility-dependent GPD, all exceedances | 0.1208 | 0.0849 |
| Chapter 6 | Volatility-dependent GPD, declustered cluster maxima | 0.1630 | 0.1153 |

The Chapter 6 shape estimates are smaller than the stationary estimates from Chapters 4 and 5.

This suggests that part of the apparent heavy-tail behavior in stationary models may be explained by changing volatility conditions. Once volatility is included in the model, less tail variation has to be absorbed by the shape parameter.

However, the Chapter 6 shape estimates remain positive, so the heavy-tail interpretation is not eliminated.

## Main Findings

1. The 97.5% loss threshold corresponds to a daily loss of about 2.42%.

2. The volatility-dependent logistic model strongly improves exceedance probability modeling.

3. A one-standard-deviation increase in lagged log-volatility multiplies the odds of a threshold exceedance by about 2.74.

4. The volatility-dependent GPD scale model strongly improves over the stationary GPD model.

5. A one-standard-deviation increase in lagged log-volatility increases the fitted GPD scale parameter by about 46%.

6. The additional time trend is not strongly supported after accounting for volatility.

7. Conditional return levels increase sharply from low-volatility to crisis-volatility states.

8. The estimated 10-year return level rises from about 3.97% in low volatility to about 10.34% in extreme crisis volatility.

9. Runs declustering with run length 5 gives an extremal index of about 0.6143, confirming clustering in extreme losses.

10. After declustering, the volatility-dependent GPD still improves over the stationary GPD.

11. The volatility coefficient remains positive after declustering, supporting the robustness of the volatility effect.

12. Crisis-volatility return levels remain similar after declustering.

13. The shape parameter remains positive across all main models, supporting the heavy-tail interpretation.

## Main Conclusion

Chapter 6 shows that extreme SPY loss behavior is non-stationary.

Volatility affects both:

1. the probability of crossing the extreme-loss threshold
2. the severity of losses once the threshold is crossed

The robustness extension shows that this volatility effect remains even after applying Chapter 5-style declustering.

The main conclusion is:

> Extreme SPY losses are heavy-tailed, clustered, and strongly dependent on market volatility. Even after accounting for clustering, volatility remains an important driver of independent extreme-loss severity.

This extends the project progression:

```text
Chapter 3: Extreme losses are heavy-tailed under a GEV block maxima model.
Chapter 4: The heavy-tail result persists under a GPD threshold model.
Chapter 5: Extreme losses cluster over time.
Chapter 6: Extreme loss risk changes with volatility, even after declustering.
```

## Main Outputs

The script saves results in:

```text
04_chapter6_nonstationary_sequences/results/
```

Important result files:

```text
threshold_summary.csv
analysis_dataset_with_covariates.csv
logistic_exceedance_model_comparison.csv
logistic_volatility_coefficients.csv
gpd_nonstationary_model_comparison.csv
likelihood_ratio_tests.csv
gpd_nonstationary_parameter_estimates.csv
volatility_scenarios.csv
conditional_return_levels_by_volatility.csv
stationary_return_levels.csv
chapter6_declustering_summary.csv
chapter6_clusters_run_length_5.csv
cluster_maxima_nonstationary_model_comparison.csv
cluster_maxima_likelihood_ratio_test.csv
cluster_maxima_parameter_estimates.csv
cluster_volatility_scenarios.csv
cluster_conditional_return_levels_by_volatility.csv
combined_conditional_return_levels_all_vs_cluster.csv
main_vs_declustered_nonstationary_summary.csv
chapter3_to_chapter6_shape_comparison.csv
chapter6_model_summary.txt
session_info.txt
```

Important figures are saved in:

```text
04_chapter6_nonstationary_sequences/figures/
```

Important figure files:

```text
daily_losses_with_threshold_and_exceedances.png
volatility_with_exceedance_days.png
exceedance_probability_vs_volatility.png
gpd_scale_vs_volatility.png
fitted_scale_at_exceedances_over_time.png
diagnostics_stationary_gpd.png
diagnostics_volatility_gpd.png
diagnostics_volatility_time_gpd.png
conditional_return_levels_by_volatility.png
stationary_vs_conditional_return_levels.png
chapter6_declustered_cluster_maxima.png
diagnostics_cluster_stationary_gpd.png
diagnostics_cluster_volatility_gpd.png
ten_year_return_levels_all_vs_cluster.png
main_vs_declustered_shape_comparison.png
main_vs_declustered_beta_vol_comparison.png
chapter3_to_chapter6_shape_comparison.png
```

## Software

The analysis is written in R.

Main packages:

```r
quantmod
xts
zoo
ismev
```
