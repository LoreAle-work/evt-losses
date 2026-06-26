# Chapter 7: Point Process Models for Extremes

## Overview

This folder applies the **point process approach to extreme value theory** to daily SPY losses.

The purpose of Chapter 7 is to connect the two main EVT approaches used earlier in the project:

1. **Block maxima models**, studied in Chapter 3 through the Generalized Extreme Value distribution.
2. **Threshold exceedance models**, studied in Chapter 4 through the Generalized Pareto distribution.

The point process framework provides a unified view of extremes by modeling exceedances as points in time-value space:

$$(t, X_t)$$

where $$t$$ is time and $$X_t$$ is the daily loss.

In this chapter, the analysis is applied to daily SPY log losses:

$$
L_t = -R_t
$$

where:

$$
R_t = \log(P_t) - \log(P_{t-1})
$$

and $$P_t$$ is the adjusted closing price of SPY.

Large positive values of $$L_t$$ correspond to large negative market returns.

---

## Main Goal

The main goal of this chapter is to estimate and compare point process models for extreme SPY losses.

The chapter has two parts:

1. A **stationary point process model**, where the parameters are constant over time.
2. Several **volatility-dependent point process models**, where lagged volatility affects the distribution of extremes.

The stationary point process model is used to show the theoretical connection between:

- Chapter 3 GEV block maxima,
- Chapter 4 GPD threshold exceedances,
- Chapter 7 point process extremes.

The volatility-dependent models extend the Chapter 6 non-stationary EVT analysis by allowing the point process parameters to change with market volatility.

---

## Data

The data consist of daily adjusted closing prices for SPY from Yahoo Finance.

The analysis uses:

- Daily log returns.
- Daily log losses.
- A high loss threshold.
- A lagged volatility covariate.

The final usable dataset contains:

| Quantity | Value |
|---|---:|
| Usable observations | 8,386 |
| Approximate years | 33.28 |
| Threshold probability | 97.5% |
| Threshold value | 0.02419849 |
| Threshold percent | 2.419849% |
| Number of exceedances | 210 |
| Exceedance fraction | 0.02504174 |

The threshold is the empirical 97.5% quantile of daily losses. Any loss above this level is treated as an extreme loss.

---

## Loss Definition

Daily log returns are computed as:

$$
R_t = \log(P_t) - \log(P_{t-1})
$$

Daily losses are then defined as:

$$
L_t = -R_t
$$

This transforms large negative returns into large positive losses.

For example, if SPY falls by about 5% in one day, the corresponding loss is approximately:

$$
L_t \approx 0.05
$$

This makes the upper tail of the loss distribution the object of interest.

---

## Volatility Covariate

To study non-stationarity, the analysis uses lagged 21-day realized volatility.

First, the rolling standard deviation of daily log returns is computed over the previous 21 trading days:

$$
RV_t = sd(R_{t-20}, \dots, R_t)
$$

Then the volatility measure is lagged by one day to avoid look-ahead bias:

$$
RV_{t-1}
$$

The covariate used in the models is the standardized log-volatility:

$$
z_t = \frac{\log(RV_{t-1}) - \overline{\log(RV)}}{sd(\log(RV))}
$$

This means that:

- $$z_t = 0$$ corresponds to average log-volatility.
- $$z_t = 1$$ corresponds to volatility one standard deviation above average.
- $$z_t = -1$$ corresponds to volatility one standard deviation below average.

The volatility covariate captures market regimes such as calm periods, high-volatility periods, and crisis periods.

---

## Point Process Theory

The point process model represents extremes as points above a high threshold.

For a high threshold $$u$$, the model focuses on observations satisfying:

$$
L_t > u
$$

The point process formulation is closely connected to the GEV distribution.

The annual tail measure is:

$$\Lambda(x)=
\left[
1 + \xi \left(\frac{x - \mu}{\sigma}\right)
\right]^{-1/\xi}
$$

where:

- $$\mu$$ is the location parameter,
- $$\sigma > 0$$ is the scale parameter,
- $$\xi$$ is the shape parameter.

The support condition is:

$$
1 + \xi \left(\frac{x - \mu}{\sigma}\right) > 0
$$

The shape parameter $$\xi$$ controls the tail behavior:

| Value of $$\xi$$ | Tail Type | Interpretation |
|---|---|---|
| $$\xi > 0$$ | Fréchet | Heavy-tailed |
| $$\xi = 0$$ | Gumbel | Exponential-type tail |
| $$\xi < 0$$ | Weibull | Bounded upper tail |

In financial loss modeling, a positive $$\xi$$ is usually interpreted as evidence of heavy-tailed extreme losses.

---

## Stationary Point Process Model

The stationary point process model assumes that the parameters are constant over time:

$$
\mu_t = \mu
$$

$$
\sigma_t = \sigma
$$

$$
\xi_t = \xi
$$

The model estimates:

$$
(\mu, \sigma, \xi)
$$

using exceedances above the 97.5% loss threshold.

### Stationary Point Process Estimates

| Parameter | Estimate | Standard Error |
|---|---:|---:|
| $$\mu$$ | 0.04461811 | 0.002122752 |
| $$\log(\sigma)$$ | -4.25149406 | 0.129167560 |
| $$\sigma$$ | 0.01424294 | |
| $$\xi$$ | 0.28411083 | 0.093000916 |

The estimated shape parameter is positive:

$$
\hat{\xi} = 0.2841
$$

This supports a heavy-tailed interpretation of SPY extreme losses.

---

## Stationary Point Process Return Levels

For the point process model, a $$T$$-year return level $$z_T$$ is defined by:

$$
\Lambda(z_T) = \frac{1}{T}
$$

Solving for $$z_T$$ gives:

$$
z_T =
\mu + \frac{\sigma}{\xi}
\left(
T^\xi - 1
\right)
$$

for $$\xi \neq 0$$.

For $$\xi = 0$$, the return level is:

$$
z_T = \mu + \sigma \log(T)
$$

### Estimated Stationary Point Process Return Levels

| Return Period | Return Level |
|---:|---:|
| 1 year | 4.461811% |
| 5 years | 7.368130% |
| 10 years | 9.091881% |
| 20 years | 11.190822% |
| 50 years | 14.682425% |
| 100 years | 17.998202% |

These return levels represent daily losses expected to be exceeded approximately once every $$T$$ years under the stationary point process model.

---

## Comparison with Chapter 3 and Chapter 4

One of the main purposes of Chapter 7 is to show that the point process model connects the GEV and GPD approaches.

The comparison uses:

- Chapter 3 monthly GEV block maxima,
- Chapter 4 stationary GPD threshold exceedances,
- Chapter 7 stationary point process exceedances.

Because Chapter 3 used monthly block maxima, the GEV return levels are adjusted using 12 monthly blocks per year.

For a $$T$$-year return level in the monthly GEV model, the relevant number of blocks is:

$$
m = 12T
$$

The GEV return level is computed using:

$$
p = 1 - \frac{1}{12T}
$$

### Return Level Comparison

| Model | 1-year | 5-year | 10-year |
|---|---:|---:|---:|
| Chapter 3 GEV monthly maxima | 3.821821% | 6.325256% | 7.664937% |
| Chapter 4 stationary GPD | 4.462900% | 7.365751% | 9.086000% |
| Chapter 7 stationary point process | 4.461811% | 7.368130% | 9.091881% |

The Chapter 4 GPD and Chapter 7 stationary point process return levels are almost identical.

This is the central Chapter 7 result:

> The stationary point process model reproduces the stationary GPD threshold model almost exactly, confirming the theoretical connection between threshold exceedances and point process extremes.

The Chapter 3 GEV estimates are lower, especially at longer horizons, because monthly block maxima use less information from the tail than threshold exceedance methods.

---

## Volatility-Dependent Point Process Models

The stationary model assumes that the extreme-loss process does not change over time.

However, earlier chapters showed that SPY extremes are not fully stationary:

- Chapter 5 showed that extremes cluster.
- Chapter 6 showed that volatility strongly affects exceedance probability and severity.

Therefore, this chapter also estimates volatility-dependent point process models.

The covariate is standardized lagged log-volatility:

$$
z_t
$$

Three non-stationary point process models are fitted:

1. Volatility-dependent scale model.
2. Volatility-dependent location model.
3. Volatility-dependent location and scale model.

---

## Model 1: Volatility-Dependent Scale Point Process

The scale-only model is:

$$
\mu_t = \mu
$$

$$
\log(\sigma_t) = \beta_0 + \beta_1 z_t
$$

$$
\xi_t = \xi
$$

This model allows volatility to affect the scale of the extreme-loss process, while keeping location fixed.

### Scale Model Estimates

| Parameter | Estimate | Standard Error |
|---|---:|---:|
| $$\mu$$ | 0.05376595 | 0.002821488 |
| $$\beta_0$$ | -3.93036198 | 0.094376178 |
| $$\beta_1$$ | -0.19367046 | 0.023850725 |
| $$\xi$$ | 0.16533822 | 0.062353671 |

The volatility coefficient is negative:

$$
\hat{\beta}_1 = -0.1937
$$

This means that, in the scale-only specification, fitted scale decreases as volatility increases.

This produced a counterintuitive result: higher-volatility states had lower conditional return levels at long horizons.

The explanation is that the model was probably too restrictive. It forced volatility to act only through scale, even though volatility likely shifts the whole extreme-loss process upward.

---

## Model 2: Volatility-Dependent Location Point Process

The location-only model is:

$$
\mu_t = \mu_0 + \mu_1 z_t
$$

$$
\sigma_t = \sigma
$$

$$
\xi_t = \xi
$$

This model allows volatility to shift the location of the extreme-loss process.

### Location Model Estimates

| Parameter | Estimate | Standard Error |
|---|---:|---:|
| $$\mu_0$$ | 0.03801090 | 0.001387859 |
| $$\mu_1$$ | 0.01341247 | 0.001366779 |
| $$\log(\sigma)$$ | -4.60413021 | 0.064662304 |
| $$\sigma$$ | 0.01001041 | |
| $$\xi$$ | -0.12785000 | 0.037226436 |

The volatility-location coefficient is positive:

$$
\hat{\mu}_1 = 0.0134
$$

This means that higher volatility shifts the fitted point process location upward.

The fitted location increases from about 2.84% in low-volatility states to about 6.02% in extreme-crisis-volatility states.

---

## Model 3: Volatility-Dependent Location and Scale Point Process

The most flexible model allows both location and scale to depend on volatility:

$$
\mu_t = \mu_0 + \mu_1 z_t
$$

$$
\log(\sigma_t) = \beta_0 + \beta_1 z_t
$$

$$
\xi_t = \xi
$$

This model allows volatility to affect both:

1. The level of the extreme-loss process.
2. The spread of the extreme-loss process.

### Location + Scale Model Estimates

| Parameter | Estimate | Standard Error |
|---|---:|---:|
| $$\mu_0$$ | 0.03591829 | 0.001393618 |
| $$\mu_1$$ | 0.01402089 | 0.001381389 |
| $$\beta_0$$ | -4.67680199 | 0.095866385 |
| $$\beta_1$$ | 0.16552245 | 0.036310732 |
| $$\xi$$ | 0.04451533 | 0.059660978 |

Both volatility effects are positive:

$$
\hat{\mu}_1 = 0.0140
$$

and

$$
\hat{\beta}_1 = 0.1655
$$

This means that higher volatility both shifts the point process upward and increases its scale.

This model gives the most intuitive non-stationary interpretation.

---

## Model Comparison

The four point process models are compared using AIC, BIC, and likelihood-ratio tests.

### AIC and BIC Comparison

| Model | Parameters | NLL | AIC | BIC | $$\xi$$ | Convergence |
|---|---:|---:|---:|---:|---:|---:|
| Stationary point process | 3 | -909.9261 | -1813.852 | -1803.811 | 0.2841 | 0 |
| Volatility scale point process | 4 | -969.4045 | -1930.809 | -1917.421 | 0.1653 | 0 |
| Volatility location point process | 4 | -1033.7484 | -2059.497 | -2046.108 | -0.1279 | 0 |
| Volatility location + scale point process | 5 | -1043.3230 | -2076.646 | -2059.910 | 0.0445 | 0 |

The best model by both AIC and BIC is:

> **Volatility-dependent location + scale point process**

This model has the lowest AIC and the lowest BIC.

All models converged successfully.

---

## Likelihood-Ratio Tests

The likelihood-ratio tests compare nested models.

| Comparison | LR Statistic | df | p-value |
|---|---:|---:|---:|
| Stationary vs volatility scale | 118.95676 | 1 | $$1.070377 \times 10^{-27}$$ |
| Stationary vs volatility location | 247.64467 | 1 | $$8.470969 \times 10^{-56}$$ |
| Volatility scale vs volatility location + scale | 147.83698 | 1 | $$5.149524 \times 10^{-34}$$ |
| Volatility location vs volatility location + scale | 19.14907 | 1 | $$1.208965 \times 10^{-5}$$ |

The tests show that:

1. The volatility scale model improves substantially over the stationary model.
2. The volatility location model improves even more over the stationary model.
3. The location + scale model improves over the scale-only model.
4. The location + scale model also improves over the location-only model.

Therefore, the evidence supports allowing volatility to affect both location and scale.

---

## Volatility Scenarios

The non-stationary return levels are evaluated at five volatility states:

| Volatility State | Quantile | $$z_t$$ |
|---|---:|---:|
| Low volatility | 25% | -0.71810142 |
| Median volatility | 50% | -0.07165893 |
| High volatility | 75% | 0.65746431 |
| Crisis volatility | 90% | 1.30156991 |
| Extreme crisis volatility | 95% | 1.65204628 |

These scenarios help interpret how extreme-loss risk changes across market regimes.

---

## Location + Scale Scenarios

For the best model, both fitted location and scale increase with volatility.

| Volatility State | $$\mu(z)$$ | $$\sigma(z)$$ | $$\xi$$ |
|---|---:|---:|---:|
| Low volatility | 0.02584987 | 0.008265510 | 0.04451533 |
| Median volatility | 0.03491357 | 0.009198976 | 0.04451533 |
| High volatility | 0.04513653 | 0.010378936 | 0.04451533 |
| Crisis volatility | 0.05416746 | 0.011546615 | 0.04451533 |
| Extreme crisis volatility | 0.05908145 | 0.012236264 | 0.04451533 |

This is the key improvement over the scale-only model.

The scale-only model forced volatility to act only through $$\sigma_t$$, which produced a negative scale coefficient.

Once $$\mu_t$$ is allowed to vary, the fitted volatility effects become more interpretable:

- Higher volatility increases location.
- Higher volatility increases scale.
- Conditional return levels increase with volatility.

---

## Conditional Return Levels: Location + Scale Model

For the best model, the conditional return levels are:

| Volatility State | 1-year | 5-year | 10-year | 20-year | 50-year | 100-year |
|---|---:|---:|---:|---:|---:|---:|
| Low volatility | 2.584987% | 3.964082% | 4.589151% | 5.233806% | 6.117102% | 6.809639% |
| Median volatility | 3.491357% | 5.026200% | 5.721861% | 6.439321% | 7.422371% | 8.193120% |
| High volatility | 4.513653% | 6.245371% | 7.030265% | 7.839754% | 8.948901% | 9.818514% |
| Crisis volatility | 5.416746% | 7.343291% | 8.216489% | 9.117049% | 10.350981% | 11.318429% |
| Extreme crisis volatility | 5.908145% | 7.949758% | 8.875109% | 9.829458% | 11.137089% | 12.162320% |

The return levels increase monotonically with volatility.

For example, the estimated 10-year return level rises from:

$$
4.59\%
$$

in low-volatility states to:

$$
8.88\%
$$

in extreme-crisis-volatility states.

This result is consistent with the financial intuition that extreme losses become more severe during volatile market regimes.

---

## Shape Parameter Comparison

The estimated shape parameters across chapters and models are:

| Model | $$\xi$$ | Standard Error |
|---|---:|---:|
| Chapter 3 GEV monthly maxima | 0.20179242 | 0.04343350 |
| Chapter 4 stationary GPD | 0.28334863 | 0.09342161 |
| Chapter 7 stationary point process | 0.28411083 | 0.09300092 |
| Chapter 7 volatility scale point process | 0.16533822 | 0.06235367 |
| Chapter 7 volatility location point process | -0.12785000 | 0.03722644 |
| Chapter 7 volatility location + scale point process | 0.04451533 | 0.05966098 |

The stationary GPD and stationary point process estimates are almost identical:

$$
\hat{\xi}_{GPD} = 0.2833
$$

$$
\hat{\xi}_{PP} = 0.2841
$$

This confirms the theoretical link between the GPD threshold model and the point process model.

The shape parameter falls after including volatility. This suggests that part of the heavy-tailed behavior estimated by stationary models was actually due to unmodeled volatility dependence.

In other words:

> Some of what looks like heavy-tailedness in a stationary model is explained by volatility regimes once non-stationarity is included.

The best non-stationary model still has a positive shape estimate:

$$
\hat{\xi} = 0.0445
$$

but it is much smaller than the stationary estimate.

---

## Diagnostics

The stationary point process diagnostics include:

1. Probability plot.
2. Exponential QQ plot.
3. Exponential residual density.
4. Residuals by exceedance order.

The probability plot follows the 45-degree line closely, indicating that the fitted CDF values are approximately uniform.

The exponential QQ plot is also close to the reference line for most of the distribution, with mild deviations in the far upper tail.

The residual density broadly follows the exponential reference density.

The residual-by-order plot shows some clustering of large residuals around crisis periods, which is consistent with the dependence found in Chapter 5.

Overall, the stationary point process model provides an acceptable baseline fit.

---

## Interpretation of the Scale-Only Issue

The scale-only model originally produced decreasing return levels as volatility increased.

This happened because the model was restricted to:

$$
\mu_t = \mu
$$

$$
\log(\sigma_t) = \beta_0 + \beta_1 z_t
$$

with $$\mu$$ fixed.

The fitted scale coefficient was negative:

$$
\hat{\beta}_1 = -0.1937
$$

This caused return levels to decrease with volatility.

This should not be interpreted as high volatility reducing financial risk.

Instead, it shows that the scale-only point process specification was too restrictive. Volatility likely shifts the entire extreme-loss process upward, which is better captured by allowing $$\mu_t$$ to vary.

Once location is allowed to move, the model finds:

$$
\hat{\mu}_1 > 0
$$

and, in the location + scale model:

$$
\hat{\beta}_1 > 0
$$

Thus, the decreasing-return-level issue disappears.

---

## Main Findings

The main findings of Chapter 7 are:

1. **The stationary point process model reproduces the Chapter 4 GPD results.**

   The stationary point process and stationary GPD return levels are almost identical. This confirms the theoretical connection between point process extremes and threshold exceedance models.

2. **SPY extreme losses are heavy-tailed under stationary EVT models.**

   The stationary point process estimate is:

   $$\hat{\xi} = 0.2841$$

   which supports a heavy-tailed, Fréchet-type interpretation.

3. **Volatility strongly improves point process modeling.**

   All volatility-dependent models improve substantially over the stationary model.

4. **The best model allows both location and scale to depend on volatility.**

   The volatility-dependent location + scale point process has the lowest AIC and BIC.

5. **Higher volatility shifts the extreme-loss process upward.**

   The fitted location parameter increases from low-volatility states to crisis-volatility states.

6. **Higher volatility also increases scale in the best model.**

   In the location + scale model, fitted scale increases with volatility.

7. **Conditional return levels increase with volatility.**

   The 10-year return level rises from 4.59% in low-volatility states to 8.88% in extreme-crisis-volatility states.

8. **Non-stationarity reduces residual tail heaviness.**

   The estimated $$\xi$$ is much smaller in the volatility-dependent location + scale model than in the stationary model. This suggests that part of the apparent heavy-tailedness in stationary EVT models is explained by volatility regimes.

---

## Final Conclusion

The Chapter 7 point process analysis confirms the theoretical bridge between block maxima and threshold exceedance methods.

The stationary point process model gives return levels almost identical to the Chapter 4 stationary GPD model, showing that both approaches describe the same underlying extreme-value structure.

The non-stationary point process analysis extends this result by showing that volatility is a major driver of extreme-loss behavior. The best model allows volatility to affect both the location and scale of the point process. In this model, higher volatility shifts the extreme-loss process upward, increases scale, and produces higher conditional return levels.

Overall, Chapter 7 strengthens the main conclusion of the project:

> SPY extreme losses are heavy-tailed, clustered, and strongly dependent on market volatility.

---

## Output Files

The script produces result files in:

```text
05_chapter7_point_processes/results/
