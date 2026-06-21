# Chapter 4: Threshold Models for SPY Daily Losses

This chapter applies the **threshold exceedance** approach from Chapter 4 of Coles' *An Introduction to Statistical Modeling of Extreme Values* to daily SPY losses.

The goal is to move beyond the block maxima method used in Chapter 3 and model all sufficiently large daily losses above a high threshold.

Chapter 3 used one maximum per month. Chapter 4 uses all losses above a threshold. In other words, Chapter 4 stops throwing away useful tail observations like a statistical aristocrat.

## Objective

The objective of this chapter is to model extreme daily SPY losses using the **Generalized Pareto Distribution**, or GPD.

The main questions are:

- What threshold should be used to define an extreme loss?
- Are threshold exceedances well described by the GPD?
- Is the estimated tail shape parameter positive?
- Are return level estimates stable across threshold choices?
- Are the Chapter 4 threshold results consistent with the Chapter 3 block maxima results?

## Data

The analysis uses daily adjusted closing prices for SPY.

Let P<sub>t</sub> denote the adjusted closing price on day t. Daily log returns are computed as:

$$
R_t = \log(P_t) - \log(P_{t-1})
$$

Because the focus is on market losses, returns are transformed into losses:

$$
L_t = -R_t
$$

This transformation turns large negative returns into large positive losses. EVT is then applied to the upper tail of the loss distribution.

## Methodology

### Threshold exceedances

Instead of taking one maximum per block, the threshold method keeps all observations above a high threshold.

Let u be a high threshold. An exceedance occurs when:

$$
L_t > u
$$

The excess over the threshold is:

$$
Y = L_t - u \mid L_t > u
$$

So if a daily loss exceeds the threshold, the model studies how far above the threshold it is.

### Generalized Pareto Distribution

For a sufficiently high threshold u, EVT suggests that threshold excesses can be approximated by the **Generalized Pareto Distribution**.

The GPD distribution function is:

$$
H(y) =
1 -
\left(
1 + \xi \frac{y}{\sigma}
\right)^{-1/\xi}
$$

with support condition:

$$
1 + \xi \frac{y}{\sigma} > 0
$$

The parameters are:

- $$\sigma$$ > 0: scale parameter
- $$\xi$$: shape parameter

The shape parameter xi controls the tail behavior:

- $$\xi$$ > 0: heavy-tailed distribution
- $$\xi$$ = 0: exponential-type tail
- $$\xi$$ < 0: bounded upper tail

For financial losses, xi is the key parameter. A positive xi suggests that extreme losses are heavy-tailed.

## Threshold Selection

Three empirical thresholds are considered:

| Threshold | Loss level | Exceedances |
|---|---:|---:|
| 95% | 1.82% | 421 |
| 97.5% | 2.42% | 211 |
| 99% | 3.25% | 85 |

The 95% threshold gives more observations but may include less extreme losses.

The 99% threshold focuses on the most extreme losses but leaves only 85 exceedances, making estimates more uncertain.

The 97.5% threshold is used as the main specification because it balances tail relevance and sample size.

## Parameter Stability

To study threshold sensitivity, GPD models are fitted across a grid of thresholds from the 90th to the 99th percentile.

The main diagnostic is the stability of the estimated shape parameter xi.

If the GPD model is appropriate above a threshold, the estimate of xi should remain reasonably stable as the threshold increases.

The results show that xi remains positive across the threshold range, although uncertainty increases at very high thresholds due to fewer exceedances.

## GPD Estimates

The fitted GPD models give the following estimates:

| Threshold | Exceedances | sigma | xi | SE(xi) |
|---|---:|---:|---:|---:|
| 95% | 421 | 0.0079 | 0.2178 | 0.0578 |
| 97.5% | 211 | 0.0083 | 0.2927 | 0.0943 |
| 99% | 85 | 0.0129 | 0.1595 | 0.1417 |

All three shape estimates are positive.

This supports the conclusion that extreme SPY losses are heavy-tailed.

The 99% estimate has a much larger standard error because there are fewer exceedances. This is the usual bias-variance tradeoff, except now it is wearing a finance costume.

## Return Levels

For daily data, return periods are measured in trading days.

The analysis uses:

$$
m = 252
$$

for a 1-year return level,

$$
m = 1260
$$

for a 5-year return level, and

$$
m = 2520
$$

for a 10-year return level.

Let zeta<sub>u</sub> be the probability of exceeding the threshold u:

$$
\zeta_u = P(L > u)
$$

For xi not equal to zero, the GPD return level is:

$$
x_m =
u +
\frac{\sigma}{\xi}
\left[
(m \zeta_u)^\xi - 1
\right]
$$

The estimated return levels are:

| Threshold | 1-year | 5-year | 10-year |
|---|---:|---:|---:|
| 95% | 4.51% | 7.16% | 8.63% |
| 97.5% | 4.45% | 7.38% | 9.13% |
| 99% | 4.56% | 7.31% | 8.73% |

The return levels are fairly stable across thresholds. This is a useful sign: the estimated extreme loss levels are not completely driven by one arbitrary threshold choice.

Using the 97.5% threshold as the main specification, the model estimates:

- 1-year return level: about 4.45%
- 5-year return level: about 7.38%
- 10-year return level: about 9.13%

These are daily loss levels, not multi-day or monthly losses.

## Diagnostics

The GPD fit is evaluated using:

- probability plots
- quantile plots
- return level plots
- density plots
- mean residual life plots
- parameter stability plots

The probability plots generally show a reasonable fit.

The quantile plots show some deviation in the most extreme observations, especially in the far upper tail. This is common in financial losses, where the largest crash days are difficult to model precisely.

The return level plots show increasing uncertainty for longer return periods, which is expected because rare-event extrapolation is hard. Apparently probability refuses to give precise answers about events we barely observe. Rude, but fair.

## Comparison with Chapter 3

Chapter 3 used the **block maxima** method.

The daily losses were divided into monthly blocks, and only the largest daily loss in each month was retained:

$$
M_j = \max(L_t \text{ in month } j)
$$

The monthly maxima were modeled using the **Generalized Extreme Value** distribution.

Chapter 4 instead uses the **threshold exceedance** method. It keeps all daily losses above a high threshold and models the excesses using the GPD.

The key comparison is:

| Chapter | Method | Model | Main object |
|---|---|---|---|
| Chapter 3 | Block maxima | GEV | Monthly maximum daily losses |
| Chapter 4 | Threshold exceedances | GPD | Daily losses above a high threshold |

### Shape parameter comparison

| Model | xi estimate |
|---|---:|
| Chapter 3 GEV monthly maxima | about 0.20 |
| Chapter 4 GPD 95% threshold | 0.2178 |
| Chapter 4 GPD 97.5% threshold | 0.2927 |
| Chapter 4 GPD 99% threshold | 0.1595 |

Both methods produce positive estimates of xi.

This means that both the block maxima approach and the threshold exceedance approach suggest heavy-tailed behavior in extreme SPY losses.

### Return level comparison

Using the Chapter 3 GEV model, the estimated return levels were approximately:

| Method | 1-year | 5-year | 10-year |
|---|---:|---:|---:|
| Chapter 3 GEV | 3.82% | 6.31% | 7.64% |
| Chapter 4 GPD, 97.5% threshold | 4.45% | 7.38% | 9.13% |

The Chapter 4 threshold model gives somewhat higher return levels than the Chapter 3 block maxima model.

This makes sense because the threshold model uses more tail information, while the block maxima model keeps only one observation per month.

The main conclusion is consistent across both chapters:

> Extreme SPY daily losses appear to be heavy-tailed.

## Main Conclusion

The threshold exceedance analysis supports and strengthens the Chapter 3 conclusion.

The GPD estimates are positive across the 95%, 97.5%, and 99% thresholds. Return levels are also reasonably stable across thresholds.

The 97.5% threshold is used as the main specification because it provides a good balance between focusing on extreme losses and retaining enough exceedances for estimation.

Overall, Chapter 4 shows that the heavy-tail conclusion is not just an artifact of the block maxima method. It also appears under a threshold-based EVT model.

## Files

This folder contains:

```text
02_chapter4_threshold_models/
├── README.md
├── chapter4_threshold_models.R
├── figures/
└── results/
```

The `figures/` folder contains diagnostic plots, threshold plots, stability plots, and return level plots.

The `results/` folder contains CSV tables with fitted parameters, return levels, threshold summaries, and comparison tables.

## Key Outputs

Important result files:

```text
results/threshold_summary.csv
results/gpd_estimates_by_threshold.csv
results/gpd_return_levels.csv
results/parameter_stability.csv
results/chapter3_chapter4_shape_comparison.csv
results/chapter3_chapter4_return_level_comparison.csv
results/chapter4_model_summary.txt
```

Important figures:

```text
figures/daily_losses_thresholds_tail_zoom.png
figures/mean_residual_life_high_thresholds.png
figures/shape_parameter_stability_with_ci.png
figures/gpd_diagnostics_975.png
figures/gpd_return_levels_comparison.png
figures/chapter3_vs_chapter4_return_levels.png
figures/chapter3_vs_chapter4_shape_comparison.png
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
