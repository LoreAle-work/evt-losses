# Chapter 5: Extremes of Dependent Sequences

This chapter extends the threshold exceedance analysis from Chapter 4 by accounting for **dependence among extreme losses**.

Chapter 4 modeled all threshold exceedances as if they were approximately independent. That is a useful starting point, but financial losses often violate this assumption. Large losses tend to cluster during crisis periods such as 2008 or 2020.

This chapter studies that clustering directly.

The main goal is to understand whether extreme SPY losses behave like independent rare events or whether they arrive in clusters.

## Main Idea

The Chapter 4 threshold model keeps all daily losses above a high threshold and fits a Generalized Pareto Distribution.

Chapter 5 asks a deeper question:

> Are those exceedances independent, or do they cluster over time?

If exceedances cluster, treating every exceedance as a separate independent event overstates the amount of independent tail information in the data.

To address this, this chapter uses **runs declustering**. Nearby exceedances are grouped into clusters, and only the maximum loss from each cluster is retained.

So the progression is:

```text
Chapter 3: Monthly block maxima -> GEV
Chapter 4: Daily threshold exceedances -> GPD
Chapter 5: Declustering threshold exceedances -> GPD on cluster maxima
```

Chapter 5 does not replace Chapter 4. It checks whether the Chapter 4 independence assumption is too simplistic.

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

The dataset contains **8403 daily observations**.

## Threshold Choice

The threshold is the 97.5% empirical quantile of daily losses.

| Quantity | Value |
|---|---:|
| Threshold probability | 97.5% |
| Threshold value | 0.0242 |
| Threshold in percent | 2.42% |
| Total observations | 8403 |
| Number of exceedances | 211 |
| Exceedance fraction | 2.51% |

So an exceedance is defined as a daily SPY loss larger than approximately **2.42%**.

Mathematically, an exceedance occurs when:

$$
L_t > u
$$

where u is the threshold.

The excess above the threshold is:

$$
Y_t = L_t - u \mid L_t > u
$$

## Why Dependence Matters

In the Chapter 4 threshold model, the 211 exceedances were treated as if they were approximately independent.

But financial extremes often arrive close together. During crisis periods, one large loss is frequently followed by another large loss soon after.

This matters because if exceedances are clustered, then the effective number of independent extreme events is smaller than the raw number of exceedances.

In simple terms:

```text
211 exceedances does not necessarily mean 211 independent extreme events.
```

That is the entire reason Chapter 5 exists. Statistics looked at financial markets and said, “No, I do not trust these observations to behave themselves.”

## Runs Declustering

Runs declustering groups nearby exceedances into clusters.

A new cluster begins only after a specified number of consecutive non-exceedance days. This number is called the **run length**.

For example, with run length 5:

```text
A new cluster starts only after at least 5 consecutive trading days without an exceedance.
```

This roughly corresponds to one trading week.

The idea is:

- exceedances close together are treated as part of the same extreme episode
- each cluster represents one independent extreme event
- only the largest loss in each cluster is retained

## Extremal Index

The **extremal index**, usually denoted theta, measures the degree of clustering in extremes.

A simple runs-based estimate is:

$$
\hat{\theta}
=
\frac{\text{number of clusters}}{\text{number of exceedances}}
$$

Interpretation:

| Extremal index | Meaning |
|---:|---|
| theta close to 1 | Extremes behave approximately independently |
| theta below 1 | Extremes cluster |
| smaller theta | stronger clustering |

For run length 5:

$$
\hat{\theta}
=
\frac{130}{211}
\approx 0.616
$$

This is clearly below 1, so the results indicate clustering in extreme SPY losses.

The reciprocal gives an approximate average cluster size:

$$
\frac{1}{\hat{\theta}}
\approx
1.62
$$

So, under the run length 5 rule, each independent extreme episode contains about 1.62 exceedances on average.

## Run Length Sensitivity

Different run lengths give different cluster definitions. A longer run length groups exceedances more aggressively, producing fewer clusters and a smaller extremal index.

| Run length | Exceedances | Clusters | Mean cluster size | Max cluster size | Extremal index | Cluster rate per year |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 211 | 187 | 1.13 | 5 | 0.886 | 5.61 |
| 3 | 211 | 148 | 1.43 | 12 | 0.701 | 4.44 |
| 5 | 211 | 130 | 1.62 | 14 | 0.616 | 3.90 |
| 10 | 211 | 97 | 2.18 | 25 | 0.460 | 2.91 |

The pattern is exactly what we expect:

- as run length increases, the number of clusters decreases
- mean cluster size increases
- the estimated extremal index decreases
- the estimated number of independent extreme episodes per year decreases

The main specification uses **run length 5**, because it corresponds roughly to one trading week and gives a moderate declustering rule.

## Model Comparison

This chapter compares two GPD models:

| Model | Data used | Interpretation |
|---|---|---|
| Chapter 4 iid GPD | All 211 threshold exceedances | Treats exceedances as approximately independent |
| Chapter 5 declustered GPD | 130 cluster maxima | Accounts for clustering by using one maximum per cluster |

The Chapter 4 model is more data-rich but assumes independence.

The Chapter 5 model uses fewer observations but is more realistic for dependent financial extremes.

## GPD Model Estimates

Both models fit a Generalized Pareto Distribution to threshold exceedances.

The GPD distribution function is:

$$
H(y)
=
1 -
\left(
1 + \xi \frac{y}{\sigma}
\right)^{-1/\xi}
$$

where:

- sigma is the scale parameter
- xi is the shape parameter

The shape parameter xi determines the tail behavior:

| xi value | Interpretation |
|---:|---|
| xi > 0 | Heavy-tailed distribution |
| xi = 0 | Exponential-type tail |
| xi < 0 | Bounded upper tail |

The fitted GPD estimates are:

| Model | Observations used | sigma | xi | SE(sigma) | SE(xi) |
|---|---:|---:|---:|---:|---:|
| Chapter 4 iid exceedances | 211 | 0.00829 | 0.2932 | 0.00092 | 0.0944 |
| Chapter 5 declustered cluster maxima | 130 | 0.00686 | 0.3577 | 0.00094 | 0.1185 |

Both xi estimates are positive.

This means that both the iid threshold model and the declustered threshold model support heavy-tailed behavior in extreme SPY losses.

The Chapter 5 estimate is higher, but it also has a larger standard error because the declustered model uses only 130 cluster maxima instead of 211 raw exceedances.

## Why the Declustered Model Uses Fewer Observations

The Chapter 4 model uses every exceedance:

```text
211 exceedances
```

The Chapter 5 model groups nearby exceedances into clusters and keeps only the largest loss from each cluster:

```text
130 cluster maxima
```

This reduces the sample size but avoids counting several losses from the same crisis episode as separate independent events.

That tradeoff is the whole point:

```text
Chapter 4: more observations, stronger independence assumption
Chapter 5: fewer observations, weaker independence assumption
```

Neither is magically perfect. Because apparently reality was not designed for clean likelihood functions.

## Return Level Comparison

Return levels estimate the size of a daily loss expected to be exceeded once over a given time horizon.

The comparison between the Chapter 4 iid model and the Chapter 5 declustered model is:

| Return period | Trading days | Chapter 4 iid GPD | Chapter 5 declustered GPD |
|---|---:|---:|---:|
| 1 year | 252 | 4.45% | 3.62% |
| 5 years | 1260 | 7.38% | 6.05% |
| 10 years | 2520 | 9.14% | 7.61% |

The declustered return levels are lower than the iid return levels.

At first, this may seem surprising because the declustered model has a larger xi estimate. However, the return level depends not only on tail shape, but also on the event rate.

Chapter 4 uses the raw exceedance rate:

$$
\hat{\zeta}_u =
\frac{211}{8403}
$$

Chapter 5 uses the cluster rate:

$$
\hat{\zeta}_{cluster}
=
\frac{130}{8403}
$$

The cluster rate is lower because clustered exceedances are treated as part of the same extreme event.

Therefore, even though the declustered xi estimate is larger, the lower effective event rate produces lower return level estimates.

This is the key Chapter 5 insight:

> Accounting for dependence changes the effective frequency of independent extreme events.

## Comparison Across Chapters 3, 4, and 5

The project progression so far is:

| Chapter | Method | Model | Data object | Dependence treatment |
|---|---|---|---|---|
| Chapter 3 | Block maxima | GEV | Monthly maximum daily losses | Implicitly assumes block maxima are approximately independent |
| Chapter 4 | Threshold exceedances | GPD | Daily losses above threshold | Treats exceedances as iid |
| Chapter 5 | Declustered threshold exceedances | GPD | Cluster maxima above threshold | Accounts for clustering |

### Shape Parameter Comparison

| Chapter | Model | xi | SE(xi) |
|---|---|---:|---:|
| Chapter 3 | GEV monthly maxima | 0.1999 | 0.0433 |
| Chapter 4 | iid GPD, 97.5% threshold | 0.2932 | 0.0944 |
| Chapter 5 | declustered GPD, run length 5 | 0.3577 | 0.1185 |

All three xi estimates are positive.

This gives a coherent result across models:

> Extreme SPY daily losses appear heavy-tailed.

The magnitude of xi increases from Chapter 3 to Chapter 5, but the uncertainty also increases. This is expected because the models use different tail samples:

- Chapter 3 uses monthly maxima
- Chapter 4 uses all threshold exceedances
- Chapter 5 uses only cluster maxima

The main conclusion is not that one xi estimate is “the true one.” The main conclusion is that all three modeling approaches point in the same direction: positive tail index, heavy-tailed losses.

## Conceptual Comparison of the Models

### Chapter 3: GEV Block Maxima

Chapter 3 is clean and classical.

It takes the maximum loss in each month:

$$
M_j = \max(L_t \text{ in month } j)
$$

Then it fits a GEV model to the monthly maxima.

Strengths:

- theoretically clean
- simple interpretation
- avoids threshold selection
- directly connected to the extremal types theorem

Weaknesses:

- throws away many large losses
- one observation per month
- can miss important within-month tail information

### Chapter 4: iid GPD Threshold Model

Chapter 4 keeps all losses above a high threshold.

Strengths:

- more data-efficient than block maxima
- uses more tail observations
- directly models exceedances
- return levels are stable across thresholds

Weaknesses:

- threshold choice matters
- assumes exceedances are approximately independent
- financial data often violate this independence assumption

### Chapter 5: Declustered GPD Model

Chapter 5 keeps the threshold idea but corrects for clustering.

Strengths:

- more realistic for financial extremes
- accounts for crisis clustering
- estimates the extremal index
- avoids overcounting clustered exceedances as independent events

Weaknesses:

- fewer observations after declustering
- run length choice matters
- larger standard errors
- return levels depend on cluster rate

## Interpretation of the Extremal Index

The main extremal index estimate is:

$$
\hat{\theta} \approx 0.616
$$

Since this is below 1, there is evidence that extreme SPY losses cluster over time.

A rough interpretation is:

```text
Only about 61.6% of raw exceedances behave like independent extreme events.
```

Equivalently, the average cluster size is approximately:

$$
\frac{1}{0.616} \approx 1.62
$$

This means that an extreme episode contains about 1.62 exceedances on average under the run length 5 rule.

## Main Findings

1. The 97.5% threshold identifies 211 extreme daily losses out of 8403 observations.

2. Runs declustering with run length 5 groups these exceedances into 130 clusters.

3. The estimated extremal index is approximately 0.616, indicating clustering in extreme SPY losses.

4. The Chapter 4 iid GPD model gives xi approximately 0.293.

5. The Chapter 5 declustered GPD model gives xi approximately 0.358.

6. Both models suggest heavy-tailed losses.

7. Declustered return levels are lower than iid return levels because the effective independent event rate is lower after clustering.

8. The heavy-tail conclusion remains consistent across Chapters 3, 4, and 5.

## Main Conclusion

Chapter 5 adds an important insight to the project:

> Extreme SPY daily losses are not only heavy-tailed, they also cluster over time.

Chapter 3 showed heavy-tailed behavior using monthly maxima.

Chapter 4 showed that the heavy-tail result survives under a threshold exceedance model.

Chapter 5 shows that threshold exceedances are temporally dependent and should not be blindly treated as independent observations.

The resulting story is:

```text
Chapter 3: Extreme losses are heavy-tailed under a GEV block maxima model.
Chapter 4: Extreme losses are heavy-tailed under a GPD threshold model.
Chapter 5: Extreme losses are heavy-tailed and clustered over time.
```

That is the key progression.

## Files

This folder contains:

```text
03_chapter5_dependent_sequences/
├── README.md
├── chapter5_dependent_sequences.R
├── figures/
└── results/
```

## Key Result Files

```text
results/threshold_summary.csv
results/declustering_results.csv
results/clusters_run_length_1.csv
results/clusters_run_length_3.csv
results/clusters_run_length_5.csv
results/clusters_run_length_10.csv
results/gpd_iid_vs_declustered_comparison.csv
results/return_level_iid_vs_declustered.csv
results/chapter3_chapter4_chapter5_shape_comparison.csv
results/chapter5_model_summary.txt
```

## Key Figures

```text
figures/daily_losses_with_exceedances.png
figures/acf_absolute_returns.png
figures/acf_exceedance_indicator.png
figures/declustered_exceedances_run_length_5.png
figures/cluster_size_distribution_run_length_5.png
figures/extremal_index_by_run_length.png
figures/gpd_diagnostics_declustered_run_length_5.png
figures/return_levels_iid_vs_declustered.png
figures/chapter3_chapter4_chapter5_shape_comparison.png
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