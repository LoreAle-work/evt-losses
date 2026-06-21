# Extreme Value Theory for Financial Losses

This repository is a progressive applied project on **Extreme Value Theory (EVT)** using financial market data.

The goal is to study how extreme financial losses can be modeled statistically, starting from classical block maxima methods and gradually moving toward more advanced EVT techniques.

The project follows the structure of Stuart Coles' *An Introduction to Statistical Modeling of Extreme Values*. Each stage introduces a new EVT method, explains the mathematical idea behind it, and applies it to financial losses using R.

The empirical focus is on SPY, an ETF tracking the S&P 500, used here as a broad proxy for the U.S. equity market.

## Project Objective

Standard financial models often describe average behavior, volatility, or typical return fluctuations. However, financial crises are not average events. They occur in the tails of the distribution.

This project focuses on the statistical modeling of rare and severe market losses.

The aim is to build a clear progression:

```text
Classical EVT and block maxima
    ↓
Threshold exceedances
    ↓
Extremes of dependent sequences
    ↓
Extremes of non-stationary sequences
    ↓
Point process characterization of extremes
    ↓
Multivariate extremes
    ↓
Further EVT topics
```

Each chapter adds a more advanced layer of extreme value modeling.

## Repository Structure

```text
evt-financial-losses/
├── README.md
├── .gitignore
├── 01_chapter3_classical_evt/
│   ├── README.md
│   ├── chapter3_classical_evt.R
│   ├── figures/
│   └── results/
├── 02_chapter4_threshold_models/
│   ├── README.md
│   ├── chapter4_threshold_models.R
│   ├── figures/
│   └── results/
├── 03_chapter5_dependent_sequences/
│   ├── README.md
│   ├── chapter5_dependent_sequences.R
│   ├── figures/
│   └── results/
├── 04_chapter6_nonstationary_sequences/
│   ├── README.md
│   ├── chapter6_nonstationary_sequences.R
│   ├── figures/
│   └── results/
├── 05_chapter7_point_processes/
│   ├── README.md
│   ├── chapter7_point_processes.R
│   ├── figures/
│   └── results/
├── 06_chapter8_multivariate_extremes/
│   ├── README.md
│   ├── chapter8_multivariate_extremes.R
│   ├── figures/
│   └── results/
└── 07_chapter9_further_topics/
    ├── README.md
    ├── chapter9_further_topics.R
    ├── figures/
    └── results/
```

Each chapter folder contains:

- an R script with the full analysis
- a chapter-specific README explaining the method
- diagnostic plots and figures
- numerical results and fitted model summaries

## Data

The project uses daily adjusted closing prices for SPY.

Let $$P\_t$$ denote the adjusted closing price on day t. Daily log returns are computed as:

$$
R_t = \log(P_t) - \log(P_{t-1})
$$

where $begin:math:text$R\_t$end:math:text$ is the daily log return.

Because the focus is on market losses, returns are transformed into losses:

$$
L_t = -R_t
$$

This transformation turns large negative returns into large positive losses. Therefore, EVT can be applied to the upper tail of the loss distribution.

In other words, instead of studying extremely negative returns directly, the project studies extremely large values of the loss variable $$L\_t$$.

## Chapter 3: Classical Extreme Value Theory and Models

The first stage uses **classical extreme value theory**, focusing on the block maxima approach.

The daily loss series is divided into monthly blocks, and the largest daily loss in each month is extracted.

If $$L\_1\, L\_2\, \\ldots\, L\_n$$ are daily losses within a block, the block maximum is:

$$
M_n = \max(L_1, L_2, \ldots, L_n)
$$

where $$M\_n$$ is the largest daily loss in that block.

The sequence of block maxima is modeled using the **Generalized Extreme Value (GEV)** distribution.

The GEV distribution is:

$$
G(z) =
\exp\left(
-\left[
1 + \xi \left( \frac{z-\mu}{\sigma} \right)
\right]^{-1/\xi}
\right)
$$

with support condition:

$$
1 + \xi \left( \frac{z-\mu}{\sigma} \right) > 0
$$

The parameters are:

- $$\(\mu\)$$: location parameter
- $$\(\sigma > 0\)$$: scale parameter
- $$\(\xi\)$$: shape parameter

The shape parameter $$\(\xi\)$$ determines the behavior of the tail:

- $$\(\xi > 0\)$$: Fréchet type, heavy-tailed distribution
- $$\(\xi = 0\)$$: Gumbel type, light-tailed distribution
- $$\(\xi < 0\)$$: Weibull type, bounded upper tail

In financial risk applications, $$\(\xi\)$$ is especially important. A positive value of $$\(\xi\)$$ suggests that extreme losses are heavy-tailed, meaning very large losses are more likely than they would be under thin-tailed models such as the normal distribution.

This chapter also estimates **return levels**, which answer questions such as:

> What size daily loss is expected to be exceeded once every 1, 5, or 10 years?

For monthly block maxima, an $$m$$-month return level $$z\_m$$ satisfies:

$$
P(M > z_m) = \frac{1}{m}
$$

Equivalently:

$$
G(z_m) = 1 - \frac{1}{m}
$$

For $$\(\xi\)$$, the return level is:

$$
z_m =
\mu - \frac{\sigma}{\xi}
\left[
1 -
\left(
-\log\left(1 - \frac{1}{m}\right)
\right)^{-\xi}
\right]
$$

This chapter reports 1-year, 5-year, and 10-year return levels using:

$$
m = 12,\quad m = 60,\quad m = 120
$$

because the model is fitted to monthly maxima.

The chapter also includes the **r-largest order statistic model**, which extends the block maxima approach by using more than just the largest observation in each block.

## Chapter 4: Threshold Models

The second stage moves from block maxima to **threshold models**, also known as the peaks-over-threshold approach.

Instead of keeping only the largest observation in each block, this method keeps all losses above a high threshold $$u$$.

The exceedances over the threshold are defined as:

$$
Y = L - u \mid L > u
$$

where $$u$$ is a high threshold and $$Y$$ is the amount by which the loss exceeds that threshold.

For a sufficiently high threshold, EVT suggests that exceedances can be modeled using the **Generalized Pareto Distribution (GPD)**.

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

- $$\(\sigma\) \> 0$$: scale parameter
- $$\(\xi\)$$: shape parameter

This approach is usually more data-efficient than block maxima because it uses more extreme observations.

The main goals of this chapter are:

- choose suitable thresholds
- fit the GPD to threshold exceedances
- compare threshold-based results with block maxima results
- estimate tail probabilities
- estimate return levels
- study sensitivity to threshold choice
- perform model checking for threshold exceedances

This chapter improves on the block maxima approach by using more information from the tail rather than discarding all but one observation per block.

## Chapter 5: Extremes of Dependent Sequences

Financial losses are not perfectly independent.

Extreme losses often cluster during crisis periods, such as 2008 or 2020. This creates dependence in the extremes. Ignoring this dependence can lead to misleading risk estimates.

This chapter studies how temporal dependence affects EVT models.

The main topics include:

- maxima of stationary sequences
- modeling stationary time series extremes
- block maxima under dependence
- threshold models under dependence
- clustering of extreme losses
- declustering methods
- the extremal index
- interpretation of return levels under dependence

The **extremal index**, usually denoted by $$\(\theta\)$$, measures the degree of clustering in extremes.

A simplified interpretation is:

- $$\(\theta\)\=1$$: extremes behave approximately independently
- $$\(\theta\)\<1$$: extremes tend to cluster

The goal is to avoid treating clustered crisis losses as if they were fully independent observations. This chapter makes the modeling framework more realistic for financial data, where market stress often persists over multiple days.

## Chapter 6: Extremes of Non-Stationary Sequences

The earlier models assume that the distribution of extremes is stable over time. In finance, this assumption can be too restrictive.

Market risk may change across regimes, especially during:

- financial crises
- volatility spikes
- monetary tightening periods
- structural breaks
- changes in market liquidity

This chapter studies **non-stationary EVT models**, where GEV or GPD parameters may depend on covariates or time.

For example, a non-stationary GEV model may allow the location parameter to vary with a covariate $$x\_t$$:

$$
\mu_t = \beta_0 + \beta_1 x_t
$$

The scale parameter may also be modeled as time-varying:

$$
\log(\sigma_t) = \gamma_0 + \gamma_1 x_t
$$

The logarithm is used to ensure that $$\(\sigma\) \> 0$$.

Possible covariates include:

- time
- realized volatility
- rolling standard deviation of returns
- market regime indicators
- crisis dummy variables

The goals of this chapter are:

- allow EVT parameters to vary over time
- compare stationary and non-stationary models
- perform model choice
- evaluate model diagnostics
- study whether extreme risk changes across market regimes
- improve interpretation of return levels under changing conditions

This chapter asks whether the probability of extreme losses is constant over time or whether it changes with market conditions.

## Chapter 7: Point Process Characterization of Extremes

Chapter 7 introduces a **point process characterization of extremes**.

This framework connects block maxima models and threshold exceedance models under a common theoretical structure. Instead of looking only at maxima or only at exceedances, the point process approach models the occurrence of extreme events over both time and magnitude.

The main idea is that extremes can be represented as points:

$$
(t_i, L_i)
$$

where $$t\_i$$ records when an extreme event occurs and $$L\_i$$ records its size.

Under suitable conditions, the limiting behavior of extreme observations can be described using a Poisson process.

The goals of this chapter are:

- introduce the basic theory of point processes
- study the Poisson process limit for extremes
- connect the point process model to the GEV distribution
- connect the point process model to threshold excess models
- estimate return levels using the point process framework
- compare point process modeling with the earlier GEV and GPD approaches

This chapter is important because it shows that block maxima and threshold exceedance methods are not separate tricks. They are connected through a broader mathematical framework.

## Chapter 8: Multivariate Extremes

Chapter 8 moves from univariate losses to **multivariate extremes**.

Financial risk is rarely isolated. Large losses in one asset often occur together with large losses in others. During crises, diversification can weaken exactly when investors need it most.

This chapter studies joint tail behavior across multiple assets.

Possible applications include:

- SPY and QQQ
- SPY and sector ETFs
- SPY and financial sector ETFs
- U.S. equity indices and volatility indices
- equity losses across multiple international markets

The goal is to understand whether extreme losses occur independently or jointly.

Key questions include:

- Do extreme losses in different assets occur together?
- How strong is tail dependence?
- Does diversification fail during extreme events?
- Can multivariate EVT describe systemic risk better than separate univariate models?

A central concept is **tail dependence**, which studies the probability that one variable is extreme given that another variable is also extreme.

This chapter extends the project from modeling isolated extreme losses to studying the dependence structure of extreme financial events.

## Chapter 9: Further Topics

Chapter 9 covers further topics in extreme value modeling.

Possible extensions include:

- Bayesian inference for extremes
- Bayesian estimation of GEV or GPD parameters
- Markov chain models for extremes
- spatial extremes
- max-stable processes
- latent spatial process models

For this project, the most relevant extension is likely **Bayesian inference**, because it allows prior information and parameter uncertainty to be incorporated directly into the model.

This chapter may be used to extend earlier models by comparing frequentist and Bayesian EVT estimates.

## Why This Project Matters

Standard financial models often focus on average behavior, volatility, or normally distributed returns.

But financial crises are not average events. They live in the tails.

EVT provides tools specifically designed for rare and severe events. This makes it useful for:

- financial risk management
- stress testing
- tail-risk estimation
- crisis modeling
- return level estimation
- understanding the probability of extreme losses
- studying whether diversification breaks down during crises

This repository is intended to show not only how EVT works mathematically, but also how it can be implemented and interpreted in a real financial setting.

The project is designed to be both technical and applied: each method is connected to financial risk questions.

## Software

The project is written in R.

Main packages include:

```r
quantmod
xts
zoo
ismev
```

Additional packages may be added in later chapters as the project becomes more advanced.

## Current Status

Completed:

- Chapter 3: classical EVT, block maxima, GEV modeling, return levels, diagnostics, and the r-largest order statistic model

Planned:

- Chapter 4: threshold models and the Generalized Pareto Distribution
- Chapter 5: extremes of dependent sequences
- Chapter 6: extremes of non-stationary sequences
- Chapter 7: point process characterization of extremes
- Chapter 8: multivariate extremes
- Chapter 9: further topics, including Bayesian inference

## Project Philosophy

This repository is designed to show progression.

Each chapter builds on the previous one:

```text
Chapter 3: Classical EVT and block maxima
    ↓
Chapter 4: Threshold models
    ↓
Chapter 5: Dependent extremes
    ↓
Chapter 6: Non-stationary extremes
    ↓
Chapter 7: Point process models
    ↓
Chapter 8: Multivariate extremes
    ↓
Chapter 9: Further topics
```

The goal is not simply to run statistical models. The goal is to understand what each EVT method adds, what assumptions it makes, and how the interpretation of financial tail risk changes as the modeling framework becomes more realistic.
