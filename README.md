# EVT on SPY Daily Losses

This project applies Extreme Value Theory to daily SPY returns using the block maxima approach from Chapter 3 of Coles' *An Introduction to Statistical Modeling of Extreme Values*.

## Objective

The goal is to model extreme negative daily stock returns using the Generalized Extreme Value distribution.

## Data

The analysis uses SPY adjusted closing prices downloaded from Yahoo Finance through the `quantmod` package in R.

## Methodology

This project applies **Extreme Value Theory (EVT)** to model unusually large daily losses of SPY. Instead of modeling the full distribution of returns, the focus is only on the extreme left tail, where the largest market losses occur.

### 1. Daily log returns and losses

Let \(P_t\) denote the adjusted closing price of SPY on day \(t\). Daily log returns are computed as

\[
R_t = \log(P_t) - \log(P_{t-1}).
\]

Since the interest is in large negative returns, returns are transformed into losses:

\[
L_t = -R_t.
\]

Under this transformation, large negative returns become large positive losses. This allows the analysis to focus on upper-tail extremes of the loss distribution.

### 2. Block maxima approach

The analysis follows the block maxima approach from Extreme Value Theory. The daily loss series is divided into monthly blocks, and from each month the maximum daily loss is extracted:

\[
M_n = \max(L_1, L_2, \ldots, L_n),
\]

where \(M_n\) represents the largest daily loss within a given block.

Monthly blocks are used because they provide more observations than annual blocks while still focusing on extreme events. Annual block maxima are also estimated as a robustness check.

### 3. Generalized Extreme Value distribution

According to the extremal types theorem, suitably normalized block maxima converge in distribution to the **Generalized Extreme Value (GEV)** distribution. The GEV distribution is given by

\[
G(z) =
\exp \left\{
-\left[
1 + \xi \left( \frac{z-\mu}{\sigma} \right)
\right]^{-1/\xi}
\right\},
\]

defined for

\[
1 + \xi \left( \frac{z-\mu}{\sigma} \right) > 0.
\]

The three parameters are:

- \(\mu\): location parameter
- \(\sigma > 0\): scale parameter
- \(\xi\): shape parameter

The shape parameter \(\xi\) is especially important because it determines the type of tail behavior:

- \(\xi > 0\): Fréchet type, heavy-tailed distribution
- \(\xi = 0\): Gumbel type, light-tailed distribution
- \(\xi < 0\): Weibull type, bounded upper tail

In the context of financial losses, a positive estimate of \(\xi\) suggests that extreme losses are heavy-tailed, meaning very large losses are more likely than they would be under a normal distribution.

### 4. Maximum likelihood estimation

The GEV parameters are estimated by maximum likelihood using the monthly block maxima. The fitted model provides estimates of

\[
\hat{\mu}, \quad \hat{\sigma}, \quad \hat{\xi}.
\]

The estimated shape parameter is then used to assess whether SPY losses show evidence of heavy-tailed extreme behavior.

### 5. Return levels

The fitted GEV model is also used to estimate return levels. An \(m\)-block return level \(z_m\) is the level expected to be exceeded once every \(m\) blocks on average.

For a return period of \(m\) months, the return level satisfies

\[
P(M > z_m) = \frac{1}{m}.
\]

Using the GEV model, the return level is computed as

\[
z_m =
\mu - \frac{\sigma}{\xi}
\left[
1 - \{-\log(1 - 1/m)\}^{-\xi}
\right],
\]

for \(\xi \neq 0\).

This project reports 1-year, 5-year, and 10-year return levels based on monthly block maxima.

### 6. Model diagnostics

The fitted GEV model is assessed using standard diagnostic plots:

- probability plot
- quantile plot
- return level plot
- density plot

These plots help evaluate whether the GEV distribution provides a reasonable fit to the observed monthly maxima, especially in the upper tail where the most severe losses occur.

### 7. Robustness checks

Two robustness checks are included.

First, the model is re-estimated using yearly block maxima. This is closer to the classical block maxima approach but provides fewer observations.

Second, the analysis uses the \(r\)-largest order statistics approach for \(r = 1, 2, 3\). Instead of keeping only the single largest loss in each month, this method also considers the second and third largest losses. The goal is to check whether the estimated shape parameter remains stable when more extreme observations are included.

If the estimated \(\xi\) remains positive and similar across these specifications, this supports the conclusion that SPY extreme losses exhibit heavy-tailed behavior.
