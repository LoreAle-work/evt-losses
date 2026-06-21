# EVT on SPY Daily Losses

This project applies Extreme Value Theory to daily SPY returns using the block maxima approach from Chapter 3 of Coles' *An Introduction to Statistical Modeling of Extreme Values*.

## Objective

The goal is to model extreme negative daily stock returns using the Generalized Extreme Value distribution.

## Data

The analysis uses SPY adjusted closing prices downloaded from Yahoo Finance through the `quantmod` package in R.

## Methodology

This project applies **Extreme Value Theory (EVT)** to model unusually large daily losses of SPY. Instead of modeling the entire return distribution, the analysis focuses on the extreme left tail, where the largest market losses occur.

### 1. Daily log returns and losses

Let \(P_t\) denote the adjusted closing price of SPY on day \(t\). Daily log returns are computed as:

$$
R_t = \log(P_t) - \log(P_{t-1})
$$

Because the objective is to study large negative returns, returns are transformed into losses:

$$
L_t = -R_t
$$

This transformation turns large negative returns into large positive losses, allowing the problem to be studied as an upper-tail extreme value problem.

### 2. Block maxima

The daily loss series is divided into monthly blocks. For each month, the largest daily loss is extracted:

$$
M_n = \max(L_1, L_2, \ldots, L_n)
$$

where \(M_n\) represents the maximum daily loss within a block of size \(n\).

Monthly blocks are used as the main specification because they provide more observations than annual blocks while still focusing on extreme losses. Annual block maxima are also used as a robustness check.

### 3. Generalized Extreme Value distribution

The block maxima approach is based on the extremal types theorem. If suitably normalized block maxima converge in distribution, their limiting distribution must belong to the **Generalized Extreme Value (GEV)** family:

$$
G(z) =
\exp \left\{
-\left[
1 + \xi \left( \frac{z-\mu}{\sigma} \right)
\right]^{-1/\xi}
\right\}
$$

defined on the support:

$$
1 + \xi \left( \frac{z-\mu}{\sigma} \right) > 0
$$

The parameters are:

- \(\mu\): location parameter
- \(\sigma > 0\): scale parameter
- \(\xi\): shape parameter

The shape parameter controls the tail behavior:

- \(\xi > 0\): Fréchet type, heavy-tailed distribution
- \(\xi = 0\): Gumbel type, light-tailed distribution
- \(\xi < 0\): Weibull type, bounded upper tail

For financial losses, the key object is \(\xi\). A positive estimate of \(\xi\) suggests heavy-tailed extreme losses, meaning very large losses are more likely than under thin-tailed models such as the normal distribution.

### 4. Maximum likelihood estimation

The GEV parameters are estimated by maximum likelihood using the monthly block maxima:

$$
\hat{\theta} = (\hat{\mu}, \hat{\sigma}, \hat{\xi})
$$

where \(\hat{\theta}\) denotes the estimated parameter vector.

The fitted model is then used to study the tail behavior of SPY losses and to estimate return levels.

### 5. Return levels

An \(m\)-block return level \(z_m\) is the loss level expected to be exceeded once every \(m\) blocks on average. It satisfies:

$$
P(M > z_m) = \frac{1}{m}
$$

Equivalently:

$$
G(z_m) = 1 - \frac{1}{m}
$$

For \(\xi \neq 0\), the GEV return level is:

$$
z_m =
\mu - \frac{\sigma}{\xi}
\left[
1 - \left\{ -\log\left(1 - \frac{1}{m}\right) \right\}^{-\xi}
\right]
$$

This project reports 1-year, 5-year, and 10-year return levels using monthly block maxima, corresponding to:

$$
m = 12,\quad m = 60,\quad m = 120
$$

### 6. Model diagnostics

The fitted GEV model is evaluated using standard diagnostic plots:

- probability plot
- quantile plot
- return level plot
- density plot

These diagnostics assess whether the GEV model provides a reasonable fit to the observed monthly maxima, especially in the upper tail.

### 7. Robustness checks

Two robustness checks are included.

First, the GEV model is re-estimated using yearly block maxima. This is closer to the classical block maxima approach but produces fewer observations, so estimates are less stable.

Second, the analysis uses the \(r\)-largest order statistics approach for:

$$
r = 1, 2, 3
$$

Instead of keeping only the largest loss in each month, this approach also considers the second and third largest losses. The aim is to check whether the estimated shape parameter remains stable when more extreme observations are included.

If \(\hat{\xi}\) remains positive and similar across monthly, yearly, and \(r\)-largest specifications, this supports the conclusion that SPY extreme losses exhibit heavy-tailed behavior.
