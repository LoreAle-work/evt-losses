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

## Methodology

This project applies **Extreme Value Theory (EVT)** to model unusually large daily losses of SPY. Instead of modeling the entire distribution of daily returns, the analysis focuses on the extreme left tail, where the largest market losses occur.

The main idea is simple: ordinary volatility models describe typical fluctuations, while EVT focuses directly on rare and severe losses.

### 1. Daily log returns and losses

Let \(P_t\) denote the adjusted closing price of SPY on day \(t\). Daily log returns are computed as:

$$
R_t = \log(P_t) - \log(P_{t-1})
$$

Since the objective is to study large negative returns, returns are transformed into losses:

$$
L_t = -R_t
$$

This transformation turns large negative returns into large positive losses. Therefore, the problem becomes an upper-tail extreme value problem applied to the loss series.

### 2. Block maxima

The daily loss series is divided into monthly blocks. For each month, the largest daily loss is extracted:

$$
M_n = \max(L_1, L_2, \ldots, L_n)
$$

where \(M_n\) represents the maximum daily loss within a block of size \(n\).

Monthly blocks are used as the main specification because they provide a larger number of block maxima than annual blocks. This gives a better sample size for estimation while still focusing on extreme observations.

Annual block maxima are also estimated as a robustness check.

### 3. Generalized Extreme Value distribution

The block maxima method is based on the extremal types theorem. If properly normalized block maxima converge in distribution, their limiting distribution belongs to the **Generalized Extreme Value (GEV)** family.

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

In financial applications, the shape parameter is especially important. A positive estimate of $$\(\xi\)$$ suggests heavy-tailed losses, meaning that very large losses are more likely than they would be under thin-tailed models such as the normal distribution.

### 4. Maximum likelihood estimation

The GEV parameters are estimated using maximum likelihood on the monthly block maxima.

The estimated parameter vector is:

$$
\hat{\theta} =
(\hat{\mu}, \hat{\sigma}, \hat{\xi})
$$

where:

- $$\(\hat{\mu}\)$$ is the estimated location parameter
- $$\(\hat{\sigma}\)$$ is the estimated scale parameter
- $$\(\hat{\xi}\)$$ is the estimated shape parameter

The estimate of $$\(\xi\)$$ is used to determine whether the extreme losses of SPY show evidence of heavy-tailed behavior.

### 5. Return levels

The fitted GEV model is used to estimate return levels. An \(m\)-block return level \(z_m\) is the level expected to be exceeded once every \(m\) blocks on average.

It satisfies:

$$
P(M > z_m) = \frac{1}{m}
$$

Equivalently:

$$
G(z_m) = 1 - \frac{1}{m}
$$

For \(\xi \neq 0\), the return level is:

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

Since the main model uses monthly maxima, the return periods are measured in months. Therefore:

$$
m = 12
$$

corresponds to a 1-year return level,

$$
m = 60
$$

corresponds to a 5-year return level, and

$$
m = 120
$$

corresponds to a 10-year return level.

These return levels estimate the size of an extreme daily loss expected to be exceeded once over the corresponding time horizon.

### 6. Diagnostic plots

The fitted GEV model is evaluated using standard diagnostic plots:

- probability plot
- quantile plot
- return level plot
- density plot

The probability and quantile plots assess whether the fitted GEV distribution matches the empirical distribution of monthly maxima.

The return level plot shows how estimated extreme losses increase with longer return periods.

The density plot compares the fitted GEV density with the observed block maxima.

These diagnostics are especially important in the far right tail of the loss distribution, where the most severe market losses appear.

### 7. Robustness checks

Two robustness checks are included.

First, the GEV model is re-estimated using yearly block maxima. This is closer to the classical block maxima approach, but it produces far fewer observations. As a result, yearly estimates are expected to have larger uncertainty.

Second, the analysis uses the \(r\)-largest order statistics approach for:

$$
r = 1, 2, 3
$$

Instead of using only the single largest loss in each month, this approach also includes the second and third largest losses. This provides more information from each block.

The purpose of the \(r\)-largest comparison is not simply to choose the model with the best numerical fit. The goal is to check whether the estimated shape parameter remains stable as more extreme observations are included.

If $$\(\hat{\xi}\)$$ remains positive and similar across the monthly GEV model, yearly robustness check, and \(r\)-largest models, this supports the conclusion that SPY extreme losses exhibit heavy-tailed behavior.
