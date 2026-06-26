# Chapter 8: Multivariate Extremes

## Joint Extreme Losses in SPY, QQQ, TLT, and VIX

This chapter extends the project from **univariate extreme-value modeling** to **multivariate extremes**. Earlier chapters studied the tail behavior of a single financial loss series, mainly SPY daily losses. Chapter 8 asks a different question:

> When one financial variable is extreme, how likely is another variable to be extreme at the same time?

This matters because financial crises are rarely isolated one-dimensional events. Extreme equity losses, bond stress, and volatility spikes can occur together, but not all pairs behave the same way. Chapter 8 is about measuring and comparing those joint tail relationships.

The analysis focuses on three pairs:

| Pair | Interpretation |
|---|---|
| SPY-QQQ | Equity-equity joint crash risk |
| SPY-TLT | Equity-bond joint stress behavior |
| SPY-VIX | Equity crash versus volatility spike dependence |

The main conclusion is that extremal dependence differs sharply across pairs:

| Rank | Pair | Main finding |
|---:|---|---|
| 1 | SPY-QQQ | Strongest joint tail dependence |
| 2 | SPY-VIX | Strong crash-volatility dependence |
| 3 | SPY-TLT | Very weak joint tail dependence |

This supports the main Chapter 8 lesson: **ordinary correlation and marginal tail risk are not enough to describe joint extremes**. The dependence structure in the tail matters.

---

## 1. Data

The analysis uses daily adjusted prices or index levels for:

| Symbol | Description | Extreme direction |
|---|---|---|
| SPY | S&P 500 ETF | Daily loss |
| QQQ | Nasdaq-100 ETF | Daily loss |
| TLT | Long-term Treasury ETF | Daily loss |
| VIX | CBOE Volatility Index | Daily increase |

The common sample begins when all four series are available.

| Quantity | Value |
|---|---:|
| Start date | 2002-07-31 |
| End date | 2026-06-25 |
| Observations | 6014 |
| Main threshold | 97.5% marginal quantile |
| Exceedances per margin at 97.5% | 151 |

The analysis uses a common sample across all assets so that cross-pair comparisons are based on the same dates.

---

## 2. Variable construction

Daily log returns are computed as:

$$
r_t = \log(P_t) - \log(P_{t-1})
$$

For ETFs, extreme risk is defined as a daily loss:

$$
X_t = -r_t
$$

For VIX, extreme risk is defined as a daily increase:

$$
X_t = r_t
$$

Therefore, all variables are oriented so that **larger positive values represent more extreme stress**.

The final extreme-direction variables are:

| Variable | Definition |
|---|---|
| SPY_extreme | negative SPY log return |
| QQQ_extreme | negative QQQ log return |
| TLT_extreme | negative TLT log return |
| VIX_extreme | positive VIX log return |

This orientation is important. A large positive SPY value means a large SPY loss, while a large positive VIX value means a large VIX increase.

---

## 3. Why multivariate extremes?

A univariate EVT model can estimate the probability of large SPY losses. But it cannot answer questions like:

- Are SPY and QQQ likely to crash together?
- Does TLT experience extreme losses on the same days as SPY?
- How often do large SPY losses coincide with large VIX increases?
- Are joint extremes balanced, or does one variable dominate the event?

Chapter 8 introduces tools for this type of joint tail analysis.

The central object is not ordinary correlation, but **extremal dependence**.

Ordinary correlation measures average linear co-movement. Extremal dependence measures whether very large observations occur together.

---

## 4. Methodology

### 4.1 Marginal GPD fits

Each marginal series is first studied using a generalized Pareto distribution above the 97.5% threshold.

For a high threshold $$u$$, exceedances are modeled as:

$$
X - u \mid X > u \sim \text{GPD}(\sigma, \xi)
$$

where:

- $$\sigma$$ is the scale parameter,
- $$\xi$$ is the shape parameter.

The estimated marginal GPD parameters are:

| Asset | Threshold | Threshold (%) | Exceedances | GPD scale | GPD shape $$\xi$$ |
|---|---:|---:|---:|---:|---:|
| SPY | 0.024386 | 2.4386% | 151 | 0.009930 | 0.2420 |
| QQQ | 0.029485 | 2.9485% | 151 | 0.009656 | 0.1431 |
| TLT | 0.018023 | 1.8023% | 151 | 0.004814 | 0.1548 |
| VIX | 0.163193 | 16.3193% | 151 | 0.068502 | 0.1352 |

All estimated marginal shape parameters are positive, although uncertainty remains. This is consistent with heavy-tailed marginal stress variables.

---

### 4.2 Empirical transformation to standard Fréchet margins

Multivariate EVT separates marginal tail behavior from dependence structure. To focus on dependence, the margins are transformed to approximately standard Fréchet scale.

First, empirical ranks are converted to pseudo-uniform variables:

$$
U_i = \frac{\text{rank}(X_i)}{n + 1}
$$

Then the standard Fréchet transformation is:

$$
Z_i = -\frac{1}{\log(U_i)}
$$

The standard Fréchet distribution has distribution function:

$$
F(z) = \exp\left(-\frac{1}{z}\right), \quad z > 0
$$

This transformation makes the marginal scales comparable and allows the analysis to focus on extremal dependence.

---

### 4.3 Empirical extremal dependence coefficient

The main tail-dependence diagnostic is:

$$
\chi(u) = P(U_Y > u \mid U_X > u)
$$

Empirically, since both marginal variables are approximately uniform,

$$
\chi(u) =
\frac{P(U_X > u, U_Y > u)}{1-u}
$$

If $$\chi(u)$$ remains high as $$u$$ increases, this suggests strong extremal dependence over the observed tail range.

If $$\chi(u)$$ moves toward zero, this suggests weak joint tail dependence or asymptotic independence.

The summary statistic used in this project is:

$$
\hat{\chi} =
\text{average of } \chi(u) \text{ for } u \geq 0.95
$$

This is an empirical high-threshold summary, not a proof of asymptotic dependence.

---

### 4.4 Residual tail dependence diagnostic

The second diagnostic is:

$$
\bar{\chi}(u)
=
\frac{2 \log(1-u)}
{\log P(U_X > u, U_Y > u)}
- 1
$$

This statistic helps diagnose the strength of residual association in the joint tail, especially when $$\chi(u)$$ is near zero.

---

### 4.5 Implied logistic dependence parameter

For the bivariate logistic extreme value model, the relationship between the tail dependence coefficient and the logistic dependence parameter is:

$$
\chi = 2 - 2^\alpha
$$

Solving for $$\alpha$$ gives:

$$
\alpha =
\frac{\log(2-\chi)}{\log(2)}
$$

This project reports an **implied logistic alpha**:

$$
\hat{\alpha}
=
\frac{\log(2-\hat{\chi})}{\log(2)}
$$

Smaller values of $$\alpha$$ indicate stronger extremal dependence. Values near 1 indicate weak extremal dependence.

Important methodological note: this is **not** a full bivariate likelihood estimate. It is an implied dependence summary derived from empirical $$\hat{\chi}$$.

---

### 4.6 Joint exceedance ratio

For each pair and quantile level, the empirical joint exceedance probability is compared with the probability expected under independence.

Let:

$$
p_X = P(X > q_X)
$$

$$
p_Y = P(Y > q_Y)
$$

$$
p_{XY} = P(X > q_X, Y > q_Y)
$$

Under independence:

$$
P(X > q_X, Y > q_Y) = p_X p_Y
$$

The joint exceedance ratio is:

$$
\text{Joint ratio}
=
\frac{p_{XY}}{p_X p_Y}
$$

A ratio greater than 1 means joint extremes occur more often than independence would predict.

---

### 4.7 Conditional co-exceedance

The conditional co-exceedance probability is:

$$
P(Y > q_Y \mid X > q_X)
$$

Because the analysis uses equal marginal quantile levels, we have:

$$
P(X > q_X) = P(Y > q_Y)
$$

Therefore:

$$
P(Y > q_Y \mid X > q_X)
=
P(X > q_X \mid Y > q_Y)
$$

This means the conditional probabilities are mechanically equal. They should be interpreted as **co-exceedance strength**, not as directional asymmetry.

Directional structure is instead analyzed using angular diagnostics.

---

### 4.8 Angular and radial decomposition

After transforming both variables to standard Fréchet scale, define:

$$
R = Z_X + Z_Y
$$

and

$$
W = \frac{Z_X}{Z_X + Z_Y}
$$

Here:

- $$R$$ measures the overall size of the bivariate extreme,
- $$W$$ measures the direction of the extreme.

Interpretation of $$W$$:

| Value of $$W$$ | Meaning |
|---|---|
| close to 0 | Y dominates the extreme |
| close to 0.5 | both variables contribute similarly |
| close to 1 | X dominates the extreme |

Radial extremes are defined as observations where $$R$$ exceeds a high threshold. The main radial threshold is the 95% quantile of $$R$$.

The angular analysis reports:

| Quantity | Interpretation |
|---|---|
| share of $$W \in [0.4, 0.6]$$ | balanced joint extremes |
| share of $$W < 0.25$$ | Y-dominated extremes |
| share of $$W > 0.75$$ | X-dominated extremes |
| endpoint-heavy share | one variable dominates the radial extreme |
| endpoint asymmetry | difference between X-dominated and Y-dominated endpoint mass |

---

### 4.9 Bootstrap uncertainty

The project reports both IID and block bootstrap intervals for $$\hat{\chi}$$.

The IID bootstrap resamples individual observations. This is simple but ignores serial dependence.

The block bootstrap resamples blocks of consecutive observations. The block length is 5 trading days, matching the declustering logic used earlier in Chapter 5.

The main interpretation uses the **block bootstrap confidence interval**, because financial extremes cluster over time.

---

### 4.10 Structure variables

Two raw structure variables are computed:

$$
Z_{\max} = \max(X, Y)
$$

and

$$
Z_{\min} = \min(X, Y)
$$

The maximum structure variable measures the worst single-variable stress event. The minimum structure variable measures joint stress because it is large only when both variables are large.

A GPD is fitted to these raw structure variables above their 97.5% thresholds.

For cross-pair comparison, especially for SPY-VIX, Fréchet-scale structure summaries are also reported. This matters because raw SPY losses and VIX increases are not directly comparable in scale.

---

## 5. Main results

### 5.1 Pair summary

| Pair | Pearson correlation | Interpretation |
|---|---:|---|
| SPY-QQQ | 0.9133 | Strong average equity-equity co-movement |
| SPY-TLT | -0.3032 | Negative average equity-bond dependence |
| SPY-VIX | 0.7356 | Strong average crash-volatility co-movement |

Correlation already suggests that the three pairs behave differently. But the EVT analysis shows that tail dependence gives a sharper picture.

---

### 5.2 Tail dependence ranking

| Rank | Pair | $$\hat{\chi}$$ | Block bootstrap 95% CI | Implied $$\alpha$$ |
|---:|---|---:|---:|---:|
| 1 | SPY-QQQ | 0.7105 | [0.5920, 0.8368] | 0.3669 |
| 2 | SPY-VIX | 0.4630 | [0.3792, 0.5567] | 0.6201 |
| 3 | SPY-TLT | 0.0451 | [0.0214, 0.0727] | 0.9671 |

Interpretation:

- SPY-QQQ has the strongest extremal dependence.
- SPY-VIX also has substantial extremal dependence.
- SPY-TLT has very weak extremal dependence.

The implied logistic alpha gives the same ranking. SPY-QQQ has the smallest implied alpha, meaning strongest tail dependence. SPY-TLT has an implied alpha close to 1, meaning weak tail dependence.

---

### 5.3 Joint exceedance ratios

At the 97.5% marginal threshold:

| Pair | Joint exceedance ratio | Conditional co-exceedance | Joint observations |
|---|---:|---:|---:|
| SPY-QQQ | 27.17 | 0.6821 | 103 |
| SPY-VIX | 16.88 | 0.4238 | 64 |
| SPY-TLT | 1.58 | 0.0397 | 6 |

This means:

- SPY-QQQ joint extremes occur about 27 times more often than independence would predict.
- SPY-VIX joint extremes occur about 17 times more often than independence would predict.
- SPY-TLT joint extremes occur only slightly more often than independence at 97.5%, with very few joint observations.

The 99% threshold results are treated as robustness checks because joint counts become small. For example, SPY-TLT has only 5 joint observations at the 99% threshold, so its high 99% joint ratio is unstable.

---

### 5.4 Conditional co-exceedance across thresholds

| Pair | 90% | 95% | 97.5% | 99% |
|---|---:|---:|---:|---:|
| SPY-QQQ | 0.7525 | 0.7209 | 0.6821 | 0.6721 |
| SPY-TLT | 0.0930 | 0.0498 | 0.0397 | 0.0820 |
| SPY-VIX | 0.6130 | 0.5249 | 0.4238 | 0.3279 |

SPY-QQQ remains strongly dependent across all thresholds. SPY-VIX also remains meaningfully dependent, though the co-exceedance probability decreases at more extreme thresholds. SPY-TLT remains weak and unstable because joint exceedance counts are small.

---

## 6. Angular results

### 6.1 Main angular summary

At the 95% radial threshold:

| Pair | Balanced share $$W \in [0.4,0.6]$$ | Endpoint-heavy share | Endpoint asymmetry |
|---|---:|---:|---:|
| SPY-QQQ | 0.3654 | 0.2126 | -0.0266 |
| SPY-VIX | 0.1661 | 0.5648 | 0.0199 |
| SPY-TLT | 0.0299 | 0.9435 | 0.0000 |

Interpretation:

- SPY-QQQ has the most balanced joint extremes.
- SPY-VIX has meaningful joint tail dependence but more endpoint dominance.
- SPY-TLT radial extremes are overwhelmingly endpoint-heavy.

The endpoint asymmetry values are small. This means the main angular difference is not strong directional asymmetry. The key difference is **balanced joint extremes versus endpoint-dominated extremes**.

---

### 6.2 Angular sensitivity

The angular analysis is repeated at radial thresholds of 90%, 95%, and 97.5%.

| Pair | Radial threshold | Balanced share | Endpoint-heavy share |
|---|---:|---:|---:|
| SPY-QQQ | 90% | 0.3937 | 0.1711 |
| SPY-QQQ | 95% | 0.3654 | 0.2126 |
| SPY-QQQ | 97.5% | 0.3377 | 0.2318 |
| SPY-TLT | 90% | 0.0365 | 0.9186 |
| SPY-TLT | 95% | 0.0299 | 0.9435 |
| SPY-TLT | 97.5% | 0.0066 | 0.9735 |
| SPY-VIX | 90% | 0.2209 | 0.4518 |
| SPY-VIX | 95% | 0.1661 | 0.5648 |
| SPY-VIX | 97.5% | 0.1523 | 0.6623 |

The ranking is stable:

```text
SPY-QQQ has the most balanced joint extremes.
SPY-VIX is intermediate.
SPY-TLT is overwhelmingly endpoint-dominated.