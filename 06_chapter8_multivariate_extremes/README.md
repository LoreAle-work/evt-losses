# Chapter 8: Multivariate Extremes

## Joint Extreme Losses in SPY, QQQ, TLT, and VIX

This chapter extends the project from **univariate extreme-value modeling** to **multivariate extremes**.

Earlier chapters studied the tail behavior of a single financial loss series, mainly SPY daily losses. Chapter 8 asks a broader question:

> When one financial variable becomes extreme, how likely is another variable to become extreme at the same time?

This matters because financial crises are rarely one-dimensional. Equity losses, bond stress, and volatility spikes can occur together, but not all pairs behave in the same way. The goal of this chapter is to measure and compare **joint tail dependence** across different financial relationships.

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
| 2 | SPY-VIX | Substantial crash-volatility dependence |
| 3 | SPY-TLT | Very weak joint tail dependence |

The chapter supports the central idea of multivariate EVT:

> Marginal tail risk is not enough. To understand financial crises, we also need to understand how extremes occur together.

---

## Summary

This chapter compares three types of joint financial extremes.

SPY-QQQ represents equity-equity crash dependence. SPY-TLT represents equity-bond joint stress. SPY-VIX represents equity crash-volatility spike dependence.

The results show that:

- SPY-QQQ has the strongest joint tail dependence.
- SPY-VIX has substantial but less balanced crash-volatility dependence.
- SPY-TLT has very weak joint tail dependence.

The ranking is stable across:

- empirical tail-dependence estimates,
- joint exceedance ratios,
- conditional co-exceedance probabilities,
- angular diagnostics,
- Fréchet-scale structure summaries.

This shows why ordinary correlation is not enough. Correlation describes average co-movement. Multivariate EVT studies dependence in the tail, where financial risk is usually most important.

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
| Main marginal threshold | 97.5% |
| Exceedances per marginal series | 151 |

A common sample is used across all pairs so that the pairwise results are comparable.

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

The final variables are:

| Variable | Definition | Interpretation |
|---|---|---|
| SPY_extreme | negative SPY log return | large SPY loss |
| QQQ_extreme | negative QQQ log return | large QQQ loss |
| TLT_extreme | negative TLT log return | large TLT loss |
| VIX_extreme | positive VIX log return | large VIX increase |

This orientation is important. A large positive value always means financial stress.

---

## 3. Why multivariate extremes?

A univariate EVT model can estimate the probability of large SPY losses. But it cannot answer questions such as:

- Are SPY and QQQ likely to crash together?
- Does TLT experience extreme losses on the same days as SPY?
- How often do large SPY losses coincide with large VIX increases?
- Are joint extremes balanced, or does one variable dominate the event?

Chapter 8 introduces tools for studying these questions.

The central object is not ordinary correlation. The central object is **extremal dependence**.

Ordinary correlation measures average linear co-movement. Extremal dependence measures whether very large observations occur together.

This distinction matters because two assets can be moderately correlated on average but strongly linked during crises. They can also be correlated in normal periods but weakly dependent in the joint tail.

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

The marginal GPD fits are:

| Asset | Threshold | Threshold (%) | Exceedances | Scale $$\sigma$$ | Shape $$\xi$$ | SE of $$\xi$$ |
|---|---:|---:|---:|---:|---:|---:|
| SPY | 0.024386 | 2.4386% | 151 | 0.009930 | 0.2420 | 0.1080 |
| QQQ | 0.029485 | 2.9485% | 151 | 0.009656 | 0.1431 | 0.0887 |
| TLT | 0.018023 | 1.8023% | 151 | 0.004814 | 0.1548 | 0.0820 |
| VIX | 0.163193 | 16.3193% | 151 | 0.068502 | 0.1352 | 0.1001 |

All estimated marginal shape parameters are positive, although the standard errors are not negligible. This is consistent with heavy-tailed marginal stress behavior.

The marginal results are not the main focus of Chapter 8, but they are useful because multivariate extremes require careful treatment of the individual margins before analyzing dependence.

---

### 4.2 Empirical transformation to standard Fréchet margins

Multivariate EVT separates marginal behavior from dependence structure. To focus on dependence, the margins are transformed to approximately standard Fréchet scale.

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

This transformation makes the marginal scales comparable. After this step, the analysis focuses on dependence rather than differences in raw units.

This is especially important for SPY-VIX, because SPY losses and VIX increases are not naturally on the same raw scale.

---

### 4.3 Empirical extremal dependence coefficient

The main tail-dependence diagnostic is:

$$
\chi(u) = P(U_Y > u \mid U_X > u)
$$

Empirically, because both variables are approximately uniform after ranking,

$$
\chi(u) =
\frac{P(U_X > u, U_Y > u)}{1-u}
$$

If $$\chi(u)$$ remains high as $$u$$ increases, this suggests strong extremal dependence over the observed tail range.

If $$\chi(u)$$ is close to zero, this suggests weak joint tail dependence.

The chapter reports:

$$\hat{\chi}=\text{average of } \chi(u) \text{ for } u \geq 0.95$$

This is an empirical high-threshold summary. It is not a proof of asymptotic dependence.

---

### 4.4 Residual tail dependence diagnostic

The second diagnostic is:

$$\bar{\chi}(u)=\frac{2 \log(1-u)}{\log P(U_X > u, U_Y > u)}- 1$$

This statistic helps diagnose residual association in the joint tail, especially when $$\chi(u)$$ is small.

Broadly:

- $$\chi(u)$$ focuses on asymptotic dependence strength.
- $$\bar{\chi}(u)$$ helps describe residual dependence when full asymptotic dependence may not hold.

The project uses both diagnostics graphically.

---

### 4.5 Implied logistic dependence parameter

For the bivariate logistic extreme-value model, the relationship between the tail dependence coefficient and the logistic dependence parameter is:

$$
\chi = 2 - 2^\alpha
$$

Solving for $$\alpha$$ gives:

$$
\alpha =
\frac{\log(2-\chi)}{\log(2)}
$$

This project reports an empirical implied logistic alpha:

$$\hat{\alpha}=
\frac{\log(2-\hat{\chi})}{\log(2)}
$$

Interpretation:

| Value of $$\alpha$$ | Interpretation |
|---|---|
| close to 0 | strong extremal dependence |
| close to 1 | weak extremal dependence |

Important note: this is **not** a full bivariate likelihood estimate. It is an implied dependence summary derived from empirical $$\hat{\chi}$$.

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

$$\text{Joint ratio}=
\frac{p_{XY}}{p_X p_Y}
$$

A ratio greater than 1 means joint extremes occur more often than independence would predict.

A ratio close to 1 means the pair behaves close to independence at that threshold.

---

### 4.7 Conditional co-exceedance

The conditional co-exceedance probability is:

$$
P(Y > q_Y \mid X > q_X)
$$

Since the same marginal quantile level is used for both variables,

$$
P(X > q_X) = P(Y > q_Y)
$$

Therefore:

$$P(Y > q_Y \mid X > q_X)=
P(X > q_X \mid Y > q_Y)
$$

These probabilities should therefore be interpreted as **co-exceedance strength**, not directional asymmetry.

Directional structure is studied separately using angular diagnostics.

This distinction matters. Otherwise, one might incorrectly interpret equal conditional probabilities as a causal or directional statement. The equality is mechanical because the same quantile level is used for both margins.

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

The angular analysis is useful because two pairs can have similar co-exceedance probabilities but different structures. One pair may have balanced joint extremes, while another may have radial extremes dominated by one variable.

---

### 4.9 Bootstrap uncertainty

The project reports both IID and block bootstrap intervals for $$\hat{\chi}$$.

The IID bootstrap resamples individual observations. This is simple but ignores serial dependence.

The block bootstrap resamples blocks of consecutive observations. The block length is 5 trading days.

The block bootstrap is preferred because financial extremes cluster over time. This connects directly to Chapter 5, where the project found clustering in extreme SPY losses.

---

### 4.10 Structure variables

Two raw structure variables are computed:

$$
Z_{\max} = \max(X, Y)
$$

and:

$$
Z_{\min} = \min(X, Y)
$$

The maximum structure variable measures worst-component stress.

The minimum structure variable measures joint stress, because it is large only when both variables are large.

A GPD is fitted to these structure variables above their 97.5% thresholds.

For cross-pair comparison, Fréchet-scale structure summaries are also reported. This is especially important for SPY-VIX, because VIX increases and ETF losses are not directly comparable in raw units.

---

## 5. Main results

### 5.1 Pair summary

| Pair | Pearson correlation | Interpretation |
|---|---:|---|
| SPY-QQQ | 0.9133 | Strong average equity-equity co-movement |
| SPY-TLT | -0.3032 | Negative average equity-bond dependence |
| SPY-VIX | 0.7356 | Strong average crash-volatility co-movement |

Correlation already suggests that the three pairs behave differently. But correlation is not enough. The EVT diagnostics show how these relationships behave specifically in the tail.

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

## 6. How to read the main diagnostics

This chapter uses several diagnostics because multivariate extremes cannot be summarized well by one number. Humanity tried doing everything with one number, and then risk management became interpretive dance.

Each diagnostic answers a different question.

| Diagnostic | Main question answered |
|---|---|
| Pearson correlation | Do the two variables move together on average? |
| Marginal GPD shape $$\xi$$ | Are the individual stress variables heavy-tailed? |
| $$\chi(u)$$ | Do very large values occur together? |
| $$\hat{\chi}$$ | What is the average high-threshold extremal dependence? |
| Implied $$\alpha$$ | How strong is tail dependence under a logistic-style summary? |
| Joint exceedance ratio | How much more frequent are joint extremes than under independence? |
| Conditional co-exceedance | When one variable is extreme, how often is the other also extreme? |
| Angular variable $$W$$ | Are joint extremes balanced or dominated by one variable? |
| Structure variables | How large are worst-component and joint-stress events? |

The key idea is:

> Correlation is not tail dependence.

A pair can have high correlation but still have a different tail structure from another highly correlated pair. A pair can also have moderate average dependence but very strong crisis dependence.

---

## 7. Correlation versus extremal dependence

The three pairs show three different dependence regimes.

| Pair | Pearson correlation | $$\hat{\chi}$$ | Interpretation |
|---|---:|---:|---|
| SPY-QQQ | 0.9133 | 0.7105 | Strong average and strong tail dependence |
| SPY-TLT | -0.3032 | 0.0451 | Negative average dependence and weak joint tail dependence |
| SPY-VIX | 0.7356 | 0.4630 | Strong crash-volatility dependence |

SPY-QQQ is the simplest case. Correlation is high, and tail dependence is also high. Both the average relationship and the extreme relationship point in the same direction.

SPY-TLT is the contrast case. The correlation is negative, and the estimated tail dependence is close to zero. Extreme losses in SPY and TLT rarely occur together relative to SPY-QQQ and SPY-VIX.

SPY-VIX is the stress-volatility case. The correlation is positive and high because VIX tends to rise when SPY falls. The tail-dependence estimate confirms that this relationship remains important in the joint tail. However, the angular results show that SPY-VIX extremes are less balanced than SPY-QQQ extremes.

This is why Chapter 8 matters. A single correlation matrix would miss much of the structure of joint extremes.

---

## 8. Interpretation of $$\chi(u)$$

The empirical $$\chi(u)$$ curves are one of the most important outputs of the chapter.

The estimate:

$$
\chi(u) = P(U_Y > u \mid U_X > u)
$$

measures the probability that one variable is extreme given that the other variable is extreme.

The plot `comparison_chi_all_pairs.png` shows a clear ranking:

```text
SPY-QQQ highest
SPY-VIX intermediate
SPY-TLT lowest
