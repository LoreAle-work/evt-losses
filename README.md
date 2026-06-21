# EVT on SPY Daily Losses

This project applies Extreme Value Theory to daily SPY returns using the block maxima approach from Chapter 3 of Coles' *An Introduction to Statistical Modeling of Extreme Values*.

## Objective

The goal is to model extreme negative daily stock returns using the Generalized Extreme Value distribution.

## Data

The analysis uses SPY adjusted closing prices downloaded from Yahoo Finance through the `quantmod` package in R.

## Method

Daily log returns are computed as:

```math
R_t = \log(P_t) - \log(P_{t-1})