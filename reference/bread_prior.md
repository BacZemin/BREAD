# Construct a BREAD prior for the summary-mode model

Default prior for
[`fit_bread_summary()`](https://baczemin.github.io/BREAD/reference/fit_bread_summary.md):
weakly-informative Normal-Inverse-Gamma. `mu0` defaults to a zero vector
(length inferred from the design) and `Lambda0` to `lambda0 * I`
(diagonal prior precision, inferred at fit time). The inverse-gamma
hyperparameters default to `a0 = b0 = 0.001`, which is approximately
Jeffreys.

## Usage

``` r
bread_prior(mu0 = NULL, Lambda0 = NULL, lambda0 = 0.01, a0 = 0.001, b0 = 0.001)
```

## Arguments

- mu0:

  Optional prior mean vector for the coefficients. If `NULL`, a zero
  vector of the correct dimension is used at fit time.

- Lambda0:

  Optional prior precision matrix. If `NULL`, `lambda0 * I` is used at
  fit time.

- lambda0:

  Scalar prior precision used when `Lambda0` is `NULL`. Default `0.01`
  (weak).

- a0, b0:

  Inverse-gamma hyperparameters for the residual variance. Defaults
  `a0 = b0 = 0.001`.

## Value

A list with class `"bread_prior"`.
