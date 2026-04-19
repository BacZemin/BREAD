# Fit BREAD summary-mode model via brms (MCMC backend)

Per-region Bayesian linear regression fitted with brms. The Stan program
is compiled **once** on the first region and reused across the remaining
regions via `update(..., recompile = FALSE)`, which makes fitting N
regions roughly "one Stan compile + N cheap sampling runs".

## Usage

``` r
fit_bread_brms(
  region_mat,
  coldata,
  design,
  contrast,
  iter = 2000L,
  chains = 4L,
  cores = 4L,
  seed = 1L,
  silent = 2L,
  refresh = 0L
)
```

## Arguments

- region_mat:

  Region-by-sample numeric matrix.

- coldata:

  Sample metadata.

- design:

  One-sided formula.

- contrast:

  Character coefficient name.

- iter, chains, cores, seed:

  MCMC settings passed to
  [`brms::brm`](https://paulbuerkner.com/brms/reference/brm.html).

- silent:

  [`brms::brm`](https://paulbuerkner.com/brms/reference/brm.html)'s
  `silent` level (default 2).

- refresh:

  [`brms::brm`](https://paulbuerkner.com/brms/reference/brm.html)'s
  `refresh` (default 0, no progress).

## Value

A list with the same shape as `fit_bread_summary`, but each per-region
fit carries `$draws` (a numeric vector of posterior draws for the
contrast coefficient) instead of `mu_n / Lambda_n_inv / a_n / b_n`.
`posterior_summary` and `posterior_draws` detect the shape
automatically.

## Priors

v1 uses brms's default priors (weakly informative). The BREAD
[`bread_prior()`](https://baczemin.github.io/BREAD/reference/bread_prior.md)
object is ignored by this backend — it only drives the conjugate
backend. Custom brms priors via the `brms_prior` argument are planned
for a future release.
