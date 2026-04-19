# Extract per-region posterior summaries

Given the output of
[`fit_bread_summary()`](https://baczemin.github.io/BREAD/reference/fit_bread_summary.md)
(conjugate backend) or
[`fit_bread_brms()`](https://baczemin.github.io/BREAD/reference/fit_bread_brms.md)
(brms backend), compute the marginal posterior location / scale /
credible interval of the contrast coefficient and the posterior
probabilities used by
[`classify_regions()`](https://baczemin.github.io/BREAD/reference/classify_regions.md).

## Usage

``` r
posterior_summary(fit, delta = 0.1, ci = 0.95)
```

## Arguments

- fit:

  Output of
  [`fit_bread_summary()`](https://baczemin.github.io/BREAD/reference/fit_bread_summary.md)
  or
  [`fit_bread_brms()`](https://baczemin.github.io/BREAD/reference/fit_bread_brms.md).

- delta:

  Effect-size threshold on the M-value scale. Default `0.10`.

- ci:

  Credible-interval mass. Default `0.95`.

## Value

A `data.frame` with one row per region and columns: `region_id`, `n`,
`mean_effect`, `median_effect`, `ci_lo`, `ci_hi`, `df`, `scale`,
`p_pos`, `p_neg`, `p_gt_delta`, `p_lt_neg_delta`, `error`. `df` is
`NA_real_` for the empirical path.

## Path selection

Each per-region fit carries either

- `$mu_n`, `$Lambda_n_inv`, `$a_n`, `$b_n` (conjugate analytical path),
  in which case the marginal posterior of the contrast coefficient is a
  location-scale Student-t with `df = 2 * a_n`, or

- `$draws` (brms empirical path), in which case posterior quantities are
  computed from the MCMC draws directly.

Columns in the returned data frame are the same in both cases.
