# Fit BREAD summary-mode Bayesian model

v1 MVP fit. For each region, fits a Normal linear regression
`y ~ design` with a conjugate Normal-Inverse-Gamma prior, yielding an
analytical posterior — no MCMC. The marginal posterior of each
coefficient is a location-scale Student-t, used downstream by
[`posterior_summary()`](https://baczemin.github.io/BREAD/reference/posterior_summary.md).

## Usage

``` r
fit_bread_summary(region_mat, coldata, design, contrast, prior = NULL)
```

## Arguments

- region_mat:

  Region-by-sample numeric matrix (from
  [`summarize_features()`](https://baczemin.github.io/BREAD/reference/summarize_features.md)).

- coldata:

  Sample metadata
  ([S4Vectors::DataFrame](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
  or `data.frame`).

- design:

  One-sided formula.

- contrast:

  Character coefficient name of interest.

- prior:

  A
  [`bread_prior()`](https://baczemin.github.io/BREAD/reference/bread_prior.md)
  object (or `NULL` for defaults).

## Value

A list with:

- `fits`: per-region list of
  `list(mu_n, Lambda_n_inv, a_n, b_n, n, error)`

- `design_matrix`: model matrix `X` actually used

- `coef_names`: coefficient names

- `contrast`, `contrast_idx`: contrast name and its column index in `X`

- `region_ids`: rownames of `region_mat`

- `prior`: the prior applied (with `mu0`/`Lambda0` filled in)

## Model

For region r and sample i: \$\$y\_{ir} \sim \mathcal{N}(X_i \beta_r,\\
\sigma_r^2),\$\$ \$\$\beta_r \mid \sigma_r^2 \sim \mathcal{N}(\mu_0,\\
\sigma_r^2 \Lambda_0^{-1}),\$\$ \$\$\sigma_r^2 \sim
\mathrm{Inv\text{-}Gamma}(a_0, b_0).\$\$ Posterior: \$\$\Lambda_n =
X^\top X + \Lambda_0,\quad \mu_n = \Lambda_n^{-1}(X^\top y + \Lambda_0
\mu_0),\$\$ \$\$a_n = a_0 + n/2,\quad b_n = b_0 + \tfrac{1}{2}(y^\top
y + \mu_0^\top \Lambda_0 \mu_0 - \mu_n^\top \Lambda_n \mu_n).\$\$
