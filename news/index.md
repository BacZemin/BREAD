# Changelog

## BREAD 0.0.0.9000 (development)

Initial MVP — milestone 1.

### New features

- [`fit_bread()`](https://baczemin.github.io/BREAD/reference/fit_bread.md)
  — end-to-end pipeline: validate → map → summarize → fit Bayesian
  region-level model → classify. Default backend is a conjugate
  Normal-Inverse-Gamma posterior computed analytically (no MCMC).
- [`bread_prior()`](https://baczemin.github.io/BREAD/reference/bread_prior.md)
  — constructor for weakly-informative priors.
- [`validate_bread_input()`](https://baczemin.github.io/BREAD/reference/validate_bread_input.md),
  [`map_probes_to_features()`](https://baczemin.github.io/BREAD/reference/map_probes_to_features.md),
  [`summarize_features()`](https://baczemin.github.io/BREAD/reference/summarize_features.md),
  [`posterior_summary()`](https://baczemin.github.io/BREAD/reference/posterior_summary.md),
  [`classify_regions()`](https://baczemin.github.io/BREAD/reference/classify_regions.md)
  — exported pipeline helpers.
- [`plot_region_posterior()`](https://baczemin.github.io/BREAD/reference/plot_region_posterior.md),
  [`plot_region_data()`](https://baczemin.github.io/BREAD/reference/plot_region_data.md),
  [`plot_feature_set()`](https://baczemin.github.io/BREAD/reference/plot_feature_set.md)
  — publication-oriented `ggplot2` helpers for BreadFit output.
- [`bread_colors()`](https://baczemin.github.io/BREAD/reference/bread_colors.md)
  — MetBrewer `Cross` palette embedded as hex (no runtime dependency on
  the MetBrewer package).
- `BreadFit` / `BreadResults` S4 classes with accessors
  [`results()`](https://baczemin.github.io/BREAD/reference/results.md),
  [`classifications()`](https://baczemin.github.io/BREAD/reference/classifications.md),
  [`posterior_draws()`](https://baczemin.github.io/BREAD/reference/posterior_draws.md),
  and [`show()`](https://rdrr.io/r/methods/show.html).

### Known scope

- `mode = "hierarchical"` and `backend` in `c("brms", "cmdstanr")` are
  scaffolded with “planned for v2” errors.
- Hand-coded conjugate backend covers summary mode only; partial pooling
  across regions and class-level priors are deferred to M3.
