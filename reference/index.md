# Package index

## Fit a BREAD model

The main entry point and the two backend fit helpers. Most users only
call
[`fit_bread()`](https://baczemin.github.io/BREAD/reference/fit_bread.md);
backend-specific functions are exposed for power users.

- [`fit_bread()`](https://baczemin.github.io/BREAD/reference/fit_bread.md)
  : Fit a Bayesian region-specific methylation model
- [`bread_prior()`](https://baczemin.github.io/BREAD/reference/bread_prior.md)
  : Construct a BREAD prior for the summary-mode model

## Pipeline building blocks

Exposed helpers for each stage of the pipeline — useful when you want to
inspect intermediate output or plug BREAD into a larger workflow.

- [`validate_bread_input()`](https://baczemin.github.io/BREAD/reference/validate_bread_input.md)
  :

  Validate inputs to
  [`fit_bread()`](https://baczemin.github.io/BREAD/reference/fit_bread.md)

- [`map_probes_to_features()`](https://baczemin.github.io/BREAD/reference/map_probes_to_features.md)
  : Map probes to user-defined regions

- [`summarize_features()`](https://baczemin.github.io/BREAD/reference/summarize_features.md)
  : Summarize probe-level methylation into region-level values

- [`posterior_summary()`](https://baczemin.github.io/BREAD/reference/posterior_summary.md)
  : Extract per-region posterior summaries

- [`classify_regions()`](https://baczemin.github.io/BREAD/reference/classify_regions.md)
  : Classify regions as hyper / hypo / inconclusive

## Plotting and palette

- [`plot_region_posterior()`](https://baczemin.github.io/BREAD/reference/plot_region_posterior.md)
  : Posterior density of the contrast coefficient per region
- [`plot_region_data()`](https://baczemin.github.io/BREAD/reference/plot_region_data.md)
  : Raw region-level values by contrast group
- [`plot_feature_set()`](https://baczemin.github.io/BREAD/reference/plot_feature_set.md)
  : Feature-set classification summary
- [`bread_colors()`](https://baczemin.github.io/BREAD/reference/bread_colors.md)
  : BREAD color palettes (MetBrewer "Cross")

## Classes and methods

- [`BreadFit`](https://baczemin.github.io/BREAD/reference/BreadFit.md)
  [`BreadFit-class`](https://baczemin.github.io/BREAD/reference/BreadFit.md)
  :

  The `BreadFit` S4 class

- [`BreadResults`](https://baczemin.github.io/BREAD/reference/BreadResults.md)
  [`BreadResults-class`](https://baczemin.github.io/BREAD/reference/BreadResults.md)
  :

  The `BreadResults` S4 class

- [`results()`](https://baczemin.github.io/BREAD/reference/results.md) :
  Extract the region-level results table from a BreadFit

- [`classifications()`](https://baczemin.github.io/BREAD/reference/classifications.md)
  : Extract region classifications from a BreadFit

- [`posterior_draws()`](https://baczemin.github.io/BREAD/reference/posterior_draws.md)
  : Draw from the posterior of a region's contrast coefficient
