# Fit a Bayesian region-specific methylation model

Main user-facing entry point. Given a
[SummarizedExperiment::SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
with a methylation assay and a
[GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
of predefined regions, BREAD maps probes to regions, summarizes them per
sample, and fits Bayesian region-level models to produce posterior
probabilities of directional methylation change under the contrast of
interest. Regions are classified as hypermethylated, hypomethylated, or
inconclusive.

## Usage

``` r
fit_bread(
  se,
  features,
  design,
  contrast = NULL,
  assay_name = "M",
  input_scale = c("M", "Beta"),
  mode = c("summary", "hierarchical"),
  summary_fun = c("mean", "median", "weighted_mean", "pc1"),
  delta = 0.1,
  prob_cutoff = 0.95,
  min_probes = 3L,
  feature_class_col = NULL,
  backend = c("conjugate", "brms", "cmdstanr"),
  prior = NULL,
  iter = 2000L,
  chains = 4L,
  cores = 4L,
  seed = 1L
)
```

## Arguments

- se:

  A
  [SummarizedExperiment::SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  with a methylation assay.

- features:

  A
  [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  of user-defined regions.

- design:

  A one-sided formula giving the model design, e.g. `~ group + sex`.

- contrast:

  Character coefficient name of interest. If `NULL` (default), the first
  non-intercept coefficient is used and a message is emitted.

- assay_name:

  Assay name in `se`. Default `"M"`.

- input_scale:

  `"M"` or `"Beta"`. Beta inputs are converted to M internally.

- mode:

  `"summary"` (v1) or `"hierarchical"` (planned).

- summary_fun:

  Region summary: `"mean"` (default), `"median"`, `"weighted_mean"`,
  `"pc1"`.

- delta:

  Effect-size threshold on the M-value scale. Default `0.10`.

- prob_cutoff:

  Posterior probability cutoff for classification. Default `0.95`.

- min_probes:

  Minimum probes per region. Default `3`.

- feature_class_col:

  Column in `mcols(features)` giving feature class (reserved for
  class-level pooling in M3).

- backend:

  One of `"conjugate"` (default, v1), `"brms"`, `"cmdstanr"`.

- prior:

  A
  [`bread_prior()`](https://baczemin.github.io/BREAD/reference/bread_prior.md)
  object or `NULL` for defaults.

- iter, chains, cores, seed:

  MCMC settings. Reserved for `"brms"` / `"cmdstanr"` backends; ignored
  under `"conjugate"`.

## Value

A [BreadFit](https://baczemin.github.io/BREAD/reference/BreadFit.md)
object. Use
[`results()`](https://baczemin.github.io/BREAD/reference/results.md),
[`classifications()`](https://baczemin.github.io/BREAD/reference/classifications.md),
or
[`posterior_draws()`](https://baczemin.github.io/BREAD/reference/posterior_draws.md)
to inspect.

## v1 backend

`fit_bread()` uses an exact conjugate Normal-Inverse-Gamma posterior
(`backend = "conjugate"`). `"brms"` and `"cmdstanr"` are planned for v2
(needed for partial pooling, non-conjugate priors, hierarchical mode)
and currently error. The MCMC arguments `iter`, `chains`, `cores`,
`seed` are reserved for those backends and are no-ops under
`"conjugate"`.

## See also

[`validate_bread_input()`](https://baczemin.github.io/BREAD/reference/validate_bread_input.md),
[`map_probes_to_features()`](https://baczemin.github.io/BREAD/reference/map_probes_to_features.md),
[`summarize_features()`](https://baczemin.github.io/BREAD/reference/summarize_features.md),
[`brms::posterior_summary()`](https://paulbuerkner.com/brms/reference/posterior_summary.html),
[`classify_regions()`](https://baczemin.github.io/BREAD/reference/classify_regions.md)
