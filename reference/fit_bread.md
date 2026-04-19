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
  delta = 0.1,
  prob_cutoff = 0.95,
  min_probes = 3L,
  feature_class_col = NULL,
  summary_fun = c("mean", "median", "weighted_mean", "pc1"),
  assay_name = NULL,
  input_scale = NULL,
  backend = c("conjugate", "brms"),
  prior = NULL,
  ...
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

- delta:

  Effect-size threshold on the M-value scale. Default `0.10`.

- prob_cutoff:

  Posterior probability cutoff for classification. Default `0.95`.

- min_probes:

  Minimum probes per region. Default `3`.

- feature_class_col:

  Column in `mcols(features)` giving feature class (used by
  [`plot_feature_set()`](https://baczemin.github.io/BREAD/reference/plot_feature_set.md);
  reserved for class-level pooling in M3).

- summary_fun:

  Region summary: `"mean"` (default), `"median"`, `"weighted_mean"`, or
  `"pc1"`.

- assay_name:

  Assay name in `se`. `NULL` auto-detects.

- input_scale:

  `"M"` or `"Beta"`. `NULL` auto-detects from value range.

- backend:

  One of `"conjugate"` (default) or `"brms"`.

- prior:

  Optional
  [`bread_prior()`](https://baczemin.github.io/BREAD/reference/bread_prior.md)
  object (conjugate backend only).

- ...:

  Additional arguments forwarded to the backend. For `backend = "brms"`,
  this accepts `iter`, `chains`, `cores`, `seed`, etc.

## Value

A [BreadFit](https://baczemin.github.io/BREAD/reference/BreadFit.md)
object. Use
[`results()`](https://baczemin.github.io/BREAD/reference/results.md),
[`classifications()`](https://baczemin.github.io/BREAD/reference/classifications.md),
or
[`posterior_draws()`](https://baczemin.github.io/BREAD/reference/posterior_draws.md)
to inspect.

## Minimal call

The typical call is:

    fit_bread(se, features, ~ condition)

Everything else has a sensible default. In particular:

- `contrast` defaults to the first non-intercept coefficient,

- `assay_name` auto-detects from the first assay that looks like
  methylation (prefers `"M"`, `"betas"`, `"Beta"`, `"beta"`),

- `input_scale` auto-detects from the assay value range (`[0,1]` →
  `"Beta"`, otherwise `"M"`).

## Backends

Default is an exact conjugate Normal-Inverse-Gamma posterior (fast, no
MCMC). `backend = "brms"` routes to
[`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html) (one
compile + per-region updates); MCMC controls `iter`, `chains`, `cores`,
`seed` can be passed through `...` to
[`fit_bread_brms()`](https://baczemin.github.io/BREAD/reference/fit_bread_brms.md).

## See also

[`validate_bread_input()`](https://baczemin.github.io/BREAD/reference/validate_bread_input.md),
[`map_probes_to_features()`](https://baczemin.github.io/BREAD/reference/map_probes_to_features.md),
[`summarize_features()`](https://baczemin.github.io/BREAD/reference/summarize_features.md),
[`posterior_summary()`](https://baczemin.github.io/BREAD/reference/posterior_summary.md),
[`classify_regions()`](https://baczemin.github.io/BREAD/reference/classify_regions.md)
