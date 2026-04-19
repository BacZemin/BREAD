# Summarize probe-level methylation into region-level values

Collapses probe × sample methylation values into a region × sample
matrix, one row per region in the mapping. BREAD models this matrix
directly.

## Usage

``` r
summarize_features(
  se,
  mapping,
  summary_fun = c("mean", "median", "weighted_mean", "pc1"),
  input_scale = c("M", "Beta"),
  assay_name = "M"
)
```

## Arguments

- se:

  A
  [SummarizedExperiment::SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html).

- mapping:

  Data frame from
  [`map_probes_to_features()`](https://baczemin.github.io/BREAD/reference/map_probes_to_features.md).

- summary_fun:

  One of `"mean"`, `"median"`, `"weighted_mean"`, `"pc1"`.

- input_scale:

  `"M"` or `"Beta"`. `"Beta"` inputs are converted to M-values before
  summarization.

- assay_name:

  Assay name in `se`. Default `"M"`.

## Value

A numeric matrix with one row per region (in the order they first appear
in `mapping`) and one column per sample (matching
`colnames(assay(se, assay_name))`). Attributes: `summary_fun`,
`input_scale`, `assay_name`.

## Summary functions

- `"mean"` (default): arithmetic mean of probes per region per sample.

- `"median"`: per-sample median across probes in the region.

- `"weighted_mean"`: inverse-variance weighting, where each probe's
  weight is `1 / max(var_across_samples, 1e-6)`. Probes with zero
  variance receive the minimum-variance weight; if all probes have zero
  variance, weights fall back to uniform (reducing to the plain mean).

- `"pc1"`: first principal component of the probes (SVD after
  row-centering). Returned scores are sign-aligned so they correlate
  positively with the per-sample mean across probes. Scale is abstract —
  `delta` loses its M-value interpretation under `"pc1"`.
