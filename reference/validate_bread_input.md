# Validate inputs to [`fit_bread()`](https://baczemin.github.io/BREAD/reference/fit_bread.md)

Cheap structural checks run before any modeling. Verifies that `se` is a
[SummarizedExperiment::SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
with row-level genomic coordinates, that `features` is a non-empty
[GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html),
that `design` is a one-sided formula whose variables all live in
`colData(se)`, that `assay_name` is present, and that `contrast` (when
non-NULL) resolves to a coefficient of the design's model matrix.

## Usage

``` r
validate_bread_input(
  se,
  features,
  design,
  contrast = NULL,
  assay_name = "M",
  input_scale = c("M", "Beta")
)
```

## Arguments

- se:

  A
  [SummarizedExperiment::SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html).

- features:

  A
  [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  of regions.

- design:

  A one-sided `formula`, e.g. `~ group + sex`.

- contrast:

  Character coefficient name, or `NULL` to defer.

- assay_name:

  Character scalar naming an assay in `se`.

- input_scale:

  `"M"` or `"Beta"`.

## Value

Invisibly `TRUE`. Errors loudly on the first violation.
