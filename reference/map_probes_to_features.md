# Map probes to user-defined regions

Computes probe-to-region overlaps, preserves feature metadata, and drops
regions with fewer than `min_probes` overlapping probes.

## Usage

``` r
map_probes_to_features(se, features, min_probes = 3L)
```

## Arguments

- se:

  A
  [SummarizedExperiment::SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  with non-empty `rowRanges()`.

- features:

  A
  [GenomicRanges::GRanges](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
  of regions. If `names(features)` is `NULL` or empty, IDs
  `region_1, region_2, ...` are generated.

- min_probes:

  Integer. Regions with fewer overlapping probes are dropped. Default
  `3L`.

## Value

A `data.frame` with (at minimum) columns `probe_id`, `probe_idx`,
`region_id`, `region_idx`, `n_probes`, plus any `mcols(features)`
columns. Attributes:

- `dropped_regions` : character vector of region IDs excluded.

- `min_probes` : the threshold applied.

- `n_features_in` : regions supplied.

- `n_features_out` : regions retained.

## Details

Probes that overlap multiple regions are emitted once per region (long
format). Probes with no region are silently excluded from the returned
mapping. Regions excluded by `min_probes` (including those with zero
overlaps) are recorded on `attr(mapping, "dropped_regions")`.
