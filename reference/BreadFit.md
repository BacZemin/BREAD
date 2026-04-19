# The `BreadFit` S4 class

Main fitted object returned by
[`fit_bread()`](https://baczemin.github.io/BREAD/reference/fit_bread.md).
Thin S4 wrapper with stable slots that downstream packages can depend
on.

## Slots

- `call`:

  The original
  [`fit_bread()`](https://baczemin.github.io/BREAD/reference/fit_bread.md)
  call.

- `params`:

  List of parameters used (delta, prob_cutoff, summary_fun, mode,
  backend, contrast, min_probes, feature_class_col, iter, chains, cores,
  seed).

- `mode`:

  `"summary"` or `"hierarchical"`.

- `assay_name`:

  Assay name used from `se`.

- `input_scale`:

  `"M"` or `"Beta"`.

- `mapping`:

  Probe-to-region data frame from
  [`map_probes_to_features()`](https://baczemin.github.io/BREAD/reference/map_probes_to_features.md).

- `features`:

  `GRanges` of regions that survived `min_probes` filtering.

- `model`:

  Internal fit object from
  [`fit_bread_summary()`](https://baczemin.github.io/BREAD/reference/fit_bread_summary.md).

- `posterior`:

  Per-region posterior summary data frame.

- `results`:

  Per-region data frame with classification column.

- `diagnostics`:

  List with backend, seed, feature counts, failure counts.
