# Feature-set level summary report

Aggregates region-level classifications into a feature-set /
feature-class summary (counts and proportions of hyper / hypo /
inconclusive).

## Usage

``` r
report_feature_set(fit, feature_class_col = NULL)
```

## Arguments

- fit:

  A [BreadFit](https://baczemin.github.io/BREAD/reference/BreadFit.md).

- feature_class_col:

  Column in `mcols(features)` defining feature class.

## Value

A data frame with one row per feature class.
