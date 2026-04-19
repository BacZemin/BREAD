# BREAD color palettes (MetBrewer "Cross")

Colorblind-reasonable palette derived from the MetBrewer `Cross` palette
(Blake Robert Mills). Embedded as hex values so there is no runtime
dependency on the `MetBrewer` package.

## Usage

``` r
bread_colors(which = c("classification", "group", "cross"))
```

## Arguments

- which:

  One of:

  - `"classification"` : named 3-vector for hyper/hypo/inconclusive

  - `"group"` : unnamed 2-vector for binary contrast plots

  - `"cross"` : full 9-color palette

## Value

A character vector of hex colors (named where applicable).
