# Posterior density of the contrast coefficient per region

Plots the analytical scaled-Student-t posterior of the contrast
coefficient for one or more regions, with vertical guides at `0` and
`+/- delta` and a color by final classification. When multiple regions
are supplied the plot facets by region.

## Usage

``` r
plot_region_posterior(fit, region_id = NULL, n_grid = 500L, show_delta = TRUE)
```

## Arguments

- fit:

  A [BreadFit](https://baczemin.github.io/BREAD/reference/BreadFit.md).

- region_id:

  Character vector of region IDs. If `NULL` (default), all fitted
  regions are plotted.

- n_grid:

  Number of grid points for the density curve. Default `500`.

- show_delta:

  Logical. Draw dashed guides at `+/- delta`? Default `TRUE`.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.
