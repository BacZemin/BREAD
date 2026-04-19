# Raw region-level values by contrast group

Boxplot + jitter of the summarized region values for a single region,
grouped by the first variable in the design formula.

## Usage

``` r
plot_region_data(fit, region_id)
```

## Arguments

- fit:

  A [BreadFit](https://baczemin.github.io/BREAD/reference/BreadFit.md).

- region_id:

  Single region ID.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.
