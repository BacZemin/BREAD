# Feature-set classification summary

Bar chart of classification counts across all fitted regions. When
`feature_class_col` is supplied, bars are stacked by feature class so
users can see, e.g., how PRC / CGI / LAD subsets partition into hyper
vs. hypo.

## Usage

``` r
plot_feature_set(fit, feature_class_col = NULL)
```

## Arguments

- fit:

  A [BreadFit](https://baczemin.github.io/BREAD/reference/BreadFit.md).

- feature_class_col:

  Optional column name in `fit@mapping` to stratify by.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.
