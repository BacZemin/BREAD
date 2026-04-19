# Draw from the posterior of a region's contrast coefficient

Samples are drawn from the marginal scaled Student-t posterior of the
contrast coefficient, `beta ~ mu_n + scale * t_{2 a_n}`.

## Usage

``` r
posterior_draws(object, region_id = NULL, ...)

# S4 method for class 'BreadFit'
posterior_draws(object, region_id = NULL, n = 4000L, seed = NULL, ...)
```

## Arguments

- object:

  A [BreadFit](https://baczemin.github.io/BREAD/reference/BreadFit.md).

- region_id:

  Character vector of region IDs (all regions if `NULL`).

- ...:

  Unused.

- n:

  Number of draws per region. Default `4000`.

- seed:

  Optional integer seed. If `NULL`, the global RNG state is used.

## Value

A long `data.frame` with columns `region_id`, `draw`, `value`.
