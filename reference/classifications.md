# Extract region classifications from a [BreadFit](https://baczemin.github.io/BREAD/reference/BreadFit.md)

Extract region classifications from a
[BreadFit](https://baczemin.github.io/BREAD/reference/BreadFit.md)

## Usage

``` r
classifications(object, ...)

# S4 method for class 'BreadFit'
classifications(object, ...)
```

## Arguments

- object:

  A [BreadFit](https://baczemin.github.io/BREAD/reference/BreadFit.md).

- ...:

  Unused.

## Value

A named character vector of classifications per region (names are region
IDs, values are one of `"hypermethylated"`, `"hypomethylated"`,
`"inconclusive"`).
