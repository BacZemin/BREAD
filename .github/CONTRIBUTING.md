# Contributing to BREAD

Thanks for your interest in BREAD. This document covers how to report a
problem and how to work on the package.

## Reporting bugs and asking questions

- **Bugs and feature requests** → open a [GitHub issue](https://github.com/BacZemin/BREAD/issues).
  Please include a [reprex](https://reprex.tidyverse.org/) and the output of
  `sessionInfo()`.
- **Usage questions** → the
  [Bioconductor support site](https://support.bioconductor.org) with the
  `BREAD` tag reaches a wider audience than the issue tracker.

A good bug report for BREAD usually needs three things: the structure of your
`SummarizedExperiment` (`assayNames()`, `dim()`, `colData()`), the `features`
`GRanges` you passed, and the exact `fit_bread()` call.

## Development setup

```r
# clone, then from the package root:
install.packages(c("devtools", "roxygen2", "testthat"))
BiocManager::install(c("SummarizedExperiment", "GenomicRanges", "S4Vectors"))

devtools::load_all()
devtools::test()
```

Optional backends and vignette dependencies (`brms`, `knowYourCG`,
`sesameData`) live in `Suggests`; tests that need them skip cleanly when they
are absent.

## Conventions

- **Documentation is roxygen2.** Never edit `NAMESPACE` or anything in `man/`
  by hand — run `devtools::document()` and commit the regenerated files.
- **Every exported object needs a runnable `@examples` block.** Bioconductor
  requires this. The packaged example data
  (`system.file("extdata", "vitc_ag06561.rds", package = "BREAD")`) is small
  and fast enough that examples can fit a real model rather than fake one.
- **Tests are testthat 3rd edition.** Shared fixtures live in
  `tests/testthat/helper-toy.R`.
- **No `library()` calls inside `R/`** — use `@importFrom`.
- **Style**: tidyverse style guide. Keep lines under 80 characters and
  functions under 50 lines where practical; BiocCheck flags both.
- Prefer `vapply()` over `sapply()`, `seq_len()`/`seq_along()` over `1:n`, and
  `TRUE`/`FALSE` over `T`/`F`.

## Before opening a pull request

```r
devtools::document()
devtools::test()
```

and a full check plus BiocCheck:

```sh
Rscript tools/run_check.R
Rscript tools/run_bioccheck.R
```

Both CI workflows (`R-CMD-check` and `bioc-check`) must be green. The
`bioc-check` workflow runs on the Bioconductor devel container and is the
authoritative gate.

Branch from `main` and use a descriptive branch name. Please do not bump the
version in a PR — that is handled at release time.

## Code of Conduct

This project follows the
[Bioconductor Code of Conduct](https://bioconductor.org/about/code-of-conduct/).
By participating you agree to abide by its terms.
