# Getting help with BREAD

## Usage questions

Ask on the **[Bioconductor support site](https://support.bioconductor.org)**
and tag your post `BREAD`. Questions there are seen by the wider Bioconductor
community and stay searchable for the next person with the same question.

Good things to include:

- what you are trying to test (the contrast, the regions)
- `assayNames(se)`, `dim(se)`, and `colData(se)`
- the `features` object (`length()`, `mcols()` column names)
- your exact `fit_bread()` call
- `sessionInfo()`

## Bugs and feature requests

Open an issue at <https://github.com/BacZemin/BREAD/issues> with a
[reprex](https://reprex.tidyverse.org/). The packaged example data is a good
basis for a self-contained reproduction:

```r
se  <- readRDS(system.file("extdata", "vitc_ag06561.rds", package = "BREAD"))
reg <- readRDS(system.file("extdata", "vitc_regions.rds", package = "BREAD"))
```

## Documentation

- Package website: <https://baczemin.github.io/BREAD/>
- `vignette("bread-intro", package = "BREAD")` — getting started
- `vignette("bread-vitc", package = "BREAD")` — real-data walkthrough
