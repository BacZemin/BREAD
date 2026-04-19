# Enrich BREAD classifications against KnowYourCG databases

For each region classified as `"hypermethylated"` or `"hypomethylated"`
(or any level you pick), extracts the constituent EPIC / HM450 / MM285
probes and runs
[`knowYourCG::testEnrichment()`](https://rdrr.io/pkg/knowYourCG/man/testEnrichment.html)
against one or more CpG annotation databases. This is how you get from a
region-level BREAD result to *"is the hyper set enriched for Polycomb
targets / AP-1 binding sites / PMDs / ..."*.

## Usage

``` r
bread_kycg(
  fit,
  which = c("hypermethylated", "hypomethylated"),
  databases = NULL,
  platform = c("EPIC", "EPICv2", "HM450", "MM285"),
  universe = NULL,
  alternative = "greater",
  include_genes = FALSE
)
```

## Arguments

- fit:

  A [BreadFit](https://baczemin.github.io/BREAD/reference/BreadFit.md).

- which:

  Classification level(s) to enrich. Default
  `c("hypermethylated", "hypomethylated")` runs both as separate
  queries. Use a single string to run one.

- databases:

  Character vector of KYCG database group titles. Default `NULL` selects
  defaults by `platform`.

- platform:

  One of `"EPIC"`, `"EPICv2"`, `"HM450"`, `"MM285"`. Default `"EPIC"`.
  Only used when `databases = NULL`.

- universe:

  Optional universe of probe IDs (character vector). Default `NULL` uses
  all probes in the BREAD mapping (i.e. probes covered by `features`),
  which is the right universe for a region-level enrichment.

- alternative:

  `"greater"` (default), `"two.sided"`, or `"less"`.

- include_genes:

  Passed through to
  [`knowYourCG::testEnrichment()`](https://rdrr.io/pkg/knowYourCG/man/testEnrichment.html).

## Value

A tidy `data.frame` with one row per tested set, containing a `query`
column naming the classification level (`"hypermethylated"` /
`"hypomethylated"` / ...) alongside standard KYCG columns (`dbname`,
`estimate`, `p.value`, `FDR`, ...).

## Databases

If `databases = NULL`, the function picks a sensible default set for the
platform you specify via `platform` (e.g. `"EPIC"`, `"EPICv2"`,
`"HM450"`, `"MM285"`). Pass a character vector of
[`knowYourCG::listDBGroups()`](https://rdrr.io/pkg/knowYourCG/man/listDBGroups.html)
titles to override.

## Examples

``` r
if (FALSE) { # \dontrun{
  fit <- fit_bread(se, features, ~ group)
  enr <- bread_kycg(fit, platform = "EPIC")
  head(enr[enr$FDR < 0.01, ])
} # }
```
