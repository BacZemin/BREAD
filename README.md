
<!-- README.md is generated from README.Rmd. Please edit that file, then run
     devtools::build_readme() and commit both files. -->

# BREAD

<!-- badges: start -->

[![R-CMD-check](https://github.com/BacZemin/BREAD/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/BacZemin/BREAD/actions/workflows/R-CMD-check.yaml)
[![bioc-check](https://github.com/BacZemin/BREAD/actions/workflows/bioc-check.yaml/badge.svg)](https://github.com/BacZemin/BREAD/actions/workflows/bioc-check.yaml)
[![pkgdown](https://github.com/BacZemin/BREAD/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/BacZemin/BREAD/actions/workflows/pkgdown.yaml)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License:
MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Docs](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://baczemin.github.io/BREAD/)
<!-- badges: end -->

**BREAD** (Bayesian Region-specific DNA methylation inference) provides
targeted Bayesian inference for **predefined** DNA methylation regions
from array data stored in `SummarizedExperiment` objects. For each
region you supply, BREAD fits a Bayesian model, computes the posterior
probability of a directional methylation change, and classifies the
region as **hypermethylated**, **hypomethylated**, or **inconclusive**
at user-configurable effect-size and probability thresholds.

Unlike genome-wide DMR callers that scan for regions, BREAD answers a
different question: *given regions I already care about (PRC2 targets,
CGIs, LADs, a chromHMM state, a custom BED), what is the posterior
evidence for methylation change in each one, and how confident am I?*
Output is a per-region posterior — effect size, credible interval, and
directional probabilities — not just a p-value.

## Installation

BREAD is in development. Install the latest version from GitHub:

``` r
# install.packages("pak")
pak::pak("BacZemin/BREAD")

# or with remotes:
# install.packages("remotes")
remotes::install_github("BacZemin/BREAD")
```

BREAD depends on Bioconductor packages (`SummarizedExperiment`,
`GenomicRanges`, `S4Vectors`). The optional `brms` backend additionally
needs `brms` plus a working Stan toolchain.

## Example

BREAD ships a small packaged dataset so you can run the whole pipeline
out of the box: 8 EPICv2 arrays from a fibroblast passage-aging ×
vitamin C experiment, plus 500 predefined regions spanning five feature
classes (PMD, PRC-CGI, bivalent, …).

``` r
library(BREAD)
suppressPackageStartupMessages(library(SummarizedExperiment))

se  <- readRDS(system.file("extdata", "vitc_ag06561.rds", package = "BREAD"))
reg <- readRDS(system.file("extdata", "vitc_regions.rds", package = "BREAD"))

# Which regions change methylation with passage in the control fibroblasts?
se_ctrl <- se[, se$condition == "ctrl"]

fit <- fit_bread(
  se_ctrl,
  features          = reg,
  design            = ~ passage,
  feature_class_col = "feature_class"
)

fit                     # summary: n regions, backend, classification counts
#> <BreadFit>
#>   mode       : summary 
#>   backend    : conjugate 
#>   input_scale: Beta 
#>   assay      : betas 
#>   contrast   : passagelate 
#>   delta      : 0.1 
#>   prob_cutoff: 0.95 
#>   n_regions  : 500 (of  500  input)
#>   classifications:
#>     hypermethylated  76
#>     hypomethylated   39
#>     inconclusive     385
```

``` r
res <- results(fit)     # one row per region
table(res$classification)
#> 
#> hypermethylated  hypomethylated    inconclusive 
#>              76              39             385
```

That is the whole pattern: a `SummarizedExperiment` of array data, a
`GRanges` of regions to test, and a model formula referencing columns of
`colData(se)`. `fit_bread()` auto-detects the assay and whether values
are on the beta or M scale, so the three arguments above are usually all
you need. `feature_class_col` is optional — supply it when your regions
carry a grouping column you want summarized.

### Reading the results

`results(fit)` returns one row per region:

| column           | meaning                                               |
|------------------|-------------------------------------------------------|
| `region_id`      | region identifier                                     |
| `n`              | number of probes summarized in the region             |
| `mean_effect`    | posterior mean methylation change (M-scale)           |
| `ci_lo`, `ci_hi` | 95% credible interval                                 |
| `p_gt_delta`     | P(effect \> +delta) — evidence for hypermethylation   |
| `p_lt_neg_delta` | P(effect \< -delta) — evidence for hypomethylation    |
| `classification` | `hypermethylated` / `hypomethylated` / `inconclusive` |

``` r
head(res[, c("region_id", "n", "mean_effect", "ci_lo", "ci_hi",
             "classification")])
#>             region_id n   mean_effect       ci_lo     ci_hi  classification
#> 1         PRC_CGI_025 4  0.1031847688 -0.07029387 0.2766634    inconclusive
#> 2 Active_promoter_080 4  0.3106513471  0.15979964 0.4615031 hypermethylated
#> 3    PMD_soloWCGW_066 4 -0.2168276139 -0.64352511 0.2098699    inconclusive
#> 4         PRC_CGI_039 4 -0.0003147594 -0.40264440 0.4020149    inconclusive
#> 5    PMD_soloWCGW_030 4  0.2989332658  0.09707228 0.5007942 hypermethylated
#> 6        Bivalent_090 4  0.0603224027 -0.31407739 0.4347222    inconclusive
```

A region is called **hypermethylated** if `p_gt_delta >= prob_cutoff`,
**hypomethylated** if `p_lt_neg_delta >= prob_cutoff`, otherwise
**inconclusive**. Defaults are `delta = 0.10` (M-scale) and
`prob_cutoff = 0.95`; both are arguments to `fit_bread()`.

`classifications(fit)` returns just the per-region calls, and
`posterior_draws(fit)` gives posterior samples for downstream summaries.

### Backends

- `backend = "conjugate"` (default) — analytic Normal-Inverse-Gamma
  posterior, no MCMC. Hundreds of regions fit in well under a second.
- `backend = "brms"` — full MCMC via Stan; compiles once, then reuses
  the compiled model across regions. Use when you need the flexibility
  of a full Bayesian fit.

### KnowYourCG enrichment

`bread_kycg()` takes the probes in your hyper- or hypo-classified
regions and runs `knowYourCG::testEnrichment()` against curated CpG
databases, so you can ask what genomic features your called regions are
enriched for.

## Vignettes

Two worked examples on real data (rendered on the [documentation
site](https://baczemin.github.io/BREAD/)):

- **Getting started** (`bread-intro`) — TCGA HM450 matched normal/tumour
  pairs over chromHMM chromatin-state regions, end-to-end through
  `bread_kycg()`.
- **Vitamin C EPICv2** (`bread-vitc`) — the packaged fibroblast
  passage-aging × vitamin C experiment, recovering classic PMD-hypo /
  PRC-CGI-hyper region signatures.

## Status

Milestone 1 (MVP) is complete: the full `fit_bread()` pipeline,
conjugate and brms backends, S4 classes with accessors, plotting
helpers, KYCG integration, two real-data vignettes, and a live pkgdown
site. Partial pooling across regions, feature-class priors,
random-effects designs, and contrast vectors are on the roadmap. See
`NEWS.md` for the changelog.

## Citation

BREAD does not have an associated publication yet. For now, cite the
package itself:

``` r
citation("BREAD")
To cite package 'BREAD' in publications use:

  Park J (2026). _BREAD: Bayesian Region-specific DNA Methylation
  Inference_. R package version 0.99.0,
  https://baczemin.github.io/BREAD,
  <https://github.com/BacZemin/BREAD>.

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {BREAD: Bayesian Region-specific DNA Methylation Inference},
    author = {Jaemin Park},
    year = {2026},
    note = {R package version 0.99.0, https://baczemin.github.io/BREAD},
    url = {https://github.com/BacZemin/BREAD},
  }
```

## License

MIT © Jaemin Park. See `LICENSE`.
