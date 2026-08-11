#' BREAD: Bayesian Region-specific DNA Methylation Inference
#'
#' Targeted Bayesian inference for predefined DNA methylation regions on array
#' data, supplied either as a [SummarizedExperiment::SummarizedExperiment] or
#' as a probe-by-sample matrix. For each user-supplied region, BREAD fits a
#' Bayesian model, computes posterior probabilities of directional methylation
#' change, and classifies regions as hypermethylated, hypomethylated,
#' unchanged, or inconclusive at user-configurable effect-size and probability
#' thresholds. The `unchanged` class reports regions whose posterior lies
#' inside the region of practical equivalence: positive evidence of no change,
#' as distinct from insufficient evidence either way.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
