#' Summarize probe-level methylation into region-level values
#'
#' Aggregates probe-level methylation per region per sample. Supports
#' `"mean"` (default), `"median"`, `"weighted_mean"`, and `"pc1"`.
#'
#' @param se A [SummarizedExperiment::SummarizedExperiment].
#' @param mapping Output of [map_probes_to_features()].
#' @param summary_fun Summary function name.
#' @param input_scale `"M"` or `"Beta"`; beta values are converted to M internally.
#' @return A numeric matrix of regions (rows) by samples (columns).
#' @keywords internal
summarize_features <- function(se,
                               mapping,
                               summary_fun = c("mean", "median", "weighted_mean", "pc1"),
                               input_scale = c("M", "Beta")) {
  summary_fun <- match.arg(summary_fun)
  input_scale <- match.arg(input_scale)
  stop("summarize_features() not yet implemented")
}
