#' Validate inputs to [fit_bread()]
#'
#' Performs cheap structural checks before any modeling: assay presence,
#' feature GRanges validity, design formula variables, contrast resolution,
#' and scale consistency. Fails early with actionable messages.
#'
#' @param se A [SummarizedExperiment::SummarizedExperiment].
#' @param features A [GenomicRanges::GRanges].
#' @param design A one-sided formula.
#' @param contrast Character coefficient name or `NULL`.
#' @param assay_name Character scalar.
#' @param input_scale `"M"` or `"Beta"`.
#' @return Invisibly `TRUE` on success; errors otherwise.
#' @keywords internal
validate_bread_input <- function(se, features, design, contrast, assay_name, input_scale) {
  stop("validate_bread_input() not yet implemented")
}
