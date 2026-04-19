#' Feature-set level summary report
#'
#' Aggregates region-level classifications into a feature-set / feature-class
#' summary (counts and proportions of hyper / hypo / inconclusive).
#'
#' @param fit A [BreadFit].
#' @param feature_class_col Column in `mcols(features)` defining feature class.
#' @return A data frame with one row per feature class.
#' @keywords internal
report_feature_set <- function(fit, feature_class_col = NULL) {
  stop("report_feature_set() not yet implemented")
}
