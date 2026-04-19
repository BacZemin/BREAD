#' BREAD plotting helpers
#'
#' - `plot_region_posterior()` : posterior interval / density per region
#' - `plot_region_data()`      : raw region-level values by group
#' - `plot_feature_set()`      : feature-set classification summary
#'
#' All return ggplot objects and avoid any lab-wide theme dependencies.
#'
#' @name BREAD-plots
#' @keywords internal
NULL

#' @rdname BREAD-plots
#' @param fit A [BreadFit].
#' @param region_id Character region ID.
#' @export
plot_region_posterior <- function(fit, region_id) {
  stop("plot_region_posterior() not yet implemented")
}

#' @rdname BREAD-plots
#' @export
plot_region_data <- function(fit, region_id) {
  stop("plot_region_data() not yet implemented")
}

#' @rdname BREAD-plots
#' @export
plot_feature_set <- function(fit, feature_class_col = NULL) {
  stop("plot_feature_set() not yet implemented")
}
