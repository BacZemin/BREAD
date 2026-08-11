#' Fit BREAD hierarchical (CpG-within-region) model
#'
#' Planned for v2. Models CpGs nested within regions with CpG-specific
#' offsets and partial pooling of region-level effects.
#'
#' @inheritParams fit_bread_summary
#' @return Currently signals an error; planned to return a list with the
#'   same shape as [fit_bread_summary()].
#' @keywords internal
fit_bread_hierarchical <- function(...) {
  stop("hierarchical mode not implemented in v1")
}
