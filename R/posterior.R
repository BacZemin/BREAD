#' Extract posterior summaries per region
#'
#' From posterior draws of each region's effect `beta_r`, computes mean,
#' median, credible interval, and posterior probabilities of positive and
#' negative effect (both signed and threshold-based on `delta`).
#'
#' @param fit Internal per-region fit object (list of `brms` fits).
#' @param delta Effect-size threshold on the M-value scale.
#' @param ci Credible-interval width. Default 0.95.
#' @return A data frame, one row per region.
#' @keywords internal
posterior_summary <- function(fit, delta = 0.10, ci = 0.95) {
  stop("posterior_summary() not yet implemented")
}
