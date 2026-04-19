#' Classify regions as hyper / hypo / inconclusive
#'
#' Applies the BREAD decision rule per region:
#' - hypermethylated if `P(beta > delta) >= prob_cutoff`
#' - hypomethylated  if `P(beta < -delta) >= prob_cutoff`
#' - inconclusive    otherwise
#'
#' @param post A data frame of posterior summaries per region (output of
#'   [posterior_summary()]).
#' @param delta Effect-size threshold. Default 0.10.
#' @param prob_cutoff Posterior probability cutoff. Default 0.95.
#' @return The input data frame with an added `classification` column.
#' @keywords internal
classify_regions <- function(post, delta = 0.10, prob_cutoff = 0.95) {
  stop("classify_regions() not yet implemented")
}
