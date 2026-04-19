#' Fit BREAD summary-mode Bayesian model
#'
#' v1 MVP. Fits one region-level Normal regression per region via `brms`,
#' with the design's main contrast as the fixed effect of interest.
#'
#' @param region_mat Region-by-sample numeric matrix.
#' @param coldata Sample metadata ([S4Vectors::DataFrame] or data frame).
#' @param design One-sided formula defining the model.
#' @param contrast Character name of the coefficient of interest.
#' @param delta Effect-size threshold on the M-value scale.
#' @param prob_cutoff Posterior probability cutoff for classification.
#' @param iter,chains,cores,seed MCMC settings.
#' @param backend `"brms"` or `"cmdstanr"`.
#' @return A list with per-region fits and posterior summaries.
#' @keywords internal
fit_bread_summary <- function(region_mat,
                              coldata,
                              design,
                              contrast,
                              delta = 0.10,
                              prob_cutoff = 0.95,
                              iter = 2000L,
                              chains = 4L,
                              cores = 4L,
                              seed = 1L,
                              backend = c("brms", "cmdstanr")) {
  backend <- match.arg(backend)
  stop("fit_bread_summary() not yet implemented")
}
