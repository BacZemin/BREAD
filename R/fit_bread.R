#' Fit a Bayesian region-specific methylation model
#'
#' Main user-facing entry point. Given a [SummarizedExperiment::SummarizedExperiment]
#' with a methylation assay and a [GenomicRanges::GRanges] of predefined regions,
#' BREAD maps probes to regions, summarizes them per sample, and fits Bayesian
#' region-level models to produce posterior probabilities of directional
#' methylation change under the contrast of interest.
#'
#' @param se A [SummarizedExperiment::SummarizedExperiment] with a methylation assay.
#' @param features A [GenomicRanges::GRanges] of user-defined regions.
#' @param design A one-sided formula giving the model design, e.g. `~ group + sex`.
#' @param contrast Character coefficient name of interest (v1). Contrast vectors planned.
#' @param assay_name Assay name in `se`. Default `"M"`.
#' @param input_scale `"M"` or `"Beta"`. Beta inputs are converted to M internally.
#' @param mode `"summary"` (v1) or `"hierarchical"` (planned).
#' @param summary_fun Region summary: `"mean"` (default), `"median"`, `"weighted_mean"`, `"pc1"`.
#' @param delta Effect-size threshold on the M-value scale. Default 0.10.
#' @param prob_cutoff Posterior probability cutoff for classification. Default 0.95.
#' @param min_probes Minimum probes per region. Default 3.
#' @param feature_class_col Column in `mcols(features)` giving feature class (future).
#' @param backend Sampler backend: `"brms"` (default) or `"cmdstanr"`.
#' @param iter,chains,cores,seed MCMC settings.
#'
#' @return A [BreadFit] object.
#' @export
fit_bread <- function(se,
                      features,
                      design,
                      contrast = NULL,
                      assay_name = "M",
                      input_scale = c("M", "Beta"),
                      mode = c("summary", "hierarchical"),
                      summary_fun = c("mean", "median", "weighted_mean", "pc1"),
                      delta = 0.10,
                      prob_cutoff = 0.95,
                      min_probes = 3L,
                      feature_class_col = NULL,
                      backend = c("brms", "cmdstanr"),
                      iter = 2000L,
                      chains = 4L,
                      cores = 4L,
                      seed = 1L) {
  input_scale <- match.arg(input_scale)
  mode        <- match.arg(mode)
  summary_fun <- match.arg(summary_fun)
  backend     <- match.arg(backend)

  stop("fit_bread() not yet implemented \u2014 scaffold only.")
}
