#' Fit a Bayesian region-specific methylation model
#'
#' Main user-facing entry point. Given a [SummarizedExperiment::SummarizedExperiment]
#' with a methylation assay and a [GenomicRanges::GRanges] of predefined regions,
#' BREAD maps probes to regions, summarizes them per sample, and fits Bayesian
#' region-level models to produce posterior probabilities of directional
#' methylation change under the contrast of interest. Regions are classified
#' as hypermethylated, hypomethylated, or inconclusive.
#'
#' @section v1 backend:
#' `fit_bread()` uses an exact conjugate Normal-Inverse-Gamma posterior
#' (`backend = "conjugate"`). `"brms"` and `"cmdstanr"` are planned for v2
#' (needed for partial pooling, non-conjugate priors, hierarchical mode) and
#' currently error. The MCMC arguments `iter`, `chains`, `cores`, `seed` are
#' reserved for those backends and are no-ops under `"conjugate"`.
#'
#' @param se A [SummarizedExperiment::SummarizedExperiment] with a methylation assay.
#' @param features A [GenomicRanges::GRanges] of user-defined regions.
#' @param design A one-sided formula giving the model design, e.g. `~ group + sex`.
#' @param contrast Character coefficient name of interest. If `NULL` (default),
#'   the first non-intercept coefficient is used and a message is emitted.
#' @param assay_name Assay name in `se`. Default `"M"`.
#' @param input_scale `"M"` or `"Beta"`. Beta inputs are converted to M internally.
#' @param mode `"summary"` (v1) or `"hierarchical"` (planned).
#' @param summary_fun Region summary: `"mean"` (default), `"median"`, `"weighted_mean"`, `"pc1"`.
#' @param delta Effect-size threshold on the M-value scale. Default `0.10`.
#' @param prob_cutoff Posterior probability cutoff for classification. Default `0.95`.
#' @param min_probes Minimum probes per region. Default `3`.
#' @param feature_class_col Column in `mcols(features)` giving feature class
#'   (reserved for class-level pooling in M3).
#' @param backend One of `"conjugate"` (default, v1), `"brms"`, `"cmdstanr"`.
#' @param prior A [bread_prior()] object or `NULL` for defaults.
#' @param iter,chains,cores,seed MCMC settings. Reserved for `"brms"` /
#'   `"cmdstanr"` backends; ignored under `"conjugate"`.
#'
#' @return A [BreadFit] object. Use [results()], [classifications()], or
#'   [posterior_draws()] to inspect.
#'
#' @seealso [validate_bread_input()], [map_probes_to_features()],
#'   [summarize_features()], [posterior_summary()], [classify_regions()]
#'
#' @importFrom methods is new
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors mcols
#' @importFrom stats model.matrix
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
                      backend = c("conjugate", "brms", "cmdstanr"),
                      prior = NULL,
                      iter = 2000L,
                      chains = 4L,
                      cores = 4L,
                      seed = 1L) {
  the_call    <- match.call()
  input_scale <- match.arg(input_scale)
  mode        <- match.arg(mode)
  summary_fun <- match.arg(summary_fun)
  backend     <- match.arg(backend)

  if (mode == "hierarchical") {
    stop("`mode = \"hierarchical\"` is not yet implemented in v1. ",
         "Use mode = \"summary\".", call. = FALSE)
  }
  if (backend == "cmdstanr") {
    stop("`backend = \"cmdstanr\"` is planned; install CmdStan and use ",
         "the brms backend instead, or the default \"conjugate\".",
         call. = FALSE)
  }

  # Phase 1: validate inputs (contrast may still be NULL at this point)
  validate_bread_input(se, features, design,
                       contrast   = contrast,
                       assay_name = assay_name,
                       input_scale = input_scale)

  # Phase 2: default contrast = first non-intercept coefficient
  if (is.null(contrast)) {
    X_preview <- stats::model.matrix(
      design,
      data = as.data.frame(SummarizedExperiment::colData(se))
    )
    if (ncol(X_preview) < 2L) {
      stop("Design has no non-intercept coefficient; supply `contrast` ",
           "explicitly.", call. = FALSE)
    }
    contrast <- colnames(X_preview)[2L]
    message("`contrast` not supplied; using '", contrast,
            "' (first non-intercept coefficient).")
  }

  # Phase 3: feature_class_col check
  if (!is.null(feature_class_col)) {
    fc <- colnames(S4Vectors::mcols(features))
    if (!feature_class_col %in% fc) {
      stop("`feature_class_col = \"", feature_class_col,
           "\"` not found in `mcols(features)`. Available: ",
           paste(shQuote(fc), collapse = ", "), ".",
           call. = FALSE)
    }
  }

  # Phase 4: probe-to-region mapping
  mapping <- map_probes_to_features(se, features, min_probes = min_probes)
  if (nrow(mapping) == 0L) {
    stop("No regions survived `min_probes = ", min_probes, "`.",
         call. = FALSE)
  }

  # Phase 5: summarize probes to region x sample matrix
  region_mat <- summarize_features(
    se, mapping,
    summary_fun = summary_fun,
    input_scale = input_scale,
    assay_name  = assay_name
  )

  # Phase 6: fit conjugate posterior per region
  model <- if (backend == "brms") {
    fit_bread_brms(
      region_mat = region_mat,
      coldata    = SummarizedExperiment::colData(se),
      design     = design,
      contrast   = contrast,
      iter       = iter,
      chains     = chains,
      cores      = cores,
      seed       = seed
    )
  } else {
    fit_bread_summary(
      region_mat = region_mat,
      coldata    = SummarizedExperiment::colData(se),
      design     = design,
      contrast   = contrast,
      prior      = prior
    )
  }

  # Phase 7: posterior + classify
  post    <- posterior_summary(model, delta = delta, ci = 0.95)
  results <- classify_regions(post, delta = delta, prob_cutoff = prob_cutoff)

  # Phase 8: subset features to regions that survived min_probes
  kept_idx <- sort(unique(mapping$region_idx))
  features_kept <- features[kept_idx]

  # Diagnostics
  errs <- vapply(model$fits, function(f) f$error, character(1L))
  diagnostics <- list(
    backend         = backend,
    seed            = seed,
    n_features_in   = attr(mapping, "n_features_in"),
    n_features_out  = attr(mapping, "n_features_out"),
    dropped_regions = attr(mapping, "dropped_regions"),
    min_probes      = min_probes,
    n_failed_fits   = sum(!is.na(errs)),
    timestamp       = Sys.time()
  )

  methods::new("BreadFit",
    call        = the_call,
    params      = list(
      contrast          = contrast,
      delta             = delta,
      prob_cutoff       = prob_cutoff,
      summary_fun       = summary_fun,
      mode              = mode,
      backend           = backend,
      min_probes        = min_probes,
      feature_class_col = if (is.null(feature_class_col)) NA_character_ else feature_class_col,
      iter              = iter,
      chains            = chains,
      cores             = cores,
      seed              = seed
    ),
    mode        = mode,
    assay_name  = assay_name,
    input_scale = input_scale,
    mapping     = mapping,
    features    = features_kept,
    model       = model,
    posterior   = post,
    results     = results,
    diagnostics = diagnostics
  )
}
