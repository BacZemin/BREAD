#' Fit a Bayesian region-specific methylation model
#'
#' Main user-facing entry point. Given a [SummarizedExperiment::SummarizedExperiment]
#' with a methylation assay and a [GenomicRanges::GRanges] of predefined
#' regions, BREAD maps probes to regions, summarizes them per sample, and
#' fits Bayesian region-level models to produce posterior probabilities of
#' directional methylation change under the contrast of interest. Regions
#' are classified as hypermethylated, hypomethylated, or inconclusive.
#'
#' @section Minimal call:
#' The typical call is:
#' ```
#' fit_bread(se, features, ~ condition)
#' ```
#' Everything else has a sensible default. In particular:
#' - `contrast` defaults to the first non-intercept coefficient,
#' - `assay_name` auto-detects from the first assay that looks like methylation
#'   (prefers `"M"`, `"betas"`, `"Beta"`, `"beta"`),
#' - `input_scale` auto-detects from the assay value range (`[0,1]` → `"Beta"`,
#'   otherwise `"M"`).
#'
#' @section Backends:
#' Default is an exact conjugate Normal-Inverse-Gamma posterior (fast, no
#' MCMC). `backend = "brms"` routes to `brms::brm()` (one compile + per-region
#' updates); MCMC controls `iter`, `chains`, `cores`, `seed` can be passed
#' through `...` to `fit_bread_brms()`.
#'
#' @param se A [SummarizedExperiment::SummarizedExperiment] with a methylation assay.
#' @param features A [GenomicRanges::GRanges] of user-defined regions.
#' @param design A one-sided formula giving the model design, e.g. `~ group + sex`.
#' @param contrast Character coefficient name of interest. If `NULL` (default),
#'   the first non-intercept coefficient is used and a message is emitted.
#' @param delta Effect-size threshold on the M-value scale. Default `0.10`.
#' @param prob_cutoff Posterior probability cutoff for classification. Default `0.95`.
#' @param min_probes Minimum probes per region. Default `3`.
#' @param feature_class_col Column in `mcols(features)` giving feature class
#'   (used by `plot_feature_set()`; reserved for class-level pooling in M3).
#' @param summary_fun Region summary: `"mean"` (default), `"median"`,
#'   `"weighted_mean"`, or `"pc1"`.
#' @param assay_name Assay name in `se`. `NULL` auto-detects.
#' @param input_scale `"M"` or `"Beta"`. `NULL` auto-detects from value range.
#' @param backend One of `"conjugate"` (default) or `"brms"`.
#' @param prior Optional [bread_prior()] object (conjugate backend only).
#' @param ... Additional arguments forwarded to the backend. For
#'   `backend = "brms"`, this accepts `iter`, `chains`, `cores`, `seed`, etc.
#'
#' @return A [BreadFit] object. Use [results()], [classifications()], or
#'   [posterior_draws()] to inspect.
#'
#' @seealso [validate_bread_input()], [map_probes_to_features()],
#'   [summarize_features()], [posterior_summary()], [classify_regions()]
#'
#' @importFrom methods is new
#' @importFrom SummarizedExperiment colData assay assayNames
#' @importFrom S4Vectors mcols
#' @importFrom stats model.matrix
#' @examples
#' suppressPackageStartupMessages(library(SummarizedExperiment))
#'
#' se  <- readRDS(system.file("extdata", "vitc_ag06561.rds", package = "BREAD"))
#' reg <- readRDS(system.file("extdata", "vitc_regions.rds", package = "BREAD"))
#' se_ctrl <- se[, se$condition == "ctrl"]
#'
#' # Which of the 500 predefined regions change methylation with passage
#' # in the control (untreated) fibroblasts?
#' fit <- fit_bread(se_ctrl, reg, ~ passage,
#'                  feature_class_col = "feature_class")
#' fit
#'
#' head(results(fit))
#' table(results(fit)$classification)
#' @export
fit_bread <- function(se, features, design,
                      contrast          = NULL,
                      delta             = 0.10,
                      prob_cutoff       = 0.95,
                      min_probes        = 3L,
                      feature_class_col = NULL,
                      summary_fun       = c("mean", "median", "weighted_mean", "pc1"),
                      assay_name        = NULL,
                      input_scale       = NULL,
                      backend           = c("conjugate", "brms"),
                      prior             = NULL,
                      ...) {
  the_call    <- match.call()
  summary_fun <- match.arg(summary_fun)
  backend     <- match.arg(backend)

  if (!methods::is(se, "SummarizedExperiment"))
    stop("`se` must be a SummarizedExperiment, not <", class(se)[1], ">.",
         call. = FALSE)

  # Auto-detect assay_name / input_scale if not given
  if (is.null(assay_name))  assay_name  <- .detect_assay_name(se)
  if (is.null(input_scale)) input_scale <- .detect_input_scale(
    SummarizedExperiment::assay(se, assay_name)
  )

  # Validate with current (possibly NULL) contrast
  validate_bread_input(se, features, design,
                       contrast    = contrast,
                       assay_name  = assay_name,
                       input_scale = input_scale)

  # Default contrast = first non-intercept coefficient
  if (is.null(contrast)) {
    X_preview <- stats::model.matrix(
      design,
      data = as.data.frame(SummarizedExperiment::colData(se))
    )
    if (ncol(X_preview) < 2L)
      stop("Design has no non-intercept coefficient; supply `contrast` ",
           "explicitly.", call. = FALSE)
    contrast <- colnames(X_preview)[2L]
    message("`contrast` not supplied; using '", contrast,
            "' (first non-intercept coefficient).")
  }

  # feature_class_col presence check
  if (!is.null(feature_class_col)) {
    fc <- colnames(S4Vectors::mcols(features))
    if (!feature_class_col %in% fc)
      stop("`feature_class_col = \"", feature_class_col,
           "\"` not found in `mcols(features)`. Available: ",
           paste(shQuote(fc), collapse = ", "), ".", call. = FALSE)
  }

  # Pipeline
  mapping <- map_probes_to_features(se, features, min_probes = min_probes)
  if (nrow(mapping) == 0L)
    stop("No regions survived `min_probes = ", min_probes, "`.",
         call. = FALSE)

  region_mat <- summarize_features(
    se, mapping,
    summary_fun = summary_fun,
    input_scale = input_scale,
    assay_name  = assay_name
  )

  model <- switch(
    backend,
    conjugate = fit_bread_summary(
      region_mat = region_mat,
      coldata    = SummarizedExperiment::colData(se),
      design     = design,
      contrast   = contrast,
      prior      = prior
    ),
    brms = fit_bread_brms(
      region_mat = region_mat,
      coldata    = SummarizedExperiment::colData(se),
      design     = design,
      contrast   = contrast,
      ...
    )
  )

  post    <- posterior_summary(model, delta = delta, ci = 0.95)
  results <- classify_regions(post, delta = delta, prob_cutoff = prob_cutoff)

  kept_idx      <- sort(unique(mapping$region_idx))
  features_kept <- features[kept_idx]

  errs <- vapply(model$fits, function(f) f$error, character(1L))
  diagnostics <- list(
    backend         = backend,
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
      backend           = backend,
      min_probes        = min_probes,
      feature_class_col = if (is.null(feature_class_col)) NA_character_ else feature_class_col
    ),
    mode        = "summary",
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
