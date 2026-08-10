#' Fit a Bayesian region-specific methylation model
#'
#' Main user-facing entry point. Given a [SummarizedExperiment::SummarizedExperiment]
#' with a methylation assay and a [GenomicRanges::GRanges] of predefined
#' regions, BREAD maps probes to regions, summarizes them per sample, and
#' fits Bayesian region-level models to produce posterior probabilities of
#' directional methylation change under the contrast of interest. Regions
#' are classified as hypermethylated, hypomethylated, unchanged (posterior
#' concentrated inside the region of practical equivalence) or inconclusive.
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
#' @param x A [SummarizedExperiment::SummarizedExperiment] with a methylation
#'   assay, or a probe-by-sample `matrix` (with `colData` and either
#'   `rowRanges` or `platform`), or a `list(betas =, sampleInfo =)` as
#'   returned by `sesameData`. See [bread_se()].
#' @param features A [GenomicRanges::GRanges] of user-defined regions. Several
#'   ranges may share a name to define one region as an exact probe set.
#' @param design A one-sided formula giving the model design, e.g. `~ group + sex`.
#' @param contrast Character coefficient name of interest. If `NULL` (default),
#'   the first non-intercept coefficient is used and a message is emitted.
#' @param colData,rowRanges,platform Only for matrix input: sample metadata,
#'   probe coordinates, and the array platform for a `sesameData` manifest
#'   lookup. Passing any of them alongside a `SummarizedExperiment` is an
#'   error. See [bread_se()].
#' @param delta Effect-size threshold on the M-value scale. Default `0.10`
#'   (a beta change of roughly 0.017 at mid-methylation, less toward the
#'   extremes -- see [bread_delta_beta()]).
#' @param prob_cutoff Posterior probability cutoff for a directional
#'   (hyper/hypo) call. Default `0.95`.
#' @param rope_cutoff Posterior probability cutoff for an `unchanged`
#'   (equivalence) call. Defaults to `prob_cutoff`; see [classify_regions()]
#'   for why it is worth setting independently.
#' @param ci Credible-interval mass reported in `ci_lo`/`ci_hi`. Default
#'   `0.95`. Independent of `prob_cutoff`.
#' @param ref_beta Reference methylation level(s) anchoring the beta-scale
#'   columns. `NULL` (default) uses each region's own mean. See
#'   [posterior_summary()].
#' @param min_probes Minimum probes per region. Default `3`.
#' @param feature_class_col Column in `mcols(features)` giving feature class
#'   (used by `plot_feature_set()`; reserved for class-level pooling in M3).
#' @param summary_fun Region summary: `"mean"` (default), `"median"`,
#'   `"weighted_mean"`, or `"pc1"`.
#' @param assay_name Assay name in `se`. `NULL` auto-detects.
#' @param input_scale `"M"` or `"Beta"`. `NULL` auto-detects from value range.
#' @param backend One of `"conjugate"` (default) or `"brms"`.
#' @param prior Optional [bread_prior()] object (conjugate backend only).
#' @param df_mode Degrees-of-freedom convention for the conjugate backend:
#'   `"conjugate"` (default, \eqn{a_n = a_0 + n/2}) or `"residual"`
#'   (\eqn{a_n = a_0 + (n-p)/2}), which reproduces the classical
#'   \eqn{t_{n-p}} marginal and matches `lm()` intervals under a weak prior.
#'   The default overstates precision by a factor \eqn{\sqrt{n/(n-p)}} on the
#'   posterior scale — negligible when \eqn{p \ll n}, material for interaction
#'   designs at small \eqn{n}. Ignored by `backend = "brms"`, which samples
#'   \eqn{\sigma^2} directly. See [fit_bread_summary()].
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
fit_bread <- function(x, features, design,
                      contrast          = NULL,
                      colData           = NULL,
                      rowRanges         = NULL,
                      platform          = NULL,
                      delta             = 0.10,
                      prob_cutoff       = 0.95,
                      rope_cutoff       = prob_cutoff,
                      ci                = 0.95,
                      ref_beta          = NULL,
                      min_probes        = 3L,
                      feature_class_col = NULL,
                      summary_fun       = c("mean", "median", "weighted_mean", "pc1"),
                      assay_name        = NULL,
                      input_scale       = NULL,
                      backend           = c("conjugate", "brms"),
                      prior             = NULL,
                      df_mode           = c("conjugate", "residual"),
                      ...) {
  the_call    <- match.call()
  summary_fun <- match.arg(summary_fun)
  backend     <- match.arg(backend)
  df_mode     <- match.arg(df_mode)
  if (identical(backend, "brms") && !missing(df_mode) &&
      identical(df_mode, "residual")) {
    message("`df_mode` is ignored for backend = \"brms\": MCMC samples ",
            "sigma^2 directly, so the residual degrees of freedom are ",
            "already accounted for.")
  }
  delta_default <- missing(delta)

  # Accept a matrix (what openSesame() hands you) as readily as a
  # SummarizedExperiment. Coerce at the boundary; everything downstream sees
  # the canonical object.
  se <- .as_bread_se(x, colData = colData, rowRanges = rowRanges,
                     platform = platform, assay_name = assay_name)

  # Auto-detect assay_name / input_scale if not given
  if (is.null(assay_name))  assay_name  <- .detect_assay_name(se)
  if (is.null(input_scale)) input_scale <- .detect_input_scale(
    SummarizedExperiment::assay(se, assay_name)
  )

  # A user handing in beta values is thinking in beta, but `delta` is on the
  # M scale. Say so once rather than letting them assume 0.10 means 10 points.
  if (identical(input_scale, "Beta") && delta_default) {
    message("`delta` = 0.10 is on the M-value scale, not beta: that is a ",
            "beta change of about ", signif(bread_delta_beta(0.10), 2),
            " at beta = 0.5, and less toward the extremes. ",
            "Use `bread_delta_m()` to pick `delta` from a target beta ",
            "change; `results()` reports a per-region `delta_beta`.")
  }

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
      prior      = prior,
      df_mode    = df_mode
    ),
    brms = fit_bread_brms(
      region_mat = region_mat,
      coldata    = SummarizedExperiment::colData(se),
      design     = design,
      contrast   = contrast,
      ...
    )
  )

  post    <- posterior_summary(model, delta = delta, ci = ci,
                               ref_beta = ref_beta)
  results <- classify_regions(post, delta = delta,
                              prob_cutoff = prob_cutoff,
                              rope_cutoff = rope_cutoff)

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
      rope_cutoff       = rope_cutoff,
      ci                = ci,
      ref_beta          = ref_beta,
      summary_fun       = summary_fun,
      backend           = backend,
      df_mode           = if (identical(backend, "conjugate")) df_mode else NA_character_,
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
