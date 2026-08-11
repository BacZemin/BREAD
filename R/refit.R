#' Re-fit or re-threshold an existing BreadFit
#'
#' Repeats the modelling step of [fit_bread()] on a fit you already have,
#' reusing its region-by-sample matrix. Probe-to-region mapping and region
#' summarization — by far the expensive parts — are never repeated.
#'
#' Every argument defaults to `NULL`, meaning "keep what the original fit
#' used". Supply only what changes.
#'
#' @section Why this exists:
#' Label-permutation calibration is the natural way to check a posterior at
#' small n: shuffle the group labels a few hundred times and see where the
#' observed effect falls in the resulting null. That needs the region matrix
#' computed once and only the fit repeated. Without a public entry point the
#' only route was `BREAD:::fit_bread_summary()`, which is exactly the sort of
#' thing users should not have to reach for.
#'
#' ```
#' nulls <- vapply(permutations, function(g) {
#'   cd <- coldata; cd$genotype <- g
#'   results(refit_bread(fit, colData = cd))$prob_hyper[i]
#' }, numeric(1))
#' ```
#'
#' @section Re-thresholding is free:
#' When only `delta`, `prob_cutoff`, `rope_cutoff`, `ci` or `ref_beta` change,
#' the model is not re-fitted at all — the stored posterior is re-summarized
#' and re-classified. So sweeping a delta x cutoff grid costs essentially
#' nothing.
#'
#' @param fit A [BreadFit] from [fit_bread()].
#' @param colData Replacement sample metadata, with one row per column of the
#'   region matrix. If it has rownames they are matched and reordered against
#'   the matrix columns.
#' @param design Replacement one-sided design formula.
#' @param contrast Replacement coefficient name.
#' @param delta,prob_cutoff,rope_cutoff,ci,ref_beta Replacement posterior and
#'   classification settings. See [fit_bread()].
#' @param prior Replacement [bread_prior()] (conjugate backend only).
#' @param backend Replacement backend. Note that a `"brms"` refit recompiles
#'   the Stan model; use `"conjugate"` for permutation work.
#' @param ... Passed to the brms backend when `backend = "brms"`.
#'
#' @return A new [BreadFit]. The `mapping`, `features`, `mode`, `assay_name`
#'   and `input_scale` slots are carried over unchanged; `diagnostics` gains a
#'   `refit_of` timestamp naming the parent fit.
#'
#' @seealso [fit_bread()], [posterior_summary()], [classify_regions()]
#' @importFrom methods is new
#' @examples
#' suppressPackageStartupMessages(library(SummarizedExperiment))
#'
#' se  <- readRDS(system.file("extdata", "vitc_ag06561.rds", package = "BREAD"))
#' reg <- readRDS(system.file("extdata", "vitc_regions.rds", package = "BREAD"))
#' se_ctrl <- se[, se$condition == "ctrl"]
#'
#' fit <- fit_bread(se_ctrl, reg, ~ passage)
#'
#' # Re-threshold without re-fitting anything
#' table(results(refit_bread(fit, delta = 0.25))$classification)
#'
#' # Relax only the equivalence bar
#' table(results(refit_bread(fit, rope_cutoff = 0.80))$classification)
#'
#' # Re-fit against shuffled labels (one draw from a permutation null)
#' cd <- as.data.frame(colData(se_ctrl))
#' cd$passage <- sample(cd$passage)
#' head(results(refit_bread(fit, colData = cd))$prob_hyper)
#' @export
refit_bread <- function(fit,
                        colData     = NULL,
                        design      = NULL,
                        contrast    = NULL,
                        delta       = NULL,
                        prob_cutoff = NULL,
                        rope_cutoff = NULL,
                        ci          = NULL,
                        ref_beta    = NULL,
                        prior       = NULL,
                        backend     = NULL,
                        ...) {
  the_call <- match.call()
  if (!methods::is(fit, "BreadFit")) {
    stop("`fit` must be a BreadFit, not <", class(fit)[1], ">.", call. = FALSE)
  }
  model <- fit@model
  if (is.null(model$region_mat)) {
    stop("This fit carries no region matrix, so it cannot be re-fitted. ",
         "Re-run fit_bread() from the SummarizedExperiment.", call. = FALSE)
  }

  p <- fit@params
  delta       <- delta       %||% p$delta
  prob_cutoff <- prob_cutoff %||% p$prob_cutoff
  rope_cutoff <- rope_cutoff %||% p$rope_cutoff %||% prob_cutoff
  ci          <- ci          %||% p$ci %||% 0.95
  if (is.null(ref_beta)) ref_beta <- p$ref_beta
  backend     <- backend     %||% p$backend

  # Anything that changes the likelihood forces a re-fit; anything else only
  # re-reads the posterior we already have.
  needs_fit <- !is.null(colData) || !is.null(design) ||
               !is.null(contrast) || !is.null(prior) ||
               !identical(backend, p$backend)

  if (needs_fit) {
    region_mat <- model$region_mat
    cd <- if (is.null(colData)) model$coldata else
            .align_refit_coldata(colData, colnames(region_mat))
    dsn      <- design   %||% model$design
    contrast <- contrast %||% p$contrast

    .validate_design_contrast(cd, dsn, contrast)

    # A stored prior is sized to the original design's coefficients, so a new
    # design invalidates it. Fall back to the default rather than erroring on
    # a length mismatch the user never asked about.
    inherited_prior <- if (is.null(design)) model$prior else NULL
    if (!is.null(design) && !is.null(model$prior) && is.null(prior)) {
      message("`design` changed; the stored prior no longer matches the ",
              "coefficients and the default is used instead. Pass `prior = ` ",
              "to set one for the new design.")
    }

    model <- switch(
      backend,
      conjugate = fit_bread_summary(
        region_mat = region_mat, coldata = cd, design = dsn,
        contrast = contrast, prior = prior %||% inherited_prior
      ),
      brms = fit_bread_brms(
        region_mat = region_mat, coldata = cd, design = dsn,
        contrast = contrast, ...
      ),
      stop("Unknown backend '", backend, "'.", call. = FALSE)
    )
  } else {
    contrast <- p$contrast
  }

  post <- posterior_summary(model, delta = delta, ci = ci,
                            ref_beta = ref_beta)
  res  <- classify_regions(post, delta = delta,
                           prob_cutoff = prob_cutoff,
                           rope_cutoff = rope_cutoff)

  errs <- vapply(model$fits, function(f) f$error, character(1L))
  diagnostics <- fit@diagnostics
  diagnostics$backend       <- backend
  diagnostics$n_failed_fits <- sum(!is.na(errs))
  diagnostics$timestamp     <- Sys.time()
  diagnostics$refit_of      <- fit@diagnostics$timestamp

  methods::new("BreadFit",
    call        = the_call,
    params      = list(
      contrast          = contrast,
      delta             = delta,
      prob_cutoff       = prob_cutoff,
      rope_cutoff       = rope_cutoff,
      ci                = ci,
      ref_beta          = ref_beta,
      summary_fun       = p$summary_fun,
      backend           = backend,
      min_probes        = p$min_probes,
      feature_class_col = p$feature_class_col
    ),
    mode        = fit@mode,
    assay_name  = fit@assay_name,
    input_scale = fit@input_scale,
    mapping     = fit@mapping,
    features    = fit@features,
    model       = model,
    posterior   = post,
    results     = res,
    diagnostics = diagnostics
  )
}

# Sample metadata for a refit is matched to the region matrix's columns.
# Positional assignment is the trap this exists to remove.
.align_refit_coldata <- function(colData, sample_ids) {
  if (!is.data.frame(colData) && !methods::is(colData, "DataFrame")) {
    stop("`colData` must be a data.frame or DataFrame, not <",
         class(colData)[1], ">.", call. = FALSE)
  }
  cd <- as.data.frame(colData, stringsAsFactors = FALSE)
  if (nrow(cd) != length(sample_ids)) {
    stop("`colData` has ", nrow(cd), " rows but the region matrix has ",
         length(sample_ids), " samples.", call. = FALSE)
  }
  rn <- rownames(colData)
  if (!is.null(rn) && !is.null(sample_ids)) {
    idx <- match(sample_ids, rn)
    if (anyNA(idx)) {
      stop("`colData` rownames do not cover every sample in the region ",
           "matrix. Missing: ",
           paste(shQuote(sample_ids[is.na(idx)][seq_len(min(5L, sum(is.na(idx))))]),
                 collapse = ", "), ".", call. = FALSE)
    }
    cd <- cd[idx, , drop = FALSE]
  }
  rownames(cd) <- sample_ids
  cd
}
