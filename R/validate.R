#' Validate inputs to [fit_bread()]
#'
#' Cheap structural checks run before any modeling. Verifies that `se` is a
#' [SummarizedExperiment::SummarizedExperiment] with row-level genomic
#' coordinates, that `features` is a non-empty [GenomicRanges::GRanges], that
#' `design` is a one-sided formula whose variables all live in `colData(se)`,
#' that `assay_name` is present, and that `contrast` (when non-NULL) resolves
#' to a coefficient of the design's model matrix.
#'
#' @param se A [SummarizedExperiment::SummarizedExperiment].
#' @param features A [GenomicRanges::GRanges] of regions.
#' @param design A one-sided `formula`, e.g. `~ group + sex`.
#' @param contrast Character coefficient name, or `NULL` to defer.
#' @param assay_name Character scalar naming an assay in `se`.
#' @param input_scale `"M"` or `"Beta"`.
#'
#' @return Invisibly `TRUE`. Errors loudly on the first violation.
#'
#' @importFrom methods is
#' @importFrom SummarizedExperiment assayNames colData rowRanges
#' @importFrom stats model.matrix
#' @examples
#' suppressPackageStartupMessages(library(SummarizedExperiment))
#'
#' se  <- readRDS(system.file("extdata", "vitc_ag06561.rds", package = "BREAD"))
#' reg <- readRDS(system.file("extdata", "vitc_regions.rds", package = "BREAD"))
#' se_ctrl <- se[, se$condition == "ctrl"]
#'
#' # The packaged data keeps beta values in an assay named "betas", so
#' # both arguments are given explicitly here. fit_bread() detects them
#' # for you; this lower-level helper does not.
#' validate_bread_input(se_ctrl, reg, ~ passage,
#'                      assay_name = "betas", input_scale = "Beta")
#' @export
validate_bread_input <- function(se,
                                 features,
                                 design,
                                 contrast = NULL,
                                 assay_name = "M",
                                 input_scale = c("M", "Beta")) {
  input_scale <- match.arg(input_scale)

  if (!methods::is(se, "SummarizedExperiment")) {
    stop("`se` must be a SummarizedExperiment, not <", class(se)[1], ">.",
         call. = FALSE)
  }
  if (!methods::is(features, "GRanges")) {
    stop("`features` must be a GRanges, not <", class(features)[1], ">.",
         call. = FALSE)
  }
  if (length(features) == 0L) {
    stop("`features` is empty \u2014 supply at least one region.", call. = FALSE)
  }
  if (!inherits(design, "formula")) {
    stop("`design` must be a formula (e.g. `~ group + sex`), not <",
         class(design)[1], ">.", call. = FALSE)
  }
  if (length(design) != 2L) {
    stop("`design` must be one-sided (e.g. `~ group`, not `y ~ group`). ",
         "BREAD supplies the response internally.", call. = FALSE)
  }

  if (!is.character(assay_name) || length(assay_name) != 1L || is.na(assay_name)) {
    stop("`assay_name` must be a single non-NA character string.", call. = FALSE)
  }
  an <- SummarizedExperiment::assayNames(se)
  if (!assay_name %in% an) {
    stop("assay `", assay_name, "` not found in se. Available: ",
         paste(shQuote(an), collapse = ", "), call. = FALSE)
  }

  # rowRanges must exist and be a GRanges
  rr <- tryCatch(
    SummarizedExperiment::rowRanges(se),
    error = function(e) NULL
  )
  if (is.null(rr) || !methods::is(rr, "GRanges") || length(rr) == 0L) {
    stop("`se` must have non-empty `rowRanges()` returning a GRanges. ",
         "Coerce row-level coordinates into rowRanges before calling BREAD.",
         call. = FALSE)
  }

  .validate_design_contrast(SummarizedExperiment::colData(se),
                            design, contrast)

  invisible(TRUE)
}

# Shared by validate_bread_input() and refit_bread(): the design variables
# must exist in the sample metadata, the model matrix must be buildable, and
# the contrast (when given) must be one of its coefficients.
#
# Also warns on a rank-deficient design. This is not cosmetic: the conjugate
# backend adds Lambda0 = 0.01 * I before the Cholesky, so a collinear design
# does not fail -- it silently returns a prior-regularised estimate that looks
# like a real fit. Permutation studies hit this routinely, and every caller
# so far has had to hand-roll its own qr() guard.
.validate_design_contrast <- function(coldata, design, contrast = NULL) {
  if (!inherits(design, "formula")) {
    stop("`design` must be a formula (e.g. `~ group + sex`), not <",
         class(design)[1], ">.", call. = FALSE)
  }
  if (length(design) != 2L) {
    stop("`design` must be one-sided (e.g. `~ group`, not `y ~ group`). ",
         "BREAD supplies the response internally.", call. = FALSE)
  }

  cd <- as.data.frame(coldata)
  missing_vars <- setdiff(all.vars(design), colnames(cd))
  if (length(missing_vars) > 0L) {
    stop("Variables in `design` not found in `colData(se)`: ",
         paste(shQuote(missing_vars), collapse = ", "), ".",
         call. = FALSE)
  }

  mm <- tryCatch(
    stats::model.matrix(design, data = cd),
    error = function(e) NULL
  )
  if (is.null(mm)) {
    stop("Could not build a model matrix from `design` and `colData(se)`. ",
         "Check for NAs or singularities in the design variables.",
         call. = FALSE)
  }

  if (!is.null(contrast)) {
    if (!is.character(contrast) || length(contrast) != 1L || is.na(contrast)) {
      stop("`contrast` must currently be a single character coefficient name (v1).",
           call. = FALSE)
    }
    if (!contrast %in% colnames(mm)) {
      stop("`contrast = \"", contrast, "\"` not found among design coefficients. ",
           "Available: ", paste(shQuote(colnames(mm)), collapse = ", "), ".",
           call. = FALSE)
    }
  }

  rk <- tryCatch(qr(mm)$rank, error = function(e) NA_integer_)
  if (!is.na(rk) && rk < ncol(mm)) {
    warning("Design matrix is rank deficient (rank ", rk, " < ", ncol(mm),
            " coefficients). The conjugate prior will absorb the deficiency ",
            "rather than error, so estimates for the collinear coefficients ",
            "are prior-driven. Columns: ",
            paste(shQuote(colnames(mm)), collapse = ", "), ".",
            call. = FALSE)
  }

  invisible(mm)
}
