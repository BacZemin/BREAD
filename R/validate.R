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

  # Design variables must live in colData
  cd <- SummarizedExperiment::colData(se)
  dvars <- all.vars(design)
  missing_vars <- setdiff(dvars, colnames(cd))
  if (length(missing_vars) > 0L) {
    stop("Variables in `design` not found in `colData(se)`: ",
         paste(shQuote(missing_vars), collapse = ", "), ".",
         call. = FALSE)
  }

  # Contrast, if given, must match a design coefficient
  if (!is.null(contrast)) {
    if (!is.character(contrast) || length(contrast) != 1L || is.na(contrast)) {
      stop("`contrast` must currently be a single character coefficient name (v1).",
           call. = FALSE)
    }
    mm <- tryCatch(
      stats::model.matrix(design, data = as.data.frame(cd)),
      error = function(e) NULL
    )
    if (is.null(mm)) {
      stop("Could not build a model matrix from `design` and `colData(se)`. ",
           "Check for NAs or singularities in the design variables.",
           call. = FALSE)
    }
    if (!contrast %in% colnames(mm)) {
      stop("`contrast = \"", contrast, "\"` not found among design coefficients. ",
           "Available: ", paste(shQuote(colnames(mm)), collapse = ", "), ".",
           call. = FALSE)
    }
  }

  invisible(TRUE)
}
