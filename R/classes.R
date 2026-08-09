#' The `BreadFit` S4 class
#'
#' Main fitted object returned by [fit_bread()]. Thin S4 wrapper with stable
#' slots that downstream packages can depend on.
#'
#' @slot call The original `fit_bread()` call.
#' @slot params List of parameters used (contrast, delta, prob_cutoff,
#'   rope_cutoff, ci, ref_beta, summary_fun, backend, min_probes,
#'   feature_class_col).
#' @slot mode `"summary"` or `"hierarchical"`.
#' @slot assay_name Assay name used from `se`.
#' @slot input_scale `"M"` or `"Beta"`.
#' @slot mapping Probe-to-region data frame from [map_probes_to_features()].
#' @slot features `GRanges` of the ranges belonging to regions that survived
#'   `min_probes` filtering. When several ranges share a `region_id` this is
#'   longer than `nrow(results(fit))`; `diagnostics$n_features_out` is the
#'   region count.
#' @slot model Internal fit object from [fit_bread_summary()] or
#'   [fit_bread_brms()]. Both backends return the same named list:
#'   `fits`, `design_matrix`, `coef_names`, `contrast`, `contrast_idx`,
#'   `region_ids`, `prior`, `region_mat`, `design`, `coldata`. This shape is
#'   relied on by [refit_bread()], [posterior_summary()] and
#'   [plot_region_data()]; treat it as part of the interface.
#' @slot posterior Per-region posterior summary data frame.
#' @slot results Per-region data frame with classification column.
#' @slot diagnostics List with backend, feature counts, dropped regions,
#'   failure counts, timestamp, and `refit_of` when produced by
#'   [refit_bread()].
#'
#' @name BreadFit
#' @aliases BreadFit-class
#' @examples
#' suppressPackageStartupMessages(library(SummarizedExperiment))
#'
#' se  <- readRDS(system.file("extdata", "vitc_ag06561.rds", package = "BREAD"))
#' reg <- readRDS(system.file("extdata", "vitc_regions.rds", package = "BREAD"))
#' se_ctrl <- se[, se$condition == "ctrl"]
#'
#' fit <- fit_bread(se_ctrl, reg, ~ passage)
#' fit
#'
#' slotNames(fit)
#' methods::slot(fit, "diagnostics")
#' @exportClass BreadFit
setClass(
  "BreadFit",
  representation(
    call        = "call",
    params      = "list",
    mode        = "character",
    assay_name  = "character",
    input_scale = "character",
    mapping     = "ANY",
    features    = "ANY",
    model       = "ANY",
    posterior   = "ANY",
    results     = "ANY",
    diagnostics = "list"
  ),
  prototype(
    params      = list(),
    mode        = NA_character_,
    assay_name  = NA_character_,
    input_scale = NA_character_,
    diagnostics = list()
  )
)

#' The `BreadResults` S4 class
#'
#' Structured result table wrapper. Reserved for downstream reporting helpers
#' that may want a dedicated class rather than a bare data.frame.
#'
#' @slot table A `data.frame` of per-region results.
#' @slot params List of classification parameters.
#'
#' @name BreadResults
#' @aliases BreadResults-class
#' @exportClass BreadResults
setClass(
  "BreadResults",
  representation(
    table  = "ANY",
    params = "list"
  ),
  prototype(params = list())
)

#' @rdname BreadResults
#'
#' @param fit A [BreadFit], as returned by [fit_bread()].
#'
#' @return `BreadResults()` returns a `BreadResults` object wrapping the
#'   region-level results table and the classification parameters used.
#'
#' @examples
#' suppressPackageStartupMessages(library(SummarizedExperiment))
#'
#' se  <- readRDS(system.file("extdata", "vitc_ag06561.rds", package = "BREAD"))
#' reg <- readRDS(system.file("extdata", "vitc_regions.rds", package = "BREAD"))
#' se_ctrl <- se[, se$condition == "ctrl"]
#'
#' fit <- fit_bread(se_ctrl, reg, ~ passage)
#'
#' br <- BreadResults(fit)
#' br
#' head(methods::slot(br, "table"))
#'
#' @importFrom methods is new
#' @export
BreadResults <- function(fit) {
  if (!methods::is(fit, "BreadFit")) {
    stop("`fit` must be a BreadFit object.", call. = FALSE)
  }
  methods::new(
    "BreadResults",
    table  = results(fit),
    params = fit@params
  )
}
