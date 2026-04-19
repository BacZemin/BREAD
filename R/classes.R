#' S4 classes for BREAD
#'
#' `BreadFit` holds a fitted Bayesian region-level methylation model together
#' with mapping metadata, posterior summaries, and diagnostics. `BreadResults`
#' is the tidy region-level result table returned by [results()].
#'
#' Both classes are intentionally thin S4 wrappers so the internals remain
#' inspectable and the slots are stable entry points for downstream packages.
#'
#' @name BREAD-classes
#' @keywords internal
NULL

#' @rdname BREAD-classes
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

#' @rdname BREAD-classes
#' @exportClass BreadResults
setClass(
  "BreadResults",
  representation(
    table  = "ANY",
    params = "list"
  ),
  prototype(params = list())
)
