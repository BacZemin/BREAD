#' Methods for [BreadFit] and [BreadResults]
#'
#' @name BREAD-methods
#' @keywords internal
NULL

#' Extract the tidy results table from a [BreadFit]
#'
#' @param object A [BreadFit].
#' @param ... Unused.
#' @return A data frame (or [BreadResults]) with one row per region.
#' @export
setGeneric("results", function(object, ...) standardGeneric("results"))

#' @rdname results
setMethod("results", "BreadFit", function(object, ...) {
  object@results
})

#' Extract region classifications from a [BreadFit]
#'
#' @param object A [BreadFit].
#' @param ... Unused.
#' @return A character vector of classifications per region.
#' @export
setGeneric("classifications", function(object, ...) standardGeneric("classifications"))

#' @rdname classifications
setMethod("classifications", "BreadFit", function(object, ...) {
  stop("classifications() not yet implemented")
})

#' Extract posterior draws for the region-level effect
#'
#' @param object A [BreadFit].
#' @param region_id Optional region ID to subset.
#' @param ... Unused.
#' @return A tidy data frame of draws.
#' @export
setGeneric("posterior_draws", function(object, region_id = NULL, ...) standardGeneric("posterior_draws"))

#' @rdname posterior_draws
setMethod("posterior_draws", "BreadFit", function(object, region_id = NULL, ...) {
  stop("posterior_draws() not yet implemented")
})

#' @export
setMethod("show", "BreadFit", function(object) {
  cat("<BreadFit>\n")
  cat("  mode       :", object@mode, "\n")
  cat("  input_scale:", object@input_scale, "\n")
  cat("  assay      :", object@assay_name, "\n")
  cat("  n_regions  :", if (!is.null(object@features)) length(object@features) else NA, "\n")
  invisible(object)
})
