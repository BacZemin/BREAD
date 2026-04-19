#' Methods for [BreadFit] and [BreadResults]
#'
#' @name BREAD-methods
#' @keywords internal
NULL

#' Extract the region-level results table from a [BreadFit]
#'
#' @param object A [BreadFit].
#' @param ... Unused.
#' @return A data frame with one row per region.
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
#' @return A named character vector of classifications per region
#'   (names are region IDs, values are one of `"hypermethylated"`,
#'   `"hypomethylated"`, `"inconclusive"`).
#' @export
setGeneric("classifications",
           function(object, ...) standardGeneric("classifications"))

#' @rdname classifications
setMethod("classifications", "BreadFit", function(object, ...) {
  res <- object@results
  if (is.null(res) || !"classification" %in% colnames(res)) return(NULL)
  out <- as.character(res$classification)
  names(out) <- as.character(res$region_id)
  out
})

#' Draw from the posterior of a region's contrast coefficient
#'
#' Samples are drawn from the marginal scaled Student-t posterior of the
#' contrast coefficient, `beta ~ mu_n + scale * t_{2 a_n}`.
#'
#' @param object A [BreadFit].
#' @param region_id Character vector of region IDs (all regions if `NULL`).
#' @param n Number of draws per region. Default `4000`.
#' @param seed Optional integer seed. If `NULL`, the global RNG state is used.
#' @param ... Unused.
#' @return A long `data.frame` with columns `region_id`, `draw`, `value`.
#'
#' @importFrom stats rt
#' @export
setGeneric("posterior_draws",
           function(object, region_id = NULL, ...) standardGeneric("posterior_draws"))

#' @rdname posterior_draws
setMethod("posterior_draws", "BreadFit",
  function(object, region_id = NULL, n = 4000L, seed = NULL, ...) {
    if (!is.null(seed)) set.seed(seed)
    fits <- object@model$fits
    if (is.null(fits))
      stop("BreadFit has no `model$fits`; was fit_bread() run successfully?",
           call. = FALSE)
    k <- object@model$contrast_idx

    target <- if (is.null(region_id)) names(fits) else region_id
    missing_rids <- setdiff(target, names(fits))
    if (length(missing_rids) > 0L) {
      stop("region_id(s) not found: ",
           paste(shQuote(missing_rids), collapse = ", "), call. = FALSE)
    }

    rows <- lapply(target, function(rid) {
      f <- fits[[rid]]
      if (!is.na(f$error)) {
        return(data.frame(region_id = rid, draw = seq_len(n),
                          value = NA_real_, stringsAsFactors = FALSE))
      }
      nu <- 2 * f$a_n
      mu <- f$mu_n[k]
      s  <- sqrt(max((f$b_n / f$a_n) * f$Lambda_n_inv[k, k],
                     .Machine$double.eps))
      data.frame(
        region_id = rid,
        draw      = seq_len(n),
        value     = mu + s * stats::rt(n, df = nu),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  }
)

#' @export
setMethod("show", "BreadFit", function(object) {
  cat("<BreadFit>\n")
  cat("  mode       :", object@mode, "\n")
  cat("  backend    :", object@diagnostics$backend %||% "unknown", "\n")
  cat("  input_scale:", object@input_scale, "\n")
  cat("  assay      :", object@assay_name, "\n")
  cat("  contrast   :",
      if (!is.null(object@params$contrast)) object@params$contrast else "unknown",
      "\n")
  cat("  delta      :", object@params$delta, "\n")
  cat("  prob_cutoff:", object@params$prob_cutoff, "\n")
  cat("  n_regions  :",
      if (length(object@features) > 0L) length(object@features) else 0L,
      "(of ",
      if (!is.null(object@diagnostics$n_features_in))
        object@diagnostics$n_features_in else NA,
      " input)\n")
  res <- object@results
  if (!is.null(res) && "classification" %in% colnames(res)) {
    tab <- table(res$classification)
    cat("  classifications:\n")
    for (nm in names(tab))
      cat("    ", format(nm, width = 17L), tab[[nm]], "\n", sep = "")
  }
  invisible(object)
})

# `%||%` — explicit because we target older R too
`%||%` <- function(a, b) if (is.null(a)) b else a
