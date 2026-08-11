#' Methods for [BreadFit] and [BreadResults]
#'
#' @return `show()` is called for its side effect of printing a summary to
#'   the console and returns its argument invisibly.
#'
#' @name BREAD-methods
#' @keywords internal
NULL

#' Extract the region-level results table from a [BreadFit]
#'
#' @param object A [BreadFit].
#' @param ... Unused.
#' @return A data frame with one row per region.
#' @examples
#' suppressPackageStartupMessages(library(SummarizedExperiment))
#'
#' se  <- readRDS(system.file("extdata", "vitc_ag06561.rds", package = "BREAD"))
#' reg <- readRDS(system.file("extdata", "vitc_regions.rds", package = "BREAD"))
#' se_ctrl <- se[, se$condition == "ctrl"]
#'
#' fit <- fit_bread(se_ctrl, reg, ~ passage)
#'
#' res <- results(fit)
#' head(res)
#' colnames(res)
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
#' @examples
#' suppressPackageStartupMessages(library(SummarizedExperiment))
#'
#' se  <- readRDS(system.file("extdata", "vitc_ag06561.rds", package = "BREAD"))
#' reg <- readRDS(system.file("extdata", "vitc_regions.rds", package = "BREAD"))
#' se_ctrl <- se[, se$condition == "ctrl"]
#'
#' fit <- fit_bread(se_ctrl, reg, ~ passage)
#'
#' cls <- classifications(fit)
#' head(cls)
#' table(cls)
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
#' @importFrom withr with_seed
#' @examples
#' suppressPackageStartupMessages(library(SummarizedExperiment))
#'
#' se  <- readRDS(system.file("extdata", "vitc_ag06561.rds", package = "BREAD"))
#' reg <- readRDS(system.file("extdata", "vitc_regions.rds", package = "BREAD"))
#' se_ctrl <- se[, se$condition == "ctrl"]
#'
#' fit <- fit_bread(se_ctrl, reg, ~ passage)
#'
#' # Always name the region(s) you want -- the default (NULL) draws from
#' # every region, which is `n` x n_regions rows.
#' rid <- results(fit)$region_id[1]
#' d <- posterior_draws(fit, region_id = rid, n = 500L, seed = 1L)
#' head(d)
#' quantile(d$value, c(0.025, 0.5, 0.975))
#' @export
setGeneric("posterior_draws",
           function(object, region_id = NULL, ...) standardGeneric("posterior_draws"))

#' @rdname posterior_draws
setMethod("posterior_draws", "BreadFit",
  function(object, region_id = NULL, n = 4000L, seed = NULL, ...) {
    fits <- object@model$fits
    if (is.null(fits))
      stop("BreadFit has no `model$fits`; was fit_bread() run successfully?",
           call. = FALSE)
    k <- object@model$contrast_idx

    target <- if (is.null(region_id)) names(fits) else as.character(region_id)
    unknown <- setdiff(target, names(fits))
    if (length(unknown) > 0L)
      stop("region_id(s) not found: ",
           paste(shQuote(unknown), collapse = ", "), call. = FALSE)

    draw_one <- function(rid) {
      f <- fits[[rid]]
      if (!is.na(f$error))
        return(data.frame(region_id = rid, draw = seq_len(n),
                          value = NA_real_, stringsAsFactors = FALSE))
      if (!is.null(f$draws) && length(f$draws) > 0L) {
        # Empirical: resample MCMC draws to n
        idx <- if (length(f$draws) >= n)
                 sample.int(length(f$draws), n, replace = FALSE)
               else
                 sample.int(length(f$draws), n, replace = TRUE)
        return(data.frame(region_id = rid, draw = seq_len(n),
                          value = f$draws[idx], stringsAsFactors = FALSE))
      }
      # Analytical: sample from scaled-t
      nu <- 2 * f$a_n
      mu <- f$mu_n[k]
      s  <- sqrt(max((f$b_n / f$a_n) * f$Lambda_n_inv[k, k],
                     .Machine$double.eps))
      data.frame(region_id = rid, draw = seq_len(n),
                 value = mu + s * stats::rt(n, df = nu),
                 stringsAsFactors = FALSE)
    }

    # with_seed() restores the caller's RNG state on exit; a bare
    # set.seed() would leak this reseed into the user's session.
    draw_all <- function() do.call(rbind, lapply(target, draw_one))
    if (is.null(seed)) draw_all() else withr::with_seed(seed, draw_all())
  }
)

#' @rdname BREAD-methods
#' @importFrom methods show
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
  cat("  rope_cutoff:", object@params$rope_cutoff %||% object@params$prob_cutoff, "\n")
  cat("  ci         :", object@params$ci %||% 0.95, "\n")
  # Count distinct regions, not ranges: `features` may hold many ranges per
  # region_id, which previously printed as e.g. "788 (of 790 input)" for what
  # was really 30 regions built from 790 probes.
  cat("  n_regions  :",
      object@diagnostics$n_features_out %||% 0L,
      "(of ",
      object@diagnostics$n_features_in %||% NA,
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

#' @rdname BREAD-methods
#' @export
setMethod("show", "BreadResults", function(object) {
  tab <- object@table
  cat("<BreadResults>\n")
  cat("  n_regions  :", if (is.null(tab)) 0L else nrow(tab), "\n")
  cat("  delta      :", object@params$delta %||% NA, "\n")
  cat("  prob_cutoff:", object@params$prob_cutoff %||% NA, "\n")
  if (!is.null(tab) && "classification" %in% colnames(tab)) {
    cat("  classifications:\n")
    counts <- table(tab$classification)
    for (nm in names(counts))
      cat("    ", format(nm, width = 17L), counts[[nm]], "\n", sep = "")
  }
  invisible(object)
})

# `%||%` — explicit because we target older R too
`%||%` <- function(a, b) if (is.null(a)) b else a
