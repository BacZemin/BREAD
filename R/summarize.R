#' Summarize probe-level methylation into region-level values
#'
#' Collapses probe × sample methylation values into a region × sample matrix,
#' one row per region in the mapping. BREAD models this matrix directly.
#'
#' @section Summary functions:
#' - `"mean"` (default): arithmetic mean of probes per region per sample.
#' - `"median"`: per-sample median across probes in the region.
#' - `"weighted_mean"`: inverse-variance weighting, where each probe's weight
#'   is `1 / max(var_across_samples, 1e-6)`. Probes with zero variance receive
#'   the minimum-variance weight; if all probes have zero variance, weights
#'   fall back to uniform (reducing to the plain mean).
#' - `"pc1"`: first principal component of the probes (SVD after row-centering).
#'   Returned scores are sign-aligned so they correlate positively with the
#'   per-sample mean across probes. Scale is abstract — `delta` loses its
#'   M-value interpretation under `"pc1"`.
#'
#' @param se A [SummarizedExperiment::SummarizedExperiment].
#' @param mapping Data frame from [map_probes_to_features()].
#' @param summary_fun One of `"mean"`, `"median"`, `"weighted_mean"`, `"pc1"`.
#' @param input_scale `"M"` or `"Beta"`. `"Beta"` inputs are converted to
#'   M-values before summarization.
#' @param assay_name Assay name in `se`. Default `"M"`.
#'
#' @return A numeric matrix with one row per region (in the order they first
#'   appear in `mapping`) and one column per sample (matching
#'   `colnames(assay(se, assay_name))`). Attributes: `summary_fun`,
#'   `input_scale`, `assay_name`.
#'
#' @importFrom methods is
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom stats median var cor
#' @export
summarize_features <- function(se,
                               mapping,
                               summary_fun = c("mean", "median", "weighted_mean", "pc1"),
                               input_scale = c("M", "Beta"),
                               assay_name  = "M") {
  summary_fun <- match.arg(summary_fun)
  input_scale <- match.arg(input_scale)

  if (!methods::is(se, "SummarizedExperiment")) {
    stop("`se` must be a SummarizedExperiment.", call. = FALSE)
  }
  if (!is.data.frame(mapping)) {
    stop("`mapping` must be a data.frame from map_probes_to_features().",
         call. = FALSE)
  }
  req <- c("probe_idx", "region_id")
  missing_cols <- setdiff(req, colnames(mapping))
  if (length(missing_cols) > 0L) {
    stop("`mapping` missing required columns: ",
         paste(shQuote(missing_cols), collapse = ", "), ".", call. = FALSE)
  }
  if (nrow(mapping) == 0L) {
    stop("`mapping` is empty \u2014 nothing to summarize.", call. = FALSE)
  }
  an <- SummarizedExperiment::assayNames(se)
  if (!assay_name %in% an) {
    stop("assay `", assay_name, "` not found in se. Available: ",
         paste(shQuote(an), collapse = ", "), ".", call. = FALSE)
  }

  X <- SummarizedExperiment::assay(se, assay_name)
  if (!is.matrix(X)) X <- as.matrix(X)
  if (input_scale == "Beta") X <- .beta_to_m(X)

  # Preserve the order regions first appear in mapping
  region_lvls <- unique(mapping$region_id)
  by_region   <- split(mapping$probe_idx,
                       factor(mapping$region_id, levels = region_lvls))
  # Drop duplicate probe indices within a region (safety net)
  by_region <- lapply(by_region, unique)

  summarise_one <- switch(
    summary_fun,
    mean = function(idx) {
      sub <- X[idx, , drop = FALSE]
      colMeans(sub, na.rm = TRUE)
    },
    median = function(idx) {
      sub <- X[idx, , drop = FALSE]
      apply(sub, 2L, stats::median, na.rm = TRUE)
    },
    weighted_mean = function(idx) {
      sub <- X[idx, , drop = FALSE]
      if (nrow(sub) == 1L) return(sub[1L, ])
      row_var <- apply(sub, 1L, stats::var, na.rm = TRUE)
      row_var[!is.finite(row_var)] <- 0
      w <- 1 / pmax(row_var, 1e-6)
      w[!is.finite(w)] <- 0
      if (sum(w) == 0) w <- rep(1, nrow(sub))
      # NA-aware weighted mean per sample
      sub0 <- sub; sub0[is.na(sub0)] <- 0
      wmat <- matrix(w, nrow = nrow(sub), ncol = ncol(sub)) * !is.na(sub)
      num <- colSums(sub0 * matrix(w, nrow = nrow(sub), ncol = ncol(sub)))
      den <- colSums(wmat)
      out <- num / den
      out[!is.finite(out)] <- NA_real_
      out
    },
    pc1 = function(idx) {
      sub <- X[idx, , drop = FALSE]
      if (nrow(sub) == 1L) return(sub[1L, ])
      if (anyNA(sub)) {
        rm_ <- rowMeans(sub, na.rm = TRUE)
        for (i in seq_len(nrow(sub))) {
          na_i <- is.na(sub[i, ])
          if (any(na_i)) sub[i, na_i] <- rm_[i]
        }
      }
      subc   <- sub - rowMeans(sub)
      s      <- svd(subc, nu = 0L, nv = 1L)
      scores <- as.numeric(s$d[1L] * s$v[, 1L])
      # sign-align so scores correlate positively with per-sample mean
      mean_across <- colMeans(sub)
      if (stats::var(mean_across) > 0 &&
          stats::cor(scores, mean_across) < 0) {
        scores <- -scores
      }
      scores
    }
  )

  mat <- do.call(rbind, lapply(by_region, summarise_one))
  colnames(mat) <- colnames(X)
  rownames(mat) <- names(by_region)

  attr(mat, "summary_fun") <- summary_fun
  attr(mat, "input_scale") <- input_scale
  attr(mat, "assay_name")  <- assay_name
  mat
}
