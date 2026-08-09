#' Classify regions as hyper / hypo / inconclusive
#'
#' Applies the BREAD decision rule to the output of [posterior_summary()]:
#' - `hypermethylated` if `p_gt_delta >= prob_cutoff`
#' - `hypomethylated`  if `p_lt_neg_delta >= prob_cutoff`
#' - `inconclusive`    otherwise
#'
#' In the rare case that both probabilities exceed the cutoff (only possible
#' for very low `prob_cutoff`), the region is assigned to whichever side has
#' the larger posterior probability.
#'
#' @param post Output of [posterior_summary()].
#' @param delta Effect-size threshold used for the rule. Default `0.10`. Stored
#'   as an attribute; does not re-evaluate the posterior probabilities (those
#'   must have been computed at this same `delta` upstream).
#' @param prob_cutoff Posterior probability cutoff. Default `0.95`.
#'
#' @return The input `data.frame` with an added `classification` factor column
#'   (levels: `hypermethylated`, `hypomethylated`, `inconclusive`). Attributes
#'   `delta` and `prob_cutoff` are updated.
#'
#' @examples
#' suppressPackageStartupMessages(library(SummarizedExperiment))
#'
#' se  <- readRDS(system.file("extdata", "vitc_ag06561.rds", package = "BREAD"))
#' reg <- readRDS(system.file("extdata", "vitc_regions.rds", package = "BREAD"))
#' se_ctrl <- se[, se$condition == "ctrl"]
#'
#' fit <- fit_bread(se_ctrl, reg, ~ passage)
#' post <- posterior_summary(fit)
#'
#' cl <- classify_regions(post)
#' table(cl$classification)
#'
#' # A stricter cutoff moves borderline regions into `inconclusive`
#' table(classify_regions(post, prob_cutoff = 0.99)$classification)
#' @export
classify_regions <- function(post, delta = 0.10, prob_cutoff = 0.95) {
  if (!is.data.frame(post)) {
    stop("`post` must be a data.frame from posterior_summary().",
         call. = FALSE)
  }
  req <- c("p_gt_delta", "p_lt_neg_delta")
  missing_cols <- setdiff(req, colnames(post))
  if (length(missing_cols) > 0L) {
    stop("`post` missing required columns: ",
         paste(shQuote(missing_cols), collapse = ", "),
         ". Run posterior_summary(fit) first.", call. = FALSE)
  }
  if (!is.numeric(prob_cutoff) || length(prob_cutoff) != 1L ||
      prob_cutoff <= 0 || prob_cutoff >= 1) {
    stop("`prob_cutoff` must be in (0, 1).", call. = FALSE)
  }

  hyper <- !is.na(post$p_gt_delta)     & post$p_gt_delta     >= prob_cutoff
  hypo  <- !is.na(post$p_lt_neg_delta) & post$p_lt_neg_delta >= prob_cutoff

  cls <- rep("inconclusive", nrow(post))
  cls[hyper] <- "hypermethylated"
  cls[hypo]  <- "hypomethylated"

  both <- hyper & hypo
  if (any(both)) {
    favor_hyper <- both & post$p_gt_delta >= post$p_lt_neg_delta
    cls[favor_hyper]           <- "hypermethylated"
    cls[both & !favor_hyper]   <- "hypomethylated"
  }

  # NA rows (failed fits) stay inconclusive
  cls[is.na(post$p_gt_delta) | is.na(post$p_lt_neg_delta)] <- "inconclusive"

  post$classification <- factor(
    cls,
    levels = c("hypermethylated", "hypomethylated", "inconclusive")
  )
  attr(post, "delta")       <- delta
  attr(post, "prob_cutoff") <- prob_cutoff
  post
}
