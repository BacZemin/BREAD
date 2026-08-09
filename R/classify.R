#' Classify regions as hyper / hypo / unchanged / inconclusive
#'
#' Applies the BREAD decision rule to the output of [posterior_summary()]:
#' - `hypermethylated` if `prob_hyper >= prob_cutoff`
#' - `hypomethylated`  if `prob_hypo  >= prob_cutoff`
#' - `unchanged`       if `prob_rope  >= rope_cutoff`
#' - `inconclusive`    otherwise
#'
#' @section Why `unchanged` is a separate class:
#' `inconclusive` used to absorb two entirely different situations: a region
#' whose posterior sits tightly inside the region of practical equivalence
#' (strong evidence of *no* change) and a region whose posterior is so diffuse
#' that nothing can be said. Collapsing them discards the one claim a p-value
#' structurally cannot make — that a region is *demonstrably* unmoved at the
#' stated `delta`. `unchanged` means "practically unchanged at this `delta`",
#' not "identical"; `inconclusive` now means only what its name says.
#'
#' @section Mutual exclusivity:
#' `prob_hyper`, `prob_hypo` and `prob_rope` partition the posterior, so they
#' sum to 1. Two of them can therefore clear their thresholds simultaneously
#' only if the two thresholds sum to no more than 1 — impossible at any
#' sensible setting (0.95 + 0.95 > 1). Should you set thresholds that low, the
#' largest of the qualifying probabilities wins, with ties resolved
#' hyper > hypo > unchanged.
#'
#' @param post Output of [posterior_summary()].
#' @param delta Effect-size threshold used for the rule. Default `0.10`. Stored
#'   as an attribute; does not re-evaluate the posterior probabilities (those
#'   must have been computed at this same `delta` upstream).
#' @param prob_cutoff Posterior probability cutoff for a *directional* call.
#'   Default `0.95`.
#' @param rope_cutoff Posterior probability cutoff for an *equivalence* call.
#'   Defaults to `prob_cutoff`. Worth setting independently: concluding
#'   equivalence requires the posterior to fit entirely inside
#'   \eqn{[-\delta, +\delta]}, a far stricter demand than a directional call,
#'   and at small n almost nothing reaches 0.95. Loosening it should not
#'   require loosening the discovery threshold too.
#'
#' @return The input `data.frame` with an added `classification` factor column
#'   (levels: `hypermethylated`, `hypomethylated`, `unchanged`,
#'   `inconclusive`). Attributes `delta`, `prob_cutoff` and `rope_cutoff` are
#'   updated. If `post` has no `prob_rope` column it is derived as
#'   `1 - prob_hyper - prob_hypo`.
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
#'
#' # Relax only the equivalence bar, leaving discovery untouched
#' table(classify_regions(post, rope_cutoff = 0.80)$classification)
#' @export
classify_regions <- function(post, delta = 0.10, prob_cutoff = 0.95,
                             rope_cutoff = prob_cutoff) {
  if (!is.data.frame(post)) {
    stop("`post` must be a data.frame from posterior_summary().",
         call. = FALSE)
  }
  req <- c("prob_hyper", "prob_hypo")
  missing_cols <- setdiff(req, colnames(post))
  if (length(missing_cols) > 0L) {
    stop("`post` missing required columns: ",
         paste(shQuote(missing_cols), collapse = ", "),
         ". Run posterior_summary(fit) first.", call. = FALSE)
  }
  .check_cutoff(prob_cutoff, "prob_cutoff")
  .check_cutoff(rope_cutoff, "rope_cutoff")

  # Hand-built posterior tables are a legitimate input; derive the ROPE mass
  # when it is absent. When present it is authoritative.
  prob_rope <- post$prob_rope
  if (is.null(prob_rope)) {
    prob_rope <- pmin(pmax(1 - post$prob_hyper - post$prob_hypo, 0), 1)
  }

  hyper <- !is.na(post$prob_hyper) & post$prob_hyper >= prob_cutoff
  hypo  <- !is.na(post$prob_hypo)  & post$prob_hypo  >= prob_cutoff
  rope  <- !is.na(prob_rope)       & prob_rope       >= rope_cutoff

  cls <- rep("inconclusive", nrow(post))
  cls[rope]  <- "unchanged"
  cls[hypo]  <- "hypomethylated"
  cls[hyper] <- "hypermethylated"

  # Only reachable when the two relevant cutoffs sum to <= 1; see @section.
  multi <- (hyper + hypo + rope) > 1L
  if (any(multi)) {
    cand <- cbind(
      ifelse(hyper[multi], post$prob_hyper[multi], -Inf),
      ifelse(hypo[multi],  post$prob_hypo[multi],  -Inf),
      ifelse(rope[multi],  prob_rope[multi],       -Inf)
    )
    # which.max takes the first maximum, so exact ties resolve
    # hyper > hypo > unchanged -- matching the pre-4-level tie-break.
    winner <- apply(cand, 1L, which.max)
    cls[multi] <- c("hypermethylated", "hypomethylated", "unchanged")[winner]
  }

  # NA rows (failed fits) are inconclusive, never unchanged: no posterior
  # means no evidence of equivalence either.
  cls[is.na(post$prob_hyper) | is.na(post$prob_hypo)] <- "inconclusive"

  post$classification <- factor(cls, levels = .BREAD_LEVELS)
  attr(post, "delta")       <- delta
  attr(post, "prob_cutoff") <- prob_cutoff
  attr(post, "rope_cutoff") <- rope_cutoff
  post
}

# Classification levels, in display order. `unchanged` sits third so that
# levels()[1:2] remain the directional pair (bread_kycg()'s `which` default
# and any positional indexing rely on it) and `inconclusive` stays last as
# the residual bucket.
.BREAD_LEVELS <- c("hypermethylated", "hypomethylated",
                   "unchanged", "inconclusive")

.check_cutoff <- function(x, nm) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x <= 0 || x >= 1) {
    stop("`", nm, "` must be in (0, 1).", call. = FALSE)
  }
  invisible(TRUE)
}
