#' Internal utilities
#'
#' @return Varies by helper: the scale transforms return numeric vectors,
#'   the detection helpers return a character scalar.
#'
#' @name BREAD-utils
#' @keywords internal
NULL

#' Beta -> M transform, clamped away from 0/1.
#'
#' @param beta Numeric vector of beta values.
#' @param eps Clamping tolerance keeping values off the 0/1 asymptotes.
#' @return Numeric vector of M-values.
#' @keywords internal
.beta_to_m <- function(beta, eps = 1e-6) {
  beta <- pmin(pmax(beta, eps), 1 - eps)
  log2(beta / (1 - beta))
}

#' M -> Beta transform.
#'
#' @param m Numeric vector of M-values.
#' @return Numeric vector of beta values in (0, 1).
#' @keywords internal
.m_to_beta <- function(m) {
  2^m / (2^m + 1)
}

# MetBrewer "Cross" palette (Blake Robert Mills). Embedded as hex to avoid
# a runtime dependency on the MetBrewer package.
.bread_cross <- c(
  "#c969a1",  # pink
  "#ce4441",  # red
  "#ee8577",  # salmon
  "#eb7926",  # orange
  "#ffbb44",  # amber
  "#859b6c",  # olive
  "#62929a",  # teal
  "#004f63",  # navy teal
  "#122451"   # deep navy
)

# Classification palette (warm = hyper gain, cool = hypo loss, neutral = inconclusive)
.col_classification <- c(
  hypermethylated = "#ce4441",
  hypomethylated  = "#004f63",
  inconclusive    = "#859b6c"
)

# Binary group palette — used by plot_region_data() for the contrast variable
.col_group <- c("#62929a", "#eb7926")

#' BREAD color palettes (MetBrewer "Cross")
#'
#' Colorblind-reasonable palette derived from the MetBrewer `Cross` palette
#' (Blake Robert Mills). Embedded as hex values so there is no runtime
#' dependency on the `MetBrewer` package.
#'
#' @param which One of:
#'   - `"classification"` : named 3-vector for hyper/hypo/inconclusive
#'   - `"group"`          : unnamed 2-vector for binary contrast plots
#'   - `"cross"`          : full 9-color palette
#'
#' @return A character vector of hex colors (named where applicable).
#' @examples
#' bread_colors("classification")
#' bread_colors("group")
#' bread_colors("cross")
#' @export
bread_colors <- function(which = c("classification", "group", "cross")) {
  which <- match.arg(which)
  switch(
    which,
    classification = .col_classification,
    group          = .col_group,
    cross          = .bread_cross
  )
}

# Internal: detect input scale from assay values.
# Returns "Beta" if all non-NA values are in [0, 1], else "M".
.detect_input_scale <- function(x) {
  rng <- suppressWarnings(range(x, na.rm = TRUE))
  if (all(is.finite(rng)) && rng[1L] >= 0 && rng[2L] <= 1) "Beta" else "M"
}

# Internal: pick a reasonable assay name when the user leaves it NULL.
# Priority: "M", "betas", "Beta", "beta", then first assay.
.detect_assay_name <- function(se) {
  an <- SummarizedExperiment::assayNames(se)
  if (is.null(an) || length(an) == 0L)
    stop("`se` has no named assays; supply `assay_name` explicitly.",
         call. = FALSE)
  for (candidate in c("M", "betas", "Beta", "beta"))
    if (candidate %in% an) return(candidate)
  an[1L]
}
