#' Translate effect sizes between the M and beta scales
#'
#' BREAD models M-values, so `delta` and every effect estimate are on the
#' M scale. Biologists generally think in beta (proportion methylated). These
#' helpers convert between the two using the local linearisation of
#' \eqn{M = \log_2(\beta / (1 - \beta))} at a reference beta:
#'
#' \deqn{d\beta \approx dM \cdot \beta(1 - \beta) \ln 2}
#'
#' @section Why the translation is not a single number:
#' The slope \eqn{\beta(1-\beta)\ln 2} is maximal at \eqn{\beta = 0.5} (0.173)
#' and shrinks toward the extremes (0.111 at \eqn{\beta = 0.2}, 0.062 at
#' \eqn{\beta = 0.1}). So the default `delta = 0.10` on the M scale means
#' \eqn{\Delta\beta \approx 0.017} at mid-methylation but only \eqn{\approx
#' 0.006} at \eqn{\beta = 0.1}. **0.017 is a ceiling, not a typical value.**
#' This is why BREAD reports a per-region `delta_beta` in [results()] rather
#' than accepting `delta` in beta units: a beta-defined threshold would
#' silently become a 20-fold wider equivalence region at the extremes.
#'
#' The default `ref_beta = 0.5` is the maximal-slope anchor. For a target
#' \eqn{\Delta\beta} it therefore returns the *smallest* `delta_m` that could
#' produce it — conservative when hunting for change, anti-conservative when
#' claiming equivalence. Anchor at your own data's methylation level when the
#' distinction matters.
#'
#' @param delta_m Effect size on the M-value scale.
#' @param delta_beta Effect size on the beta scale.
#' @param ref_beta Reference methylation level at which to linearise, in
#'   (0, 1). Default `0.5`.
#'
#' @return A numeric vector the length of the recycled inputs.
#'
#' @examples
#' # The default BREAD threshold, in beta units, at mid-methylation
#' bread_delta_beta(0.10)
#'
#' # ... and how much smaller it is toward the extremes
#' bread_delta_beta(0.10, ref_beta = c(0.5, 0.2, 0.1))
#'
#' # Going the other way: what delta_m gives a 2-percentage-point window?
#' bread_delta_m(0.02)
#'
#' # Round trip
#' bread_delta_m(bread_delta_beta(0.10, 0.3), 0.3)
#' @name bread_scale
NULL

# Local slope d(beta)/d(M) at a reference beta.
.dbeta_per_dm <- function(ref_beta) ref_beta * (1 - ref_beta) * log(2)

.check_ref_beta <- function(ref_beta) {
  if (!is.numeric(ref_beta) || length(ref_beta) == 0L) {
    stop("`ref_beta` must be a non-empty numeric vector.", call. = FALSE)
  }
  bad <- !is.na(ref_beta) & (ref_beta <= 0 | ref_beta >= 1)
  if (any(bad)) {
    offending <- ref_beta[bad]
    stop("`ref_beta` must be in (0, 1); got ",
         paste(offending[seq_len(min(3L, length(offending)))],
               collapse = ", "), ".",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' @rdname bread_scale
#' @export
bread_delta_beta <- function(delta_m, ref_beta = 0.5) {
  if (!is.numeric(delta_m)) stop("`delta_m` must be numeric.", call. = FALSE)
  .check_ref_beta(ref_beta)
  delta_m * .dbeta_per_dm(ref_beta)
}

#' @rdname bread_scale
#' @export
bread_delta_m <- function(delta_beta, ref_beta = 0.5) {
  if (!is.numeric(delta_beta)) {
    stop("`delta_beta` must be numeric.", call. = FALSE)
  }
  .check_ref_beta(ref_beta)
  delta_beta / .dbeta_per_dm(ref_beta)
}
