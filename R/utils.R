#' Internal utilities
#'
#' @name BREAD-utils
#' @keywords internal
NULL

#' Beta -> M transform, clamped away from 0/1.
#' @keywords internal
.beta_to_m <- function(beta, eps = 1e-6) {
  beta <- pmin(pmax(beta, eps), 1 - eps)
  log2(beta / (1 - beta))
}

#' M -> Beta transform.
#' @keywords internal
.m_to_beta <- function(m) {
  2^m / (2^m + 1)
}
