#' Extract per-region posterior summaries
#'
#' Given the output of [fit_bread_summary()], compute the marginal posterior
#' location/scale/df of the contrast coefficient for each region, along with
#' directional and threshold probabilities used by [classify_regions()].
#'
#' @section Math:
#' Under the Normal-Inverse-Gamma conjugate update, the marginal posterior of
#' the contrast coefficient is a location-scale Student-t:
#' \deqn{\beta_{r,k} \mid y \sim t_{2 a_n}\!\left(\mu_{n,k},\ \tfrac{b_n}{a_n} [\Lambda_n^{-1}]_{kk}\right).}
#' Posterior probabilities are computed as `P(beta > c) = 1 - pt((c - mu)/s, df)`.
#'
#' @param fit Output of [fit_bread_summary()].
#' @param delta Effect-size threshold on the M-value scale. Default `0.10`.
#' @param ci Credible-interval mass. Default `0.95`.
#'
#' @return A `data.frame` with one row per region and columns:
#'   `region_id`, `n`, `mean_effect`, `median_effect`, `ci_lo`, `ci_hi`,
#'   `df`, `scale`, `p_pos`, `p_neg`, `p_gt_delta`, `p_lt_neg_delta`, `error`.
#'   Attributes: `delta`, `ci`, `contrast`.
#'
#' @importFrom stats pt qt
#' @export
posterior_summary <- function(fit, delta = 0.10, ci = 0.95) {
  if (!is.list(fit) || is.null(fit$fits) || is.null(fit$contrast_idx)) {
    stop("`fit` must be the output of fit_bread_summary().", call. = FALSE)
  }
  if (!is.numeric(delta) || length(delta) != 1L || delta < 0) {
    stop("`delta` must be a non-negative scalar.", call. = FALSE)
  }
  if (!is.numeric(ci) || length(ci) != 1L || ci <= 0 || ci >= 1) {
    stop("`ci` must be in (0, 1).", call. = FALSE)
  }

  k     <- fit$contrast_idx
  alpha <- (1 - ci) / 2

  rows <- lapply(seq_along(fit$fits), function(i) {
    f   <- fit$fits[[i]]
    rid <- fit$region_ids[i]
    if (!is.na(f$error)) {
      return(data.frame(
        region_id      = rid,
        n              = f$n,
        mean_effect    = NA_real_,
        median_effect  = NA_real_,
        ci_lo          = NA_real_,
        ci_hi          = NA_real_,
        df             = NA_real_,
        scale          = NA_real_,
        p_pos          = NA_real_,
        p_neg          = NA_real_,
        p_gt_delta     = NA_real_,
        p_lt_neg_delta = NA_real_,
        error          = f$error,
        stringsAsFactors = FALSE
      ))
    }
    nu  <- 2 * f$a_n
    mu  <- f$mu_n[k]
    s2  <- (f$b_n / f$a_n) * f$Lambda_n_inv[k, k]
    s   <- sqrt(max(s2, .Machine$double.eps))

    q_lo <- mu + s * stats::qt(alpha,     df = nu)
    q_hi <- mu + s * stats::qt(1 - alpha, df = nu)

    p_pos  <- 1 - stats::pt((0      - mu) / s, df = nu)
    p_neg  <-     stats::pt((0      - mu) / s, df = nu)
    p_gt_d <- 1 - stats::pt(( delta - mu) / s, df = nu)
    p_lt_d <-     stats::pt((-delta - mu) / s, df = nu)

    data.frame(
      region_id      = rid,
      n              = f$n,
      mean_effect    = mu,
      median_effect  = mu,
      ci_lo          = q_lo,
      ci_hi          = q_hi,
      df             = nu,
      scale          = s,
      p_pos          = p_pos,
      p_neg          = p_neg,
      p_gt_delta     = p_gt_d,
      p_lt_neg_delta = p_lt_d,
      error          = NA_character_,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "delta")    <- delta
  attr(out, "ci")       <- ci
  attr(out, "contrast") <- fit$contrast
  out
}
