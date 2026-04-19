#' Extract per-region posterior summaries
#'
#' Given the output of [fit_bread_summary()] (conjugate backend) or
#' [fit_bread_brms()] (brms backend), compute the marginal posterior
#' location / scale / credible interval of the contrast coefficient and the
#' posterior probabilities used by [classify_regions()].
#'
#' @section Path selection:
#' Each per-region fit carries either
#' - `$mu_n`, `$Lambda_n_inv`, `$a_n`, `$b_n` (conjugate analytical path), in
#'   which case the marginal posterior of the contrast coefficient is a
#'   location-scale Student-t with `df = 2 * a_n`, or
#' - `$draws` (brms empirical path), in which case posterior quantities are
#'   computed from the MCMC draws directly.
#'
#' Columns in the returned data frame are the same in both cases.
#'
#' @param fit Output of [fit_bread_summary()] or [fit_bread_brms()].
#' @param delta Effect-size threshold on the M-value scale. Default `0.10`.
#' @param ci Credible-interval mass. Default `0.95`.
#'
#' @return A `data.frame` with one row per region and columns:
#'   `region_id`, `n`, `mean_effect`, `median_effect`, `ci_lo`, `ci_hi`,
#'   `df`, `scale`, `p_pos`, `p_neg`, `p_gt_delta`, `p_lt_neg_delta`, `error`.
#'   `df` is `NA_real_` for the empirical path.
#'
#' @importFrom stats pt qt median quantile sd
#' @export
posterior_summary <- function(fit, delta = 0.10, ci = 0.95) {
  if (!is.list(fit) || is.null(fit$fits) || is.null(fit$contrast_idx)) {
    stop("`fit` must be the output of fit_bread_summary() or fit_bread_brms().",
         call. = FALSE)
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

    na_row <- data.frame(
      region_id      = rid,
      n              = if (is.null(f$n)) NA_integer_ else f$n,
      mean_effect    = NA_real_, median_effect = NA_real_,
      ci_lo          = NA_real_, ci_hi         = NA_real_,
      df             = NA_real_, scale         = NA_real_,
      p_pos          = NA_real_, p_neg         = NA_real_,
      p_gt_delta     = NA_real_, p_lt_neg_delta = NA_real_,
      error          = if (is.null(f$error)) NA_character_ else f$error,
      stringsAsFactors = FALSE
    )
    if (!is.null(f$error) && !is.na(f$error)) return(na_row)

    # Empirical path (brms)
    if (!is.null(f$draws) && length(f$draws) > 0L) {
      d <- f$draws
      q <- stats::quantile(d, c(alpha, 1 - alpha), names = FALSE, type = 7L)
      return(data.frame(
        region_id      = rid,
        n              = f$n,
        mean_effect    = mean(d),
        median_effect  = stats::median(d),
        ci_lo          = q[1L],
        ci_hi          = q[2L],
        df             = NA_real_,
        scale          = stats::sd(d),
        p_pos          = mean(d > 0),
        p_neg          = mean(d < 0),
        p_gt_delta     = mean(d >  delta),
        p_lt_neg_delta = mean(d < -delta),
        error          = NA_character_,
        stringsAsFactors = FALSE
      ))
    }

    # Analytical path (conjugate NIG)
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
