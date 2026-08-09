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
#' @section Equivalence (`prob_rope`):
#' `prob_hyper`, `prob_hypo` and `prob_rope` are mutually exclusive and
#' exhaustive: they are the posterior mass above `+delta`, below `-delta`, and
#' inside the region of practical equivalence \eqn{[-\delta, +\delta]}, and
#' they sum to 1. `prob_rope` is what lets BREAD state that a region is
#' *confidently unchanged* rather than merely undetected — a claim no p-value
#' can make. See [classify_regions()].
#'
#' @section Beta-scale columns:
#' BREAD models M-values, but reports a beta-scale translation of the effect
#' and the ROPE half-width via the local linearisation
#' \eqn{d\beta \approx dM \cdot \beta(1-\beta)\ln 2}, anchored per region at
#' `ref_beta`. By default `ref_beta` is the region's own mean methylation,
#' back-transformed from the mean M-value — well defined for every design and
#' contrast type, unlike the reference level of a factor. The same multiplier
#' is applied to the effect, both interval bounds and `delta`, so the
#' beta-scale comparison can never contradict the M-scale classification
#' beside it. See [bread_delta_beta()].
#'
#' @param fit A [BreadFit] (as returned by [fit_bread()]), or the internal
#'   model list from [fit_bread_summary()] / [fit_bread_brms()].
#' @param delta Effect-size threshold on the M-value scale. Default `0.10`.
#' @param ci Credible-interval mass. Default `0.95`.
#' @param ref_beta Reference methylation level for the beta-scale columns.
#'   `NULL` (default) derives it per region from the fitted region matrix.
#'   Otherwise a single value applied to every region, or a numeric vector
#'   named by `region_id`. Values must lie in (0, 1).
#'
#' @return A `data.frame` with one row per region and columns:
#'   `region_id`, `n`, `mean_effect`, `median_effect`, `ci_lo`, `ci_hi`,
#'   `df`, `scale`, `prob_pos`, `prob_neg`, `prob_hyper`, `prob_hypo`,
#'   `prob_rope`, `ref_beta`, `mean_dbeta`, `dbeta_lo`, `dbeta_hi`,
#'   `delta_beta`, `error`.
#'   `df` is `NA_real_` for the empirical path. The beta-scale columns are
#'   `NA_real_` when no region matrix is available, or when
#'   `summary_fun = "pc1"` (PC1 scores are not M-values).
#'
#' @importFrom stats pt qt median quantile sd
#' @importFrom methods is
#' @examples
#' suppressPackageStartupMessages(library(SummarizedExperiment))
#'
#' se  <- readRDS(system.file("extdata", "vitc_ag06561.rds", package = "BREAD"))
#' reg <- readRDS(system.file("extdata", "vitc_regions.rds", package = "BREAD"))
#' se_ctrl <- se[, se$condition == "ctrl"]
#'
#' fit <- fit_bread(se_ctrl, reg, ~ passage)
#'
#' # A BreadFit is accepted directly; the internal model list also works.
#' post <- posterior_summary(fit)
#' head(post)
#'
#' # A wider credible interval
#' head(posterior_summary(fit, ci = 0.99))
#'
#' # Posterior mass inside the region of practical equivalence
#' summary(posterior_summary(fit)$prob_rope)
#' @export
posterior_summary <- function(fit, delta = 0.10, ci = 0.95, ref_beta = NULL) {
  # Accept the user-facing object as well as the internal model list, so
  # callers never have to reach into the BreadFit's slots.
  if (methods::is(fit, "BreadFit")) fit <- fit@model
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
      prob_pos       = NA_real_, prob_neg  = NA_real_,
      prob_hyper     = NA_real_, prob_hypo = NA_real_,
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
        prob_pos       = mean(d > 0),
        prob_neg       = mean(d < 0),
        prob_hyper     = mean(d >  delta),
        prob_hypo      = mean(d < -delta),
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

    prob_pos   <- 1 - stats::pt((0      - mu) / s, df = nu)
    prob_neg   <-     stats::pt((0      - mu) / s, df = nu)
    prob_hyper <- 1 - stats::pt(( delta - mu) / s, df = nu)
    prob_hypo  <-     stats::pt((-delta - mu) / s, df = nu)

    data.frame(
      region_id      = rid,
      n              = f$n,
      mean_effect    = mu,
      median_effect  = mu,
      ci_lo          = q_lo,
      ci_hi          = q_hi,
      df             = nu,
      scale          = s,
      prob_pos       = prob_pos,
      prob_neg       = prob_neg,
      prob_hyper     = prob_hyper,
      prob_hypo      = prob_hypo,
      error          = NA_character_,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  # Posterior mass inside the ROPE. Computed here rather than inside each of
  # the three data.frame branches for two reasons: NA propagates through
  # pmin/pmax, so failed fits get NA for free; and every future backend
  # inherits it. The complement is used deliberately -- it loses *relative*
  # accuracy only where prob_rope is near zero, i.e. where the answer is "not
  # unchanged" either way, and the decision rule compares against ~0.95, where
  # absolute accuracy is what matters. Do not "fix" this into a direct
  # integral. The clamp catches the analytical path rounding to ~-1e-17.
  out$prob_rope <- pmin(pmax(1 - out$prob_hyper - out$prob_hypo, 0), 1)

  # Beta-scale translation, anchored per region.
  rb <- .resolve_ref_beta(fit, out$region_id, ref_beta)
  k  <- .dbeta_per_dm(rb)
  out$ref_beta   <- rb
  out$mean_dbeta <- out$mean_effect * k
  out$dbeta_lo   <- out$ci_lo * k
  out$dbeta_hi   <- out$ci_hi * k
  out$delta_beta <- delta * k

  out <- out[, .POST_COLS, drop = FALSE]
  attr(out, "delta")    <- delta
  attr(out, "ci")       <- ci
  attr(out, "contrast") <- fit$contrast
  out
}

# Canonical column order for posterior_summary(). `error` stays last.
.POST_COLS <- c(
  "region_id", "n", "mean_effect", "median_effect", "ci_lo", "ci_hi",
  "df", "scale",
  "prob_pos", "prob_neg", "prob_hyper", "prob_hypo", "prob_rope",
  "ref_beta", "mean_dbeta", "dbeta_lo", "dbeta_hi", "delta_beta",
  "error"
)

# Resolve the per-region beta anchor for the M -> beta linearisation.
#
# Default: the region's own mean methylation, as the back-transform of its
# mean M-value. This is defined for every design (a factor's reference level
# is not -- consider `~ passage`, an interaction contrast, or `~ 0 + group`),
# needs nothing beyond state the model already carries, and sits where the
# linearisation error is smallest on average. Note it is m_to_beta(mean(M)),
# not mean(beta); those differ by Jensen, and the former is the right anchor
# for a linearisation of an M-scale model.
.resolve_ref_beta <- function(model, region_ids, ref_beta = NULL) {
  n <- length(region_ids)

  if (!is.null(ref_beta)) {
    .check_ref_beta(ref_beta)
    if (length(ref_beta) == 1L) return(rep(as.numeric(ref_beta), n))
    if (!is.null(names(ref_beta))) {
      return(unname(as.numeric(ref_beta[match(region_ids, names(ref_beta))])))
    }
    if (length(ref_beta) != n) {
      stop("`ref_beta` must be length 1, named by region_id, or length ",
           n, " (one per region); got ", length(ref_beta),
           " unnamed values.", call. = FALSE)
    }
    return(as.numeric(ref_beta))
  }

  rm_ <- model$region_mat
  if (is.null(rm_) || !is.matrix(rm_) || nrow(rm_) == 0L) {
    return(rep(NA_real_, n))
  }
  if (identical(attr(rm_, "summary_fun"), "pc1")) {
    message("Beta-scale columns are NA under `summary_fun = \"pc1\"`: ",
            "PC1 scores are not M-values, so no beta translation exists. ",
            "Pass `ref_beta` explicitly to override.")
    return(rep(NA_real_, n))
  }

  idx <- match(region_ids, rownames(rm_))
  mM  <- rep(NA_real_, n)
  ok  <- !is.na(idx)
  if (any(ok)) mM[ok] <- rowMeans(rm_[idx[ok], , drop = FALSE], na.rm = TRUE)
  out <- .m_to_beta(mM)
  out[!is.finite(out)] <- NA_real_
  out
}
