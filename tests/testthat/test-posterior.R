# posterior_summary() is the bridge between a backend fit and the
# classification rule. Its column contract is what downstream code and the
# BreadFit results table both depend on.

# Spelled out independently of the package's own .POST_COLS -- that is the
# point of a contract test.
POST_COLS <- c("region_id", "n", "mean_effect", "median_effect",
               "ci_lo", "ci_hi", "df", "scale", "prob_pos", "prob_neg",
               "prob_hyper", "prob_hypo", "prob_rope",
               "ref_beta", "mean_dbeta", "dbeta_lo", "dbeta_hi", "delta_beta",
               "error")

.toy_fit <- function() {
  fit_bread(.make_toy_signal_se(), .make_toy_features(), ~ group)
}

test_that("the returned column contract is stable", {
  post <- posterior_summary(.toy_fit())
  expect_s3_class(post, "data.frame")
  expect_identical(colnames(post), POST_COLS)
})

test_that("a BreadFit is accepted and matches the internal model list", {
  fit <- .toy_fit()
  expect_equal(posterior_summary(fit), posterior_summary(fit@model))
})

test_that("a wider ci widens every interval but moves no point estimate", {
  fit <- .toy_fit()
  p95 <- posterior_summary(fit, ci = 0.95)
  p99 <- posterior_summary(fit, ci = 0.99)
  ok  <- !is.na(p95$ci_lo) & !is.na(p99$ci_lo)
  expect_true(all(p99$ci_lo[ok] <= p95$ci_lo[ok]))
  expect_true(all(p99$ci_hi[ok] >= p95$ci_hi[ok]))
  expect_equal(p95$mean_effect, p99$mean_effect)
})

test_that("directional probabilities are coherent", {
  post <- posterior_summary(.toy_fit())
  ok   <- !is.na(post$prob_pos)
  expect_equal(post$prob_pos[ok] + post$prob_neg[ok], rep(1, sum(ok)),
               tolerance = 1e-8)
  expect_true(all(post$prob_hyper[ok]     <= post$prob_pos[ok] + 1e-8))
  expect_true(all(post$prob_hypo[ok] <= post$prob_neg[ok] + 1e-8))
  for (col in c("prob_pos", "prob_neg", "prob_hyper", "prob_hypo")) {
    expect_true(all(post[[col]][ok] >= 0 & post[[col]][ok] <= 1))
  }
})

test_that("a larger delta cannot increase the directional probabilities", {
  fit   <- .toy_fit()
  small <- posterior_summary(fit, delta = 0.05)
  large <- posterior_summary(fit, delta = 0.50)
  ok    <- !is.na(small$prob_hyper)
  expect_true(all(large$prob_hyper[ok]     <= small$prob_hyper[ok] + 1e-12))
  expect_true(all(large$prob_hypo[ok] <= small$prob_hypo[ok] + 1e-12))
})

test_that("delta, ci and contrast are recorded as attributes", {
  post <- posterior_summary(.toy_fit(), delta = 0.2, ci = 0.9)
  expect_equal(attr(post, "delta"), 0.2)
  expect_equal(attr(post, "ci"), 0.9)
  expect_true(is.character(attr(post, "contrast")))
})

test_that("bad input is rejected", {
  fit <- .toy_fit()
  expect_error(posterior_summary(list(a = 1)), "fit_bread_summary")
  expect_error(posterior_summary(fit, delta = -1), "non-negative")
  expect_error(posterior_summary(fit, ci = 0), "must be in \\(0, 1\\)")
  expect_error(posterior_summary(fit, ci = 1), "must be in \\(0, 1\\)")
})


# ---- prob_rope -------------------------------------------------------------

test_that("the three posterior masses partition the line", {
  post <- posterior_summary(.toy_fit())
  ok   <- !is.na(post$prob_hyper)
  expect_equal(post$prob_hyper[ok] + post$prob_hypo[ok] + post$prob_rope[ok],
               rep(1, sum(ok)), tolerance = 1e-12)
  expect_true(all(post$prob_rope[ok] >= 0 & post$prob_rope[ok] <= 1))
  expect_identical(is.na(post$prob_rope), is.na(post$prob_hyper))
})

test_that("a wider ROPE can only absorb more posterior mass", {
  fit   <- .toy_fit()
  small <- posterior_summary(fit, delta = 0.05)
  large <- posterior_summary(fit, delta = 0.50)
  ok    <- !is.na(small$prob_rope)
  expect_true(all(large$prob_rope[ok] >= small$prob_rope[ok] - 1e-12))
})

test_that("a zero-width ROPE holds no mass", {
  post <- posterior_summary(.toy_fit(), delta = 0)
  ok   <- !is.na(post$prob_rope)
  expect_equal(post$prob_rope[ok], rep(0, sum(ok)), tolerance = 1e-12)
})


# ---- beta-scale columns ----------------------------------------------------

test_that("beta columns are the linearisation of the M-scale columns", {
  post <- posterior_summary(.toy_fit(), delta = 0.10)
  ok   <- !is.na(post$ref_beta)
  expect_true(any(ok))
  expect_true(all(post$ref_beta[ok] > 0 & post$ref_beta[ok] < 1))

  k <- post$ref_beta[ok] * (1 - post$ref_beta[ok]) * log(2)
  expect_equal(post$mean_dbeta[ok], post$mean_effect[ok] * k)
  expect_equal(post$dbeta_lo[ok],   post$ci_lo[ok] * k)
  expect_equal(post$dbeta_hi[ok],   post$ci_hi[ok] * k)
  expect_equal(post$delta_beta[ok], 0.10 * k)
  expect_true(all(post$dbeta_lo[ok] <= post$dbeta_hi[ok]))
})

test_that("one multiplier keeps the beta scale consistent with the M scale", {
  # The whole reason for a single linearisation rather than an exact secant:
  # the beta comparison must never contradict the classification beside it.
  post <- posterior_summary(.toy_fit(), delta = 0.10)
  ok   <- !is.na(post$ref_beta)
  expect_identical(post$mean_effect[ok] > 0.10,
                   post$mean_dbeta[ok]  > post$delta_beta[ok])
})

test_that("ref_beta accepts a scalar and a named vector", {
  fit  <- .toy_fit()
  ids  <- posterior_summary(fit)$region_id

  flat <- posterior_summary(fit, ref_beta = 0.3)
  expect_equal(flat$ref_beta, rep(0.3, nrow(flat)))

  named <- stats::setNames(seq(0.2, 0.4, length.out = length(ids)), ids)
  byid  <- posterior_summary(fit, ref_beta = named)
  expect_equal(byid$ref_beta, unname(named[byid$region_id]))
})

test_that("ref_beta rejects impossible values and ambiguous lengths", {
  fit <- .toy_fit()
  expect_error(posterior_summary(fit, ref_beta = 0),   "must be in \\(0, 1\\)")
  expect_error(posterior_summary(fit, ref_beta = 1.2), "must be in \\(0, 1\\)")
  expect_error(posterior_summary(fit, ref_beta = c(0.3, 0.4, 0.5)),
               "named by region_id")
})

test_that("beta columns are NA when no region matrix is available", {
  m <- .toy_fit()@model
  m$region_mat <- NULL
  post <- posterior_summary(m)
  expect_true(all(is.na(post$ref_beta)))
  expect_true(all(is.na(post$mean_dbeta)))
  # ... but the M-scale results are untouched
  expect_false(all(is.na(post$mean_effect)))
})

test_that("pc1 scores get no beta translation", {
  fit <- suppressMessages(
    fit_bread(.make_toy_signal_se(), .make_toy_features(), ~ group,
              summary_fun = "pc1")
  )
  expect_message(post <- posterior_summary(fit), "PC1 scores are not M-values")
  expect_true(all(is.na(post$ref_beta)))
  expect_true(all(is.na(post$delta_beta)))
  # An explicit anchor overrides the refusal
  post2 <- posterior_summary(fit, ref_beta = 0.5)
  expect_true(all(!is.na(post2$delta_beta)))
})


# ---- exported converters ---------------------------------------------------

test_that("bread_delta_beta matches the documented arithmetic", {
  # The corrected value: 0.10 M-units is ~1.7 percentage points at beta = 0.5,
  # NOT 3.5 (that is the width of the full +/-delta window).
  expect_equal(bread_delta_beta(0.10), 0.10 * 0.25 * log(2))
  expect_equal(round(bread_delta_beta(0.10), 3), 0.017)

  # Scale dependence: the translation shrinks toward the extremes
  v <- bread_delta_beta(0.10, ref_beta = c(0.5, 0.2, 0.1))
  expect_true(all(diff(v) < 0))

  # Round trip
  expect_equal(bread_delta_m(bread_delta_beta(0.10, 0.3), 0.3), 0.10)
  expect_equal(bread_delta_beta(bread_delta_m(0.02, 0.4), 0.4), 0.02)

  expect_error(bread_delta_beta(0.1, ref_beta = 0), "must be in \\(0, 1\\)")
  expect_error(bread_delta_m(0.1, ref_beta = 1),    "must be in \\(0, 1\\)")
})
