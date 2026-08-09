# posterior_summary() is the bridge between a backend fit and the
# classification rule. Its column contract is what downstream code and the
# BreadFit results table both depend on.

POST_COLS <- c("region_id", "n", "mean_effect", "median_effect",
               "ci_lo", "ci_hi", "df", "scale", "p_pos", "p_neg",
               "p_gt_delta", "p_lt_neg_delta", "error")

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
  ok   <- !is.na(post$p_pos)
  expect_equal(post$p_pos[ok] + post$p_neg[ok], rep(1, sum(ok)),
               tolerance = 1e-8)
  expect_true(all(post$p_gt_delta[ok]     <= post$p_pos[ok] + 1e-8))
  expect_true(all(post$p_lt_neg_delta[ok] <= post$p_neg[ok] + 1e-8))
  for (col in c("p_pos", "p_neg", "p_gt_delta", "p_lt_neg_delta")) {
    expect_true(all(post[[col]][ok] >= 0 & post[[col]][ok] <= 1))
  }
})

test_that("a larger delta cannot increase the directional probabilities", {
  fit   <- .toy_fit()
  small <- posterior_summary(fit, delta = 0.05)
  large <- posterior_summary(fit, delta = 0.50)
  ok    <- !is.na(small$p_gt_delta)
  expect_true(all(large$p_gt_delta[ok]     <= small$p_gt_delta[ok] + 1e-12))
  expect_true(all(large$p_lt_neg_delta[ok] <= small$p_lt_neg_delta[ok] + 1e-12))
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
