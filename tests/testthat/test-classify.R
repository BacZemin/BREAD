# classify_regions() applies the decision rule; it does not recompute
# posterior probabilities. These tests drive it with hand-built posterior
# tables so the boundaries are exact rather than approximate.

`%|NULL|%` <- function(a, b) if (is.null(a)) b else a

.fake_post <- function(p_gt, p_lt, ids = NULL) {
  n <- length(p_gt)
  data.frame(
    region_id      = ids %|NULL|% sprintf("r%02d", seq_len(n)),
    p_gt_delta     = p_gt,
    p_lt_neg_delta = p_lt,
    stringsAsFactors = FALSE
  )
}

test_that("classification levels are exact and in a fixed order", {
  out <- classify_regions(.fake_post(c(0.99, 0.01, 0.5), c(0.01, 0.99, 0.5)))
  expect_s3_class(out$classification, "factor")
  expect_identical(
    levels(out$classification),
    c("hypermethylated", "hypomethylated", "inconclusive")
  )
  expect_identical(
    as.character(out$classification),
    c("hypermethylated", "hypomethylated", "inconclusive")
  )
})

test_that("the cutoff comparison is inclusive (>=), not strict", {
  at    <- classify_regions(.fake_post(0.95, 0.0), prob_cutoff = 0.95)
  below <- classify_regions(.fake_post(0.95 - 1e-12, 0.0), prob_cutoff = 0.95)
  expect_identical(as.character(at$classification),    "hypermethylated")
  expect_identical(as.character(below$classification), "inconclusive")
})

test_that("a stricter cutoff only ever moves regions to inconclusive", {
  set.seed(42)
  post   <- .fake_post(runif(200), runif(200))
  loose  <- classify_regions(post, prob_cutoff = 0.80)$classification
  strict <- classify_regions(post, prob_cutoff = 0.99)$classification
  moved  <- loose != strict
  expect_true(all(strict[moved] == "inconclusive"))
})

test_that("failed fits (NA probabilities) become inconclusive, never NA", {
  out <- classify_regions(.fake_post(c(NA, 0.99), c(NA, 0.001)))
  expect_false(any(is.na(out$classification)))
  expect_identical(as.character(out$classification)[1], "inconclusive")
})

test_that("when both tails clear the cutoff the larger probability wins", {
  out <- classify_regions(.fake_post(c(0.60, 0.40), c(0.40, 0.60)),
                          prob_cutoff = 0.30)
  expect_identical(as.character(out$classification),
                   c("hypermethylated", "hypomethylated"))
})

test_that("delta and prob_cutoff are recorded as attributes", {
  out <- classify_regions(.fake_post(0.99, 0.0), delta = 0.25,
                          prob_cutoff = 0.9)
  expect_equal(attr(out, "delta"), 0.25)
  expect_equal(attr(out, "prob_cutoff"), 0.9)
})

test_that("bad input is rejected with an informative error", {
  expect_error(classify_regions("nope"), "data.frame")
  expect_error(classify_regions(data.frame(x = 1)), "missing required columns")
  expect_error(classify_regions(.fake_post(0.9, 0.1), prob_cutoff = 1),
               "must be in \\(0, 1\\)")
  expect_error(classify_regions(.fake_post(0.9, 0.1), prob_cutoff = 0),
               "must be in \\(0, 1\\)")
})
