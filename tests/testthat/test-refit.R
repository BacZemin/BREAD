# refit_bread() exists so that label-permutation calibration does not have to
# recompute the region matrix, and does not have to reach into the namespace.
# The load-bearing property is therefore: the matrix is reused, never rebuilt.

.refit_fit <- function() {
  suppressMessages(fit_bread(.make_toy_signal_se(), .make_toy_features(),
                             ~ group))
}

test_that("a no-op refit reproduces the original exactly", {
  fit <- .refit_fit()
  re  <- refit_bread(fit)
  expect_equal(re@posterior, fit@posterior)
  expect_equal(results(re),  results(fit))
  expect_identical(re@params$contrast, fit@params$contrast)
})

test_that("nothing is re-summarized", {
  fit <- .refit_fit()
  re  <- refit_bread(fit, delta = 0.5)
  # The proof that mapping/summarization did not run again
  expect_identical(re@model$region_mat, fit@model$region_mat)
  expect_identical(re@mapping,  fit@mapping)
  expect_identical(re@features, fit@features)
})

test_that("re-thresholding matches doing it by hand", {
  fit  <- .refit_fit()
  re   <- refit_bread(fit, delta = 0.05, prob_cutoff = 0.8, rope_cutoff = 0.6,
                      ci = 0.9)
  hand <- classify_regions(
    posterior_summary(fit@model, delta = 0.05, ci = 0.9),
    delta = 0.05, prob_cutoff = 0.8, rope_cutoff = 0.6
  )
  expect_equal(results(re), hand)
  expect_equal(re@params$delta, 0.05)
  expect_equal(re@params$ci, 0.9)
  expect_equal(re@params$rope_cutoff, 0.6)
})

test_that("unspecified settings are inherited from the parent fit", {
  fit <- suppressMessages(
    fit_bread(.make_toy_signal_se(), .make_toy_features(), ~ group,
              delta = 0.3, prob_cutoff = 0.9, rope_cutoff = 0.7, ci = 0.8)
  )
  re <- refit_bread(fit)
  expect_equal(re@params$delta,       0.3)
  expect_equal(re@params$prob_cutoff, 0.9)
  expect_equal(re@params$rope_cutoff, 0.7)
  expect_equal(re@params$ci,          0.8)
})

test_that("a permuted colData changes the fit but not the region matrix", {
  fit <- .refit_fit()
  cd  <- as.data.frame(fit@model$coldata)
  cd$group <- rev(cd$group)

  re <- refit_bread(fit, colData = cd)
  expect_identical(re@model$region_mat, fit@model$region_mat)
  expect_false(isTRUE(all.equal(results(re)$mean_effect,
                                results(fit)$mean_effect)))
  expect_identical(re@diagnostics$refit_of, fit@diagnostics$timestamp)
})

test_that("colData is matched by rowname, not position", {
  fit <- .refit_fit()
  cd  <- as.data.frame(fit@model$coldata)
  shuffled <- cd[sample(nrow(cd)), , drop = FALSE]

  # Same information, different row order -- must give the same answer
  expect_equal(results(refit_bread(fit, colData = shuffled)),
               results(refit_bread(fit, colData = cd)))
})

test_that("a mis-sized or disjoint colData is rejected", {
  fit <- .refit_fit()
  cd  <- as.data.frame(fit@model$coldata)
  expect_error(refit_bread(fit, colData = cd[1:3, , drop = FALSE]),
               "rows but the region matrix has")

  bad <- cd; rownames(bad) <- paste0("X", seq_len(nrow(bad)))
  expect_error(refit_bread(fit, colData = bad), "do not cover every sample")
})

test_that("an unknown contrast lists the available coefficients", {
  fit <- .refit_fit()
  expect_error(refit_bread(fit, contrast = "groupNOPE"),
               "not found among design coefficients")
})

test_that("a rank-deficient design warns instead of silently regularising", {
  fit <- .refit_fit()
  cd  <- as.data.frame(fit@model$coldata)
  cd$dupe <- cd$group           # perfectly collinear with group
  expect_warning(refit_bread(fit, colData = cd, design = ~ group + dupe),
                 "rank deficient")
})

test_that("refit_bread rejects non-BreadFit input", {
  expect_error(refit_bread(list()), "must be a BreadFit")
})

test_that("a label-permutation null runs end to end through the public API", {
  # The workflow this function exists for, in miniature.
  fit <- .refit_fit()
  cd  <- as.data.frame(fit@model$coldata)
  rid <- results(fit)$region_id[1]
  obs <- results(fit)$mean_effect[results(fit)$region_id == rid]

  null <- withr::with_seed(7L, vapply(seq_len(24L), function(i) {
    cdp <- cd
    cdp$group <- sample(cdp$group)
    r <- results(refit_bread(fit, colData = cdp))
    r$mean_effect[r$region_id == rid]
  }, numeric(1)))

  expect_length(null, 24L)
  expect_false(anyNA(null))
  p <- (sum(abs(null) >= abs(obs)) + 1) / (length(null) + 1)
  expect_true(p >= 0 && p <= 1)
})
