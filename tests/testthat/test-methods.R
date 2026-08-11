# Accessors and show methods for BreadFit / BreadResults.

.toy_fit <- function() {
  fit_bread(.make_toy_signal_se(), .make_toy_features(), ~ group)
}

test_that("results() returns the region table with a classification column", {
  res <- results(.toy_fit())
  expect_s3_class(res, "data.frame")
  expect_true(all(c("region_id", "classification") %in% colnames(res)))
  expect_gt(nrow(res), 0L)
})

test_that("classifications() names line up with results() region ids", {
  fit <- .toy_fit()
  cls <- classifications(fit)
  res <- results(fit)
  expect_type(cls, "character")
  expect_identical(names(cls), as.character(res$region_id))
  expect_identical(unname(cls), as.character(res$classification))
  expect_true(all(cls %in% c("hypermethylated", "hypomethylated",
                             "unchanged", "inconclusive")))
})

test_that("posterior_draws() is reproducible for a fixed seed", {
  fit <- .toy_fit()
  rid <- results(fit)$region_id[1]
  a <- posterior_draws(fit, region_id = rid, n = 200L, seed = 7L)
  b <- posterior_draws(fit, region_id = rid, n = 200L, seed = 7L)
  expect_equal(a, b)
  expect_identical(colnames(a), c("region_id", "draw", "value"))
  expect_identical(nrow(a), 200L)
  expect_identical(unique(a$region_id), rid)
})

test_that("posterior_draws() differs across seeds", {
  fit <- .toy_fit()
  rid <- results(fit)$region_id[1]
  a <- posterior_draws(fit, region_id = rid, n = 200L, seed = 1L)
  b <- posterior_draws(fit, region_id = rid, n = 200L, seed = 2L)
  expect_false(isTRUE(all.equal(a$value, b$value)))
})

test_that("posterior_draws() restores the caller's RNG state", {
  # local_seed() must not leak its reseed into the calling session.
  fit <- .toy_fit()
  rid <- results(fit)$region_id[1]
  set.seed(99)
  before <- runif(1)
  set.seed(99)
  invisible(posterior_draws(fit, region_id = rid, n = 10L, seed = 123L))
  after <- runif(1)
  expect_equal(before, after)
})

test_that("posterior_draws() defaults to every region", {
  fit <- .toy_fit()
  n_regions <- nrow(results(fit))
  d <- posterior_draws(fit, n = 10L, seed = 1L)
  expect_identical(nrow(d), as.integer(n_regions * 10L))
})

test_that("an unknown region_id errors and names the offender", {
  expect_error(posterior_draws(.toy_fit(), region_id = "no_such_region"),
               "not found")
})

test_that("show() prints the expected BreadFit header", {
  fit <- .toy_fit()
  expect_output(show(fit), "<BreadFit>")
  expect_output(show(fit), "classifications:")
  expect_output(show(fit), "backend")
})

test_that("show() prints the expected BreadResults header", {
  br <- BreadResults(.toy_fit())
  expect_output(show(br), "<BreadResults>")
  expect_output(show(br), "n_regions")
})
