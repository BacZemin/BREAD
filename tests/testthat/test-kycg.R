# bread_kycg() reaches out to KnowYourCG reference databases, so these tests
# deliberately cover only the validation that happens before any network or
# annotation-hub access.

.toy_fit <- function() {
  fit_bread(.make_toy_signal_se(), .make_toy_features(), ~ group)
}

test_that("a non-BreadFit is rejected before any database access", {
  expect_error(bread_kycg(data.frame(x = 1)), "BreadFit")
  expect_error(bread_kycg(NULL), "BreadFit")
})

test_that("`which` and `platform` are matched against their allowed values", {
  skip_if_not_installed("knowYourCG")
  fit <- .toy_fit()
  expect_error(bread_kycg(fit, platform = "NotAPlatform"))
  expect_error(bread_kycg(fit, which = "sideways"))
})

test_that("the fit exposes the mapping columns bread_kycg() reads", {
  # Not a network test: pins the columns bread_kycg() depends on, so a
  # refactor of map_probes_to_features() cannot silently break it.
  fit <- .toy_fit()
  expect_true(all(c("probe_id", "region_id") %in% colnames(fit@mapping)))
  expect_type(fit@mapping$probe_id, "character")
})
