test_that("package loads and exported symbols exist", {
  expect_true(exists("fit_bread", mode = "function"))
  expect_true(isClass("BreadFit"))
  expect_true(isClass("BreadResults"))
})

test_that("utility transforms round-trip", {
  expect_equal(BREAD:::.m_to_beta(BREAD:::.beta_to_m(0.3)), 0.3, tolerance = 1e-6)
  expect_equal(BREAD:::.m_to_beta(BREAD:::.beta_to_m(0.8)), 0.8, tolerance = 1e-6)
})

test_that("fit_bread() errors with scaffold-only message", {
  expect_error(
    fit_bread(se = NULL, features = NULL, design = ~1),
    "not yet implemented"
  )
})
