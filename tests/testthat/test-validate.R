test_that("validate passes with clean inputs", {
  se <- .make_toy_se()
  gr <- .make_toy_features()
  expect_invisible(
    expect_true(
      validate_bread_input(se, gr, design = ~ group, contrast = "groupold",
                           assay_name = "M", input_scale = "M")
    )
  )
})

test_that("validate accepts NULL contrast", {
  se <- .make_toy_se()
  gr <- .make_toy_features()
  expect_true(
    validate_bread_input(se, gr, design = ~ group + sex, contrast = NULL,
                         assay_name = "M", input_scale = "M")
  )
})

test_that("validate rejects non-SE", {
  expect_error(
    validate_bread_input(se = list(), features = .make_toy_features(),
                         design = ~ group, assay_name = "M"),
    "SummarizedExperiment"
  )
})

test_that("validate rejects non-GRanges features", {
  se <- .make_toy_se()
  expect_error(
    validate_bread_input(se, features = data.frame(), design = ~ group,
                         assay_name = "M"),
    "GRanges"
  )
})

test_that("validate rejects two-sided design", {
  se <- .make_toy_se(); gr <- .make_toy_features()
  expect_error(
    validate_bread_input(se, gr, design = y ~ group, assay_name = "M"),
    "one-sided"
  )
})

test_that("validate rejects missing assay", {
  se <- .make_toy_se(); gr <- .make_toy_features()
  expect_error(
    validate_bread_input(se, gr, design = ~ group, assay_name = "beta"),
    "not found in se"
  )
})

test_that("validate rejects design variables missing from colData", {
  se <- .make_toy_se(); gr <- .make_toy_features()
  expect_error(
    validate_bread_input(se, gr, design = ~ treatment, assay_name = "M"),
    "not found in .colData.se.."
  )
})

test_that("validate rejects unknown contrast", {
  se <- .make_toy_se(); gr <- .make_toy_features()
  expect_error(
    validate_bread_input(se, gr, design = ~ group, contrast = "groupYOUNG",
                         assay_name = "M"),
    "not found among design coefficients"
  )
})

test_that("validate rejects empty features", {
  se <- .make_toy_se()
  gr <- .make_toy_features()[0]
  expect_error(
    validate_bread_input(se, gr, design = ~ group, assay_name = "M"),
    "empty"
  )
})
