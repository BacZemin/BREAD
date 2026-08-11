# The S4 slot layout is a public contract: the plotting helpers and the
# accessors reach into these slots by name. Pin them.

test_that("BreadFit slots are the documented set, in order", {
  expect_identical(
    methods::slotNames("BreadFit"),
    c("call", "params", "mode", "assay_name", "input_scale",
      "mapping", "features", "model", "posterior", "results", "diagnostics")
  )
})

test_that("BreadResults slots are the documented set", {
  expect_identical(methods::slotNames("BreadResults"), c("table", "params"))
})

test_that("a fitted BreadFit populates every slot it promises", {
  fit <- fit_bread(.make_toy_signal_se(), .make_toy_features(), ~ group)
  expect_s4_class(fit, "BreadFit")
  expect_true(is.call(fit@call))
  expect_identical(fit@mode, "summary")
  expect_type(fit@params, "list")
  expect_type(fit@diagnostics, "list")
  expect_s3_class(fit@mapping, "data.frame")
  expect_s3_class(fit@results, "data.frame")
  expect_s3_class(fit@posterior, "data.frame")
  expect_type(fit@model, "list")
  expect_identical(fit@diagnostics$backend, "conjugate")
})

test_that("BreadResults() round-trips the results table", {
  fit <- fit_bread(.make_toy_signal_se(), .make_toy_features(), ~ group)
  br  <- BreadResults(fit)
  expect_s4_class(br, "BreadResults")
  expect_equal(methods::slot(br, "table"), results(fit))
  expect_equal(methods::slot(br, "params"), fit@params)
})

test_that("BreadResults() rejects non-BreadFit input", {
  expect_error(BreadResults(data.frame(x = 1)), "must be a BreadFit")
  expect_error(BreadResults(NULL), "must be a BreadFit")
})

test_that("the model list shape is a stable contract", {
  # refit_bread(), posterior_summary() and plot_region_data() all reach into
  # these names. Pin them so a backend refactor cannot quietly break them.
  fit <- suppressMessages(
    fit_bread(.make_toy_signal_se(), .make_toy_features(), ~ group)
  )
  expect_identical(
    names(fit@model),
    c("fits", "design_matrix", "coef_names", "contrast", "contrast_idx",
      "region_ids", "prior", "df_mode", "region_mat", "design", "coldata")
  )
  expect_identical(fit@model$df_mode, "conjugate")
  expect_true(is.matrix(fit@model$region_mat))
  expect_identical(rownames(fit@model$region_mat), fit@model$region_ids)
})
