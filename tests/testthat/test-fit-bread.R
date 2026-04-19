test_that("fit_bread end-to-end produces a populated BreadFit", {
  se  <- .make_toy_signal_se()
  gr  <- .make_toy_features()
  fit <- suppressMessages(
    fit_bread(se, gr, design = ~ group, contrast = "groupold",
              delta = 0.10, prob_cutoff = 0.95, min_probes = 3L)
  )
  expect_s4_class(fit, "BreadFit")
  expect_equal(fit@mode,        "summary")
  expect_equal(fit@input_scale, "M")
  expect_equal(fit@params$contrast, "groupold")
  expect_equal(fit@diagnostics$backend, "conjugate")
  expect_equal(fit@diagnostics$n_features_in,  3L)
  expect_equal(fit@diagnostics$n_features_out, 2L)
  expect_true("regB" %in% fit@diagnostics$dropped_regions)
  expect_equal(fit@diagnostics$n_failed_fits, 0L)
  expect_equal(length(fit@features), 2L)
  expect_s3_class(fit@results, "data.frame")
  expect_true("classification" %in% colnames(fit@results))
})

test_that("fit_bread recovers injected hyper/hypo signal", {
  se <- .make_toy_signal_se()
  gr <- .make_toy_features()
  fit <- suppressMessages(
    fit_bread(se, gr, ~ group, "groupold",
              delta = 0.10, prob_cutoff = 0.95, min_probes = 3L)
  )
  cls <- classifications(fit)
  expect_equal(cls[["regA"]], "hypermethylated")
  expect_equal(cls[["regC"]], "hypomethylated")
})

test_that("fit_bread defaults contrast with a message", {
  se <- .make_toy_signal_se(); gr <- .make_toy_features()
  msgs <- testthat::capture_messages(
    fit <- fit_bread(se, gr, design = ~ group, contrast = NULL,
                     min_probes = 3L)
  )
  expect_match(paste(msgs, collapse = " "),
               "first non-intercept coefficient")
  expect_equal(fit@params$contrast, "groupold")
})

test_that("results() and classifications() accessors return expected objects", {
  se <- .make_toy_signal_se(); gr <- .make_toy_features()
  fit <- suppressMessages(fit_bread(se, gr, ~ group, "groupold",
                                    min_probes = 3L))
  r <- results(fit)
  expect_s3_class(r, "data.frame")
  expect_true(all(c("region_id","classification","p_gt_delta","p_lt_neg_delta")
                  %in% colnames(r)))
  cls <- classifications(fit)
  expect_type(cls, "character")
  expect_setequal(names(cls), c("regA","regC"))
})

test_that("posterior_draws samples from scaled-t around recovered effect", {
  se <- .make_toy_signal_se(); gr <- .make_toy_features()
  fit <- suppressMessages(fit_bread(se, gr, ~ group, "groupold",
                                    min_probes = 3L))
  d <- posterior_draws(fit, region_id = "regA", n = 2000L, seed = 1L)
  expect_s3_class(d, "data.frame")
  expect_equal(nrow(d), 2000L)
  expect_named(d, c("region_id","draw","value"))
  # Injected hyper signal should push draws mean clearly above 0
  expect_gt(mean(d$value), 0.5)
  # All regions
  d_all <- posterior_draws(fit, n = 500L, seed = 1L)
  expect_setequal(unique(d_all$region_id), c("regA","regC"))
})

test_that("posterior_draws errors on unknown region_id", {
  se <- .make_toy_signal_se(); gr <- .make_toy_features()
  fit <- suppressMessages(fit_bread(se, gr, ~ group, "groupold",
                                    min_probes = 3L))
  expect_error(posterior_draws(fit, region_id = "nope"),
               "region_id\\(s\\) not found")
})

test_that("fit_bread guards hierarchical mode and brms backend", {
  se <- .make_toy_signal_se(); gr <- .make_toy_features()
  expect_error(
    fit_bread(se, gr, ~ group, "groupold", mode = "hierarchical"),
    "not yet implemented"
  )
  expect_error(
    fit_bread(se, gr, ~ group, "groupold", backend = "brms"),
    "planned for v2"
  )
  expect_error(
    fit_bread(se, gr, ~ group, "groupold", backend = "cmdstanr"),
    "planned for v2"
  )
})

test_that("fit_bread rejects unknown feature_class_col", {
  se <- .make_toy_signal_se(); gr <- .make_toy_features()
  expect_error(
    suppressMessages(fit_bread(se, gr, ~ group, "groupold",
                               feature_class_col = "not_a_col",
                               min_probes = 3L)),
    "not found in"
  )
})

test_that("fit_bread errors when all regions drop below min_probes", {
  se <- .make_toy_signal_se(); gr <- .make_toy_features()
  expect_error(
    suppressMessages(fit_bread(se, gr, ~ group, "groupold",
                               min_probes = 100L)),
    "No regions survived"
  )
})

test_that("BreadFit show() prints classification counts", {
  se <- .make_toy_signal_se(); gr <- .make_toy_features()
  fit <- suppressMessages(fit_bread(se, gr, ~ group, "groupold",
                                    min_probes = 3L))
  out <- capture.output(show(fit))
  expect_true(any(grepl("classifications:", out)))
  expect_true(any(grepl("hypermethylated", out)))
  expect_true(any(grepl("hypomethylated",  out)))
})
