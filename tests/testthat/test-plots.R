test_that("bread_colors() returns expected palettes", {
  cl <- bread_colors("classification")
  expect_named(cl, c("hypermethylated", "hypomethylated", "inconclusive"))
  expect_match(cl, "^#[0-9a-fA-F]{6}$")

  gr <- bread_colors("group")
  expect_length(gr, 2L)

  cx <- bread_colors("cross")
  expect_length(cx, 9L)
})

test_that("plot_region_posterior returns a ggplot with expected layers", {
  se  <- .make_toy_signal_se(); gr <- .make_toy_features()
  fit <- suppressMessages(fit_bread(se, gr, ~ group, "groupold",
                                    min_probes = 3L))
  p <- plot_region_posterior(fit, region_id = "regA")
  expect_s3_class(p, "ggplot")
  layers <- vapply(p$layers,
                   function(l) class(l$geom)[1L], character(1))
  expect_true(any(layers == "GeomLine"))
  expect_true(any(layers == "GeomVline"))
})

test_that("plot_region_posterior facets when multiple regions", {
  se  <- .make_toy_signal_se(); gr <- .make_toy_features()
  fit <- suppressMessages(fit_bread(se, gr, ~ group, "groupold",
                                    min_probes = 3L))
  p <- plot_region_posterior(fit)
  expect_s3_class(p, "ggplot")
  expect_s3_class(p$facet, "FacetWrap")
})

test_that("plot_region_posterior errors on unknown region", {
  se  <- .make_toy_signal_se(); gr <- .make_toy_features()
  fit <- suppressMessages(fit_bread(se, gr, ~ group, "groupold",
                                    min_probes = 3L))
  expect_error(plot_region_posterior(fit, "nope"),
               "region_id\\(s\\) not found")
})

test_that("plot_region_data returns a ggplot with boxplot + jitter", {
  se  <- .make_toy_signal_se(); gr <- .make_toy_features()
  fit <- suppressMessages(fit_bread(se, gr, ~ group, "groupold",
                                    min_probes = 3L))
  p <- plot_region_data(fit, "regA")
  expect_s3_class(p, "ggplot")
  layers <- vapply(p$layers,
                   function(l) class(l$geom)[1L], character(1))
  expect_true(any(layers == "GeomBoxplot"))
  expect_true(any(layers == "GeomPoint"))  # geom_jitter -> GeomPoint
})

test_that("plot_region_data errors on bad region", {
  se  <- .make_toy_signal_se(); gr <- .make_toy_features()
  fit <- suppressMessages(fit_bread(se, gr, ~ group, "groupold",
                                    min_probes = 3L))
  expect_error(plot_region_data(fit, "nope"),        "not found in fit")
  expect_error(plot_region_data(fit, c("regA","regC")), "single")
})

test_that("plot_feature_set returns a ggplot bar chart", {
  se  <- .make_toy_signal_se(); gr <- .make_toy_features()
  fit <- suppressMessages(fit_bread(se, gr, ~ group, "groupold",
                                    min_probes = 3L))
  p <- plot_feature_set(fit)
  expect_s3_class(p, "ggplot")
  layers <- vapply(p$layers,
                   function(l) class(l$geom)[1L], character(1))
  expect_true(any(layers == "GeomBar"))
})

test_that("plot_feature_set stacks by feature class when provided", {
  se  <- .make_toy_signal_se(); gr <- .make_toy_features()
  fit <- suppressMessages(fit_bread(se, gr, ~ group, "groupold",
                                    min_probes = 3L))
  p <- plot_feature_set(fit, feature_class_col = "feature_class")
  expect_s3_class(p, "ggplot")
  # Built data should have a feature_class column as the x aesthetic
  pb <- ggplot2::ggplot_build(p)
  expect_true(any(grepl("feature_class",
                        as.character(p$mapping),
                        fixed = FALSE)))
})

test_that("plot_feature_set errors on unknown class column", {
  se  <- .make_toy_signal_se(); gr <- .make_toy_features()
  fit <- suppressMessages(fit_bread(se, gr, ~ group, "groupold",
                                    min_probes = 3L))
  expect_error(plot_feature_set(fit, feature_class_col = "wat"),
               "not found in fit@mapping")
})
