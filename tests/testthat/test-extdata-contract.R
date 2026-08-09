# Every @examples block in this package is written against the two packaged
# extdata objects. If their shape drifts -- an assay rename, a dropped
# colData column, a change in the feature classes -- the examples break at
# R CMD check time with an opaque error. These tests pin the contract so the
# failure lands here instead, with a message that says what changed.

.se  <- function() {
  readRDS(system.file("extdata", "vitc_ag06561.rds", package = "BREAD"))
}
.reg <- function() {
  readRDS(system.file("extdata", "vitc_regions.rds", package = "BREAD"))
}

test_that("both extdata files are installed and loadable", {
  expect_true(nzchar(system.file("extdata", "vitc_ag06561.rds",
                                 package = "BREAD")))
  expect_true(nzchar(system.file("extdata", "vitc_regions.rds",
                                 package = "BREAD")))
})

test_that("the packaged SE has the assay name and scale the examples assume", {
  se <- .se()
  expect_s4_class(se, "RangedSummarizedExperiment")
  expect_identical(SummarizedExperiment::assayNames(se), "betas")

  x <- SummarizedExperiment::assay(se, "betas")
  expect_true(all(x >= 0 & x <= 1, na.rm = TRUE))
  # Examples rely on fit_bread() auto-detecting both of these.
  expect_identical(BREAD:::.detect_input_scale(x), "Beta")
  expect_identical(BREAD:::.detect_assay_name(se), "betas")
})

test_that("the colData columns the examples subset on are present", {
  cd <- SummarizedExperiment::colData(.se())

  expect_true(all(c("condition", "passage") %in% colnames(cd)))
  expect_s3_class(cd$condition, "factor")
  expect_s3_class(cd$passage, "factor")
  expect_identical(levels(cd$condition), c("ctrl", "aa57"))
  expect_identical(levels(cd$passage), c("early", "late"))

  # `se[, se$condition == "ctrl"]` must leave both passage levels populated,
  # otherwise `~ passage` is not estimable.
  ctrl <- cd[cd$condition == "ctrl", , drop = FALSE]
  expect_gt(nrow(ctrl), 0L)
  expect_setequal(as.character(unique(ctrl$passage)), c("early", "late"))
})

test_that("the packaged regions carry the feature_class column", {
  reg <- .reg()
  expect_s4_class(reg, "GRanges")
  expect_gt(length(reg), 0L)
  expect_true("feature_class" %in% colnames(S4Vectors::mcols(reg)))
  expect_false(is.null(names(reg)))
  expect_identical(anyDuplicated(names(reg)), 0L)
})

test_that("regions and probes overlap at the example threshold", {
  # If this fails, every example calling fit_bread() on the packaged data
  # errors with "no regions retained".
  mapping <- map_probes_to_features(.se(), .reg(), min_probes = 3L)
  expect_gt(nrow(mapping), 0L)
  expect_gt(length(unique(mapping$region_id)), 0L)
  expect_true("feature_class" %in% colnames(mapping))
})

test_that("the documented example fit runs end to end", {
  se  <- .se()
  reg <- .reg()
  se_ctrl <- se[, se$condition == "ctrl"]

  fit <- fit_bread(se_ctrl, reg, ~ passage,
                   feature_class_col = "feature_class")
  expect_s4_class(fit, "BreadFit")

  res <- results(fit)
  expect_gt(nrow(res), 0L)
  expect_true(all(c("region_id", "classification") %in% colnames(res)))
  expect_true(any(!is.na(res$p_gt_delta)))
})
