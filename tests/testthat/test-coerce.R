# Matrix input has to reach exactly the same answer as SummarizedExperiment
# input -- otherwise the sesame workflow is a second, subtly different code
# path rather than a convenience.

.toy_parts <- function() {
  se <- .make_toy_signal_se()
  list(
    se  = se,
    mat = SummarizedExperiment::assay(se, "M"),
    cd  = as.data.frame(SummarizedExperiment::colData(se)),
    gr  = SummarizedExperiment::rowRanges(se)
  )
}

test_that("matrix input and SE input give identical results", {
  p <- .toy_parts()
  gr_feat <- .make_toy_features()

  fit_se <- suppressMessages(fit_bread(p$se, gr_feat, ~ group))
  fit_m  <- suppressMessages(
    fit_bread(p$mat, gr_feat, ~ group, colData = p$cd, rowRanges = p$gr)
  )

  expect_equal(results(fit_m),   results(fit_se))
  expect_equal(fit_m@posterior,  fit_se@posterior)
  expect_equal(fit_m@mapping,    fit_se@mapping)
  expect_identical(fit_m@input_scale, fit_se@input_scale)
})

test_that("bread_se() round-trips a decomposed SummarizedExperiment", {
  p  <- .toy_parts()
  se2 <- bread_se(p$mat, colData = p$cd, rowRanges = p$gr)
  expect_s4_class(se2, "SummarizedExperiment")
  expect_equal(SummarizedExperiment::assay(se2, "M"), p$mat)
  expect_equal(nrow(SummarizedExperiment::colData(se2)), ncol(p$mat))
  expect_identical(names(SummarizedExperiment::rowRanges(se2)), rownames(p$mat))
})

test_that("a sesameData-style list is unwrapped", {
  p  <- .toy_parts()
  se2 <- bread_se(list(betas = p$mat, sampleInfo = p$cd), rowRanges = p$gr)
  expect_s4_class(se2, "SummarizedExperiment")
  expect_true("group" %in% colnames(SummarizedExperiment::colData(se2)))
})

test_that("the assay is named from the value range, as fit_bread expects", {
  p <- .toy_parts()
  se_m <- bread_se(p$mat, colData = p$cd, rowRanges = p$gr)
  expect_identical(SummarizedExperiment::assayNames(se_m), "M")

  betas <- 2^p$mat / (2^p$mat + 1)
  se_b  <- bread_se(betas, colData = p$cd, rowRanges = p$gr)
  expect_identical(SummarizedExperiment::assayNames(se_b), "betas")
})


# ---- rowRanges alignment ---------------------------------------------------

test_that("a named manifest is subset and reordered to the matrix rows", {
  p <- .toy_parts()
  shuffled <- p$gr[sample(length(p$gr))]
  extra    <- p$gr[1:3]
  names(extra) <- paste0("zz", 1:3)
  manifest <- c(shuffled, extra)          # longer, out of order

  se2 <- bread_se(p$mat, colData = p$cd, rowRanges = manifest)
  expect_identical(names(SummarizedExperiment::rowRanges(se2)), rownames(p$mat))
  expect_equal(SummarizedExperiment::rowRanges(se2), p$gr)
})

test_that("probes absent from the manifest are dropped with a message", {
  p <- .toy_parts()
  partial <- p$gr[1:15]                    # 5 probes have no coordinates
  expect_message(
    se2 <- bread_se(p$mat, colData = p$cd, rowRanges = partial),
    "Dropping 5 of 20 probes"
  )
  expect_equal(nrow(se2), 15L)
})

test_that("an entirely mismatched manifest errors rather than dropping all", {
  p <- .toy_parts()
  wrong <- p$gr
  names(wrong) <- paste0("nope", seq_along(wrong))
  expect_error(bread_se(p$mat, colData = p$cd, rowRanges = wrong),
               "None of `rownames\\(x\\)`")
})

test_that("an unnamed rowRanges must match the row count exactly", {
  p <- .toy_parts()
  unnamed <- p$gr; names(unnamed) <- NULL
  se2 <- bread_se(p$mat, colData = p$cd, rowRanges = unnamed)
  expect_identical(names(SummarizedExperiment::rowRanges(se2)), rownames(p$mat))

  expect_error(bread_se(p$mat, colData = p$cd, rowRanges = unnamed[1:5]),
               "one range per row")
})


# ---- colData alignment -----------------------------------------------------

test_that("colData is reordered by rownames, not trusted positionally", {
  p <- .toy_parts()
  shuffled <- p$cd[sample(nrow(p$cd)), , drop = FALSE]

  ordered  <- bread_se(p$mat, colData = p$cd,     rowRanges = p$gr)
  reorder  <- bread_se(p$mat, colData = shuffled, rowRanges = p$gr)
  expect_equal(SummarizedExperiment::colData(reorder),
               SummarizedExperiment::colData(ordered))
})

test_that("colData that does not cover every sample errors", {
  p <- .toy_parts()
  expect_error(
    bread_se(p$mat, colData = p$cd[1:4, , drop = FALSE], rowRanges = p$gr),
    "no row for"
  )
})

test_that("colData without rownames is accepted but warns", {
  p <- .toy_parts()
  bare <- p$cd; rownames(bare) <- NULL
  expect_warning(bread_se(p$mat, colData = bare, rowRanges = p$gr),
                 "assuming its rows are in the same order")
})


# ---- refusals --------------------------------------------------------------

test_that("coordinates are never guessed from probe IDs", {
  p <- .toy_parts()
  expect_error(bread_se(p$mat, colData = p$cd), "does not guess the platform")
})

test_that("matrix-only arguments are rejected alongside an SE", {
  p <- .toy_parts()
  expect_error(bread_se(p$se, colData = p$cd),   "supplied alongside")
  expect_error(bread_se(p$se, rowRanges = p$gr), "supplied alongside")
  expect_error(bread_se(p$se, platform = "EPIC"),"supplied alongside")
  expect_s4_class(bread_se(p$se), "SummarizedExperiment")
})

test_that("a matrix without dimnames is rejected", {
  p <- .toy_parts()
  m <- p$mat; rownames(m) <- NULL
  expect_error(bread_se(m, colData = p$cd, rowRanges = p$gr), "rownames")

  m2 <- p$mat; colnames(m2) <- NULL
  expect_error(bread_se(m2, colData = p$cd, rowRanges = p$gr), "colnames")
})

test_that("unsupported input still names SummarizedExperiment in the error", {
  # test-smoke.R relies on this string
  expect_error(bread_se(NULL), "SummarizedExperiment")
  expect_error(fit_bread(NULL, NULL, ~ 1), "SummarizedExperiment")
})

test_that("the platform route reaches sesameData", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("sesameData")
  # Probe IDs are taken from the manifest itself, so this tests the lookup
  # rather than whether some other dataset happens to share its ID
  # convention. (The packaged vitc subset does not: it carries stripped
  # EPICv2 IDs while the manifest keeps the replicate suffix.)
  man <- sesameData::sesameData_getManifestGRanges("EPICv2")
  skip_if(length(man) == 0L, "EPICv2 manifest unavailable offline")

  ids <- names(man)[seq_len(50L)]
  mat <- matrix(stats::runif(50L * 4L), nrow = 50L,
                dimnames = list(ids, sprintf("S%d", 1:4)))
  cd  <- data.frame(group = rep(c("a", "b"), 2), row.names = colnames(mat))

  se2 <- suppressMessages(bread_se(mat, colData = cd, platform = "EPICv2"))
  expect_s4_class(se2, "SummarizedExperiment")
  expect_equal(nrow(se2), 50L)
  expect_identical(names(SummarizedExperiment::rowRanges(se2)), ids)
  expect_identical(SummarizedExperiment::assayNames(se2), "betas")
})

test_that("stripped EPICv2 suffixes get a specific diagnosis", {
  p <- .toy_parts()
  manifest <- p$gr
  names(manifest) <- paste0(names(manifest), "_BC11")
  expect_error(bread_se(p$mat, colData = p$cd, rowRanges = manifest),
               "EPICv2 replicate suffixes")
})
