test_that("mapping returns expected columns and probe counts per region", {
  se <- .make_toy_se(); gr <- .make_toy_features()
  m  <- suppressMessages(map_probes_to_features(se, gr, min_probes = 1L))
  expect_s3_class(m, "data.frame")
  expect_named(m,
    c("probe_id","probe_idx","region_id","region_idx","feature_class","n_probes"),
    ignore.order = TRUE)

  counts <- tapply(m$probe_id, m$region_id, length)
  # Region boundaries: regA 1-5500 (probes @ 1,1001,2001,3001,4001,5001 = 6);
  # regB 10001-11500 (probes @ 10001, 11001 = 2); regC 15001-20000 (probes @
  # 15001..19001 = 5).
  expect_equal(counts[["regA"]], 6L)
  expect_equal(counts[["regB"]], 2L)
  expect_equal(counts[["regC"]], 5L)
})

test_that("mapping drops regions below min_probes and records them", {
  se <- .make_toy_se(); gr <- .make_toy_features()
  m  <- suppressMessages(map_probes_to_features(se, gr, min_probes = 3L))
  expect_false("regB" %in% m$region_id)
  expect_equal(attr(m, "dropped_regions"), "regB")
  expect_equal(attr(m, "n_features_in"),  3L)
  expect_equal(attr(m, "n_features_out"), 2L)
  expect_equal(attr(m, "min_probes"),     3L)
})

test_that("mapping preserves mcols feature metadata", {
  se <- .make_toy_se(); gr <- .make_toy_features()
  m  <- suppressMessages(map_probes_to_features(se, gr, min_probes = 1L))
  expect_true("feature_class" %in% colnames(m))
  # regA -> PRC on all its rows
  expect_true(all(m$feature_class[m$region_id == "regA"] == "PRC"))
})

test_that("mapping generates region IDs when features are unnamed", {
  se <- .make_toy_se(); gr <- .make_toy_features()
  names(gr) <- NULL
  m  <- suppressMessages(map_probes_to_features(se, gr, min_probes = 1L))
  expect_true(all(grepl("^region_[0-9]+$", unique(m$region_id))))
})

test_that("mapping errors on zero overlaps", {
  se <- .make_toy_se()
  gr <- GenomicRanges::GRanges("chr99", IRanges::IRanges(1, 100))
  expect_error(map_probes_to_features(se, gr, min_probes = 1L),
               "No probes overlap")
})

test_that("mapping rejects non-positive min_probes", {
  se <- .make_toy_se(); gr <- .make_toy_features()
  expect_error(map_probes_to_features(se, gr, min_probes = 0L),
               "positive integer")
})
