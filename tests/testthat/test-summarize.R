test_that("summarize returns matrix with expected shape and names", {
  se <- .make_toy_se()
  gr <- .make_toy_features()
  mp <- suppressMessages(map_probes_to_features(se, gr, min_probes = 1L))
  mat <- summarize_features(se, mp, summary_fun = "mean", input_scale = "M")

  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), 3L)
  expect_equal(ncol(mat), ncol(se))
  expect_setequal(rownames(mat), c("regA","regB","regC"))
  expect_equal(colnames(mat), colnames(se))
  expect_identical(attr(mat, "summary_fun"), "mean")
  expect_identical(attr(mat, "input_scale"), "M")
})

test_that("mean summary matches colMeans of the probe submatrix", {
  se <- .make_toy_se(); gr <- .make_toy_features()
  mp <- suppressMessages(map_probes_to_features(se, gr, min_probes = 1L))
  mat <- summarize_features(se, mp, "mean", "M")

  X <- SummarizedExperiment::assay(se, "M")
  for (rid in rownames(mat)) {
    idx <- unique(mp$probe_idx[mp$region_id == rid])
    expect_equal(unname(mat[rid, ]), unname(colMeans(X[idx, , drop = FALSE])))
  }
})

test_that("median summary matches per-sample median of the probe submatrix", {
  se <- .make_toy_se(); gr <- .make_toy_features()
  mp <- suppressMessages(map_probes_to_features(se, gr, min_probes = 1L))
  mat <- summarize_features(se, mp, "median", "M")

  X <- SummarizedExperiment::assay(se, "M")
  for (rid in rownames(mat)) {
    idx <- unique(mp$probe_idx[mp$region_id == rid])
    expected <- apply(X[idx, , drop = FALSE], 2L, stats::median)
    expect_equal(unname(mat[rid, ]), unname(expected))
  }
})

test_that("weighted_mean equals mean when probe variances are equal", {
  # Equal-variance probes: each probe has identical across-sample variance
  set.seed(42)
  n_probes <- 4L; n_samples <- 6L
  probes <- GenomicRanges::GRanges("chr1",
    IRanges::IRanges(start = seq(1, by = 100, length.out = n_probes), width = 1))
  names(probes) <- paste0("cg", seq_len(n_probes))
  base <- rnorm(n_samples)  # same variance pattern replicated
  X <- t(sapply(seq_len(n_probes), function(i) base + rnorm(1, sd = 0.01)))
  colnames(X) <- paste0("S", seq_len(n_samples))
  rownames(X) <- names(probes)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(M = X), rowRanges = probes,
    colData = S4Vectors::DataFrame(group = factor(rep(c("a","b"), length.out = n_samples)),
                                   row.names = colnames(X))
  )
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 10000),
                               feature_class = "X")
  names(gr) <- "r1"
  mp <- suppressMessages(map_probes_to_features(se, gr, min_probes = 1L))
  m_mean <- summarize_features(se, mp, "mean",          "M")
  m_wtd  <- summarize_features(se, mp, "weighted_mean", "M")
  expect_equal(as.numeric(m_wtd), as.numeric(m_mean), tolerance = 1e-10)
})

test_that("pc1 sign-aligns with per-sample mean across probes", {
  se <- .make_toy_se(); gr <- .make_toy_features()
  mp <- suppressMessages(map_probes_to_features(se, gr, min_probes = 1L))
  mat_mean <- summarize_features(se, mp, "mean", "M")
  mat_pc1  <- summarize_features(se, mp, "pc1",  "M")
  for (rid in rownames(mat_pc1)) {
    # pc1 should correlate non-negatively with the plain mean
    expect_gte(stats::cor(mat_pc1[rid, ], mat_mean[rid, ]), 0)
  }
})

test_that("single-probe region collapses to that probe's values (all summaries)", {
  se <- .make_toy_se()
  gr <- GenomicRanges::GRanges("chr1",
    IRanges::IRanges(start = 1L, end = 10L), feature_class = "S")
  names(gr) <- "rsolo"
  mp <- suppressMessages(map_probes_to_features(se, gr, min_probes = 1L))
  expect_equal(nrow(mp), 1L)
  X <- SummarizedExperiment::assay(se, "M")
  ref <- X[mp$probe_idx, ]
  for (sf in c("mean","median","weighted_mean","pc1")) {
    mat <- summarize_features(se, mp, sf, "M")
    expect_equal(unname(mat["rsolo", ]), unname(ref), info = sf)
  }
})

test_that("Beta input is converted to M before summarization", {
  # Build SE whose assay is beta-valued; compare to pre-converted M summary
  se <- .make_toy_se()
  X  <- SummarizedExperiment::assay(se, "M")
  B  <- 2^X / (2^X + 1)  # M -> Beta
  se_beta <- SummarizedExperiment::SummarizedExperiment(
    assays = list(M = B),
    rowRanges = SummarizedExperiment::rowRanges(se),
    colData = SummarizedExperiment::colData(se)
  )
  gr <- .make_toy_features()
  mp <- suppressMessages(map_probes_to_features(se_beta, gr, min_probes = 1L))

  m_from_beta <- summarize_features(se_beta, mp, "mean", "Beta")
  m_from_m    <- summarize_features(se,      mp, "mean", "M")
  expect_equal(as.numeric(m_from_beta), as.numeric(m_from_m), tolerance = 1e-10)
  expect_identical(attr(m_from_beta, "input_scale"), "Beta")
})

test_that("NAs in a probe row are tolerated by mean/median/weighted_mean", {
  se <- .make_toy_se()
  X  <- SummarizedExperiment::assay(se, "M")
  X[1L, 1L] <- NA_real_
  SummarizedExperiment::assay(se, "M") <- X
  gr <- .make_toy_features()
  mp <- suppressMessages(map_probes_to_features(se, gr, min_probes = 1L))
  for (sf in c("mean","median","weighted_mean")) {
    mat <- summarize_features(se, mp, sf, "M")
    expect_false(anyNA(mat), info = sf)
  }
})

test_that("mapping missing required columns errors", {
  se <- .make_toy_se()
  bad <- data.frame(probe_id = "x", region_id = "r")  # no probe_idx
  expect_error(summarize_features(se, bad, "mean", "M"),
               "missing required columns")
})
