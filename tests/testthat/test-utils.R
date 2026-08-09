# Internal transforms, palettes, and the auto-detection helpers that make
# fit_bread()'s three-argument form work.

test_that("beta <-> M round-trips across the usable range", {
  betas <- c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)
  expect_equal(BREAD:::.m_to_beta(BREAD:::.beta_to_m(betas)), betas,
               tolerance = 1e-8)
})

test_that("beta -> M clamps at the boundaries instead of returning Inf", {
  out <- BREAD:::.beta_to_m(c(0, 1))
  expect_true(all(is.finite(out)))
  expect_lt(out[1], 0)
  expect_gt(out[2], 0)
})

test_that("beta -> M is monotone increasing and centred at 0.5", {
  m <- BREAD:::.beta_to_m(seq(0.05, 0.95, by = 0.05))
  expect_true(all(diff(m) > 0))
  expect_equal(BREAD:::.beta_to_m(0.5), 0, tolerance = 1e-6)
})

test_that("input scale detection keys off the [0, 1] range", {
  expect_identical(BREAD:::.detect_input_scale(c(0, 0.5, 1)),  "Beta")
  expect_identical(BREAD:::.detect_input_scale(c(0.2, 0.8)),   "Beta")
  expect_identical(BREAD:::.detect_input_scale(c(-2, 0.5, 3)), "M")
  expect_identical(BREAD:::.detect_input_scale(c(0.5, 1.5)),   "M")
})

test_that("assay-name detection follows the documented priority", {
  mk <- function(nms) {
    m <- matrix(0.5, nrow = 2, ncol = 2,
                dimnames = list(c("p1", "p2"), c("s1", "s2")))
    a <- stats::setNames(replicate(length(nms), m, simplify = FALSE), nms)
    SummarizedExperiment::SummarizedExperiment(assays = a)
  }
  expect_identical(BREAD:::.detect_assay_name(mk(c("betas", "M"))),    "M")
  expect_identical(BREAD:::.detect_assay_name(mk(c("Beta", "betas"))), "betas")
  expect_identical(BREAD:::.detect_assay_name(mk(c("beta", "Beta"))),  "Beta")
  expect_identical(BREAD:::.detect_assay_name(mk("beta")),             "beta")
  expect_identical(BREAD:::.detect_assay_name(mk(c("weird", "other"))), "weird")
})

test_that("bread_colors() returns the documented shapes", {
  cl <- bread_colors("classification")
  expect_length(cl, 3L)
  expect_identical(names(cl),
                   c("hypermethylated", "hypomethylated", "inconclusive"))

  grp <- bread_colors("group")
  expect_length(grp, 2L)
  expect_null(names(grp))

  cross <- bread_colors("cross")
  expect_length(cross, 9L)

  for (pal in list(cl, grp, cross)) {
    expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", pal)))
  }
})

test_that("bread_colors() defaults to the classification palette", {
  expect_identical(bread_colors(), bread_colors("classification"))
  expect_error(bread_colors("nonsense"))
})
