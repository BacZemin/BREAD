# Degrees-of-freedom handling: the n <= p guard and the "residual" mode.
#
# Motivation: `a_n = a0 + n/2` makes nu = 2*a_n a function of n alone, never of
# p. Under the weak default prior that overstates precision by exactly
# sqrt(n / (n - p)) on the posterior scale, and at n == p the residuals are
# identically zero so the scale collapses to the prior floor.

# --- helpers ---------------------------------------------------------------

# Region x sample matrix with a known design; y is pure noise unless `beta` set.
sim_mat <- function(n_per_cell, n_regions = 5L, sd = 1, beta = 0, seed = 1L) {
  set.seed(seed)
  g <- rep(c("a", "b"), each = 2L * n_per_cell)
  t <- rep(rep(c("x", "y"), each = n_per_cell), 2L)
  cd <- data.frame(g = factor(g), t = factor(t))
  X  <- stats::model.matrix(~ g * t, cd)
  k  <- which(colnames(X) == "gb:ty")
  eta <- as.numeric(X[, k]) * beta
  m <- matrix(rep(eta, each = n_regions) + stats::rnorm(n_regions * nrow(cd), sd = sd),
              nrow = n_regions, dimnames = list(paste0("R", seq_len(n_regions)), rownames(cd)))
  list(mat = m, cd = cd, p = ncol(X), contrast = "gb:ty")
}

fit_one <- function(s, df_mode = "conjugate", prior = NULL) {
  BREAD:::fit_bread_summary(s$mat, s$cd, ~ g * t, s$contrast,
                            prior = prior, df_mode = df_mode)
}

# --- the n <= p guard ------------------------------------------------------

test_that("regions with n == p are dropped rather than fitted", {
  s <- sim_mat(n_per_cell = 1L)          # n = 4, p = 4
  expect_identical(s$p, 4L)
  f <- suppressWarnings(fit_one(s))
  errs <- vapply(f$fits, function(z) z$error, character(1))
  expect_true(all(errs == "n <= number of coefficients"))
  expect_true(all(vapply(f$fits, function(z) is.na(z$a_n), logical(1))))
})

test_that("the n <= p guard applies under both df_mode settings", {
  s <- sim_mat(n_per_cell = 1L)
  for (dm in c("conjugate", "residual")) {
    f <- suppressWarnings(fit_one(s, df_mode = dm))
    expect_true(all(vapply(f$fits, function(z) z$error, character(1)) ==
                      "n <= number of coefficients"))
  }
})

test_that("n < 2 still reports the pre-existing reason", {
  s <- sim_mat(n_per_cell = 2L)
  s$mat[1, ] <- NA_real_
  s$mat[1, 1] <- 0.5                      # a single non-NA sample
  f <- suppressWarnings(fit_one(s))
  expect_identical(f$fits[[1]]$error, "too few non-NA samples")
})

test_that("n > p fits normally and carries no error", {
  s <- sim_mat(n_per_cell = 3L)          # n = 12, p = 4
  f <- fit_one(s)
  expect_true(all(is.na(vapply(f$fits, function(z) z$error, character(1)))))
})

# --- low residual df warns once, not per region ---------------------------

test_that("fewer than 3 residual df warns exactly once for the whole fit", {
  s <- sim_mat(n_per_cell = 2L, n_regions = 10L)   # n = 8, p = 4 -> n - p = 4
  expect_silent(fit_one(s))

  s2 <- sim_mat(n_per_cell = 2L, n_regions = 10L)
  s2$cd$z <- factor(rep(c("u", "v"), length.out = nrow(s2$cd)))
  # ~ g * t + z  -> p = 5, n = 8, n - p = 3 -> still silent
  f5 <- BREAD:::fit_bread_summary(s2$mat, s2$cd, ~ g * t + z, "gb:ty")
  expect_true(is.list(f5$fits))

  s3 <- sim_mat(n_per_cell = 2L, n_regions = 10L)
  s3$mat[, 1:3] <- NA_real_              # n drops to 5, p = 4 -> n - p = 1
  w <- capture_warnings(fit_one(s3))
  expect_length(w, 1L)
  expect_match(w, "residual degrees of freedom", fixed = FALSE)
})

test_that("the warning names df_mode = residual only in conjugate mode", {
  s <- sim_mat(n_per_cell = 2L, n_regions = 4L)
  s$mat[, 1:3] <- NA_real_
  expect_match(capture_warnings(fit_one(s, df_mode = "conjugate")),
               "df_mode")
  expect_false(any(grepl("df_mode",
                         capture_warnings(fit_one(s, df_mode = "residual")))))
})

# --- residual mode reproduces the classical answer ------------------------

test_that("df_mode = 'residual' matches lm() df and standard error", {
  s <- sim_mat(n_per_cell = 4L, n_regions = 3L, seed = 7L)   # n = 16, p = 4
  # weak prior so the posterior should collapse onto OLS
  pr <- bread_prior(lambda0 = 1e-8, a0 = 1e-8, b0 = 1e-8)
  f  <- fit_one(s, df_mode = "residual", prior = pr)
  k  <- f$contrast_idx

  for (i in seq_len(nrow(s$mat))) {
    ml <- stats::lm(s$mat[i, ] ~ g * t, data = s$cd)
    fo <- f$fits[[i]]
    scale_b <- sqrt((fo$b_n / fo$a_n) * fo$Lambda_n_inv[k, k])
    expect_equal(2 * fo$a_n, ml$df.residual, tolerance = 1e-5)
    expect_equal(scale_b, summary(ml)$coefficients[k, 2], tolerance = 1e-4)
    expect_equal(fo$mu_n[k], unname(coef(ml)[k]), tolerance = 1e-5)
  }
})

test_that("conjugate vs residual differ by exactly sqrt(n / (n - p))", {
  s  <- sim_mat(n_per_cell = 4L, n_regions = 3L, seed = 11L)  # n = 16, p = 4
  pr <- bread_prior(lambda0 = 1e-8, a0 = 1e-8, b0 = 1e-8)
  fc <- fit_one(s, df_mode = "conjugate", prior = pr)
  fr <- fit_one(s, df_mode = "residual",  prior = pr)
  k  <- fc$contrast_idx

  sc <- vapply(fc$fits, function(z) sqrt((z$b_n / z$a_n) * z$Lambda_n_inv[k, k]), 0)
  sr <- vapply(fr$fits, function(z) sqrt((z$b_n / z$a_n) * z$Lambda_n_inv[k, k]), 0)
  n <- 16L; p <- 4L
  expect_equal(unname(sr / sc), rep(sqrt(n / (n - p)), length(sc)), tolerance = 1e-5)
})

test_that("residual mode widens credible intervals", {
  s  <- sim_mat(n_per_cell = 4L, n_regions = 5L, seed = 3L)
  pr <- bread_prior(lambda0 = 1e-8, a0 = 1e-8, b0 = 1e-8)
  pc <- posterior_summary(fit_one(s, "conjugate", pr))
  prr <- posterior_summary(fit_one(s, "residual",  pr))
  expect_true(all((prr$ci_hi - prr$ci_lo) > (pc$ci_hi - pc$ci_lo)))
  expect_equal(prr$mean_effect, pc$mean_effect, tolerance = 1e-6)
})

# --- end-to-end through fit_bread() ---------------------------------------

test_that("fit_bread() accepts df_mode and records it in params", {
  skip_if_not_installed("SummarizedExperiment")
  se <- .make_toy_signal_se(); gr <- .make_toy_features()
  fit_c <- fit_bread(se, gr, ~ group)
  fit_r <- fit_bread(se, gr, ~ group, df_mode = "residual")

  expect_identical(fit_c@params$df_mode, "conjugate")
  expect_identical(fit_r@params$df_mode, "residual")

  rc <- results(fit_c); rr <- results(fit_r)
  ok <- is.na(rc$error) & is.na(rr$error)
  expect_true(any(ok))
  expect_true(all(rr$df[ok] < rc$df[ok]))
  expect_true(all(rr$scale[ok] >= rc$scale[ok]))
  expect_equal(rr$mean_effect[ok], rc$mean_effect[ok], tolerance = 1e-8)
})

test_that("df_mode is rejected when misspelled", {
  se <- .make_toy_signal_se(); gr <- .make_toy_features()
  expect_error(fit_bread(se, gr, ~ group, df_mode = "residuals"))
})
