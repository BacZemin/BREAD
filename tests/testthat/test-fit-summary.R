# Build a small SE with known ground truth per region
.sim_region_mat <- function(n_samples = 20L, true_betas, seed = 1L, sigma = 0.2) {
  set.seed(seed)
  n_regions <- length(true_betas)
  G <- rep(c(0, 1), length.out = n_samples)  # "young" = 0, "old" = 1
  mat <- t(sapply(true_betas, function(b) rnorm(n_samples, mean = b * G, sd = sigma)))
  rownames(mat) <- paste0("r", seq_len(n_regions))
  colnames(mat) <- sprintf("S%02d", seq_len(n_samples))
  coldata <- S4Vectors::DataFrame(
    group = factor(ifelse(G == 0, "young", "old"), levels = c("young","old")),
    row.names = colnames(mat)
  )
  list(mat = mat, coldata = coldata, truth = true_betas, G = G)
}

test_that("posterior mean approximates OLS under weak prior", {
  sim <- .sim_region_mat(n_samples = 30L,
                         true_betas = c(0, 0.4, -0.4), seed = 1L)
  fit <- fit_bread_summary(sim$mat, sim$coldata, design = ~ group,
                           contrast = "groupold",
                           prior = bread_prior(lambda0 = 1e-6,
                                               a0 = 1e-4, b0 = 1e-4))
  # Compare to per-region lm() estimates
  for (i in seq_len(nrow(sim$mat))) {
    df  <- data.frame(y = sim$mat[i, ], group = sim$coldata$group)
    ols <- coef(lm(y ~ group, data = df))["groupold"]
    bay <- fit$fits[[i]]$mu_n[fit$contrast_idx]
    expect_equal(unname(bay), unname(ols), tolerance = 1e-3,
                 info = paste("region", i))
  }
})

test_that("posterior_summary recovers sign and orders true effects", {
  sim <- .sim_region_mat(n_samples = 40L,
                         true_betas = c(0, 0, 0.5, 0.5, -0.5, -0.5),
                         seed = 7L)
  fit  <- fit_bread_summary(sim$mat, sim$coldata, design = ~ group,
                            contrast = "groupold")
  post <- posterior_summary(fit, delta = 0.10, ci = 0.95)
  expect_equal(nrow(post), 6L)
  expect_true(all(!is.na(post$mean_effect)))
  # Hyper regions have positive effect, hypo negative, nulls near zero
  expect_true(all(post$mean_effect[3:4]  > 0.25))
  expect_true(all(post$mean_effect[5:6]  < -0.25))
  expect_true(all(abs(post$mean_effect[1:2]) < 0.25))
  # Probs in [0,1]
  for (col in c("p_pos","p_neg","p_gt_delta","p_lt_neg_delta"))
    expect_true(all(post[[col]] >= 0 & post[[col]] <= 1))
  # p_pos + p_neg == 1 (within tolerance)
  expect_equal(post$p_pos + post$p_neg, rep(1, nrow(post)),
               tolerance = 1e-8)
})

test_that("classify_regions recovers truth with strong signal", {
  sim <- .sim_region_mat(
    n_samples  = 60L,
    true_betas = c(rep(0, 3),         # nulls
                   rep(0.8, 3),       # hyper
                   rep(-0.8, 3),      # hypo
                   rep(0.02, 3)),     # weak (below delta)
    seed = 42L, sigma = 0.2
  )
  fit  <- fit_bread_summary(sim$mat, sim$coldata,
                            design = ~ group, contrast = "groupold")
  post <- posterior_summary(fit, delta = 0.10)
  cls  <- classify_regions(post, delta = 0.10, prob_cutoff = 0.95)

  got <- as.character(cls$classification)
  expect_true(all(got[1:3]   == "inconclusive"),     info = paste(got[1:3],  collapse=","))
  expect_true(all(got[4:6]   == "hypermethylated"),  info = paste(got[4:6],  collapse=","))
  expect_true(all(got[7:9]   == "hypomethylated"),   info = paste(got[7:9],  collapse=","))
  expect_true(all(got[10:12] == "inconclusive"),     info = paste(got[10:12],collapse=","))
  expect_identical(attr(cls, "delta"), 0.10)
  expect_identical(attr(cls, "prob_cutoff"), 0.95)
})

test_that("CI coverage is near nominal over simulated regions", {
  # Many null regions at same sigma; 95% CIs should cover 0 ~95% of the time
  set.seed(2026)
  n_samples <- 30L; n_regions <- 500L
  G <- rep(c(0,1), length.out = n_samples)
  mat <- matrix(rnorm(n_regions * n_samples, sd = 0.3),
                n_regions, n_samples)
  rownames(mat) <- paste0("r", seq_len(n_regions))
  colnames(mat) <- sprintf("S%02d", seq_len(n_samples))
  cd <- S4Vectors::DataFrame(
    group = factor(ifelse(G == 0, "young", "old"), levels = c("young","old")),
    row.names = colnames(mat)
  )
  fit  <- fit_bread_summary(mat, cd, ~ group, "groupold")
  post <- posterior_summary(fit, delta = 0, ci = 0.95)
  covered <- post$ci_lo <= 0 & post$ci_hi >= 0
  expect_gte(mean(covered), 0.92)
  expect_lte(mean(covered), 0.98)
})

test_that("NA handling in region row: n < 2 flags error; summary marks NA", {
  sim <- .sim_region_mat(n_samples = 10L, true_betas = c(0.4),
                         seed = 3L)
  # Knock out all but one sample
  bad <- sim$mat
  bad[1L, 2:10] <- NA_real_
  fit  <- fit_bread_summary(bad, sim$coldata, ~ group, "groupold")
  expect_equal(fit$fits[[1]]$error, "too few non-NA samples")
  post <- posterior_summary(fit, delta = 0.10)
  expect_true(is.na(post$mean_effect))
  cls <- classify_regions(post)
  expect_equal(as.character(cls$classification), "inconclusive")
})

test_that("fit_bread_summary errors on mismatched sample count", {
  sim <- .sim_region_mat(n_samples = 10L, true_betas = c(0.4))
  cd_bad <- sim$coldata[1:5, , drop = FALSE]
  expect_error(
    fit_bread_summary(sim$mat, cd_bad, ~ group, "groupold"),
    "rows but"
  )
})

test_that("fit_bread_summary errors on unknown contrast", {
  sim <- .sim_region_mat(n_samples = 10L, true_betas = c(0.4))
  expect_error(
    fit_bread_summary(sim$mat, sim$coldata, ~ group, "groupYoung"),
    "not found among design coefficients"
  )
})

test_that("bread_prior() validates inputs", {
  expect_s3_class(bread_prior(), "bread_prior")
  expect_error(bread_prior(lambda0 = -1), "lambda0")
  expect_error(bread_prior(a0 = 0),       "a0")
  expect_error(bread_prior(b0 = -1),      "b0")
})
