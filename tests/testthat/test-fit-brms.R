# Slow: Stan compile + sampling. ~60s on HPC. Skipped on CRAN, on systems
# without brms, and when $_R_CHECK_FORCE_SUGGESTS_ is FALSE and brms missing.
test_that("fit_bread(backend = 'brms') recovers injected signal", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("rstan")

  se <- .make_toy_signal_se(n_samples = 20L,
                            hyper_effect = 1.0, hypo_effect = -1.0)
  gr <- .make_toy_features()

  fit <- suppressMessages(
    fit_bread(
      se, gr, design = ~ group, contrast = "groupold",
      backend     = "brms",
      min_probes  = 3L,
      delta       = 0.10,
      prob_cutoff = 0.95,
      iter        = 400L,
      chains      = 2L,
      cores       = 2L,
      seed        = 1L
    )
  )

  expect_s4_class(fit, "BreadFit")
  expect_equal(fit@diagnostics$backend, "brms")
  expect_equal(fit@diagnostics$n_features_in,  3L)
  expect_equal(fit@diagnostics$n_features_out, 2L)

  # Empirical path: fits should have $draws, not mu_n
  f_regA <- fit@model$fits[["regA"]]
  expect_true(!is.null(f_regA$draws))
  expect_gt(length(f_regA$draws), 100L)

  cls <- classifications(fit)
  expect_equal(cls[["regA"]], "hypermethylated")
  expect_equal(cls[["regC"]], "hypomethylated")

  # posterior_summary works on the brms fit shape
  res <- results(fit)
  expect_true(is.na(res$df[1L]))  # df undefined for empirical path
  expect_true(all(!is.na(res$mean_effect)))
  expect_true(all(res$p_pos + res$p_neg >= 0.999))

  # posterior_draws returns actual MCMC draws (subsampled to n)
  d <- posterior_draws(fit, region_id = "regA", n = 200L, seed = 1L)
  expect_equal(nrow(d), 200L)
  expect_gt(mean(d$value), 0.5)
})
