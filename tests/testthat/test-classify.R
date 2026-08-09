# classify_regions() applies the decision rule; it does not recompute
# posterior probabilities. These tests drive it with hand-built posterior
# tables so the boundaries are exact rather than approximate.
#
# prob_hyper, prob_hypo and prob_rope partition the posterior, so fixtures
# must sum to 1 -- an incoherent triple can exercise states the sampler can
# never reach and would let a real bug hide.

`%|NULL|%` <- function(a, b) if (is.null(a)) b else a

.fake_post <- function(p_gt, p_lt, p_rope = NULL, ids = NULL) {
  n <- length(p_gt)
  out <- data.frame(
    region_id  = ids %|NULL|% sprintf("r%02d", seq_len(n)),
    prob_hyper = p_gt,
    prob_hypo  = p_lt,
    stringsAsFactors = FALSE
  )
  if (!is.null(p_rope)) out$prob_rope <- p_rope
  out
}

# Uniform draws from the 2-simplex (Dirichlet(1,1,1)).
.coherent_post <- function(n, seed = 42L) {
  withr::with_seed(seed, {
    e <- matrix(stats::rexp(3L * n), ncol = 3L)
    p <- e / rowSums(e)
  })
  .fake_post(p[, 1L], p[, 2L], p[, 3L])
}


test_that("classification levels are exact and in a fixed order", {
  out <- classify_regions(
    .fake_post(c(0.99, 0.01, 0.02, 0.34),
               c(0.01, 0.99, 0.02, 0.33),
               c(0.00, 0.00, 0.96, 0.33))
  )
  expect_s3_class(out$classification, "factor")
  expect_identical(
    levels(out$classification),
    c("hypermethylated", "hypomethylated", "unchanged", "inconclusive")
  )
  expect_identical(
    as.character(out$classification),
    c("hypermethylated", "hypomethylated", "unchanged", "inconclusive")
  )
})

test_that("the motivating case is called unchanged, not inconclusive", {
  # REGULON_promoter from the mousearray_609G CArG analysis: the posterior
  # sits almost entirely inside the ROPE. Before the fourth level existed
  # this was labelled `inconclusive`, indistinguishable from a region whose
  # posterior spanned the whole line.
  out <- classify_regions(.fake_post(0.021, 0.021, 0.958), prob_cutoff = 0.95)
  expect_identical(as.character(out$classification), "unchanged")
})

test_that("the cutoff comparison is inclusive (>=), not strict", {
  at    <- classify_regions(.fake_post(0.95, 0.0, 0.05), prob_cutoff = 0.95)
  below <- classify_regions(.fake_post(0.95 - 1e-12, 0.0, 0.05),
                            prob_cutoff = 0.95)
  expect_identical(as.character(at$classification),    "hypermethylated")
  expect_identical(as.character(below$classification), "inconclusive")

  rope_at    <- classify_regions(.fake_post(0.02, 0.03, 0.95),
                                 rope_cutoff = 0.95)
  rope_below <- classify_regions(.fake_post(0.02, 0.03, 0.95 - 1e-12),
                                 rope_cutoff = 0.95)
  expect_identical(as.character(rope_at$classification),    "unchanged")
  expect_identical(as.character(rope_below$classification), "inconclusive")
})

test_that("at sane cutoffs the four classes are mutually exclusive", {
  post <- .coherent_post(10000L)
  cls  <- classify_regions(post, prob_cutoff = 0.95,
                           rope_cutoff = 0.95)$classification
  n_clearing <- (post$prob_hyper >= 0.95) +
                (post$prob_hypo  >= 0.95) +
                (post$prob_rope  >= 0.95)
  expect_true(all(n_clearing <= 1L))
  expect_false(any(is.na(cls)))
  # Every row that clears something is labelled, and nothing else is
  expect_identical(cls != "inconclusive", n_clearing == 1L)
})

test_that("tightening both bars only ever moves regions to inconclusive", {
  post   <- .coherent_post(500L)
  loose  <- classify_regions(post, prob_cutoff = 0.60)$classification
  strict <- classify_regions(post, prob_cutoff = 0.99)$classification
  moved  <- loose != strict
  expect_true(all(strict[moved] == "inconclusive"))
})

test_that("rope_cutoff moves only the equivalence bar", {
  post   <- .coherent_post(500L)
  base   <- classify_regions(post, prob_cutoff = 0.95)$classification
  looser <- classify_regions(post, prob_cutoff = 0.95,
                             rope_cutoff = 0.70)$classification
  # Directional calls are untouched ...
  dir_base <- base   %in% c("hypermethylated", "hypomethylated")
  dir_loose<- looser %in% c("hypermethylated", "hypomethylated")
  expect_identical(dir_base, dir_loose)
  # ... and the only movement is inconclusive -> unchanged
  moved <- base != looser
  expect_true(all(base[moved]   == "inconclusive"))
  expect_true(all(looser[moved] == "unchanged"))
  expect_true(any(moved))

  # It defaults to prob_cutoff
  expect_identical(
    classify_regions(post, prob_cutoff = 0.8)$classification,
    classify_regions(post, prob_cutoff = 0.8, rope_cutoff = 0.8)$classification
  )
})

test_that("failed fits (NA probabilities) become inconclusive, never NA", {
  out <- classify_regions(.fake_post(c(NA, 0.99), c(NA, 0.001)))
  expect_false(any(is.na(out$classification)))
  expect_identical(as.character(out$classification)[1], "inconclusive")
})

test_that("an NA posterior is never called unchanged", {
  # A missing posterior is absence of evidence, not evidence of equivalence.
  out <- classify_regions(.fake_post(NA_real_, NA_real_, 1.0))
  expect_identical(as.character(out$classification), "inconclusive")
})

test_that("when several classes clear a low cutoff the largest wins", {
  out <- classify_regions(
    .fake_post(c(0.40, 0.31, 0.31), c(0.31, 0.40, 0.31), c(0.29, 0.29, 0.38)),
    prob_cutoff = 0.30, rope_cutoff = 0.30
  )
  expect_identical(as.character(out$classification),
                   c("hypermethylated", "hypomethylated", "unchanged"))

  # Exact ties resolve hyper > hypo > unchanged, as before the fourth level
  tie <- classify_regions(.fake_post(1/3, 1/3, 1/3), prob_cutoff = 0.30)
  expect_identical(as.character(tie$classification), "hypermethylated")
})

test_that("prob_rope is derived when absent and clamped when incoherent", {
  # Hand-built two-column tables remain a legitimate input
  derived <- classify_regions(.fake_post(0.01, 0.01), rope_cutoff = 0.95)
  expect_identical(as.character(derived$classification), "unchanged")

  # hyper + hypo > 1 cannot happen from a real posterior, but must not
  # produce a negative ROPE mass or an NA class
  bad <- classify_regions(.fake_post(0.8, 0.8), prob_cutoff = 0.95)
  expect_identical(as.character(bad$classification), "inconclusive")

  # A supplied prob_rope is authoritative over the complement
  supplied <- classify_regions(.fake_post(0.01, 0.01, 0.10),
                               rope_cutoff = 0.95)
  expect_identical(as.character(supplied$classification), "inconclusive")
})

test_that("delta and both cutoffs are recorded as attributes", {
  out <- classify_regions(.fake_post(0.99, 0.0), delta = 0.25,
                          prob_cutoff = 0.9, rope_cutoff = 0.7)
  expect_equal(attr(out, "delta"), 0.25)
  expect_equal(attr(out, "prob_cutoff"), 0.9)
  expect_equal(attr(out, "rope_cutoff"), 0.7)
})

test_that("bad input is rejected with an informative error", {
  expect_error(classify_regions("nope"), "data.frame")
  expect_error(classify_regions(data.frame(x = 1)), "missing required columns")
  expect_error(classify_regions(.fake_post(0.9, 0.1), prob_cutoff = 1),
               "`prob_cutoff` must be in \\(0, 1\\)")
  expect_error(classify_regions(.fake_post(0.9, 0.1), prob_cutoff = 0),
               "`prob_cutoff` must be in \\(0, 1\\)")
  expect_error(classify_regions(.fake_post(0.9, 0.1), rope_cutoff = 1),
               "`rope_cutoff` must be in \\(0, 1\\)")
})
