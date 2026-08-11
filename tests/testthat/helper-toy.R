# Toy SE: 20 probes on chr1 spaced every 1000bp, 8 samples, assay "M"
.make_toy_se <- function(n_probes = 20L, n_samples = 8L, seed = 1L) {
  set.seed(seed)
  probes <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges   = IRanges::IRanges(
      start = seq(1L, by = 1000L, length.out = n_probes),
      width = 1L
    )
  )
  names(probes) <- sprintf("cg%04d", seq_len(n_probes))
  m <- matrix(rnorm(n_probes * n_samples, mean = 0, sd = 1),
              nrow = n_probes, ncol = n_samples,
              dimnames = list(names(probes), sprintf("S%02d", seq_len(n_samples))))
  cd <- S4Vectors::DataFrame(
    group = factor(rep(c("young", "old"), length.out = n_samples),
                   levels = c("young", "old")),
    # Period 4 against group's period 2, so `~ group + sex` is full rank.
    # (When both alternated, the two were perfectly collinear and any design
    # using both silently fell back on the prior.)
    sex   = factor(rep(c("F", "F", "M", "M"), length.out = n_samples),
                   levels = c("F", "M")),
    row.names = colnames(m)
  )
  SummarizedExperiment::SummarizedExperiment(
    assays    = list(M = m),
    rowRanges = probes,
    colData   = cd
  )
}

# Toy features: regA covers probes 1-6 (6), regB probes 11-12 (2), regC probes 16-20 (5)
# With min_probes = 3L, regB is dropped.
.make_toy_features <- function() {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges   = IRanges::IRanges(
      start = c(     1L, 10001L, 15001L),
      end   = c(  5500L, 11500L, 20000L)
    ),
    feature_class = c("PRC", "CGI", "LAD")
  )
  names(gr) <- c("regA", "regB", "regC")
  gr
}

# Toy features where ONE region is defined by several disjoint ranges sharing
# a name -- the only way to pin a region to an exact probe set, since a single
# bounding interval would sweep in the probes between them.
# regD = probes 1-3 (1..2500) + probes 8-10 (7001..9500) = 6 probes, 2 ranges.
# regE = probes 16-20, 1 range. So: 2 distinct regions from 3 ranges.
.make_toy_features_dup <- function() {
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges   = IRanges::IRanges(
      start = c(     1L,  7001L, 15001L),
      end   = c(  2500L,  9500L, 20000L)
    ),
    feature_class = c("PRC", "PRC", "LAD")
  )
  names(gr) <- c("regD", "regD", "regE")
  gr
}

# Toy SE with injected signal: regA becomes hyper, regC becomes hypo,
# in "old" vs "young" under ~ group. Sample size large enough for
# prob_cutoff = 0.95 classification to recover truth.
.make_toy_signal_se <- function(n_samples = 24L, seed = 11L,
                                hyper_effect = 1.0, hypo_effect = -1.0) {
  se <- .make_toy_se(n_probes = 20L, n_samples = n_samples, seed = seed)
  old_mask <- SummarizedExperiment::colData(se)$group == "old"
  X <- SummarizedExperiment::assay(se, "M")
  # regA covers probes 1-6; regC covers probes 16-20 (see .make_toy_features())
  X[1:6,   old_mask] <- X[1:6,   old_mask] + hyper_effect
  X[16:20, old_mask] <- X[16:20, old_mask] + hypo_effect
  SummarizedExperiment::assay(se, "M") <- X
  se
}
