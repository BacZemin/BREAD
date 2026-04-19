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
    sex   = factor(rep(c("F", "M"),       length.out = n_samples),
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
