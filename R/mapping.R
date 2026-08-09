#' Map probes to user-defined regions
#'
#' Computes probe-to-region overlaps, preserves feature metadata, and drops
#' regions with fewer than `min_probes` overlapping probes.
#'
#' Probes that overlap multiple regions are emitted once per region (long
#' format). Probes with no region are silently excluded from the returned
#' mapping. Regions excluded by `min_probes` (including those with zero
#' overlaps) are recorded on `attr(mapping, "dropped_regions")`.
#'
#' @param se A [SummarizedExperiment::SummarizedExperiment] with
#'   non-empty `rowRanges()`.
#' @param features A [GenomicRanges::GRanges] of regions. If `names(features)`
#'   is `NULL` or empty, IDs `region_1, region_2, ...` are generated.
#' @param min_probes Integer. Regions with fewer overlapping probes are
#'   dropped. Default `3L`.
#'
#' @return A `data.frame` with (at minimum) columns
#'   `probe_id`, `probe_idx`, `region_id`, `region_idx`, `n_probes`, plus any
#'   `mcols(features)` columns. Attributes:
#'   - `dropped_regions` : character vector of region IDs excluded.
#'   - `min_probes`      : the threshold applied.
#'   - `n_features_in`   : regions supplied.
#'   - `n_features_out`  : regions retained.
#'
#' @importFrom methods is
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom GenomicRanges findOverlaps
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom S4Vectors queryHits subjectHits mcols
#' @examples
#' suppressPackageStartupMessages(library(SummarizedExperiment))
#'
#' se  <- readRDS(system.file("extdata", "vitc_ag06561.rds", package = "BREAD"))
#' reg <- readRDS(system.file("extdata", "vitc_regions.rds", package = "BREAD"))
#' se_ctrl <- se[, se$condition == "ctrl"]
#'
#' mapping <- map_probes_to_features(se, reg)
#' head(mapping)
#' nrow(mapping)
#'
#' # Regions carrying fewer than `min_probes` probes are dropped and
#' # recorded in an attribute rather than silently disappearing.
#' attr(mapping, "dropped_regions")
#' @export
map_probes_to_features <- function(se, features, min_probes = 3L) {
  if (!methods::is(se, "SummarizedExperiment")) {
    stop("`se` must be a SummarizedExperiment.", call. = FALSE)
  }
  if (!methods::is(features, "GRanges")) {
    stop("`features` must be a GRanges.", call. = FALSE)
  }
  if (!is.numeric(min_probes) || length(min_probes) != 1L ||
      is.na(min_probes) || min_probes < 1L) {
    stop("`min_probes` must be a positive integer.", call. = FALSE)
  }
  min_probes <- as.integer(min_probes)

  probes <- SummarizedExperiment::rowRanges(se)
  if (is.null(probes) || length(probes) == 0L) {
    stop("`rowRanges(se)` is empty.", call. = FALSE)
  }

  # Probe IDs: GRanges names -> rownames(se) -> generated
  probe_ids <- names(probes)
  if (is.null(probe_ids) || any(probe_ids == "") || anyNA(probe_ids)) {
    probe_ids <- rownames(se)
  }
  if (is.null(probe_ids) || any(probe_ids == "") || anyNA(probe_ids)) {
    probe_ids <- paste0("probe_", seq_along(probes))
  }

  # Region IDs: GRanges names or generated
  region_ids <- names(features)
  if (is.null(region_ids) || any(region_ids == "") || anyNA(region_ids)) {
    region_ids <- paste0("region_", seq_along(features))
  }

  # findOverlaps() warns about non-overlapping seqlevels, which is not
  # actionable on its own -- the error below reports them instead.
  hits <- suppressWarnings(GenomicRanges::findOverlaps(probes, features))
  if (length(hits) == 0L) {
    .fmt_seqlevels <- function(x) {
      s <- GenomeInfoDb::seqlevels(x)
      if (length(s) == 0L) return("<none>")
      paste(c(s[seq_len(min(5L, length(s)))],
              if (length(s) > 5L) "..."), collapse = ", ")
    }
    stop("No probes overlap any of the ", length(features),
         " features. Check that `seqlevels()` and genome builds agree.\n",
         "  probe seqlevels  : ", .fmt_seqlevels(probes), "\n",
         "  feature seqlevels: ", .fmt_seqlevels(features),
         call. = FALSE)
  }

  qh <- S4Vectors::queryHits(hits)
  sh <- S4Vectors::subjectHits(hits)

  mapping <- data.frame(
    probe_id   = probe_ids[qh],
    probe_idx  = qh,
    region_id  = region_ids[sh],
    region_idx = sh,
    stringsAsFactors = FALSE
  )

  fmcols <- as.data.frame(S4Vectors::mcols(features), stringsAsFactors = FALSE)
  if (ncol(fmcols) > 0L) {
    mapping <- cbind(mapping, fmcols[sh, , drop = FALSE])
  }

  # Per-region probe counts
  counts <- table(mapping$region_id)
  mapping$n_probes <- as.integer(counts[mapping$region_id])

  # Regions with fewer than min_probes (including zero-overlap regions) are dropped
  kept <- names(counts)[counts >= min_probes]
  dropped <- setdiff(region_ids, kept)
  mapping <- mapping[mapping$region_id %in% kept, , drop = FALSE]
  rownames(mapping) <- NULL

  if (length(dropped) > 0L) {
    message("Dropped ", length(dropped), " of ", length(features),
            " regions with < ", min_probes, " probes.")
  }

  attr(mapping, "dropped_regions") <- dropped
  attr(mapping, "min_probes")      <- min_probes
  attr(mapping, "n_features_in")   <- length(features)
  attr(mapping, "n_features_out")  <- length(unique(mapping$region_id))
  mapping
}
