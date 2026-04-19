#' Map probes to user-defined regions
#'
#' Computes probe-region overlap, filters regions with fewer than `min_probes`,
#' and returns a tidy mapping with preserved feature metadata.
#'
#' @param se A [SummarizedExperiment::SummarizedExperiment] with `rowRanges()`.
#' @param features A [GenomicRanges::GRanges] of regions.
#' @param min_probes Minimum probes per region. Default 3.
#' @return A data frame with columns `probe_id`, `region_id`, feature metadata,
#'   and per-region overlap counts.
#' @keywords internal
map_probes_to_features <- function(se, features, min_probes = 3L) {
  stop("map_probes_to_features() not yet implemented")
}
