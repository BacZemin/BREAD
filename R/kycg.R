#' Enrich BREAD classifications against KnowYourCG databases
#'
#' For each region classified as `"hypermethylated"` or `"hypomethylated"`
#' (or any level you pick), extracts the constituent EPIC / HM450 / MM285
#' probes and runs [knowYourCG::testEnrichment()] against one or more CpG
#' annotation databases. This is how you get from a region-level BREAD
#' result to *"is the hyper set enriched for Polycomb targets / AP-1 binding
#' sites / PMDs / ..."*.
#'
#' @section Databases:
#' If `databases = NULL`, the function picks a sensible default set for the
#' platform you specify via `platform` (e.g. `"EPIC"`, `"EPICv2"`, `"HM450"`,
#' `"MM285"`). Pass a character vector of [knowYourCG::listDBGroups()] titles
#' to override.
#'
#' @param fit A [BreadFit].
#' @param which Classification level(s) to enrich. Default
#'   `c("hypermethylated", "hypomethylated")` runs both as separate queries.
#'   Use a single string to run one.
#' @param databases Character vector of KYCG database group titles. Default
#'   `NULL` selects defaults by `platform`.
#' @param platform One of `"EPIC"`, `"EPICv2"`, `"HM450"`, `"MM285"`. Default
#'   `"EPIC"`. Only used when `databases = NULL`.
#' @param universe Optional universe of probe IDs (character vector). Default
#'   `NULL` uses all probes in the BREAD mapping (i.e. probes covered by
#'   `features`), which is the right universe for a region-level enrichment.
#' @param alternative `"greater"` (default), `"two.sided"`, or `"less"`.
#' @param include_genes Passed through to [knowYourCG::testEnrichment()].
#'
#' @return A tidy `data.frame` with one row per tested set, containing a
#'   `query` column naming the classification level (`"hypermethylated"` /
#'   `"hypomethylated"` / ...) alongside standard KYCG columns (`dbname`,
#'   `estimate`, `p.value`, `FDR`, ...).
#'
#' @examples
#' \dontrun{
#'   fit <- fit_bread(se, features, ~ group)
#'   enr <- bread_kycg(fit, platform = "EPIC")
#'   head(enr[enr$FDR < 0.01, ])
#' }
#'
#' @importFrom methods is
#' @export
bread_kycg <- function(fit,
                       which        = c("hypermethylated", "hypomethylated"),
                       databases    = NULL,
                       platform     = c("EPIC", "EPICv2", "HM450", "MM285"),
                       universe     = NULL,
                       alternative  = "greater",
                       include_genes = FALSE) {
  if (!methods::is(fit, "BreadFit"))
    stop("`fit` must be a BreadFit.", call. = FALSE)
  if (!requireNamespace("knowYourCG", quietly = TRUE))
    stop("Package `knowYourCG` is required. Install via ",
         "BiocManager::install(\"knowYourCG\") and try again.",
         call. = FALSE)

  platform <- match.arg(platform)

  res <- fit@results
  if (is.null(res) || !"classification" %in% colnames(res))
    stop("`fit@results` has no classification column.", call. = FALSE)

  which_levels <- intersect(as.character(which),
                            levels(res$classification))
  if (length(which_levels) == 0L)
    stop("None of `which` matched the classification levels. ",
         "Available: ",
         paste(shQuote(levels(res$classification)), collapse = ", "),
         call. = FALSE)

  mapping <- fit@mapping
  if (is.null(mapping) || !all(c("region_id", "probe_id") %in% colnames(mapping)))
    stop("`fit@mapping` is missing or malformed.", call. = FALSE)

  # Universe: all probes that BREAD actually mapped to any surviving region
  if (is.null(universe))
    universe <- unique(as.character(mapping$probe_id))

  # Default DB selection by platform — TFBS + chromHMM + CGI
  dbs <- databases
  if (is.null(dbs)) {
    groups <- knowYourCG::listDBGroups()
    want   <- paste0("^KYCG\\.", platform, "\\.",
                     "(TFBS|ChromHMM|chromHMM|CGI)", "\\.")
    dbs <- groups$Title[grepl(want, groups$Title)]
    if (length(dbs) == 0L) {
      message("bread_kycg(): no default KYCG databases found for platform '",
              platform, "'; pass `databases = ...` explicitly.")
      return(data.frame())
    }
  }

  # One enrichment run per requested classification level
  do.call(rbind, lapply(which_levels, function(lvl) {
    rids <- as.character(res$region_id[res$classification == lvl])
    pids <- unique(as.character(mapping$probe_id[mapping$region_id %in% rids]))
    if (length(pids) == 0L) {
      return(data.frame(query = lvl, n_query = 0L,
                        stringsAsFactors = FALSE))
    }
    enr <- tryCatch(
      knowYourCG::testEnrichment(
        query         = pids,
        databases     = dbs,
        universe      = universe,
        alternative   = alternative,
        include_genes = include_genes,
        platform      = platform,
        silent        = TRUE
      ),
      error = function(e) {
        warning("knowYourCG::testEnrichment failed for '", lvl, "': ",
                conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    if (is.null(enr) || nrow(enr) == 0L)
      return(data.frame(query = lvl, n_query = length(pids),
                        stringsAsFactors = FALSE))
    enr$query   <- lvl
    enr$n_query <- length(pids)
    enr
  }))
}
