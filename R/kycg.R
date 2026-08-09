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
#' suppressPackageStartupMessages(library(SummarizedExperiment))
#'
#' se  <- readRDS(system.file("extdata", "vitc_ag06561.rds", package = "BREAD"))
#' reg <- readRDS(system.file("extdata", "vitc_regions.rds", package = "BREAD"))
#' se_ctrl <- se[, se$condition == "ctrl"]
#'
#' fit <- fit_bread(se_ctrl, reg, ~ passage)
#' table(results(fit)$classification)
#'
#' # The enrichment call downloads KnowYourCG reference databases, so it is
#' # not run here.
#' \dontrun{
#' enr <- bread_kycg(fit, which = "hypermethylated",
#'                   platform = "EPICv2")
#' head(enr[enr$FDR < 0.01, ])
#' }
#'
#' @param mtc_by_group Correct for multiple testing within each knowledgebase
#'   group rather than across all of them. Passed to
#'   `knowYourCG::testEnrichment()` when the installed version supports it
#'   (added after Bioconductor 3.20) and ignored with a message when it does
#'   not.
#' @param mtc_method Multiple-testing correction method, as for
#'   [stats::p.adjust()]. Same version caveat as `mtc_by_group`.
#'
#' @importFrom methods is
#' @importFrom utils packageVersion
#' @export
bread_kycg <- function(fit,
                       which        = c("hypermethylated", "hypomethylated"),
                       databases    = NULL,
                       platform     = c("EPIC", "EPICv2", "HM450", "MM285"),
                       universe     = NULL,
                       alternative  = "greater",
                       include_genes = FALSE,
                       mtc_by_group = TRUE,
                       mtc_method   = "fdr") {
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

  # Default DB selection by platform
  dbs <- databases
  if (is.null(dbs)) {
    groups <- tryCatch(knowYourCG::listDBGroups(),
                       error = function(e) {
                         stop("knowYourCG::listDBGroups() failed (it reaches ",
                              "ExperimentHub): ", conditionMessage(e),
                              "\nPass `databases = ` to work offline.",
                              call. = FALSE)
                       })
    dbs <- .kycg_default_dbs(platform, groups$Title)
    if (length(dbs) == 0L) {
      avail <- grep(paste0("^KYCG\\.", platform, "\\."), groups$Title,
                    value = TRUE)
      warning("bread_kycg(): no default KYCG groups matched platform '",
              platform, "'. ",
              if (length(avail))
                paste0("Available for this platform:\n  ",
                       paste(avail, collapse = "\n  "), "\n")
              else "No groups at all are registered for this platform. ",
              "Pass `databases = ` explicitly.", call. = FALSE)
      return(data.frame())
    }
  }

  # One enrichment run per requested classification level
  parts <- lapply(which_levels, function(lvl) {
    rids <- as.character(res$region_id[res$classification == lvl])
    pids <- unique(as.character(mapping$probe_id[mapping$region_id %in% rids]))
    if (length(pids) == 0L) return(NULL)

    args <- list(
      query         = pids,
      databases     = dbs,
      universe      = universe,
      alternative   = alternative,
      include_genes = include_genes,
      platform      = platform,
      silent        = TRUE,
      mtc_by_group  = mtc_by_group,
      mtc_method    = mtc_method
    )
    # knowYourCG gained mtc_by_group / mtc_method after Bioc 3.20. Filter
    # against the installed signature rather than testing a version number:
    # BREAD is developed against one release and CI runs against devel.
    keep <- names(args) %in% names(formals(knowYourCG::testEnrichment))
    if (!all(keep) && (!missing(mtc_by_group) || !missing(mtc_method))) {
      message("bread_kycg(): knowYourCG ", utils::packageVersion("knowYourCG"),
              " has no ", paste(names(args)[!keep], collapse = ", "),
              " argument; ignoring.")
    }

    enr <- tryCatch(
      do.call(knowYourCG::testEnrichment, args[keep]),
      error = function(e) {
        warning("knowYourCG::testEnrichment failed for '", lvl, "': ",
                conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    if (is.null(enr) || nrow(enr) == 0L) return(NULL)
    enr$query   <- lvl
    enr$n_query <- length(pids)
    enr
  })

  # Drop empties before rbind: a level with no probes used to contribute a
  # 2-column stub, which rbind refuses to combine with a real result table.
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0L) {
    warning("bread_kycg(): no classification level in `which` yielded any ",
            "probes to test. Requested: ",
            paste(shQuote(which_levels), collapse = ", "), ".",
            call. = FALSE)
    return(data.frame())
  }
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

# Default knowledgebase groups per platform, as name fragments matched right
# after "KYCG.<platform>.". Deliberately no trailing "\\.": the real MM285
# titles are KYCG.MM285.TFBSconsensus.20220116 and
# KYCG.MM285.HMconsensus.20220116, which an anchored "TFBS\\." cannot match --
# that is why mouse users silently got an empty data.frame.
#
# Technical annotations (Mask, chromosome, probeType, seqContext) are excluded
# on purpose: enrichment against them is not biologically interpretable.
.KYCG_DEFAULT_GROUPS <- list(
  EPIC   = c("TFBS", "ChromHMM", "chromHMM", "CGI"),
  EPICv2 = c("TFBS", "ChromHMM", "chromHMM", "CGI"),
  HM450  = c("TFBS", "ChromHMM", "chromHMM", "CGI"),
  MM285  = c("chromHMM", "TFBSconsensus", "HMconsensus",
             "designGroup", "tissueSignature", "metagene")
)

.kycg_default_dbs <- function(platform, titles) {
  frags <- .KYCG_DEFAULT_GROUPS[[platform]]
  if (is.null(frags) || length(titles) == 0L) return(character(0L))
  want <- paste0("^KYCG\\.", platform, "\\.(",
                 paste(frags, collapse = "|"), ")")
  titles[grepl(want, titles)]
}
