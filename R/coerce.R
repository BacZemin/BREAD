#' Assemble a SummarizedExperiment for BREAD from a matrix
#'
#' BREAD models a `SummarizedExperiment` carrying row-level genomic
#' coordinates. Most methylation pipelines do not hand you one: `sesame`'s
#' `openSesame()` returns a plain beta matrix, and its packaged example data
#' are `list(betas = <matrix>, sampleInfo = <data.frame>)`. This helper builds
#' the object BREAD needs, so a matrix workflow does not stall at the first
#' step. [fit_bread()] calls it for you; use it directly when you want to
#' coerce once and reuse the result.
#'
#' @section Where coordinates come from:
#' A bare matrix has no coordinates, so you must supply them one of two ways:
#' pass `rowRanges` (a [GenomicRanges::GRanges], ideally named by probe ID),
#' or pass `platform` to look the manifest up through `sesameData`.
#'
#' **The platform is never guessed.** `cg########` identifiers are shared
#' across HM450, EPIC and MM285, so inferring the array from probe names would
#' silently return the wrong coordinates for a substantial fraction of probes,
#' assign them to the wrong regions, and produce confident, wrong biology with
#' no error anywhere. One word from you removes that entire failure mode.
#'
#' @param x A `SummarizedExperiment` (returned unchanged), a probe-by-sample
#'   `matrix`, or a `list` with a `betas` element and sample metadata under
#'   `sampleInfo`, `meta` or `pd`.
#' @param colData Sample metadata: a `data.frame` or `DataFrame` with one row
#'   per column of `x`. If it has rownames they are matched against
#'   `colnames(x)` and reordered; otherwise rows are assumed to be in column
#'   order and a warning is emitted. Required unless `design` has no variables.
#' @param rowRanges A [GenomicRanges::GRanges] of probe coordinates. If named,
#'   it is subset and reordered to `rownames(x)`; unnamed, it must already be
#'   in row order.
#' @param platform Array platform for the `sesameData` manifest lookup, e.g.
#'   `"EPIC"`, `"EPICv2"`, `"HM450"`, `"MM285"`. Requires the `sesameData`
#'   package. Ignored when `rowRanges` is supplied.
#' @param assay_name Name for the assay. Defaults to `"betas"` when the values
#'   all fall in \[0, 1\] and `"M"` otherwise, matching what
#'   [fit_bread()] auto-detects.
#'
#' @return A [SummarizedExperiment::SummarizedExperiment].
#'
#' @importFrom methods is
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors DataFrame
#' @seealso [fit_bread()]
#' @examples
#' suppressPackageStartupMessages(library(SummarizedExperiment))
#'
#' # Take a packaged SE apart, then put it back together the matrix way
#' se  <- readRDS(system.file("extdata", "vitc_ag06561.rds", package = "BREAD"))
#' mat <- assay(se, "betas")
#' cd  <- as.data.frame(colData(se))
#' gr  <- rowRanges(se)
#'
#' se2 <- bread_se(mat, colData = cd, rowRanges = gr)
#' se2
#' @export
bread_se <- function(x, colData = NULL, rowRanges = NULL, platform = NULL,
                     assay_name = NULL) {
  .as_bread_se(x, colData = colData, rowRanges = rowRanges,
               platform = platform, assay_name = assay_name)
}

# Internal workhorse. Kept separate so fit_bread() can distinguish "the user
# passed these" from "these are defaults", which drives the SE-plus-extras
# error below.
.as_bread_se <- function(x, colData = NULL, rowRanges = NULL, platform = NULL,
                         assay_name = NULL) {

  # --- already canonical -----------------------------------------------
  if (methods::is(x, "SummarizedExperiment")) {
    extras <- c(colData = !is.null(colData),
                rowRanges = !is.null(rowRanges),
                platform = !is.null(platform))
    if (any(extras)) {
      stop("`", paste(names(extras)[extras], collapse = "`, `"),
           "` supplied alongside a SummarizedExperiment. These arguments ",
           "describe matrix input only; edit the object itself instead.",
           call. = FALSE)
    }
    return(x)
  }

  # --- sesameData-style list -------------------------------------------
  if (is.list(x) && !is.data.frame(x) && "betas" %in% names(x)) {
    meta_nm <- intersect(c("sampleInfo", "meta", "pd"), names(x))
    if (is.null(colData) && length(meta_nm) > 0L) colData <- x[[meta_nm[1L]]]
    x <- x[["betas"]]
  }

  # --- matrix ------------------------------------------------------------
  if (is.null(dim(x)) || length(dim(x)) != 2L || is.data.frame(x)) {
    stop("`x` must be a SummarizedExperiment or a probe-by-sample matrix, ",
         "not <", class(x)[1], ">.", call. = FALSE)
  }
  if (!is.matrix(x)) x <- as.matrix(x)

  if (is.null(rownames(x))) {
    stop("`x` needs rownames giving probe IDs: they are what ties the ",
         "matrix to genomic coordinates.", call. = FALSE)
  }
  if (is.null(colnames(x))) {
    stop("`x` needs colnames giving sample IDs.", call. = FALSE)
  }

  cd <- .align_coldata(colData, colnames(x))
  rr <- .resolve_rowranges(rownames(x), rowRanges, platform)

  # Drop probes without coordinates rather than failing the whole call
  keep <- !is.na(rr$idx)
  if (!all(keep)) {
    message("Dropping ", sum(!keep), " of ", nrow(x),
            " probes with no coordinates in the supplied manifest.")
    x <- x[keep, , drop = FALSE]
    rr$gr <- rr$gr[keep]
  }
  if (nrow(x) == 0L) {
    stop("No probes in `x` have coordinates.", call. = FALSE)
  }

  if (is.null(assay_name)) {
    assay_name <- if (identical(.detect_input_scale(x), "Beta")) "betas" else "M"
  }
  assays <- list(x); names(assays) <- assay_name

  SummarizedExperiment::SummarizedExperiment(
    assays    = assays,
    rowRanges = rr$gr,
    colData   = cd
  )
}

# A data.frame always reports rownames, even when the user never set any --
# base R substitutes "1", "2", ... and stores them as a compact integer
# attribute. Treating those as sample IDs would make every rownameless
# colData look like a total mismatch, so distinguish real labels from the
# automatic ones.
.explicit_rownames <- function(df) {
  if (methods::is(df, "DataFrame")) return(rownames(df))
  rn <- attr(df, "row.names")
  if (is.null(rn) || is.integer(rn)) return(NULL)
  as.character(rn)
}

# Match sample metadata to the matrix columns, reordering when it is safe to.
.align_coldata <- function(colData, sample_ids) {
  n <- length(sample_ids)
  if (is.null(colData)) {
    # A design with no variables (~ 1) is legal; let validate_bread_input()
    # be the one to complain if the design actually needs columns.
    return(S4Vectors::DataFrame(row.names = sample_ids))
  }
  if (!is.data.frame(colData) && !methods::is(colData, "DataFrame")) {
    stop("`colData` must be a data.frame or DataFrame, not <",
         class(colData)[1], ">.", call. = FALSE)
  }
  cd <- S4Vectors::DataFrame(as.data.frame(colData, stringsAsFactors = FALSE))

  rn <- .explicit_rownames(colData)
  if (!is.null(rn)) {
    idx <- match(sample_ids, rn)
    if (anyNA(idx)) {
      miss <- sample_ids[is.na(idx)]
      stop("`colData` has no row for ", sum(is.na(idx)), " sample(s) in ",
           "`colnames(x)`: ",
           paste(shQuote(miss[seq_len(min(5L, length(miss)))]),
                 collapse = ", "),
           if (length(miss) > 5L) ", ..." else "", ".", call. = FALSE)
    }
    cd <- cd[idx, , drop = FALSE]
  } else {
    if (nrow(cd) != n) {
      stop("`colData` has ", nrow(cd), " rows but `x` has ", n,
           " columns.", call. = FALSE)
    }
    warning("`colData` has no rownames; assuming its rows are in the same ",
            "order as `colnames(x)`.", call. = FALSE)
  }
  rownames(cd) <- sample_ids
  cd
}

# Return list(gr = <GRanges in probe order>, idx = <match index, NA = absent>).
.resolve_rowranges <- function(probe_ids, rowRanges, platform) {
  n <- length(probe_ids)

  if (!is.null(rowRanges)) {
    if (!methods::is(rowRanges, "GRanges")) {
      stop("`rowRanges` must be a GRanges, not <", class(rowRanges)[1], ">.",
           call. = FALSE)
    }
    if (is.null(names(rowRanges))) {
      if (length(rowRanges) != n) {
        stop("Unnamed `rowRanges` has ", length(rowRanges), " ranges but ",
             "`x` has ", n, " rows. Name it by probe ID, or supply one ",
             "range per row in order.", call. = FALSE)
      }
      gr <- rowRanges
      names(gr) <- probe_ids
      return(list(gr = gr, idx = seq_len(n)))
    }
    idx <- match(probe_ids, names(rowRanges))
    if (all(is.na(idx))) {
      # A very common cause on EPICv2: the manifest keeps the replicate
      # suffix (cg00381604_BC11) while the analysis pipeline stripped it.
      suffixed  <- any(grepl("_[A-Z]{2}[0-9]{2}$", names(rowRanges)))
      bare      <- !any(grepl("_[A-Z]{2}[0-9]{2}$", probe_ids))
      stop("None of `rownames(x)` appear in `names(rowRanges)`. ",
           "Are these the same platform?",
           if (suffixed && bare)
             paste0("\n  The coordinates carry EPICv2 replicate suffixes ",
                    "(e.g. '", names(rowRanges)[1], "') but `x` does not ",
                    "(e.g. '", probe_ids[1], "'). Match them before ",
                    "calling BREAD.")
           else "",
           call. = FALSE)
    }
    gr <- rowRanges[ifelse(is.na(idx), 1L, idx)]
    names(gr) <- probe_ids
    return(list(gr = gr, idx = idx))
  }

  if (!is.null(platform)) {
    if (!requireNamespace("sesameData", quietly = TRUE)) {
      stop("`platform` lookup needs the 'sesameData' package. Install it, ",
           "or pass `rowRanges` directly.", call. = FALSE)
    }
    man <- sesameData::sesameData_getManifestGRanges(platform)
    return(.resolve_rowranges(probe_ids, man, NULL))
  }

  hint <- if (any(grepl("_[A-Z]{2}[0-9]{2}$", probe_ids))) {
    " (the probe ID suffixes look like EPICv2)"
  } else ""
  stop("No probe coordinates. Supply `rowRanges` (a GRanges named by probe ",
       "ID) or `platform` (looked up via sesameData)", hint, ". BREAD does ",
       "not guess the platform from probe IDs: cg-numbers are shared across ",
       "arrays and a wrong guess gives wrong coordinates silently.",
       call. = FALSE)
}
