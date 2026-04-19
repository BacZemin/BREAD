# ---------------------------------------------------------------------------
# Build BREAD's packaged example:
#   inst/extdata/vitc_ag06561.rds  — SummarizedExperiment, 8 samples
#   inst/extdata/vitc_regions.rds  — GRanges of predefined feature-class regions
#
# Source:
#   SE   = /varidata/.../aging_serial_culture/vitc_fibroblasts/data/se_vitc_fibroblasts.rds
#   MNFT = /varidata/.../EPIC_manifest/EPICv2_manifestV3.rds
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
})
set.seed(2026L)

BP  <- "/varidata/research/projects/laird/jaemin.park/projects/BREAD"
SRC <- "/varidata/research/projects/laird/jaemin.park/projects/aging_serial_culture/vitc_fibroblasts/data/se_vitc_fibroblasts.rds"
MNF <- "/varidata/research/projects/laird/jaemin.park/EPIC_manifest/EPICv2_manifestV3.rds"

# 1. Load + subset to AG06561 only (8 samples)
se_full <- readRDS(SRC)
keep <- grepl("^AG06561_", as.character(colData(se_full)$group))
se <- se_full[, keep]
stopifnot(ncol(se) == 8L)

# Tidy colData: condition + passage as factors with reference levels
cd <- colData(se)
cd$condition <- factor(as.character(cd$condition), levels = c("ctrl", "aa57"))
cd$passage   <- factor(as.character(cd$passage),   levels = c("early", "late"))
colData(se)  <- cd

# 2. Load manifest, restrict to autosomal probes with coordinates
m <- readRDS(MNF)
m <- m[!is.na(m$CpG_chrm) & grepl("^chr[0-9]+$", m$CpG_chrm), ]

# 3. Attach rowRanges + selected mcols from manifest
common <- intersect(rownames(se), m$Probe_ID)
se <- se[common, ]
mm <- m[match(common, m$Probe_ID), ]
stopifnot(identical(mm$Probe_ID, rownames(se)))

rr <- GRanges(
  seqnames = mm$CpG_chrm,
  ranges   = IRanges(start = mm$CpG_beg, end = mm$CpG_end)
)
mcols(rr) <- DataFrame(
  Probe_ID      = mm$Probe_ID,
  CpG_Island    = mm$CpG_Island,
  REMC_ChromHMM = mm$REMC.ChromHMM,
  PMDsoloWCGW   = mm$PMDsoloWCGW,
  H3K27ME3      = mm$H3K27ME3,
  H3K27AC       = mm$H3K27AC,
  H3K4ME3       = mm$H3K4ME3
)
names(rr) <- mm$Probe_ID
rowRanges(se) <- rr

# 4. Build feature-class regions
make_class_regions <- function(gr_probes, tag, flank_bp = 2000L,
                               min_probes = 3L, max_regions = 100L) {
  if (length(gr_probes) == 0L) return(GRanges())
  ext <- resize(gr_probes, width = 2L * flank_bp + 1L, fix = "center")
  red <- reduce(ext)
  ov  <- countOverlaps(red, gr_probes)
  red <- red[ov >= min_probes]
  if (length(red) > max_regions) {
    red <- red[sample.int(length(red), max_regions)]
    red <- sort(red)
  }
  mcols(red)$feature_class <- tag
  names(red) <- sprintf("%s_%03d", tag, seq_along(red))
  red
}

rr_probes <- rowRanges(se)

regs <- list(
  PRC_CGI          = make_class_regions(
    rr_probes[rr_probes$H3K27ME3 & rr_probes$CpG_Island], "PRC_CGI"),
  Bivalent         = make_class_regions(
    rr_probes[rr_probes$H3K4ME3 & rr_probes$H3K27ME3], "Bivalent"),
  Active_promoter  = make_class_regions(
    rr_probes[rr_probes$CpG_Island & rr_probes$H3K4ME3 & !rr_probes$H3K27ME3],
    "Active_promoter"),
  Active_enhancer  = make_class_regions(
    rr_probes[rr_probes$H3K27AC & !rr_probes$CpG_Island & !rr_probes$H3K4ME3],
    "Active_enhancer"),
  PMD_soloWCGW     = make_class_regions(
    rr_probes[!is.na(rr_probes$PMDsoloWCGW) & rr_probes$PMDsoloWCGW],
    "PMD_soloWCGW")
)
regions <- do.call(c, unname(regs))
regions <- sort(regions)
message("regions per class:")
print(table(regions$feature_class))

# 5. Subset SE to probes overlapping our regions
ov <- findOverlaps(rowRanges(se), regions)
probe_idx <- sort(unique(queryHits(ov)))
se_sub <- se[probe_idx, ]

message("final SE: ", nrow(se_sub), " probes x ", ncol(se_sub), " samples")

# 6. Serialize to inst/extdata
dir.create(file.path(BP, "inst/extdata"), recursive = TRUE, showWarnings = FALSE)
saveRDS(se_sub,  file.path(BP, "inst/extdata/vitc_ag06561.rds"), compress = "xz")
saveRDS(regions, file.path(BP, "inst/extdata/vitc_regions.rds"), compress = "xz")

message("SE   RDS size: ",
        format(file.info(file.path(BP, "inst/extdata/vitc_ag06561.rds"))$size,
               big.mark = ","))
message("regs RDS size: ",
        format(file.info(file.path(BP, "inst/extdata/vitc_regions.rds"))$size,
               big.mark = ","))
