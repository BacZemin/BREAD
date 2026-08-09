## Run BiocCheck against a freshly built BREAD tarball.
##
## Preflight only -- the authoritative BiocCheck runs on the Bioconductor
## devel container in CI (.github/workflows/bioc-check.yaml), because the
## HPC R is older than Bioc devel. This wrapper exists to catch the cheap
## problems before spending a CI cycle.
##
## Submit via sbatch (never the login node):
##   sbatch --partition=laird --mem=32G --cpus-per-task=4 --time=2:00:00 \
##     --wrap='cd <pkg> && Rscript tools/run_bioccheck.R'

Sys.setenv(
  LANG = "C.UTF-8", LC_ALL = "C.UTF-8",
  RSTUDIO_PANDOC = "/varidata/research/projects/laird/jaemin.park/quarto/quarto-1.6.40/bin/tools/x86_64"
)
Sys.setenv("_R_CHECK_FORCE_SUGGESTS_" = "false")

suppressPackageStartupMessages({
  library(BiocCheck)
})

pkg_dir   <- "/varidata/research/projects/laird/jaemin.park/projects/BREAD"
check_dir <- file.path(pkg_dir, "docs", "check")
dir.create(check_dir, recursive = TRUE, showWarnings = FALSE)

## ---- 1. Git-clone-level checks (run on the source dir, not the tarball) ----
message("== BiocCheckGitClone ==")
gitres <- try(BiocCheck::BiocCheckGitClone(pkg_dir), silent = TRUE)
if (inherits(gitres, "try-error")) {
  message("BiocCheckGitClone failed: ", conditionMessage(attr(gitres, "condition")))
}

## ---- 2. Build a tarball ----------------------------------------------------
## BiocCheck's `new-package` checks want the built tarball, not the source dir.
message("== R CMD build ==")
old <- setwd(check_dir)
on.exit(setwd(old), add = TRUE)

build_log <- system2(
  file.path(R.home("bin"), "R"),
  c("CMD", "build", "--no-resave-data", shQuote(pkg_dir)),
  stdout = TRUE, stderr = TRUE
)
cat(build_log, sep = "\n")

tarballs <- list.files(check_dir, pattern = "^BREAD_.*\\.tar\\.gz$", full.names = TRUE)
if (!length(tarballs)) stop("R CMD build produced no tarball; see log above.")
tarball <- tarballs[order(file.mtime(tarballs), decreasing = TRUE)][1]
message("Using tarball: ", tarball)

## ---- 3. BiocCheck ----------------------------------------------------------
message("== BiocCheck (new-package = TRUE) ==")
res <- BiocCheck::BiocCheck(tarball, `new-package` = TRUE)

saveRDS(res, file.path(check_dir, "bioccheck_result.rds"))

## ---- 4. Machine-readable summary ------------------------------------------
## The printed BiocCheck output is long; this block is what gets read back
## over SSH to build the fix list.
summarise <- function(res) {
  for (sev in c("error", "warning", "note")) {
    items <- tryCatch(res[[sev]], error = function(e) NULL)
    cat("\n########## ", toupper(sev), " (", length(items), ") ##########\n", sep = "")
    if (!length(items)) next
    for (nm in names(items)) {
      cat("- ", nm, "\n", sep = "")
      det <- items[[nm]]
      if (length(det)) cat(paste0("    ", unlist(det), collapse = "\n"), "\n", sep = "")
    }
  }
}
cat("\n\n================ BIOCCHECK SUMMARY ================\n")
try(summarise(res))
cat("\n=================== END SUMMARY ===================\n")

cat("\nBiocCheck artifacts:\n")
print(list.files(check_dir, pattern = "BiocCheck", full.names = TRUE))
