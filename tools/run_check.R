Sys.setenv(
  LANG = "C.UTF-8", LC_ALL = "C.UTF-8",
  RSTUDIO_PANDOC = "/varidata/research/projects/laird/jaemin.park/quarto/quarto-1.6.40/bin/tools/x86_64"
)
Sys.setenv("_R_CHECK_CRAN_INCOMING_"   = "false",
           "_R_CHECK_FORCE_SUGGESTS_"  = "false")
suppressPackageStartupMessages({ library(rcmdcheck) })

res <- rcmdcheck::rcmdcheck(
  "/varidata/research/projects/laird/jaemin.park/projects/BREAD",
  args      = c("--no-manual", "--as-cran"),
  check_dir = "/varidata/research/projects/laird/jaemin.park/projects/BREAD/docs/check",
  error_on  = "never",
  quiet     = FALSE
)
saveRDS(res, "/varidata/research/projects/laird/jaemin.park/projects/BREAD/docs/check/rcmdcheck_result.rds")
