## Fast inner loop: regenerate docs, run every example, run the test suite.
## Skips vignettes and the full R CMD check entirely (~2 min vs ~14 min).
Sys.setenv(LANG = "C.UTF-8", LC_ALL = "C.UTF-8")
pkg <- "/varidata/research/projects/laird/jaemin.park/projects/BREAD"

cat("\n########## document() ##########\n")
suppressPackageStartupMessages(library(roxygen2))
roxygen2::roxygenise(pkg, clean = TRUE)

cat("\n########## NAMESPACE ##########\n")
cat(readLines(file.path(pkg, "NAMESPACE")), sep = "\n")

cat("\n########## run_examples() ##########\n")
suppressPackageStartupMessages(library(devtools))
ok <- TRUE
res <- tryCatch(
  devtools::run_examples(pkg, document = FALSE, run_donttest = TRUE),
  error = function(e) { ok <<- FALSE; message("EXAMPLES FAILED: ",
                                              conditionMessage(e)); NULL }
)
cat("\nexamples_ok:", ok, "\n")

cat("\n########## test() ##########\n")
tr <- tryCatch(devtools::test(pkg, stop_on_failure = FALSE),
               error = function(e) { message("TESTS ERRORED: ",
                                             conditionMessage(e)); NULL })
cat("\n########## DONE ##########\n")
