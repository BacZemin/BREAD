suppressPackageStartupMessages({ library(pkgload); library(rmarkdown) })
BP <- "/varidata/research/projects/laird/jaemin.park/projects/BREAD"
pkgload::load_all(BP, quiet = TRUE)
rmarkdown::render(
  file.path(BP, "vignettes/bread-intro.Rmd"),
  output_dir  = file.path(BP, "docs"),
  output_file = "bread-intro.html",
  output_options = list(self_contained = TRUE),
  quiet = TRUE
)
cat("OK ->", file.path(BP, "docs/bread-intro.html"), "\n")
