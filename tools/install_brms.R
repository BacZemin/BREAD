# brms + stan stack installer — Posit PPM (rhel9) binaries
options(
  Ncpus = 4L,
  install.packages.check.source = "no",
  repos = c(
    PPM  = "https://packagemanager.posit.co/cran/__linux__/rhel9/latest",
    CRAN = "https://cran.rstudio.com/"
  ),
  # PPM serves binaries only when User-Agent reports the platform
  HTTPUserAgent = sprintf(
    "R/%s R (%s)", getRversion(),
    paste(getRversion(), R.version$platform, R.version$arch, R.version$os)
  )
)
cat(R.version.string, "\n")
print(.libPaths())

pkgs <- c("StanHeaders", "rstan", "posterior", "brms", "tidybayes")
installed <- rownames(installed.packages())
need <- setdiff(pkgs, intersect(installed, pkgs))
cat("\nneed:", if (length(need)) paste(need, collapse = ", ") else "nothing", "\n\n")

if (length(need) > 0L) install.packages(need, Ncpus = 4L)

cat("\n=== installed versions ===\n")
for (p in pkgs) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) "MISSING")
  cat(sprintf("%-14s %s\n", p, v))
}

cat("\n=== brms load test ===\n")
ok <- tryCatch({
  suppressPackageStartupMessages(library(brms))
  "OK"
}, error = function(e) conditionMessage(e))
cat(ok, "\n")

cat("\n=== rstan smoke test (trivial Normal fit) ===\n")
suppressPackageStartupMessages(library(rstan))
code <- "data { int<lower=1> N; vector[N] y; }\nparameters { real mu; real<lower=0> sigma; }\nmodel { y ~ normal(mu, sigma); }\n"
set.seed(1)
fit <- tryCatch(
  rstan::stan(model_code = code,
              data = list(N = 30L, y = rnorm(30L, 2, 0.5)),
              iter = 500, chains = 1, refresh = 0, verbose = FALSE),
  error = function(e) { cat("FAIL:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(fit)) {
  s <- summary(fit)$summary[1:2, c("mean","sd","Rhat")]
  print(s)
  cat("\nstan compile+sample: OK\n")
}
