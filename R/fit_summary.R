#' Construct a BREAD prior for the summary-mode model
#'
#' Default prior for `fit_bread_summary()`: weakly-informative Normal-Inverse-Gamma.
#' `mu0` defaults to a zero vector (length inferred from the design) and
#' `Lambda0` to `lambda0 * I` (diagonal prior precision, inferred at fit time).
#' The inverse-gamma hyperparameters default to `a0 = b0 = 0.001`, which is
#' approximately Jeffreys.
#'
#' @param mu0 Optional prior mean vector for the coefficients. If `NULL`, a
#'   zero vector of the correct dimension is used at fit time.
#' @param Lambda0 Optional prior precision matrix. If `NULL`, `lambda0 * I`
#'   is used at fit time.
#' @param lambda0 Scalar prior precision used when `Lambda0` is `NULL`.
#'   Default `0.01` (weak).
#' @param a0,b0 Inverse-gamma hyperparameters for the residual variance.
#'   Defaults `a0 = b0 = 0.001`.
#'
#' @return A list with class `"bread_prior"`.
#' @examples
#' # Defaults: weak coefficient precision, near-flat inverse-gamma on the
#' # residual variance.
#' bread_prior()
#'
#' # A tighter prior, e.g. when regions are small and n is low
#' bread_prior(lambda0 = 0.05, a0 = 0.01, b0 = 0.01)
#' @export
bread_prior <- function(mu0 = NULL, Lambda0 = NULL,
                        lambda0 = 0.01, a0 = 0.001, b0 = 0.001) {
  stopifnot(
    is.null(mu0)     || is.numeric(mu0),
    is.null(Lambda0) || is.matrix(Lambda0),
    is.numeric(lambda0), length(lambda0) == 1L, lambda0 > 0,
    is.numeric(a0),      length(a0)      == 1L, a0 > 0,
    is.numeric(b0),      length(b0)      == 1L, b0 > 0
  )
  structure(
    list(mu0 = mu0, Lambda0 = Lambda0,
         lambda0 = lambda0, a0 = a0, b0 = b0),
    class = "bread_prior"
  )
}

#' Fit BREAD summary-mode Bayesian model
#'
#' v1 MVP fit. For each region, fits a Normal linear regression
#' `y ~ design` with a conjugate Normal-Inverse-Gamma prior, yielding an
#' analytical posterior — no MCMC. The marginal posterior of each coefficient
#' is a location-scale Student-t, used downstream by [posterior_summary()].
#'
#' @section Model:
#' For region r and sample i:
#' \deqn{y_{ir} \sim \mathcal{N}(X_i \beta_r,\, \sigma_r^2),}
#' \deqn{\beta_r \mid \sigma_r^2 \sim \mathcal{N}(\mu_0,\, \sigma_r^2 \Lambda_0^{-1}),}
#' \deqn{\sigma_r^2 \sim \mathrm{Inv\text{-}Gamma}(a_0, b_0).}
#' Posterior:
#' \deqn{\Lambda_n = X^\top X + \Lambda_0,\quad \mu_n = \Lambda_n^{-1}(X^\top y + \Lambda_0 \mu_0),}
#' \deqn{a_n = a_0 + n/2,\quad b_n = b_0 + \tfrac{1}{2}(y^\top y + \mu_0^\top \Lambda_0 \mu_0 - \mu_n^\top \Lambda_n \mu_n).}
#'
#' @param region_mat Region-by-sample numeric matrix (from [summarize_features()]).
#' @param coldata Sample metadata ([S4Vectors::DataFrame] or `data.frame`).
#' @param design One-sided formula.
#' @param contrast Character coefficient name of interest.
#' @param prior A [bread_prior()] object (or `NULL` for defaults).
#'
#' @return A list with:
#'   - `fits`: per-region list of `list(mu_n, Lambda_n_inv, a_n, b_n, n, error)`
#'   - `design_matrix`: model matrix `X` actually used
#'   - `coef_names`: coefficient names
#'   - `contrast`, `contrast_idx`: contrast name and its column index in `X`
#'   - `region_ids`: rownames of `region_mat`
#'   - `prior`: the prior applied (with `mu0`/`Lambda0` filled in)
#'
#' @importFrom methods is
#' @importFrom stats model.matrix
#' @keywords internal
fit_bread_summary <- function(region_mat, coldata, design, contrast,
                              prior = NULL) {
  if (!is.matrix(region_mat)) {
    stop("`region_mat` must be a numeric matrix (regions \u00d7 samples).",
         call. = FALSE)
  }
  if (!inherits(design, "formula") || length(design) != 2L) {
    stop("`design` must be a one-sided formula.", call. = FALSE)
  }
  if (!is.character(contrast) || length(contrast) != 1L || is.na(contrast)) {
    stop("`contrast` must be a single character coefficient name.",
         call. = FALSE)
  }
  cd_df <- as.data.frame(coldata)
  if (nrow(cd_df) != ncol(region_mat)) {
    stop("`coldata` has ", nrow(cd_df), " rows but `region_mat` has ",
         ncol(region_mat), " samples.", call. = FALSE)
  }

  X <- stats::model.matrix(design, data = cd_df)
  coef_names <- colnames(X)
  if (!contrast %in% coef_names) {
    stop("`contrast = \"", contrast, "\"` not found among design coefficients. ",
         "Available: ", paste(shQuote(coef_names), collapse = ", "), ".",
         call. = FALSE)
  }
  contrast_idx <- which(coef_names == contrast)

  if (is.null(prior)) prior <- bread_prior()
  p <- ncol(X)
  if (is.null(prior$mu0))     prior$mu0     <- rep(0, p)
  if (is.null(prior$Lambda0)) prior$Lambda0 <- diag(prior$lambda0, p)
  if (length(prior$mu0) != p)
    stop("`prior$mu0` length (", length(prior$mu0),
         ") does not match n coefficients (", p, ").", call. = FALSE)
  if (!identical(dim(prior$Lambda0), c(p, p)))
    stop("`prior$Lambda0` must be a ", p, "x", p, " matrix.", call. = FALSE)

  fits <- lapply(seq_len(nrow(region_mat)), function(i) {
    .fit_one_region(y       = region_mat[i, ],
                    X       = X,
                    mu0     = prior$mu0,
                    Lambda0 = prior$Lambda0,
                    a0      = prior$a0,
                    b0      = prior$b0)
  })
  names(fits) <- rownames(region_mat)

  list(
    fits          = fits,
    design_matrix = X,
    coef_names    = coef_names,
    contrast      = contrast,
    contrast_idx  = contrast_idx,
    region_ids    = rownames(region_mat),
    prior         = prior,
    region_mat    = region_mat,
    design        = design,
    coldata       = cd_df
  )
}

# Internal: conjugate NIG posterior for one region
.fit_one_region <- function(y, X, mu0, Lambda0, a0, b0) {
  p  <- length(mu0)
  ok <- !is.na(y)
  y  <- y[ok]
  X  <- X[ok, , drop = FALSE]
  n  <- length(y)

  na_fit <- function(reason) list(
    mu_n         = rep(NA_real_, p),
    Lambda_n_inv = matrix(NA_real_, p, p),
    a_n          = NA_real_,
    b_n          = NA_real_,
    n            = n,
    error        = reason
  )
  if (n < 2L) return(na_fit("too few non-NA samples"))

  XtX <- crossprod(X)
  Xty <- drop(crossprod(X, y))

  Lambda_n <- XtX + Lambda0
  chol_Ln  <- tryCatch(chol(Lambda_n), error = function(e) NULL)
  if (is.null(chol_Ln)) return(na_fit("singular posterior precision"))

  Lambda_n_inv <- chol2inv(chol_Ln)
  rhs  <- Xty + drop(Lambda0 %*% mu0)
  mu_n <- drop(Lambda_n_inv %*% rhs)

  a_n    <- a0 + n / 2
  qprior <- drop(crossprod(mu0, Lambda0 %*% mu0))
  qpost  <- drop(crossprod(mu_n, Lambda_n %*% mu_n))
  b_n    <- b0 + 0.5 * (sum(y^2) + qprior - qpost)
  b_n    <- max(b_n, .Machine$double.eps)

  list(mu_n = mu_n, Lambda_n_inv = Lambda_n_inv,
       a_n = a_n, b_n = b_n, n = n, error = NA_character_)
}
