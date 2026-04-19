#' Fit BREAD summary-mode model via brms (MCMC backend)
#'
#' Per-region Bayesian linear regression fitted with \pkg{brms}. The Stan
#' program is compiled **once** on the first region and reused across the
#' remaining regions via \code{update(..., recompile = FALSE)}, which makes
#' fitting N regions roughly "one Stan compile + N cheap sampling runs".
#'
#' @section Priors:
#' v1 uses brms's default priors (weakly informative). The \pkg{BREAD}
#' \code{bread_prior()} object is ignored by this backend — it only drives
#' the conjugate backend. Custom brms priors via the \code{brms_prior}
#' argument are planned for a future release.
#'
#' @param region_mat Region-by-sample numeric matrix.
#' @param coldata Sample metadata.
#' @param design One-sided formula.
#' @param contrast Character coefficient name.
#' @param iter,chains,cores,seed MCMC settings passed to \code{brms::brm}.
#' @param silent \code{brms::brm}'s \code{silent} level (default 2).
#' @param refresh \code{brms::brm}'s \code{refresh} (default 0, no progress).
#' @return A list with the same shape as \code{fit_bread_summary}, but each
#'   per-region fit carries \code{$draws} (a numeric vector of posterior
#'   draws for the contrast coefficient) instead of \code{mu_n / Lambda_n_inv
#'   / a_n / b_n}. \code{posterior_summary} and \code{posterior_draws} detect
#'   the shape automatically.
#'
#' @importFrom methods is
#' @importFrom stats as.formula model.matrix
#' @keywords internal
fit_bread_brms <- function(region_mat, coldata, design, contrast,
                           iter = 2000L, chains = 4L, cores = 4L,
                           seed = 1L, silent = 2L, refresh = 0L) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    stop("Package \"brms\" is required for backend = \"brms\". ",
         "Install via install.packages(\"brms\") or switch to backend = \"conjugate\".",
         call. = FALSE)
  }
  if (!is.matrix(region_mat)) stop("`region_mat` must be a matrix.", call. = FALSE)
  if (!inherits(design, "formula") || length(design) != 2L) {
    stop("`design` must be a one-sided formula.", call. = FALSE)
  }
  if (!is.character(contrast) || length(contrast) != 1L) {
    stop("`contrast` must be a single character coefficient name.", call. = FALSE)
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
  region_ids <- rownames(region_mat)

  # Build y ~ <design rhs> formula
  rhs    <- paste(deparse(design[[2L]]), collapse = " ")
  bform  <- stats::as.formula(paste("y ~", rhs))

  # brms prefixes fixed-effect coefficients with "b_"; our contrast arrives unprefixed
  draw_var <- paste0("b_", contrast)

  extract_draws <- function(brmsfit) {
    mat <- as.matrix(brmsfit, variable = draw_var)
    as.numeric(mat)
  }

  # Fit first region: this compiles Stan once
  first_rid <- region_ids[1L]
  df1 <- data.frame(y = region_mat[first_rid, ], cd_df, check.names = FALSE)
  fit_first <- tryCatch(
    suppressMessages(suppressWarnings(
      brms::brm(
        formula = bform, data = df1, family = stats::gaussian(),
        chains = chains, iter = iter, cores = cores, seed = seed,
        silent = silent, refresh = refresh
      )
    )),
    error = function(e) list(brms_error = conditionMessage(e))
  )
  if (!is.null(fit_first$brms_error)) {
    stop("brms failed to compile/sample the first region (", first_rid, "): ",
         fit_first$brms_error, call. = FALSE)
  }

  fits <- vector("list", length(region_ids))
  names(fits) <- region_ids
  fits[[first_rid]] <- list(
    n     = sum(!is.na(region_mat[first_rid, ])),
    error = NA_character_,
    draws = extract_draws(fit_first)
  )

  # Remaining regions: reuse compiled Stan model via update()
  if (length(region_ids) > 1L) {
    for (i in seq.int(2L, length(region_ids))) {
      rid  <- region_ids[i]
      df_i <- data.frame(y = region_mat[rid, ], cd_df, check.names = FALSE)
      res  <- tryCatch(
        suppressMessages(suppressWarnings(
          stats::update(fit_first, newdata = df_i, recompile = FALSE,
                        iter = iter, chains = chains, seed = seed + i,
                        silent = silent, refresh = refresh)
        )),
        error = function(e) list(brms_error = conditionMessage(e))
      )
      if (!is.null(res$brms_error)) {
        fits[[rid]] <- list(n = NA_integer_, error = res$brms_error, draws = NULL)
      } else {
        fits[[rid]] <- list(
          n     = sum(!is.na(region_mat[rid, ])),
          error = NA_character_,
          draws = extract_draws(res)
        )
      }
    }
  }

  list(
    fits          = fits,
    design_matrix = X,
    coef_names    = coef_names,
    contrast      = contrast,
    contrast_idx  = contrast_idx,
    region_ids    = region_ids,
    prior         = NULL,
    region_mat    = region_mat,
    design        = design,
    coldata       = cd_df
  )
}
