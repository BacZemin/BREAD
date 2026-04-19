#' BREAD plotting helpers
#'
#' Three publication-oriented plot functions for [BreadFit] output:
#' [plot_region_posterior()], [plot_region_data()], and [plot_feature_set()].
#' All return [ggplot2::ggplot] objects and use the MetBrewer `Cross` palette
#' (see [bread_colors()]).
#'
#' @name BREAD-plots
#' @keywords internal
NULL

#' Posterior density of the contrast coefficient per region
#'
#' Plots the analytical scaled-Student-t posterior of the contrast coefficient
#' for one or more regions, with vertical guides at `0` and `+/- delta` and a
#' color by final classification. When multiple regions are supplied the plot
#' facets by region.
#'
#' @param fit A [BreadFit].
#' @param region_id Character vector of region IDs. If `NULL` (default),
#'   all fitted regions are plotted.
#' @param n_grid Number of grid points for the density curve. Default `500`.
#' @param show_delta Logical. Draw dashed guides at `+/- delta`? Default `TRUE`.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @importFrom methods is
#' @importFrom stats dt
#' @importFrom rlang .data
#' @export
plot_region_posterior <- function(fit, region_id = NULL,
                                  n_grid = 500L, show_delta = TRUE) {
  if (!methods::is(fit, "BreadFit"))
    stop("`fit` must be a BreadFit.", call. = FALSE)
  fits <- fit@model$fits
  if (is.null(fits))
    stop("BreadFit has no `model$fits`; was fit_bread() run successfully?",
         call. = FALSE)
  if (is.null(region_id)) region_id <- names(fits)
  region_id <- as.character(region_id)

  unknown <- setdiff(region_id, names(fits))
  if (length(unknown) > 0L)
    stop("region_id(s) not found: ",
         paste(shQuote(unknown), collapse = ", "), call. = FALSE)

  k     <- fit@model$contrast_idx
  delta <- fit@params$delta
  ci    <- if (!is.null(attr(fit@results, "ci"))) attr(fit@results, "ci") else 0.95

  rows <- lapply(region_id, function(rid) {
    f <- fits[[rid]]
    if (!is.na(f$error)) return(NULL)
    nu <- 2 * f$a_n
    mu <- f$mu_n[k]
    s  <- sqrt(max((f$b_n / f$a_n) * f$Lambda_n_inv[k, k], .Machine$double.eps))
    x  <- seq(mu - 5 * s, mu + 5 * s, length.out = n_grid)
    y  <- stats::dt((x - mu) / s, df = nu) / s
    data.frame(region_id = rid, x = x, y = y,
               stringsAsFactors = FALSE)
  })
  dens <- do.call(rbind, rows)
  if (is.null(dens))
    stop("No plottable regions (all fits failed).", call. = FALSE)

  res_cols <- c("region_id", "classification", "mean_effect", "ci_lo", "ci_hi")
  meta <- fit@results[, intersect(res_cols, colnames(fit@results)), drop = FALSE]
  dens <- merge(dens, meta, by = "region_id", sort = FALSE)

  p <- ggplot2::ggplot(dens,
         ggplot2::aes(x = .data$x, y = .data$y,
                      color = .data$classification)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = 0,
                        linetype = "dotted", color = "gray40") +
    ggplot2::scale_color_manual(values = .col_classification, drop = FALSE,
                                name = "classification") +
    ggplot2::labs(
      x = "region effect (M-value units)",
      y = "posterior density",
      title = sprintf("Posterior of contrast: %s",
                      fit@params$contrast %||% "")
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )

  if (isTRUE(show_delta))
    p <- p + ggplot2::geom_vline(xintercept = c(-delta, delta),
                                 linetype = "dashed", color = "gray60")

  if (length(unique(dens$region_id)) > 1L)
    p <- p + ggplot2::facet_wrap(~ .data$region_id, scales = "free_y")

  p
}

#' Raw region-level values by contrast group
#'
#' Boxplot + jitter of the summarized region values for a single region,
#' grouped by the first variable in the design formula.
#'
#' @param fit A [BreadFit].
#' @param region_id Single region ID.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @importFrom methods is
#' @importFrom rlang .data
#' @export
plot_region_data <- function(fit, region_id) {
  if (!methods::is(fit, "BreadFit"))
    stop("`fit` must be a BreadFit.", call. = FALSE)
  if (length(region_id) != 1L)
    stop("Supply a single `region_id`.", call. = FALSE)

  region_mat <- fit@model$region_mat
  if (is.null(region_mat))
    stop("`fit@model$region_mat` missing; refit with current BREAD.",
         call. = FALSE)
  if (!region_id %in% rownames(region_mat))
    stop("region_id `", region_id, "` not found in fit.", call. = FALSE)

  coldata <- fit@model$coldata
  design  <- fit@model$design
  dvars   <- all.vars(design)
  if (length(dvars) == 0L)
    stop("Design has no variables; nothing to group by.", call. = FALSE)
  group_var <- dvars[1L]

  y <- region_mat[region_id, ]
  df <- data.frame(
    sample = colnames(region_mat),
    y      = unname(y),
    group  = coldata[[group_var]],
    stringsAsFactors = FALSE
  )

  # Reuse the two-color group palette; wider categorical fallback via Cross
  group_levels <- unique(df$group)
  if (length(group_levels) <= 2L) {
    pal <- stats::setNames(.col_group[seq_along(group_levels)], group_levels)
  } else {
    pal <- stats::setNames(.bread_cross[seq_along(group_levels)], group_levels)
  }

  eff  <- fit@results$mean_effect[fit@results$region_id == region_id]
  cls  <- as.character(fit@results$classification[fit@results$region_id == region_id])

  ggplot2::ggplot(df,
    ggplot2::aes(x = .data$group, y = .data$y,
                 color = .data$group, fill = .data$group)) +
    ggplot2::geom_boxplot(alpha = 0.15, outlier.shape = NA,
                          width = 0.5, color = "gray30") +
    ggplot2::geom_jitter(width = 0.15, size = 2.2, alpha = 0.9) +
    ggplot2::scale_color_manual(values = pal, guide = "none") +
    ggplot2::scale_fill_manual(values  = pal, guide = "none") +
    ggplot2::labs(
      x = group_var,
      y = "region summary (M-value)",
      title    = region_id,
      subtitle = sprintf("%s | effect = %.2f", cls, eff)
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}

#' Feature-set classification summary
#'
#' Bar chart of classification counts across all fitted regions. When
#' `feature_class_col` is supplied, bars are stacked by feature class so users
#' can see, e.g., how PRC / CGI / LAD subsets partition into hyper vs. hypo.
#'
#' @param fit A [BreadFit].
#' @param feature_class_col Optional column name in `fit@mapping` to stratify by.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @importFrom methods is
#' @importFrom rlang .data
#' @export
plot_feature_set <- function(fit, feature_class_col = NULL) {
  if (!methods::is(fit, "BreadFit"))
    stop("`fit` must be a BreadFit.", call. = FALSE)
  res <- fit@results
  if (is.null(res) || !"classification" %in% colnames(res))
    stop("`fit@results` has no classification column.", call. = FALSE)

  if (!is.null(feature_class_col)) {
    mp <- fit@mapping
    if (is.null(mp) || !feature_class_col %in% colnames(mp)) {
      stop("`feature_class_col = \"", feature_class_col,
           "\"` not found in fit@mapping.", call. = FALSE)
    }
    lut <- unique(mp[, c("region_id", feature_class_col)])
    names(lut)[2L] <- "feature_class"
    res <- merge(res, lut, by = "region_id", all.x = TRUE)
    p <- ggplot2::ggplot(res,
           ggplot2::aes(x = .data$feature_class, fill = .data$classification)) +
      ggplot2::geom_bar(position = "stack") +
      ggplot2::labs(x = feature_class_col, y = "n regions")
  } else {
    p <- ggplot2::ggplot(res,
           ggplot2::aes(x = .data$classification, fill = .data$classification)) +
      ggplot2::geom_bar() +
      ggplot2::labs(x = NULL, y = "n regions")
  }

  p +
    ggplot2::scale_fill_manual(values = .col_classification, drop = FALSE,
                               name = "classification") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "right",
                   plot.title = ggplot2::element_text(face = "bold"))
}
