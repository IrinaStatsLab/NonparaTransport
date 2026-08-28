#' Plot Fitted Nonparanormal Frechet Regression Bivariate Densities
#'
#' Visualizes fitted bivariate marginal densities from an \code{\link{npt_frechetreg}}
#' object across selected predictor evaluation points.
#'
#' @details
#' For each requested predictor evaluation point \eqn{z}, the function reconstructs the continuous
#' bivariate distributional response from the fitted nonparanormal parameters:
#' \enumerate{
#'   \item \strong{Latent Gaussian sampling}: Draws Monte Carlo samples from the latent
#'     bivariate Gaussian copula:
#'     \deqn{L \sim \mathcal{N}_2(0, \hat{R}(z))}
#'   \item \strong{Probability integral transform}: Maps latent Gaussian draws to uniform
#'     probabilities on \eqn{(0, 1)}:
#'     \deqn{U_j = \Phi(L_j), \quad j = 1, 2}
#'   \item \strong{Quantile transformation}: Maps uniform coordinates to the observational
#'     domain via the fitted marginal quantile functions:
#'     \deqn{Y_j = \hat{Q}_j(U_j; z), \quad j = 1, 2}
#'   \item \strong{Shared 2D Kernel Density Estimation (KDE)}: Evaluates a smooth bivariate
#'     density on a common coordinate grid with pooled bandwidths and a shared color scale,
#'     enabling direct visual comparison across different predictor values.
#' }
#'
#' @param fit An \code{npt_frechetreg} object created by \code{\link{npt_frechetreg}}
#'   (evaluated either at training points or at target points \code{Z}).
#' @param variables A character vector of two variable names or integer indices specifying
#'   the bivariate pair to plot. Defaults to the first two variables \code{c(1, 2)}.
#' @param predictions Prediction points to plot, supplied as integer row indices or
#'   character names. If \code{NULL} (default), selects up to five approximately equally
#'   spaced prediction rows across the range of \eqn{Z}.
#' @param n_samples Number of Monte Carlo draws per prediction panel to reconstruct the
#'   continuous bivariate distribution (default is `4000L`).
#' @param seed Optional random seed integer for reproducible Monte Carlo sampling.
#'   The caller's random seed state is preserved and restored after execution.
#' @param grid_size Number of evaluation points along each axis for the 2D density grid (default is `100L`).
#' @param trim A length-2 numeric vector of lower and upper probabilities in \eqn{[0, 1]}
#'   defining the robust shared display range across panels (default is `c(0.01, 0.99)`).
#' @param bandwidth_adjust Multiplier applied to the pooled Gaussian KDE bandwidths.
#'   Values above one produce smoother density surfaces; the default is `1.4`.
#' @param panel_columns Optional integer number of columns in the multi-panel plot layout.
#' @param col Color palette used by the continuous density scale. The default
#'   progresses from near-white at zero density to deep blue at high density,
#'   so visually bold color indicates greater concentration.
#'
#' @return A \code{ggplot} object. The reconstruction inputs and evaluated KDEs
#'   are stored in \code{attr(plot, "npt_density_data")} as a list containing:
#' \describe{
#'   \item{\code{variables}}{Names of the two visualized response variables.}
#'   \item{\code{prediction_indices}}{Row indices of the plotted prediction points.}
#'   \item{\code{prediction_names}}{Labels of the plotted prediction points.}
#'   \item{\code{predictor_values}}{Predictor values corresponding to each panel.}
#'   \item{\code{pair_correlations}}{Latent correlation values \eqn{\hat{\rho}(z)} for each panel.}
#'   \item{\code{samples}}{List of simulated bivariate Monte Carlo samples per panel.}
#'   \item{\code{grid}}{List with \code{x} and \code{y} coordinate evaluation vectors.}
#'   \item{\code{densities}}{List of evaluated \eqn{\text{grid\_size} \times \text{grid\_size}} density matrices.}
#'   \item{\code{bandwidth}}{Length-2 numeric vector of pooled KDE bandwidths.}
#'   \item{\code{limits}}{\eqn{2 \times 2} matrix holding the shared \code{x} and \code{y} axis limits.}
#'   \item{\code{density_limit}}{Length-2 vector \code{c(0, max_density)} defining the shared color scale.}
#' }
#'
#' @seealso \code{\link{npt_frechetreg}}, \code{\link{plot_npt_correlation}}
#' @importFrom rlang .data
#' @export
plot_npt_density <- function(
  fit,
  variables = NULL,
  predictions = NULL,
  n_samples = 4000L,
  seed = NULL,
  grid_size = 100L,
  trim = c(0.01, 0.99),
  bandwidth_adjust = 1.4,
  panel_columns = NULL,
  col = c("#F7FBFF", "#CFE8F3", "#73A9C2", "#2C6C9B", "#08306B")
) {
  # 1. Validate inputs and resolve variable and prediction selections
  .require_npt_fit(fit)

  variable_indices <- .resolve_density_variables(variables, fit$variable_names)
  prediction_indices <- .resolve_density_predictions(
    predictions,
    fit$prediction_names
  )
  n_samples <- .as_integer_count(n_samples, "n_samples", minimum = 2L)
  grid_size <- .as_integer_count(grid_size, "grid_size", minimum = 2L)
  seed <- .validate_optional_seed(seed)

  if (!is.numeric(trim) || length(trim) != 2L || any(!is.finite(trim)) ||
    trim[1L] < 0 || trim[2L] > 1 || trim[1L] >= trim[2L]) {
    stop("`trim` must contain two increasing probabilities in [0, 1].", call. = FALSE)
  }
  if (!is.numeric(bandwidth_adjust) || length(bandwidth_adjust) != 1L ||
    !is.finite(bandwidth_adjust) || bandwidth_adjust <= 0) {
    stop("`bandwidth_adjust` must be one finite positive number.", call. = FALSE)
  }
  if (!is.null(panel_columns)) {
    panel_columns <- .as_integer_count(
      panel_columns,
      "panel_columns",
      minimum = 1L
    )
  }

  # 2. Manage random seed to ensure reproducibility without side effects
  if (!is.null(seed)) {
    had_random_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_random_seed) {
      old_random_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    on.exit(
      {
        if (had_random_seed) {
          assign(".Random.seed", old_random_seed, envir = .GlobalEnv)
        } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      },
      add = TRUE
    )
    set.seed(seed)
  }

  variable_names <- fit$variable_names[variable_indices]
  prediction_names <- fit$prediction_names[prediction_indices]

  # 3. Simulate bivariate observational samples for each prediction point
  samples <- lapply(prediction_indices, function(prediction) {
    pair_correlation <- fit$correlations[[prediction]][
      variable_indices,
      variable_indices,
      drop = FALSE
    ]
    root <- .symmetric_psd_root(pair_correlation)
    latent_standard_normal <- matrix(
      stats::rnorm(n_samples * 2L),
      nrow = n_samples,
      ncol = 2L
    )
    latent_sample <- latent_standard_normal %*% t(root)
    probabilities <- stats::pnorm(latent_sample)

    # Transform uniform probabilities to observational values via fitted marginal quantiles
    observed_sample <- vapply(seq_len(2L), function(variable) {
      stats::approx(
        x = fit$probabilities,
        y = fit$quantiles[[variable_indices[variable]]][prediction, ],
        xout = probabilities[, variable],
        rule = 2,
        ties = "ordered"
      )$y
    }, numeric(n_samples))
    colnames(observed_sample) <- variable_names
    observed_sample
  })
  names(samples) <- prediction_names

  pooled_sample <- do.call(rbind, samples)
  if (any(!is.finite(pooled_sample))) {
    stop("Fitted quantiles produced non-finite visualization samples.", call. = FALSE)
  }

  # 4. Compute pooled display limits, bandwidths, and 2D KDE densities
  limits <- vapply(seq_len(2L), function(variable) {
    display_range <- stats::quantile(
      pooled_sample[, variable],
      probs = trim,
      names = FALSE,
      type = 8
    )
    span <- diff(display_range)
    if (!is.finite(span) || span <= 0) {
      stop(
        sprintf("Variable `%s` has no spread to display.", variable_names[variable]),
        call. = FALSE
      )
    }
    display_range + c(-0.04, 0.04) * span
  }, numeric(2L))
  rownames(limits) <- c("lower", "upper")
  colnames(limits) <- variable_names

  bandwidth <- vapply(
    seq_len(2L),
    function(variable) stats::bw.nrd0(pooled_sample[, variable]),
    numeric(1L)
  ) * bandwidth_adjust
  names(bandwidth) <- variable_names
  if (any(!is.finite(bandwidth)) || any(bandwidth <= 0)) {
    stop("The pooled samples have insufficient spread for a bivariate KDE.", call. = FALSE)
  }

  grid <- list(
    x = seq(limits[1L, 1L], limits[2L, 1L], length.out = grid_size),
    y = seq(limits[1L, 2L], limits[2L, 2L], length.out = grid_size)
  )
  densities <- lapply(samples, .bivariate_gaussian_kde,
    x = grid$x,
    y = grid$y,
    bandwidth = bandwidth
  )
  density_limit <- c(0, max(vapply(densities, max, numeric(1L))))

  # 5. Assemble tidy grid data and a customizable ggplot visualization
  number_panels <- length(prediction_indices)
  if (is.null(panel_columns)) {
    panel_columns <- min(3L, number_panels)
  }
  pair_correlations <- vapply(seq_len(number_panels), function(panel) {
    prediction <- prediction_indices[panel]
    fit$correlations[[prediction]][
      variable_indices[1L],
      variable_indices[2L]
    ]
  }, numeric(1L))
  names(pair_correlations) <- prediction_names

  panel_labels <- vapply(seq_len(number_panels), function(panel) {
    .density_panel_title(
      fit,
      prediction_indices[panel],
      pair_correlations[panel]
    )
  }, character(1L))

  density_data <- do.call(rbind, lapply(seq_len(number_panels), function(panel) {
    data.frame(
      x = rep(grid$x, times = length(grid$y)),
      y = rep(grid$y, each = length(grid$x)),
      density = as.vector(densities[[panel]]),
      panel = factor(
        panel_labels[panel],
        levels = panel_labels,
        ordered = TRUE
      )
    )
  }))
  rownames(density_data) <- NULL

  density_plot <- ggplot2::ggplot(
    density_data,
    ggplot2::aes(x = .data$x, y = .data$y, fill = .data$density)
  ) +
    ggplot2::geom_raster(interpolate = TRUE) +
    ggplot2::geom_contour(
      mapping = ggplot2::aes(x = .data$x, y = .data$y, z = .data$density),
      breaks = seq(
        0.15 * density_limit[2L],
        0.9 * density_limit[2L],
        length.out = 6L
      ),
      colour = "white",
      linewidth = 0.25,
      alpha = 0.55,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    ggplot2::facet_wrap(stats::as.formula("~ panel"), ncol = panel_columns) +
    ggplot2::scale_fill_gradientn(
      colours = col,
      limits = density_limit,
      name = "Density"
    ) +
    ggplot2::coord_cartesian(
      xlim = limits[, 1L],
      ylim = limits[, 2L],
      expand = FALSE
    ) +
    ggplot2::labs(
      title = "Predicted bivariate distributions",
      subtitle = "Shared coordinate limits, bandwidths, and density scale"
      # x = variable_names[1L],
      # y = variable_names[2L]
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      aspect.ratio = 1,
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        colour = "grey70",
        fill = NA,
        linewidth = 0.4
      ),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )

  # 6. Keep the numerical reconstruction available without obscuring ggplot use
  density_result <- list(
    variables = variable_names,
    prediction_indices = prediction_indices,
    prediction_names = prediction_names,
    predictor_values = fit$Z[prediction_indices, , drop = FALSE],
    pair_correlations = pair_correlations,
    samples = samples,
    grid = grid,
    densities = densities,
    bandwidth = bandwidth,
    limits = limits,
    density_limit = density_limit,
    plot_data = density_data
  )
  attr(density_plot, "npt_density_data") <- density_result
  density_plot
}

#' Plot a Predictor-Indexed In-Sample Latent Correlation Fit
#'
#' Plots observed latent Gaussian correlations for one pair of response variables
#' against a selected predictor and overlays the fitted nonparanormal Frechet
#' regression trajectory.
#'
#' @details
#' The points are latent correlations estimated from the empirical distributional
#' responses supplied to \code{\link{npt_frechetreg}}. The fitted curve evaluates
#' \eqn{\hat R(x)} on an evenly spaced grid over the selected predictor while
#' holding all remaining predictors at their sample means.
#'
#' This diagnostic requires a fit created with \code{Z = NULL}, because it uses
#' the aligned training predictors and observed correlation matrices retained by
#' the in-sample fit.
#'
#' @param fit An in-sample \code{npt_frechetreg} object created by calling
#'   \code{\link{npt_frechetreg}} with \code{Z = NULL}.
#' @param variables A character vector of two variable names or integer indices
#'   specifying the correlation to plot. Defaults to the first two variables
#'   \code{c(1, 2)}.
#' @param predictor A predictor name or integer column index for the horizontal
#'   axis. Defaults to the first predictor.
#' @param grid_size Number of predictor values used to draw the fitted curve
#'   (default is `100L`).
#' @param point_alpha Opacity of the observed correlation points in \eqn{[0, 1]}
#'   (default is `0.25`).
#'
#' @return A \code{ggplot} object. The values are stored in
#'   \code{attr(plot, "npt_correlation_data")} as a list containing the selected
#'   variables, predictor, observed point data, and fitted curve data.
#'
#' @seealso \code{\link{plot_npt_density}}, \code{\link{npt_frechetreg}}
#' @export
plot_npt_correlation <- function(
  fit,
  variables = NULL,
  predictor = 1L,
  grid_size = 100L,
  point_alpha = 0.25
) {
  .require_in_sample_npt_fit(fit)
  if (is.null(fit$observed_correlations) ||
    length(fit$observed_correlations) != length(fit$correlations)) {
    stop("`fit` does not contain aligned observed correlations.", call. = FALSE)
  }
  variable_indices <- .resolve_density_variables(variables, fit$variable_names)
  predictor_index <- .resolve_plot_predictor(predictor, colnames(fit$Z))
  grid_size <- .as_integer_count(grid_size, "grid_size", minimum = 2L)
  if (!is.numeric(point_alpha) || length(point_alpha) != 1L ||
    !is.finite(point_alpha) || point_alpha < 0 || point_alpha > 1) {
    stop("`point_alpha` must be one finite number in [0, 1].", call. = FALSE)
  }

  variable_names <- fit$variable_names[variable_indices]
  predictor_name <- colnames(fit$Z)[predictor_index]
  observed <- vapply(fit$observed_correlations, function(correlation) {
    correlation[variable_indices[1L], variable_indices[2L]]
  }, numeric(1L))
  observed_data <- data.frame(
    distribution = fit$prediction_names,
    predictor = fit$Z[, predictor_index],
    correlation = observed,
    check.names = FALSE
  )

  predictor_grid <- seq(
    min(fit$Z[, predictor_index]),
    max(fit$Z[, predictor_index]),
    length.out = grid_size
  )
  curve_predictors <- matrix(
    rep(colMeans(fit$Z), each = grid_size),
    nrow = grid_size,
    ncol = ncol(fit$Z),
    dimnames = list(NULL, colnames(fit$Z))
  )
  curve_predictors[, predictor_index] <- predictor_grid

  curve_fit <- correlation_regression_cpp(
    fit$observed_correlations,
    fit$Z,
    curve_predictors,
    fit$max_iter,
    fit$tol,
    1e-12
  )
  fitted_curve <- vapply(curve_fit$correlations, function(correlation) {
    correlation[variable_indices[1L], variable_indices[2L]]
  }, numeric(1L))
  curve_data <- data.frame(
    predictor = predictor_grid,
    correlation = fitted_curve
  )

  correlation_plot <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = observed_data,
      mapping = ggplot2::aes(
        x = .data$predictor,
        y = .data$correlation,
        colour = "Observed"
      ),
      size = 1.8,
      alpha = point_alpha
    ) +
    ggplot2::geom_line(
      data = curve_data,
      mapping = ggplot2::aes(
        x = .data$predictor,
        y = .data$correlation,
        colour = "Fitted"
      ),
      linewidth = 1.2
    ) +
    ggplot2::scale_colour_manual(
      name = NULL,
      values = c("Observed" = "grey35", "Fitted" = "#D55E00"),
      breaks = c("Observed", "Fitted"),
      guide = ggplot2::guide_legend(
        override.aes = list(
          shape = c(16, NA),
          linetype = c("blank", "solid"),
          alpha = c(0.7, 1),
          linewidth = c(0, 1.2)
        )
      )
    ) +
    ggplot2::labs(
      subtitle = sprintf("%s vs. %s", variable_names[1L], variable_names[2L]),
      x = predictor_name,
      y = "Latent correlation"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        colour = "grey75",
        fill = NA,
        linewidth = 0.4
      )
    )

  attr(correlation_plot, "npt_correlation_data") <- list(
    variables = variable_names,
    predictor = predictor_name,
    observed = observed_data,
    fitted_curve = curve_data,
    fixed_predictors = colMeans(fit$Z)[-predictor_index]
  )
  correlation_plot
}

# ------------------------------------------------------------------------------
# Internal Helpers for Regression Visualization
# ------------------------------------------------------------------------------

# Validates that an object was produced by npt_frechetreg.
.require_npt_fit <- function(fit) {
  if (!inherits(fit, "npt_frechetreg")) {
    stop("`fit` must be created by `npt_frechetreg()`.", call. = FALSE)
  }
  invisible(fit)
}

# In-sample diagnostics require training observations aligned with fitted responses.
.require_in_sample_npt_fit <- function(fit) {
  .require_npt_fit(fit)
  if (!isTRUE(fit$in_sample)) {
    stop(
      "`fit` must be an in-sample fit created with `npt_frechetreg(..., Z = NULL)`.",
      call. = FALSE
    )
  }
  invisible(fit)
}

# Resolves one predictor used for the horizontal axis of a fitted trajectory.
.resolve_plot_predictor <- function(predictor, predictor_names) {
  if (is.character(predictor) && length(predictor) == 1L &&
    !is.na(predictor) && predictor %in% predictor_names) {
    return(match(predictor, predictor_names))
  }
  if (is.numeric(predictor) && length(predictor) == 1L &&
    is.finite(predictor) && predictor == trunc(predictor) &&
    predictor >= 1L && predictor <= length(predictor_names)) {
    return(as.integer(predictor))
  }
  stop("`predictor` must be one fitted predictor name or column index.", call. = FALSE)
}

# Resolves and validates the pair of response variables to plot
.resolve_density_variables <- function(variables, variable_names) {
  if (length(variable_names) < 2L) {
    stop("At least two fitted response variables are required.", call. = FALSE)
  }
  if (is.null(variables)) {
    return(seq_len(2L))
  }

  if (is.character(variables)) {
    if (length(variables) != 2L || anyNA(variables) ||
      anyDuplicated(variables) || !all(variables %in% variable_names)) {
      stop("`variables` must name two distinct fitted response variables.", call. = FALSE)
    }
    return(match(variables, variable_names))
  }

  if (is.numeric(variables) && length(variables) == 2L &&
    all(is.finite(variables)) && all(variables == trunc(variables)) &&
    !anyDuplicated(variables) && all(variables >= 1L) &&
    all(variables <= length(variable_names))) {
    return(as.integer(variables))
  }

  stop("`variables` must contain two distinct fitted names or indices.", call. = FALSE)
}

# Resolves and validates the prediction points to plot
.resolve_density_predictions <- function(predictions, prediction_names) {
  number_predictions <- length(prediction_names)
  if (number_predictions < 1L) {
    stop("The fitted object contains no prediction rows.", call. = FALSE)
  }
  if (is.null(predictions)) {
    number_selected <- min(5L, number_predictions)
    return(unique(as.integer(round(seq(
      1,
      number_predictions,
      length.out = number_selected
    )))))
  }

  if (is.character(predictions)) {
    if (length(predictions) < 1L || anyNA(predictions) ||
      !all(predictions %in% prediction_names)) {
      stop("Unknown prediction name in `predictions`.", call. = FALSE)
    }
    return(match(predictions, prediction_names))
  }

  if (is.numeric(predictions) && length(predictions) >= 1L &&
    all(is.finite(predictions)) && all(predictions == trunc(predictions)) &&
    all(predictions >= 1L) && all(predictions <= number_predictions)) {
    return(as.integer(predictions))
  }

  stop("`predictions` must contain fitted prediction names or indices.", call. = FALSE)
}

# Computes a symmetric PSD square root factor V * Lambda^{1/2}
.symmetric_psd_root <- function(matrix) {
  symmetric <- (matrix + t(matrix)) / 2
  decomposition <- eigen(symmetric, symmetric = TRUE)
  eigen_scale <- max(1, max(abs(decomposition$values)))
  if (min(decomposition$values) < -1e-8 * eigen_scale) {
    stop("A fitted latent correlation is not positive semidefinite.", call. = FALSE)
  }

  root_values <- sqrt(pmax(decomposition$values, 0))
  sweep(decomposition$vectors, 2L, root_values, `*`)
}

# Evaluates 2D Gaussian kernel density estimate on a regular grid
.bivariate_gaussian_kde <- function(sample, x, y, bandwidth) {
  standardized_x <- outer(x, sample[, 1L], `-`) / bandwidth[1L]
  standardized_y <- outer(y, sample[, 2L], `-`) / bandwidth[2L]
  kernel_x <- exp(-0.5 * standardized_x^2)
  kernel_y <- exp(-0.5 * standardized_y^2)

  (kernel_x %*% t(kernel_y)) /
    (nrow(sample) * 2 * pi * bandwidth[1L] * bandwidth[2L])
}

# Generates informative panel title with predictor value and latent correlation
.density_panel_title <- function(fit, prediction, pair_correlation) {
  if (base::ncol(fit$Z) == 1L) {
    predictor_name <- colnames(fit$Z)[1L]
    if (is.null(predictor_name) || is.na(predictor_name) || predictor_name == "") {
      predictor_name <- "predictor"
    }
    predictor_value <- format(fit$Z[prediction, 1L], digits = 3, trim = TRUE)
    first_line <- sprintf("%s = %s", predictor_name, predictor_value)
  } else if (base::ncol(fit$Z) == 2L) {
    predictor_values <- format(
      fit$Z[prediction, ],
      digits = 2,
      trim = TRUE
    )
    first_line <- sprintf(
      "X = (%s, %s)",
      predictor_values[1L],
      predictor_values[2L]
    )
  } else {
    first_line <- fit$prediction_names[prediction]
  }
  display_correlation <- if (abs(pair_correlation) < 0.005) 0 else pair_correlation
  sprintf("%s\nlatent rho = %.2f", first_line, display_correlation)
}
