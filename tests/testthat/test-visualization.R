make_density_plot_fit <- function() {
  probabilities <- (seq_len(31L) - 0.5) / 31L
  predictor_values <- c(-1, 0, 1)
  prediction_names <- c("low", "middle", "high")
  variable_names <- c("location", "positive")

  location_quantiles <- vapply(predictor_values, function(value) {
    stats::qnorm(probabilities) + 0.8 * value
  }, numeric(length(probabilities)))
  positive_quantiles <- vapply(predictor_values, function(value) {
    exp(0.35 * stats::qnorm(probabilities) + 0.25 * value)
  }, numeric(length(probabilities)))

  correlations <- lapply(c(-0.55, 0, 0.55), function(rho) {
    matrix(
      c(1, rho, rho, 1),
      nrow = 2L,
      dimnames = list(variable_names, variable_names)
    )
  })
  names(correlations) <- prediction_names
  observed_correlations <- lapply(c(-0.7, 0.1, 0.4), function(rho) {
    matrix(
      c(1, rho, rho, 1),
      nrow = 2L,
      dimnames = list(variable_names, variable_names)
    )
  })
  names(observed_correlations) <- prediction_names

  structure(
    list(
      quantiles = list(
        location = t(location_quantiles),
        positive = t(positive_quantiles)
      ),
      correlations = correlations,
      probabilities = probabilities,
      Z = matrix(
        predictor_values,
        ncol = 1L,
        dimnames = list(prediction_names, "exposure")
      ),
      prediction_names = prediction_names,
      variable_names = variable_names,
      in_sample = TRUE,
      observed_correlations = observed_correlations,
      max_iter = 1000L,
      tol = 1e-6
    ),
    class = "npt_frechetreg"
  )
}

test_that("plot_npt_density returns a reproducible ggplot density visualization", {
  fit <- make_density_plot_fit()
  set.seed(90210)
  random_state_before_plot <- .Random.seed

  first_plot <- plot_npt_density(
    fit,
    variables = c("location", "positive"),
    predictions = c("low", "high"),
    n_samples = 250L,
    seed = 17L,
    grid_size = 30L,
    panel_columns = 2L
  )
  first <- attr(first_plot, "npt_density_data")

  expect_s3_class(first_plot, "ggplot")
  expect_silent(ggplot2::ggplot_build(first_plot))
  density_scale <- first_plot$scales$get_scales("fill")
  expect_identical(density_scale$palette(0), "#F7FBFF")
  expect_identical(density_scale$palette(1), "#08306B")
  expect_identical(.Random.seed, random_state_before_plot)
  expect_identical(first$variables, c("location", "positive"))
  expect_identical(first$prediction_names, c("low", "high"))
  expect_equal(unname(vapply(first$samples, dim, integer(2L))), matrix(
    c(250L, 2L, 250L, 2L),
    nrow = 2L
  ))
  expect_true(all(vapply(first$samples, function(x) all(is.finite(x)), logical(1L))))
  expect_true(all(vapply(first$densities, function(x) all(is.finite(x)), logical(1L))))
  expect_true(all(vapply(first$densities, function(x) all(x >= 0), logical(1L))))
  expect_equal(unname(vapply(first$densities, dim, integer(2L))), matrix(
    rep(c(30L, 30L), 2L),
    nrow = 2L
  ))

  second_plot <- plot_npt_density(
    fit,
    variables = c(1L, 2L),
    predictions = c(1L, 3L),
    n_samples = 250L,
    seed = 17L,
    grid_size = 30L,
    panel_columns = 2L
  )
  second <- attr(second_plot, "npt_density_data")
  expect_equal(second$samples, first$samples)
  expect_equal(second$densities, first$densities)
})

test_that("plot_npt_density validates pair and prediction selection", {
  fit <- make_density_plot_fit()

  expect_error(
    plot_npt_density(fit, variables = c("location", "location")),
    "two distinct"
  )
  expect_error(
    plot_npt_density(fit, predictions = "unknown"),
    "Unknown prediction"
  )
  expect_error(
    plot_npt_density(list(), variables = c(1L, 2L)),
    "npt_frechetreg"
  )

  out_of_sample_fit <- fit
  out_of_sample_fit$in_sample <- FALSE
  out_of_sample_fit$observed_correlations <- NULL

  # plot_npt_density supports out-of-sample prediction objects
  out_of_sample_plot <- plot_npt_density(out_of_sample_fit, predictions = 1L)
  expect_s3_class(out_of_sample_plot, "ggplot")

  # plot_npt_correlation requires in-sample fits to overlay observed correlations
  expect_error(
    plot_npt_correlation(out_of_sample_fit),
    "Z = NULL",
    fixed = TRUE
  )
})

test_that("plot_npt_density accepts a positive-semidefinite boundary fit", {
  fit <- make_density_plot_fit()
  fit$correlations[[1L]][1L, 2L] <- 1
  fit$correlations[[1L]][2L, 1L] <- 1

  density_plot <- plot_npt_density(
    fit,
    predictions = 1L,
    n_samples = 100L,
    seed = 5L,
    grid_size = 25L
  )
  result <- attr(density_plot, "npt_density_data")

  # At rho = 1 both latent coordinates are identical. Monotone quantile maps
  # preserve their ranks even though the two observed margins have different
  # shapes and scales.
  expect_equal(
    stats::cor(result$samples[[1L]], method = "spearman")[1L, 2L],
    1,
    tolerance = 1e-12
  )
})

test_that("plot_npt_correlation overlays a predictor-indexed fitted curve", {
  fit <- make_density_plot_fit()
  correlation_plot <- plot_npt_correlation(
    fit,
    variables = c("location", "positive"),
    predictor = "exposure",
    grid_size = 25L,
    point_alpha = 0.2
  )
  result <- attr(correlation_plot, "npt_correlation_data")

  expect_s3_class(correlation_plot, "ggplot")
  expect_silent(ggplot2::ggplot_build(correlation_plot))
  expect_identical(result$variables, c("location", "positive"))
  expect_identical(result$predictor, "exposure")
  expect_equal(result$observed$predictor, c(-1, 0, 1))
  expect_equal(result$observed$correlation, c(-0.7, 0.1, 0.4))
  expect_equal(nrow(result$fitted_curve), 25L)
  expect_equal(range(result$fitted_curve$predictor), c(-1, 1))
  expect_true(all(is.finite(result$fitted_curve$correlation)))
})

test_that("plot_npt_correlation requires an in-sample fit and two variables", {
  fit <- make_density_plot_fit()
  expect_error(
    plot_npt_correlation(fit, variables = c("location", "location")),
    "two distinct"
  )
  expect_error(
    plot_npt_correlation(fit, predictor = "unknown"),
    "fitted predictor"
  )
  expect_error(
    plot_npt_correlation(fit, point_alpha = 2),
    "point_alpha"
  )

  fit$in_sample <- FALSE
  fit$observed_correlations <- NULL
  expect_error(
    plot_npt_correlation(fit),
    "Z = NULL",
    fixed = TRUE
  )
})
