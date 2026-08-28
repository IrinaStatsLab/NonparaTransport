test_that("component-wise R-squared is separate and has no scalar total", {
  example <- make_regression_example(d = 2L)
  result <- npt_component_r2(example$X, example$Y)

  expect_identical(
    result$component,
    c("response_1", "response_2", "latent_correlation")
  )
  expect_false("total" %in% result$component)
  expect_false("adjusted_r_squared" %in% names(result))
  expect_equal(
    result$r_squared,
    1 - result$residual_sum / result$total_sum,
    tolerance = 1e-12
  )
})

test_that("marginal component losses use equal midpoint-cell weights", {
  component_losses <- getFromNamespace(
    "matched_component_losses_cpp",
    "NonparaTransport"
  )
  observed <- matrix(c(0, 2, 5), nrow = 1L)
  predicted <- matrix(c(1, 0, 2), nrow = 1L)

  result <- component_losses(
    observed_quantiles = list(observed),
    predicted_quantiles = list(predicted),
    observed_correlations = list(matrix(1, 1L, 1L)),
    predicted_correlations = list(matrix(1, 1L, 1L)),
    eigen_floor = 1e-15
  )

  expect_equal(unname(result$marginal), mean((observed - predicted)^2))
  expect_equal(result$correlation, 0)
})

test_that("component-wise R-squared is undefined for a constant response", {
  set.seed(91)
  one_distribution <- matrix(rnorm(90), ncol = 3L)
  colnames(one_distribution) <- c("A", "B", "C")
  distributions <- rep(list(one_distribution), 6L)
  names(distributions) <- paste0("distribution_", seq_len(6L))
  X <- matrix(
    seq_len(6L),
    ncol = 1L,
    dimnames = list(names(distributions), "x")
  )
  summaries <- as_nonparanormal(distributions, M = 20L)

  result <- npt_component_r2(X, summaries)
  expect_true(all(result$total_sum < 1e-12))
  expect_true(all(is.na(result$r_squared)))
})

test_that("permutation inference returns component-wise null distributions", {
  example <- make_regression_example(d = 2L, n = 6L, sample_size = 25L)
  result <- npt_permutation_test(
    example$X,
    example$Y,
    B = 50L,
    workers = 1L,
    seed = 17L
  )

  expect_s3_class(result, "npt_permutation_test")
  expect_equal(dim(result$null_distribution), c(50L, 3L))
  expect_identical(colnames(result$null_distribution), result$results$component)
  expect_false("total" %in% result$results$component)
  expect_true(all(result$results$p_value >= 1 / 51))
  expect_true(all(result$results$p_value_adjusted >= result$results$p_value))
  expect_true(all(result$results$p_value_adjusted <= 1))
})

test_that("single-step min-p uses one common reference set with conservative ties", {
  null <- cbind(
    marginal = c(0.4, 0.4, 0.2, 0.1),
    latent_correlation = c(0.5, 0.3, 0.3, 0.1)
  )
  observed <- c(marginal = 0.35, latent_correlation = 0.25)
  min_p <- getFromNamespace(".westfall_young_min_p", "NonparaTransport")

  result <- min_p(null, observed)

  expect_equal(result$unadjusted, c(marginal = 0.6, latent_correlation = 0.8))
  expect_equal(result$adjusted, c(marginal = 0.8, latent_correlation = 0.8))
})

test_that("an extreme component cannot bypass the min-p family adjustment", {
  # The observed row is largest for component_a, while the final permutation
  # row is largest for component_b. Both rows therefore attain the reference
  # minimum 1 / 21, so component_a must have adjusted p-value 2 / 21.
  null <- cbind(
    component_a = 20:1,
    component_b = 1:20
  )
  observed <- c(component_a = 21, component_b = 0)
  min_p <- getFromNamespace(".westfall_young_min_p", "NonparaTransport")

  result <- min_p(null, observed)

  expect_equal(result$unadjusted, c(component_a = 1 / 21, component_b = 1))
  expect_equal(result$adjusted, c(component_a = 2 / 21, component_b = 1))
})

test_that("pre-generated permutations are invariant to worker scheduling", {
  example <- make_regression_example(d = 2L, n = 6L, sample_size = 20L)

  sequential <- npt_permutation_test(
    example$X,
    example$Y,
    B = 6L,
    workers = 1L,
    seed = 73L
  )
  parallel <- npt_permutation_test(
    example$X,
    example$Y,
    B = 6L,
    workers = 2L,
    seed = 73L
  )

  expect_identical(parallel$null_distribution, sequential$null_distribution)
  expect_identical(parallel$results, sequential$results)
})
