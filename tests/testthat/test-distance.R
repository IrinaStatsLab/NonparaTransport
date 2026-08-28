# These tests verify the distance pipeline from cache handling through the two
# mathematical components and finally against a direct R implementation.

test_that("cached and on-demand matrix square roots give the same distances", {
  data <- example_distributions()
  cached <- as_nonparanormal(data, M = 21, cache_sqrt = TRUE)
  uncached <- as_nonparanormal(data, M = 21, cache_sqrt = FALSE)

  expect_equal(
    pairwise_npt_distance(cached, decompose = TRUE),
    pairwise_npt_distance(uncached, decompose = TRUE),
    tolerance = 1e-12
  )
})

test_that("cached roots agree with on-demand roots above dimension two", {
  # The bivariate implementation uses a scalar formula, so d=3 is required to
  # exercise the general matrix-square-root cache path.
  set.seed(42)
  data <- list(
    first = matrix(rnorm(36), ncol = 3),
    second = matrix(rnorm(45), ncol = 3),
    third = matrix(rnorm(54), ncol = 3)
  )
  cached <- as_nonparanormal(data, M = 19, cache_sqrt = TRUE)
  uncached <- as_nonparanormal(data, M = 19, cache_sqrt = FALSE)

  expect_equal(
    pairwise_npt_distance(cached, decompose = TRUE),
    pairwise_npt_distance(uncached, decompose = TRUE),
    tolerance = 1e-11
  )
})

test_that("distance decomposition is symmetric and additive", {
  summaries <- as_nonparanormal(
    example_distributions(),
    M = 21,
    cache_sqrt = TRUE
  )
  distances <- pairwise_npt_distance(summaries, decompose = TRUE)

  expect_equal(distances$marginal, t(distances$marginal), tolerance = 1e-14)
  expect_equal(distances$correlation, t(distances$correlation), tolerance = 1e-14)
  expect_equal(unname(diag(distances$total)), c(0, 0), tolerance = 1e-14)
  expect_equal(
    distances$total,
    distances$marginal + distances$correlation,
    tolerance = 1e-14
  )
  expect_equal(pairwise_npt_distance(summaries), distances$total, tolerance = 1e-14)
})

test_that("BW weight changes only the correlation component", {
  summaries <- as_nonparanormal(
    example_distributions(),
    M = 21L,
    standardize = TRUE
  )
  unweighted <- pairwise_npt_distance(
    summaries,
    decompose = TRUE,
    bw_weight = 1
  )
  weighted <- pairwise_npt_distance(
    summaries,
    decompose = TRUE,
    bw_weight = 2.5
  )

  expect_equal(weighted$marginal, unweighted$marginal, tolerance = 1e-14)
  expect_equal(
    weighted$correlation,
    2.5 * unweighted$correlation,
    tolerance = 1e-14
  )
  expect_equal(
    weighted$total,
    weighted$marginal + weighted$correlation,
    tolerance = 1e-14
  )
  expect_error(pairwise_npt_distance(summaries, bw_weight = 0), "positive")
  expect_error(pairwise_npt_distance(summaries, bw_weight = Inf), "positive")
})

test_that("single distance matches the corresponding pairwise matrix entry", {
  data <- c(
    example_distributions(),
    list(third = example_distributions()[[1L]] + 0.25)
  )
  names(data)[1:2] <- c("first", "second")
  summaries <- as_nonparanormal(data, M = 19L, standardize = TRUE)
  pairwise <- pairwise_npt_distance(
    summaries,
    decompose = TRUE,
    bw_weight = 1.75
  )

  single <- npt_distance(
    summaries,
    "first",
    "third",
    decompose = TRUE,
    bw_weight = 1.75
  )
  expect_equal(single$marginal, pairwise$marginal["first", "third"])
  expect_equal(single$correlation, pairwise$correlation["first", "third"])
  expect_equal(single$total, pairwise$total["first", "third"])
  expect_equal(
    npt_distance(summaries, 1L, 3L, bw_weight = 1.75),
    pairwise$total[1L, 3L]
  )
  expect_error(npt_distance(summaries, "missing", "third"), "does not match")
  expect_error(npt_distance(summaries, 0L, 2L), "between 1 and 3")
})

test_that("bivariate BW term agrees with its closed form", {
  summaries <- as_nonparanormal(
    example_distributions(),
    M = 15,
    pd_shrinkage = 0.001,
    cache_sqrt = TRUE
  )
  correlations <- vapply(
    summaries$correlations,
    function(matrix) matrix[1, 2],
    numeric(1)
  )
  # For [[1,rho],[rho,1]], both matrices share eigenvectors and the BW formula
  # reduces to products of square roots of the eigenvalues 1+rho and 1-rho.
  expected <- unname(4 - 2 * (
    sqrt((1 + correlations[1]) * (1 + correlations[2])) +
      sqrt((1 - correlations[1]) * (1 - correlations[2]))
  ))

  observed <- pairwise_npt_distance(summaries, decompose = TRUE)$correlation[1, 2]
  expect_equal(observed, expected, tolerance = 1e-10)
})

test_that("identical summaries have zero pairwise distance", {
  distribution <- example_distributions()[[1L]]
  summaries <- as_nonparanormal(
    list(first = distribution, copy = distribution),
    M = 17,
    cache_sqrt = TRUE
  )

  expect_equal(pairwise_npt_distance(summaries)[1, 2], 0, tolerance = 1e-12)
})

test_that("marginal distances remain accurate after a large common shift", {
  set.seed(91)
  shifts <- c(-0.4, 0, 0.2, 0.7)
  common_auxiliary <- rnorm(60)
  data <- lapply(shifts, function(shift) {
    cbind(
      signal = 1e8 + shift + rnorm(60),
      auxiliary = common_auxiliary
    )
  })
  names(data) <- paste0("distribution_", seq_along(data))
  summaries <- as_nonparanormal(data, M = 37L, cache_sqrt = FALSE)
  observed <- pairwise_npt_distance(summaries, decompose = TRUE)$marginal

  # Direct row differences avoid the cancellation that an uncentered expansion
  # of squared norms would incur at this common location.
  expected <- matrix(0, length(data), length(data))
  for (i in seq_len(length(data) - 1L)) {
    for (k in (i + 1L):length(data)) {
      expected[i, k] <- sum(vapply(
        summaries$quantiles,
        function(quantiles) mean((quantiles[i, ] - quantiles[k, ])^2),
        numeric(1)
      ))
      expected[k, i] <- expected[i, k]
    }
  }

  expect_equal(unname(observed), expected, tolerance = 1e-12)
})

test_that("compiled full matrix agrees with a direct R definition", {
  # This is the end-to-end numerical oracle: unlike the focused tests above, it
  # independently reconstructs every unordered pair and both NPT components.
  set.seed(2026)
  data <- list(
    first = matrix(rnorm(30), ncol = 3),
    second = matrix(rnorm(36), ncol = 3),
    third = matrix(rnorm(42), ncol = 3)
  )
  summaries <- as_nonparanormal(data, M = 17, cache_sqrt = TRUE)
  observed <- pairwise_npt_distance(summaries, decompose = TRUE)

  # Direct R versions favor transparency over speed.
  midpoint_integral <- function(values) mean(values)
  symmetric_sqrt <- function(matrix) {
    decomposition <- eigen((matrix + t(matrix)) / 2, symmetric = TRUE)
    decomposition$vectors %*%
      diag(sqrt(pmax(decomposition$values, 0)), nrow(matrix)) %*%
      t(decomposition$vectors)
  }

  expected_marginal <- matrix(0, 3, 3)
  expected_correlation <- matrix(0, 3, 3)

  # Compute one triangle using the defining formulas, then mirror it exactly as
  # the compiled implementation does.
  for (i in 1:2) {
    for (j in (i + 1L):3) {
      expected_marginal[i, j] <- sum(vapply(
        summaries$quantiles,
        function(quantiles) midpoint_integral((quantiles[i, ] - quantiles[j, ])^2),
        numeric(1)
      ))
      first <- summaries$correlations[[i]]
      second <- summaries$correlations[[j]]
      middle_root <- symmetric_sqrt(
        summaries$correlation_sqrts[[i]] %*% second %*%
          summaries$correlation_sqrts[[i]]
      )
      expected_correlation[i, j] <-
        sum(diag(first)) + sum(diag(second)) - 2 * sum(diag(middle_root))
    }
  }
  expected_marginal <- expected_marginal + t(expected_marginal)
  expected_correlation <- expected_correlation + t(expected_correlation)

  expect_equal(unname(observed$marginal), expected_marginal, tolerance = 1e-11)
  expect_equal(unname(observed$correlation), expected_correlation, tolerance = 1e-10)
  expect_equal(
    unname(observed$total),
    expected_marginal + expected_correlation,
    tolerance = 1e-10
  )
})
