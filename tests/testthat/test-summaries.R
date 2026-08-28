# These tests follow the summary pipeline in mathematical order: tau-a,
# quantiles, correlation constraints, and finally the public input contract.

test_that("compiled Kendall tau-a agrees with its sign definition", {
  examples <- list(
    list(x = c(1, 2, 3), y = c(0, 0, 1)),
    list(x = c(0, 0, 1), y = c(2, 1, 3)),
    list(x = c(1, 1, 2, 2, 3, 3), y = c(4, 4, 1, 2, 3, 3)),
    list(x = c(1, 2, 3, 4, 5), y = c(5, 4, 3, 2, 1))
  )

  for (example in examples) {
    observed <- NonparaTransport:::kendall_tau_a_cpp(example$x, example$y)
    expect_equal(observed, brute_tau_a(example$x, example$y), tolerance = 1e-14)
  }
})

test_that("Kendall tau-a is symmetric and invariant to row order", {
  x <- c(0, 0, 1, 2, 2, 3, 4)
  y <- c(3, 1, 2, 2, 5, 4, 6)
  permutation <- c(6, 2, 7, 1, 5, 3, 4)

  tau_xy <- NonparaTransport:::kendall_tau_a_cpp(x, y)
  expect_equal(
    tau_xy,
    NonparaTransport:::kendall_tau_a_cpp(y, x),
    tolerance = 1e-14
  )
  expect_equal(
    tau_xy,
    NonparaTransport:::kendall_tau_a_cpp(x[permutation], y[permutation]),
    tolerance = 1e-14
  )
})

test_that("compiled Kendall tau-a matches random tied samples", {
  set.seed(2718)
  for (iteration in seq_len(50)) {
    sample_size <- sample(4:15, 1)
    x <- sample(-2:3, sample_size, replace = TRUE)
    y <- sample(-1:4, sample_size, replace = TRUE)
    expect_equal(
      NonparaTransport:::kendall_tau_a_cpp(x, y),
      brute_tau_a(x, y),
      tolerance = 1e-14
    )
  }
})

test_that("as_nonparanormal uses type-7 quantiles on equal-interval midpoints", {
  data <- example_distributions()
  summaries <- as_nonparanormal(data, M = 11, cache_sqrt = TRUE)
  probabilities <- (seq_len(11) - 0.5) / 11

  expect_s3_class(summaries, "nonparanormal")
  expect_equal(summaries$probabilities, probabilities, tolerance = 1e-15)
  expect_equal(names(summaries$sample_sizes), names(data))
  expect_equal(
    dimnames(summaries$correlations[[1L]]),
    list(colnames(data[[1L]]), colnames(data[[1L]]))
  )
  # Compare every row of the stored n-by-M matrices against R's defining
  # type-7 implementation,
  # not against another approximation written specifically for the test.
  for (distribution in seq_along(data)) {
    for (variable in seq_len(ncol(data[[distribution]]))) {
      expected <- unname(quantile(
        data[[distribution]][, variable],
        probs = probabilities,
        type = 7
      ))
      expect_equal(
        unname(summaries$quantiles[[variable]][distribution, ]),
        expected,
        tolerance = 1e-14
      )
    }
  }
})

test_that("pooled median and raw MAD standardization matches direct preprocessing", {
  data <- example_distributions()
  pooled <- do.call(rbind, data)
  centers <- apply(pooled, 2L, median)
  scales <- vapply(
    seq_len(ncol(pooled)),
    function(variable) median(abs(pooled[, variable] - centers[variable])),
    numeric(1)
  )
  names(centers) <- colnames(pooled)
  names(scales) <- colnames(pooled)

  manually_standardized <- lapply(data, function(distribution) {
    sweep(sweep(distribution, 2L, centers, "-"), 2L, scales, "/")
  })
  observed <- as_nonparanormal(data, M = 13L, standardize = TRUE)
  expected <- as_nonparanormal(
    manually_standardized,
    M = 13L,
    standardize = FALSE
  )

  expect_true(observed$standardization$applied)
  expect_identical(observed$standardization$method, "pooled_median_mad")
  expect_equal(observed$standardization$center, centers, tolerance = 1e-14)
  expect_equal(observed$standardization$scale, scales, tolerance = 1e-14)
  expect_equal(observed$quantiles, expected$quantiles, tolerance = 1e-13)
  expect_equal(observed$correlations, expected$correlations, tolerance = 1e-13)

  unscaled <- as_nonparanormal(data, M = 13L)
  expect_false(unscaled$standardization$applied)
  expect_identical(unscaled$standardization$method, "none")
})

test_that("standardization rejects a zero pooled MAD without affecting raw summaries", {
  data <- list(
    first = cbind(x = c(0, 0, 0, 1), y = c(0, 0, 0, 2)),
    second = cbind(x = c(0, 0, 0, -1), y = c(0, 0, 0, -2))
  )

  expect_s3_class(
    as_nonparanormal(data, M = 7L, standardize = FALSE),
    "nonparanormal"
  )
  expect_error(
    as_nonparanormal(data, M = 7L, standardize = TRUE),
    "pooled MAD is zero"
  )
})

test_that("tied type-7 quantiles remain monotone at machine precision", {
  data <- list(tied = cbind(
    first = c(0, 1 / 3, 1 / 3, 2 / 3),
    second = c(0, 1, 2, 3)
  ))
  summaries <- as_nonparanormal(data, M = 16, cache_sqrt = FALSE)

  expect_true(all(diff(summaries$quantiles$first[1, ]) >= 0))
})

test_that("summary consumers reject a non-midpoint probability grid", {
  summaries <- as_nonparanormal(example_distributions(), M = 7)
  summaries$probabilities <- seq(0, 1, length.out = 7)

  expect_error(
    pairwise_npt_distance(summaries),
    "equal-interval midpoint grid"
  )
})

test_that("projected latent correlations are strict correlation matrices", {
  shrinkage <- 0.01
  summaries <- as_nonparanormal(
    example_distributions(),
    M = 9,
    pd_shrinkage = shrinkage,
    cache_sqrt = TRUE
  )

  # The projection/shrinkage contract has three parts: symmetry, unit diagonal,
  # and eigenvalues bounded below by the requested shrinkage. Cached roots must
  # also square back to the same matrices.
  for (i in seq_along(summaries$correlations)) {
    correlation <- summaries$correlations[[i]]
    square_root <- summaries$correlation_sqrts[[i]]

    expect_equal(correlation, t(correlation), tolerance = 1e-12)
    expect_equal(unname(diag(correlation)), rep(1, nrow(correlation)), tolerance = 1e-12)
    expect_gte(min(eigen(correlation, symmetric = TRUE, only.values = TRUE)$values), shrinkage - 1e-10)
    expect_equal(square_root %*% square_root, correlation, tolerance = 1e-10)
  }
})

test_that("an indefinite pairwise estimate is projected to the correlation cone", {
  # This fixed example was chosen so the independently assembled pairwise
  # bridge estimates are indefinite; it forces the nontrivial Dykstra branch.
  data <- matrix(
    c(
      3, 1, 2, 4, 1, 3,
      2, 1, 4, 1, 1, 3,
      3, 4, 4, 4, 3, 2,
      4, 3, 4, 4, 4, 2,
      3, 4, 1, 2, 4, 4,
      3, 3, 2, 3, 4, 4
    ),
    nrow = 6,
    ncol = 6
  )

  # Reconstruct the raw sine-bridge matrix in transparent R code and verify the
  # negative eigenvalue before asking the package to repair it.
  raw <- diag(ncol(data))
  for (i in 2:ncol(data)) {
    for (j in seq_len(i - 1L)) {
      tau <- brute_tau_a(data[, i], data[, j])
      raw[i, j] <- raw[j, i] <- sin(pi * tau / 2)
    }
  }
  expect_lt(min(eigen(raw, symmetric = TRUE, only.values = TRUE)$values), 0)

  shrinkage <- 0.005
  projected <- as_nonparanormal(
    list(example = data),
    M = 7,
    pd_shrinkage = shrinkage
  )$correlations[[1L]]
  expect_equal(unname(diag(projected)), rep(1, ncol(data)), tolerance = 1e-12)
  expect_gte(
    min(eigen(projected, symmetric = TRUE, only.values = TRUE)$values),
    shrinkage - 1e-9
  )
})

test_that("R validation rejects malformed or unidentified inputs", {
  expect_error(as_nonparanormal(matrix(1:4, 2)), "non-empty list")
  expect_error(
    as_nonparanormal(list(matrix(c(1, NA, 2, 3), ncol = 2))),
    "finite"
  )
  expect_error(
    as_nonparanormal(list(a = cbind(x = 1:3, y = rep(1, 3)))),
    "not identifiable"
  )
  expect_error(
    as_nonparanormal(list(matrix(1:2, nrow = 1))),
    "at least two rows"
  )
  expect_error(
    as_nonparanormal(list(matrix(1:4, ncol = 1))),
    "at least 2"
  )
  expect_error(
    as_nonparanormal(example_distributions(), standardize = NA),
    "`standardize` must be TRUE or FALSE"
  )
})
