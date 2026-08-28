test_that("npt_frechetreg follows the fastfrechet X-Y-Z contract", {
  example <- make_regression_example(d = 2L)
  Z <- matrix(
    c(-0.75, 0, 0.75),
    ncol = 1L,
    dimnames = list(c("low", "middle", "high"), "biomarker")
  )
  second_range <- range(example$Y$quantiles[[2L]]) + c(-0.1, 0.1)
  bounds <- list(response_1 = c(-Inf, Inf), response_2 = second_range)

  fit <- npt_frechetreg(example$X, example$Y, Z, bounds = bounds)

  expect_s3_class(fit, "npt_frechetreg")
  expect_identical(names(fit$quantiles), example$Y$variable_names)
  expect_equal(dim(fit$quantiles[[1L]]), c(3L, 25L))
  expect_identical(rownames(fit$quantiles[[1L]]), rownames(Z))
  expect_identical(names(fit$correlations), rownames(Z))
  expect_equal(fit$Z, Z)
  expect_false(fit$in_sample)
  expect_null(fit$observed_correlations)
  expect_true(all(fit$correlation_diagnostics$status == "closed_form"))

  direct <- fastfrechet::frechetreg_univar2wass(
    X = example$X,
    Y = example$Y$quantiles[[1L]],
    Z = Z
  )
  expect_equal(unname(fit$quantiles[[1L]]), unname(direct$Qhat), tolerance = 1e-12)
  expect_true(all(apply(fit$quantiles[[2L]], 1L, function(x) all(diff(x) >= -1e-12))))
  expect_true(all(fit$quantiles[[2L]] >= second_range[1L] - 1e-12))
  expect_true(all(fit$quantiles[[2L]] <= second_range[2L] + 1e-12))
})

test_that("Z = NULL returns aligned fitted values", {
  example <- make_regression_example(d = 3L)
  fit <- npt_frechetreg(example$X, example$Y, max_iter = 100L)

  expect_equal(dim(fit$quantiles[[1L]]), c(7L, 25L))
  expect_identical(fit$prediction_names, example$Y$distribution_names)
  expect_identical(rownames(fit$Z), example$Y$distribution_names)
  expect_true(fit$in_sample)
  expect_equal(fit$observed_correlations, example$Y$correlations)
  expect_true(all(vapply(
    fit$correlations,
    function(matrix) {
      isTRUE(all.equal(matrix, t(matrix), tolerance = 1e-10)) &&
        max(abs(diag(matrix) - 1)) < 1e-10 &&
        min(eigen(matrix, symmetric = TRUE, only.values = TRUE)$values) > -1e-8
    },
    logical(1)
  )))
  expect_true(all(fit$correlation_diagnostics$converged))
})

test_that("regression validates predictor alignment and dimensions", {
  example <- make_regression_example(d = 2L)
  reversed <- example$X[nrow(example$X):1L, , drop = FALSE]

  aligned_fit <- npt_frechetreg(reversed, example$Y)
  reference_fit <- npt_frechetreg(example$X, example$Y)
  expect_equal(aligned_fit$quantiles, reference_fit$quantiles)

  bad_names <- example$X
  rownames(bad_names)[1L] <- "unknown_subject"
  expect_error(
    npt_frechetreg(bad_names, example$Y),
    "row names of `X`"
  )
  expect_error(
    npt_frechetreg(example$X[-1L, , drop = FALSE], example$Y),
    "must have 7 rows"
  )
  expect_error(
    npt_frechetreg(example$X, example$Y, Z = matrix(1:4, ncol = 2L)),
    "must have 1 predictor columns"
  )
})

test_that("data-frame X and Z follow the same contract as matrices", {
  example <- make_regression_example(d = 2L)
  X_data_frame <- data.frame(biomarker = example$X[, 1L])
  Z_data_frame <- data.frame(biomarker = c(-0.5, 0.5))

  data_frame_fit <- npt_frechetreg(X_data_frame, example$Y, Z_data_frame)
  matrix_fit <- npt_frechetreg(
    unname(example$X),
    example$Y,
    matrix(c(-0.5, 0.5), ncol = 1L)
  )

  expect_equal(
    unname(data_frame_fit$quantiles[[1L]]),
    unname(matrix_fit$quantiles[[1L]])
  )
  expect_identical(colnames(data_frame_fit$Z), "biomarker")
  expect_identical(
    rownames(data_frame_fit$Z),
    c("prediction_1", "prediction_2")
  )
})

test_that("bivariate extreme extrapolation selects the correct boundary", {
  rho <- c(-0.99, 0, 0.99)
  correlations <- lapply(rho, function(value) {
    matrix(c(1, value, value, 1), nrow = 2L)
  })
  names(correlations) <- paste0("distribution_", 1:3)
  common_quantiles <- matrix(rep(seq(-1, 1, length.out = 5L), 3L), nrow = 3L, byrow = TRUE)

  summaries <- structure(
    list(
      probabilities = (seq_len(5L) - 0.5) / 5L,
      quantiles = list(A = common_quantiles, B = common_quantiles),
      correlations = correlations,
      distribution_names = names(correlations),
      variable_names = c("A", "B")
    ),
    class = "nonparanormal"
  )
  X <- matrix(
    c(0, 0, 1, 0, 0, 1),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(names(correlations), c("x1", "x2"))
  )

  # At z = (-4, 2.6), the saturated global-regression weights are
  # (2.4, -4, 2.6). Both bivariate objective coefficients are negative,
  # with the sqrt(1 + rho) coefficient larger, so rho = +1 minimizes the
  # objective among the two boundary endpoints.
  fit <- npt_frechetreg(X, summaries, Z = matrix(c(-4, 2.6), nrow = 1L))
  expect_equal(fit$correlations[[1L]][1L, 2L], 1, tolerance = 1e-12)
})

test_that("correlation regression agrees with the Python core reference", {
  X <- matrix(
    c(-1, -0.5, -0.5, 0.75, 0, -1, 0.5, 0.25, 1, 0.5),
    ncol = 2L,
    byrow = TRUE
  )
  Z <- matrix(c(-0.75, 0.1, 0.2, -0.4, 0.8, 0.9), ncol = 2L, byrow = TRUE)
  correlations <- list(
    matrix(c(1, 0.2, -0.1, 0.2, 1, 0.25, -0.1, 0.25, 1), 3L, byrow = TRUE),
    matrix(c(1, -0.15, 0.3, -0.15, 1, 0.1, 0.3, 0.1, 1), 3L, byrow = TRUE),
    matrix(c(1, 0.35, 0.05, 0.35, 1, -0.2, 0.05, -0.2, 1), 3L, byrow = TRUE),
    matrix(c(1, 0.05, -0.25, 0.05, 1, 0.4, -0.25, 0.4, 1), 3L, byrow = TRUE),
    matrix(c(1, -0.3, -0.05, -0.3, 1, 0.2, -0.05, 0.2, 1), 3L, byrow = TRUE)
  )
  expected <- list(
    matrix(
      c(
        1, 0.05708402093331598, 0.09255648318579689,
        0.05708402093331598, 1, 0.18933701456210145,
        0.09255648318579689, 0.18933701456210145, 1
      ),
      3L,
      byrow = TRUE
    ),
    matrix(
      c(
        1, 0.1522435202047135, -0.06954343943812133,
        0.1522435202047135, 1, 0.07939956052826905,
        -0.06954343943812133, 0.07939956052826905, 1
      ),
      3L,
      byrow = TRUE
    ),
    matrix(
      c(
        1, -0.29276893400126913, -0.026238636882333348,
        -0.29276893400126913, 1, 0.30022827876329217,
        -0.026238636882333348, 0.30022827876329217, 1
      ),
      3L,
      byrow = TRUE
    )
  )

  # The reference was generated from NPTFrechet/functions/corr_barycenter.py
  # using these exact precomputed correlations, X, and Z. Bypassing raw samples
  # keeps quantile conventions outside this regression-kernel comparison.
  core <- getFromNamespace("correlation_regression_cpp", "NonparaTransport")
  actual <- core(correlations, X, Z, 1000L, 1e-6, 1e-10)$correlations
  expect_equal(actual, expected, tolerance = 1e-12)
})
