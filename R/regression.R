#' Nonparanormal Fréchet Regression for Multivariate Distributional Responses
#'
#' Fits global Fréchet regression where the response is a collection of empirically
#' sampled multivariate distributions endowed with nonparanormal transport (NPT).
#'
#' @details
#' Under the nonparanormal model, the Fréchet regression of multivariate
#' distributional responses decomposes into independent marginal regressions and a
#' latent Gaussian correlation regression:
#' \enumerate{
#'   \item \strong{Marginal quantile regressions}: For each variable \eqn{j = 1, \dots, d},
#'     the 1D 2-Wasserstein Fréchet regression of marginal quantile functions
#'     \eqn{Q_{P_i, j}} against predictors \eqn{X} is computed via
#'     \code{fastfrechet::frechetreg_univar2wass}.
#'   \item \strong{Latent correlation regression}: The latent Gaussian correlation
#'     matrices \eqn{R_{P_i}} are regressed using projected Bures--Wasserstein (BW)
#'     Riemannian gradient descent (Park and Gaynanova 2026). For bivariate responses
#'     (\eqn{d = 2}), an exact closed-form solution is used.
#' }
#'
#' Passing precomputed response representations \code{Y} (from \code{\link{as_nonparanormal}})
#' avoids repeated sorting and correlation estimation during model fitting.
#'
#' @param X A numeric vector, matrix, or data frame of training predictors (\eqn{n \times p}).
#'   Rows correspond to the \eqn{n} distributions in \code{Y}.
#' @param Y A \code{nonparanormal} object returned by \code{\link{as_nonparanormal}}.
#' @param Z Optional numeric vector, matrix, or data frame of predictor values for prediction.
#'   If \code{NULL} (default), fitted values at the training points \code{X} are computed.
#' @param bounds Optional marginal support bounds \eqn{[a_j, b_j]} for each response variable.
#'   Can be a single length-2 numeric vector applied to all variables, a named list of length-2 vectors,
#'   or a \eqn{d \times 2} matrix with column names \code{"lower"} and \code{"upper"}.
#'   Default is \code{c(-Inf, Inf)} for all variables.
#' @param max_iter Maximum number of projected Riemannian iterations for latent
#'   correlation estimation when \eqn{d > 2} (default is `1000L`).
#' @param tol Convergence tolerance on the Frobenius norm between successive
#'   correlation iterates when \eqn{d > 2} (default is `1e-6`).
#'
#' @return An object of class \code{"npt_frechetreg"}, containing:
#' \describe{
#'   \item{\code{quantiles}}{A named list of \eqn{d} matrices (one per variable), each of dimension
#'     \eqn{k \times M}, containing predicted marginal quantiles across \eqn{k} prediction points.}
#'   \item{\code{correlations}}{A named list of \eqn{k} predicted \eqn{d \times d} latent correlation matrices.}
#'   \item{\code{probabilities}}{Numeric vector of length \eqn{M} holding the probability evaluation grid.}
#'   \item{\code{Z}}{The \eqn{k \times p} matrix of prediction predictor values.}
#'   \item{\code{prediction_names}}{Character vector of prediction point labels.}
#'   \item{\code{variable_names}}{Character vector of response variable names.}
#'   \item{\code{in_sample}}{Logical indicator that the fit was evaluated at the
#'     training predictors by calling the function with \code{Z = NULL}.}
#'   \item{\code{observed_correlations}}{For an in-sample fit, the observed latent
#'     correlation matrices aligned with the fitted correlations; otherwise \code{NULL}.}
#'   \item{\code{max_iter}, \code{tol}}{Correlation-regression controls retained
#'     so in-sample visualization methods can evaluate fitted trajectories.}
#'   \item{\code{correlation_diagnostics}}{A data frame with convergence diagnostics for each correlation fit.}
#' }
#'
#' @seealso \code{\link{as_nonparanormal}}, \code{\link{pairwise_npt_distance}}
#' @export
npt_frechetreg <- function(
    X,
    Y,
    Z = NULL,
    bounds = NULL,
    max_iter = 1000L,
    tol = 1e-6) {
  # Validate inputs and prepare predictor/response structures
  prepared <- .prepare_npt_regression(X, Y, Z, bounds, max_iter, tol)
  .fit_npt_regression(prepared)
}

#' @export
print.npt_frechetreg <- function(x, ...) {
  cat("<npt_frechetreg>\n")
  cat("Prediction points:", nrow(x$Z), "\n")
  cat("Predictors:", ncol(x$Z), "\n")
  cat("Marginal components:", length(x$quantiles), "\n")
  if (length(x$quantiles) > 2L) {
    cat(
      "Correlation fits converged:",
      sum(x$correlation_diagnostics$converged),
      "of",
      nrow(x$correlation_diagnostics),
      "\n"
    )
  }
  invisible(x)
}

# Input validation helper for npt_frechetreg
.prepare_npt_regression <- function(X, Y, Z, bounds, max_iter, tol) {
  # 1. Validate response summaries and align training predictors
  Y <- .require_npt_summaries(Y)
  number_distributions <- length(Y$correlations)
  X_info <- .as_training_predictor_matrix(
    X,
    expected_rows = number_distributions,
    distribution_names = Y$distribution_names
  )

  # 2. Format prediction points (Z = NULL defaults to in-sample training points)
  at_training_points <- is.null(Z)
  if (at_training_points) {
    Z_matrix <- X_info$matrix
    prediction_names <- Y$distribution_names
  } else {
    Z_info <- .as_prediction_matrix(
      Z,
      predictor_names = X_info$predictor_names
    )
    Z_matrix <- Z_info$matrix
    prediction_names <- Z_info$row_names
  }

  # 3. Validate control parameters and normalize marginal bounds
  max_iter <- .as_integer_count(max_iter, "max_iter", minimum = 1L)
  tol <- .as_positive_scalar(tol, "tol")

  list(
    X = X_info$matrix,
    Y = Y,
    Z = Z_matrix,
    fastfrechet_Z = if (at_training_points) NULL else Z_matrix,
    in_sample = at_training_points,
    prediction_names = prediction_names,
    predictor_names = X_info$predictor_names,
    bounds = .normalize_marginal_bounds(bounds, Y$variable_names),
    max_iter = max_iter,
    tol = tol
  )
}

# Actual regression fit logic
.fit_npt_regression <- function(prepared) {
  Y <- prepared$Y
  variable_names <- Y$variable_names
  number_variables <- length(variable_names)

  # 1. Fit univariate Frechet regression for each marginal quantile function
  marginal_predictions <- vector("list", number_variables)
  for (variable in seq_len(number_variables)) {
    marginal_fit <- fastfrechet::frechetreg_univar2wass(
      X = prepared$X,
      Y = Y$quantiles[[variable]],
      Z = prepared$fastfrechet_Z,
      lower = prepared$bounds[variable, "lower"],
      upper = prepared$bounds[variable, "upper"]
    )
    prediction <- as.matrix(marginal_fit$Qhat)
    storage.mode(prediction) <- "double"
    rownames(prediction) <- prepared$prediction_names
    marginal_predictions[[variable]] <- prediction
  }
  names(marginal_predictions) <- variable_names

  # 2. Fit latent correlation regression across all prediction points (in C++)
  correlation_fit <- correlation_regression_cpp(
    Y$correlations,
    prepared$X,
    prepared$Z,
    prepared$max_iter,
    prepared$tol,
    1e-12     # eigenvalue numerical stability threshold
  )
  correlations <- correlation_fit$correlations
  names(correlations) <- prepared$prediction_names
  for (prediction in seq_along(correlations)) {
    dimnames(correlations[[prediction]]) <- list(variable_names, variable_names)
  }

  # 3. Assemble convergence diagnostics for correlation estimation
  diagnostics <- data.frame(
    prediction = prepared$prediction_names,
    iterations = as.integer(correlation_fit$iterations),
    converged = as.logical(correlation_fit$converged),
    final_change = as.numeric(correlation_fit$final_change),
    status = as.character(correlation_fit$status),
    row.names = prepared$prediction_names,
    check.names = FALSE
  )

  Z_output <- prepared$Z
  rownames(Z_output) <- prepared$prediction_names
  colnames(Z_output) <- prepared$predictor_names

  # 4. Package and return the fitted regression object
  structure(
    list(
      quantiles = marginal_predictions,
      correlations = correlations,
      probabilities = Y$probabilities,
      Z = Z_output,
      prediction_names = prepared$prediction_names,
      variable_names = variable_names,
      in_sample = prepared$in_sample,
      observed_correlations = if (prepared$in_sample) Y$correlations else NULL,
      max_iter = prepared$max_iter,
      tol = prepared$tol,
      correlation_diagnostics = diagnostics
    ),
    class = "npt_frechetreg"
  )
}

.as_training_predictor_matrix <- function(
    X,
    expected_rows,
    distribution_names) {
  explicit_rows <- .explicit_row_names(X)
  matrix <- .coerce_numeric_predictors(X, "X", vector_as_row = FALSE)

  if (nrow(matrix) != expected_rows) {
    stop(
      sprintf("`X` must have %d rows, one per distribution in `Y`.", expected_rows),
      call. = FALSE
    )
  }
  if (ncol(matrix) < 1L) {
    stop("`X` must contain at least one predictor column.", call. = FALSE)
  }

  # Align rows of X if explicit distribution row names are provided
  if (!is.null(explicit_rows)) {
    if (anyDuplicated(explicit_rows) || !setequal(explicit_rows, distribution_names)) {
      stop(
        "Explicit row names of `X` must match the distribution names in `Y`.",
        call. = FALSE
      )
    }
    matrix <- matrix[match(distribution_names, explicit_rows), , drop = FALSE]
  }

  # Ensure valid predictor column names
  predictor_names <- colnames(matrix)
  if (is.null(predictor_names)) {
    predictor_names <- paste0("X", seq_len(ncol(matrix)))
  }
  if (anyNA(predictor_names) || any(!nzchar(predictor_names)) || anyDuplicated(predictor_names)) {
    stop("Predictor column names must be non-empty and unique.", call. = FALSE)
  }
  colnames(matrix) <- predictor_names
  rownames(matrix) <- distribution_names

  list(matrix = matrix, predictor_names = predictor_names)
}

.as_prediction_matrix <- function(Z, predictor_names) {
  expected_columns <- length(predictor_names)
  vector_as_row <- expected_columns > 1L
  explicit_rows <- .explicit_row_names(Z)

  # Handle vector inputs: k points when p = 1; single point of length p when p > 1
  if (is.numeric(Z) && is.null(dim(Z))) {
    if (expected_columns == 1L) {
      matrix <- matrix(as.numeric(Z), ncol = 1L)
    } else if (length(Z) == expected_columns) {
      matrix <- matrix(as.numeric(Z), nrow = 1L)
    } else {
      stop(
        sprintf("A vector `Z` must have length %d for this model.", expected_columns),
        call. = FALSE
      )
    }
  } else {
    matrix <- .coerce_numeric_predictors(Z, "Z", vector_as_row = vector_as_row)
  }

  if (nrow(matrix) < 1L || ncol(matrix) != expected_columns) {
    stop(
      sprintf("`Z` must have %d predictor columns and at least one row.", expected_columns),
      call. = FALSE
    )
  }

  # Match column order to X if column names are supplied
  supplied_names <- colnames(matrix)
  if (!is.null(supplied_names)) {
    if (anyDuplicated(supplied_names) || !setequal(supplied_names, predictor_names)) {
      stop("The predictor columns of `Z` must match those of `X`.", call. = FALSE)
    }
    matrix <- matrix[, match(predictor_names, supplied_names), drop = FALSE]
  }
  colnames(matrix) <- predictor_names

  # Assign row names for prediction points
  row_names <- explicit_rows
  if (is.null(row_names)) {
    row_names <- paste0("prediction_", seq_len(nrow(matrix)))
  }
  if (length(row_names) != nrow(matrix) || anyNA(row_names) ||
      any(!nzchar(row_names)) || anyDuplicated(row_names)) {
    stop("Prediction-point row names must be non-empty and unique.", call. = FALSE)
  }
  rownames(matrix) <- row_names

  list(matrix = matrix, row_names = row_names)
}

.coerce_numeric_predictors <- function(x, name, vector_as_row) {
  if (is.numeric(x) && is.null(dim(x))) {
    matrix <- if (vector_as_row) matrix(x, nrow = 1L) else matrix(x, ncol = 1L)
  } else if (is.data.frame(x)) {
    numeric_columns <- vapply(x, is.numeric, logical(1))
    if (!all(numeric_columns)) {
      stop(sprintf("Every column of `%s` must be numeric.", name), call. = FALSE)
    }
    matrix <- as.matrix(x)
  } else if (is.matrix(x) && is.numeric(x)) {
    matrix <- x
  } else {
    stop(sprintf("`%s` must be a numeric vector, matrix, or data frame.", name), call. = FALSE)
  }

  if (any(!is.finite(matrix))) {
    stop(sprintf("`%s` must contain only finite values.", name), call. = FALSE)
  }
  storage.mode(matrix) <- "double"
  matrix
}

.explicit_row_names <- function(x) {
  # Ignore automatic integer row names from data frames
  if (is.data.frame(x) && .row_names_info(x, type = 1L) < 0L) {
    return(NULL)
  }
  rownames(x)
}

.normalize_marginal_bounds <- function(bounds, variable_names) {
  number_variables <- length(variable_names)

  # Standardize bounds to a d-by-2 numeric matrix with columns "lower" and "upper"
  if (is.null(bounds)) {
    result <- cbind(
      lower = rep(-Inf, number_variables),
      upper = rep(Inf, number_variables)
    )
  } else if (is.numeric(bounds) && is.null(dim(bounds))) {
    if (length(bounds) != 2L) {
      stop("A numeric `bounds` vector must have length two.", call. = FALSE)
    }
    result <- matrix(
      rep(as.numeric(bounds), each = number_variables),
      nrow = number_variables,
      dimnames = list(NULL, c("lower", "upper"))
    )
  } else if (is.list(bounds) && !is.data.frame(bounds)) {
    if (length(bounds) != number_variables) {
      stop("A `bounds` list must have one element per response variable.", call. = FALSE)
    }
    if (!is.null(names(bounds))) {
      if (anyDuplicated(names(bounds)) || !setequal(names(bounds), variable_names)) {
        stop("Named `bounds` must match the response variable names.", call. = FALSE)
      }
      bounds <- bounds[variable_names]
    }
    valid <- vapply(
      bounds,
      function(value) is.numeric(value) && length(value) == 2L,
      logical(1)
    )
    if (!all(valid)) {
      stop("Every element of `bounds` must be a length-two numeric vector.", call. = FALSE)
    }
    result <- do.call(rbind, lapply(bounds, as.numeric))
    colnames(result) <- c("lower", "upper")
  } else if ((is.matrix(bounds) || is.data.frame(bounds)) &&
             all(dim(bounds) == c(number_variables, 2L))) {
    if (is.data.frame(bounds) && !all(vapply(bounds, is.numeric, logical(1)))) {
      stop("Both columns of `bounds` must be numeric.", call. = FALSE)
    }
    result <- as.matrix(bounds)
    if (!is.numeric(result)) {
      stop("`bounds` must be numeric.", call. = FALSE)
    }
    if (!is.null(rownames(result))) {
      if (anyDuplicated(rownames(result)) || !setequal(rownames(result), variable_names)) {
        stop("The row names of `bounds` must match the response variables.", call. = FALSE)
      }
      result <- result[match(variable_names, rownames(result)), , drop = FALSE]
    }
    colnames(result) <- c("lower", "upper")
  } else {
    stop(
      "`bounds` must be NULL, a length-two numeric vector, a named list, or a d-by-2 numeric table.",
      call. = FALSE
    )
  }

  storage.mode(result) <- "double"
  rownames(result) <- variable_names
  if (anyNA(result) || any(result[, "lower"] >= result[, "upper"])) {
    stop("Each lower bound must be strictly smaller than its upper bound.", call. = FALSE)
  }
  result
}
