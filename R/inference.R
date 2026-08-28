#' Component-wise R-squared for Nonparanormal Fréchet Regression
#'
#' Computes component-wise generalized \eqn{R^2} goodness-of-fit statistics for
#' nonparanormal Fréchet regression: one \eqn{R^2} value for each marginal quantile
#' function and one for the latent Gaussian correlation structure.
#'
#' @details
#' Under the nonparanormal transport framework, total and residual variation
#' decompose additively into marginal quantile components and a latent correlation
#' component. For each component, the generalized \eqn{R^2} is defined as:
#' \deqn{R^2 = 1 - \frac{\sum_{i=1}^n d^2(P_i, \hat{P}(X_i))}{\sum_{i=1}^n d^2(P_i, \hat{P}(\bar{X}))}}
#' where:
#' \itemize{
#'   \item \strong{Marginal components (\eqn{j = 1, \dots, d})}: Evaluates squared
#'     \eqn{L^2} Wasserstein loss between observed and predicted marginal quantile functions:
#'     \deqn{R_j^2 = 1 - \frac{\sum_{i=1}^n \|Q_i^{(j)} - \hat{Q}^{(j)}(X_i)\|_{L^2}^2}{\sum_{i=1}^n \|Q_i^{(j)} - \hat{Q}^{(j)}(\bar{X})\|_{L^2}^2}}
#'   \item \strong{Latent correlation component}: Evaluates squared Bures--Wasserstein
#'     loss between observed and predicted latent correlation matrices:
#'     \deqn{R_{\mathrm{corr}}^2 = 1 - \frac{\sum_{i=1}^n \mathcal{B}^2(\Sigma_i, \hat{\Sigma}(X_i))}{\sum_{i=1}^n \mathcal{B}^2(\Sigma_i, \hat{\Sigma}(\bar{X}))}}
#' }
#' Here \eqn{\bar{X} = \frac{1}{n} \sum_{i=1}^n X_i} denotes the sample mean predictor vector.
#' If a component exhibits zero total variation around \eqn{\bar{X}}, its \eqn{R^2} is
#' undefined and returned as \code{NA}.
#'
#' @inheritParams npt_frechetreg
#'
#' @return A data frame with one row per component and columns:
#' \describe{
#'   \item{\code{component}}{Name of the response component (marginal variable names and \code{"latent_correlation"}).}
#'   \item{\code{r_squared}}{Estimated component-wise \eqn{R^2} value in \eqn{(-\infty, 1]}.}
#'   \item{\code{residual_sum}}{Sum of squared residual Fréchet losses around fitted values \eqn{\hat{P}(X_i)}.}
#'   \item{\code{total_sum}}{Sum of squared total Fréchet losses around the baseline fit \eqn{\hat{P}(\bar{X})}.}
#' }
#'
#' @seealso \code{\link{npt_frechetreg}}, \code{\link{npt_permutation_test}}
#' @export
npt_component_r2 <- function(
    X,
    Y,
    bounds = NULL,
    max_iter = 1000L,
    tol = 1e-6) {
  # Fit in-sample and baseline regressions to compute component-wise R^2
  prepared <- .prepare_npt_regression(X, Y, NULL, bounds, max_iter, tol)
  calculation <- .component_r2_from_prepared(prepared)
  calculation$table
}

#' Permutation Inference for Component-wise R-squared
#'
#' Tests the joint significance of predictors \eqn{X} on each response component
#' (marginal quantile functions and latent Gaussian correlation matrix) via
#' permutation testing.
#'
#' @details
#' Permutation tests evaluate the null hypothesis of independence between the
#' predictors \eqn{X} and the distributional responses \eqn{Y}.
#' In each permutation replicate \eqn{b = 1, \dots, B}:
#' \enumerate{
#'   \item Rows of \eqn{X} are randomly permuted together across distributions,
#'     preserving the correlation structure among predictors.
#'   \item The Fréchet regression model is refit to the permuted data.
#'   \item Permuted component-wise \eqn{R^2} values are computed using the fixed
#'     baseline total sum of squares.
#' }
#'
#' Raw permutation p-values use the standard add-one correction:
#' \deqn{p_\ell = \frac{1 + \sum_{b=1}^B \mathbf{1}\left( R_{\ell, (b)}^2 \ge R_{\ell, \mathrm{obs}}^2 \right)}{B + 1},}
#' where \eqn{\ell} indexes the \eqn{d} marginal components and the latent correlation component.
#' Family-wise error rate (FWER) across all \eqn{d + 1} components is controlled
#' using the single-step Westfall--Young min-\eqn{p} procedure.
#'
#' @inheritParams npt_component_r2
#' @param B Number of random permutations (default is `999L`).
#' @param workers Number of parallel CPU worker processes (default is `1L` for
#'   sequential execution).
#' @param seed Optional random seed integer for reproducible permutation tests.
#'
#' @return An object of class \code{"npt_permutation_test"}, containing:
#' \describe{
#'   \item{\code{results}}{A data frame summarizing observed \eqn{R^2}, unadjusted p-values,
#'     and multiplicity-adjusted p-values per component.}
#'   \item{\code{observed}}{The observed component-wise \eqn{R^2} summary table.}
#'   \item{\code{null_distribution}}{A \eqn{B \times (d + 1)} matrix of permuted \eqn{R^2} statistics.}
#'   \item{\code{B}}{Number of permutations performed.}
#'   \item{\code{workers}}{Number of parallel workers used.}
#'   \item{\code{seed}}{Random seed used (if specified).}
#'   \item{\code{adjustment}}{Multiple testing correction method name.}
#' }
#'
#' @seealso \code{\link{npt_component_r2}}, \code{\link{npt_frechetreg}}
#' @export
npt_permutation_test <- function(
    X,
    Y,
    B = 999L,
    workers = 1L,
    seed = NULL,
    bounds = NULL,
    max_iter = 1000L,
    tol = 1e-6) {
  B <- .as_integer_count(B, "B", minimum = 1L)
  workers <- .as_integer_count(workers, "workers", minimum = 1L)
  workers <- min(workers, B)
  seed <- .validate_optional_seed(seed)

  # 1. Compute observed R^2 statistics and baseline denominator
  prepared <- .prepare_npt_regression(X, Y, NULL, bounds, max_iter, tol)
  observed <- .component_r2_from_prepared(prepared)
  observed_values <- observed$table$r_squared
  names(observed_values) <- observed$table$component

  # 2. Pre-generate permutation indices for reproducible parallel execution
  if (!is.null(seed)) {
    set.seed(seed)
  }
  permutations <- replicate(
    B,
    sample.int(nrow(prepared$X), replace = FALSE),
    simplify = FALSE
  )

  # 3. Evaluate permuted R^2 across replicates (sequential or parallel)
  if (workers == 1L) {
    null_rows <- lapply(
      permutations,
      .permutation_component_r2,
      prepared = prepared,
      denominator = observed$denominator
    )
  } else {
    cluster <- parallel::makeCluster(workers)
    on.exit(parallel::stopCluster(cluster), add = TRUE)

    library_paths <- .libPaths()
    parallel::clusterCall(
      cluster,
      function(paths) {
        .libPaths(paths)
        library(NonparaTransport)
        NULL
      },
      library_paths
    )

    worker_context <- list(
      prepared = prepared,
      denominator = observed$denominator
    )
    worker_function <- .permutation_component_r2
    parallel::clusterExport(
      cluster,
      c("worker_context", "worker_function"),
      envir = environment()
    )
    null_rows <- parallel::parLapply(
      cluster,
      permutations,
      function(index) {
        worker_function(
          index,
          worker_context$prepared,
          worker_context$denominator
        )
      }
    )
  }

  # 4. Assemble null distribution matrix and compute Westfall-Young p-values
  null_distribution <- do.call(rbind, null_rows)
  colnames(null_distribution) <- names(observed_values)
  p_values <- .westfall_young_min_p(null_distribution, observed_values)

  result_table <- data.frame(
    component = names(observed_values),
    r_squared = as.numeric(observed_values),
    p_value = unname(p_values$unadjusted),
    p_value_adjusted = unname(p_values$adjusted),
    row.names = names(observed_values),
    check.names = FALSE
  )

  structure(
    list(
      results = result_table,
      observed = observed$table,
      null_distribution = null_distribution,
      B = B,
      workers = workers,
      seed = seed,
      adjustment = "Westfall-Young single-step min-p"
    ),
    class = "npt_permutation_test"
  )
}

#' @export
print.npt_permutation_test <- function(x, ...) {
  cat("<npt_permutation_test>\n")
  cat("Permutations:", x$B, "\n")
  cat("Workers:", x$workers, "\n")
  print(x$results, row.names = FALSE)
  invisible(x)
}

# ------------------------------------------------------------------------------
# Internal Calculation Helpers
# ------------------------------------------------------------------------------

# Computes observed component-wise R^2 table and total variation denominator
.component_r2_from_prepared <- function(prepared) {
  # 1. Residual variation: losses around in-sample fitted values hat{P}(X_i)
  fitted <- .fit_npt_regression(prepared)
  numerator <- .component_losses(prepared$Y, fitted)

  # 2. Total variation: losses around baseline prediction at sample mean predictor mean(X)
  baseline_prepared <- prepared
  baseline_prepared$Z <- matrix(
    colMeans(prepared$X),
    nrow = 1L,
    dimnames = list("mean_predictor", prepared$predictor_names)
  )
  baseline_prepared$fastfrechet_Z <- baseline_prepared$Z
  baseline_prepared$prediction_names <- "mean_predictor"
  baseline <- .fit_npt_regression(baseline_prepared)
  denominator <- .component_losses(prepared$Y, baseline)

  # 3. Generalized R^2 = 1 - (residual_sum / total_sum)
  r_squared <- rep(NA_real_, length(denominator))
  positive_variation <- .has_component_variation(denominator)
  r_squared[positive_variation] <-
    1 - numerator[positive_variation] / denominator[positive_variation]

  component_names <- names(denominator)
  table <- data.frame(
    component = component_names,
    r_squared = unname(r_squared),
    residual_sum = unname(numerator),
    total_sum = unname(denominator),
    row.names = component_names,
    check.names = FALSE
  )

  list(table = table, denominator = denominator)
}

# Computes sum of squared Fréchet losses across marginal and correlation components
.component_losses <- function(Y, predictions) {
  losses <- matched_component_losses_cpp(
    Y$quantiles,
    predictions$quantiles,
    Y$correlations,
    predictions$correlations,
    1e-15
  )
  result <- stats::setNames(
    as.numeric(losses$marginal),
    Y$variable_names
  )
  if (length(Y$variable_names) > 1L) {
    result <- c(result, latent_correlation = as.numeric(losses$correlation))
  }
  result
}

# Computes permuted component-wise R^2 for a single permutation replicate
.permutation_component_r2 <- function(index, prepared, denominator) {
  permuted <- prepared
  # Permute full rows of X jointly to preserve predictor covariance
  permuted$X <- prepared$X[index, , drop = FALSE]
  rownames(permuted$X) <- prepared$Y$distribution_names
  permuted$Z <- permuted$X
  permuted$fastfrechet_Z <- NULL
  permuted$prediction_names <- prepared$Y$distribution_names

  # Fit regression on permuted predictors and evaluate residual losses
  fitted <- .fit_npt_regression(permuted)
  numerator <- .component_losses(prepared$Y, fitted)

  # Baseline denominator is unchanged under predictor permutation (mean(X) is invariant)
  result <- rep(NA_real_, length(denominator))
  positive_variation <- .has_component_variation(denominator)
  result[positive_variation] <-
    1 - numerator[positive_variation] / denominator[positive_variation]
  names(result) <- names(denominator)
  result
}

# Computes raw and Westfall-Young min-p adjusted permutation p-values
.westfall_young_min_p <- function(null_distribution, observed) {
  B <- nrow(null_distribution)
  component_names <- colnames(null_distribution)
  unadjusted <- adjusted <- stats::setNames(
    rep(NA_real_, length(component_names)),
    component_names
  )

  valid <- is.finite(observed) & vapply(
    seq_along(component_names),
    function(column) all(is.finite(null_distribution[, column])),
    logical(1)
  )
  if (!any(valid)) {
    return(list(unadjusted = unadjusted, adjusted = adjusted))
  }

  null_valid <- null_distribution[, valid, drop = FALSE]
  observed_valid <- observed[valid]

  # 1. Rank observed and permuted statistics in a single (B + 1) reference set
  reference <- rbind(observed_valid, null_valid)
  reference_size <- B + 1L
  reference_p <- apply(
    reference,
    2L,
    function(values) rank(-values, ties.method = "max") / reference_size
  )
  if (is.null(dim(reference_p))) {
    reference_p <- matrix(reference_p, ncol = 1L)
  }

  # 2. Extract observed unadjusted p-values (first row) and component-wise minima
  raw <- reference_p[1L, ]
  minimum_p <- apply(reference_p, 1L, min)

  # 3. Compute Westfall--Young single-step min-p adjusted p-values
  adjusted_valid <- vapply(
    raw,
    function(value) mean(minimum_p <= value),
    numeric(1)
  )

  unadjusted[valid] <- raw
  adjusted[valid] <- adjusted_valid
  list(unadjusted = unadjusted, adjusted = adjusted)
}

# Validates an optional integer seed
.validate_optional_seed <- function(seed) {
  if (is.null(seed)) {
    return(NULL)
  }
  integer_minimum <- -.Machine$integer.max - 1
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != trunc(seed) || seed < integer_minimum ||
      seed > .Machine$integer.max) {
    stop("`seed` must be NULL or one finite integer.", call. = FALSE)
  }
  as.integer(seed)
}

# Checks if a total variation sum is strictly positive above roundoff threshold
.has_component_variation <- function(denominator) {
  denominator > 100 * .Machine$double.eps
}
