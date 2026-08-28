# ==============================================================================
# Input Validation and Preprocessing
#
# Functions in this file validate user inputs in R before passing them to the
# C++ computation engine. Handling checks in R allows for informative error
# messages and guarantees key invariants (sample sizes, shared variable names,
# finiteness, and statistical identifiability).
# ==============================================================================

#' Validate and normalize a list of empirical distribution matrices
#'
#' @param data A list of numeric matrices representing empirical distributions.
#' @return A list containing sanitized data matrices and metadata.
#' @noRd
.validate_distributional_data <- function(data) {
  # 1. Outer structure: data must be a non-empty list
  if (!is.list(data) || !length(data)) {
    stop("`data` must be a non-empty list of numeric matrices.", call. = FALSE)
  }

  # 2. Distribution names: used to label rows and columns in pairwise outputs.
  # If names are missing or blank, assign default 'distribution_1', 'distribution_2', etc.
  distribution_names <- names(data)
  if (is.null(distribution_names)) {
    distribution_names <- paste0("distribution_", seq_along(data))
  } else {
    blank <- is.na(distribution_names) | !nzchar(distribution_names)
    distribution_names[blank] <- paste0("distribution_", which(blank))
  }

  if (anyDuplicated(distribution_names)) {
    stop("Distribution names must be unique.", call. = FALSE)
  }
  names(data) <- distribution_names

  # 3. Matrix format and numeric type check
  valid_matrix <- vapply(
    data,
    function(x) is.matrix(x) && is.numeric(x),
    logical(1)
  )
  if (!all(valid_matrix)) {
    stop("Every element of `data` must be a numeric matrix.", call. = FALSE)
  }

  # 4. Dimension and sample size requirements
  # - Sample sizes N_i (rows) can vary, but each must have at least 2 observations
  #   to compute pairwise concordance (Kendall's tau).
  # - The number of variables d (columns) must be at least 2 and identical across distributions.
  dimensions <- lapply(data, dim)
  sample_sizes <- vapply(dimensions, `[[`, integer(1), 1L)
  variable_counts <- vapply(dimensions, `[[`, integer(1), 2L)

  if (any(sample_sizes < 2L)) {
    bad <- names(data)[which(sample_sizes < 2L)[1L]]
    stop(sprintf("Distribution `%s` must contain at least two rows.", bad), call. = FALSE)
  }
  if (any(variable_counts < 2L) || length(unique(variable_counts)) != 1L) {
    stop("All distribution matrices must have the same number of columns (at least 2).", call. = FALSE)
  }
  if (any(vapply(data, function(x) any(!is.finite(x)), logical(1)))) {
    stop("Distribution matrices must contain only finite values.", call. = FALSE)
  }

  # 5. Variable (column) consistency
  # All distributions must share the same variable identities in the exact same column order.
  # If column names are provided, they must match across all distributions;
  # otherwise, default to 'V1', 'V2', ..., 'Vd'.
  supplied_names <- lapply(data, colnames)
  any_names <- any(vapply(supplied_names, Negate(is.null), logical(1)))
  if (any_names) {
    if (any(vapply(supplied_names, is.null, logical(1)))) {
      stop("Column names must be supplied for every distribution or for none of them.", call. = FALSE)
    }
    reference_names <- supplied_names[[1L]]
    same_names <- vapply(supplied_names, identical, logical(1), reference_names)
    if (!all(same_names) || anyDuplicated(reference_names) || any(!nzchar(reference_names))) {
      stop("Column names must be non-empty, unique, and identical across distributions.", call. = FALSE)
    }
  } else {
    reference_names <- paste0("V", seq_len(variable_counts[[1L]]))
  }

  # 6. Statistical identifiability check
  # A constant variable within a sample has zero variance, making Kendall's tau
  # and latent copula correlations undefined.
  for (i in seq_along(data)) {
    constant <- vapply(
      seq_len(ncol(data[[i]])),
      function(j) all(data[[i]][, j] == data[[i]][1L, j]),
      logical(1)
    )
    if (any(constant)) {
      bad_variable <- reference_names[which(constant)[1L]]
      stop(
        sprintf(
          "Variable `%s` is constant in distribution `%s`; its latent correlation is not identifiable.",
          bad_variable,
          names(data)[i]
        ),
        call. = FALSE
      )
    }
  }

  # 7. Format data for C++ consumption
  # Ensure double precision and strip dimension names to avoid unnecessary copies.
  # Semantic names are re-attached to the summary object after C++ execution.
  data <- lapply(data, function(x) {
    storage.mode(x) <- "double"
    unname(x)
  })
  names(data) <- distribution_names

  list(
    data = data,
    distribution_names = distribution_names,
    variable_names = reference_names,
    sample_sizes = sample_sizes
  )
}

#' Standardize all distributions with pooled robust location and scale
#'
#' Distance preprocessing pools observations across distributions within each
#' variable. Its MAD is the literal sample median absolute deviation and does
#' not use R's normal-consistency multiplier.
#' @noRd
.standardize_pooled_median_mad <- function(data, variable_names) {
  number_variables <- ncol(data[[1L]])
  centers <- numeric(number_variables)
  scales <- numeric(number_variables)

  for (variable in seq_len(number_variables)) {
    pooled <- unlist(
      lapply(data, function(distribution) distribution[, variable]),
      use.names = FALSE
    )
    centers[variable] <- stats::median(pooled)
    scales[variable] <- stats::median(abs(pooled - centers[variable]))
  }
  names(centers) <- variable_names
  names(scales) <- variable_names

  zero_scale <- which(!is.finite(scales) | scales <= 0)
  if (length(zero_scale)) {
    stop(
      sprintf(
        paste0(
          "Pooled median/MAD standardization is undefined for variable `%s` ",
          "because its pooled MAD is zero."
        ),
        variable_names[zero_scale[1L]]
      ),
      call. = FALSE
    )
  }

  standardized <- lapply(data, function(distribution) {
    centered <- sweep(distribution, 2L, centers, "-")
    sweep(centered, 2L, scales, "/")
  })
  names(standardized) <- names(data)

  list(data = standardized, center = centers, scale = scales)
}

# ------------------------------------------------------------------------------
# Scalar and Control Option Validators
# ------------------------------------------------------------------------------

#' Validate an integer count parameter (e.g., quantile-grid size M)
#' @noRd
.as_integer_count <- function(x, name, minimum) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x != as.integer(x) || x < minimum) {
    stop(sprintf("`%s` must be a single integer >= %d.", name, minimum), call. = FALSE)
  }
  as.integer(x)
}

#' Validate a logical flag (e.g., cache_sqrt, decompose)
#' @noRd
.as_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  x
}

#' Validate a strictly positive finite scalar
#' @noRd
.as_positive_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop(sprintf("`%s` must be one positive finite number.", name), call. = FALSE)
  }
  as.numeric(x)
}

#' Validate the positive-definiteness shrinkage parameter lambda in (0, 1)
#' @noRd
.validate_pd_shrinkage <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x >= 1) {
    stop("`pd_shrinkage` must be a finite number strictly between 0 and 1.", call. = FALSE)
  }
  as.numeric(x)
}

#' Verify that an object is a precomputed 'nonparanormal' structure
#' @noRd
.require_npt_summaries <- function(x) {
  if (!inherits(x, "nonparanormal")) {
    stop("`summaries` must be created by `as_nonparanormal()`.", call. = FALSE)
  }

  # Quantile columns have a fixed numerical meaning: column k represents the
  # midpoint of probability cell k. Checking that contract here prevents an
  # endpoint or arbitrary grid from being combined with midpoint integration.
  probabilities <- x$probabilities
  if (!is.numeric(probabilities) || length(probabilities) < 2L ||
      any(!is.finite(probabilities))) {
    stop("`summaries$probabilities` must be a finite numeric grid.", call. = FALSE)
  }
  M <- length(probabilities)
  expected <- (seq_len(M) - 0.5) / M
  if (!isTRUE(all.equal(
    as.numeric(probabilities),
    expected,
    tolerance = 100 * .Machine$double.eps,
    check.attributes = FALSE
  ))) {
    stop(
      "`summaries$probabilities` must use the equal-interval midpoint grid.",
      call. = FALSE
    )
  }
  if (!is.list(x$quantiles) || !length(x$quantiles) ||
      any(vapply(x$quantiles, function(value) {
        !is.matrix(value) || !is.numeric(value) || ncol(value) != M
      }, logical(1)))) {
    stop("Every summary quantile matrix must have one column per midpoint.", call. = FALSE)
  }
  x
}
