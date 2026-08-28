#' Pairwise Squared Nonparanormal Transport Distances
#'
#' Computes the matrix of pairwise squared nonparanormal transport (NPT) distances
#' across distributions using precomputed representations from [as_nonparanormal()].
#'
#' @details
#' Under the nonparanormal model, the squared distance between two continuous 
#' multivariate distributions \eqn{P} and \eqn{Q} decomposes additively into
#' marginal quantile differences and latent Gaussian correlation distances:
#' \deqn{\mathrm{NPT}_c^2(P, Q) = \sum_{j=1}^d \|Q_{P,j} - Q_{Q,j}\|_{L^2}^2 + c\,\mathrm{BW}^2(R_P, R_Q)}
#'
#' \describe{
#'   \item{\strong{Marginal distance (\eqn{L^2} between quantiles)}}{
#'     The sum of squared \eqn{L^2} distances between marginal quantile functions:
#'     \deqn{\sum_{j=1}^d \|Q_{P,j} - Q_{Q,j}\|_{L^2}^2 = \sum_{j=1}^d \int_0^1 \left( Q_{P,j}(p) - Q_{Q,j}(p) \right)^2 dp}
#'     approximated via midpoint numerical integration on the probability grid.
#'   }
#'   \item{\strong{Weighted correlation distance (Bures--Wasserstein)}}{
#'     The squared Bures--Wasserstein distance between the latent Gaussian correlation matrices
#'     is multiplied by the positive weight \eqn{c=\code{bw_weight}}:
#'     \deqn{c\,\mathrm{BW}^2(R_P, R_Q) = c\left[\operatorname{tr}(R_P) + \operatorname{tr}(R_Q) - 2 \operatorname{tr}\left( \left( R_P^{1/2} R_Q R_P^{1/2} \right)^{1/2} \right)\right]}
#'     For bivariate distributions (\eqn{d = 2}), an exact closed-form scalar formula is used for maximum speed.
#'   }
#' }
#'
#' @param summaries A `nonparanormal` object returned by [as_nonparanormal()].
#' @param decompose Logical; if `TRUE`, returns a named list with the separate
#'   `marginal`, `correlation`, and `total` distance matrices. If `FALSE` (default),
#'   returns only the `total` squared-distance matrix.
#' @param bw_weight Positive finite scalar \eqn{c} multiplying the squared
#'   Bures--Wasserstein correlation term (default is `1`). Changing this weight
#'   does not alter the marginal term.
#'
#' @return A symmetric \eqn{n \times n} numeric matrix of total squared pairwise
#'   distances, where \eqn{n} is the number of distributions, or a list containing:
#'   \describe{
#'     \item{\code{marginal}}{\eqn{n \times n} matrix of squared marginal quantile distances.}
#'     \item{\code{correlation}}{\eqn{n \times n} matrix of weighted squared
#'       Bures--Wasserstein correlation distances \eqn{c\,\mathrm{BW}^2}.}
#'     \item{\code{total}}{\eqn{n \times n} matrix of total squared NPT distances
#'       (\code{marginal + correlation}).}
#'   }
#'
#' @seealso [as_nonparanormal()], [npt_distance()]
#' @export
pairwise_npt_distance <- function(summaries, decompose = FALSE, bw_weight = 1) {
  # 1. Validate inputs
  summaries <- .require_npt_summaries(summaries)
  decompose <- .as_flag(decompose, "decompose")
  bw_weight <- .as_positive_scalar(bw_weight, "bw_weight")

  # 2. Compute and label the complete pairwise matrix
  distances <- .summary_distances(summaries, bw_weight = bw_weight)
  for (name in names(distances)) {
    dimnames(distances[[name]]) <- list(
      summaries$distribution_names,
      summaries$distribution_names
    )
  }

  # 3. Return full decomposition or total matrix
  if (decompose) distances else distances$total
}

#' Squared Nonparanormal Transport Distance Between Two Distributions
#'
#' Computes one squared NPT distance between two distributions selected from a
#' common `nonparanormal` representation object. Using the common representation object is important
#' when [as_nonparanormal()] applied pooled median/MAD preprocessing:
#' both the single distance and the corresponding pairwise-matrix entry then use
#' exactly the same global transformation.
#'
#' @inheritParams pairwise_npt_distance
#' @param first,second A distribution name or a one-based integer index selecting
#'   the two distributions in `summaries`.
#'
#' @return A numeric scalar containing the total squared NPT distance. If
#'   `decompose = TRUE`, returns a list with scalar components `marginal`,
#'   `correlation` (the weighted term \eqn{c\,\mathrm{BW}^2}), and `total`.
#'
#' @seealso [pairwise_npt_distance()], [as_nonparanormal()]
#' @export
npt_distance <- function(
    summaries,
    first,
    second,
    decompose = FALSE,
    bw_weight = 1) {
  summaries <- .require_npt_summaries(summaries)
  decompose <- .as_flag(decompose, "decompose")
  bw_weight <- .as_positive_scalar(bw_weight, "bw_weight")

  first_index <- .distribution_index(
    first,
    summaries$distribution_names,
    "first"
  )
  second_index <- .distribution_index(
    second,
    summaries$distribution_names,
    "second"
  )

  distances <- .summary_distances(
    summaries,
    indices = c(first_index, second_index),
    bw_weight = bw_weight
  )
  components <- lapply(distances, function(matrix) unname(matrix[1L, 2L]))

  if (decompose) components else components$total
}

# Compute distance components for the full summary object or a two-row subset.
.summary_distances <- function(summaries, indices = NULL, bw_weight) {
  quantiles <- summaries$quantiles
  correlations <- summaries$correlations

  # Reuse cached R^{1/2} matrices when available. For a single requested pair,
  # subset the cache together with the quantiles and correlations.
  square_roots <- summaries$correlation_sqrts

  if (!is.null(indices)) {
    quantiles <- lapply(
      quantiles,
      function(values) values[indices, , drop = FALSE]
    )
    correlations <- correlations[indices]
    if (!is.null(square_roots)) {
      square_roots <- square_roots[indices]
    }
  }

  if (is.null(square_roots)) {
    square_roots <- list()
  }

  distances <- pairwise_summary_distance_cpp(
    quantiles,
    correlations,
    square_roots
  )
  distances$correlation <- bw_weight * distances$correlation
  distances$total <- distances$marginal + distances$correlation
  distances
}

# Resolve a name or one-based integer without partial matching.
.distribution_index <- function(value, distribution_names, argument) {
  if (is.character(value) && length(value) == 1L && !is.na(value)) {
    index <- match(value, distribution_names)
    if (is.na(index)) {
      stop(
        sprintf("`%s` does not match a distribution name.", argument),
        call. = FALSE
      )
    }
    return(index)
  }

  if (is.numeric(value) && length(value) == 1L && is.finite(value) &&
      value == as.integer(value) && value >= 1L &&
      value <= length(distribution_names)) {
    return(as.integer(value))
  }

  stop(
    sprintf(
      "`%s` must be one distribution name or an integer index between 1 and %d.",
      argument,
      length(distribution_names)
    ),
    call. = FALSE
  )
}
