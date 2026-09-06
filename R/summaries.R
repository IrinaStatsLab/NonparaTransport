#' Fit and Represent Empirical Distributions as Nonparanormal Distributions
#'
#' Estimates marginal empirical quantiles and latent Gaussian correlation
#' matrices from a collection of multivariate empirical distributions.
#' These nonparanormal representations serve as the core input for computing
#' nonparanormal transport metrics and fitting Fréchet regressions.
#'
#' @details
#' Under the nonparanormal model, the nonparanormal transport (NPT) metric between
#' two continuous multivariate distributions requires:
#' \enumerate{
#'   \item \strong{Marginal quantile functions}: Empirical quantiles evaluated
#'     at the midpoints \eqn{p_k=(k-1/2)/M}, \eqn{k=1,\ldots,M}, of \eqn{M}
#'     equal-width probability intervals.
#'   \item \strong{Latent correlation matrices}: Dependence structures estimated
#'     via Kendall's \eqn{\tau_a} and transformed through the sine bridge
#'     \eqn{\rho = \sin(\pi \tau_a / 2)}. The computation follows the default
#'     computation of the R package `latentcor`.
#' }
#'
#' Precomputing these representations once for each of the \eqn{n} distributions avoids
#' repeatedly sorting data and estimating correlations across all \eqn{O(n^2)}
#' pairs in \code{\link{pairwise_npt_distance}} or downstream regression models.
#'
#' @param data A list of \eqn{n} numeric matrices, where the \eqn{i}th
#'   \eqn{N_i \times d} matrix represents an empirical distribution (sample).
#'   Rows are observations, so sample sizes \eqn{N_i} may vary; columns are the
#'   \eqn{d} common variables.
#' @param M Number of equal-width probability intervals whose midpoints
#'   are used to evaluate marginal quantiles (default is `100L`).
#' @param pd_shrinkage Shrinkage intensity in \eqn{(0, 1)} ensuring the
#'   estimated latent correlation matrices are strictly positive definite via
#'   \eqn{(1 - \lambda) R + \lambda I} (default is `0.001`).
#' @param cache_sqrt Logical; if `TRUE` (default), precomputes and caches
#'   matrix square roots \eqn{R^{1/2}} to speed up pairwise distance calculations
#'   when dimension \eqn{d > 2}. For \eqn{d = 2}, roots are not computed or
#'   cached, even when `TRUE`, because bivariate distances use scalar correlations.
#' @param standardize Logical; if `TRUE`, preprocesses the pooled collection for
#'   distance computation by subtracting each variable's pooled median
#'   and dividing by its pooled MAD from `stats::mad()`. The default is `FALSE`;
#'   regression or other methods therefore use original measurement units
#'   unless standardization is explicitly requested.
#'
#' @return An object of class \code{"nonparanormal"}, containing:
#' \describe{
#'   \item{\code{quantiles}}{A list of \eqn{d} matrices (one per variable), each
#'     of dimension \eqn{n \times M}, holding the evaluated marginal quantile functions.}
#'   \item{\code{correlations}}{A list of \eqn{n} positive-definite \eqn{d \times d}
#'     latent correlation matrices.}
#'   \item{\code{correlation_sqrts}}{A list of \eqn{n} matrix square roots \eqn{R^{1/2}}
#'     if \code{cache_sqrt = TRUE} and \eqn{d > 2}, or \code{NULL} otherwise.}
#'   \item{\code{cache_sqrt}}{The effective caching setting; always \code{FALSE}
#'     for bivariate data.}
#'   \item{\code{probabilities}}{Numeric vector containing the \eqn{M}
#'     probability-interval midpoints in \eqn{(0, 1)}.}
#'   \item{\code{sample_sizes}}{Named integer vector of sample sizes \eqn{N_i} across distributions.}
#'   \item{\code{tied_pair_fractions}}{An \eqn{n \times d} matrix containing the fraction
#'     of tied pairs per variable and distribution.}
#'   \item{\code{distribution_names}}{Character vector of distribution names.}
#'   \item{\code{variable_names}}{Character vector of variable names.}
#'   \item{\code{standardization}}{Preprocessing metadata recording whether
#'     pooled median/MAD standardization was applied and, when applicable, its
#'     named center and scale vectors.}
#' }
#'
#' @seealso \code{\link{pairwise_npt_distance}}, \code{\link{npt_frechetreg}}
#' @export
as_nonparanormal <- function(
  data,
  M = 100L,
  pd_shrinkage = 0.001,
  cache_sqrt = TRUE,
  standardize = FALSE
) {
  # Validate input matrices and normalize parameters
  validated <- .validate_distributional_data(data)
  M <- .as_integer_count(M, "M", minimum = 2L)
  pd_shrinkage <- .validate_pd_shrinkage(pd_shrinkage)
  cache_sqrt <- .as_flag(cache_sqrt, "cache_sqrt")
  standardize <- .as_flag(standardize, "standardize")

  # Distance computation can use one common robust transformation across all
  # distributions; the default retains the original measurement units.
  standardization <- list(
    applied = FALSE,
    method = "none",
    center = NULL,
    scale = NULL
  )
  if (standardize) {
    transformed <- .standardize_pooled_median_mad(
      validated$data,
      validated$variable_names
    )
    validated$data <- transformed$data
    standardization <- list(
      applied = TRUE,
      method = "pooled_median_mad",
      center = transformed$center,
      scale = transformed$scale
    )
  }

  # Precompute marginal quantiles and latent correlations in C++
  summaries <- distribution_summaries_cpp(
    validated$data,
    M,
    pd_shrinkage,
    cache_sqrt
  )

  # Attach variable and distribution names to output components
  names(summaries$quantiles) <- validated$variable_names
  for (j in seq_along(summaries$quantiles)) {
    rownames(summaries$quantiles[[j]]) <- validated$distribution_names
  }

  names(summaries$correlations) <- validated$distribution_names
  for (i in seq_along(summaries$correlations)) {
    dimnames(summaries$correlations[[i]]) <- list(
      validated$variable_names,
      validated$variable_names
    )
  }

  if (!is.null(summaries$correlation_sqrts)) {
    names(summaries$correlation_sqrts) <- validated$distribution_names
    for (i in seq_along(summaries$correlation_sqrts)) {
      dimnames(summaries$correlation_sqrts[[i]]) <- list(
        validated$variable_names,
        validated$variable_names
      )
    }
  }

  names(summaries$sample_sizes) <- validated$distribution_names
  rownames(summaries$tied_pair_fractions) <- validated$distribution_names
  colnames(summaries$tied_pair_fractions) <- validated$variable_names

  # Attach metadata and assign S3 class
  summaries$distribution_names <- validated$distribution_names
  summaries$variable_names <- validated$variable_names
  summaries$pd_shrinkage <- pd_shrinkage
  summaries$cache_sqrt <- !is.null(summaries$correlation_sqrts)
  summaries$standardization <- standardization

  structure(summaries, class = "nonparanormal")
}

#' @export
print.nonparanormal <- function(x, ...) {
  cat("<nonparanormal distributions>\n")
  cat("Distributions:", length(x$correlations), "\n")
  cat("Variables:", length(x$quantiles), "\n")
  cat("Quantile grid points:", length(x$probabilities), "\n")
  cat(
    "Pooled median/MAD standardization: ",
    if (isTRUE(x$standardization$applied)) "yes" else "no",
    "\n",
    sep = ""
  )
  cat("Cached correlation square roots:", if (x$cache_sqrt) "yes" else "no", "\n")
  invisible(x)
}
