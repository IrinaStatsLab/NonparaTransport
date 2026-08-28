// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <vector>

namespace {

// =============================================================================
// Pairwise Nonparanormal Transport (NPT) Distance Matrix (C++ Engine)
//
// Computes the pairwise squared NPT distance matrix between n distributions
// P_1, ..., P_n using precomputed summaries:
//
//   d_{NPT}^2(P_i, P_k) = sum_{j=1}^d ||Q_i^{(j)} - Q_k^{(j)}||^2_{L^2} +
//   \mathcal{B}^2(Sigma_i, Sigma_k)
//
// Components:
//   1. Marginal L^2 distance:
//        sum_{j=1}^d ||Q_i^{(j)} - Q_k^{(j)}||^2_{L^2}
//      evaluated via midpoint numerical integration over M grid points.
//   2. Correlation distance (Bures-Wasserstein):
//        \mathcal{B}^2(Sigma_i, Sigma_k) = tr(Sigma_i) + tr(Sigma_k) - 2 tr{
//        (Sigma_i^{1/2} Sigma_k Sigma_i^{1/2})^{1/2} }
//      For d = 2, simplified to an exact closed-form scalar expression.
//
// Reading Guide (Execution Flow):
//   The main exported function is `pairwise_summary_distance_cpp()` at the
//   bottom of this file. It evaluates unordered pairs (i, k) for i < k and
//   mirrors the results to construct symmetric n x n distance matrices:
//     1. Marginal loop: evaluates L^2 quantile differences across all d
//     variables.
//     2. Dependence loop:
//        - If d = 2: calls squared_bivariate_bures() for scalar fast evaluation.
//        - If d > 2: calls squared_bures_wasserstein() using cached Sigma^{1/2} roots
//          and trace_symmetric_sqrt() for eigenvalue-only trace evaluation.
// =============================================================================

// Principal symmetric square root used when summaries do not cache Sigma^{1/2}
arma::mat symmetric_sqrt(const arma::mat &matrix) {
  arma::mat result;
  const arma::mat symmetric = 0.5 * (matrix + matrix.t());
  if (!arma::sqrtmat_sympd(result, symmetric)) {
    Rcpp::stop("Unable to compute a positive-definite matrix square root.");
  }
  return 0.5 * (result + result.t());
}

// For the Bures--Wasserstein middle matrix A, tr(A^{1/2}) = sum_j sqrt(lambda_j(A)); this
// scalar term needs eigenvalues, not the reconstructed square-root matrix.
double trace_symmetric_sqrt(const arma::mat &matrix) {
  arma::vec eigenvalues;
  const arma::mat symmetric = 0.5 * (matrix + matrix.t());
  if (!arma::eig_sym(eigenvalues, symmetric)) {
    Rcpp::stop("Unable to eigendecompose a positive-definite matrix.");
  }
  eigenvalues.transform([](double value) { return std::max(0.0, value); });
  return arma::sum(arma::sqrt(eigenvalues));
}

// Squared Bures-Wasserstein distance between two d x d latent correlation
// matrices:
//   \mathcal{B}^2(Sigma_P, Sigma_Q) = tr(Sigma_P) + tr(Sigma_Q) - 2 * tr{
//   (Sigma_P^{1/2} Sigma_Q Sigma_P^{1/2})^{1/2} }
// `first_sqrt` (Sigma_P^{1/2}) is precomputed to avoid redundant O(d^3) roots.
double squared_bures_wasserstein(const arma::mat &first,
                                 const arma::mat &second,
                                 const arma::mat &first_sqrt) {
  const arma::mat middle = first_sqrt * second * first_sqrt;
  const double distance =
      arma::trace(first) + arma::trace(second) -
      2.0 * trace_symmetric_sqrt(middle);

  // Numerical precision guard: clamp tiny negative roundoff errors to zero
  return std::max(0.0, distance);
}

// Closed-form squared Bures-Wasserstein distance for bivariate correlations (d
// = 2). For 2x2 correlation matrices with parameters rho_P and rho_Q,
// eigenvalues are (1 + rho) and (1 - rho) with shared eigenvectors, simplifying
// the distance to:
//   4 - 2 * [ sqrt((1 + rho_P)(1 + rho_Q)) + sqrt((1 - rho_P)(1 - rho_Q)) ]
double squared_bivariate_bures(double first_correlation,
                               double second_correlation) {
  const double positive =
      std::sqrt((1.0 + first_correlation) * (1.0 + second_correlation));
  const double negative =
      std::sqrt((1.0 - first_correlation) * (1.0 - second_correlation));
  return std::max(0.0, 4.0 - 2.0 * (positive + negative));
}

} // namespace

// [[Rcpp::export]]
Rcpp::List pairwise_summary_distance_cpp(Rcpp::List quantiles,
                                         Rcpp::List correlations,
                                         Rcpp::List cached_sqrts) {
  const int d = quantiles.size();
  const int n = correlations.size();
  if (d < 1 || n < 1) {
    Rcpp::stop("Compiled distance input is invalid.");
  }

  // Unpack R lists into C++ Armadillo structures:
  // - quantile_matrices: d matrices of dimension (n x M)
  // - correlation_matrices: n matrices of dimension (d x d)
  std::vector<arma::mat> quantile_matrices(d);
  for (int variable = 0; variable < d; ++variable) {
    quantile_matrices[variable] = Rcpp::as<arma::mat>(quantiles[variable]);
  }

  std::vector<arma::mat> correlation_matrices(n);
  for (int i = 0; i < n; ++i) {
    correlation_matrices[i] = Rcpp::as<arma::mat>(correlations[i]);
  }

  // For d > 2, prepare one matrix square root Sigma^{1/2} per distribution.
  // Use cached roots if available, or compute them on demand.
  std::vector<arma::mat> square_roots;
  if (d > 2) {
    square_roots.resize(n);
    const bool use_cache = (cached_sqrts.size() == n);
    for (int i = 0; i < n; ++i) {
      square_roots[i] = use_cache ? Rcpp::as<arma::mat>(cached_sqrts[i])
                                  : symmetric_sqrt(correlation_matrices[i]);
    }
  }

  // Allocate symmetric n x n distance matrices
  arma::mat marginal(n, n, arma::fill::zeros);
  arma::mat correlation(n, n, arma::fill::zeros);

  // 1. Marginal component: sum of L^2 squared quantile differences across all d
  // variables. Since quantiles are evaluated on M equal-width midpoint cells,
  // the L^2 integral is the mean squared difference across the M columns.
  for (int variable = 0; variable < d; ++variable) {
    arma::mat &current = quantile_matrices[variable];

    // Subtract the same probability-wise mean curve from every distribution.
    // Pairwise differences are unchanged while controlling numerical stability.
    const arma::rowvec probability_means = arma::mean(current, 0);
    current.each_row() -= probability_means;

    // For rows q_i and q_k, ||q_i-q_k||^2/M equals
    // (||q_i||^2 + ||q_k||^2 - 2 q_i^T q_k)/M. One matrix multiplication
    // computes q_i^T q_k for every pair, replacing the nested pair loop.
    const double inverse_M = 1.0 / static_cast<double>(current.n_cols);
    const arma::vec squared_norms =
        inverse_M * arma::sum(arma::square(current), 1);
    marginal += (-2.0 * inverse_M) * (current * current.t());
    marginal.each_col() += squared_norms;
    marginal.each_row() += squared_norms.t();
  }

  // Copy one computed triangle to the other and remove negative roundoff from
  // the expanded squared-distance identity; true squared distances are >= 0.
  marginal = arma::symmatu(marginal);
  marginal.transform([](double value) { return std::max(0.0, value); });
  marginal.diag().zeros();

  // 2. Correlation component: Bures-Wasserstein distance between latent
  // correlation matrices.
  // - d = 2: exact closed-form scalar formula
  // - d > 2: general matrix Bures-Wasserstein formula using cached Sigma^{1/2}
  if (d == 2) {
    for (int i = 0; i < n; ++i) {
      const double rho_i = correlation_matrices[i](0, 1);
      for (int k = i + 1; k < n; ++k) {
        const double value =
            squared_bivariate_bures(rho_i, correlation_matrices[k](0, 1));
        correlation(i, k) = value;
        correlation(k, i) = value;
      }
    }
  } else if (d > 2) {
    for (int i = 0; i < n; ++i) {
      for (int k = i + 1; k < n; ++k) {
        const double value = squared_bures_wasserstein(
            correlation_matrices[i], correlation_matrices[k], square_roots[i]);
        correlation(i, k) = value;
        correlation(k, i) = value;
      }
    }
  }

  return Rcpp::List::create(Rcpp::Named("marginal") = marginal,
                            Rcpp::Named("correlation") = correlation,
                            Rcpp::Named("total") = marginal + correlation);
}
