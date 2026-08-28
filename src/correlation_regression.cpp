// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace {

// =============================================================================
// Latent Gaussian Correlation Regression & Component Losses (C++ Engine)
//
// Computes global Fréchet regression for latent Gaussian correlation matrices
// Sigma_1, ..., Sigma_n under the Bures--Wasserstein metric, and evaluates
// component-wise sum of squared Fréchet losses for goodness-of-fit (R^2).
//
// Theoretical Formulation:
//   Given training predictors X (n x p), prediction points Z (k x p), and
//   training latent correlation matrices Sigma_1, ..., Sigma_n (d x d):
//
//   1. Global Fréchet Weights:
//        w(z) = (1/n) 1 + (Z_c (X_c^T X_c)^+ X_c^T)^T
//      where X_c, Z_c are centered and standardized by training sample moments.
//
//   2. Correlation Barycenter / Regression:
//      - For d = 2: Closed-form scalar solution for rho(z) in [-1, 1].
//      - For d > 2: Projected Bures--Wasserstein Riemannian gradient descent:
//          T_i(Sigma) = Sigma^{-1/2}
//                       (Sigma^{1/2} Sigma_i Sigma^{1/2})^{1/2} Sigma^{-1/2}
//          G         = sum_{i=1}^n w_i(z) T_i(Sigma)
//          C         = G Sigma G^T
//          Sigma_new = diag(C)^{-1/2} C diag(C)^{-1/2}  (projected to correlation
//          cone)
//
//   3. Component Losses:
//      - Marginal: sum_{i=1}^n ||Q_i^{(j)} - Qhat^{(j)}(z_i)||^2_{L^2}
//      - Latent Correlation: sum_{i=1}^n \mathcal{B}^2(Sigma_i, Sigmahat(z_i))
//
//   1. Helper Linear Algebra / Manifold Operations:
//      - symmetric_roots()        -> Simultaneous eigendecomposition of
//                                    Sigma^{1/2} & Sigma^{-1/2}
//      - symmetric_square_root()  -> Full symmetric PSD square root A^{1/2}
//      - trace_symmetric_sqrt()   -> Direct scalar trace tr(A^{1/2}) = sum sqrt(lambda_j)
//      - project_to_correlation() -> Normalizes covariance to unit diagonal
//      - global_frechet_weights() -> Computes normalized Fréchet weights w(z)
//      - bivariate_barycenter()   -> Exact scalar closed form for d = 2
//      - integrated_squared_diff()-> Midpoint L^2 quantile loss
//      - squared_bures_wasserst() -> \mathcal{B}^2 distance between correlation
//                                    matrices
//   2. Exported Entry Points:
//      - correlation_regression_cpp()   -> Predicts Sigmahat(z) for all prediction
//      points
//      - matched_component_losses_cpp() -> Computes in-sample / baseline loss
//      sums
// =============================================================================

// Holds both Sigma^{1/2} and Sigma^{-1/2} computed from a single eigendecomposition
struct SymmetricRoots {
  arma::mat square_root;
  arma::mat inverse_square_root;
};

// Computes Sigma^{1/2} and Sigma^{-1/2} with eigenvalue thresholding for numerical
// stability
SymmetricRoots symmetric_roots(const arma::mat &matrix, double eigen_floor) {
  arma::vec eigenvalues;
  arma::mat eigenvectors;
  const arma::mat symmetric = 0.5 * (matrix + matrix.t());
  if (!arma::eig_sym(eigenvalues, eigenvectors, symmetric)) {
    Rcpp::stop("Unable to eigendecompose a correlation iterate.");
  }

  eigenvalues.transform(
      [eigen_floor](double value) { return std::max(value, eigen_floor); });
  const arma::vec square_roots = arma::sqrt(eigenvalues);

  SymmetricRoots result;
  result.square_root =
      eigenvectors * arma::diagmat(square_roots) * eigenvectors.t();
  result.inverse_square_root =
      eigenvectors * arma::diagmat(1.0 / square_roots) * eigenvectors.t();
  return result;
}

// Computes the symmetric positive semidefinite square root A^{1/2}
arma::mat symmetric_square_root(const arma::mat &matrix, double eigen_floor) {
  arma::vec eigenvalues;
  arma::mat eigenvectors;
  const arma::mat symmetric = 0.5 * (matrix + matrix.t());
  if (!arma::eig_sym(eigenvalues, eigenvectors, symmetric)) {
    Rcpp::stop("Unable to eigendecompose a positive-semidefinite matrix.");
  }
  eigenvalues.transform(
      [eigen_floor](double value) { return std::max(value, eigen_floor); });
  return eigenvectors * arma::diagmat(arma::sqrt(eigenvalues)) *
         eigenvectors.t();
}

// Computes tr(A^{1/2}) = tr(V Lambda^{1/2} V^T) = sum_j sqrt(lambda_j(A)).
// Avoids reconstructing the full O(d^3) square-root matrix when only its trace is needed.
double trace_symmetric_sqrt(const arma::mat &matrix, double eigen_floor) {
  arma::vec eigenvalues;
  const arma::mat symmetric = 0.5 * (matrix + matrix.t());
  if (!arma::eig_sym(eigenvalues, symmetric)) {
    Rcpp::stop("Unable to eigendecompose a positive-semidefinite matrix.");
  }
  eigenvalues.transform(
      [eigen_floor](double value) { return std::max(value, eigen_floor); });
  return arma::sum(arma::sqrt(eigenvalues));
}

// Projects a covariance matrix C to the correlation manifold: D^{-1/2} C
// D^{-1/2}
arma::mat project_to_correlation(const arma::mat &covariance,
                                 double diagonal_floor) {
  const arma::mat symmetric = 0.5 * (covariance + covariance.t());
  const arma::vec diagonal = symmetric.diag();
  if (arma::any(diagonal <= diagonal_floor) || !diagonal.is_finite()) {
    Rcpp::stop(
        "A projected BW update has a non-positive diagonal; the weighted "
        "correlation objective is numerically degenerate at this prediction "
        "point.");
  }

  const arma::vec scales = arma::sqrt(diagonal);
  arma::mat correlation = symmetric;
  correlation.each_col() /= scales;
  correlation.each_row() /= scales.t();
  correlation = 0.5 * (correlation + correlation.t());
  correlation.diag().ones();
  return correlation;
}

// Computes global Fréchet regression weights w(z) across all prediction points:
//   w(z)^T = (1/n) 1^T + z_tilde^T (X_tilde^T X_tilde)^+ X_tilde^T
arma::mat global_frechet_weights(const arma::mat &X, const arma::mat &Z,
                                 double scale_floor) {
  // Standardize predictors by training sample mean and variance
  const arma::rowvec center = arma::mean(X, 0);
  arma::mat scaled_X = X.each_row() - center;
  arma::rowvec scales = arma::sqrt(arma::mean(arma::square(scaled_X), 0));
  scales.transform(
      [scale_floor](double value) { return std::max(value, scale_floor); });
  scaled_X.each_row() /= scales;

  arma::mat scaled_Z = Z.each_row() - center;
  scaled_Z.each_row() /= scales;

  // Compute k x n weight matrix via Moore-Penrose pseudoinverse
  const arma::mat gram_inverse = arma::pinv(scaled_X.t() * scaled_X);
  arma::mat weights = scaled_Z * gram_inverse * scaled_X.t();
  weights += 1.0 / static_cast<double>(X.n_rows);

  // Normalize row sums to exactly 1.0 to eliminate floating-point drift
  const arma::vec row_sums = arma::sum(weights, 1);
  if (arma::any(arma::abs(row_sums) <= scale_floor) || !row_sums.is_finite()) {
    Rcpp::stop("Unable to normalize global Frechet regression weights.");
  }
  weights.each_col() /= row_sums;
  return weights;
}

// Closed-form correlation barycenter for bivariate distributions (d = 2):
// Maximizes a * sqrt(1 + rho) + b * sqrt(1 - rho) over rho in [-1, 1]
arma::mat bivariate_barycenter(const std::vector<arma::mat> &samples,
                               const arma::rowvec &weights,
                               bool &constant_objective) {
  double positive = 0.0;
  double negative = 0.0;
  for (arma::uword i = 0; i < weights.n_elem; ++i) {
    const double rho = std::max(-1.0, std::min(1.0, samples[i](0, 1)));
    positive += weights[i] * std::sqrt(1.0 + rho);
    negative += weights[i] * std::sqrt(1.0 - rho);
  }

  double rho;
  constant_objective = false;
  if (positive < 0.0 && negative >= 0.0) {
    // Strictly decreasing on [-1, 1] -> minimizer at -1
    rho = -1.0;
  } else if (negative < 0.0 && positive >= 0.0) {
    // Strictly increasing on [-1, 1] -> minimizer at +1
    rho = 1.0;
  } else if (positive < 0.0 && negative < 0.0) {
    // Both negative -> compare boundary endpoints
    rho = positive > negative ? 1.0 : -1.0;
  } else {
    // Nonnegative coefficients -> stationary point rho = (a^2 - b^2) / (a^2 +
    // b^2)
    const double denominator = positive * positive + negative * negative;
    if (denominator <= 1e-30) {
      rho = samples.front()(0, 1);
      constant_objective = true;
    } else {
      rho = (positive * positive - negative * negative) / denominator;
    }
  }

  arma::mat result(2, 2, arma::fill::eye);
  result(0, 1) = rho;
  result(1, 0) = rho;
  return result;
}

// Computes integrated squared difference between two quantile curves via
// midpoint grid
double integrated_squared_difference(const arma::rowvec &first,
                                     const arma::rowvec &second) {
  const arma::uword M = first.n_elem;
  if (M < 2 || second.n_elem != M) {
    Rcpp::stop("Quantile matrices must share at least two grid points.");
  }
  const arma::rowvec squared = arma::square(first - second);
  return arma::mean(squared);
}

// Computes squared Bures--Wasserstein distance between two d x d correlation
// matrices
double squared_bures_wasserstein(const arma::mat &first,
                                 const arma::mat &second, double eigen_floor) {
  const arma::uword d = first.n_rows;
  if (d == 2) {
    const double first_rho = std::max(-1.0, std::min(1.0, first(0, 1)));
    const double second_rho = std::max(-1.0, std::min(1.0, second(0, 1)));
    const double value =
        4.0 - 2.0 * (std::sqrt((1.0 + first_rho) * (1.0 + second_rho)) +
                     std::sqrt((1.0 - first_rho) * (1.0 - second_rho)));
    return std::max(0.0, value);
  }

  // General d > 2: tr(Sigma_1) + tr(Sigma_2) -
  // 2 tr((Sigma_1^{1/2} Sigma_2 Sigma_1^{1/2})^{1/2})
  const arma::mat first_sqrt = symmetric_square_root(first, eigen_floor);
  const double middle_sqrt_trace =
      trace_symmetric_sqrt(first_sqrt * second * first_sqrt, eigen_floor);
  const double value =
      arma::trace(first) + arma::trace(second) - 2.0 * middle_sqrt_trace;
  return std::max(0.0, value);
}

} // namespace

// =============================================================================
// Exported Functions
// =============================================================================

// [[Rcpp::export]]
Rcpp::List correlation_regression_cpp(Rcpp::List correlations,
                                      const arma::mat &X, const arma::mat &Z,
                                      int max_iter, double tolerance,
                                      double eigen_floor) {
  const int n = correlations.size();
  if (n < 1 || X.n_rows != static_cast<arma::uword>(n) ||
      X.n_cols != Z.n_cols || Z.n_rows < 1) {
    Rcpp::stop(
        "Compiled correlation-regression input has inconsistent dimensions.");
  }

  // Convert R list to standard vector of correlation matrices
  std::vector<arma::mat> samples(n);
  arma::uword d = 0;
  for (int i = 0; i < n; ++i) {
    samples[i] = Rcpp::as<arma::mat>(correlations[i]);
    if (i == 0) {
      d = samples[i].n_rows;
    }
    if (samples[i].n_rows != d || samples[i].n_cols != d) {
      Rcpp::stop("All latent correlation matrices must have the same square "
                 "dimension.");
    }
  }

  if (d < 2) {
    Rcpp::stop("Latent correlation regression requires at least two variables (d >= 2).");
  }

  // Precompute global Fréchet weights w(z) for all k prediction points (k x n)
  const arma::mat weights = global_frechet_weights(X, Z, 1e-10);
  const int number_predictions = Z.n_rows;
  Rcpp::List predictions(number_predictions);
  Rcpp::IntegerVector iterations(number_predictions);
  Rcpp::LogicalVector converged(number_predictions);
  Rcpp::NumericVector final_change(number_predictions);
  Rcpp::CharacterVector status(number_predictions);

  for (int prediction = 0; prediction < number_predictions; ++prediction) {
    // Case 1: Bivariate (d = 2) -> closed-form scalar barycenter
    if (d == 2) {
      bool constant_objective = false;
      const arma::mat result = bivariate_barycenter(
          samples, weights.row(prediction), constant_objective);
      predictions[prediction] = result;
      iterations[prediction] = 1;
      converged[prediction] = true;
      final_change[prediction] = arma::norm(result - samples.front(), "fro");
      status[prediction] =
          constant_objective ? "constant_objective" : "closed_form";
      continue;
    }

    // Case 2: Multivariate (d > 2) -> Projected Riemannian Gradient Descent
    arma::mat current = samples.front();
    bool did_converge = false;
    double change = NA_REAL;
    int completed_iterations = 0;

    for (int iteration = 1; iteration <= max_iter; ++iteration) {
      const SymmetricRoots current_roots =
          symmetric_roots(current, eigen_floor);
      arma::mat middle_average(d, d, arma::fill::zeros);

      // For U = Sigma^{-1/2} and
      // H_i = (Sigma^{1/2} Sigma_i Sigma^{1/2})^{1/2}, U is common,
      // so sum_i w_i U H_i U = U (sum_i w_i H_i) U. This cuts 2n matrix
      // products to 2 per iteration; each nonlinear H_i stays inside the loop.
      for (int i = 0; i < n; ++i) {
        const arma::mat middle_sqrt = symmetric_square_root(
            current_roots.square_root * samples[i] * current_roots.square_root,
            eigen_floor);
        middle_average += weights(prediction, i) * middle_sqrt;
      }
      const arma::mat transport_average = current_roots.inverse_square_root *
                                          middle_average *
                                          current_roots.inverse_square_root;

      // Covariance update with step size 1 and projection back to correlation
      // manifold
      const arma::mat covariance =
          transport_average * current * transport_average.t();
      const arma::mat updated = project_to_correlation(covariance, eigen_floor);
      change = arma::norm(updated - current, "fro");
      current = updated;
      completed_iterations = iteration;

      if (change < tolerance) {
        did_converge = true;
        break;
      }
    }

    predictions[prediction] = current;
    iterations[prediction] = completed_iterations;
    converged[prediction] = did_converge;
    final_change[prediction] = change;
    status[prediction] = did_converge ? "converged" : "max_iterations";
  }

  return Rcpp::List::create(Rcpp::Named("correlations") = predictions,
                            Rcpp::Named("iterations") = iterations,
                            Rcpp::Named("converged") = converged,
                            Rcpp::Named("final_change") = final_change,
                            Rcpp::Named("status") = status);
}

// [[Rcpp::export]]
Rcpp::List matched_component_losses_cpp(Rcpp::List observed_quantiles,
                                        Rcpp::List predicted_quantiles,
                                        Rcpp::List observed_correlations,
                                        Rcpp::List predicted_correlations,
                                        double eigen_floor) {
  const int d = observed_quantiles.size();
  const int n = observed_correlations.size();
  if (d < 1 || n < 1 || predicted_quantiles.size() != d) {
    Rcpp::stop("Compiled component-loss input has inconsistent dimensions.");
  }

  // 1. Marginal quantile losses across all d variables
  Rcpp::NumericVector marginal(d);
  for (int variable = 0; variable < d; ++variable) {
    const arma::mat observed =
        Rcpp::as<arma::mat>(observed_quantiles[variable]);
    const arma::mat predicted =
        Rcpp::as<arma::mat>(predicted_quantiles[variable]);
    if (observed.n_rows != static_cast<arma::uword>(n) ||
        observed.n_cols != predicted.n_cols ||
        (predicted.n_rows != 1 && predicted.n_rows != observed.n_rows)) {
      Rcpp::stop("Observed and predicted quantile matrices are not aligned.");
    }

    double loss = 0.0;
    for (int i = 0; i < n; ++i) {
      const arma::uword prediction_row = predicted.n_rows == 1 ? 0 : i;
      loss += integrated_squared_difference(observed.row(i),
                                            predicted.row(prediction_row));
    }
    marginal[variable] = loss;
  }

  if (predicted_correlations.size() != 1 &&
      predicted_correlations.size() != n) {
    Rcpp::stop("Observed and predicted latent correlations are not aligned.");
  }

  // 2. Latent correlation Bures--Wasserstein loss
  double correlation_loss = 0.0;
  for (int i = 0; i < n; ++i) {
    const arma::mat observed = Rcpp::as<arma::mat>(observed_correlations[i]);
    const int prediction_index = predicted_correlations.size() == 1 ? 0 : i;
    const arma::mat predicted =
        Rcpp::as<arma::mat>(predicted_correlations[prediction_index]);
    correlation_loss +=
        squared_bures_wasserstein(observed, predicted, eigen_floor);
  }

  return Rcpp::List::create(Rcpp::Named("marginal") = marginal,
                            Rcpp::Named("correlation") = correlation_loss);
}
