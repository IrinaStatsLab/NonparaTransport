// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

namespace {

// =============================================================================
// Nonparanormal Distribution Summaries (C++ Engine)
//
// For each N_i x d empirical sample matrix, this file precomputes:
//   1. d marginal quantile functions evaluated on an M-point grid.
//   2. A d x d latent Gaussian copula correlation matrix (via Kendall's tau-a,
//      the sine bridge, and projection to the correlation cone).
//   3. (Optional) Matrix square roots R^{1/2} for fast distance calculations.
//
// Key optimization:
//   Each column is sorted once in O(N_i log N_i). The sorted values and
//   permutation orders are reused for both quantile interpolation and fast
//   Kendall's tau-a across all pairs of variables.
//
// Reading Guide (Execution Flow):
//   The main exported entry point is `distribution_summaries_cpp()` at the
//   bottom of this file. It orchestrates the internal helper functions in
//   order:
//     1. sort_column()            -> Single-pass sorting and tie-block
//     detection
//     2. quantiles_type7()        -> Type-7 empirical quantile evaluation
//     3. correlation_from_sorted()-> Pairwise Kendall's tau-a and sine bridge
//        `- merge_and_count()     -> O(N_i log N_i) discordant pair counting
//     4. nearest_correlation()    -> Dykstra PSD projection + pd_shrinkage
//     5. symmetric_sqrt()         -> Optional R^{1/2} caching
// =============================================================================

// 64-bit integer for pair counts (N_i*(N_i-1)/2 grows quadratically)
using count_t = std::int64_t;

// Half-open interval [begin, end) representing a run of identical values (ties)
struct TieBlock {
  arma::uword begin;
  arma::uword end; // Index one past the last tied element
};

// Cached sort information for a single column/variable
struct SortedColumn {
  arma::vec values; // Sorted observations (for quantile interpolation)
  arma::uvec order; // Permutation index from original to sorted order
  std::vector<TieBlock>
      tie_blocks;         // Contiguous intervals with >= 2 identical values
  count_t tied_pairs = 0; // Total pairs tied in this column: sum choose(len, 2)
};

// Binomial coefficient choose(N_i, 2) = N_i * (N_i - 1) / 2
inline count_t choose_two(arma::uword count) {
  return static_cast<count_t>(count) * static_cast<count_t>(count - 1) / 2;
}

// Sort a column and record runs of tied values in a single O(N_i log N_i) pass
SortedColumn sort_column(const arma::vec &x) {
  SortedColumn result;

  // Use stable sort to preserve a deterministic ordering
  result.order = arma::stable_sort_index(x);
  result.values = x.elem(result.order);

  // Identify contiguous runs of identical values (ties) and count tied pairs
  arma::uword begin = 0;
  while (begin < result.values.n_elem) {
    arma::uword end = begin + 1;
    while (end < result.values.n_elem &&
           result.values[end] == result.values[begin]) {
      ++end;
    }
    const arma::uword block_size = end - begin;
    if (block_size > 1) {
      result.tie_blocks.push_back({begin, end});
      result.tied_pairs += choose_two(block_size);
    }
    begin = end;
  }

  return result;
}

// Count inversions in O(N_i log N_i) time during merge-sort.
// When observations are ordered by x, an inversion in y corresponds to a
// strictly discordant pair (x_i < x_j and y_i > y_j).
count_t merge_and_count(arma::vec &values, arma::vec &work, arma::uword left,
                        arma::uword right) {
  if (left >= right) {
    return 0;
  }

  const arma::uword middle = left + (right - left) / 2;
  count_t inversions = merge_and_count(values, work, left, middle) +
                       merge_and_count(values, work, middle + 1, right);

  arma::uword i = left;
  arma::uword j = middle + 1;
  arma::uword k = left;
  while (i <= middle && j <= right) {
    if (values[i] <= values[j]) {
      work[k++] = values[i++];
    } else {
      work[k++] = values[j++];
      inversions += static_cast<count_t>(middle - i + 1);
    }
  }
  while (i <= middle) {
    work[k++] = values[i++];
  }
  while (j <= right) {
    work[k++] = values[j++];
  }
  for (arma::uword index = left; index <= right; ++index) {
    values[index] = work[index];
  }

  return inversions;
}

// Fast Kendall's tau-a between two pre-sorted columns.
// Definition: tau_a = (Concordant - Discordant) / choose(N_i, 2)
//
// Pairs tied in x, tied in y, or tied in both contribute 0 to concordance.
// Using inclusion-exclusion:
//   Comparable pairs = Total - Tied_x - Tied_y + Tied_both
//   Concordant = Comparable - Discordant
//   Concordant - Discordant = Comparable - 2 * Discordant
double kendall_tau_a_from_sorted(const SortedColumn &x, const SortedColumn &y,
                                 const arma::vec &y_original) {
  const arma::uword sample_size = x.values.n_elem;
  const count_t total_pairs = choose_two(sample_size);
  arma::vec y_in_x_order = y_original.elem(x.order);
  count_t tied_in_both = 0;

  // Reorder y according to x's permutation.
  // Within blocks where x is tied, sort y to prevent false discordance counts
  // and simultaneously count pairs that are tied in both variables.
  for (const TieBlock &block : x.tie_blocks) {
    std::sort(y_in_x_order.memptr() + block.begin,
              y_in_x_order.memptr() + block.end);

    arma::uword begin = block.begin;
    while (begin < block.end) {
      arma::uword end = begin + 1;
      while (end < block.end && y_in_x_order[end] == y_in_x_order[begin]) {
        ++end;
      }
      tied_in_both += choose_two(end - begin);
      begin = end;
    }
  }

  // Count strictly discordant pairs via O(N_i log N_i) merge-sort inversion
  // count
  arma::vec work(sample_size);
  const count_t discordant =
      merge_and_count(y_in_x_order, work, 0, sample_size - 1);

  // Net concordance score
  const count_t score =
      total_pairs - x.tied_pairs - y.tied_pairs + tied_in_both - 2 * discordant;

  const double tau =
      static_cast<double>(score) / static_cast<double>(total_pairs);
  return std::max(-1.0, std::min(1.0, tau));
}

// Midpoint of the k-th equal-width probability cell: (k + 0.5) / M.
// Evaluates quantiles at cell centers in (0, 1), avoiding boundary extremes (0
// and 1).
inline double midpoint_probability(int index, int M) {
  return (static_cast<double>(index) + 0.5) / static_cast<double>(M);
}

// Compute empirical quantiles using R's default type-7 definition:
// h = (N_i - 1) * p, with linear interpolation between floor(h) and ceil(h).
arma::vec quantiles_type7(const arma::vec &sorted_values, int M) {
  const arma::uword sample_size = sorted_values.n_elem;
  arma::vec result(M);
  for (int k = 0; k < M; ++k) {
    const double probability = midpoint_probability(k, M);
    const double index = probability * static_cast<double>(sample_size - 1);
    const arma::uword lower = static_cast<arma::uword>(std::floor(index));
    const arma::uword upper = static_cast<arma::uword>(std::ceil(index));
    const double weight = index - static_cast<double>(lower);
    // Difference-form interpolation preserves exact ties and avoids roundoff
    // that can make theoretically nondecreasing type-7 quantiles decrease.
    result[k] = sorted_values[lower] +
                weight * (sorted_values[upper] - sorted_values[lower]);
  }
  return result;
}

// Orthogonal projection onto the positive semi-definite (PSD) cone under
// Frobenius norm: Performs eigenvalue decomposition and clamps negative
// eigenvalues to zero.
arma::mat project_psd(const arma::mat &matrix) {
  arma::vec eigenvalues;
  arma::mat eigenvectors;
  if (!arma::eig_sym(eigenvalues, eigenvectors, 0.5 * (matrix + matrix.t()))) {
    Rcpp::stop(
        "Eigenvalue decomposition failed during correlation projection.");
  }
  eigenvalues.transform([](double value) { return std::max(0.0, value); });
  return eigenvectors * arma::diagmat(eigenvalues) * eigenvectors.t();
}

// Nearest positive-definite correlation matrix projection (Higham 2002 /
// Dykstra). Alternates projections between the PSD cone and the unit-diagonal
// affine space, followed by shrinkage (1 - lambda) * R + lambda * I to
// guarantee strict PD.
arma::mat nearest_correlation(const arma::mat &raw, double shrinkage,
                              double tolerance = 1e-8,
                              int max_iterations = 100) {
  arma::mat y = 0.5 * (raw + raw.t());
  y.diag().ones();

  // Fast path: if the matrix is already PSD (or within numerical tolerance),
  // apply shrinkage directly and skip iterative projection.
  arma::vec raw_eigenvalues;
  if (!arma::eig_sym(raw_eigenvalues, y)) {
    Rcpp::stop(
        "Eigenvalue decomposition failed during correlation projection.");
  }
  if (raw_eigenvalues.min() >= -tolerance) {
    arma::mat correlation =
        (1.0 - shrinkage) * y + shrinkage * arma::eye(y.n_rows, y.n_cols);
    correlation.diag().ones();
    return 0.5 * (correlation + correlation.t());
  }

  // Dykstra's alternating projection algorithm
  arma::mat dykstra(y.n_rows, y.n_cols, arma::fill::zeros);

  for (int iteration = 0; iteration < max_iterations; ++iteration) {
    const arma::mat previous = y;

    // Project onto PSD cone and update Dykstra correction
    const arma::mat residual = y - dykstra;
    const arma::mat psd = project_psd(residual);
    dykstra = psd - residual;

    // Project onto unit-diagonal space
    y = psd;
    y.diag().ones();

    // Check convergence in relative Frobenius norm
    const double denominator = std::max(1.0, arma::norm(previous, "fro"));
    if (arma::norm(y - previous, "fro") / denominator <= tolerance) {
      break;
    }
  }

  // Final PSD projection and diagonal rescaling: D^{-1/2} R D^{-1/2}
  arma::mat correlation = project_psd(y);
  const arma::vec diagonal = correlation.diag();
  if (arma::any(diagonal <= std::numeric_limits<double>::epsilon())) {
    Rcpp::stop("Correlation projection produced a non-positive diagonal.");
  }
  const arma::vec scale = 1.0 / arma::sqrt(diagonal);
  correlation = arma::diagmat(scale) * correlation * arma::diagmat(scale);
  correlation = 0.5 * (correlation + correlation.t());
  correlation.diag().ones();

  // Shrink toward identity: guarantees all eigenvalues >= shrinkage
  correlation = (1.0 - shrinkage) * correlation +
                shrinkage * arma::eye(correlation.n_rows, correlation.n_cols);
  correlation.diag().ones();
  return 0.5 * (correlation + correlation.t());
}

// Estimate pairwise latent correlations from pre-sorted columns and project to
// correlation cone. Uses Gaussian copula sine bridge: rho_{ij} = sin(pi * tau_a
// / 2).
arma::mat correlation_from_sorted(const arma::mat &data,
                                  const std::vector<SortedColumn> &columns,
                                  double shrinkage) {
  const arma::uword d = data.n_cols;
  arma::mat correlation(d, d, arma::fill::eye);

  for (arma::uword i = 1; i < d; ++i) {
    for (arma::uword j = 0; j < i; ++j) {
      const double tau =
          kendall_tau_a_from_sorted(columns[i], columns[j], data.col(j));
      correlation(i, j) = std::sin(arma::datum::pi * tau / 2.0);
    }
  }

  // Symmetrize and project to the nearest strictly positive-definite
  // correlation matrix
  correlation = arma::symmatl(correlation);
  return nearest_correlation(correlation, shrinkage);
}

// Symmetric matrix square root R^{1/2} using Armadillo's sqrtmat_sympd
arma::mat symmetric_sqrt(const arma::mat &matrix) {
  arma::mat result;
  const arma::mat symmetric = 0.5 * (matrix + matrix.t());
  if (!arma::sqrtmat_sympd(result, symmetric)) {
    Rcpp::stop("Unable to compute a correlation-matrix square root.");
  }
  return 0.5 * (result + result.t());
}

} // namespace

// Standalone exported entry point for Kendall's tau-a (useful for unit testing)
// [[Rcpp::export]]
double kendall_tau_a_cpp(const arma::vec &x, const arma::vec &y) {
  if (x.n_elem != y.n_elem || x.n_elem < 2) {
    Rcpp::stop("`x` and `y` must have the same length of at least two.");
  }
  if (!x.is_finite() || !y.is_finite()) {
    Rcpp::stop("`x` and `y` must contain only finite values.");
  }
  const SortedColumn sorted_x = sort_column(x);
  const SortedColumn sorted_y = sort_column(y);
  return kendall_tau_a_from_sorted(sorted_x, sorted_y, y);
}

// Main C++ representation precomputation kernel called by as_nonparanormal()
// [[Rcpp::export]]
Rcpp::List distribution_summaries_cpp(Rcpp::List data, int M,
                                      double pd_shrinkage, bool cache_sqrt) {
  const int number_distributions = data.size();
  if (number_distributions < 1 || M < 2) {
    Rcpp::stop("Compiled summary input is invalid.");
  }

  // Unpack R list into C++ Armadillo matrices
  std::vector<arma::mat> distributions(number_distributions);
  for (int i = 0; i < number_distributions; ++i) {
    distributions[i] = Rcpp::as<arma::mat>(data[i]);
  }

  const arma::uword d = distributions[0].n_cols;

  // Allocate storage for outputs:
  // - Quantiles: list of d matrices, each of size (n x M)
  // - Correlations: list of n matrices, each of size (d x d)
  // - Correlation square roots: optional list of n matrices (d x d)
  std::vector<arma::mat> quantile_storage(
      d, arma::mat(number_distributions, M, arma::fill::none));
  Rcpp::List correlations(number_distributions);
  Rcpp::List correlation_sqrts(cache_sqrt ? number_distributions : 0);
  Rcpp::IntegerVector sample_sizes(number_distributions);
  arma::mat tied_pair_fractions(number_distributions, d, arma::fill::zeros);

  // Compute summaries for each empirical distribution
  for (int distribution_index = 0; distribution_index < number_distributions;
       ++distribution_index) {
    const arma::mat &current = distributions[distribution_index];
    const arma::uword sample_size = current.n_rows;
    const count_t total_pairs = choose_two(sample_size);
    sample_sizes[distribution_index] = static_cast<int>(sample_size);

    // 1. Sort each variable once and compute marginal quantiles & tie
    // diagnostics
    std::vector<SortedColumn> columns;
    columns.reserve(d);
    for (arma::uword variable = 0; variable < d; ++variable) {
      columns.push_back(sort_column(current.col(variable)));
      quantile_storage[variable].row(distribution_index) =
          quantiles_type7(columns.back().values, M).t();
      tied_pair_fractions(distribution_index, variable) =
          static_cast<double>(columns.back().tied_pairs) /
          static_cast<double>(total_pairs);
    }

    // 2. Compute pairwise Kendall's tau-a, sine-bridge, and nearest correlation
    // matrix
    const arma::mat correlation =
        correlation_from_sorted(current, columns, pd_shrinkage);
    correlations[distribution_index] = correlation;

    // 3. (Optional) Precompute symmetric square root R^{1/2} for downstream
    // distance calculations
    if (cache_sqrt) {
      correlation_sqrts[distribution_index] = symmetric_sqrt(correlation);
    }
  }

  // Convert quantile matrices into an R list of length d
  Rcpp::List quantiles(d);
  for (arma::uword variable = 0; variable < d; ++variable) {
    quantiles[variable] = quantile_storage[variable];
  }

  // Construct probability grid vector
  Rcpp::NumericVector probabilities(M);
  for (int k = 0; k < M; ++k) {
    probabilities[k] = midpoint_probability(k, M);
  }

  return Rcpp::List::create(
      Rcpp::Named("probabilities") = probabilities,
      Rcpp::Named("quantiles") = quantiles,
      Rcpp::Named("correlations") = correlations,
      Rcpp::Named("correlation_sqrts") =
          cache_sqrt ? static_cast<SEXP>(correlation_sqrts) : R_NilValue,
      Rcpp::Named("sample_sizes") = sample_sizes,
      Rcpp::Named("tied_pair_fractions") = tied_pair_fractions);
}
