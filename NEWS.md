# NonparaTransport 1.0.1

## Behavior changes

* `as_nonparanormal(standardize = TRUE)` now uses R's `mad()` scaling
  (1.4826 times raw MAD). Compared with previous releases, standardized marginal
  squared distances are divided by 1.4826^2. Latent correlations are unchanged.

## Performance improvements

* `as_nonparanormal()` avoids redundant matrix decompositions and uses a looser
  stopping tolerance for iterative correlation projection. Results may differ
  slightly when projection is needed; the requested identity shrinkage is retained.
* `pairwise_npt_distance()` and `npt_distance()` avoid redundant matrix
  allocations when assembling distance components and applying `bw_weight`.
* `npt_frechetreg()` speeds up bivariate correlation fitting at multiple
  prediction points by reusing scalar roots of observed correlations.
* Bivariate summaries now return `correlation_sqrts = NULL` and
  `cache_sqrt = FALSE`, even when caching is requested, because bivariate
  distances do not use matrix square roots.
* Minor efficiency improvements in latent correlation preprocessing.

# NonparaTransport 1.0.0

* First release.
