make_regression_example <- function(d = 3L, n = 7L, sample_size = 35L) {
  set.seed(2026 + d)
  predictor <- matrix(
    seq(-1, 1, length.out = n),
    ncol = 1L,
    dimnames = list(paste0("subject_", seq_len(n)), "biomarker")
  )

  distributions <- lapply(seq_len(n), function(i) {
    latent <- matrix(rnorm(sample_size * d), ncol = d)
    latent[, 1L] <- latent[, 1L] + 0.4 * predictor[i, 1L]
    if (d > 1L) {
      latent[, 2L] <- latent[, 2L] +
        0.35 * predictor[i, 1L] * latent[, 1L]
    }
    colnames(latent) <- paste0("response_", seq_len(d))
    latent
  })
  names(distributions) <- rownames(predictor)

  list(
    X = predictor,
    Y = as_nonparanormal(distributions, M = 25L)
  )
}
