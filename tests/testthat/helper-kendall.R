# Transparent O(N_i^2) reference used only to verify the compiled O(N_i log N_i)
# implementation. Ties contribute zero through base R's sign(0) behavior.
brute_tau_a <- function(x, y) {
  sample_size <- length(x)
  score <- 0
  for (i in seq_len(sample_size - 1L)) {
    for (j in seq.int(i + 1L, sample_size)) {
      score <- score + sign(x[i] - x[j]) * sign(y[i] - y[j])
    }
  }
  score / choose(sample_size, 2)
}

example_distributions <- function() {
  list(
    first = cbind(
      depth = c(1, 1, 2, 3, 3, 4),
      duration = c(3, 1, 4, 2, 5, 6)
    ),
    second = cbind(
      depth = c(2, 2, 2, 3, 4, 5, 5),
      duration = c(1, 3, 2, 5, 4, 7, 6)
    )
  )
}
