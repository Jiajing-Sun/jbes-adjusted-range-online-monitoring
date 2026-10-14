# Check transition indexing independently of random number generation.
# Mock draws return conditional means, making the first changed observation
# identifiable exactly. The production generator itself is not modified.
root <- normalizePath(getwd())
model <- new.env(parent = globalenv())
sys.source(file.path(root, '20_dgps_par1.R'), model)
model$rpois <- function(n, lambda) rep(lambda, n)
model$rbinom <- function(n, size, prob) rep(size * prob, n)
for (m in c(4L, 10L)) for (k_star in c(1L, m - 1L)) for (type in c('abrupt', 'smooth')) {
  observations <- model$sim_par1_break(m, 1, k_star, 0, 1, 0, 9, type, burnin = 0)
  first_changed <- m + k_star + 2L
  stopifnot(
    length(observations) == 2 * m + 1L,
    all(observations[seq_len(first_changed - 1L)] == 1)
  )
  expected <- if (type == 'abrupt') 9 else 1 + 8 * (1 - exp(-.05))
  stopifnot(abs(observations[first_changed] - expected) < 1e-12)
}
# A score at row j uses the pair (X_{j-1}, X_j). This check verifies
# the m training contributions and the first monitored pair explicitly.
model$par1_score_transition <- function(x, y, alpha, lambda) c(y, x)
observations <- seq.int(0L, 20L)
scores <- model$par1_score_matrix(observations, .2, 2)
stopifnot(nrow(scores) == 20L, identical(unname(scores[11L, ]), c(10, 11)))
cat('PAR first changed transition, smooth timing and score-pair indexing passed.\n')
