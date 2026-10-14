# Exact containment diagnostic for the theoretical (unregularized) SSMS.
# G includes the time-zero bridge value. D uses the original m training times.
training_containment_kappa <- function(G, D) {
  G <- as.matrix(G); q <- ncol(G)
  H <- chol(D)
  W <- G %*% solve(H)
  if (q == 1L) return(diff(range(W[, 1L]))^2)
  if (q == 2L) W <- W[chull(W[, 1L], W[, 2L]), , drop = FALSE]
  # Maximum of a convex quadratic over the increment hull occurs at an increment.
  max(as.matrix(dist(W))^2)
}

# Reuse a two-dimensional training hull across all monitoring scenarios.
# The polygon formula is the same exact finite-grid gauge used by ae_gauge_values.
persistent_training_objects <- function(X, m) {
  stopifnot(ncol(X) == 2L)
  psi <- sweep(X, 2L, colMeans(X[seq_len(m), , drop = FALSE]))
  fit <- aersn_training_fit(psi[seq_len(m), , drop = FALSE])
  G <- fit$path$G
  D <- aersn::aersn_shao_normalizer(fit, integration = 'calendar')$matrix
  gauge <- aersn_gauge_evaluator(fit$hull)
  list(
    G = G, D = D, psi = psi, gauge = gauge,
    Sinv = inv_pd(D, ridge = 1e-8),
    Cinv = inv_pd(hac_lrv(psi[seq_len(m), , drop = FALSE]), ridge = 1e-8),
    kappa = training_containment_kappa(G, D)
  )
}
persistent_base_statistics <- function(obj, Z) {
  cbind(
    RSMS = obj$gauge(Z)^2,
    SSMS = rowSums((Z %*% obj$Sinv) * Z),
    CSMS = rowSums((Z %*% obj$Cinv) * Z)
  )
}
