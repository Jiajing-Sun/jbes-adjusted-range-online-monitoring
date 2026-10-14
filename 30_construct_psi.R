# ============================================================
# 30_construct_psi.R  -- helper constructors for simulation inputs
# ============================================================

# Apply an abrupt or smooth mean change after t_star = m + k_star (observation index)
apply_mean_change <- function(X, m, T, k_star, Delta, shift_vec,
                              type = c("abrupt", "smooth")) {
  type <- match.arg(type)
  mT <- as.integer(round(m * T))
  n <- m + mT
  stopifnot(nrow(X) >= n)
  Y <- X[1:n, , drop = FALSE]
  t_star <- m + k_star
  stopifnot(t_star < n)

  if (type == "abrupt") {
    idx <- (t_star + 1):n
    Y[idx, ] <- Y[idx, ] + matrix(Delta, nrow = length(idx), ncol = length(shift_vec), byrow = TRUE) *
      matrix(shift_vec, nrow = length(idx), ncol = length(shift_vec), byrow = TRUE)
  } else {
    idx <- (t_star + 1):n
    # linear ramp reaching Delta at horizon end
    w <- (idx - t_star) / (n - t_star)
    Y[idx, ] <- Y[idx, ] + (w * Delta) %*% t(shift_vec)
  }
  Y
}

# Construct psi for mean monitoring: psi_t = Y_t - mu_hat (mu_hat from training window)
make_mean_psi <- function(Y, m, T) {
  mT <- as.integer(round(m * T))
  n <- m + mT
  stopifnot(nrow(Y) >= n)
  mu_hat <- colMeans(Y[1:m, , drop = FALSE])
  psi <- sweep(Y[1:n, , drop = FALSE], 2, mu_hat, FUN = "-")
  list(psi = psi, mu_hat = mu_hat)
}

# Dense vs sparse direction vectors (unit l2 norm)
make_shift_vec <- function(q, type = c("dense", "sparse"), s = 1L) {
  type <- match.arg(type)
  if (type == "dense") {
    v <- rep(1, q)
  } else {
    s <- as.integer(s)
    if (s < 1 || s > q) stop("s must be in {1,...,q}")
    v <- rep(0, q)
    v[1:s] <- 1
  }
  v / sqrt(sum(v^2))
}
