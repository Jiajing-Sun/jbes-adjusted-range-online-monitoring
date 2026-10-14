# ============================================================
# 10_dgps_mean.R  -- DGPs for mean-monitoring simulations
# ============================================================

# Helper: Toeplitz covariance with parameter rho
toeplitz_cov <- function(q, rho) {
  mat <- outer(1:q, 1:q, function(i, j) rho^(abs(i - j)))
  mat
}

# Helper: banded adjacency matrix (1 on immediate off-diagonals)
band_adj <- function(q) {
  A <- matrix(0, q, q)
  if (q >= 2) {
    for (i in 1:(q - 1)) {
      A[i, i + 1] <- 1
      A[i + 1, i] <- 1
    }
  }
  A
}

# Default VAR(1) coefficient used in multiple DGPs:
# Psi = phi * I_q + cross * A, where A is band adjacency (stable for small cross)
make_Psi <- function(q, phi = 0.5, cross = 0.1) {
  phi * diag(q) + cross * band_adj(q)
}

# ------------------------------------------------------------
# Simulate VAR(1): X_t = Psi X_{t-1} + eps_t
# dgp in {"DGP1","DGP2","DGP3","DGP4"}
# - DGP1: Psi = 0.5 I, eps ~ N(0, I)
# - DGP2: Psi = 0.5 I + 0.1 A, eps ~ N(0, Toeplitz(rho=0.1))
# - DGP3: Psi as DGP2, eps has componentwise GARCH(1,1) with e_t ~ N(0,I)
# - DGP4: Psi as DGP2, eps variance shift at t > floor(m_train/2)
# ------------------------------------------------------------
sim_var1_dgp <- function(dgp, n, m_train, q = 2,
                         burnin = 200,
                         phi = 0.5, cross = 0.1, rho = 0.1,
                         garch_alpha = 0.1, garch_beta = 0.2,
                         vol_shift_delta = 0.2,
                         seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  dgp <- toupper(dgp)
  Psi <- switch(dgp,
    "DGP1" = phi * diag(q),
    "DGP2" = make_Psi(q, phi = phi, cross = cross),
    "DGP3" = make_Psi(q, phi = phi, cross = cross),
    "DGP4" = make_Psi(q, phi = phi, cross = cross),
    stop("Unknown dgp: ", dgp)
  )

  # error covariance baseline
  Sigma_eps <- switch(dgp,
    "DGP1" = diag(q),
    "DGP2" = toeplitz_cov(q, rho),
    "DGP3" = diag(q),
    "DGP4" = toeplitz_cov(q, rho)
  )

  # R returns an upper triangular factor U with t(U) %*% U = Sigma_eps.
  # Innovations are column vectors, so use the lower factor t(U).
  # Using U directly would instead generate covariance U %*% t(U).
  chol_Sigma <- t(chol(Sigma_eps))

  Ntot <- n + burnin
  X <- matrix(0, Ntot, q)

  if (dgp %in% c("DGP1", "DGP2", "DGP4")) {
    # Draw in the same time-major order as repeated rnorm(q) calls. The first
    # state is zero, so it consumes no innovation. Row-vector multiplication
    # by t(chol_Sigma) is equivalent to chol_Sigma %*% z for each column draw.
    innovations <- matrix(rnorm((Ntot - 1L) * q), Ntot - 1L, q, byrow = TRUE) %*% t(chol_Sigma)
    if (dgp == "DGP4") {
      effective_time <- (2:Ntot) - burnin
      innovations <- innovations * ifelse(effective_time > floor(m_train / 2),
        sqrt(1 + vol_shift_delta), 1
      )
    }
    # Psi is symmetric in all listed designs. Rotate to its orthogonal
    # eigenbasis and apply scalar AR recursions in compiled stats::filter.
    # This preserves the model and random draws while avoiding a slow R loop.
    # A direct-recursion comparison checks numerical equivalence separately.
    decomposition <- eigen(Psi, symmetric = TRUE)
    rotated <- rbind(rep(0, q), innovations %*% decomposition$vectors)
    for (j in seq_len(q)) rotated[, j] <- as.numeric(stats::filter(
      rotated[, j],
      filter = decomposition$values[j], method = "recursive"
    ))
    X <- rotated %*% t(decomposition$vectors)
  }

  if (dgp == "DGP3") {
    # componentwise GARCH(1,1) on eps (diagonal Sigma_t)
    sig2 <- matrix(1, Ntot, q)
    eps_mat <- matrix(0, Ntot, q)
    c0 <- 1 - garch_alpha - garch_beta
    for (t in 2:Ntot) {
      z <- rnorm(q) # independent innovations
      sig2[t, ] <- c0 + garch_alpha * (eps_mat[t - 1, ]^2) + garch_beta * sig2[t - 1, ]
      eps <- sqrt(sig2[t, ]) * z
      eps_mat[t, ] <- eps
      X[t, ] <- as.numeric(Psi %*% X[t - 1, ]) + eps
    }
  }


  # drop burnin
  X_use <- X[(burnin + 1):Ntot, , drop = FALSE]
  stopifnot(nrow(X_use) == n)
  X_use
}
