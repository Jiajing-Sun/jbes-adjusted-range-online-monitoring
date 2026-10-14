# ============================================================
# 20_dgps_par1.R  -- PAR(1) DGP + MLE + score for monitoring
# ============================================================

# Simulate a PAR(1): X_t = Binomial(X_{t-1}, alpha) + Poisson(lambda)
sim_par1 <- function(n, alpha, lambda, burnin = 200, x0 = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(x0)) {
    mu <- lambda / (1 - alpha)
    x0 <- rpois(1, mu)
  }
  Ntot <- n + burnin
  X <- integer(Ntot)
  X[1] <- x0
  for (t in 2:Ntot) {
    thin <- rbinom(1, X[t-1], alpha)
    innov <- rpois(1, lambda)
    X[t] <- thin + innov
  }
  X_use <- X[(burnin+1):Ntot]
  stopifnot(length(X_use) == n)
  X_use
}

# Time-varying PAR(1) with a single break at t_star (observation index)
sim_par1_break <- function(m, T, k_star,
                           alpha1, lambda1, alpha2, lambda2,
                           type = c("abrupt", "smooth"),
                           delta_smooth = 0.05,
                           burnin = 400,
                           seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  type <- match.arg(type)

  mT <- as.integer(round(m * T))
  n_obs <- m + mT + 1L # +1 so we have m+mT transitions

  # break location in observation time:
  t_star <- m + k_star

  # initial from stationary approx of pre-change regime
  mu <- lambda1 / (1 - alpha1)
  X <- integer(n_obs + burnin)
  X[1] <- rpois(1, mu)

  for (t in 2:(n_obs + burnin)) {
    # Map 't' to post-burnin observation index
    t_eff <- t - burnin
    # Use pre-change params during burn-in
    if (t_eff <= 0) {
      a <- alpha1; lam <- lambda1
    } else {
      if (t_eff <= t_star) {
        a <- alpha1; lam <- lambda1
      } else {
        if (type == "abrupt") {
          a <- alpha2; lam <- lambda2
        } else {
          # smooth exponential transition
          s <- t_eff - t_star
          a <- alpha1 + (alpha2 - alpha1) * (1 - exp(-delta_smooth * s))
          lam <- lambda1 + (lambda2 - lambda1) * (1 - exp(-delta_smooth * s))
        }
      }
    }

    thin <- rbinom(1, X[t-1], a)
    innov <- rpois(1, lam)
    X[t] <- thin + innov
  }

  X_use <- X[(burnin+1):(burnin + n_obs)]
  stopifnot(length(X_use) == n_obs)
  X_use
}

# Conditional probability P(X_t = x | X_{t-1} = y) under PAR(1)
par1_cond_prob <- function(x, y, alpha, lambda) {
  jmax <- min(x, y)
  j <- 0:jmax
  sum(dbinom(j, y, alpha) * dpois(x - j, lambda))
}

# Negative log-likelihood (conditional on X1), using transformed params:
# par[1] = logit(alpha), par[2] = log(lambda)
par1_negloglik_trans <- function(par, X) {
  alpha <- exp(par[1]) / (1 + exp(par[1]))
  lambda <- exp(par[2])

  n <- length(X)
  ll <- 0
  for (t in 2:n) {
    p <- par1_cond_prob(X[t], X[t-1], alpha, lambda)
    p <- max(p, 1e-300)
    ll <- ll + log(p)
  }
  -ll
}

# Fit PAR(1) by MLE on training observations X_train (length >= 3)
par1_fit_mle <- function(X_train) {
  n <- length(X_train)
  stopifnot(n >= 3)

  # crude initial values
  acf1 <- suppressWarnings(cor(X_train[2:n], X_train[1:(n-1)]))
  if (!is.finite(acf1)) acf1 <- 0.2
  alpha0 <- min(max(acf1, 0.01), 0.95)

  mu_hat <- mean(X_train)
  lambda0 <- max(mu_hat * (1 - alpha0), 0.05)

  par0 <- c(qlogis(alpha0), log(lambda0))
  opt <- optim(par0, par1_negloglik_trans, X = X_train,
               method = "BFGS",
               control = list(maxit = 300))
  alpha_hat <- plogis(opt$par[1])
  lambda_hat <- exp(opt$par[2])
  list(alpha = alpha_hat, lambda = lambda_hat, opt = opt)
}

# Exact score for one transition (y -> x)
par1_score_transition <- function(x, y, alpha, lambda) {
  jmax <- min(x, y)
  j <- 0:jmax

  w <- dbinom(j, y, alpha) * dpois(x - j, lambda)
  P <- sum(w)
  P <- max(P, 1e-300)

  dP_da <- sum(w * (j / alpha - (y - j) / (1 - alpha)))
  dP_dl <- sum(w * ((x - j) / lambda - 1))

  c(dP_da / P, dP_dl / P)
}

# Score matrix for a full series X (length n_obs):
# returns (n_obs - 1) x 2 matrix of scores for transitions t=1..n_obs-1.
par1_score_matrix <- function(X, alpha, lambda) {
  n <- length(X)
  S <- matrix(0, n - 1L, 2L)
  for (t in 1:(n-1L)) {
    S[t, ] <- par1_score_transition(X[t+1], X[t], alpha, lambda)
  }
  colnames(S) <- c("score_alpha", "score_lambda")
  S
}
