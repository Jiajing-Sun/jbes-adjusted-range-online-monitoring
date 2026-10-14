# Exact finite-grid maximum and first crossing, with certified LP screening.
# These functions preserve every original monitoring time and strict > boundary.
rsms_grid_objects <- function(psi, m, T, gamma) {
  psi <- as.matrix(psi);q <- ncol(psi);N <- as.integer(round(m * T))
  stopifnot(m > q, N >= 1, nrow(psi) >= m + N, gamma >= 0, gamma < .5)
  train <- psi[seq_len(m), , drop = FALSE];mon <- psi[m + seq_len(N), , drop = FALSE]
  G <- aersn::aersn_path(train)$G
  Z <- matrix(apply(mon, 2, cumsum), ncol = q) / sqrt(m)
  k <- seq_len(N);w <- 1 / ((1 + k / m) * (k / (m + k))^gamma)
  list(G = G, Z = Z, w = w, N = N)
}
rsms_maximum <- function(psi, m, T, gamma) {
  x <- rsms_grid_objects(psi, m, T, gamma)
  as.numeric(ae_gauge_maxima(x$G, x$Z, x$N, matrix(x$w))^2)
}
ae_gauge_first_crossing <- function(G, Z, weights, boundary) {
  G <- unique(as.matrix(G));Z <- as.matrix(Z);q <- ncol(G);N <- nrow(Z);n <- nrow(G)
  stopifnot(length(weights) == N, all(weights > 0), boundary > 0)
  if (q <= 2L) {
    hit <- which(ae_gauge_values(G, Z) * weights > boundary)
    return(if (length(hit)) hit[1L] else Inf)
  }
  sv <- svd(G, nu = 0, nv = q)
  if (length(sv$d) < q || min(sv$d) <= max(sv$d) * 1e-12) stop('Rank deficient training hull')
  H <- sweep(sv$v, 2, sv$d, '/');G <- G %*% H;Z <- Z %*% H
  gauge_value <- ae_lp_engine(G)
  V <- crossprod(scale(G, scale = FALSE)) / n
  E <- rbind(rep(0, q), Z %*% solve(chol(V)) / 2 * weights)
  cache <- c(0, rep(NA_real_, N));count <- 0L
  value <- function(i) {
    if (!is.na(cache[i])) return(cache[i])
    j <- i - 1L
    v <- gauge_value(Z[j, ]) * weights[j]
    count <<- count + 1L;cache[i] <<- v;v
  }
  value(N + 1L)
  visit <- function(l, r) {
    if (cache[l] > boundary) return(l - 1L)
    if (r - l <= 1L) return(if (cache[r] > boundary) r - 1L else Inf)
    ids <- seq.int(l, r);EE <- E[ids, , drop = FALSE]
    dl <- sqrt(rowSums((EE - rep(E[l, ], each = nrow(EE)))^2));dr <- sqrt(rowSums((EE - rep(E[r, ], each = nrow(EE)))^2))
    upper <- max(pmin(cache[l] + dl, cache[r] + dr))
    if (upper + 1e-7 * (1 + boundary) < boundary) return(Inf)
    mid <- as.integer((l + r) %/% 2L);value(mid)
    a <- visit(l, mid)
    if (is.finite(a)) return(a)
    visit(mid, r)
  }
  ans <- visit(1L, N + 1L);attr(ans, 'n_lp') <- count;ans
}
rsms_first_crossing <- function(psi, m, T, gamma, cval) {
  x <- rsms_grid_objects(psi, m, T, gamma)
  out <- ae_gauge_first_crossing(x$G, x$Z, x$w, sqrt(cval))
  if (is.finite(out)) out else x$N + 1L
}
