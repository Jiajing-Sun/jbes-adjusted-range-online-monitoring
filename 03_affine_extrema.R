# Exact maxima on the original simulation grid. No grid points are removed.
# Popoviciu's inequality implies K(G) contains the ellipsoid
# 2 * V^(1/2) * unit_ball, where V is the covariance of the grid points.
# Therefore p_K(x-y) <= ||(x-y) V^(-1/2)||/2. Triangle inequalities
# certify intervals that cannot exceed the largest already evaluated gauge.
ae_gauge_maxima <- function(G, Z, horizons, weights) {
  G <- unique(as.matrix(G));Z <- as.matrix(Z);weights <- as.matrix(weights)
  q <- ncol(G);n <- nrow(G);N <- nrow(Z)
  stopifnot(ncol(Z) == q, nrow(weights) == N, all(weights > 0), all(horizons >= 1), all(horizons <= N))
  out <- matrix(NA_real_, length(horizons), ncol(weights))
  if (q <= 2L) {
    vals <- ae_gauge_values(G, Z)
    for (h in seq_along(horizons)) for (w in seq_len(ncol(weights)))
      out[h, w] <- max(vals[seq_len(horizons[h])] * weights[seq_len(horizons[h]), w])
    return(out)
  }
  sv <- svd(G, nu = 0, nv = q)
  if (length(sv$d) < q || min(sv$d) <= max(sv$d) * 1e-12) stop('Rank deficient training hull')
  H <- sweep(sv$v, 2, sv$d, '/');G <- G %*% H;Z <- Z %*% H
  gauge_value <- ae_lp_engine(G)
  V <- crossprod(scale(G, scale = FALSE)) / n
  Eraw <- Z %*% solve(chol(V)) / 2
  cache <- rep(NA_real_, N);count <- 0L
  evalone <- function(i) {
    if (!is.na(cache[i])) return(cache[i])
    v <- gauge_value(Z[i, ])
    cache[i] <<- max(v, 0);count <<- count + 1L;cache[i]
  }
  for (h in seq_along(horizons)) for (w in seq_len(ncol(weights))) {
    nH <- as.integer(horizons[h]); ww <- weights[seq_len(nH), w];E <- Eraw[seq_len(nH), , drop = FALSE] * ww
    best <- max(evalone(1L) * ww[1L], evalone(nH) * ww[nH])
    jseed <- which.max(rowSums(E^2));best <- max(best, evalone(jseed) * ww[jseed])
    known <- which(!is.na(cache[seq_len(nH)]));best <- max(best, cache[known] * ww[known])
    visit <- function(l, r) {
      if (r - l <= 1L) return()
      ids <- seq.int(l, r);EE <- E[ids, , drop = FALSE]
      dl <- sqrt(rowSums((EE - rep(E[l, ], each = nrow(EE)))^2))
      dr <- sqrt(rowSums((EE - rep(E[r, ], each = nrow(EE)))^2))
      upper <- max(pmin(cache[l] * ww[l] + dl, cache[r] * ww[r] + dr))
      # Conservative numerical cushion: never screen close to the incumbent.
      if (upper + 1e-7 * (1 + abs(best)) < best) return()
      j <- as.integer((l + r) %/% 2L);v <- evalone(j) * ww[j];best <<- max(best, v)
      visit(l, j);visit(j, r)
    }
    visit(1L, nH);out[h, w] <- best
  }
  attr(out, 'n_lp') <- count
  out
}
