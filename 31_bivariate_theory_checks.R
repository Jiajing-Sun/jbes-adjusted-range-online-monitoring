# Bivariate full-rank training profile: interpretation and output guide.
#
# theta indexes the interior vertex of the piecewise-linear training path.
# Training means are m times its discrete increments, not ordinary iid means.
# The last training regime (-4,0) persists into monitoring. The additional
# local shift is (-eta,delta)/sqrt(m): eta controls the predicted favorable or
# adverse direction; delta is transverse in the first-order limit.
# rho is the correlation in a positive-definite covariance matrix with
# diagonal entries 1 and 1.5. Both methods receive identical oracle DGP
# information and independent calibration/null/alternative sample blocks.
#
# This design checks the specified full-rank proposition, not the different
# rank-one single-training-break law. Saved stopping times retain all signs
# and transverse shifts. Paired errors condition on the estimated thresholds.

# Full-rank bivariate profile from the manuscript's directional proposition.
# JBES_RESUME_VERIFIED=1 is only for checkpoints with externally verified identical designs, seeds and thresholds.
# All signs and transverse directions are retained. This is not a rank-one model.
root <- normalizePath(Sys.getenv('JBES_PROJECT_ROOT', getwd()))
for (f in c('00_utils.R', '02_methods.R', '08_mc_audit.R')) source(file.path(root, f))
jbes_set_root(root);out <- file.path(root, Sys.getenv('JBES_BIVARIATE_OUT', 'results/bivariate_theory'));dir.create(out, recursive = TRUE, showWarnings = FALSE)
N_CAL <- env_int('JBES_THEORY_CAL', 10000L);N_EVAL <- env_int('JBES_THEORY_EVAL', 5000L)
CORES <- env_int('JBES_N_CORES', 4L);MS <- env_num_vec('JBES_THEORY_M', c(500, 2000))
THETAS <- c(.25, .75);RHOS <- c(0, .6);T0 <- 2;TAU <- .4;ALPHA <- .05
SCENES <- expand.grid(eta = c(-4, 0, 2, 4), delta = c(-2, 0, 2));METHODS <- c('RSMS', 'SSMS', 'CSMS')
bivariate_record <- function(i, phase, m, theta, rho, critical = NULL) {
  config <- (match(theta, THETAS) - 1) * length(RHOS) * length(MS) + (match(rho, RHOS) - 1) * length(MS) + match(m, MS)
  set.seed(20261015L + config * 100000L + phase * 20000L + i)
  N <- m * T0;k <- seq_len(N);ks <- floor(m * TAU);r <- (0:m) / m
  f <- approx(c(0, .25, .5, .75, 1), c(0, 1, theta, 1, 0), xout = r)$y
  h <- approx(c(0, .125, .25, .5, .75, .875, 1), c(0, 0, 1, 0, -1, 0, 0), xout = r)$y
  g <- cbind(f, h);drift <- m * (g[-1, , drop = FALSE] - g[-nrow(g), , drop = FALSE])
  Sigma <- matrix(c(1, rho * sqrt(1.5), rho * sqrt(1.5), 1.5), 2)
  eps <- matrix(rnorm((m + N) * 2), m + N, 2) %*% chol(Sigma)
  x <- eps + rbind(drift, cbind(rep(-4, N), rep(0, N)))
  baseline <- colMeans(x[1:m, , drop = FALSE]);train <- sweep(x[1:m, , drop = FALSE], 2, baseline)
  fit <- aersn_training_fit(train);gauge <- aersn_gauge_evaluator(fit$hull)
  Sinv <- inv_pd(aersn::aersn_shao_normalizer(fit, integration = 'calendar')$matrix, 1e-8)
  Cinv <- inv_pd(hac_lrv(train), 1e-8)
  Z0 <- apply(sweep(x[m + k, , drop = FALSE], 2, baseline), 2, cumsum) / sqrt(m);w <- (1 + k / m)^2
  path <- function(eta, delta) {
    Z <- Z0 + (pmax(k - ks, 0) / m) %o% c(-eta, delta)
    cbind(RSMS = gauge(Z)^2, SSMS = rowSums((Z %*% Sinv) * Z), CSMS = rowSums((Z %*% Cinv) * Z)) / w
  }
  if (phase == 1) return(apply(path(0, 0), 2, max))
  scenes <- if (phase == 2) data.frame(eta = 0, delta = 0) else SCENES
  stops <- matrix(N + 1, nrow(scenes), 3)
  for (s in seq_len(nrow(scenes))) {
    p <- path(scenes$eta[s], scenes$delta[s]);for (j in 1:3) stops[s, j] <- first_exceed(p[, j], critical[j])
  }
  stops
}
summary_rows <- pair_rows <- cv_rows <- curve_rows <- list()
for (theta in THETAS) for (rho in RHOS) for (m in MS) {
  tag <- sprintf('theta%.2f_rho%.2f_m%d', theta, rho, m);cat('Bivariate theory:', tag, '\n');flush.console()
  cache_cal <- file.path(out, paste0(tag, '_calibration.rds'))
  if (identical(Sys.getenv('JBES_RESUME_VERIFIED'), '1') && file.exists(cache_cal)) {
    maxima <- readRDS(cache_cal);stopifnot(nrow(maxima) == N_CAL, ncol(maxima) == 3, all(is.finite(maxima)))
  } else {
    cal <- par_lapply_safe(seq_len(N_CAL), function(i) bivariate_record(i, 1, m, theta, rho), CORES)
    maxima <- do.call(rbind, cal)
  }
  critical <- apply(maxima, 2, quantile, 1 - ALPHA, type = 8)
  saveRDS(maxima, file.path(out, paste0(tag, '_calibration.rds')))
  prob <- pmin(1, pmax(0, (1 - ALPHA) + c(-1, 1) * qnorm(.975) * sqrt(ALPHA * (1 - ALPHA) / N_CAL)))
  for (j in 1:3) cv_rows[[length(cv_rows) + 1]] <- data.frame(theta, rho, m, method = METHODS[j], critical = critical[j], lo = quantile(maxima[, j], prob[1]), hi = quantile(maxima[, j], prob[2]), n_cal = N_CAL)
  for (phase in 2:3) {
    cache_eval <- file.path(out, paste0(tag, '_phase', phase, '.rds'))
    if (identical(Sys.getenv('JBES_RESUME_VERIFIED'), '1') && file.exists(cache_eval)) {
      ans <- readRDS(cache_eval);stopifnot(length(ans) == N_EVAL, all(vapply(ans, function(x) all(is.finite(x)), FALSE)))
    } else ans <- par_lapply_safe(seq_len(N_EVAL), function(i) bivariate_record(i, phase, m, theta, rho, critical), CORES)
    saveRDS(ans, file.path(out, paste0(tag, '_phase', phase, '.rds')))
    scenes <- if (phase == 2) data.frame(eta = 0, delta = 0) else SCENES
    vs <- T0 + T0^2 + 16 * T0^2 * 3 * (23 * theta^2 + 59 * theta + 46) / (160 * (theta^2 + theta + 2)^2);vr <- T0 + 3 * T0^2
    for (s in seq_len(nrow(scenes))) {
      st <- t(vapply(ans, function(z) z[s, ], numeric(3)));colnames(st) <- METHODS;a <- mc_event_audit(st, floor(m * TAU), m * T0)
      meta <- data.frame(theta, rho, m, phase, eta = scenes$eta[s], delta = scenes$delta[s])
      row <- cbind(meta[rep(1, 3), , drop = FALSE], a$events);row$theory_power <- c(pnorm(qnorm(1 - ALPHA) - scenes$eta[s] * (T0 - TAU) / sqrt(c(vr, vs)), lower.tail = FALSE), NA_real_)
      summary_rows[[length(summary_rows) + 1]] <- row
      pair_rows[[length(pair_rows) + 1]] <- cbind(meta[rep(1, nrow(a$paired)), , drop = FALSE], a$paired)
      for (u in c(.1, .2, .4, .6, .8, .9, 1) * T0) for (j in 1:3) curve_rows[[length(curve_rows) + 1]] <- data.frame(meta, method = METHODS[j], u, first_crossing_probability = mean(st[, j] <= floor(m * u)), n_eval = N_EVAL)
    }
  }
  write.csv(do.call(rbind, summary_rows), file.path(out, 'bivariate_summary.csv'), row.names = FALSE)
  write.csv(do.call(rbind, pair_rows), file.path(out, 'bivariate_paired.csv'), row.names = FALSE)
  write.csv(do.call(rbind, cv_rows), file.path(out, 'bivariate_calibration.csv'), row.names = FALSE)
  write.csv(do.call(rbind, curve_rows), file.path(out, 'bivariate_first_crossing.csv'), row.names = FALSE)
}
writeLines(capture.output(sessionInfo()), file.path(out, 'session.txt'))
