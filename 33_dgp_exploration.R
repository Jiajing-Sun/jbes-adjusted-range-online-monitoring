# Exploratory DGP comparison followed by independent confirmation.
#
# Purpose
#   Identify conditions in which RSMS has higher post-change detection than
#   SSMS and feasible HAC-based CUSUM. The complete screening grid is retained.
#   A favorable screening result is not treated as confirmatory evidence.
#
# Modes
#   JBES_DGP_PHASE=screen  : all listed DGPs, m=250, three baseline methods.
#   JBES_DGP_PHASE=confirm : the top three screening DGPs, m=250 and 500,
#                           independent streams and two automatic HAC controls.
#
# Statistical conventions
#   Every stationary noise process below has long-run covariance I_2.
#   JBES_DGP_DIM=1 uses its first-coordinate marginal (long-run variance 1);
#   the default is dimension 2. Separate output and seed blocks keep the
#   additional scalar screening batch distinct from the original grid.
#   JBES_DGP_SET=profile selects ten separate full-rank training profiles,
#   with iid noise, gamma=0, break fraction 0.2 and unit-length angle shifts.
#   This follow-up batch uses its own seed offset and output directory.
#   Thus an additional shift h*v/sqrt(m), ||v||=1, has the same long-run
#   signal scale across DGPs. Persistent training shifts are identified
#   explicitly and reproduced in their no-additional-change calibration.
#   Calibration, independent null evaluation, and alternative evaluation use
#   separate streams. Methods share each evaluation path. Pre-change alarms
#   are not detections of the additional change and are excluded from CADD.
#   All standard errors condition on the estimated Monte Carlo thresholds.
#
# Outputs
#   dgp_definitions.csv records every model and its parameters.
#   summary.csv and paired.csv retain all settings and ranking reversals.
#   ranking.csv is exploratory; selected_dgps.csv freezes the confirmation set.
#   Raw maxima/stopping records stay in local results, outside the Git repo.

root <- normalizePath(Sys.getenv('JBES_PROJECT_ROOT', getwd()))
for (filename in c('00_utils.R', '01_critical_values.R', '02_methods.R', '08_mc_audit.R')) {
  source(file.path(root, filename))
}
jbes_set_root(root)
phase_name <- Sys.getenv('JBES_DGP_PHASE', 'screen')
stopifnot(phase_name %in% c('screen', 'confirm'))
is_confirmation <- phase_name == 'confirm'
grid_type <- Sys.getenv('JBES_DGP_SET', 'stationary')
stopifnot(grid_type %in% c('stationary', 'profile'))
is_profile <- grid_type == 'profile'
dimension <- env_int('JBES_DGP_DIM', 2L)
stopifnot(dimension %in% c(1L, 2L))
stopifnot(!is_profile || dimension == 2L)
output_suffix <- if (is_profile) '_profile' else if (dimension == 1L) '_q1' else ''
n_cal <- env_int('JBES_DGP_CAL', if (is_confirmation) 10000L else 1000L)
n_eval <- env_int('JBES_DGP_EVAL', if (is_confirmation) 5000L else 1000L)
n_cores <- env_int('JBES_N_CORES', 2L)
training_lengths <- env_num_vec('JBES_DGP_M', if (is_confirmation) c(250, 500) else 250)
horizon <- 2
weight_exponent <- if (is_profile) 0 else 0.15
alpha <- 0.05
burn_in <- 2000L
out_dir <- file.path(root, 'results', paste0('dgp_', phase_name, output_suffix))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# The grid is specified before viewing screening outcomes. Parameter meanings:
# phi: AR coefficient; theta: MA coefficient; lag: seasonal lag;
# arch/garch: conditional-variance recursion; shift/fraction: training change.
specification <- function(id, family, phi = 0, theta = 0, lag = 1, innovation = 'normal',
                          arch = 0, garch = 0, cross = 0, shift = 0, fraction = 0.25,
                          left = 1, right = 1, middle = 0.25) {
  data.frame(id, family, phi, theta, lag, innovation, arch, garch, cross, shift, fraction, left, right, middle)
}
dgp_grid <- do.call(rbind, list(
  specification('iid_normal', 'ar'),
  specification('iid_t5', 'ar', innovation = 't5'),
  specification('iid_mixture', 'ar', innovation = 'mixture'),
  specification('ar_pos05', 'ar', phi = 0.5),
  specification('ar_pos08', 'ar', phi = 0.8),
  specification('ar_pos095', 'ar', phi = 0.95),
  specification('ar_neg05', 'ar', phi = -0.5),
  specification('ar_neg08', 'ar', phi = -0.8),
  specification('ar_pos08_t5', 'ar', phi = 0.8, innovation = 't5'),
  specification('ar_neg08_t5', 'ar', phi = -0.8, innovation = 't5'),
  specification('arma_pos08_neg05', 'arma', phi = 0.8, theta = -0.5),
  specification('arma_pos08_neg09', 'arma', phi = 0.8, theta = -0.9),
  specification('arma_neg05_pos08', 'arma', phi = -0.5, theta = 0.8),
  specification('seasonal_ma5', 'ma', theta = 0.8, lag = 5),
  specification('seasonal_ma12', 'ma', theta = 0.8, lag = 12),
  specification('seasonal_ma20', 'ma', theta = 0.8, lag = 20),
  specification('seasonal_ar5', 'ar', phi = 0.6, lag = 5),
  specification('seasonal_ar12', 'ar', phi = 0.6, lag = 12),
  specification('var_cross08', 'var', phi = 0.7, cross = 0.8),
  specification('var_cross15', 'var', phi = 0.7, cross = 1.5),
  specification('garch095', 'garch', arch = 0.05, garch = 0.90),
  specification('garch098_t5', 'garch', arch = 0.03, garch = 0.95, innovation = 't5'),
  specification('ar_garch', 'garch', phi = 0.5, arch = 0.05, garch = 0.90),
  specification('persistent_b05_pi25', 'ar', phi = 0.4, shift = 0.5),
  specification('persistent_b10_pi25', 'ar', phi = 0.4, shift = 1),
  specification('persistent_b10_pi50', 'ar', phi = 0.4, shift = 1, fraction = 0.5)
))
# The profile extension is a separate exploratory batch, defined after the
# stationary screen. It perturbs the two peaks of the paper's full-rank
# training path; it is not a variance-profile or information-clock experiment.
# All ten profiles and both signs of the transverse component are retained.
if (is_profile) {
  profile_rows <- list()
  peaks <- rbind(c(.7, 1), c(1, 1), c(1.3, 1), c(1, .7), c(1, 1.3))
  for (i in seq_len(nrow(peaks))) for (middle in c(.25, .75)) {
    left <- peaks[i, 1];right <- peaks[i, 2]
    id <- sprintf('profile_L%.1f_R%.1f_M%.2f', left, right, middle)
    profile_rows[[length(profile_rows) + 1L]] <- specification(
      id, 'profile',
      left = left, right = right, middle = middle
    )
  }
  dgp_grid <- do.call(rbind, profile_rows)
}
dgp_grid$index <- seq_len(nrow(dgp_grid))
write.csv(dgp_grid, file.path(out_dir, 'dgp_definitions.csv'), row.names = FALSE)
if (is_confirmation) {
  selected <- read.csv(file.path(root, 'results', paste0('dgp_screen', output_suffix), 'selected_dgps.csv'))
  dgp_grid <- dgp_grid[dgp_grid$id %in% selected$id, ]
  stopifnot(nrow(dgp_grid) > 0L)
}
filter_ids <- Sys.getenv('JBES_DGP_IDS', '')
if (nzchar(filter_ids)) dgp_grid <- dgp_grid[dgp_grid$id %in% strsplit(filter_ids, ',')[[1]], ]
stopifnot(nrow(dgp_grid) > 0L)

# Standardized innovations have variance one. The Gaussian mixture has a
# 5% component with standard deviation 5, so its variance divisor is 2.2.
innovation_sample <- function(n, distribution) {
  if (distribution == 't5') return(matrix(rt(2 * n, 5) / sqrt(5 / 3), n, 2))
  if (distribution == 'mixture') {
    scale <- ifelse(runif(2 * n) < 0.05, 5, 1)
    return(matrix(rnorm(2 * n) * scale / sqrt(2.2), n, 2))
  }
  matrix(rnorm(2 * n), n, 2)
}

stationary_noise <- function(n, design) {
  total <- n + burn_in
  innovations <- innovation_sample(total, design$innovation)
  if (design$family == 'garch') {
    variance <- rep(1, 2)
    previous <- rep(0, 2)
    shocks <- innovations
    for (t in seq_len(total)) {
      variance <- 1 - design$arch - design$garch +
        design$arch * previous^2 + design$garch * variance
      shocks[t, ] <- sqrt(variance) * innovations[t, ]
      previous <- shocks[t, ]
    }
    innovations <- shocks
  }
  if (design$family == 'var') {
    # Eigenvalues are phi and -0.2. Innovation loading I-A gives LRV I_2,
    # including for the nonnormal triangular transition matrix.
    transition <- matrix(c(design$phi, 0, design$cross, -0.2), 2, 2)
    loading <- diag(2) - transition
    observations <- matrix(0, total, 2)
    for (t in 2:total) observations[t, ] <-
      as.numeric(transition %*% observations[t - 1, ] + loading %*% innovations[t, ])
  } else {
    observations <- innovations
    if (design$family %in% c('ma', 'arma')) {
      lag <- design$lag
      previous <- rbind(matrix(0, lag, 2), innovations[seq_len(total - lag), , drop = FALSE])
      observations <- (innovations + design$theta * previous) / (1 + design$theta)
    }
    if (design$phi != 0) {
      coefficients <- c(rep(0, design$lag - 1), design$phi)
      observations <- apply(observations, 2, function(z)
        as.numeric(stats::filter((1 - design$phi) * z, coefficients, method = 'recursive')))
    }
  }
  observations[(burn_in + 1L):total, , drop = FALSE]
}

# Each alternative is fixed in advance. h is a local effect in LRV units;
# the actual per-observation mean change is h*v/sqrt(m).
scenarios <- expand.grid(
  break_fraction = c(0.1, 0.8), h = c(2, 4, 6),
  direction = if (dimension == 1L) 'first' else c('dense', 'first'), stringsAsFactors = FALSE
)
if (is_profile) scenarios <- expand.grid(
  break_fraction = .2, h = c(2, 4),
  direction = paste0('angle_', c(-60, -30, 0, 30, 60)), stringsAsFactors = FALSE
)
method_names <- c('RSMS', 'SSMS', 'CSMS_default')
if (is_confirmation) method_names <- c(method_names, 'CSMS_Andrews_QS', 'CSMS_NeweyWest')

one_path <- function(replication, sample_phase, design, m, thresholds = NULL) {
  configuration <- 10L * design$index + match(m, c(100, 250, 500))
  # Distinct blocks separate designs, lengths, phases and confirmation.
  seed <- 20261016 + as.integer(is_confirmation) * 1000000000 +
    configuration * 1000000 + sample_phase * 100000 + replication +
    as.integer(dimension == 1L) * 400000000 + as.integer(is_profile) * 200000000
  set.seed(seed)
  n_monitor <- as.integer(m * horizon)
  monitoring_index <- seq_len(n_monitor)
  observations <- stationary_noise(m + n_monitor, design)[, seq_len(dimension), drop = FALSE]
  if (is_profile) {
    # g is the cumulative centered training-mean profile. Its discrete
    # derivative gives each training mean; the final slope persists online.
    r <- (0:m) / m
    f <- approx(c(0, .25, .5, .75, 1),
      c(0, design$left, design$middle, design$right, 0),
      xout = r
    )$y
    transverse <- approx(c(0, .125, .25, .5, .75, .875, 1),
      c(0, 0, 1, 0, -1, 0, 0),
      xout = r
    )$y
    g <- cbind(f, transverse)
    training_means <- m * (g[-1, , drop = FALSE] - g[-nrow(g), , drop = FALSE])
    monitoring_means <- cbind(rep(-4 * design$right, n_monitor), rep(0, n_monitor))
    observations <- observations + rbind(training_means, monitoring_means)
  }
  if (design$shift != 0) {
    start <- floor(m * (1 - design$fraction)) + 1L
    observations[start:nrow(observations), 1] <-
      observations[start:nrow(observations), 1] + design$shift
  }
  training_mean <- colMeans(observations[seq_len(m), , drop = FALSE])
  scores <- sweep(observations, 2, training_mean)
  fit <- aersn_training_fit(scores[seq_len(m), , drop = FALSE])
  gauge <- aersn_gauge_evaluator(fit$hull)
  normalizers <- list(
    aersn::aersn_shao_normalizer(fit, integration = 'calendar')$matrix,
    hac_lrv(scores[seq_len(m), , drop = FALSE])
  )
  if (is_confirmation) {
    normalizers <- c(normalizers, list(
      aersn::aersn_hac_lrv(fit,
        kernel = 'Quadratic Spectral', bandwidth = 'andrews',
        prewhite = FALSE, adjust = FALSE, center = TRUE
      )$matrix,
      aersn::aersn_hac_lrv(fit,
        kernel = 'Bartlett', bandwidth = 'newey-west',
        prewhite = FALSE, adjust = FALSE, center = TRUE
      )$matrix
    ))
  }
  inverses <- lapply(normalizers, inv_pd, ridge = 1e-8)
  base_sum <- apply(scores[m + monitoring_index, , drop = FALSE], 2, cumsum) / sqrt(m)
  weight <- (1 + monitoring_index / m)^2 *
    (monitoring_index / (m + monitoring_index))^(2 * weight_exponent)
  statistic_path <- function(change_index = 0L, h = 0, direction = 'dense') {
    vector <- if (direction == 'dense') rep(1, dimension) / sqrt(dimension) else c(1, rep(0, dimension - 1L))
    if (grepl('^angle_', direction)) {
      angle <- as.numeric(sub('angle_', '', direction)) * pi / 180
      vector <- c(-cos(angle), sin(angle))
    }
    cumulative <- base_sum + (pmax(monitoring_index - change_index, 0) / m) %o% (h * vector)
    values <- cbind(gauge(cumulative)^2, vapply(inverses, function(inverse)
      rowSums((cumulative %*% inverse) * cumulative), numeric(n_monitor))) / weight
    stopifnot(all(is.finite(values)))
    values
  }
  if (sample_phase == 1L) return(apply(statistic_path(), 2, max))
  settings <- if (sample_phase == 2L) data.frame(break_fraction = 0, h = 0, direction = 'dense') else scenarios
  stops <- matrix(n_monitor + 1L, nrow(settings), length(method_names))
  for (setting in seq_len(nrow(settings))) {
    k_star <- floor(m * horizon * settings$break_fraction[setting])
    paths <- statistic_path(k_star, settings$h[setting], settings$direction[setting])
    for (method in seq_along(method_names))
      stops[setting, method] <- first_exceed(paths[, method], thresholds[method])
  }
  stops
}

summary_rows <- paired_rows <- calibration_rows <- list()
for (design_index in seq_len(nrow(dgp_grid))) for (m in training_lengths) {
  design <- dgp_grid[design_index, ]
  tag <- paste0(design$id, '_m', m)
  cat(phase_name, tag, 'calibration=', n_cal, 'evaluation=', n_eval, '\n');flush.console()
  calibration <- par_lapply_safe(seq_len(n_cal), function(i) one_path(i, 1L, design, m), n_cores)
  maxima <- do.call(rbind, calibration)
  thresholds <- apply(maxima, 2, quantile, probs = 1 - alpha, type = 8)
  saveRDS(maxima, file.path(out_dir, paste0(tag, '_calibration.rds')))
  calibration_rows[[length(calibration_rows) + 1L]] <- data.frame(
    dgp = design$id, m, method = method_names, threshold = thresholds, n_cal
  )
  for (sample_phase in 2:3) {
    records <- par_lapply_safe(
      seq_len(n_eval),
      function(i) one_path(i, sample_phase, design, m, thresholds), n_cores
    )
    saveRDS(records, file.path(out_dir, paste0(tag, '_phase', sample_phase, '.rds')))
    settings <- if (sample_phase == 2) data.frame(break_fraction = 0, h = 0, direction = 'dense') else scenarios
    for (setting in seq_len(nrow(settings))) {
      stops <- t(vapply(records, function(z) z[setting, ], numeric(length(method_names))))
      colnames(stops) <- method_names
      k_star <- floor(m * horizon * settings$break_fraction[setting])
      audit <- mc_event_audit(stops, k_star, m * horizon)
      metadata <- data.frame(
        dgp = design$id, m, sample_phase, settings[setting, ],
        training_shift = design$shift, training_profile = design$family,
        left = design$left, right = design$right, middle = design$middle
      )
      summary_rows[[length(summary_rows) + 1L]] <- cbind(
        metadata[rep(1, nrow(audit$events)), , drop = FALSE], audit$events
      )
      paired_rows[[length(paired_rows) + 1L]] <- cbind(
        metadata[rep(1, nrow(audit$paired)), , drop = FALSE], audit$paired
      )
    }
  }
  write.csv(do.call(rbind, summary_rows), file.path(out_dir, 'summary.csv'), row.names = FALSE)
  write.csv(do.call(rbind, paired_rows), file.path(out_dir, 'paired.csv'), row.names = FALSE)
  write.csv(do.call(rbind, calibration_rows), file.path(out_dir, 'calibration.csv'), row.names = FALSE)
}

if (!is_confirmation) {
  # Rank the minimum paired lower bound against BOTH baseline competitors.
  # This is a screening score, not a multiplicity-adjusted significance test.
  all_pairs <- do.call(rbind, paired_rows)
  all_events <- do.call(rbind, summary_rows)
  ranking <- list()
  for (id in dgp_grid$id) {
    null <- all_events[all_events$dgp == id & all_events$sample_phase == 2, ]
    size_ok <- all(abs(null$post - alpha) <= 2.58 * sqrt(alpha * (1 - alpha) / n_eval))
    pairs <- all_pairs[all_pairs$dgp == id & all_pairs$sample_phase == 3 & all_pairs$first == 'RSMS', ]
    keys <- unique(pairs[c('m', 'break_fraction', 'h', 'direction')])
    for (j in seq_len(nrow(keys))) {
      key <- keys[j, ]
      group <- pairs[pairs$m == key$m & pairs$break_fraction == key$break_fraction &
        pairs$h == key$h & pairs$direction == key$direction, ]
      ranking[[length(ranking) + 1L]] <- data.frame(id, key, size_ok,
        minimum_gap = min(group$post_difference),
        screening_score = min(group$post_difference - 1.96 * group$paired_se)
      )
    }
  }
  ranking <- do.call(rbind, ranking)
  ranking <- ranking[order(!ranking$size_ok, -ranking$screening_score), ]
  write.csv(ranking, file.path(out_dir, 'ranking.csv'), row.names = FALSE)
  selected <- head(ranking[!duplicated(ranking$id), ], 3)
  write.csv(selected, file.path(out_dir, 'selected_dgps.csv'), row.names = FALSE)
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, 'session.txt'))
cat('Finished', phase_name, '; all designs retained.\n')
