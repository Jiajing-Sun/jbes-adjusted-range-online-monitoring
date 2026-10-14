# Focused follow-up to the completed DGP exploration.
#
# This grid was chosen after inspecting the original results. It is a NEW
# exploratory batch, not a confirmatory reanalysis of those same draws.
# Seven models cross two common monitoring weights. All five methods are
# calibrated on independent null paths under the SAME model and weight.
# Model-specific null calibration is an oracle simulation benchmark; it is
# not a proposed calibration method for unknown real-world training changes.
#
# Screening keeps every direction, effect and change location. The minimum
# paired lower bound against all four competitors ranks the configurations.
# The top three distinct model/weight settings are frozen, with their exact
# alternative, before fresh confirmation at m=250 and m=500. Screening
# lower bounds are ranking devices, not multiplicity-adjusted inference.
# Independent calibration, null and alternative records are saved locally.
# Conditional Monte Carlo errors exclude calibration-threshold uncertainty.
#
# Run with JBES_DGP_PHASE=screen, then JBES_DGP_PHASE=confirm.
# JBES_FOCUS_OUTPUT can name a separate temporary directory for smoke tests.
# Cached phases are reused only when the full run manifest is identical.

root <- normalizePath(Sys.getenv('JBES_PROJECT_ROOT', getwd()))
for (filename in c('00_utils.R', '01_critical_values.R', '02_methods.R', '08_mc_audit.R'))
  source(file.path(root, filename))
jbes_set_root(root)
phase_name <- Sys.getenv('JBES_DGP_PHASE', 'screen')
stopifnot(phase_name %in% c('screen', 'confirm'))
is_confirmation <- phase_name == 'confirm'
n_cal <- env_int('JBES_DGP_CAL', if (is_confirmation) 10000L else 1000L)
n_eval <- env_int('JBES_DGP_EVAL', if (is_confirmation) 5000L else 1000L)
n_cores <- env_int('JBES_N_CORES', 2L)
training_lengths <- env_num_vec('JBES_DGP_M', if (is_confirmation) c(250, 500) else 250)
stopifnot(n_cal < 100000L, n_eval < 100000L, all(training_lengths %in% c(100, 250, 500)))
dimension <- 2L
horizon <- 2
alpha <- .05
burn_in <- 2000L
output_base <- Sys.getenv('JBES_FOCUS_OUTPUT', file.path(root, 'results', 'dgp_focused'))
out_dir <- paste0(output_base, '_', phase_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

specification <- function(id, family, phi = 0, cross = 0, left = 1) {
  data.frame(id, family, phi,
    theta = 0, lag = 1, innovation = 'normal', arch = 0, garch = 0,
    cross, shift = 0, fraction = .25, left, right = 1, middle = .75
  )
}
# The profile is the centered cumulative training-mean curve used in the
# previous exploration. Its final slope persists during monitoring. The
# left peak changes while the right peak and transverse curve stay fixed.
base_models <- do.call(rbind, list(
  specification('iid_reference', 'ar'),
  specification('var_cross15', 'var', phi = .7, cross = 1.5),
  specification('var_cross30', 'var', phi = .7, cross = 3),
  specification('profile_L10', 'profile', left = 1),
  specification('profile_L13', 'profile', left = 1.3),
  specification('profile_L16', 'profile', left = 1.6),
  specification('profile_L20', 'profile', left = 2)
))
dgp_grid <- do.call(rbind, lapply(c(0, .15), function(gamma) {
  z <- base_models;z$gamma <- gamma
  z$id <- paste0(z$id, sprintf('_g%.2f', gamma));z
}))
dgp_grid$index <- seq_len(nrow(dgp_grid))
all_scenarios <- expand.grid(
  break_fraction = c(.1, .4), h = c(2, 4, 6),
  direction = paste0('angle_', c(-90, -60, -30, 0, 30, 45, 60, 75, 90)),
  stringsAsFactors = FALSE
)
method_names <- c('RSMS', 'SSMS', 'CSMS_default', 'CSMS_Andrews_QS', 'CSMS_NeweyWest')
selected <- NULL
if (is_confirmation) {
  selected <- read.csv(paste0(output_base, '_screen/selected_dgps.csv'))
  stopifnot(nrow(selected) == 3L, !anyDuplicated(selected$id))
  dgp_grid <- dgp_grid[dgp_grid$id %in% selected$id, ]
}
filter_ids <- Sys.getenv('JBES_DGP_IDS', '')
if (nzchar(filter_ids)) dgp_grid <- dgp_grid[dgp_grid$id %in% strsplit(filter_ids, ',')[[1]], ]
stopifnot(nrow(dgp_grid) > 0L)
# The manifest prevents a smoke run, changed design, or changed backend from
# being mistaken for a completed production phase. RDS writes are atomic.
source_files <- c(
  file.path(root, c(
    '00_utils.R', '01_critical_values.R',
    '02_methods.R', '07_aersn_backend.R', '08_mc_audit.R'
  )),
  file.path(root, '50_focused_dgp_search.R')
)
manifest <- list(
  phase = phase_name, models = dgp_grid, scenarios = all_scenarios,
  selected = selected, n_cal = n_cal, n_eval = n_eval, m = training_lengths, workers = n_cores,
  alpha = alpha, horizon = horizon, burn_in = burn_in, methods = method_names,
  source_md5 = tools::md5sum(source_files), R = R.version.string,
  aersn = as.character(packageVersion('aersn')),
  highs = as.character(packageVersion('highs'))
)
manifest_path <- file.path(out_dir, 'manifest.rds')
if (file.exists(manifest_path)) {
  if (!identical(manifest, readRDS(manifest_path))) stop('Checkpoint manifest mismatch: ', out_dir)
} else saveRDS(manifest, manifest_path)
write.csv(dgp_grid, file.path(out_dir, 'dgp_definitions.csv'), row.names = FALSE)
write.csv(all_scenarios, file.path(out_dir, 'screening_scenarios.csv'), row.names = FALSE)
checkpoint <- function(name, compute) {
  path <- file.path(out_dir, paste0(name, '.rds'))
  if (file.exists(path)) return(readRDS(path))
  value <- compute()
  temporary <- paste0(path, '.tmp');saveRDS(value, temporary)
  stopifnot(file.rename(temporary, path));value
}

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
    # Triangular recursion: filtering the second component first is
    # algebraically equivalent to A %*% X[t-1,] + (I-A) %*% e[t,].
    forcing <- innovations %*% t(loading)
    forcing[1, ] <- 0
    second <- as.numeric(stats::filter(forcing[, 2], -.2, method = 'recursive'))
    first <- as.numeric(stats::filter(forcing[, 1] +
      design$cross * c(0, head(second, -1)), design$phi, method = 'recursive'))
    observations <- cbind(first, second)
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

one_path <- function(replication, sample_phase, design, m, thresholds = NULL) {
  # Each replication has its own seed; worker scheduling cannot change it.
  # The follow-up block is disjoint from the original three screening batches.
  seed <- 700000000 + as.integer(is_confirmation) * 1000000000 +
    design$index * 10000000 + match(m, c(100, 250, 500)) * 1000000 +
    sample_phase * 100000 + replication
  RNGkind('Mersenne-Twister', 'Inversion', 'Rejection')
  set.seed(seed)
  n_monitor <- as.integer(m * horizon)
  monitoring_index <- seq_len(n_monitor)
  observations <- stationary_noise(m + n_monitor, design)[, seq_len(dimension), drop = FALSE]
  if (design$family == 'profile') {
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
  inverses <- lapply(normalizers, inv_pd, ridge = 1e-8)
  base_sum <- apply(scores[m + monitoring_index, , drop = FALSE], 2, cumsum) / sqrt(m)
  weight <- (1 + monitoring_index / m)^2 *
    (monitoring_index / (m + monitoring_index))^(2 * design$gamma)
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
  scenarios <- all_scenarios
  if (is_confirmation) scenarios <- selected[selected$id == design$id,
    c('break_fraction', 'h', 'direction'),
    drop = FALSE
  ]
  cat(phase_name, tag, 'calibration=', n_cal, 'evaluation=', n_eval, '\n');flush.console()
  maxima <- checkpoint(paste0(tag, '_calibration'), function()
    do.call(rbind, par_lapply_safe(
      seq_len(n_cal),
      function(i) one_path(i, 1L, design, m), n_cores
    )))
  thresholds <- apply(maxima, 2, quantile, probs = 1 - alpha, type = 8)
  calibration_rows[[length(calibration_rows) + 1L]] <- data.frame(
    dgp = design$id, m, method = method_names, threshold = thresholds, n_cal
  )
  for (sample_phase in 2:3) {
    cat(tag, 'phase', sample_phase, '\n');flush.console()
    records <- checkpoint(paste0(tag, '_phase', sample_phase), function()
      par_lapply_safe(
        seq_len(n_eval),
        function(i) one_path(i, sample_phase, design, m, thresholds), n_cores
      ))
    settings <- if (sample_phase == 2) data.frame(break_fraction = 0, h = 0, direction = 'dense') else scenarios
    for (setting in seq_len(nrow(settings))) {
      stops <- t(vapply(records, function(z) z[setting, ], numeric(length(method_names))))
      colnames(stops) <- method_names
      k_star <- floor(m * horizon * settings$break_fraction[setting])
      audit <- mc_event_audit(stops, k_star, m * horizon)
      metadata <- data.frame(
        dgp = design$id, m, sample_phase, settings[setting, ],
        gamma = design$gamma, training_shift = design$shift, training_profile = design$family,
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
  # Rank the minimum paired lower bound against ALL FOUR competitors.
  # This is a screening score, not a multiplicity-adjusted significance test.
  all_pairs <- do.call(rbind, paired_rows)
  all_events <- do.call(rbind, summary_rows)
  ranking <- list()
  for (id in dgp_grid$id) {
    null <- all_events[all_events$dgp == id & all_events$sample_phase == 2, ]
    size_ok <- all(abs(null$post - alpha) <= 2.58 * sqrt(alpha * (1 - alpha) * (1 / n_eval + 1 / n_cal)))
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
  ranking <- ranking[order(!ranking$size_ok, -ranking$screening_score, ranking$id), ]
  write.csv(ranking, file.path(out_dir, 'ranking.csv'), row.names = FALSE)
  selected <- head(ranking[!duplicated(ranking$id), ], 3)
  write.csv(selected, file.path(out_dir, 'selected_dgps.csv'), row.names = FALSE)
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, 'session.txt'))
cat('Finished', phase_name, '; all designs retained.\n')
