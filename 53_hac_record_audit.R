# Audit saved fixed and automatic HAC sensitivity experiments.
# Run after the mean, bandwidth, kernel and automatic-bandwidth stages.
# This script does not simulate new data or change any calibrated threshold.
#
# The RSMS benchmark comes from the separate baseline mean experiment.
# Consequently, RSMS-minus-CSMS uncertainty uses independent-sample formulas,
# not a paired-path formula. All standard errors condition on fitted thresholds.

results <- 'results'
event_rows <- comparison_rows <- list()
reference_count <- 0L
for (module in c('bandwidth', 'kernel')) {
  directory <- file.path(results, paste0('csms_', module, '_records'))
  reference_files <- list.files(directory, pattern = '_reference\\.rds$', full.names = TRUE)
  stopifnot(length(reference_files) == if (module == 'bandwidth') 12L else 6L)
  reported <- read.csv(file.path(results, paste0('csms_', module, '_sensitivity.csv')))
  for (file in reference_files) {
    reference <- readRDS(file)
    stopifnot(
      reference$n_null == 1000L, reference$n_calibration == 1000L,
      nrow(reference$null_stops) == 1000L,
      length(reference$calibration_maxima) == 1000L,
      all(is.finite(reference$calibration_maxima))
    )
    fitted <- unname(quantile(reference$calibration_maxima, .95, type = 8))
    stopifnot(isTRUE(all.equal(fitted, reference$estimated_threshold)))
    if (module == 'bandwidth') {
      index <- as.numeric(sub('.*_eta', '', reference$tag))
      row <- reported[reported$method == 'CSMS' & reported$T == reference$T &
        !is.na(reported$eta) & reported$eta == index, ]
    } else {
      index <- sub('.*_', '', reference$tag)
      row <- reported[reported$method == 'CSMS' & reported$T == reference$T &
        !is.na(reported$kernel) & reported$kernel == index, ]
    }
    stopifnot(nrow(row) == 1L)
    horizon <- reference$m * reference$T
    stopifnot(
      abs(mean(reference$null_stops[, 1] <= horizon) - row$size_05) < 1e-12,
      abs(mean(reference$null_stops[, 2] <= horizon) - row$size_10) < 1e-12
    )
    reference_count <- reference_count + 1L
    for (timing in c('early', 'late')) {
      alternative <- readRDS(sub('_reference.rds$', paste0('_', timing, '.rds'), file))
      stops <- alternative$stops
      n <- length(stops)
      k_star <- alternative$k_star
      stopifnot(
        n == 1000L, alternative$n_evaluation == n,
        all(is.finite(stops)), all(stops >= 1L), all(stops <= horizon + 1L)
      )
      pre <- stops <= k_star
      post <- stops > k_star & stops <= horizon
      delay <- stops[post] - k_star
      probability <- mean(post)
      cadd <- if (length(delay)) mean(delay) else NA_real_
      cadd_se <- if (length(delay) > 1L) sd(delay) / sqrt(length(delay)) else NA_real_
      stopifnot(
        sum(pre) + sum(post) + sum(stops > horizon) == n,
        abs(probability - row[[paste0(timing, '_DP_post_adj')]]) < 1e-12,
        abs(mean(pre) - row[[paste0(timing, '_PFA_pre_adj')]]) < 1e-12,
        isTRUE(all.equal(cadd, row[[paste0(timing, '_CADD_adj')]]))
      )
      event_rows[[length(event_rows) + 1L]] <- data.frame(
        module,
        tag = reference$tag, timing, n,
        pre_count = sum(pre), post_count = sum(post), no_crossing_count = sum(stops > horizon),
        post = probability, post_se = sqrt(probability * (1 - probability) / n), cadd, cadd_se
      )
      fraction <- if (timing == 'early') .1 else .8
      magnitude <- if (timing == 'early') .5 else 1.5
      baseline_tag <- sprintf(
        'mean_alt_DGP1_m250_T%d_g0.15_bf%.1f_abrupt_D%.2f',
        reference$T, fraction, magnitude
      )
      baseline <- read.csv(file.path(
        results, 'event_audit',
        paste0(baseline_tag, '_stop_adj_events.csv')
      ))
      baseline <- baseline[baseline$method == 'RSMS', ]
      stopifnot(nrow(baseline) == 1L, baseline$n == 1000L)
      # These independent-sample errors include evaluation noise in both
      # procedures but exclude uncertainty from either calibration batch.
      comparison_rows[[length(comparison_rows) + 1L]] <- data.frame(
        module,
        tag = reference$tag, timing,
        post_gap_R_minus_C = baseline$post - probability,
        independent_se = sqrt(baseline$post_se^2 + probability * (1 - probability) / n),
        cadd_gap_R_minus_C = baseline$cadd - cadd,
        cadd_gap_se = sqrt(baseline$cadd_se^2 + cadd_se^2)
      )
    }
  }
}
stopifnot(reference_count == 18L, length(event_rows) == 36L)
write.csv(do.call(rbind, event_rows), file.path(results, 'fixed_hac_event_audit.csv'), row.names = FALSE)
write.csv(do.call(rbind, comparison_rows), file.path(results, 'fixed_hac_independent_comparisons.csv'), row.names = FALSE)
cat('Verified 18 independent reference batches and 36 alternative batches, each with 1,000 paths.\n')

# Automatic selectors use common evaluation paths, so these comparisons ARE
# paired. Keep their uncertainty separate from the fixed-tuning rows above.
script_argument <- grep('^--file=', commandArgs(), value = TRUE)
stopifnot(length(script_argument) == 1L)
script_directory <- dirname(normalizePath(sub('^--file=', '', script_argument)))
source(file.path(script_directory, '08_mc_audit.R'))
automatic <- file.path(results, 'automatic_bandwidth')
critical_values <- read.csv(file.path(automatic, 'selector_critical_values.csv'))
reported <- read.csv(file.path(automatic, 'selector_results.csv'))
automatic_events <- automatic_pairs <- list()
for (horizon_units in c(1, 5)) {
  critical <- critical_values[critical_values$T == horizon_units, ]
  methods <- critical$method
  calibration <- readRDS(file.path(automatic, paste0('calibration_T', horizon_units, '.rds')))
  maxima <- do.call(rbind, lapply(calibration, `[[`, 'maxima'))
  stopifnot(
    length(calibration) == 1000L, ncol(maxima) == 8L,
    max(abs(apply(maxima, 2, quantile, probs = .95, type = 8) - critical$chat)) < 1e-10
  )
  for (phase in 2:4) {
    records <- readRDS(file.path(automatic, paste0('evaluation_T', horizon_units, '_phase', phase, '.rds')))
    stopifnot(length(records) == 1000L)
    horizon <- 250L * horizon_units
    k_star <- if (phase == 2L) 0L else floor(horizon * if (phase == 3L) .1 else .8)
    for (boundary in 1:2) {
      stops <- t(vapply(records, function(record) record$stops[, boundary], numeric(8)))
      colnames(stops) <- methods
      audit <- mc_event_audit(stops, k_star, horizon)
      metadata <- data.frame(
        T = horizon_units, phase = c('null', 'early', 'late')[phase - 1L],
        boundary = c('asymptotic', 'independent_calibration')[boundary]
      )
      for (method in methods) {
        displayed <- reported[reported$T == horizon_units & reported$phase == metadata$phase &
          reported$boundary == metadata$boundary & reported$method == method, ]
        event <- audit$events[audit$events$method == method, ]
        stopifnot(
          nrow(displayed) == 1L,
          abs(event$post - displayed$DP_post) < 1e-12,
          abs(event$pre_count / event$n - displayed$PFA_pre) < 1e-12
        )
      }
      automatic_events[[length(automatic_events) + 1L]] <- cbind(metadata, audit$events)
      automatic_pairs[[length(automatic_pairs) + 1L]] <- cbind(metadata, audit$paired)
    }
  }
}
write.csv(do.call(rbind, automatic_events), file.path(results, 'automatic_hac_event_audit.csv'), row.names = FALSE)
write.csv(do.call(rbind, automatic_pairs), file.path(results, 'automatic_hac_paired_comparisons.csv'), row.names = FALSE)
cat('Verified two automatic-selector calibration batches and six evaluation batches, with paired comparisons.\n')
