# Reconcile the completed PAR experiment with its saved simulation paths.
# Run from the replication root after par_null and par_alternative.
# This script reads existing paths; it never resamples observations or changes
# a threshold. All Monte Carlo errors condition on the fitted thresholds.

arguments <- commandArgs(trailingOnly = FALSE)
script <- sub('^--file=', '', arguments[grepl('^--file=', arguments)])
source(file.path(dirname(normalizePath(script[1])), '08_mc_audit.R'))

directory <- 'results'
reported <- read.csv(file.path(directory, 'par_alt_performance_grid.csv'))
methods <- c('RSMS', 'SSMS', 'CSMS')
design_columns <- c('m', 'T', 'gamma', 'break_frac', 'change_type', 'alpha2', 'lambda2')
designs <- unique(reported[design_columns])
stopifnot(nrow(designs) == 1728L, nrow(reported) == 5184L, all(reported$N == 1000L))

same <- function(actual, expected) {
  stopifnot(isTRUE(all.equal(actual, expected, tolerance = 1e-11, check.attributes = FALSE)))
}

events <- pairs <- list()
for (i in seq_len(nrow(designs))) {
  design <- designs[i, ]
  tag <- with(design, sprintf(
    'par_alt_m%d_T%d_g%.2f_bf%.1f_%s_a2%.2f_l2%.1f',
    m, T, gamma, break_frac, change_type, alpha2, lambda2
  ))
  saved <- readRDS(file.path(directory, 'event_audit', paste0(tag, '.rds')))
  stopifnot(length(saved$replications) == 1000L, saved$horizon == design$m * design$T)
  selected <- rep(TRUE, nrow(reported))
  for (column in design_columns) selected <- selected & reported[[column]] == design[[column]]
  row <- reported[selected, ]
  row <- row[match(methods, row$method), ]
  stopifnot(nrow(row) == 3L, !anyNA(row$method), all(row$k_star == saved$k_star))

  for (mode in c('asym', 'adj')) {
    field <- paste0('stop_', mode)
    stops <- sapply(methods, function(method) {
      vapply(saved$replications, function(path) unname(path[[method]][field]), numeric(1))
    })
    stopifnot(all(is.finite(stops)), all(stops >= 1L), all(stops <= saved$horizon + 1L))
    recalculated <- mc_event_audit(stops, saved$k_star, saved$horizon)
    event_file <- file.path(directory, 'event_audit', paste0(tag, '_', field, '_events.csv'))
    pair_file <- file.path(directory, 'event_audit', paste0(tag, '_', field, '_paired.csv'))
    same(recalculated$events, read.csv(event_file))
    same(recalculated$paired, read.csv(pair_file))

    counts <- recalculated$events
    stopifnot(all(counts$pre_count + counts$post_count + counts$no_crossing_count == 1000L))
    same(counts$post, row[[paste0('DP_post_', mode)]])
    same(counts$pre_count / 1000, row[[paste0('PFA_pre_', mode)]])
    same((counts$pre_count + counts$post_count) / 1000, row[[paste0('power_', mode)]])
    same(counts$cadd, row[[paste0('CADD_', mode)]])
    # Each delay mean uses that method's own post-change detections.
    # A shared-detection subset would estimate a different conditional mean.
    stopifnot(all(is.na(counts$cadd) == (counts$post_count == 0L)))
    event_design <- design[rep(1L, nrow(counts)), , drop = FALSE]
    pair_design <- design[rep(1L, nrow(recalculated$paired)), , drop = FALSE]
    rownames(event_design) <- rownames(pair_design) <- NULL
    events[[length(events) + 1L]] <- cbind(event_design, mode, counts)
    pairs[[length(pairs) + 1L]] <- cbind(pair_design, mode, recalculated$paired)
  }
}

events <- do.call(rbind, events)
pairs <- do.call(rbind, pairs)
stopifnot(nrow(events) == 10368L, nrow(pairs) == 10368L)
write.csv(events, file.path(directory, 'par_all_events_verified.csv'), row.names = FALSE)
write.csv(pairs, file.path(directory, 'par_all_pairs_verified.csv'), row.names = FALSE)

# The no-change parameter pair is retained in the full output but excluded
# from detection-power rankings. Comparisons use unrounded paired estimates.
genuine <- pairs[pairs$mode == 'adj' & (pairs$alpha2 != .2 | pairs$lambda2 != 2), ]
summary_rows <- list()
for (view in c('full_grid', 'displayed_m250_T1_abrupt')) {
  sample <- if (view == 'full_grid') genuine else genuine[
    genuine$m == 250 & genuine$T == 1 & genuine$change_type == 'abrupt',
  ]
  for (second in c('SSMS', 'CSMS')) {
    comparison <- sample[sample$first == 'RSMS' & sample$second == second, ]
    # Integer count differences avoid turning floating-point roundoff into a ranking.
    difference <- round(comparison$post_difference * 1000)
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      view,
      first = 'RSMS', second, configurations = length(difference),
      positive = sum(difference > 0), equal = sum(difference == 0),
      negative = sum(difference < 0), mean_difference_pp = mean(difference) / 10,
      minimum_paired_se_pp = min(comparison$paired_se) * 100,
      maximum_paired_se_pp = max(comparison$paired_se) * 100
    )
  }
}
write.csv(do.call(rbind, summary_rows), file.path(directory, 'par_detection_comparison.csv'), row.names = FALSE)
cat('Verified 1,728 designs, 5,184 method rows and 10,368 event and paired-comparison rows.\n')
