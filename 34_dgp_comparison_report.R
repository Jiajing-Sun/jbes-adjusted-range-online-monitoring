# Report independent confirmation at the settings frozen by screening.
#
# Each screen selects three DGPs and one scenario per DGP. Confirmation uses
# new random streams at m=250 and 500, with two automatic HAC comparators.
# This script does not select a different scenario after seeing confirmation.
# All full-grid outcomes remain in their original summary and paired CSVs.
# Paired standard errors condition on Monte Carlo calibration thresholds;
# the displayed intervals do not include threshold-estimation uncertainty.

root <- normalizePath(Sys.getenv('JBES_PROJECT_ROOT', getwd()))
sets <- c(bivariate = '', scalar = '_q1', training_profile = '_profile')
rows <- contrasts <- list()
for (batch in names(sets)) {
  suffix <- sets[[batch]]
  screen_dir <- file.path(root, 'results', paste0('dgp_screen', suffix))
  confirm_dir <- file.path(root, 'results', paste0('dgp_confirm', suffix))
  if (!file.exists(file.path(confirm_dir, 'session.txt')))
    stop('Confirmation is incomplete for ', batch, '; no final report was generated.')
  selected <- read.csv(file.path(screen_dir, 'selected_dgps.csv'))
  events <- read.csv(file.path(confirm_dir, 'summary.csv'))
  paired <- read.csv(file.path(confirm_dir, 'paired.csv'))
  for (i in seq_len(nrow(selected))) {
    scenario <- selected[i, ]
    keep <- function(x) x$dgp == scenario$id & x$sample_phase == 3 &
      x$break_fraction == scenario$break_fraction & x$h == scenario$h &
      x$direction == scenario$direction
    values <- events[keep(events), ]
    stopifnot(nrow(values) == 10L, all(values$n == 5000L))
    values$batch <- batch
    values$screening_score <- scenario$screening_score
    values$null_alarm <- values$null_standard_error <- NA_real_
    for (j in seq_len(nrow(values))) {
      null <- events[events$dgp == values$dgp[j] & events$m == values$m[j] &
        events$method == values$method[j] & events$sample_phase == 2, ]
      stopifnot(nrow(null) == 1L)
      values$null_alarm[j] <- null$post
      values$null_standard_error[j] <- null$post_se
    }
    # Keep a common schema across the original and follow-up batches.
    # Detailed profile parameters remain in each batch's dgp_definitions.csv.
    rows[[length(rows) + 1L]] <- values[c(
      'batch', 'dgp', 'm', 'break_fraction', 'h',
      'direction', 'method', 'n', 'pre_count', 'post_count', 'no_crossing_count',
      'post', 'post_se', 'survival_count', 'conditional_post', 'cadd', 'cadd_se',
      'null_alarm', 'null_standard_error', 'screening_score'
    )]
    comparison <- paired[keep(paired) & paired$first == 'RSMS', ]
    stopifnot(nrow(comparison) == 8L)
    comparison$batch <- batch
    comparison$conditional_lower <- comparison$post_difference - 1.96 * comparison$paired_se
    comparison$conditional_upper <- comparison$post_difference + 1.96 * comparison$paired_se
    contrasts[[length(contrasts) + 1L]] <- comparison[c(
      'batch', 'dgp', 'm',
      'break_fraction', 'h', 'direction', 'first', 'second', 'post_difference',
      'paired_se', 'cadd_difference', 'cadd_difference_se',
      'conditional_lower', 'conditional_upper'
    )]
  }
}
output <- file.path(root, 'results', 'dgp_report')
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write.csv(do.call(rbind, rows), file.path(output, 'frozen_scenario_confirmation.csv'), row.names = FALSE)
write.csv(do.call(rbind, contrasts), file.path(output, 'frozen_scenario_paired.csv'), row.names = FALSE)
writeLines(c(
  'The scenarios were frozen by screening; confirmation uses independent paths.',
  'Each of 18 settings has 10,000 calibration paths and separate sets of 5,000 null and alternative paths.',
  'Three batches are reported: 26 bivariate noise/training-contamination designs, their 26 scalar marginals, and ten follow-up full-rank training profiles.',
  'The profile batch was defined after the initial stationary screen. It is a separate exploratory extension.',
  'All outcomes and ranking reversals are retained. No confirmatory scenario is substituted for a screening choice.',
  'Intervals condition on estimated thresholds and are not simultaneous confidence bands.',
  'Compare detection with pre-change alarms, independent null size and CADD detection counts; smaller CADD alone does not establish an advantage.'
), file.path(output, 'interpretation.txt'))
cat('Frozen-scenario confirmation report written.\n')
