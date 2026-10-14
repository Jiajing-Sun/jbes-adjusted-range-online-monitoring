# Report the frozen follow-up candidates without changing their alternatives.
# All five methods, both training lengths and unfavorable contrasts remain.
# Reported intervals condition on the independent calibration thresholds.
# They are descriptive pointwise Monte Carlo intervals, not familywise
# inference after screening or uncertainty over estimated calibration draws.
root <- normalizePath(Sys.getenv('JBES_PROJECT_ROOT', getwd()))
base <- file.path(root, 'results')
screen <- file.path(base, 'dgp_focused_screen')
confirm <- file.path(base, 'dgp_focused_confirm')
selected <- read.csv(file.path(screen, 'selected_dgps.csv'))
events <- read.csv(file.path(confirm, 'summary.csv'))
paired <- read.csv(file.path(confirm, 'paired.csv'))
calibration <- read.csv(file.path(confirm, 'calibration.csv'))
stopifnot(
  nrow(selected) == 3L, nrow(events) == 60L, nrow(paired) == 120L,
  nrow(calibration) == 30L, all(events$n == 5000),
  all(calibration$n_cal == 10000), all(events$m %in% c(250, 500))
)
stopifnot(all(events$pre_count + events$post_count + events$no_crossing_count == events$n))
for (i in seq_len(nrow(selected))) {
  a <- subset(events, dgp == selected$id[i] & sample_phase == 3)
  stopifnot(
    all(a$break_fraction == selected$break_fraction[i]),
    all(a$h == selected$h[i]), all(a$direction == selected$direction[i])
  )
}
out <- file.path(base, 'dgp_focused_report')
dir.create(out, showWarnings = FALSE)
write.csv(events, file.path(out, 'all_confirmation_events.csv'), row.names = FALSE)
contrasts <- subset(paired, sample_phase == 3 & first == 'RSMS')
contrasts$conditional_lower <- contrasts$post_difference - 1.96 * contrasts$paired_se
contrasts$conditional_upper <- contrasts$post_difference + 1.96 * contrasts$paired_se
write.csv(contrasts, file.path(out, 'rsms_paired_comparisons.csv'), row.names = FALSE)
summary <- list()
for (id in selected$id) for (m0 in c(250, 500)) {
  a <- subset(events, dgp == id & m == m0 & sample_phase == 3)
  z <- subset(events, dgp == id & m == m0 & sample_phase == 2)
  p <- subset(contrasts, dgp == id & m == m0)
  original <- subset(p, second %in% c('SSMS', 'CSMS_default'))
  one <- data.frame(
    dgp = id, m = m0, break_fraction = a$break_fraction[1], h = a$h[1],
    direction = a$direction[1], gamma = a$gamma[1],
    minimum_gap_original = min(original$post_difference),
    minimum_gap_all = min(p$post_difference),
    minimum_conditional_lower_original = min(original$conditional_lower),
    minimum_conditional_lower_all = min(p$conditional_lower),
    largest_comparator = p$second[which.min(p$post_difference)]
  )
  for (method in a$method) {
    aa <- a[a$method == method, ];zz <- z[z$method == method, ]
    one[[paste0('post_', method)]] <- aa$post
    one[[paste0('pre_', method)]] <- aa$pre_count / aa$n
    one[[paste0('null_', method)]] <- zz$post
    one[[paste0('cadd_', method)]] <- aa$cadd
  }
  summary[[length(summary) + 1L]] <- one
}
summary <- do.call(rbind, summary)
write.csv(summary, file.path(out, 'frozen_candidate_summary.csv'), row.names = FALSE)
writeLines(
  c(
    'These are the exact three alternatives frozen after the complete focused screen.',
    'The grid was designed after an earlier search; it is an exploratory follow-up.',
    'The independent confirmation uses new calibration, null and alternative paths.',
    'Original contrasts compare RSMS with SSMS and the prespecified default CSMS.',
    'All-method contrasts also include Andrews-QS and Newey-West-Bartlett CSMS.',
    'Model-specific calibration is an oracle benchmark, not feasible unknown-DGP calibration.',
    'Intervals condition on the estimated thresholds and exclude calibration uncertainty.',
    'A higher detection probability does not establish a smaller conditional detection delay.',
    'No results from this follow-up have been inserted into the manuscript.'
  ),
  file.path(out, 'interpretation.txt')
)
cat('Reported all six frozen candidate/length combinations and 24 RSMS contrasts.\n')
