# Combine completed production runs without changing the q=1,2 reference
# used by the finite-sample experiments and the empirical application.
root <- normalizePath(Sys.getenv('JBES_PROJECT_ROOT', unset = getwd()))
keys <- c('method', 'd', 'T', 'gamma', 'alpha')
for (kind in c('finite', 'open_end')) {
  name <- paste0('affine_', kind, '_critical_values.csv')
  current <- read.csv(file.path(root, 'results', name))
  high <- read.csv(file.path(root, 'results', 'full_dimensions', name))
  keys <- intersect(c('method', 'd', 'T', 'gamma', 'alpha'), names(current))
  low <- current[current$d %in% 1:2, , drop = FALSE]
  stopifnot(
    setequal(unique(high$d), 3:20), all(high$n_rep == 10000),
    all(low$n_rep == 10000), all(is.finite(high$critical_value)),
    all(high$critical_value > 0)
  )
  expected <- if (kind == 'finite') 48L else 12L
  stopifnot(
    nrow(low) == 2L * expected, nrow(high) == 18L * expected,
    identical(names(low), names(high))
  )
  merged <- rbind(low, high)
  stopifnot(!anyDuplicated(merged[keys]))
  merged <- merged[do.call(order, unname(as.list(merged[keys]))), , drop = FALSE]
  target <- file.path(root, 'results', name)
  temporary <- paste0(target, '.tmp')
  write.csv(merged, temporary, row.names = FALSE)
  stopifnot(file.rename(temporary, target))
  cat('Merged complete reference:', name, ';', nrow(merged), 'rows\n')
}
