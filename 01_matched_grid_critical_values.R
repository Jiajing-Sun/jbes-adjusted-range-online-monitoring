# Same Gaussian reference construction, with the actual training/monitoring grids.
# Independent of every DGP-specific size/power evaluation stream.
root <- normalizePath(getwd())
stopifnot(file.exists(file.path(root, '02_methods.R')))
Sys.setenv(JBES_PROJECT_ROOT = root)
allrows <- list()
for (m in c(50, 75, 100, 250, 500)) {
  cat('Matched grid m=', m, '; 10000 reference draws\n');flush.console()
  env <- c(
    'JBES_CV_D_GRID=1,2',
    paste0('JBES_CV_N_GRID_PER_UNIT=', m), 'JBES_CV_N_REP=10000', 'JBES_N_CORES=2',
    paste0('JBES_CV_SUBDIR=matched_grid/m', m)
  )
  log <- file.path(root, paste0('logs_matched_grid_m', m, '.txt'))
  status <- system2(file.path(R.home('bin'), 'Rscript'),
    shQuote(file.path(root, '00_finite_critical_values.R')),
    env = env, stdout = log, stderr = log
  )
  if (status != 0) stop('Matched reference failed for m=', m)
  x <- read.csv(file.path(root, 'results', 'matched_grid', paste0('m', m), 'affine_finite_critical_values.csv'))
  x$m <- m;x$calibration <- 'matched_grid';allrows[[length(allrows) + 1L]] <- x
}
write.csv(do.call(rbind, allrows), file.path(root, 'results', 'affine_matched_grid_critical_values.csv'), row.names = FALSE)
cat('Matched-grid critical values complete.\n')
