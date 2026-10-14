# End-to-end check in a temporary workspace, leaving production output untouched.
repo <- normalizePath(getwd())
# Resolve library paths before changing directories. Child R sessions must
# still find a package installed through a relative R_LIBS_USER setting.
Sys.setenv(R_LIBS_USER = paste(.libPaths(), collapse = .Platform$path.sep))
# Keep a failed workspace available for inspection; remove it on success.
work <- tempfile('jbes-smoke-', tmpdir = repo);dir.create(work)
files <- list.files(repo, pattern = '\\.R$', full.names = TRUE);invisible(file.copy(files, work))
oldwd <- getwd();setwd(work)
Sys.setenv(
  JBES_PROJECT_ROOT = work, JBES_N_CORES = 2, JBES_CV_D_GRID = '1,2', JBES_CV_T_GRID = '1,2',
  JBES_CV_N_REP = 100, JBES_CV_N_REP_OPEN = 100, JBES_CV_N_GRID_PER_UNIT = 80,
  JBES_CV_N_GRID_TRAIN_OPEN = 80, JBES_CV_N_GRID_MON_OPEN = 120,
  JBES_N_MC = 20, JBES_N_CAL = 20, JBES_N_MC_ALT = 20, JBES_DGP_GRID = 'DGP1',
  JBES_M_GRID = 40, JBES_T_GRID = 1, JBES_GAMMA_GRID = 0
)
run <- function(file) {
  status <- system2(file.path(R.home('bin'), 'Rscript'), shQuote(file));if (status != 0) stop('Smoke stage failed: ', file, '; retained work: ', work)
}
for (file in c('00_finite_critical_values.R', '00_open_end_critical_values.R', '01_mean_null_calibrate_and_size.R', '02_mean_alternatives_power_delay.R')) run(file)
stopifnot(file.exists('results/mean_null_size.csv'), file.exists('results/mean_alt_performance.csv'))
a <- read.csv('results/mean_null_size.csv');b <- read.csv('results/mean_alt_performance.csv');stopifnot(nrow(a) > 0, nrow(b) > 0)
setwd(oldwd);unlink(work, recursive = TRUE)
cat('Finite/open references, independent null calibration and alternative pipeline passed.\n')
