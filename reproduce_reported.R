# Fresh R processes preserve module isolation and the production RNG streams.
# Default: list the plan. --run executes all; --stage NAME executes one stage.
root <- normalizePath(Sys.getenv('JBES_PROJECT_ROOT', unset=getwd()))
stopifnot(file.exists(file.path(root,'00_utils.R')))
args <- commandArgs(trailingOnly=TRUE)
step <- function(name, script, workers=1L, settings=character(), seed='none') {
  list(name=name,script=script,workers=workers,settings=settings,seed=seed)
}
plan <- list(
 step('finite_q1q2','00_finite_critical_values.R',8L,c(JBES_CV_D_GRID='1,2'),13579),
 step('open_q1q2','00_open_end_critical_values.R',4L,c(JBES_CV_D_GRID='1,2'),975246),
 step('finite_q3q20','00_finite_critical_values.R',6L,c(JBES_CV_D_GRID=paste(3:20,collapse=','),JBES_CV_CHECKPOINT='1',JBES_CV_SUBDIR='full_dimensions'),13579),
 step('open_q3q20','00_open_end_critical_values.R',4L,c(JBES_CV_D_GRID=paste(3:20,collapse=','),JBES_CV_CHECKPOINT='1',JBES_CV_SUBDIR='full_dimensions'),975246),
 step('merge_reference','02_merge_reported_reference.R'),
 step('reference_assets','39_critical_value_assets.R'),
 step('matched_reference','01_matched_grid_critical_values.R',2L,seed=13579),
 step('mean_null','01_mean_null_calibrate_and_size.R',6L,seed=12345),
 step('mean_alternative','02_mean_alternatives_power_delay.R',6L,seed=24680),
 step('mean_summary','03_mean_summary_tables.R'),
 step('mean_plots','04_mean_plots.R'),
 step('mean_tables','35_mean_alt_latex_tables.R'),
 step('par_null','11_par_null_calibrate_and_size.R',4L,seed=54321),
 step('null_tables','34_null_latex_tables.R'),
 step('par_alternative','12_par_alternatives_power_delay.R',2L,seed=13579),
 step('par_plots','13_par_plots.R'),
 step('par_tables','40_par_grid_tables.R'),
 step('moderate','21_moderate_dim_sparse.R',2L,seed=11223),
 step('moderate_assets','36_moderate_assets.R'),
 step('dependence','22_dependence_sensitivity.R',2L,seed=33445),
 step('dependence_assets','37_dependence_assets.R'),
 step('contamination','23_training_contamination.R',2L,seed=77889),
 step('contamination_assets','38_contamination_assets.R'),
 step('persistent_training','28_persistent_training.R',8L,seed=91427000),
 step('bivariate_theory','31_bivariate_theory_checks.R',4L,seed=20261015),
 step('scalar_theory','30_scalar_theory_checks.R',4L,seed=20261014),
 step('persistent_assets','49_persistent_assets.R'),
 step('bandwidth','24_csms_bandwidth_sensitivity.R',1L,seed=97531),
 step('bandwidth_assets','41_csms_bandwidth_assets.R'),
 step('kernel','25_csms_kernel_sensitivity.R',1L,seed=86420),
 step('kernel_assets','42_csms_kernel_assets.R'),
 step('automatic_bandwidth','29_csms_automatic_bandwidth.R',4L,seed=20260914),
 step('matched_size','27_matched_grid_size_comparison.R',6L,seed=12345),
 step('empirical_plots','31_empirical_plots.R'),
 step('empirical_surfaces','43_empirical_3d_snapshots.R'),
 step('empirical_monitoring','32_empirical_monitoring_table.R'),
 step('empirical_tables','33_empirical_latex_tables.R'),
 step('empirical_bootstrap','34_empirical_bootstrap_check.R',1L,seed=20260913),
 step('multimedia_mesh','47_fx_multimedia_mesh.R')
)
if(!('--run' %in% args) && !('--stage' %in% args)) {
  for(p in plan) cat(sprintf('%-23s workers=%d  seed=%s  %s\n',p$name,p$workers,p$seed,p$script))
  cat('\nUse --run for all stages, or --stage NAME for one stage. Run from this replication directory.\n')
} else {
  if(requireNamespace('highs',quietly=TRUE) && as.character(packageVersion('highs')) != '1.14.0.2')
    stop('Use highs 1.14.0.2 for the recorded native LP interface; run install_dependencies.R.')
  if('--stage' %in% args) {
    index <- match('--stage',args)
    stopifnot(index<length(args))
    wanted <- args[index+1L]
    plan <- Filter(function(p)identical(p$name,wanted),plan)
    if(length(plan)!=1L) stop('Unknown stage: ',wanted)
  }
  # Remove inherited smoke-test and checkpoint switches. Each stage below
  # explicitly sets its production values; empirical data-path settings stay.
  managed <- c('JBES_N_CORES','JBES_N_MC','JBES_N_MC_ALT','JBES_N_MC_SIZE',
    'JBES_N_CAL','JBES_N_SIZE','JBES_CV_N_REP','JBES_CV_N_REP_OPEN',
    'JBES_CV_N_GRID_PER_UNIT','JBES_CV_N_GRID_TRAIN_OPEN','JBES_CV_N_GRID_MON_OPEN',
    'JBES_CV_D_GRID','JBES_CV_SUBDIR','JBES_CV_CHECKPOINT')
  # Clear any additional inherited design or smoke-test overrides, including
  # M/T/gamma grids and plotting filters, while retaining explicit file paths.
  managed <- unique(c(managed,setdiff(
    grep('^JBES_',names(Sys.getenv()),value=TRUE),
    c('JBES_PROJECT_ROOT','JBES_GBPUSD_PATH','JBES_LP_DIAGNOSTICS'))))
  saved <- Sys.getenv(managed,unset=NA_character_)
  dir.create(file.path(root,'logs'),showWarnings=FALSE)
  oldwd <- getwd(); setwd(root)
  for(p in plan) {
    Sys.unsetenv(managed)
    settings <- c(JBES_PROJECT_ROOT=root,JBES_N_CORES=as.character(p$workers),
      JBES_N_MC='1000',JBES_N_MC_ALT='1000',JBES_N_MC_SIZE='1000',
      JBES_N_CAL='1000',JBES_N_SIZE='1000',JBES_CV_N_REP='10000',
      JBES_CV_N_REP_OPEN='10000',JBES_CV_N_GRID_PER_UNIT='400',
      JBES_CV_N_GRID_TRAIN_OPEN='400',JBES_CV_N_GRID_MON_OPEN='5000',
      p$settings)
    do.call(Sys.setenv,as.list(settings))
    log <- file.path(root,'logs',paste0('reproduce_',p$name,'.txt'))
    cat('Running',p$name,'; log:',log,'\n'); flush.console()
    status <- system2(file.path(R.home('bin'),'Rscript'),
      shQuote(file.path(root,p$script)),stdout=log,stderr=log)
    Sys.unsetenv(managed)
    restore <- saved[!is.na(saved)]
    if(length(restore)) do.call(Sys.setenv,as.list(restore))
    if(status!=0L) stop('Stage failed: ',p$name,'; see ',log)
  }
  setwd(oldwd)
}
