# ============================================================
# 32_empirical_monitoring_table.R
# - Clean empirical monitoring outputs for 2016 Brexit and 2020 COVID
# - Produces a tidy table with statistic / critical value / reject / first date
# ============================================================

rm(list = ls())

AUTO_INSTALL <- TRUE
ALPHA <- 0.05
M_VALUES_2016 <- c(50, 75, 100)
M_VALUES_2020 <- c(50)
T_VALUES <- c(1, 2, 5)
GAMMA_VALUES <- c(0, 0.15)
FVE_THRESHOLD <- 0.95
MAX_DIM <- 8L

bootstrap_find_root <- function() {
  env_root <- Sys.getenv("JBES_PROJECT_ROOT", unset = "")
  if (nzchar(env_root)) {
    return(normalizePath(env_root, winslash = "/", mustWork = FALSE))
  }
  candidates <- unique(c(getwd()))
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    ctx_path <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(ctx_path)) candidates <- unique(c(dirname(ctx_path), candidates))
  }
  for (st in candidates) {
    st <- normalizePath(st, winslash = "/", mustWork = FALSE)
    cur <- st
    for (i in 0:25) {
      if (file.exists(file.path(cur, "00_utils.R"))) return(cur)
      parent <- dirname(cur)
      if (parent == cur) break
      cur <- parent
    }
  }
  stop("Could not locate project root.")
}

ROOT <- bootstrap_find_root()
source(file.path(ROOT, "00_utils.R"))
source(file.path(ROOT, "01_critical_values.R"))
source(file.path(ROOT, "02_methods.R"))
source(file.path(ROOT, "50_empirical_utils.R"))
jbes_set_root(ROOT)
msg("Project root: ", ROOT)

ensure_pkgs(c("dplyr", "readr", "lubridate", "fda"), auto_install = AUTO_INSTALL)

data_path <- env_string(
  "JBES_GBPUSD_PATH",
  file.path(ROOT, "data", "GBPUSD.csv")
)

RES_DIR <- file.path(ROOT, "results")
if (!dir.exists(RES_DIR)) dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

M_VALUES_2016 <- env_num_vec("JBES_EMP_M_2016", M_VALUES_2016)
M_VALUES_2020 <- env_num_vec("JBES_EMP_M_2020", M_VALUES_2020)
T_VALUES <- env_num_vec("JBES_EMP_T", T_VALUES)
GAMMA_VALUES <- env_num_vec("JBES_EMP_GAMMA", GAMMA_VALUES)
FVE_THRESHOLD <- env_num("JBES_EMP_FVE", FVE_THRESHOLD)
MAX_DIM <- env_int("JBES_EMP_MAX_DIM", MAX_DIM)

raw_df <- load_gbpusd_data(data_path)

event_specs <- list(
  list(year = 2016L, event = "Brexit", event_date = as.Date("2016-06-23"), m_values = M_VALUES_2016),
  list(year = 2020L, event = "COVID", event_date = as.Date("2020-03-17"), m_values = M_VALUES_2020)
)

rows <- list()
curve_meta <- list()

for (spec in event_specs) {
  msg("Building daily curves for year ", spec$year)
  year_df <- prepare_gbpusd_year(raw_df, spec$year)
  fd_bundle <- build_daily_fd(year_df)
  eval_bundle <- evaluate_fd_matrix(fd_bundle)

  curve_meta[[length(curve_meta) + 1L]] <- data.frame(
    year = spec$year,
    event = spec$event,
    valid_days = fd_bundle$valid_days,
    nbasis = fd_bundle$nbasis,
    min_obs = fd_bundle$min_obs,
    max_obs = fd_bundle$max_obs,
    stringsAsFactors = FALSE
  )

  cfg <- feasible_empirical_configs(
    dates = fd_bundle$dates,
    event_date = spec$event_date,
    m_values = spec$m_values,
    T_values = T_VALUES,
    gamma_values = GAMMA_VALUES
  )

  if (nrow(cfg) == 0L) {
    warning("No feasible monitoring configurations for year ", spec$year)
    next
  }

  for (i in seq_len(nrow(cfg))) {
    m <- cfg$m[i]
    T <- cfg$T[i]
    gamma <- cfg$gamma[i]
    msg(sprintf("Empirical monitoring: year=%d, m=%d, T=%d, gamma=%.2f", spec$year, m, T, gamma))

    score_bundle <- compute_empirical_score_matrix(
      value_matrix = eval_bundle$values,
      dates = fd_bundle$dates,
      m = m,
      T = T,
      fve_threshold = FVE_THRESHOLD,
      max_dim = MAX_DIM
    )

    res <- run_empirical_monitoring(
      score_mat = score_bundle$scores,
      dates = score_bundle$dates,
      m = m,
      T = T,
      gamma = gamma,
      alpha = ALPHA
    )

    res$year <- spec$year
    res$event <- spec$event
    res$event_date <- as.character(spec$event_date)
    res$pre_event_detection <- ifelse(
      is.na(res$first_detection_date),
      NA,
      as.Date(res$first_detection_date) < spec$event_date
    )
    res$valid_days <- fd_bundle$valid_days
    rows[[length(rows) + 1L]] <- res
  }
}

res_df <- do.call(rbind, rows)
res_df <- res_df[, c(
  "year", "event", "event_date", "valid_days",
  "m", "T", "gamma", "method", "d_eff",
  "statistic", "critical_value", "reject", "stop_k", "first_detection_date", "pre_event_detection"
)]
res_df <- res_df[order(res_df$year, res_df$m, res_df$T, res_df$gamma, res_df$method), ]

meta_df <- do.call(rbind, curve_meta)

csv_path <- file.path(RES_DIR, "empirical_monitoring_results.csv")
reject_path <- file.path(RES_DIR, "empirical_monitoring_results_reject_only.csv")
meta_path <- file.path(RES_DIR, "empirical_curve_metadata.csv")
write.csv(res_df, csv_path, row.names = FALSE)
write.csv(subset(res_df, reject), reject_path, row.names = FALSE)
write.csv(meta_df, meta_path, row.names = FALSE)

msg("Saved: ", csv_path)
msg("Saved: ", reject_path)
msg("Saved: ", meta_path)
msg("DONE.")
