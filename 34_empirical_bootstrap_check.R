# ============================================================
# 34_empirical_bootstrap_check.R
# - Centered circular-block bootstrap calibration check
#   for the USD/GBP empirical application.
# ============================================================

rm(list = ls())

AUTO_INSTALL <- TRUE
ALPHA <- 0.05
N_BOOT <- 1000L
set.seed(20260913L)
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

N_BOOT <- env_int("JBES_EMP_BOOT", N_BOOT)
M_VALUES_2016 <- env_num_vec("JBES_EMP_M_2016", M_VALUES_2016)
M_VALUES_2020 <- env_num_vec("JBES_EMP_M_2020", M_VALUES_2020)
T_VALUES <- env_num_vec("JBES_EMP_T", T_VALUES)
GAMMA_VALUES <- env_num_vec("JBES_EMP_GAMMA", GAMMA_VALUES)
FVE_THRESHOLD <- env_num("JBES_EMP_FVE", FVE_THRESHOLD)
MAX_DIM <- env_int("JBES_EMP_MAX_DIM", MAX_DIM)

sample_circular_blocks <- function(train_scores, total_len, block_len) {
  m <- nrow(train_scores)
  starts <- sample.int(m, ceiling(total_len / block_len), replace = TRUE)
  idx <- integer(0)
  for (st in starts) {
    idx <- c(idx, ((st - 1L + seq_len(block_len) - 1L) %% m) + 1L)
  }
  train_scores[idx[seq_len(total_len)], , drop = FALSE]
}

bootstrap_cvalue <- function(score_mat, m, T, gamma, method, n_boot, alpha) {
  total_len <- as.integer(m + m * T)
  train <- score_mat[seq_len(m), , drop = FALSE]
  train_centered <- scale(train, center = TRUE, scale = FALSE)
  block_len <- max(2L, ceiling(1.5 * m^(1 / 3)))
  boot_max <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    boot_scores <- sample_circular_blocks(train_centered, total_len, block_len)
    # Refit the bootstrap training mean, as in the original estimator.
    boot_scores <- sweep(boot_scores, 2L, colMeans(boot_scores[seq_len(m),,drop=FALSE]), "-")
    boot_max[b] <- max(compute_path(method, boot_scores, m, T, gamma))
  }
  as.numeric(stats::quantile(boot_max, probs = 1 - alpha, type = 8, names = FALSE))
}

raw_df <- load_gbpusd_data(data_path)
event_specs <- list(
  list(year = 2016L, event = "Brexit", event_date = as.Date("2016-06-23"), m_values = M_VALUES_2016),
  list(year = 2020L, event = "COVID", event_date = as.Date("2020-03-17"), m_values = M_VALUES_2020)
)

rows <- list()

for (spec in event_specs) {
  msg("Building daily curves for year ", spec$year)
  year_df <- prepare_gbpusd_year(raw_df, spec$year)
  fd_bundle <- build_daily_fd(year_df)
  eval_bundle <- evaluate_fd_matrix(fd_bundle)

  cfg <- feasible_empirical_configs(
    dates = fd_bundle$dates,
    event_date = spec$event_date,
    m_values = spec$m_values,
    T_values = T_VALUES,
    gamma_values = GAMMA_VALUES
  )

  for (i in seq_len(nrow(cfg))) {
    m <- cfg$m[i]
    T <- cfg$T[i]
    gamma <- cfg$gamma[i]
    score_bundle <- compute_empirical_score_matrix(
      value_matrix = eval_bundle$values,
      dates = fd_bundle$dates,
      m = m,
      T = T,
      fve_threshold = FVE_THRESHOLD,
      max_dim = MAX_DIM
    )
    asym <- run_empirical_monitoring(
      score_mat = score_bundle$scores,
      dates = score_bundle$dates,
      m = m,
      T = T,
      gamma = gamma,
      alpha = ALPHA
    )
    for (method in c("RSMS", "SSMS", "CSMS")) {
      msg(sprintf("Bootstrap: year=%d, m=%d, T=%d, gamma=%.2f, method=%s", spec$year, m, T, gamma, method))
      stat <- asym$statistic[asym$method == method]
      asym_reject <- asym$reject[asym$method == method]
      boot_cv <- bootstrap_cvalue(score_bundle$scores, m, T, gamma, method, N_BOOT, ALPHA)
      rows[[length(rows) + 1L]] <- data.frame(
        year = spec$year,
        event = spec$event,
        m = m,
        T = T,
        gamma = gamma,
        method = method,
        statistic = stat,
        asymptotic_reject = asym_reject,
        bootstrap_critical_value = boot_cv,
        bootstrap_reject = stat > boot_cv,
        n_boot = N_BOOT,
        stringsAsFactors = FALSE
      )
    }
  }
}

boot_df <- do.call(rbind, rows)
out_csv <- file.path(RES_DIR, "empirical_bootstrap_check.csv")
write.csv(boot_df, out_csv, row.names = FALSE)

summary_df <- boot_df |>
  dplyr::group_by(year, method) |>
  dplyr::summarise(
    configurations = dplyr::n(),
    asymptotic_rejects = sum(asymptotic_reject),
    bootstrap_rejects = sum(bootstrap_reject),
    disagreements = sum(asymptotic_reject != bootstrap_reject),
    .groups = "drop"
  )
out_summary <- file.path(RES_DIR, "empirical_bootstrap_check_summary.csv")
write.csv(summary_df, out_summary, row.names = FALSE)

msg("Saved: ", out_csv)
msg("Saved: ", out_summary)
msg("DONE.")
