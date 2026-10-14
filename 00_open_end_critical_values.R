# ============================================================
# 00_open_end_critical_values.R
# - Simulate open-end (T = Inf) asymptotic critical values.
#
# Purpose:
#   Use the same training-grid discretization as the finite-horizon
#   critical-value script (N_GRID_TRAIN = 400), while simulating the
#   open-end monitoring Brownian motion on u = s/(1+s) with a finer
#   grid (N_GRID_MON = 5000).
#
# Output:
#   results/affine_open_end_critical_values.csv
#   results/open_end_critical_values_table_rows.tex
#   results/open_end_critical_values_vs_finite_T10.csv
# ============================================================

rm(list = ls())

AUTO_INSTALL <- TRUE
SEED_BASE <- 975246L
N_REP <- 10000L
N_GRID_TRAIN <- 400L
N_GRID_MON <- 5000L
N_CORES <- min(8L, max(1L, parallel::detectCores() - 2L))

D_GRID <- 1:20
GAMMA_GRID <- c(0, 0.15)
ALPHA_GRID <- c(0.05, 0.10)
METHODS <- c("RSMS", "CSMS", "SSMS")

bootstrap_find_root <- function() {
  env_root <- Sys.getenv("JBES_PROJECT_ROOT", unset = "")
  if (nzchar(env_root)) {
    return(normalizePath(env_root, winslash = "/", mustWork = FALSE))
  }

  candidates <- character(0)
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  mi <- grep(file_arg, args)
  if (length(mi) > 0) {
    script_path <- sub(file_arg, "", args[mi[1]])
    candidates <- c(candidates, dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE)))
  }

  ofile <- tryCatch(attr(sys.frame(1), "ofile"), error = function(e) NULL)
  if (!is.null(ofile) && nzchar(ofile)) {
    candidates <- c(candidates, dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE)))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    ctx_path <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(ctx_path)) {
      candidates <- c(candidates, dirname(normalizePath(ctx_path, winslash = "/", mustWork = FALSE)))
    }
  }

  candidates <- unique(c(candidates, getwd()))
  for (st in candidates) {
    st <- normalizePath(st, winslash = "/", mustWork = FALSE)
    cur <- st
    for (i in 0:25) {
      if (file.exists(file.path(cur, "00_utils.R"))) {
        return(cur)
      }
      parent <- dirname(cur)
      if (parent == cur) break
      cur <- parent
    }
  }

  stop(
    "Could not locate project root. Set setwd('<project_root>') or ",
    "Sys.setenv(JBES_PROJECT_ROOT='<project_root>')."
  )
}

ROOT <- bootstrap_find_root()
source(file.path(ROOT, "00_utils.R"))
source(file.path(ROOT, "02_methods.R"))
source(file.path(ROOT, "03_affine_extrema.R"))
jbes_set_root(ROOT)

N_REP <- env_int("JBES_CV_N_REP_OPEN", N_REP)
SEED_BASE <- env_int("JBES_CV_SEED_OPEN", SEED_BASE)
N_GRID_TRAIN <- env_int("JBES_CV_N_GRID_TRAIN_OPEN", N_GRID_TRAIN)
N_GRID_MON <- env_int("JBES_CV_N_GRID_MON_OPEN", N_GRID_MON)
N_CORES <- env_int("JBES_N_CORES", N_CORES)
D_GRID <- env_num_vec("JBES_CV_D_GRID", D_GRID)

RES_DIR <- file.path(ROOT, "results", Sys.getenv("JBES_CV_SUBDIR", unset=""))
if (!dir.exists(RES_DIR)) dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

ensure_pkgs("parallel", auto_install = AUTO_INSTALL)

msg("Project root: ", ROOT)
msg(
  "Simulating open-end critical values with ",
  "N_REP=", N_REP,
  ", N_GRID_TRAIN=", N_GRID_TRAIN,
  ", N_GRID_MON=", N_GRID_MON,
  ", N_CORES=", N_CORES
)

simulate_open_end_one <- function(rep_id) {
  d_max <- max(D_GRID)

  u_train <- seq_len(N_GRID_TRAIN) / N_GRID_TRAIN
  dt_train <- 1 / N_GRID_TRAIN
  dW_train <- matrix(rnorm(N_GRID_TRAIN * d_max, sd = sqrt(dt_train)),
                     nrow = N_GRID_TRAIN, ncol = d_max)
  W_train <- apply(dW_train, 2, cumsum)
  B_train <- W_train - u_train %o% W_train[N_GRID_TRAIN, ]

  V_train <- crossprod(B_train) * dt_train

  u_mon <- seq_len(N_GRID_MON) / N_GRID_MON
  dt_mon <- 1 / N_GRID_MON
  dW_mon <- matrix(rnorm(N_GRID_MON * d_max, sd = sqrt(dt_mon)),
                   nrow = N_GRID_MON, ncol = d_max)
  W_mon <- apply(dW_mon, 2, cumsum)

  mon_sq_cum <- t(apply(W_mon^2, 1, cumsum))

  out_dim <- c(length(D_GRID), length(GAMMA_GRID))
  out_names <- list(paste0("d", D_GRID), paste0("g", GAMMA_GRID))
  csms_stat <- matrix(NA_real_, nrow = out_dim[1], ncol = out_dim[2], dimnames = out_names)
  rsms_stat <- csms_stat
  ssms_stat <- csms_stat

  gauge_weights <- vapply(GAMMA_GRID, function(g) u_mon^(-g), numeric(N_GRID_MON))
  for (di in seq_along(D_GRID)) {
    d <- D_GRID[di]
    rsms_stat[di,] <- ae_gauge_maxima(rbind(0, B_train[,1:d,drop=FALSE]),
      W_mon[,1:d,drop=FALSE], N_GRID_MON, gauge_weights)^2
  }

  ssms_quad_list <- vector("list", length(D_GRID))
  for (d_idx in seq_along(D_GRID)) {
    d <- D_GRID[d_idx]
    V_inv <- inv_pd(V_train[1:d, 1:d, drop = FALSE], ridge = 1e-10)
    W_d <- W_mon[, 1:d, drop = FALSE]
    ssms_quad_list[[d_idx]] <- rowSums((W_d %*% V_inv) * W_d)
  }

  for (g_idx in seq_along(GAMMA_GRID)) {
    weight <- u_mon^(-2 * GAMMA_GRID[g_idx])
    csms_stat[, g_idx] <- apply(sweep(mon_sq_cum[, D_GRID, drop = FALSE], 1, weight, `*`), 2, max)
    ssms_stat[, g_idx] <- vapply(ssms_quad_list, function(x) max(x * weight), numeric(1))
  }

  list(RSMS = rsms_stat, CSMS = csms_stat, SSMS = ssms_stat)
}

if (identical(Sys.getenv("JBES_CV_CHECKPOINT"), "1")) {
 source(file.path(ROOT,"04_reference_replicates.R"))
 reps <- reference_replicates(simulate_open_end_one, N_REP, N_CORES, SEED_BASE,
  file.path(RES_DIR,"checkpoints_open"), list(mode="open",d=D_GRID,gamma=GAMMA_GRID,grid_train=N_GRID_TRAIN,grid_mon=N_GRID_MON))
} else {
reps <- par_lapply_safe(
  X = seq_len(N_REP),
  fun = simulate_open_end_one,
  n_cores = N_CORES,
  seed = SEED_BASE,
  packages = "parallel"
)
}

collect_method_matrix <- function(method) {
  out <- array(NA_real_, dim = c(N_REP, length(D_GRID), length(GAMMA_GRID)))
  for (i in seq_len(N_REP)) out[i, , ] <- reps[[i]][[method]]
  out
}

rows <- list()
for (method in METHODS) {
  arr <- collect_method_matrix(method)
  for (d_idx in seq_along(D_GRID)) {
    for (g_idx in seq_along(GAMMA_GRID)) {
      vals <- arr[, d_idx, g_idx]
      for (alpha in ALPHA_GRID) {
        rows[[length(rows) + 1L]] <- data.frame(
          method = method,
          d = D_GRID[d_idx],
          gamma = GAMMA_GRID[g_idx],
          alpha = alpha,
          critical_value = as.numeric(stats::quantile(vals, probs = 1 - alpha, names = FALSE, type = 8)),
          n_rep = N_REP,
          n_grid_train = N_GRID_TRAIN,
          n_grid_monitoring = N_GRID_MON,
          seed_base = SEED_BASE,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

cv_df <- do.call(rbind, rows)
cv_df <- cv_df[order(match(cv_df$method, METHODS), cv_df$d, cv_df$gamma, cv_df$alpha), ]

csv_path <- file.path(RES_DIR, "affine_open_end_critical_values.csv")
write.csv(cv_df, csv_path, row.names = FALSE)
msg("Saved: ", csv_path)

fmt_cv <- function(method, d, gamma, alpha) {
  x <- cv_df$critical_value[
    cv_df$method == method &
      cv_df$d == d &
      abs(cv_df$gamma - gamma) < 1e-12 &
      abs(cv_df$alpha - alpha) < 1e-12
  ]
  if (length(x) != 1L) stop("Lookup failure for ", method, ", d=", d)
  formatC(x, format = "f", digits = 1)
}

panel_lines <- function(method, label) {
  lines <- c(
    sprintf("\\multicolumn{5}{l}{\\textbf{%s}} \\\\", label),
    "\\midrule"
  )
  for (d in D_GRID) {
    lines <- c(
      lines,
      sprintf(
        "%d & %s & %s & %s & %s \\\\",
        d,
        fmt_cv(method, d, 0, 0.05),
        fmt_cv(method, d, 0, 0.10),
        fmt_cv(method, d, 0.15, 0.05),
        fmt_cv(method, d, 0.15, 0.10)
      )
    )
  }
  lines
}

tex_lines <- c(
  panel_lines("RSMS", "Panel A: RSMS"),
  "\\addlinespace[0.4em]",
  panel_lines("CSMS", "Panel B: CSMS"),
  "\\addlinespace[0.4em]",
  panel_lines("SSMS", "Panel C: SSMS")
)

tex_path <- file.path(RES_DIR, "open_end_critical_values_table_rows.tex")
writeLines(tex_lines, tex_path)
msg("Saved: ", tex_path)

finite_path <- file.path(RES_DIR, "affine_finite_critical_values.csv")
if (file.exists(finite_path)) {
  finite_df <- read.csv(finite_path, stringsAsFactors = FALSE)
  finite_t10 <- finite_df[finite_df$T == 10, c("method", "d", "gamma", "alpha", "critical_value")]
  names(finite_t10)[names(finite_t10) == "critical_value"] <- "finite_T10"

  cmp <- merge(
    finite_t10,
    cv_df[, c("method", "d", "gamma", "alpha", "critical_value")],
    by = c("method", "d", "gamma", "alpha"),
    all = TRUE
  )
  names(cmp)[names(cmp) == "critical_value"] <- "open_end"
  cmp$difference_open_minus_T10 <- cmp$open_end - cmp$finite_T10
  cmp_path <- file.path(RES_DIR, "open_end_critical_values_vs_finite_T10.csv")
  write.csv(cmp, cmp_path, row.names = FALSE)
  msg("Saved: ", cmp_path)
}

msg("DONE.")
