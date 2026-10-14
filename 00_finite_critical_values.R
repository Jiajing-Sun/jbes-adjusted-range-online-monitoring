# ============================================================
# 00_finite_critical_values.R
# - Simulate finite-horizon asymptotic critical values for
#   RSMS / CSMS / SSMS under T = 1,2,5,10
# - Outputs CSV for lookup and manuscript tables
# ============================================================

rm(list = ls())

AUTO_INSTALL <- TRUE
SEED_BASE <- 13579
N_REP <- 10000L
N_GRID_PER_UNIT <- 400L
N_CORES <- max(1L, parallel::detectCores() - 2L)

D_GRID <- 1:20
T_GRID <- c(1, 2, 5, 10)
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

RES_DIR <- file.path(ROOT, "results", Sys.getenv("JBES_CV_SUBDIR", unset=""))
if (!dir.exists(RES_DIR)) dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

ensure_pkgs("parallel", auto_install = AUTO_INSTALL)

N_REP <- env_int("JBES_CV_N_REP", N_REP)
N_GRID_PER_UNIT <- env_int("JBES_CV_N_GRID_PER_UNIT", N_GRID_PER_UNIT)
N_CORES <- env_int("JBES_N_CORES", N_CORES)
D_GRID <- env_num_vec("JBES_CV_D_GRID", D_GRID)
T_GRID <- env_num_vec("JBES_CV_T_GRID", T_GRID)

msg("Project root: ", ROOT)
msg(
  "Simulating finite-horizon asymptotic critical values with ",
  "N_REP=", N_REP,
  ", N_GRID_PER_UNIT=", N_GRID_PER_UNIT,
  ", N_CORES=", N_CORES
)

simulate_finite_one <- function(rep_id) {
  d_max <- max(D_GRID)
  T_max <- max(T_GRID)
  dt <- 1 / N_GRID_PER_UNIT
  total_steps <- as.integer((T_max + 1) * N_GRID_PER_UNIT)

  dW <- matrix(rnorm(total_steps * d_max, sd = sqrt(dt)),
               nrow = total_steps, ncol = d_max)
  W <- rbind(0, apply(dW, 2, cumsum))

  idx_t1 <- N_GRID_PER_UNIT + 1L
  W1 <- W[idx_t1, ]

  train_times <- seq_len(N_GRID_PER_UNIT) / N_GRID_PER_UNIT
  B_train <- W[2:(N_GRID_PER_UNIT + 1L), , drop = FALSE] - train_times %o% W1

  V_full <- crossprod(B_train) * dt
  V_inv_list <- lapply(D_GRID, function(d) inv_pd(V_full[1:d, 1:d, drop = FALSE], ridge = 1e-10))

  arr_dims <- c(length(D_GRID), length(T_GRID), length(GAMMA_GRID))
  arr_names <- list(paste0("d", D_GRID), paste0("T", T_GRID), paste0("g", GAMMA_GRID))
  rsms_stat <- array(NA_real_, dim = arr_dims, dimnames = arr_names)
  csms_stat <- array(NA_real_, dim = arr_dims, dimnames = arr_names)
  ssms_stat <- array(NA_real_, dim = arr_dims, dimnames = arr_names)

  n_full <- as.integer(T_max * N_GRID_PER_UNIT)
  s_full <- seq_len(n_full) / N_GRID_PER_UNIT
  U_full <- W[(idx_t1 + 1L):(idx_t1 + n_full), , drop = FALSE] - (1 + s_full) %o% W1
  gauge_weights <- vapply(GAMMA_GRID, function(g) 1 / ((1 + s_full) * (s_full / (1 + s_full))^g), numeric(n_full))
  for (di in seq_along(D_GRID)) {
    d <- D_GRID[di]
    rsms_stat[di,,] <- ae_gauge_maxima(rbind(0, B_train[,1:d,drop=FALSE]),
      U_full[,1:d,drop=FALSE], as.integer(T_GRID * N_GRID_PER_UNIT), gauge_weights)^2
  }

  for (t_idx in seq_along(T_GRID)) {
    T0 <- T_GRID[t_idx]
    n_mon <- as.integer(T0 * N_GRID_PER_UNIT)
    s <- seq_len(n_mon) / N_GRID_PER_UNIT
    idx_mon <- (idx_t1 + 1L):(idx_t1 + n_mon)
    U <- W[idx_mon, , drop = FALSE] - (1 + s) %o% W1

    U2 <- U^2
    csms_cum <- t(apply(U2, 1, cumsum))

    ssms_quad_list <- vector("list", length(D_GRID))
    for (d_idx in seq_along(D_GRID)) {
      d <- D_GRID[d_idx]
      U_d <- U[, 1:d, drop = FALSE]
      ssms_quad_list[[d_idx]] <- rowSums((U_d %*% V_inv_list[[d_idx]]) * U_d)
    }

    for (g_idx in seq_along(GAMMA_GRID)) {
      gamma0 <- GAMMA_GRID[g_idx]
      weight <- 1 / ((1 + s)^2 * (s / (1 + s))^(2 * gamma0))
      csms_stat[, t_idx, g_idx] <- apply(sweep(csms_cum[, D_GRID, drop = FALSE], 1, weight, `*`), 2, max)
      ssms_stat[, t_idx, g_idx] <- vapply(ssms_quad_list, function(x) max(x * weight), numeric(1))
    }
  }

  list(RSMS = rsms_stat, CSMS = csms_stat, SSMS = ssms_stat)
}

if (identical(Sys.getenv("JBES_CV_CHECKPOINT"), "1")) {
 source(file.path(ROOT,"04_reference_replicates.R"))
 reps <- reference_replicates(simulate_finite_one, N_REP, N_CORES, SEED_BASE,
  file.path(RES_DIR,"checkpoints_finite"), list(mode="finite",d=D_GRID,T=T_GRID,gamma=GAMMA_GRID,grid=N_GRID_PER_UNIT))
} else {
reps <- par_lapply_safe(
  X = seq_len(N_REP),
  fun = simulate_finite_one,
  n_cores = N_CORES,
  seed = SEED_BASE,
  packages = "parallel"
)
}

collect_method_array <- function(method) {
  out <- array(NA_real_, dim = c(N_REP, length(D_GRID), length(T_GRID), length(GAMMA_GRID)))
  for (i in seq_len(N_REP)) out[i, , , ] <- reps[[i]][[method]]
  out
}

rows <- list()
for (method in METHODS) {
  arr <- collect_method_array(method)
  for (d_idx in seq_along(D_GRID)) {
    for (t_idx in seq_along(T_GRID)) {
      for (g_idx in seq_along(GAMMA_GRID)) {
        vals <- arr[, d_idx, t_idx, g_idx]
        for (alpha in ALPHA_GRID) {
          rows[[length(rows) + 1L]] <- data.frame(
            method = method,
            d = D_GRID[d_idx],
            T = T_GRID[t_idx],
            gamma = GAMMA_GRID[g_idx],
            alpha = alpha,
            critical_value = as.numeric(stats::quantile(vals, probs = 1 - alpha, names = FALSE)),
            n_rep = N_REP,
            n_grid_per_unit = N_GRID_PER_UNIT,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}

cv_df <- do.call(rbind, rows)
cv_df <- cv_df[order(match(cv_df$method, METHODS), cv_df$d, cv_df$T, cv_df$gamma, cv_df$alpha), ]

csv_path <- file.path(RES_DIR, "affine_finite_critical_values.csv")
write.csv(cv_df, csv_path, row.names = FALSE)
msg("Saved: ", csv_path)
msg("DONE.")
