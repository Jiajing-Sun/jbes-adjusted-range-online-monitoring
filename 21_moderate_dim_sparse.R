# ============================================================
# 21_moderate_dim_sparse.R
# - Appendix experiment: moderate dimension + dense vs sparse shifts
# - Design: q in {5,8,12,20}, early/late breaks, abrupt changes
# - Reports a compact size-adjusted size check plus strict post-break detection metrics
# ============================================================

rm(list = ls())

AUTO_INSTALL <- TRUE
N_CORES <- max(1L, parallel::detectCores() - 2L)
SEED_BASE <- 11223

N_CAL <- 1000
N_SIZE <- 1000
N_MC <- 1000

METHODS <- c("RSMS", "SSMS", "CSMS")
ALPHA <- 0.05

# Design
Q_GRID <- c(5, 8, 12, 20)
m <- 500
T <- 2
GAMMA_GRID <- c(0, 0.15)

BREAK_FRAC <- c(0.1, 0.8)
CHANGE_TYPE <- "abrupt"
DELTA_GRID <- c(0.25, 0.5, 0.75)

SHIFT_TYPES <- list(
  dense = list(type = "dense", s = NA),
  sparse = list(type = "sparse", s = 1)
)

# Moderate-dim DGP parameters (as in appendix text)
phi <- 0.5
cross <- 0.1
rho <- 0.3


# --- bootstrap project root and source utilities ----------------------
bootstrap_find_root <- function() {
  env_root <- Sys.getenv("JBES_PROJECT_ROOT", unset = "")
  if (nzchar(env_root)) {
    return(normalizePath(env_root, winslash = "/", mustWork = FALSE))
  }

  candidates <- character(0)

  # 1) Rscript --file=
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  mi <- grep(file_arg, args)
  if (length(mi) > 0) {
    script_path <- sub(file_arg, "", args[mi[1]])
    candidates <- c(candidates, dirname(normalizePath(script_path, winslash = "/", mustWork = FALSE)))
  }

  # 2) source() / knitr (may be NA)
  ofile <- tryCatch(attr(sys.frame(1), "ofile"), error = function(e) NULL)
  if (!is.null(ofile) && nzchar(ofile)) {
    candidates <- c(candidates, dirname(normalizePath(ofile, winslash = "/", mustWork = FALSE)))
  }

  # 3) RStudio active document
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

  stop(paste0(
    "Could not locate project root.\\n",
    "Fix options: (1) setwd('<project_root>') before running; or\\n",
    "            (2) set Sys.setenv(JBES_PROJECT_ROOT='<project_root>') once."
  ))
}

ROOT <- bootstrap_find_root()
source(file.path(ROOT, "00_utils.R"))
jbes_set_root(ROOT)
msg("Project root: ", ROOT)

N_CORES <- env_int("JBES_N_CORES", N_CORES)
N_CAL <- env_int("JBES_N_CAL", N_CAL)
N_SIZE <- env_int("JBES_N_SIZE", N_SIZE)
N_MC <- env_int("JBES_N_MC", N_MC)
Q_GRID <- env_num_vec("JBES_Q_GRID", Q_GRID)
GAMMA_GRID <- env_num_vec("JBES_GAMMA_GRID", GAMMA_GRID)
BREAK_FRAC <- env_num_vec("JBES_BREAK_FRAC", BREAK_FRAC)
DELTA_GRID <- env_num_vec("JBES_DELTA_GRID", DELTA_GRID)

# Paths
RES_DIR <- file.path(ROOT, "results")
LOG_DIR <- file.path(ROOT, "logs")
if (!dir.exists(RES_DIR)) dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
ensure_pkgs(c("MASS", "Matrix"), auto_install = AUTO_INSTALL)

source(file.path(ROOT, "02_methods.R"))
source(file.path(ROOT, "01_critical_values.R"))
source(file.path(ROOT, "03_affine_extrema.R"))
source(file.path(ROOT, "05_affine_stopping.R"))
source(file.path(ROOT, "10_dgps_mean.R"))
source(file.path(ROOT, "30_construct_psi.R"))
source(file.path(ROOT, "40_calibration.R"))
source(file.path(ROOT, "08_mc_audit.R"))

CHECK_DIR <- file.path(RES_DIR, "moderate_checkpoints")
dir.create(CHECK_DIR, recursive = TRUE, showWarnings = FALSE)
# A checkpoint is valid only for the recorded design, generator, package and
# reference values. Stop on a mismatch instead of mixing numerical versions.
manifest <- list(
  schema = 2L, q = Q_GRID, m = m, T = T, gamma = GAMMA_GRID, break_fraction = BREAK_FRAC,
  phi = phi, cross = cross, rho = rho, shift_types = SHIFT_TYPES, alpha = ALPHA,
  delta = DELTA_GRID, n_cal = N_CAL, n_size = N_SIZE, n_eval = N_MC, workers = N_CORES,
  seed = SEED_BASE, package = as.character(packageVersion('aersn')),
  source = unname(tools::md5sum(file.path(ROOT, c(
    '10_dgps_mean.R',
    '02_methods.R', '07_aersn_backend.R', '03_affine_extrema.R', '05_affine_stopping.R'
  )))),
  reference = unname(tools::md5sum(file.path(RES_DIR, 'affine_finite_critical_values.csv')))
)
manifest_path <- file.path(CHECK_DIR, 'manifest.rds')
if (file.exists(manifest_path)) {
  if (!identical(readRDS(manifest_path), manifest))
    stop('Moderate-dimension checkpoints differ from this run; use a fresh results directory.')
} else {
  if (length(list.files(CHECK_DIR, pattern = '\\.rds$')))
    stop('Unverified moderate-dimension checkpoints found; archive them before running.')
  saveRDS(manifest, manifest_path)
}

checkpoint_design <- function(name, expr) {
  path <- file.path(CHECK_DIR, paste0(name, ".rds"))
  if (file.exists(path)) return(readRDS(path))
  value <- force(expr);saveRDS(value, path);value
}
msg("Moderate-dim + sparsity experiment starting...")

rows <- operational_rows <- event_rows <- pair_rows <- list()
cnt <- 0L

for (q in Q_GRID) {
  for (gamma in GAMMA_GRID) {
    mT <- as.integer(round(m * T))
    # --- calibrate critical values under H0 for this (q,m,T,gamma)
    tag_cal <- sprintf("moderate_q%d_m%d_T%d_g%.2f", q, m, T, gamma)
    msg("Calibrating c_hat under H0: ", tag_cal)

    gen_psi_null <- function(i) {
      X <- sim_var1_dgp(
        dgp = "DGP2", n = m + mT, m_train = m, q = q,
        burnin = 300, phi = phi, cross = cross, rho = rho
      )
      make_mean_psi(X, m, T)$psi
    }
    c_hat_df <- checkpoint_design(paste0(tag_cal, "_cal"), calibrate_critical_values(
      METHODS = METHODS, gen_psi = gen_psi_null,
      m = m, T = T, gamma = gamma,
      alpha = ALPHA, N_CAL = N_CAL,
      n_cores = N_CORES,
      seed = SEED_BASE + 5000L + q * 10L + as.integer(100 * gamma),
      tag = tag_cal
    ))

    # helper to pull method-specific c_hat
    c_hat_map <- setNames(c_hat_df$c_hat, c_hat_df$method)

    # Evaluate operational reference size on the same independent null paths.
    # These are separate thresholds from the DGP-matched empirical calibration.
    c_operational <- vapply(METHODS, function(method) switch(method,
      RSMS = get_critical_value_rsms(q, T, gamma, ALPHA),
      SSMS = get_critical_value_ssms(q, T, gamma, ALPHA),
      CSMS = get_critical_value_csms(q, T, gamma, ALPHA)
    ), numeric(1))

    # --- compact size check under the calibrated c_hat
    msg("Running compact size check: ", tag_cal)
    reps_size <- checkpoint_design(paste0(tag_cal, "_size"), par_lapply_safe(
      1:N_SIZE,
      function(i) {
        psi <- gen_psi_null(i)
        maxima <- sapply(METHODS, function(method) {
          if (method == "RSMS") rsms_maximum(psi, m, T, gamma) else
            max(compute_path(method, psi, m, T, gamma))
        })
        c(
          as.integer(maxima > c_hat_map[METHODS]),
          as.integer(maxima > c_operational[METHODS])
        )
      },
      n_cores = N_CORES,
      seed = SEED_BASE + 90000L + q * 10L + as.integer(100 * gamma)
    ))
    size_mat <- do.call(rbind, reps_size)
    colnames(size_mat) <- c(METHODS, paste0(METHODS, '_operational'))
    for (method in METHODS) {
      count <- sum(size_mat[, paste0(method, '_operational')]); probability <- count / N_SIZE
      operational_rows[[length(operational_rows) + 1L]] <- data.frame(
        q, m, T, gamma, method,
        alpha = ALPHA, critical = c_operational[[method]],
        rejection_count = count, n_eval = N_SIZE, size = probability,
        standard_error = sqrt(probability * (1 - probability) / N_SIZE)
      )
    }

    for (method in METHODS) {
      rows[[length(rows) + 1L]] <- data.frame(
        exp = "size_adj",
        q = q, m = m, T = T, gamma = gamma,
        break_frac = NA_real_, k_star = NA_real_,
        shift = NA_character_, Delta = NA_real_,
        method = method,
        size_adj = mean(size_mat[, method]),
        power_adj = NA_real_,
        PFA_pre_adj = NA_real_,
        DP_post_adj = NA_real_,
        CADD_adj = NA_real_,
        ARL_adj = NA_real_,
        ADD_adj = NA_real_,
        N = N_SIZE, N_CAL = N_CAL,
        stringsAsFactors = FALSE
      )
    }

    for (bf in BREAK_FRAC) {
      k_star <- as.integer(floor(bf * mT))
      k_star <- max(1L, min(k_star, mT - 1L))

      for (shift_name in names(SHIFT_TYPES)) {
        shift_spec <- SHIFT_TYPES[[shift_name]]
        v <- make_shift_vec(q, type = shift_spec$type, s = shift_spec$s)

        for (Delta in DELTA_GRID) {
          cnt <- cnt + 1L
          tag <- sprintf("moderate_q%d_%s_bf%.1f_D%.2f_g%.2f", q, shift_name, bf, Delta, gamma)
          msg("Running power: ", tag)

          one_rep <- function(i) {
            X <- sim_var1_dgp(
              dgp = "DGP2", n = m + mT, m_train = m, q = q,
              burnin = 300, phi = phi, cross = cross, rho = rho
            )
            Y <- apply_mean_change(X, m, T, k_star = k_star, Delta = Delta, shift_vec = v, type = CHANGE_TYPE)
            psi <- make_mean_psi(Y, m, T)$psi

            out <- sapply(METHODS, function(method) {
              if (method == "RSMS") rsms_first_crossing(psi, m, T, gamma, c_hat_map[[method]]) else
                first_exceed(compute_path(method, psi, m, T, gamma), c_hat_map[[method]])
            })
            out
          }

          reps <- checkpoint_design(tag, par_lapply_safe(1:N_MC, one_rep,
            n_cores = N_CORES, seed = SEED_BASE + 1000L * cnt
          ))

          stop_mat <- do.call(rbind, reps) # N_MC x length(METHODS)
          colnames(stop_mat) <- METHODS
          audit <- mc_event_audit(stop_mat, k_star, mT)
          metadata <- data.frame(q, m, T, gamma, break_fraction = bf, shift = shift_name, Delta)
          event_rows[[length(event_rows) + 1L]] <- cbind(
            metadata[rep(1, nrow(audit$events)), , drop = FALSE], audit$events
          )
          pair_rows[[length(pair_rows) + 1L]] <- cbind(
            metadata[rep(1, nrow(audit$paired)), , drop = FALSE], audit$paired
          )

          for (method in METHODS) {
            stops <- stop_mat[, method]
            metrics <- detection_metrics(stops, k_star, mT)

            rows[[length(rows) + 1L]] <- data.frame(
              exp = "power_delay",
              q = q, m = m, T = T, gamma = gamma,
              break_frac = bf, k_star = k_star,
              shift = shift_name, Delta = Delta,
              method = method,
              size_adj = NA_real_,
              power_adj = metrics["DP_horizon"],
              PFA_pre_adj = metrics["PFA_pre"],
              DP_post_adj = metrics["DP_post"],
              CADD_adj = metrics["CADD"],
              ARL_adj = metrics["ARL_excess"],
              ADD_adj = metrics["ADD_excess"],
              N = N_MC, N_CAL = N_CAL,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
}

df <- do.call(rbind, rows)
out_path <- file.path(RES_DIR, "moderate_dim_sparse_results.csv")
write.csv(df, out_path, row.names = FALSE)
write.csv(do.call(rbind, operational_rows), file.path(RES_DIR, 'moderate_operational_size.csv'), row.names = FALSE)
write.csv(do.call(rbind, event_rows), file.path(RES_DIR, 'moderate_events.csv'), row.names = FALSE)
write.csv(do.call(rbind, pair_rows), file.path(RES_DIR, 'moderate_paired.csv'), row.names = FALSE)

msg("Saved: ", out_path)
msg("DONE.")
