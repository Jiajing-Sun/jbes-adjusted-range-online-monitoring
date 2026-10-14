# ============================================================
# 02_mean_alternatives_power_delay.R
# - Alternative simulations for mean monitoring
# - Reports raw horizon rejection and strict post-break detection metrics
#   using asymptotic and size-adjusted critical values.
# ============================================================

rm(list = ls())

AUTO_INSTALL <- TRUE
N_CORES <- 1L
SEED_BASE <- 24680

N_MC_ALT <- 1000

# Design (baseline q=2)
DGP_GRID <- c("DGP1", "DGP2", "DGP3", "DGP4")
M_GRID <- c(100, 250, 500)
T_GRID <- c(1, 2, 5) # alternatives: focus on shorter horizons
GAMMA_GRID <- c(0, 0.15)

BREAK_FRAC <- c(0.1, 0.8) # early vs late (close to horizon end)
CHANGE_TYPES <- c("abrupt", "smooth")
EARLY_DELTA_GRID <- c(0.25, 0.50, 0.75, 1.00)
LATE_DELTA_GRID <- c(0.50, 1.00, 1.50, 2.00)

METHODS <- c("RSMS", "SSMS", "CSMS")
ALPHA <- 0.05

# Paths

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
source(file.path(ROOT, "08_mc_audit.R"))
msg("Project root: ", ROOT)

N_CORES <- env_int("JBES_N_CORES", N_CORES)
N_MC_ALT <- env_int("JBES_N_MC_ALT", N_MC_ALT)
DGP_GRID <- env_chr_vec("JBES_DGP_GRID", DGP_GRID)
M_GRID <- env_num_vec("JBES_M_GRID", M_GRID)
T_GRID <- env_num_vec("JBES_T_GRID", T_GRID)
GAMMA_GRID <- env_num_vec("JBES_GAMMA_GRID", GAMMA_GRID)
BREAK_FRAC <- env_num_vec("JBES_BREAK_FRAC", BREAK_FRAC)
CHANGE_TYPES <- env_chr_vec("JBES_CHANGE_TYPES", CHANGE_TYPES)
EARLY_DELTA_GRID <- env_num_vec("JBES_DELTA_GRID_EARLY", EARLY_DELTA_GRID)
LATE_DELTA_GRID <- env_num_vec("JBES_DELTA_GRID_LATE", LATE_DELTA_GRID)

# Paths
RES_DIR <- file.path(ROOT, "results")
LOG_DIR <- file.path(ROOT, "logs")
if (!dir.exists(RES_DIR)) dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
ensure_pkgs(c("MASS", "Matrix"), auto_install = AUTO_INSTALL)
source(file.path(ROOT, "01_critical_values.R"))
source(file.path(ROOT, "02_methods.R"))
source(file.path(ROOT, "10_dgps_mean.R"))
source(file.path(ROOT, "30_construct_psi.R"))

# Load calibration critical values
calib_path <- file.path(RES_DIR, "mean_null_calibration_cHat.csv")
if (!file.exists(calib_path)) stop("Calibration file not found: ", calib_path)

calib_df <- read.csv(calib_path, stringsAsFactors = FALSE)
calib_df <- subset(calib_df, abs(alpha - ALPHA) < 1e-12)

get_c_hat <- function(dgp0, m0, T0, gamma0, method0) {
  sub <- calib_df[calib_df$dgp == dgp0 &
    calib_df$m == m0 &
    calib_df$T == T0 &
    abs(calib_df$gamma - gamma0) < 1e-12 &
    calib_df$method == method0, ]
  if (nrow(sub) != 1) {
    stop(
      "Missing or duplicate calibration for: ",
      paste(dgp0, m0, T0, gamma0, method0, sep = " / ")
    )
  }
  sub$c_hat[1]
}


msg("Starting ALT simulations (mean monitoring)")

perf_rows <- list()
cnt <- 0L

for (dgp in DGP_GRID) {
  for (m in M_GRID) {
    for (T in T_GRID) {
      for (gamma in GAMMA_GRID) {
        mT <- as.integer(round(m * T))
        # Thresholds are fixed within a design. Resolve them once here,
        # avoiding repeated CSV reads inside every Monte Carlo replication.
        c_asym_map <- c(
          RSMS = get_critical_value_rsms(2, T, gamma, ALPHA),
          SSMS = get_critical_value_ssms(2, T, gamma, ALPHA),
          CSMS = get_critical_value_csms(2, T, gamma, ALPHA)
        )
        c_hat_map <- setNames(vapply(METHODS, function(method)
          get_c_hat(dgp, m, T, gamma, method), numeric(1)), METHODS)

        for (bf in BREAK_FRAC) {
          k_star <- as.integer(floor(bf * mT))
          k_star <- max(1L, min(k_star, mT - 1L))
          delta_grid <- if (bf < 0.5) EARLY_DELTA_GRID else LATE_DELTA_GRID

          for (chg in CHANGE_TYPES) {
            for (Delta in delta_grid) {
              cnt <- cnt + 1L
              tag <- sprintf(
                "mean_alt_%s_m%d_T%d_g%.2f_bf%.1f_%s_D%.2f",
                dgp, m, T, gamma, bf, chg, Delta
              )
              msg("Design: ", tag)

              # One replication under H1
              one_rep <- function(i) {
                X <- sim_var1_dgp(dgp = dgp, n = m + mT, m_train = m, q = 2, burnin = 250)
                shift_vec <- rep(1, 2) # matches main text
                Y <- apply_mean_change(X, m, T,
                  k_star = k_star, Delta = Delta,
                  shift_vec = shift_vec, type = chg
                )
                psi <- make_mean_psi(Y, m, T)$psi

                # store per method
                out <- list()
                for (method in METHODS) {
                  Mvec <- compute_path(method, psi, m, T, gamma)
                  max_stat <- max(Mvec)

                  c_asym <- c_asym_map[[method]]
                  stop_asym <- first_exceed(Mvec, c_asym)

                  c_hat <- c_hat_map[[method]]
                  stop_adj <- first_exceed(Mvec, c_hat)

                  out[[method]] <- c(
                    max_stat = max_stat,
                    stop_asym = stop_asym,
                    stop_adj = stop_adj
                  )
                }
                out
              }

              # Resume only after independently verifying unchanged inputs and streams.
              reused <- FALSE
              checkpoint <- file.path(RES_DIR, 'event_audit', paste0(tag, '.rds'))
              if (identical(Sys.getenv('JBES_RESUME_VERIFIED'), '1') && file.exists(checkpoint)) {
                saved <- readRDS(checkpoint)
                stopifnot(
                  identical(saved$tag, tag), saved$k_star == k_star,
                  saved$horizon == mT, length(saved$replications) == N_MC_ALT
                )
                reps <- saved$replications
                reused <- TRUE
              } else reps <- par_lapply_safe(1:N_MC_ALT, one_rep,
                n_cores = N_CORES, seed = SEED_BASE + 1000L * cnt
              )

              # Existing verified records already contain the event audit.
              # Rebuild missing sidecars, without rewriting complete caches.
              audit_names <- as.vector(outer(c('stop_asym', 'stop_adj'),
                c('events', 'paired'), paste,
                sep = '_'
              ))
              audit_files <- file.path(
                RES_DIR, 'event_audit',
                paste0(tag, '_', audit_names, '.csv')
              )
              if (!reused || !all(file.exists(audit_files)))
                save_mc_event_audit(reps, METHODS, k_star, mT, tag)

              for (method in METHODS) {
                stop_asym <- sapply(reps, function(rr) rr[[method]]["stop_asym"])
                stop_adj <- sapply(reps, function(rr) rr[[method]]["stop_adj"])

                metrics_asym <- detection_metrics(stop_asym, k_star, mT)
                metrics_adj <- detection_metrics(stop_adj, k_star, mT)

                perf_rows[[length(perf_rows) + 1L]] <- data.frame(
                  dgp = dgp, m = m, T = T, gamma = gamma,
                  break_frac = bf, k_star = k_star,
                  change_type = chg, Delta = Delta,
                  method = method,
                  power_asym = metrics_asym["DP_horizon"],
                  power_adj = metrics_adj["DP_horizon"],
                  PFA_pre_asym = metrics_asym["PFA_pre"],
                  PFA_pre_adj = metrics_adj["PFA_pre"],
                  DP_post_asym = metrics_asym["DP_post"],
                  DP_post_adj = metrics_adj["DP_post"],
                  CADD_asym = metrics_asym["CADD"],
                  CADD_adj = metrics_adj["CADD"],
                  ARL_asym = metrics_asym["ARL_excess"],
                  ARL_adj = metrics_adj["ARL_excess"],
                  ADD_asym = metrics_asym["ADD_excess"],
                  ADD_adj = metrics_adj["ADD_excess"],
                  N = N_MC_ALT,
                  stringsAsFactors = FALSE
                )
              }
            }
          }
        }
      }
    }
  }
}

perf_df <- do.call(rbind, perf_rows)
out_path <- file.path(RES_DIR, "mean_alt_performance.csv")
write.csv(perf_df, out_path, row.names = FALSE)

msg("Saved: ", out_path)
msg("DONE.")
