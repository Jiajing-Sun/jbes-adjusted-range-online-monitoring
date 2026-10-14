# ============================================================
# 22_dependence_sensitivity.R
# - Appendix experiment: sensitivity to dependence strength (VAR persistence / GARCH persistence)
# - Uses size-adjusted post-break detection metrics (MC calibrated CVs)
# ============================================================

rm(list = ls())

AUTO_INSTALL <- TRUE
N_CORES <- 1L
SEED_BASE <- 33445

N_CAL <- 1000
N_MC <- 1000
ALPHA <- 0.05

METHODS <- c("RSMS", "SSMS", "CSMS")

# Base monitoring design
m <- 500
T <- 2
mT <- as.integer(round(m * T))
GAMMA_GRID <- c(0, 0.15)
BREAK_FRAC <- c(0.1, 0.8)
CHANGE_TYPE <- "abrupt"
DELTA <- 0.5
q <- 2
shift_vec <- rep(1, q)

# Dependence configurations
PHI_GRID <- round(seq(0.1, 0.9, by = 0.1), 2) # VAR persistence
GARCH_PERSIST <- round(seq(0.1, 0.9, by = 0.1), 2) # alpha+beta (keep alpha fixed)


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
N_MC <- env_int("JBES_N_MC", N_MC)
GAMMA_GRID <- env_num_vec("JBES_GAMMA_GRID", GAMMA_GRID)
BREAK_FRAC <- env_num_vec("JBES_BREAK_FRAC", BREAK_FRAC)
PHI_GRID <- env_num_vec("JBES_PHI_GRID", PHI_GRID)
GARCH_PERSIST <- env_num_vec("JBES_GARCH_PERSIST", GARCH_PERSIST)
DELTA <- env_num("JBES_DELTA", DELTA)

# Paths
RES_DIR <- file.path(ROOT, "results")
LOG_DIR <- file.path(ROOT, "logs")
if (!dir.exists(RES_DIR)) dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
ensure_pkgs(c("MASS", "Matrix"), auto_install = AUTO_INSTALL)

source(file.path(ROOT, "02_methods.R"))
source(file.path(ROOT, "10_dgps_mean.R"))
source(file.path(ROOT, "30_construct_psi.R"))
source(file.path(ROOT, "40_calibration.R"))

msg("Dependence sensitivity experiment starting...")

rows <- list()
cnt <- 0L

# --- (A) VAR persistence sensitivity (use DGP2 structure)
for (phi in PHI_GRID) {
  for (gamma in GAMMA_GRID) {
    for (bf in BREAK_FRAC) {
      k_star <- as.integer(floor(bf * mT))
      k_star <- max(1L, min(k_star, mT - 1L))

      tag <- sprintf("dep_VAR_phi%.2f_g%.2f_bf%.1f", phi, gamma, bf)
      msg("Calibrating (VAR): ", tag)

      gen_psi_null <- function(i) {
        X <- sim_var1_dgp(
          dgp = "DGP2", n = m + mT, m_train = m, q = q,
          burnin = 300, phi = phi, cross = 0.1, rho = 0.1
        )
        make_mean_psi(X, m, T)$psi
      }
      c_hat_df <- calibrate_critical_values(
        METHODS = METHODS, gen_psi = gen_psi_null,
        m = m, T = T, gamma = gamma,
        alpha = ALPHA, N_CAL = N_CAL,
        n_cores = N_CORES,
        seed = SEED_BASE + 1000L + as.integer(phi * 100) + as.integer(100 * gamma),
        tag = tag
      )
      c_hat_map <- setNames(c_hat_df$c_hat, c_hat_df$method)

      msg("Power (VAR): ", tag)
      one_rep <- function(i) {
        X <- sim_var1_dgp(
          dgp = "DGP2", n = m + mT, m_train = m, q = q,
          burnin = 300, phi = phi, cross = 0.1, rho = 0.1
        )
        Y <- apply_mean_change(X, m, T,
          k_star = k_star, Delta = DELTA,
          shift_vec = shift_vec, type = CHANGE_TYPE
        )
        psi <- make_mean_psi(Y, m, T)$psi
        sapply(METHODS, function(method) {
          Mvec <- compute_path(method, psi, m, T, gamma)
          first_exceed(Mvec, c_hat_map[[method]])
        })
      }

      reps <- par_lapply_safe(1:N_MC, one_rep,
        n_cores = N_CORES, seed = SEED_BASE + 5000L + as.integer(phi * 100) + cnt
      )
      stop_mat <- do.call(rbind, reps)
      colnames(stop_mat) <- METHODS

      for (method in METHODS) {
        stops <- stop_mat[, method]
        metrics <- detection_metrics(stops, k_star, mT)
        rows[[length(rows) + 1L]] <- data.frame(
          exp = "VAR_phi", phi = phi, garch_persist = NA,
          m = m, T = T, gamma = gamma, break_frac = bf, k_star = k_star,
          method = method,
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
      cnt <- cnt + 1L
    }
  }
}

# --- (B) GARCH persistence sensitivity (use DGP3 structure)
garch_alpha <- 0.1
for (persist in GARCH_PERSIST) {
  garch_beta <- max(0, persist - garch_alpha)
  for (gamma in GAMMA_GRID) {
    for (bf in BREAK_FRAC) {
      k_star <- as.integer(floor(bf * mT))
      k_star <- max(1L, min(k_star, mT - 1L))

      tag <- sprintf("dep_GARCH_p%.2f_g%.2f_bf%.1f", persist, gamma, bf)
      msg("Calibrating (GARCH): ", tag)

      gen_psi_null <- function(i) {
        X <- sim_var1_dgp(
          dgp = "DGP3", n = m + mT, m_train = m, q = q,
          burnin = 400, phi = 0.5, cross = 0.1, rho = 0.1,
          garch_alpha = garch_alpha, garch_beta = garch_beta
        )
        make_mean_psi(X, m, T)$psi
      }
      c_hat_df <- calibrate_critical_values(
        METHODS = METHODS, gen_psi = gen_psi_null,
        m = m, T = T, gamma = gamma,
        alpha = ALPHA, N_CAL = N_CAL,
        n_cores = N_CORES,
        seed = SEED_BASE + 20000L + as.integer(persist * 100) + as.integer(100 * gamma),
        tag = tag
      )
      c_hat_map <- setNames(c_hat_df$c_hat, c_hat_df$method)

      msg("Power (GARCH): ", tag)
      one_rep <- function(i) {
        X <- sim_var1_dgp(
          dgp = "DGP3", n = m + mT, m_train = m, q = q,
          burnin = 400, phi = 0.5, cross = 0.1, rho = 0.1,
          garch_alpha = garch_alpha, garch_beta = garch_beta
        )
        Y <- apply_mean_change(X, m, T,
          k_star = k_star, Delta = DELTA,
          shift_vec = shift_vec, type = CHANGE_TYPE
        )
        psi <- make_mean_psi(Y, m, T)$psi
        sapply(METHODS, function(method) {
          Mvec <- compute_path(method, psi, m, T, gamma)
          first_exceed(Mvec, c_hat_map[[method]])
        })
      }

      reps <- par_lapply_safe(1:N_MC, one_rep,
        n_cores = N_CORES, seed = SEED_BASE + 30000L + as.integer(persist * 100) + cnt
      )
      stop_mat <- do.call(rbind, reps)
      colnames(stop_mat) <- METHODS

      for (method in METHODS) {
        stops <- stop_mat[, method]
        metrics <- detection_metrics(stops, k_star, mT)
        rows[[length(rows) + 1L]] <- data.frame(
          exp = "GARCH_persist", phi = NA, garch_persist = persist,
          m = m, T = T, gamma = gamma, break_frac = bf, k_star = k_star,
          method = method,
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
      cnt <- cnt + 1L
    }
  }
}

df <- do.call(rbind, rows)
out_path <- file.path(RES_DIR, "dependence_sensitivity_results.csv")
write.csv(df, out_path, row.names = FALSE)

msg("Saved: ", out_path)
msg("DONE.")
