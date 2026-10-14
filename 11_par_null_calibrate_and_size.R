# ============================================================
# 11_par_null_calibrate_and_size.R
# - Null simulations for PAR(1) monitoring (Appendix)
# - Reports size at 5% and 10% using asymptotic critical values
# - Computes empirical critical values (size-adjusted) at alpha=0.05
# ============================================================

rm(list = ls())

AUTO_INSTALL <- TRUE
N_CORES <- 1L
SEED_BASE <- 54321

N_MC  <- 1000
N_CAL <- 1000

M_GRID     <- c(100, 250, 500)
T_GRID     <- c(1, 2, 5, 10)
GAMMA_GRID <- c(0, 0.15)
ALPHA_SIZE <- c(0.05, 0.10)

METHODS <- c("RSMS","SSMS","CSMS")

# Null parameters
ALPHA1 <- 0.2
LAMBDA1 <- 2.0


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
N_MC <- env_int("JBES_N_MC", N_MC)
N_CAL <- env_int("JBES_N_CAL", N_CAL)
M_GRID <- env_num_vec("JBES_M_GRID", M_GRID)
T_GRID <- env_num_vec("JBES_T_GRID", T_GRID)
GAMMA_GRID <- env_num_vec("JBES_GAMMA_GRID", GAMMA_GRID)

# Paths
RES_DIR <- file.path(ROOT, "results")
LOG_DIR <- file.path(ROOT, "logs")
if (!dir.exists(RES_DIR)) dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
ensure_pkgs(c("MASS", "Matrix"), auto_install = AUTO_INSTALL)

source(file.path(ROOT, "01_critical_values.R"))
source(file.path(ROOT, "02_methods.R"))
source(file.path(ROOT, "20_dgps_par1.R"))
source(file.path(ROOT, "40_calibration.R"))

msg("Starting NULL simulations (PAR(1))")

size_rows <- list()
calib_rows <- list()
cnt <- 0L

for (m in M_GRID) {
  for (T in T_GRID) {
    for (gamma in GAMMA_GRID) {
      cnt <- cnt + 1L
      tag <- sprintf("par_null_m%d_T%d_g%.2f", m, T, gamma)
      msg("Design point: ", tag)

      mT <- as.integer(round(m * T))

      # Precompute asymptotic CVs
      c_asym <- list()
      for (method in METHODS) {
        c_asym[[method]] <- sapply(ALPHA_SIZE, function(a) {
          switch(method,
                 "RSMS" = get_critical_value_rsms(d = 2, T = T, gamma = gamma, alpha_level = a),
                 "SSMS" = get_critical_value_ssms(d = 2, T = T, gamma = gamma, alpha_level = a),
                 "CSMS" = get_critical_value_csms(d = 2, T = T, gamma = gamma, alpha_level = a))
        })
        names(c_asym[[method]]) <- as.character(ALPHA_SIZE)
      }

      # NULL size MC (store only stop times)
      one_rep_size <- function(i) {
        X <- sim_par1(n = m + mT + 1L, alpha = ALPHA1, lambda = LAMBDA1, burnin = 400)

        fit <- par1_fit_mle(X_train = X[1:(m+1L)])
        psi <- par1_score_matrix(X, alpha = fit$alpha, lambda = fit$lambda)  # (m+mT) x 2

        out <- numeric(length(METHODS) * length(ALPHA_SIZE))
        nm <- character(length(out))
        idx <- 0L

        for (method in METHODS) {
          Mvec <- compute_path(method, psi, m, T, gamma)
          for (a in ALPHA_SIZE) {
            idx <- idx + 1L
            nm[idx] <- paste0(method, "_a", a)
            out[idx] <- first_exceed(Mvec, c_asym[[method]][as.character(a)])
          }
        }
        names(out) <- nm
        out
      }

      reps <- par_lapply_safe(1:N_MC, one_rep_size,
                              n_cores = N_CORES, seed = SEED_BASE + 1000L * cnt)
      stop_mat <- do.call(rbind, reps)

      for (a in ALPHA_SIZE) {
        for (method in METHODS) {
          stops <- stop_mat[, paste0(method, "_a", a)]
          size_hat <- mean(stops <= mT)

          size_rows[[length(size_rows) + 1L]] <- data.frame(
            m = m, T = T, gamma = gamma, method = method,
            alpha = a, size = size_hat, N = N_MC,
            stringsAsFactors = FALSE
          )
        }
      }

      # Calibration for size-adjusted power (alpha=0.05)
      gen_psi <- function(i) {
        X <- sim_par1(n = m + mT + 1L, alpha = ALPHA1, lambda = LAMBDA1, burnin = 400)
        fit <- par1_fit_mle(X_train = X[1:(m+1L)])
        par1_score_matrix(X, alpha = fit$alpha, lambda = fit$lambda)
      }

      calib <- calibrate_critical_values(METHODS = METHODS, gen_psi = gen_psi,
                                         m = m, T = T, gamma = gamma,
                                         alpha = 0.05, N_CAL = N_CAL,
                                         n_cores = N_CORES, seed = SEED_BASE + 9999L + 1000L * cnt,
                                         tag = tag)
      calib_rows[[length(calib_rows) + 1L]] <- calib
    }
  }
}

size_df <- do.call(rbind, size_rows)
calib_df <- do.call(rbind, calib_rows)

write.csv(size_df,  file.path(RES_DIR, "par_null_size.csv"), row.names = FALSE)
write.csv(calib_df, file.path(RES_DIR, "par_null_calibration_cHat.csv"), row.names = FALSE)

msg("Saved PAR null results to results/.")
msg("DONE.")
