# ============================================================
# 12_par_alternatives_power_delay.R
# - Alternative simulations for PAR(1) monitoring (Appendix)
# - Produces grid-based horizon rejection and strict post-break detection metrics
# ============================================================

rm(list = ls())

AUTO_INSTALL <- TRUE
N_CORES <- 1L
SEED_BASE <- 13579

N_MC_ALT <- 1000

METHODS <- c("RSMS", "SSMS", "CSMS")
ALPHA <- 0.05

# Design
M_GRID <- c(100, 250, 500)
T_GRID <- c(1, 2, 5)
GAMMA_GRID <- c(0, 0.15)
BREAK_FRAC <- c(0.1, 0.8)
CHANGE_TYPES <- c("abrupt", "smooth")

# Null / base parameters
ALPHA1 <- 0.2
LAMBDA1 <- 2.0

# Post-change grid
ALPHA2_GRID <- seq(0.2, 0.35, by = 0.05)
LAMBDA2_GRID <- seq(2.0, 2.5, by = 0.1)


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
M_GRID <- env_num_vec("JBES_M_GRID", M_GRID)
T_GRID <- env_num_vec("JBES_T_GRID", T_GRID)
GAMMA_GRID <- env_num_vec("JBES_GAMMA_GRID", GAMMA_GRID)
BREAK_FRAC <- env_num_vec("JBES_BREAK_FRAC", BREAK_FRAC)
CHANGE_TYPES <- env_chr_vec("JBES_CHANGE_TYPES", CHANGE_TYPES)
ALPHA2_GRID <- env_num_vec("JBES_ALPHA2_GRID", ALPHA2_GRID)
LAMBDA2_GRID <- env_num_vec("JBES_LAMBDA2_GRID", LAMBDA2_GRID)

# Paths
RES_DIR <- file.path(ROOT, "results")
LOG_DIR <- file.path(ROOT, "logs")
if (!dir.exists(RES_DIR)) dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
ensure_pkgs(c("MASS", "Matrix"), auto_install = AUTO_INSTALL)

source(file.path(ROOT, "01_critical_values.R"))
source(file.path(ROOT, "02_methods.R"))
source(file.path(ROOT, "20_dgps_par1.R"))

# Load calibration c_hat under null
calib_path <- file.path(RES_DIR, "par_null_calibration_cHat.csv")
if (!file.exists(calib_path)) stop("Missing: ", calib_path)
calib_df <- read.csv(calib_path, stringsAsFactors = FALSE)
calib_df <- subset(calib_df, abs(alpha - ALPHA) < 1e-12)

get_c_hat <- function(m0, T0, gamma0, method0) {
  sub <- calib_df[calib_df$m == m0 &
    calib_df$T == T0 &
    abs(calib_df$gamma - gamma0) < 1e-12 &
    calib_df$method == method0, ]
  if (nrow(sub) != 1) {
    stop(
      "Missing/duplicate calib for: ",
      paste(m0, T0, gamma0, method0, sep = " / ")
    )
  }
  sub$c_hat[1]
}


msg("Starting ALT simulations (PAR(1))")

# A versioned cache keeps earlier unverified files intact. Every design,
# replication count and random-stream choice must match before reuse.
CHECK_DIR <- file.path(RES_DIR, "par_alternative_checkpoints_v2")
dir.create(CHECK_DIR, recursive = TRUE, showWarnings = FALSE)
manifest <- list(
  m = M_GRID, T = T_GRID, gamma = GAMMA_GRID, break_fraction = BREAK_FRAC,
  change_type = CHANGE_TYPES, alpha2 = ALPHA2_GRID, lambda2 = LAMBDA2_GRID,
  base_parameter = c(ALPHA1, LAMBDA1), alpha = ALPHA, methods = METHODS,
  N = N_MC_ALT, workers = N_CORES, seed = SEED_BASE, R = R.version.string,
  aersn = as.character(packageVersion("aersn")),
  highs = as.character(packageVersion("highs")),
  source = unname(tools::md5sum(c(
    file.path(ROOT, c(
      "00_utils.R",
      "01_critical_values.R", "02_methods.R", "07_aersn_backend.R",
      "08_mc_audit.R", "20_dgps_par1.R"
    )),
    file.path(ROOT, "12_par_alternatives_power_delay.R")
  ))),
  reference = unname(tools::md5sum(file.path(
    RES_DIR,
    c("affine_finite_critical_values.csv", "par_null_calibration_cHat.csv")
  )))
)
stopifnot(!anyNA(manifest$source), !anyNA(manifest$reference))
manifest_path <- file.path(CHECK_DIR, "manifest.rds")
if (file.exists(manifest_path)) {
  if (!identical(readRDS(manifest_path), manifest))
    stop("PAR alternative checkpoint manifest differs; use a fresh results directory.")
} else {
  if (length(list.files(CHECK_DIR, pattern = "\\.rds$")))
    stop("Unverified PAR alternative checkpoints found; preserve them before restarting.")
  saveRDS(manifest, manifest_path)
}
checkpoint_par <- function(tag, expr) {
  path <- file.path(CHECK_DIR, paste0(tag, ".rds"))
  if (file.exists(path)) {
    value <- readRDS(path)
    stopifnot(is.list(value), length(value) == N_MC_ALT)
    return(value)
  }
  value <- force(expr)
  stopifnot(is.list(value), length(value) == N_MC_ALT)
  temporary <- paste0(path, ".tmp");saveRDS(value, temporary)
  if (!file.rename(temporary, path)) stop("Cannot save PAR alternative checkpoint: ", path)
  value
}
perf_rows <- list()
cnt <- 0L

for (m in M_GRID) {
  for (T in T_GRID) {
    mT <- as.integer(round(m * T))
    for (gamma in GAMMA_GRID) {
      for (bf in BREAK_FRAC) {
        k_star <- as.integer(floor(bf * mT))
        k_star <- max(1L, min(k_star, mT - 1L))

        for (chg in CHANGE_TYPES) {
          for (a2 in ALPHA2_GRID) {
            for (l2 in LAMBDA2_GRID) {
              cnt <- cnt + 1L
              tag <- sprintf(
                "par_alt_m%d_T%d_g%.2f_bf%.1f_%s_a2%.2f_l2%.1f",
                m, T, gamma, bf, chg, a2, l2
              )
              msg("Design: ", tag)

              c_asym_map <- c(
                RSMS = get_critical_value_rsms(2, T, gamma, ALPHA),
                SSMS = get_critical_value_ssms(2, T, gamma, ALPHA), CSMS = get_critical_value_csms(2, T, gamma, ALPHA)
              )
              c_hat_map <- setNames(vapply(METHODS, function(method) get_c_hat(m, T, gamma, method), numeric(1)), METHODS)

              one_rep <- function(i) {
                X <- sim_par1_break(
                  m = m, T = T, k_star = k_star,
                  alpha1 = ALPHA1, lambda1 = LAMBDA1,
                  alpha2 = a2, lambda2 = l2,
                  type = chg, delta_smooth = 0.05,
                  burnin = 500
                )

                # Fit on training obs (length m+1)
                fit <- par1_fit_mle(X_train = X[1:(m + 1L)])
                psi <- par1_score_matrix(X, alpha = fit$alpha, lambda = fit$lambda) # (m+mT) x 2

                out <- list()
                for (method in METHODS) {
                  Mvec <- compute_path(method, psi, m, T, gamma)
                  c_asym <- c_asym_map[[method]]
                  stop_asym <- first_exceed(Mvec, c_asym)

                  c_hat <- c_hat_map[[method]]
                  stop_adj <- first_exceed(Mvec, c_hat)

                  out[[method]] <- c(stop_asym = stop_asym, stop_adj = stop_adj)
                }
                out
              }

              reps <- checkpoint_par(tag, par_lapply_safe(1:N_MC_ALT, one_rep,
                n_cores = N_CORES, seed = SEED_BASE + 1000L * cnt
              ))

              save_mc_event_audit(reps, METHODS, k_star, mT, tag)

              for (method in METHODS) {
                stop_asym <- sapply(reps, function(rr) rr[[method]]["stop_asym"])
                stop_adj <- sapply(reps, function(rr) rr[[method]]["stop_adj"])

                metrics_asym <- detection_metrics(stop_asym, k_star, mT)
                metrics_adj <- detection_metrics(stop_adj, k_star, mT)

                perf_rows[[length(perf_rows) + 1L]] <- data.frame(
                  m = m, T = T, gamma = gamma, break_frac = bf, k_star = k_star,
                  change_type = chg,
                  alpha2 = a2, lambda2 = l2,
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
out_path <- file.path(RES_DIR, "par_alt_performance_grid.csv")
write.csv(perf_df, out_path, row.names = FALSE)

msg("Saved: ", out_path)
msg("DONE.")
