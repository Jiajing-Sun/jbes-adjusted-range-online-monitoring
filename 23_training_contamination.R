# ============================================================
# 23_training_contamination.R
# - Appendix experiment: contaminated training window
#   (a) mild contamination: l_m = o(sqrt(m)) last training obs shifted
#   (b) fixed-fraction contamination: pi*m training obs shifted
# - Reports:
#   * size under contaminated training (using asymptotic CVs under clean null)
#   * size-adjusted post-break detection metrics for a monitoring break when training is contaminated
# ============================================================

rm(list = ls())

AUTO_INSTALL <- TRUE
N_CORES <- 1L
SEED_BASE <- 77889

N_SIZE <- 1000 # size reps
N_CAL <- 1000 # calibration reps for size-adjusted power under contamination
N_MC <- 1000 # power reps under contamination

METHODS <- c("RSMS", "SSMS", "CSMS")
ALPHA <- 0.05

# Base mean monitoring design
dgp <- "DGP2"
q <- 2
m <- 500
T <- 2
mT <- as.integer(round(m * T))
GAMMA_GRID <- c(0, 0.15)

BREAK_FRAC <- c(0.1, 0.8) # break in monitoring window
DELTA <- 0.5 # monitoring break magnitude (target)
CHANGE_TYPE <- "abrupt"
shift_vec <- rep(1, q)

# Contamination regimes
l_mild <- ceiling(m^0.4) # o(sqrt(m))
pi_frac <- 0.2
l_frac <- ceiling(pi_frac * m)

B_GRID <- round(seq(0, 1, by = 0.05), 2) # contamination magnitude


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
N_SIZE <- env_int("JBES_N_SIZE", N_SIZE)
N_CAL <- env_int("JBES_N_CAL", N_CAL)
N_MC <- env_int("JBES_N_MC", N_MC)
GAMMA_GRID <- env_num_vec("JBES_GAMMA_GRID", GAMMA_GRID)
BREAK_FRAC <- env_num_vec("JBES_BREAK_FRAC", BREAK_FRAC)
B_GRID <- env_num_vec("JBES_B_GRID", B_GRID)
DELTA <- env_num("JBES_DELTA", DELTA)

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
source(file.path(ROOT, "40_calibration.R"))

msg("Training contamination experiment starting...")

# Apply contamination only to the TRAINING window (monitoring is clean under H0)
# Apply contamination only to the TRAINING window (monitoring is clean under H0)
apply_training_contamination <- function(X, m, l_contam, B, shift_vec) {
  Y <- X
  if (l_contam <= 0) return(Y)
  if (length(B) == 1L && isTRUE(all.equal(B, 0))) return(Y)

  idx <- (m - l_contam + 1L):m

  # Ensure matrix shape (avoid accidental vector drops when q=1)
  Ysub <- Y[idx, , drop = FALSE]

  shift_vec <- as.numeric(shift_vec)
  q <- ncol(Ysub)
  if (length(shift_vec) != q) {
    stop(sprintf("shift_vec has length %d but data has q=%d columns.", length(shift_vec), q))
  }

  # Allow B to be scalar or a length-|idx| vector
  if (length(B) == 1L) {
    Bvec <- rep(B, length(idx))
  } else if (length(B) == length(idx)) {
    Bvec <- as.numeric(B)
  } else {
    stop(sprintf("B must be scalar or length %d (got length %d).", length(idx), length(B)))
  }

  Y[idx, ] <- Ysub + outer(Bvec, shift_vec)
  Y
}

rows <- list()
cnt <- 0L

for (gamma in GAMMA_GRID) {
  for (regime in c("mild", "fraction")) {
    l_contam <- if (regime == "mild") l_mild else l_frac

    for (B in B_GRID) {
      cnt <- cnt + 1L
      tag <- sprintf("contam_%s_l%d_B%.2f_g%.2f", regime, l_contam, B, gamma)
      msg("Regime: ", tag)

      # --------------------------------------------------------
      # (1) SIZE under contaminated training (H0)
      # Use asymptotic CVs (clean-null) as the benchmark.
      # --------------------------------------------------------
      one_rep_size <- function(i) {
        X <- sim_var1_dgp(dgp = dgp, n = m + mT, m_train = m, q = q, burnin = 300)
        Y <- apply_training_contamination(X, m, l_contam, B, shift_vec)
        psi <- make_mean_psi(Y, m, T)$psi
        sapply(METHODS, function(method) {
          Mvec <- compute_path(method, psi, m, T, gamma)
          c_asym <- switch(method,
            "RSMS" = get_critical_value_rsms(d = q, T = T, gamma = gamma, alpha_level = ALPHA),
            "SSMS" = get_critical_value_ssms(d = q, T = T, gamma = gamma, alpha_level = ALPHA),
            "CSMS" = get_critical_value_csms(d = q, T = T, gamma = gamma, alpha_level = ALPHA)
          )
          as.integer(first_exceed(Mvec, c_asym) <= mT)
        })
      }
      reps_size <- par_lapply_safe(1:N_SIZE, one_rep_size,
        n_cores = N_CORES, seed = SEED_BASE + 1000L * cnt
      )
      mat_size <- do.call(rbind, reps_size)
      colnames(mat_size) <- METHODS

      for (method in METHODS) {
        rows[[length(rows) + 1L]] <- data.frame(
          exp = "size_contaminated_training",
          regime = regime, l_contam = l_contam, B = B,
          gamma = gamma, break_frac = NA, k_star = NA,
          method = method,
          size = mean(mat_size[, method]),
          power_adj = NA, PFA_pre_adj = NA, DP_post_adj = NA, CADD_adj = NA,
          ARL_adj = NA, ADD_adj = NA,
          N = N_SIZE,
          stringsAsFactors = FALSE
        )
      }

      # --------------------------------------------------------
      # (2) Calibrate critical values under this contaminated H0 (for size-adjusted power)
      # --------------------------------------------------------
      gen_psi_null <- function(i) {
        X <- sim_var1_dgp(dgp = dgp, n = m + mT, m_train = m, q = q, burnin = 300)
        Y <- apply_training_contamination(X, m, l_contam, B, shift_vec)
        make_mean_psi(Y, m, T)$psi
      }
      c_hat_df <- calibrate_critical_values(
        METHODS = METHODS, gen_psi = gen_psi_null,
        m = m, T = T, gamma = gamma,
        alpha = ALPHA, N_CAL = N_CAL,
        n_cores = N_CORES,
        seed = SEED_BASE + 50000L + 1000L * cnt,
        tag = tag
      )
      c_hat_map <- setNames(c_hat_df$c_hat, c_hat_df$method)

      # --------------------------------------------------------
      # (3) POWER under monitoring break with the same contaminated training regime
      # --------------------------------------------------------
      for (bf in BREAK_FRAC) {
        k_star <- as.integer(floor(bf * mT))
        k_star <- max(1L, min(k_star, mT - 1L))

        one_rep_pow <- function(i) {
          X <- sim_var1_dgp(dgp = dgp, n = m + mT, m_train = m, q = q, burnin = 300)
          Y0 <- apply_training_contamination(X, m, l_contam, B, shift_vec)
          Y <- apply_mean_change(Y0, m, T,
            k_star = k_star, Delta = DELTA,
            shift_vec = shift_vec, type = CHANGE_TYPE
          )
          psi <- make_mean_psi(Y, m, T)$psi

          sapply(METHODS, function(method) {
            Mvec <- compute_path(method, psi, m, T, gamma)
            first_exceed(Mvec, c_hat_map[[method]])
          })
        }
        reps_pow <- par_lapply_safe(1:N_MC, one_rep_pow,
          n_cores = N_CORES, seed = SEED_BASE + 80000L + 1000L * cnt + as.integer(100 * bf)
        )
        mat_stop <- do.call(rbind, reps_pow)
        colnames(mat_stop) <- METHODS

        for (method in METHODS) {
          stops <- mat_stop[, method]
          metrics <- detection_metrics(stops, k_star, mT)

          rows[[length(rows) + 1L]] <- data.frame(
            exp = "power_with_contaminated_training",
            regime = regime, l_contam = l_contam, B = B,
            gamma = gamma, break_frac = bf, k_star = k_star,
            method = method,
            size = NA,
            power_adj = metrics["DP_horizon"],
            PFA_pre_adj = metrics["PFA_pre"],
            DP_post_adj = metrics["DP_post"],
            CADD_adj = metrics["CADD"],
            ARL_adj = metrics["ARL_excess"],
            ADD_adj = metrics["ADD_excess"],
            N = N_MC,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}

df <- do.call(rbind, rows)
out_path <- file.path(RES_DIR, "training_contamination_results.csv")
write.csv(df, out_path, row.names = FALSE)

msg("Saved: ", out_path)
msg("DONE.")
