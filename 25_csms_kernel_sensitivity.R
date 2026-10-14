# ============================================================
# 25_csms_kernel_sensitivity.R
# - Kernel-sensitivity analysis for CSMS
# - Focal design: DGP1, m = 250, gamma = 0.15, abrupt changes
# - Holds the benchmark bandwidth fixed and varies the HAC kernel
# ============================================================

rm(list = ls())

AUTO_INSTALL <- TRUE
N_CORES <- 1L
SEED_BASE <- 86420

N_MC_SIZE <- 1000
N_CAL <- 1000
N_MC_ALT <- 1000

DGP <- "DGP1"
M <- 250
GAMMA <- 0.15
T_GRID <- c(1, 5)
KERNEL_GRID <- c("bartlett", "parzen", "qs")
ALPHA_SIZE <- c(0.05, 0.10)

EARLY_BREAK_FRAC <- 0.10
LATE_BREAK_FRAC <- 0.80
EARLY_DELTA <- 0.50
LATE_DELTA <- 1.50

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
  stop("Could not locate project root.")
}

ROOT <- bootstrap_find_root()
source(file.path(ROOT, "00_utils.R"))
jbes_set_root(ROOT)
msg("Project root: ", ROOT)

N_CORES <- env_int("JBES_N_CORES", N_CORES)
N_MC_SIZE <- env_int("JBES_N_MC_SIZE", N_MC_SIZE)
N_CAL <- env_int("JBES_N_CAL", N_CAL)
N_MC_ALT <- env_int("JBES_N_MC_ALT", N_MC_ALT)
SAFE_N_CORES <- 1L

RES_DIR <- file.path(ROOT, "results")
if (!dir.exists(RES_DIR)) dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

ensure_pkgs(c("MASS", "Matrix"), auto_install = AUTO_INSTALL)
source(file.path(ROOT, "01_critical_values.R"))
source(file.path(ROOT, "02_methods.R"))
source(file.path(ROOT, "10_dgps_mean.R"))
source(file.path(ROOT, "30_construct_psi.R"))

mean_size_df <- read.csv(file.path(ROOT, "results", "mean_null_size.csv"), stringsAsFactors = FALSE)
mean_alt_df <- read.csv(Sys.getenv("JBES_MEAN_ALT_PATH", file.path(ROOT, "results", "mean_alt_performance.csv")), stringsAsFactors = FALSE)

base_bandwidth <- ceiling(1.06 * (M^(1 / 5)))

kernel_label <- function(kernel) {
  k <- normalize_hac_kernel(kernel)
  if (k == "bartlett") return("Bartlett")
  if (k == "parzen") return("Parzen")
  if (k == "qs") return("Quadratic spectral")
  k
}

res_rows <- list()
counter <- 0L

for (T in T_GRID) {
  T0 <- T
  mT <- as.integer(round(M * T0))

  for (kernel in KERNEL_GRID) {
    counter <- counter + 1L
    klabel <- kernel_label(kernel)
    tag <- sprintf("csms_kernel_DGP1_m%d_T%d_g%.2f_%s", M, T0, GAMMA, normalize_hac_kernel(kernel))
    msg("Kernel sensitivity design: ", tag)

    null_reps <- par_lapply_safe(
      1:N_MC_SIZE,
      function(i) {
        X <- sim_var1_dgp(dgp = DGP, n = M + mT, m_train = M, q = 2, burnin = 250)
        psi <- make_mean_psi(X, M, T0)$psi
        Mvec <- csms_path(psi, M, T0, GAMMA, L = base_bandwidth, kernel = kernel)
        cvals <- sapply(ALPHA_SIZE, function(a) {
          get_critical_value_csms(d = 2, T = T0, gamma = GAMMA, alpha_level = a)
        })
        sapply(cvals, function(cv) first_exceed(Mvec, cv))
      },
      n_cores = SAFE_N_CORES,
      seed = SEED_BASE + 1000L * counter
    )
    null_mat <- do.call(rbind, null_reps)
    size5 <- mean(null_mat[, 1] <= mT)
    size10 <- mean(null_mat[, 2] <= mT)

    cal_reps <- par_lapply_safe(
      1:N_CAL,
      function(i) {
        X <- sim_var1_dgp(dgp = DGP, n = M + mT, m_train = M, q = 2, burnin = 250)
        psi <- make_mean_psi(X, M, T0)$psi
        max(csms_path(psi, M, T0, GAMMA, L = base_bandwidth, kernel = kernel))
      },
      n_cores = SAFE_N_CORES,
      seed = SEED_BASE + 50000L + 1000L * counter
    )
    c_hat <- as.numeric(quantile(unlist(cal_reps), probs = 0.95, type = 8))

    alt_summary <- function(break_frac, Delta, suffix) {
      k_star <- as.integer(floor(break_frac * mT))
      k_star <- max(1L, min(k_star, mT - 1L))
      stopv <- par_lapply_safe(
        1:N_MC_ALT,
        function(i) {
          X <- sim_var1_dgp(dgp = DGP, n = M + mT, m_train = M, q = 2, burnin = 250)
          Y <- apply_mean_change(X, M, T0,
            k_star = k_star, Delta = Delta,
            shift_vec = rep(1, 2), type = "abrupt"
          )
          psi <- make_mean_psi(Y, M, T0)$psi
          Mvec <- csms_path(psi, M, T0, GAMMA, L = base_bandwidth, kernel = kernel)
          first_exceed(Mvec, c_hat)
        },
        n_cores = SAFE_N_CORES,
        seed = SEED_BASE + 90000L + 1000L * counter + ifelse(suffix == "early", 10L, 20L)
      )
      stopv <- unlist(stopv)
      metrics <- detection_metrics(stopv, k_star, mT)
      c(
        power_adj = unname(metrics["DP_horizon"]),
        PFA_pre_adj = unname(metrics["PFA_pre"]),
        DP_post_adj = unname(metrics["DP_post"]),
        CADD_adj = unname(metrics["CADD"]),
        ADD_adj = unname(metrics["ADD_excess"])
      )
    }

    early <- alt_summary(EARLY_BREAK_FRAC, EARLY_DELTA, "early")
    late <- alt_summary(LATE_BREAK_FRAC, LATE_DELTA, "late")

    res_rows[[length(res_rows) + 1L]] <- data.frame(
      method = "CSMS",
      T = T0,
      kernel = normalize_hac_kernel(kernel),
      kernel_label = klabel,
      L = base_bandwidth,
      size_05 = size5,
      size_10 = size10,
      early_power_adj = early["power_adj"],
      early_PFA_pre_adj = early["PFA_pre_adj"],
      early_DP_post_adj = early["DP_post_adj"],
      early_CADD_adj = early["CADD_adj"],
      early_ADD_adj = early["ADD_adj"],
      late_power_adj = late["power_adj"],
      late_PFA_pre_adj = late["PFA_pre_adj"],
      late_DP_post_adj = late["DP_post_adj"],
      late_CADD_adj = late["CADD_adj"],
      late_ADD_adj = late["ADD_adj"],
      stringsAsFactors = FALSE
    )
  }

  rsms_size_05 <- mean_size_df$size[
    mean_size_df$dgp == DGP &
      mean_size_df$m == M &
      mean_size_df$T == T0 &
      abs(mean_size_df$gamma - GAMMA) < 1e-12 &
      mean_size_df$method == "RSMS" &
      abs(mean_size_df$alpha - 0.05) < 1e-12
  ]
  rsms_size_10 <- mean_size_df$size[
    mean_size_df$dgp == DGP &
      mean_size_df$m == M &
      mean_size_df$T == T0 &
      abs(mean_size_df$gamma - GAMMA) < 1e-12 &
      mean_size_df$method == "RSMS" &
      abs(mean_size_df$alpha - 0.10) < 1e-12
  ]
  rsms_early <- mean_alt_df[
    mean_alt_df$dgp == DGP &
      mean_alt_df$m == M &
      mean_alt_df$T == T0 &
      abs(mean_alt_df$gamma - GAMMA) < 1e-12 &
      mean_alt_df$change_type == "abrupt" &
      abs(mean_alt_df$break_frac - EARLY_BREAK_FRAC) < 1e-12 &
      abs(mean_alt_df$Delta - EARLY_DELTA) < 1e-12 &
      mean_alt_df$method == "RSMS",
  ]
  rsms_late <- mean_alt_df[
    mean_alt_df$dgp == DGP &
      mean_alt_df$m == M &
      mean_alt_df$T == T0 &
      abs(mean_alt_df$gamma - GAMMA) < 1e-12 &
      mean_alt_df$change_type == "abrupt" &
      abs(mean_alt_df$break_frac - LATE_BREAK_FRAC) < 1e-12 &
      abs(mean_alt_df$Delta - LATE_DELTA) < 1e-12 &
      mean_alt_df$method == "RSMS",
  ]

  res_rows[[length(res_rows) + 1L]] <- data.frame(
    method = "RSMS",
    T = T0,
    kernel = NA_character_,
    kernel_label = NA_character_,
    L = NA_integer_,
    size_05 = rsms_size_05[1],
    size_10 = rsms_size_10[1],
    early_power_adj = rsms_early$power_adj[1],
    early_PFA_pre_adj = rsms_early$PFA_pre_adj[1],
    early_DP_post_adj = rsms_early$DP_post_adj[1],
    early_CADD_adj = rsms_early$CADD_adj[1],
    early_ADD_adj = rsms_early$ADD_adj[1],
    late_power_adj = rsms_late$power_adj[1],
    late_PFA_pre_adj = rsms_late$PFA_pre_adj[1],
    late_DP_post_adj = rsms_late$DP_post_adj[1],
    late_CADD_adj = rsms_late$CADD_adj[1],
    late_ADD_adj = rsms_late$ADD_adj[1],
    stringsAsFactors = FALSE
  )
}

out_df <- do.call(rbind, res_rows)
out_path <- file.path(RES_DIR, "csms_kernel_sensitivity.csv")
write.csv(out_df, out_path, row.names = FALSE)
msg("Saved: ", out_path)
msg("DONE.")
