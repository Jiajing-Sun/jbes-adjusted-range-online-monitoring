# ============================================================
# 03_mean_summary_tables.R
# - Reads results produced by scripts 01 & 02
# - Produces wide-format CSV summaries and optional LaTeX tables
# ============================================================

rm(list = ls())

AUTO_INSTALL <- TRUE
WRITE_LATEX <- TRUE


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

# Paths
RES_DIR <- file.path(ROOT, "results")
LOG_DIR <- file.path(ROOT, "logs")
if (!dir.exists(RES_DIR)) dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
ensure_pkgs(c("data.table"), auto_install = AUTO_INSTALL)

library(data.table)

size_path <- file.path(RES_DIR, "mean_null_size.csv")
alt_path <- file.path(RES_DIR, "mean_alt_performance.csv")

if (!file.exists(size_path)) stop("Missing: ", size_path)
if (!file.exists(alt_path)) stop("Missing: ", alt_path)

size_dt <- fread(size_path)
alt_dt <- fread(alt_path)

# --------
# Size summary (wide): one row per (dgp,m,gamma,alpha,method), columns T
# --------
size_wide <- dcast(size_dt, dgp + m + gamma + alpha + method ~ T, value.var = "size")
fwrite(size_wide, file.path(RES_DIR, "mean_null_size_wide.csv"))

# --------
# Alt summary: keep both power_adj and ARL_adj (size-adjusted)
# One row per (dgp,m,T,gamma,break_frac,change_type,Delta,method)
# --------
alt_keep <- alt_dt[, .(
  dgp, m, T, gamma, break_frac, change_type, Delta, method,
  PFA_pre_asym, PFA_pre_adj,
  DP_post_asym, DP_post_adj, CADD_asym, CADD_adj,
  power_asym, power_adj, ARL_asym, ARL_adj, ADD_asym, ADD_adj
)]

fwrite(alt_keep, file.path(RES_DIR, "mean_alt_performance_clean.csv"))

# Optional: LaTeX tables (requires xtable)
if (WRITE_LATEX) {
  ensure_pkgs(c("xtable"), auto_install = AUTO_INSTALL)
  library(xtable)

  # Example: size table for DGP1, gamma=0, alpha=0.05
  ex <- size_dt[dgp == "DGP1" & abs(gamma - 0) < 1e-12 & abs(alpha - 0.05) < 1e-12]
  ex_w <- dcast(ex, m + method ~ T, value.var = "size")
  tex <- print(xtable(ex_w, digits = 3),
    include.rownames = FALSE, sanitize.text.function = identity
  )
  writeLines(tex, con = file.path(RES_DIR, "latex_example_mean_size_DGP1_g0_a05.tex"))
}

msg("Wrote summaries to results/.")
msg("DONE.")
