# ============================================================
# 13_par_plots.R
# - Heatmaps/contours for PAR(1) appendix simulations
# ============================================================

rm(list = ls())

AUTO_INSTALL <- TRUE


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
FIG_DIR <- file.path(ROOT, "figures")
LOG_DIR <- file.path(ROOT, "logs")
if (!dir.exists(RES_DIR)) dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
ensure_pkgs(c("data.table", "ggplot2"), auto_install = AUTO_INSTALL)

library(data.table)
library(ggplot2)

path <- file.path(RES_DIR, "par_alt_performance_grid.csv")
if (!file.exists(path)) stop("Missing: ", path)
dt <- fread(path)

M_FOCUS <- env_int("JBES_M_FOCUS", 250L)
T_FOCUS <- env_num_vec("JBES_T_FOCUS", c(1, 2, 5))
GAMMA_FOCUS <- env_num_vec("JBES_GAMMA_FOCUS", c(0, 0.15))
BREAK_FOCUS <- env_num_vec("JBES_BREAK_FOCUS", c(0.1, 0.8))
CHANGE_FOCUS <- env_chr_vec("JBES_CHANGE_FOCUS", c("abrupt", "smooth"))

slice_tag <- function(T0, gamma0, bf0, chg0) {
  bf_lab <- if (abs(bf0 - 0.1) < 1e-12) "early" else if (abs(bf0 - 0.8) < 1e-12) "late" else paste0("bf", gsub("\\.", "p", sprintf("%.2f", bf0)))
  g_lab <- gsub("\\.", "p", sprintf("%.2f", gamma0))
  paste0("m", M_FOCUS, "_T", T0, "_g", g_lab, "_", bf_lab, "_", chg0)
}

par_theme <- theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold"),
    panel.grid = element_blank()
  )

for (T0 in T_FOCUS) {
  for (gamma0 in GAMMA_FOCUS) {
    for (bf0 in BREAK_FOCUS) {
      for (chg0 in CHANGE_FOCUS) {
        dt1 <- dt[
          m == M_FOCUS &
          T == T0 &
          abs(gamma - gamma0) < 1e-12 &
          abs(break_frac - bf0) < 1e-12 &
          change_type == chg0
        ]

        if (nrow(dt1) == 0L) next

        tag <- slice_tag(T0, gamma0, bf0, chg0)
        title_stub <- sprintf("PAR(1), m=%d, T=%s, gamma=%.2f, break=%.1f, %s",
                              M_FOCUS, T0, gamma0, bf0, chg0)

        p <- ggplot(dt1, aes(x = lambda2, y = alpha2, fill = DP_post_adj)) +
          geom_tile() +
          facet_wrap(~ method, nrow = 1) +
          labs(title = paste0("Size-adjusted post-break detection probability (", title_stub, ")"),
               x = expression(lambda[2]), y = expression(alpha[2]), fill = expression(DP[post])) +
          scale_fill_gradient(low = "white", high = "black") +
          par_theme

        ggsave(file.path(FIG_DIR, paste0("par_power_adj_heatmap_", tag, ".pdf")),
               p, width = 12.5, height = 4.2)

        p2 <- ggplot(dt1, aes(x = lambda2, y = alpha2, z = CADD_adj)) +
          geom_contour(aes(color = after_stat(level)), linewidth = 0.65, bins = 8, show.legend = TRUE) +
          facet_wrap(~ method, nrow = 1) +
          labs(title = paste0("CADD contours (", title_stub, ")"),
               x = expression(lambda[2]), y = expression(alpha[2]), color = "CADD") +
          scale_color_gradient(low = "grey80", high = "black") +
          par_theme

        ggsave(file.path(FIG_DIR, paste0("par_ADD_adj_contour_", tag, ".pdf")),
               p2, width = 12.5, height = 4.2)

        p3 <- ggplot(dt1, aes(x = lambda2, y = alpha2, z = ARL_adj)) +
          geom_contour(aes(color = after_stat(level)), linewidth = 0.65, bins = 8, show.legend = TRUE) +
          facet_wrap(~ method, nrow = 1) +
          labs(title = paste0("ARL contours (", title_stub, ")"),
               x = expression(lambda[2]), y = expression(alpha[2]), color = "ARL") +
          scale_color_gradient(low = "grey80", high = "black") +
          par_theme

        ggsave(file.path(FIG_DIR, paste0("par_ARL_adj_contour_", tag, ".pdf")),
               p3, width = 12.5, height = 4.2)
      }
    }
  }
}

msg("Saved PAR plots to: ", FIG_DIR)
msg("DONE.")
