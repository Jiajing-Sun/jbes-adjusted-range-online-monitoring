# ============================================================
# 37_dependence_assets.R
# - Paper-ready figures for dependence-sensitivity appendix results
# ============================================================

rm(list = ls())

AUTO_INSTALL <- TRUE

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
      if (file.exists(file.path(cur, "00_utils.R"))) return(cur)
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

RES_DIR <- file.path(ROOT, "results")
FIG_DIR <- file.path(ROOT, "figures")
if (!dir.exists(RES_DIR)) dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

ensure_pkgs(c("data.table", "ggplot2", "patchwork"), auto_install = AUTO_INSTALL)
library(data.table)
library(ggplot2)
library(patchwork)

path <- file.path(RES_DIR, "dependence_sensitivity_results.csv")
if (!file.exists(path)) stop("Missing: ", path)
dt <- fread(path)

method_linetypes <- c(RSMS = "solid", SSMS = "longdash", CSMS = "dotdash")
method_shapes <- c(RSMS = 16, SSMS = 2, CSMS = 0)

plot_dep <- function(sub_dt, x_var, x_lab, metric, ylab, file_stub, zoom_late = FALSE) {
  gamma_vals <- sort(unique(sub_dt$gamma))
  break_vals <- sort(unique(sub_dt$break_frac))

  make_panel <- function(gamma0, break0, zoom_y = FALSE) {
    panel_dt <- sub_dt[abs(gamma - gamma0) < 1e-12 & abs(break_frac - break0) < 1e-12]
    p <- ggplot(
      panel_dt,
      aes(x = .data[[x_var]], y = .data[[metric]], linetype = method, shape = method, group = method)
    ) +
      geom_line(linewidth = 0.75, colour = "black") +
      geom_point(size = 2, colour = "black") +
      scale_x_continuous(breaks = scales::breaks_pretty(n = 4)) +
      scale_y_continuous(breaks = scales::breaks_pretty(n = 3), expand = expansion(mult = c(0.08, 0.10))) +
      scale_linetype_manual(values = method_linetypes) +
      scale_shape_manual(values = method_shapes) +
      labs(
        title = bquote(gamma == .(gamma0) * "," ~ "break" == .(break0)),
        x = x_lab,
        y = ylab,
        linetype = "Method",
        shape = "Method"
      ) +
      theme_bw(base_size = 10) +
      theme(
        legend.position = "bottom",
        legend.key.width = grid::unit(1.4, "cm"),
        plot.title = element_text(size = 10, face = "bold"),
        panel.grid.minor = element_blank()
      )
    p
  }

  plot_list <- list()
  for (gamma0 in gamma_vals) {
    for (break0 in break_vals) {
      zoom_y <- zoom_late && break0 >= 0.8
      plot_list[[length(plot_list) + 1L]] <- make_panel(gamma0, break0, zoom_y)
    }
  }

  p <- wrap_plots(plotlist = plot_list, ncol = 2, byrow = TRUE, guides = "collect") &
    theme(legend.position = "bottom")

  out <- file.path(FIG_DIR, paste0(file_stub, ".pdf"))
  ggsave(out, p, width = 8, height = 3.4)
  msg("Wrote: ", out)
}

plot_dep(dt[exp == "VAR_phi"], "phi", expression(phi), "DP_post_adj",
         expression(DP[post]), "dependence_VAR_power_adj", zoom_late = TRUE)
plot_dep(dt[exp == "VAR_phi"], "phi", expression(phi), "CADD_adj",
         "CADD", "dependence_VAR_ADD_adj", zoom_late = TRUE)
plot_dep(dt[exp == "GARCH_persist"], "garch_persist", expression(alpha[1] + beta[1]), "DP_post_adj",
         expression(DP[post]), "dependence_GARCH_power_adj", zoom_late = TRUE)
plot_dep(dt[exp == "GARCH_persist"], "garch_persist", expression(alpha[1] + beta[1]), "CADD_adj",
         "CADD", "dependence_GARCH_ADD_adj", zoom_late = TRUE)

msg("DONE.")
