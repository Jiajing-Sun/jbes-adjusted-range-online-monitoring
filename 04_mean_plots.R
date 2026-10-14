# ============================================================
# 04_mean_plots.R
# - Example plots for the main mean-monitoring simulations
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
ensure_pkgs(c("data.table", "ggplot2", "patchwork", "showtext"), auto_install = AUTO_INSTALL)

library(data.table)
library(ggplot2)
library(patchwork)
library(showtext)
showtext_auto()

alt_path <- file.path(RES_DIR, "mean_alt_performance.csv")
if (!file.exists(alt_path)) stop("Missing: ", alt_path)
dt <- fread(alt_path)

DGP_FOCUS <- env_chr_vec("JBES_DGP_FOCUS", c("DGP1", "DGP2", "DGP3", "DGP4"))
M_FOCUS <- env_num_vec("JBES_M_FOCUS", c(250))
GAMMA_FOCUS <- env_num_vec("JBES_GAMMA_FOCUS", c(0, 0.15))

save_tag <- function(dgp0, m0, gamma0) {
  paste0(dgp0, "_m", m0, "_g", gsub("\\.", "p", sprintf("%.2f", gamma0)))
}

method_linetypes <- c(RSMS = "solid", SSMS = "longdash", CSMS = "dotdash")
method_shapes <- c(RSMS = 16, SSMS = 2, CSMS = 0)
method_colours <- c(RSMS = "black", SSMS = "black", CSMS = "black")

mean_theme <- theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom",
        legend.key.width = grid::unit(1.4, "cm"),
    legend.title = element_blank(),
    strip.text = element_text(size = 10),
    plot.title = element_text(size = 13, face = "bold"),
    panel.grid.minor = element_blank()
  )

for (dgp0 in DGP_FOCUS) {
  for (m0 in M_FOCUS) {
    for (gamma0 in GAMMA_FOCUS) {
      dt1 <- dt[
        dgp == dgp0 &
        m == m0 &
        abs(gamma - gamma0) < 1e-12
      ]

      if (nrow(dt1) == 0L) next

      dt1[, panel := ifelse(
        break_frac < 0.5,
        paste0("Early, T=", T, ", ", change_type),
        paste0("Late, T=", T, ", ", change_type)
      )]
      dt1[, timing := ifelse(break_frac < 0.5, "Early break", "Late break")]
      dt1[, spec := paste0("T=", T, ", ", change_type)]
      panel_levels <- c(
        as.vector(outer(c(1, 2, 5), c("abrupt", "smooth"), function(tt, cc) paste0("Early, T=", tt, ", ", cc))),
        as.vector(outer(c(1, 2, 5), c("abrupt", "smooth"), function(tt, cc) paste0("Late, T=", tt, ", ", cc)))
      )
      dt1[, panel := factor(panel, levels = panel_levels)]
      spec_levels <- as.vector(outer(c(1, 2, 5), c("abrupt", "smooth"), function(tt, cc) paste0("T=", tt, ", ", cc)))
      dt1[, spec := factor(spec, levels = spec_levels)]

      tag <- save_tag(dgp0, m0, gamma0)

      title_suffix <- sprintf("(%s, m = %d, \u03b3 = %.2f)", dgp0, m0, gamma0)

      p1 <- ggplot(dt1, aes(x = Delta, y = DP_post_adj, group = method, linetype = method, shape = method)) +
        geom_line(linewidth = 0.7, colour = "black") +
        geom_point(size = 1.9, colour = "black") +
        facet_wrap(~ panel, ncol = 3, scales = "free") +
        labs(title = paste("Size-adjusted post-break detection probability vs \u0394", title_suffix),
             x = "\u0394", y = "Post-break detection probability (size-adjusted)") +
        scale_linetype_manual(values = method_linetypes) +
        scale_shape_manual(values = method_shapes) +
        mean_theme

      ggsave(filename = file.path(FIG_DIR, paste0("mean_power_adj_", tag, ".pdf")),
             plot = p1, width = 13, height = 7.4)

      dt_early <- dt1[timing == "Early break"]
      dt_late <- copy(dt1[timing == "Late break"])

      p2_early <- ggplot(dt_early, aes(x = Delta, y = CADD_adj, group = method, linetype = method, shape = method, colour = method)) +
        geom_line(linewidth = 0.7) +
        geom_point(size = 1.9) +
        facet_wrap(~ spec, ncol = 3, scales = "free") +
        labs(
          title = "Early break",
          x = "\u0394",
          y = "CADD"
        ) +
        scale_colour_manual(values = method_colours) +
        scale_y_continuous(breaks = scales::breaks_pretty(n = 3),
                           expand = expansion(mult = c(0.08, 0.10))) +
        scale_linetype_manual(values = method_linetypes) +
        scale_shape_manual(values = method_shapes) +
        mean_theme +
        theme(
          plot.title = element_text(size = 11, face = "bold"),
          text = element_text(size = 10),
          axis.text = element_text(size = 9),
          axis.title = element_text(size = 10),
          strip.text = element_text(size = 10),
          panel.spacing = grid::unit(0.7, "lines"),
          legend.position = "bottom"
        )

      p2_late <- ggplot(dt_late, aes(x = Delta, y = CADD_adj, group = method, linetype = method, shape = method, colour = method)) +
        geom_line(linewidth = 0.7) +
        geom_point(size = 1.9) +
        facet_wrap(~ spec, ncol = 3, scales = "free") +
        labs(
          title = "Late break",
          x = "\u0394",
          y = "CADD"
        ) +
        scale_colour_manual(values = method_colours) +
        scale_y_continuous(breaks = scales::breaks_pretty(n = 3),
                           expand = expansion(mult = c(0.08, 0.10))) +
        scale_linetype_manual(values = method_linetypes) +
        scale_shape_manual(values = method_shapes) +
        mean_theme +
        theme(
          plot.title = element_text(size = 11, face = "bold"),
          text = element_text(size = 10),
          axis.text = element_text(size = 9),
          axis.title = element_text(size = 10),
          strip.text = element_text(size = 10),
          panel.spacing = grid::unit(0.7, "lines"),
          legend.position = "bottom"
        )

      p2 <- (p2_early / p2_late) +
        plot_layout(guides = "collect") &
        theme(legend.position = "bottom")

      ggsave(filename = file.path(FIG_DIR, paste0("mean_ADD_adj_", tag, ".pdf")),
             plot = p2, width = 8, height = 7.4)

      p3 <- ggplot(dt1, aes(x = Delta, y = ARL_adj, group = method, linetype = method, shape = method)) +
        geom_line(linewidth = 0.7, colour = "black") +
        geom_point(size = 1.9, colour = "black") +
        facet_wrap(~ panel, ncol = 3, scales = "free") +
        labs(title = paste("ARL (size-adjusted) vs \u0394", title_suffix),
             x = "\u0394", y = "ARL (unconditional delay)") +
        scale_linetype_manual(values = method_linetypes) +
        scale_shape_manual(values = method_shapes) +
        mean_theme

      ggsave(filename = file.path(FIG_DIR, paste0("mean_ARL_adj_", tag, ".pdf")),
             plot = p3, width = 13, height = 7.4)
    }
  }
}

msg("Saved plots to: ", FIG_DIR)
msg("DONE.")
