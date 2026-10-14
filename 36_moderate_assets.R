# ============================================================
# 36_moderate_assets.R
# - Paper-ready appendix assets for the moderate-dimension experiment
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

RES_DIR <- file.path(ROOT, "results")
FIG_DIR <- file.path(ROOT, "figures")
if (!dir.exists(RES_DIR)) dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

ensure_pkgs(c("data.table", "ggplot2", "patchwork", "showtext"), auto_install = AUTO_INSTALL)
library(data.table)
library(ggplot2)
library(patchwork)
library(showtext)
showtext_auto()

path <- file.path(RES_DIR, "moderate_dim_sparse_results.csv")
if (!file.exists(path)) stop("Missing: ", path)
dt <- fread(path)

dt_size <- dt[exp == "size_adj"]
dt_pow <- dt[exp == "power_delay"]

fmt_num <- function(x) ifelse(is.na(x), "--", sprintf("%.3f", x))

method_greys <- c(RSMS = "black", SSMS = "black", CSMS = "black")
method_linetypes <- c(RSMS = "solid", SSMS = "longdash", CSMS = "dotdash")
shift_shapes <- c(dense = 16, sparse = 17)

make_size_table <- function() {
  lines <- c(
    "{\\color{blue}",
    "\\begin{table}[htbp]",
    "\\caption{Moderate-dimension size check under Monte Carlo calibrated critical values. The design uses $m=500$, $T=2$, and reports the empirical rejection rate under the null for $q\\in\\{5,8,12,20\\}$ at the 5\\% level.}",
    "\\label{table-moderate-size-blue}",
    "\\centering",
    "\\renewcommand{\\arraystretch}{1.05}",
    "\\begin{tabular}{ccccc}",
    "\\toprule",
    "$q$ & $\\gamma$ & R & S & C \\\\",
    "\\midrule"
  )

  for (q0 in sort(unique(dt_size$q))) {
    for (gamma0 in sort(unique(dt_size$gamma))) {
      row_vals <- sapply(c("RSMS", "SSMS", "CSMS"), function(method0) {
        sub <- dt_size[q == q0 & abs(gamma - gamma0) < 1e-12 & method == method0]
        if (nrow(sub) == 0L) return("--")
        fmt_num(sub$size_adj[1])
      })
      lines <- c(lines, sprintf("%d & %.2f & %s & %s & %s \\\\", q0, gamma0, row_vals[1], row_vals[2], row_vals[3]))
    }
    lines <- c(lines, "\\midrule")
  }

  lines[length(lines)] <- "\\bottomrule"
  lines <- c(lines, "\\end{tabular}", "\\end{table}", "}")
  out <- file.path(RES_DIR, "moderate_size_table_blue.tex")
  writeLines(lines, out)
  msg("Wrote: ", out)
}

plot_metric <- function(metric, ylab, file_stub) {
  dt_plot <- copy(dt_pow)
  dt_plot[, row_label := fifelse(
    abs(gamma) < 1e-12 & break_frac < 0.5, "g0_early",
    fifelse(
      abs(gamma) < 1e-12 & break_frac >= 0.5, "g0_late",
      fifelse(
        abs(gamma - 0.15) < 1e-12 & break_frac < 0.5, "g015_early",
        "g015_late"
      )
    )
  )]
  row_levels <- c("g0_early", "g0_late", "g015_early", "g015_late")
  dt_plot[, row_id := factor(row_label, levels = row_levels)]
  dt_plot[, q_lab := factor(paste0("q=", q), levels = paste0("q=", sort(unique(q))))]

  base_plot <- function(subdt, row_title, ylim = NULL) {
    p_row <- ggplot(
      subdt,
      aes(x = Delta, y = .data[[metric]], colour = method, linetype = method, shape = shift, group = interaction(method, shift))
    ) +
      geom_line(linewidth = 0.75) +
      geom_point(size = 2.1) +
      facet_wrap(~q_lab, nrow = 1, scales = "free_y") +
      scale_y_continuous(breaks = scales::breaks_pretty(n = 3), expand = expansion(mult = c(0.08, 0.10))) +
      scale_x_continuous(breaks = sort(unique(dt_pow$Delta))) +
      scale_colour_manual(values = method_greys) +
      scale_linetype_manual(values = method_linetypes) +
      scale_shape_manual(values = shift_shapes) +
      guides(colour = guide_legend(override.aes = list(shape = NA))) +
      labs(
        title = row_title,
        x = "\u0394",
        y = ylab,
        colour = "Method",
        linetype = "Method",
        shape = "Shift type"
      ) +
      theme_bw(base_size = 12) +
      theme(
        legend.position = "bottom",
        legend.key.width = grid::unit(1.4, "cm"),
        strip.text = element_text(size = 10),
        plot.title = element_text(size = 11, face = "bold"),
        panel.grid.minor = element_blank()
      )
    p_row
  }

  late_limits <- if (metric == "DP_post_adj") {
    list(c(0, 0.15), c(0, 0.10))
  } else {
    list(c(70, 150), c(70, 160))
  }
  row_specs <- list(
    list(key = "g0_early", title = "\u03b3 = 0, early break", ylim = NULL),
    list(key = "g0_late", title = "\u03b3 = 0, late break", ylim = late_limits[[1]]),
    list(key = "g015_early", title = "\u03b3 = 0.15, early break", ylim = NULL),
    list(key = "g015_late", title = "\u03b3 = 0.15, late break", ylim = late_limits[[2]])
  )

  row_plots <- lapply(row_specs, function(spec) {
    base_plot(
      dt_plot[row_id == spec$key],
      row_title = spec$title,
      ylim = spec$ylim
    )
  })

  p <- wrap_plots(row_plots, ncol = 1, guides = "collect") &
    theme(legend.position = "bottom")

  out <- file.path(FIG_DIR, paste0(file_stub, ".pdf"))
  ggsave(out, p, width = 10.8, height = 9.2)
  msg("Wrote: ", out)
}

make_size_table()
plot_metric("DP_post_adj", expression(DP[post]), "moderate_power_adj")
plot_metric("CADD_adj", "CADD", "moderate_ADD_adj")

msg("DONE.")
