# ============================================================
# 38_contamination_assets.R
# - Paper-ready appendix assets for contaminated-training experiments
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

path <- file.path(RES_DIR, "training_contamination_results.csv")
if (!file.exists(path)) stop("Missing: ", path)
dt <- fread(path)

dt_size <- dt[exp == "size_contaminated_training"]
dt_pow <- dt[exp == "power_with_contaminated_training"]
plot_B_grid <- sort(unique(dt_pow$B))

fmt_num <- function(x) ifelse(is.na(x), "--", sprintf("%.3f", x))
method_linetypes <- c(RSMS = "solid", SSMS = "longdash", CSMS = "dotdash")
method_shapes <- c(RSMS = 16, SSMS = 2, CSMS = 0)

make_size_table <- function() {
  lines <- c(
    "{\\color{blue}",
    "\\begin{table}[htbp]",
    "\\caption{Size under contaminated training when asymptotic critical values from the clean design are used. Regime ``mild'' corresponds to $\\ell_m=o(\\sqrt{m})$, while ``fixed'' denotes a fixed-fraction contamination stress test. The column $B$ records the magnitude of the contamination shift applied to the affected training observations.}",
    "\\label{table-contamination-size-blue}",
    "\\centering",
    "\\renewcommand{\\arraystretch}{1.05}",
    "\\resizebox{\\textwidth}{!}{%",
    "\\begin{tabular}{cccccc}",
    "\\toprule",
    "Regime & Contamination magnitude $B$ & $\\gamma$ & R & S & C \\\\",
    "\\midrule"
  )

  regimes <- unique(dt_size$regime)
  regimes <- regimes[order(regimes)]
  for (reg in regimes) {
    sub_reg <- dt_size[regime == reg]
    for (B0 in sort(unique(sub_reg$B))) {
      for (gamma0 in sort(unique(sub_reg$gamma))) {
        vals <- sapply(c("RSMS", "SSMS", "CSMS"), function(method0) {
          sub <- sub_reg[abs(B - B0) < 1e-12 & abs(gamma - gamma0) < 1e-12 & method == method0]
          if (nrow(sub) == 0L) return("--")
          fmt_num(sub$size[1])
        })
        lines <- c(lines, sprintf("%s & %.2f & %.2f & %s & %s & %s \\\\", reg, B0, gamma0, vals[1], vals[2], vals[3]))
      }
    }
    lines <- c(lines, "\\midrule")
  }

  lines[length(lines)] <- "\\bottomrule"
  lines <- c(lines, "\\end{tabular}", "}", "\\end{table}", "}")
  out <- file.path(RES_DIR, "contamination_size_table_blue.tex")
  writeLines(lines, out)
  msg("Wrote: ", out)
}

plot_metric <- function(metric, ylab, file_stub, zoom_late = FALSE) {
  regime_vals <- unique(dt_pow$regime)
  regime_vals <- regime_vals[order(regime_vals)]
  gamma_vals <- sort(unique(dt_pow$gamma))
  break_vals <- sort(unique(dt_pow$break_frac))

  make_panel <- function(regime0, gamma0, break0, zoom_y = FALSE) {
    panel_dt <- dt_pow[
      regime == regime0 &
        abs(gamma - gamma0) < 1e-12 &
        abs(break_frac - break0) < 1e-12 &
        B %in% plot_B_grid
    ][order(B, method)]
    p <- ggplot(
      panel_dt,
      aes(x = B, y = .data[[metric]], linetype = method, shape = method, group = method)
    ) +
      geom_line(linewidth = 0.75, colour = "black", na.rm = TRUE) +
      geom_point(size = 1.5, colour = "black", na.rm = TRUE) +
      scale_x_continuous(breaks = seq(0, 1, by = 0.25)) +
      scale_y_continuous(breaks = scales::breaks_pretty(n = 3), expand = expansion(mult = c(0.08, 0.10))) +
      scale_linetype_manual(values = method_linetypes) +
      scale_shape_manual(values = method_shapes) +
      labs(
        title = bquote(.(regime0) * "," ~ gamma == .(gamma0) * "," ~ "break" == .(break0)),
        x = expression("Contamination magnitude " * B),
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
  for (regime0 in regime_vals) {
    for (gamma0 in gamma_vals) {
      for (break0 in break_vals) {
        zoom_y <- zoom_late && break0 >= 0.8
        plot_list[[length(plot_list) + 1L]] <- make_panel(regime0, gamma0, break0, zoom_y)
      }
    }
  }

  p <- wrap_plots(plotlist = plot_list, ncol = 2, byrow = TRUE, guides = "collect") &
    theme(legend.position = "bottom")

  out <- file.path(FIG_DIR, paste0(file_stub, ".pdf"))
  ggsave(out, p, width = 8, height = 8.2)
  msg("Wrote: ", out)
}

make_size_table()
plot_metric("DP_post_adj", expression(DP[post]), "contamination_power_adj", zoom_late = TRUE)
plot_metric("CADD_adj", "CADD", "contamination_ADD_adj", zoom_late = TRUE)

msg("DONE.")
