rm(list = ls())

AUTO_INSTALL <- TRUE

bootstrap_find_root <- function() {
  env_root <- Sys.getenv("JBES_PROJECT_ROOT", unset = "")
  if (nzchar(env_root)) {
    return(normalizePath(env_root, winslash = "/", mustWork = FALSE))
  }

  candidates <- unique(c(getwd()))
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    ctx_path <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(ctx_path)) candidates <- unique(c(dirname(ctx_path), candidates))
  }

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

ensure_pkgs(c("readr"), auto_install = AUTO_INSTALL)

RES_DIR <- file.path(ROOT, "results")
path <- file.path(RES_DIR, "par_alt_performance_grid.csv")
if (!file.exists(path)) stop("Missing: ", path)

df <- read.csv(path, stringsAsFactors = FALSE)
df <- subset(df, m == 250 & T == 1 & change_type == "abrupt")

alpha_vals <- c(0.20, 0.25, 0.30, 0.35)
lambda_vals <- c(2.0, 2.1, 2.2, 2.3, 2.4, 2.5)
method_order <- c("RSMS", "SSMS", "CSMS")
break_panels <- list("Panel A: Early break" = 0.1, "Panel B: Late break" = 0.8)

fmt_power <- function(x) sprintf("%.0f", 100 * x)
fmt_add <- function(x) sprintf("%.1f", x)

cell_entry <- function(gamma0, break0, lambda0, alpha0) {
  hit <- subset(
    df,
    abs(gamma - gamma0) < 1e-12 &
      abs(break_frac - break0) < 1e-12 &
      abs(lambda2 - lambda0) < 1e-12 &
      abs(alpha2 - alpha0) < 1e-12
  )
  if (nrow(hit) != 3L) stop("Missing or duplicate PAR grid entry.")
  hit <- hit[match(method_order, hit$method), ]

  power_line <- paste(fmt_power(hit$DP_post_adj), collapse = "/")
  add_line <- paste(fmt_add(hit$CADD_adj), collapse = "/")
  sprintf("\\shortstack[c]{%s\\\\%s}", power_line, add_line)
}

make_panel <- function(gamma0, panel_title, break0) {
  col_spec <- paste(rep(">{\\centering\\arraybackslash}m{0.17\\linewidth}", 5), collapse = "")
  lines <- c(
    sprintf("\\begin{tabular}{%s}", col_spec),
    "\\toprule",
    sprintf("\\multicolumn{5}{c}{%s} \\\\", panel_title),
    "\\midrule",
    "$\\lambda_2 \\backslash \\alpha_2$ & 0.20 & 0.25 & 0.30 & 0.35 \\\\",
    "\\midrule"
  )

  for (lambda0 in lambda_vals) {
    cells <- vapply(alpha_vals, function(alpha0) cell_entry(gamma0, break0, lambda0, alpha0), character(1))
    row <- paste(c(sprintf("%.1f", lambda0), cells), collapse = " & ")
    lines <- c(lines, paste0(row, " \\\\"))
  }

  c(lines, "\\bottomrule", "\\end{tabular}")
}

write_gamma_table <- function(gamma0, file_name, label_name) {
  gamma_label <- if (abs(gamma0) < 1e-12) "0" else "0.15"
  caption_text <- paste0(
    "PAR(1) grid experiment with $m=250$, $T=1$, and $\\gamma=",
    gamma_label,
    "$. Separate panels report early and late abrupt breaks over the $(\\lambda_2,\\alpha_2)$ grid. ",
    "Within each cell, the first line gives the size-adjusted post-break detection probability (\\%) and the second line gives CADD, both ordered as R/S/C = RSMS/SSMS/CSMS."
  )

  lines <- c(
    "{\\color{blue}",
    "\\begin{table}[p]",
    sprintf("\\caption{%s}", caption_text),
    sprintf("\\label{%s}", label_name),
    "\\centering",
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{2pt}",
    "\\renewcommand{\\arraystretch}{1.05}"
  )

  first_panel <- TRUE
  for (panel_name in names(break_panels)) {
    if (!first_panel) {
      lines <- c(lines, "\\par\\medskip")
    }
    lines <- c(lines, make_panel(gamma0, panel_name, break_panels[[panel_name]]))
    first_panel <- FALSE
  }

  lines <- c(
    lines,
    "\\par\\medskip",
    "\\begin{minipage}{0.94\\linewidth}",
    "\\footnotesize \\emph{Notes:} In each cell, the first line reports the size-adjusted post-break detection probability (in \\%) and the second line reports the conditional average detection delay (CADD). Entries are ordered as R/S/C = RSMS/SSMS/CSMS.",
    "\\end{minipage}",
    "\\end{table}",
    "}"
  )

  writeLines(lines, file.path(RES_DIR, file_name))
}

write_gamma_table(0, "par_grid_table_g0_blue.tex", "table-par-grid-g0-blue")
write_gamma_table(0.15, "par_grid_table_g015_blue.tex", "table-par-grid-g015-blue")

msg("Saved PAR appendix tables to: ", RES_DIR)
