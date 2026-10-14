# ============================================================
# 39_critical_value_assets.R
# - Generate LaTeX row fragments for critical-value tables
#   using finite- and open-end simulated asymptotic values
# ============================================================

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
ensure_pkgs("data.table", auto_install = AUTO_INSTALL)
library(data.table)

RES_DIR <- file.path(ROOT, "results")
finite_path <- file.path(RES_DIR, "affine_finite_critical_values.csv")
open_path <- file.path(RES_DIR, "affine_open_end_critical_values.csv")
if (!file.exists(finite_path)) stop("Missing: ", finite_path)
if (!file.exists(open_path)) stop("Missing: ", open_path)

finite_dt <- fread(finite_path)
open_dt <- fread(open_path)

T_GRID <- c(1, 2, 5, 10)
GAMMA_GRID <- c(0, 0.15)
ALPHA_GRID <- c(0.05, 0.10)
METHODS <- c("RSMS", "CSMS", "SSMS")

q_all <- sort(unique(finite_dt$d))
q_all <- q_all[q_all %in% sort(unique(open_dt$d))]

fmt1 <- function(x) formatC(as.numeric(x), format = "f", digits = 1)
blue_cell <- function(x) sprintf("{\\color{blue}%s}", x)

get_finite <- function(method_name, d_val, T_val, gamma_val, alpha_val) {
  hit <- finite_dt[
    method == method_name &
      d == d_val &
      abs(T - T_val) < 1e-12 &
      abs(gamma - gamma_val) < 1e-12 &
      abs(alpha - alpha_val) < 1e-12
  ]
  if (nrow(hit) != 1L) stop("Missing finite critical value for ", method_name, ", d=", d_val, ", T=", T_val)
  blue_cell(fmt1(hit$critical_value[[1]]))
}

get_open <- function(method_name, d_val, gamma_val, alpha_val) {
  hit <- open_dt[
    method == method_name &
      d == d_val &
      abs(gamma - gamma_val) < 1e-12 &
      abs(alpha - alpha_val) < 1e-12
  ]
  if (nrow(hit) != 1L) stop("Missing open-end critical value for ", method_name, ", d=", d_val)
  blue_cell(fmt1(hit$critical_value[[1]]))
}

make_row <- function(method, d) {
  vals <- c(blue_cell(as.character(d)))
  for (gamma in GAMMA_GRID) {
    for (T in T_GRID) {
      for (alpha in ALPHA_GRID) {
        vals <- c(vals, get_finite(method, d, T, gamma, alpha))
      }
    }
    for (alpha in ALPHA_GRID) {
      vals <- c(vals, get_open(method, d, gamma, alpha))
    }
  }
  sprintf("%s \\tabularnewline", paste(vals, collapse = " & "))
}

write_rows <- function(method, filename) {
  lines <- vapply(q_all, function(d) make_row(method, d), character(1))
  out <- file.path(RES_DIR, filename)
  writeLines(lines, out)
  msg("Saved: ", out)
}

write_full_table <- function(lines, filename, caption, label, note_lines = NULL, placement = "[htbp]") {
  tab_lines <- c(
    sprintf("\\begin{table}%s", placement),
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    "\\centering",
    "\\setlength{\\tabcolsep}{4pt}",
    "\\renewcommand{\\arraystretch}{1.05}",
    "\\resizebox{\\textwidth}{!}{",
    "\\begin{tabular}{c cc cc cc cc cc cc cc cc cc cc}",
    "\\toprule",
    "\\multicolumn{1}{c}{} & \\multicolumn{10}{c}{$\\gamma=0$} & \\multicolumn{10}{c}{$\\gamma=0.15$} \\\\",
    "\\cmidrule(lr){2-11} \\cmidrule(lr){12-21}",
    " & \\multicolumn{2}{c}{$T=1$} & \\multicolumn{2}{c}{$T=2$} & \\multicolumn{2}{c}{$T=5$} & \\multicolumn{2}{c}{$T=10$} & \\multicolumn{2}{c}{{\\color{blue}$T=\\infty$}}",
    " & \\multicolumn{2}{c}{$T=1$} & \\multicolumn{2}{c}{$T=2$} & \\multicolumn{2}{c}{$T=5$} & \\multicolumn{2}{c}{$T=10$} & \\multicolumn{2}{c}{{\\color{blue}$T=\\infty$}} \\\\",
    "\\midrule",
    "{\\color{blue}$q\\backslash\\alpha$} & 5\\% & 10\\% & 5\\% & 10\\% & 5\\% & 10\\% & 5\\% & 10\\% & 5\\% & 10\\%",
    "                   & 5\\% & 10\\% & 5\\% & 10\\% & 5\\% & 10\\% & 5\\% & 10\\% & 5\\% & 10\\% \\\\",
    "\\midrule",
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    "}"
  )
  if (!is.null(note_lines) && length(note_lines) > 0L) {
    tab_lines <- c(tab_lines, note_lines)
  }
  tab_lines <- c(tab_lines, "\\end{table}")
  out <- file.path(RES_DIR, filename)
  writeLines(tab_lines, out)
  msg("Saved: ", out)
}

rsms_lines <- vapply(q_all, function(d) make_row("RSMS", d), character(1))
csms_lines <- vapply(q_all, function(d) make_row("CSMS", d), character(1))
ssms_lines <- vapply(q_all, function(d) make_row("SSMS", d), character(1))
q_main <- c(1L, 2L, 5L, 10L, 15L, 20L)
rsms_lines_main <- vapply(q_main, function(d) make_row("RSMS", d), character(1))

write_rows("RSMS", "rsms_cv_rows_blue.tex")
write_rows("CSMS", "csms_cv_rows_blue.tex")
write_rows("SSMS", "ssms_cv_rows_blue.tex")

write_full_table(
  rsms_lines_main,
  "rsms_cv_table_main_blue.tex",
  "Simulated asymptotic critical values for the RSMS at representative dimensions with $T = 1, 2, 5, 10, {\\color{blue}\\infty}$ and $\\gamma = 0$ or $\\gamma = 0.15$",
  "table:cv-RSMS",
  note_lines = c(
    "\\begin{flushleft}",
    "\\footnotesize",
    "  \\emph{Notes:} Table~\\ref{table:cv-RSMS} reports representative values for $q\\in\\{1,2,5,10,15,20\\}$. {\\color{blue}The complete grid for $q=1,\\ldots,20$ is reported in Appendix Table~\\ref{tab:critical_values_RSMS_full}.} \\textcolor{darkred}{Asymptotic decision boundaries $c_R$ are computed from 10{,}000 Monte Carlo replications. For finite horizons \\(T\\in\\{1,2,5,10\\}\\), Brownian paths are simulated on a regular grid with 400 increments per unit interval over \\([0,T+1]\\). For the open-end case \\(T=\\infty\\), after the transformation \\(u=s/(1+s)\\), the monitoring numerator is simulated as a standard Brownian motion on \\(u\\in(0,1]\\), independently of the Brownian bridge used to construct the training normalizer, on a grid of 5{,}000 points.}",
    "\\end{flushleft}"
  )
)

write_full_table(
  rsms_lines,
  "rsms_cv_table_appendix_blue.tex",
  "Simulated asymptotic critical values for the RSMS on the full dimension grid $q=1,\\ldots,20$ with $T = 1, 2, 5, 10, {\\color{blue}\\infty}$ and $\\gamma = 0$ or $\\gamma = 0.15$",
  "tab:critical_values_RSMS_full",
  note_lines = c(
    "\\begin{flushleft}",
    "\\footnotesize",
    "  \\emph{Notes:} \\textcolor{darkred}{Asymptotic decision boundaries $c_R$ are computed from 10{,}000 Monte Carlo replications. For finite horizons \\(T\\in\\{1,2,5,10\\}\\), Brownian paths are simulated on a regular grid with 400 increments per unit interval over \\([0,T+1]\\). For the open-end case \\(T=\\infty\\), after the transformation \\(u=s/(1+s)\\), the monitoring numerator is simulated as a standard Brownian motion on \\(u\\in(0,1]\\), independently of the Brownian bridge used to construct the training normalizer, on a grid of 5{,}000 points.}",
    "\\end{flushleft}"
  )
)

write_full_table(
  csms_lines,
  "csms_cv_table_blue.tex",
  "Simulated asymptotic critical values for the CSMS under horizons ($T = 1, 2, 5, 10, {\\color{blue}\\infty}$) and sensitivity parameters ($\\gamma = 0$ or $0.15$).",
  "tab:critical_values_CSMS",
  placement = "[H]"
)

write_full_table(
  ssms_lines,
  "ssms_cv_table_blue.tex",
  "Simulated asymptotic critical values for the SSMS with $T = 1, 2, 5, 10, {\\color{blue}\\infty}$ and $\\gamma = 0$ or $\\gamma = 0.15$",
  "tab:critical_values_SSMS",
  placement = "[hptb]"
)

msg("DONE.")
