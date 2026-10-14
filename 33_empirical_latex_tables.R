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

ensure_pkgs(c("dplyr"), auto_install = AUTO_INSTALL)

RES_DIR <- file.path(ROOT, "results")
in_path <- file.path(RES_DIR, "empirical_monitoring_results.csv")
if (!file.exists(in_path)) stop("Missing: ", in_path)

df <- read.csv(in_path, stringsAsFactors = FALSE)
df$ratio <- df$statistic / df$critical_value

format_date <- function(x) {
  ifelse(is.na(x) | !nzchar(x), "", format(as.Date(x), "%m-%d"))
}

format_gamma <- function(x) sprintf("%.2f", as.numeric(x))
format_ratio <- function(x) sprintf("%.2f", as.numeric(x))

all_methods <- c("RSMS", "SSMS", "CSMS")

key_df <- df |>
  dplyr::group_by(year, event, event_date, m, T, gamma) |>
  dplyr::summarise(
    has_all_methods = all(all_methods %in% method),
    all_three_reject = all(reject[match(all_methods, method)]),
    .groups = "drop"
  ) |>
  dplyr::filter(has_all_methods, all_three_reject) |>
  dplyr::arrange(year, m, T, gamma)

df_keep <- dplyr::inner_join(
  df,
  key_df[, c("year", "event", "event_date", "m", "T", "gamma")],
  by = c("year", "event", "event_date", "m", "T", "gamma")
)

format_cell <- function(block, method_name) {
  x <- block[block$method == method_name, ]
  if (nrow(x) == 0L) return("")
  paste0(format_ratio(x$ratio[[1]]), "/", format_date(x$first_detection_date[[1]]))
}

render_panel <- function(year_value, panel_title) {
  keys <- key_df[key_df$year == year_value, ]
  keys <- keys[order(keys$m, keys$T, keys$gamma), ]
  if (nrow(keys) == 0L) {
    return(c(
      sprintf("\\multicolumn{6}{c}{\\textbf{%s}} \\\\", panel_title),
      "\\midrule",
      "\\multicolumn{6}{c}{No configurations in which all three methods reject.} \\\\"
    ))
  }

  rows <- apply(keys, 1, function(row) {
    m_i <- as.numeric(row[["m"]])
    T_i <- as.numeric(row[["T"]])
    gamma_i <- as.numeric(row[["gamma"]])
    block <- df_keep[
      df_keep$year == year_value &
        df_keep$m == m_i &
        df_keep$T == T_i &
        abs(df_keep$gamma - gamma_i) < 1e-12,
    ]
    paste(
      m_i,
      T_i,
      format_gamma(gamma_i),
      format_cell(block, "RSMS"),
      format_cell(block, "SSMS"),
      format_cell(block, "CSMS"),
      sep = " & "
    )
  })

  c(
    sprintf("\\multicolumn{6}{c}{\\textbf{%s}} \\\\", panel_title),
    "\\cmidrule(lr){1-6}",
    "$m$ & $T$ & $\\gamma$ & R & S & C \\\\",
    "\\midrule",
    paste0(rows, " \\\\")
  )
}

lines <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{USD/GBP: statistic-to-boundary ratios and first detection dates for RSMS (R), SSMS (S), and CSMS (C)}",
  "\\label{tab:exchange-rate-rejection}",
  "\\begin{threeparttable}",
  "\\scriptsize",
  "\\setlength{\\tabcolsep}{2pt}",
  "\\renewcommand{\\arraystretch}{1.03}",
  "\\begin{tabular}{@{}ccc ccc@{}}",
  "\\toprule",
  render_panel(2016, "Panel A: 2016 Brexit referendum"),
  "\\midrule",
  render_panel(2020, "Panel B: 2020 COVID-19 shock"),
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{threeparttable}",
  "\\begin{flushleft}",
  "\\footnotesize",
  "\\emph{Notes:} Rows are configurations in which all three methods reject \\(H_0\\). Each method cell reports the maximal monitoring statistic divided by its asymptotic critical value, followed by the first detection date. Ratios above one indicate rejection; reporting ratios rather than raw statistics puts the three methods on a comparable boundary scale. Scheme abbreviations: R = RSMS, S = SSMS, C = CSMS. Dates are in month--day format within the corresponding calendar year.",
  "\\end{flushleft}",
  "\\end{table}"
)

out_path <- file.path(RES_DIR, "empirical_reject_ratio_table.tex")
writeLines(lines, out_path)

msg("Saved: ", out_path)
msg("DONE.")
