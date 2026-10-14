# ============================================================
# 35_mean_alt_latex_tables.R
# - Build paper-ready LaTeX tables for size-adjusted power
#   in the baseline mean-monitoring experiments
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
if (!dir.exists(RES_DIR)) dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)
ensure_pkgs(c("data.table"), auto_install = AUTO_INSTALL)

library(data.table)

M_MAIN <- env_int("JBES_M_MAIN", 250L)
DGP_ORDER <- c("DGP1", "DGP2", "DGP3", "DGP4")
METHOD_ORDER <- c("RSMS", "SSMS", "CSMS")
METHOD_LABEL <- c(RSMS = "R", SSMS = "S", CSMS = "C")
BREAK_LABEL <- c(`0.1` = "$k^*=0.1mT$", `0.8` = "$k^*=0.8mT$")
EARLY_DELTA_GRID <- env_num_vec("JBES_DELTA_GRID_EARLY", c(0.25, 0.50, 0.75, 1.00))
LATE_DELTA_GRID <- env_num_vec("JBES_DELTA_GRID_LATE", c(0.50, 1.00, 1.50, 2.00))

if (length(EARLY_DELTA_GRID) != length(LATE_DELTA_GRID)) {
  stop("EARLY_DELTA_GRID and LATE_DELTA_GRID must have the same length for table pairing.")
}

alt_path <- file.path(RES_DIR, "mean_alt_performance.csv")
if (!file.exists(alt_path)) stop("Missing: ", alt_path)
dt <- fread(alt_path)

fmt_pct <- function(x) {
  ifelse(is.na(x), "--", sprintf("%.1f", 100 * x))
}

value_lookup <- function(dgp0, bf0, T0, Delta0, method0, chg0, gamma0) {
  sub <- dt[
    dgp == dgp0 &
      m == M_MAIN &
      abs(gamma - gamma0) < 1e-12 &
      abs(break_frac - bf0) < 1e-12 &
      T == T0 &
      abs(Delta - Delta0) < 1e-12 &
      change_type == chg0 &
      method == method0
  ]
  if (nrow(sub) == 0L) return("--")
  fmt_pct(sub$DP_post_adj[1])
}

build_rows <- function(gamma0, chg0) {
  rows <- character(0)
  for (T0 in c(1, 2, 5)) {
    for (ii in seq_along(EARLY_DELTA_GRID)) {
      DeltaE <- EARLY_DELTA_GRID[ii]
      DeltaL <- LATE_DELTA_GRID[ii]
      early_vals <- character(0)
      late_vals <- character(0)
      for (dgp0 in DGP_ORDER) {
        for (method0 in METHOD_ORDER) {
          early_vals <- c(early_vals, value_lookup(dgp0, 0.1, T0, DeltaE, method0, chg0, gamma0))
          late_vals <- c(late_vals, value_lookup(dgp0, 0.8, T0, DeltaL, method0, chg0, gamma0))
        }
      }
      row <- paste(
        paste0("$T=", T0, "$"),
        sprintf("%.2f", DeltaE),
        paste(early_vals, collapse = " & "),
        "",
        sprintf("%.2f", DeltaL),
        paste(late_vals, collapse = " & "),
        sep = " & "
      )
      rows <- c(rows, paste0(row, " \\\\"))
    }
    if (T0 < 5) rows <- c(rows, "\\hline")
  }
  rows
}

build_rows_side <- function(gamma0, chg0, bf0, delta_grid) {
  rows <- character(0)
  for (T0 in c(1, 2, 5)) {
    for (Delta0 in delta_grid) {
      vals <- character(0)
      for (dgp0 in DGP_ORDER) {
        for (method0 in METHOD_ORDER) {
          vals <- c(vals, value_lookup(dgp0, bf0, T0, Delta0, method0, chg0, gamma0))
        }
      }
      row <- paste(
        paste0("$T=", T0, "$"),
        sprintf("%.2f", Delta0),
        paste(vals, collapse = " & "),
        sep = " & "
      )
      rows <- c(rows, paste0(row, " \\\\"))
    }
    if (T0 < 5) rows <- c(rows, "\\hline")
  }
  rows
}

power_col_spec <- function() {
  paste0(
    "@{}",
    ">{\\centering\\arraybackslash}m{0.075\\textwidth}",
    ">{\\centering\\arraybackslash}m{0.060\\textwidth}",
    paste(rep(">{\\centering\\arraybackslash}m{0.049\\textwidth}", 12), collapse = ""),
    "@{}"
  )
}

build_block <- function(title, bf0, delta_grid, gamma0, change_types = c("abrupt", "smooth"),
                        top_rule = TRUE, bottom_rule = TRUE) {
  tab_spec <- power_col_spec()
  out <- c(
    paste0("\\begin{tabular}{", tab_spec, "}")
  )
  if (top_rule) out <- c(out, "\\toprule")
  out <- c(
    out,
    paste0("\\multicolumn{14}{c}{", title, "} \\\\"),
    "\\midrule",
    "\\multicolumn{1}{c}{} & $\\Delta$ & \\multicolumn{3}{c}{DGP1} & \\multicolumn{3}{c}{DGP2} & \\multicolumn{3}{c}{DGP3} & \\multicolumn{3}{c}{DGP4} \\\\",
    "\\cmidrule(lr){3-5} \\cmidrule(lr){6-8} \\cmidrule(lr){9-11} \\cmidrule(lr){12-14}",
    "$T$ & $\\Delta$ & R & S & C & R & S & C & R & S & C & R & S & C \\\\"
  )

  first_block <- TRUE
  for (chg0 in change_types) {
    if (length(change_types) == 1L) {
      out <- c(out, "\\midrule", build_rows_side(gamma0, chg0, bf0, delta_grid))
    } else {
      out <- c(
        out,
        "\\midrule",
        paste0("\\multicolumn{14}{c}{", if (chg0 == "abrupt") "Abrupt change" else "Smooth change", "} \\\\"),
        "\\midrule",
        build_rows_side(gamma0, chg0, bf0, delta_grid)
      )
    }
    first_block <- FALSE
  }

  if (bottom_rule) out <- c(out, "\\bottomrule")
  c(out, "\\end{tabular}")
}

write_single_table <- function(gamma0, break_tag, break_title, bf0, delta_grid,
                               change_types = c("abrupt", "smooth"),
                               file_suffix = "") {
  gamma_lab <- if (abs(gamma0) < 1e-12) "0" else "0.15"
  gamma_stub <- gsub("\\.", "p", gamma_lab)
  if (length(change_types) == 1L) {
    ct_lab <- if (change_types == "abrupt") "abrupt" else "smooth"
    caption_prefix <- if (change_types == "abrupt") "Abrupt-change" else "Smooth-change"
    suffix_piece <- paste0("_", ct_lab)
    label_suffix <- paste0("-", ct_lab)
  } else {
    caption_prefix <- "Size-adjusted"
    suffix_piece <- ""
    label_suffix <- ""
  }
  if (nzchar(file_suffix)) {
    suffix_piece <- paste0(suffix_piece, "_", file_suffix)
    label_suffix <- paste0(label_suffix, "-", gsub("_", "-", file_suffix))
  }
  out_file <- file.path(
    RES_DIR,
    paste0("mean_power_adj_m", M_MAIN, "_g", gamma_stub, "_", break_tag, suffix_piece, "_blue.tex")
  )
  label_name <- paste0("table-mean-poweradj-m", M_MAIN, "-g", gsub("\\.", "", gamma_lab), "-", break_tag, label_suffix, "-blue")
  caption_text <- paste0(
    caption_prefix,
    " detection probabilities (\\%) for detecting a mean shift using RSMS (R), SSMS (S), and CSMS (C). Training size $m=",
    M_MAIN,
    "$; $\\gamma=",
    gamma_lab,
    "$; ",
    break_title,
    ". The ",
    break_tag,
    "-break grid uses $\\Delta\\in\\{",
    paste(sprintf("%.2f", delta_grid), collapse = ","),
    "\\}$."
  )

  lines <- c(
    "{\\color{blue}",
    "\\begin{table}[htbp]",
    paste0("\\caption{", caption_text, "}"),
    paste0("\\label{", label_name, "}"),
    "\\centering",
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{3pt}",
    "\\renewcommand{\\arraystretch}{1.05}",
    "\\begin{adjustbox}{max width=\\textwidth,center}",
    build_block(
      if (break_tag == "early") "Early break ($k^*=0.1mT$)" else "Late break ($k^*=0.8mT$)",
      bf0,
      delta_grid,
      gamma0,
      change_types = change_types
    ),
    "\\end{adjustbox}",
    "\\end{table}",
    "}"
  )

  writeLines(lines, out_file)
  msg("Wrote: ", out_file)
}

write_combined_smooth_table <- function() {
  gamma_lab <- "0.15"
  out_file <- file.path(
    RES_DIR,
    paste0("mean_power_adj_m", M_MAIN, "_g0p15_smooth_blue.tex")
  )
  caption_text <- paste0(
    "Smooth-change detection probabilities (\\%) for detecting a mean shift using RSMS (R), SSMS (S), and CSMS (C). ",
    "Training size $m=", M_MAIN, "$; $\\gamma=0.15$. ",
    "Panel A reports early breaks with $\\Delta\\in\\{",
    paste(sprintf("%.2f", EARLY_DELTA_GRID), collapse = ","),
    "\\}$, and Panel B reports late breaks with $\\Delta\\in\\{",
    paste(sprintf("%.2f", LATE_DELTA_GRID), collapse = ","),
    "\\}$."
  )

  lines <- c(
    "{\\color{blue}",
    "\\begin{table}[htbp]",
    paste0("\\caption{", caption_text, "}"),
    "\\label{table-mean-poweradj-m250-g015-smooth-blue}",
    "\\centering",
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{3pt}",
    "\\renewcommand{\\arraystretch}{1.05}",
    "\\begin{adjustbox}{max totalsize={\\textwidth}{0.82\\textheight},center}",
    "\\begin{minipage}{\\textwidth}",
    "\\centering",
    build_block(
      "Panel A. Early break ($k^*=0.1mT$)",
      0.1,
      EARLY_DELTA_GRID,
      0.15,
      change_types = "smooth",
      bottom_rule = FALSE
    ),
    "\\par\\smallskip",
    build_block(
      "Panel B. Late break ($k^*=0.8mT$)",
      0.8,
      LATE_DELTA_GRID,
      0.15,
      change_types = "smooth",
      top_rule = TRUE
    ),
    "\\end{minipage}",
    "\\end{adjustbox}",
    "\\end{table}",
    "}"
  )

  writeLines(lines, out_file)
  msg("Wrote: ", out_file)
}

for (gamma0 in c(0, 0.15)) {
  write_single_table(gamma0, "early", "early break ($k^*=0.1mT$)", 0.1, EARLY_DELTA_GRID)
  write_single_table(gamma0, "late", "late break ($k^*=0.8mT$)", 0.8, LATE_DELTA_GRID)
}

# Main-text / appendix companion tables for gamma = 0.15
write_single_table(0.15, "early", "early break ($k^*=0.1mT$)", 0.1, EARLY_DELTA_GRID,
  change_types = "abrupt"
)
write_single_table(0.15, "late", "late break ($k^*=0.8mT$)", 0.8, LATE_DELTA_GRID,
  change_types = "abrupt"
)
write_single_table(0.15, "early", "early break ($k^*=0.1mT$)", 0.1, EARLY_DELTA_GRID,
  change_types = "smooth"
)
write_single_table(0.15, "late", "late break ($k^*=0.8mT$)", 0.8, LATE_DELTA_GRID,
  change_types = "smooth"
)
write_combined_smooth_table()

msg("DONE.")
