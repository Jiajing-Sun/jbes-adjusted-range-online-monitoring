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

ensure_pkgs(c("dplyr", "tidyr", "readr"), auto_install = AUTO_INSTALL)

RES_DIR <- file.path(ROOT, "results")

mean_path <- file.path(RES_DIR, "mean_null_size.csv")
par_path <- file.path(RES_DIR, "par_null_size.csv")
if (!file.exists(mean_path)) stop("Missing: ", mean_path)
if (!file.exists(par_path)) stop("Missing: ", par_path)

mean_df <- read.csv(mean_path, stringsAsFactors = FALSE)
par_df <- read.csv(par_path, stringsAsFactors = FALSE)

method_map <- c(RSMS = "R", SSMS = "S", CSMS = "C")
t_grid <- c(1, 2, 5, 10)
alpha_grid <- c(0.05, 0.10)
gamma_grid <- c(0, 0.15)

fmt <- function(x) sprintf("%.1f", 100 * as.numeric(x))

mean_value <- function(df, dgp0, m0, gamma0, T0, method0, alpha0) {
  hit <- subset(df,
                dgp == dgp0 &
                  m == m0 &
                  abs(gamma - gamma0) < 1e-12 &
                  T == T0 &
                  method == method0 &
                  abs(alpha - alpha0) < 1e-12)
  if (nrow(hit) != 1L) stop("Missing or duplicate mean-null entry.")
  fmt(hit$size[1])
}

par_value <- function(df, m0, gamma0, T0, method0, alpha0) {
  hit <- subset(df,
                m == m0 &
                  abs(gamma - gamma0) < 1e-12 &
                  T == T0 &
                  method == method0 &
                  abs(alpha - alpha0) < 1e-12)
  if (nrow(hit) != 1L) stop("Missing or duplicate PAR-null entry.")
  fmt(hit$size[1])
}

make_mean_row <- function(df, dgp, m, gamma, t_subset, dgp_label = "") {
  vals <- c()
  for (T in t_subset) {
    for (method in c("RSMS", "SSMS", "CSMS")) {
      for (alpha in alpha_grid) {
        vals <- c(vals, mean_value(df, dgp, m, gamma, T, method, alpha))
      }
    }
  }
  paste(c(dgp_label, paste0("$m=", m, "$"), vals), collapse = " & ")
}

make_par_row <- function(df, m, gamma) {
  vals <- c()
  for (T in t_grid) {
    for (method in c("RSMS", "SSMS", "CSMS")) {
      for (alpha in alpha_grid) {
        vals <- c(vals, par_value(df, m, gamma, T, method, alpha))
      }
    }
  }
  paste(c(paste0("$m=", m, "$"), vals), collapse = " & ")
}

null_col_spec <- function() {
  paste0(
    "@{}",
    ">{\\raggedright\\arraybackslash}m{0.080\\textwidth}",
    ">{\\centering\\arraybackslash}m{0.065\\textwidth}",
    paste(rep(">{\\centering\\arraybackslash}m{0.047\\textwidth}", 12), collapse = ""),
    "@{}"
  )
}

build_mean_null_block <- function(t_subset, title) {
  lines <- c(
    paste0("\\begin{tabular}{", null_col_spec(), "}"),
    "\\toprule",
    paste0("\\multicolumn{14}{c}{", title, "} \\\\"),
    "\\midrule",
    "&  & \\multicolumn{2}{c}{R} & \\multicolumn{2}{c}{S} & \\multicolumn{2}{c}{C} & \\multicolumn{2}{c}{R} & \\multicolumn{2}{c}{S} & \\multicolumn{2}{c}{C} \\\\",
    paste0(
      "&  & \\multicolumn{2}{c}{$T=", t_subset[1], "$} & \\multicolumn{2}{c}{$T=", t_subset[1], "$} & \\multicolumn{2}{c}{$T=", t_subset[1], "$} & ",
      "\\multicolumn{2}{c}{$T=", t_subset[2], "$} & \\multicolumn{2}{c}{$T=", t_subset[2], "$} & \\multicolumn{2}{c}{$T=", t_subset[2], "$} \\\\"
    ),
    "\\cmidrule(lr){3-8} \\cmidrule(lr){9-14}",
    "&  & 5\\% & 10\\% & 5\\% & 10\\% & 5\\% & 10\\% & 5\\% & 10\\% & 5\\% & 10\\% & 5\\% & 10\\% \\\\",
    "\\midrule"
  )
  for (g_idx in seq_along(gamma_grid)) {
    gamma <- gamma_grid[g_idx]
    if (g_idx > 1L) lines <- c(lines, "\\midrule")
    lines <- c(lines, sprintf("&  & \\multicolumn{12}{c}{{\\color{blue}$\\gamma = %.2f$}} \\\\", gamma), "\\midrule")
    dgp_vals <- unique(mean_df$dgp)
    for (dgp in dgp_vals) {
      m_vals <- sort(unique(mean_df$m[mean_df$dgp == dgp]))
      for (i in seq_along(m_vals)) {
        dgp_lab <- if (i == 1L) dgp else ""
        line <- make_mean_row(mean_df, dgp, m_vals[i], gamma, t_subset, dgp_lab)
        suffix <- if (i == length(m_vals) && dgp != tail(dgp_vals, 1)) " \\\\ \\hline" else " \\\\"
        lines <- c(lines, paste0(line, suffix))
      }
    }
  }
  c(lines, "\\bottomrule", "\\end{tabular}")
}

write_mean_null_table <- function(t_subset, tag, label_name) {
  title <- paste0("Horizon block $T\\in\\{", paste(t_subset, collapse = ","), "\\}$")
  lines <- c(
    "\\begin{table}[htbp]",
    paste0("\\caption{{\\color{blue}Simulated type I error rates (\\%) for RSMS (R), SSMS (S), and CSMS (C) under mean monitoring with training sizes $m\\in\\{100,250,500\\}$ and horizon block $T\\in\\{", paste(t_subset, collapse = ","), "\\}$.}}"),
    paste0("\\label{", label_name, "}"),
    "\\centering",
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{3pt}",
    "\\renewcommand{\\arraystretch}{1.05}",
    "\\begin{adjustbox}{max width=\\textwidth,center}",
    build_mean_null_block(t_subset, title),
    "\\end{adjustbox}",
    "\\end{table}"
  )
  out_file <- file.path(RES_DIR, paste0("mean_null_table_", tag, "_blue.tex"))
  writeLines(lines, out_file)
  out_file
}

par_lines <- c(
  NULL
)

par_col_spec <- function() {
  paste0(
    "@{}",
    ">{\\centering\\arraybackslash}m{0.09\\textwidth}",
    paste(rep(">{\\centering\\arraybackslash}m{0.065\\textwidth}", 12), collapse = ""),
    "@{}"
  )
}

build_par_null_block <- function(t_subset, title) {
  lines <- c(
    paste0("\\begin{tabular}{", par_col_spec(), "}"),
    "\\toprule",
    paste0("\\multicolumn{13}{c}{", title, "} \\\\"),
    "\\midrule",
    " & \\multicolumn{2}{c}{R} & \\multicolumn{2}{c}{S} & \\multicolumn{2}{c}{C} & \\multicolumn{2}{c}{R} & \\multicolumn{2}{c}{S} & \\multicolumn{2}{c}{C} \\\\",
    paste0(
      " & \\multicolumn{2}{c}{$T=", t_subset[1], "$} & \\multicolumn{2}{c}{$T=", t_subset[1], "$} & \\multicolumn{2}{c}{$T=", t_subset[1], "$} & ",
      "\\multicolumn{2}{c}{$T=", t_subset[2], "$} & \\multicolumn{2}{c}{$T=", t_subset[2], "$} & \\multicolumn{2}{c}{$T=", t_subset[2], "$} \\\\"
    ),
    "\\cmidrule(lr){2-7} \\cmidrule(lr){8-13}",
    " & 5\\% & 10\\% & 5\\% & 10\\% & 5\\% & 10\\% & 5\\% & 10\\% & 5\\% & 10\\% & 5\\% & 10\\% \\\\",
    "\\midrule"
  )

  for (gamma in gamma_grid) {
    lines <- c(lines, sprintf("\\multicolumn{13}{c}{{\\color{blue}$\\gamma = %.2f$}} \\\\ \\midrule", gamma))
    m_vals <- sort(unique(par_df$m))
    for (m in m_vals) {
      vals <- c()
      for (T in t_subset) {
        for (method in c("RSMS", "SSMS", "CSMS")) {
          for (alpha in alpha_grid) {
            vals <- c(vals, par_value(par_df, m, gamma, T, method, alpha))
          }
        }
      }
      lines <- c(lines, paste0(paste(c(paste0("$m=", m, "$"), vals), collapse = " & "), " \\\\"))
    }
    if (gamma != tail(gamma_grid, 1)) lines <- c(lines, "\\midrule")
  }

  c(lines, "\\bottomrule", "\\end{tabular}")
}

write_par_null_table <- function(t_subset, tag, label_name) {
  title <- paste0("Horizon block $T\\in\\{", paste(t_subset, collapse = ","), "\\}$")
  lines <- c(
    "{\\color{blue}",
    "\\begin{table}[htbp]",
    paste0("\\caption{Type I error rates (\\%) for RSMS (R), SSMS (S), and CSMS (C) under the PAR(1) null with training sizes $m\\in\\{100,250,500\\}$ and horizon block $T\\in\\{", paste(t_subset, collapse = ","), "\\}$.}"),
    paste0("\\label{", label_name, "}"),
    "\\centering",
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{2pt}",
    "\\renewcommand{\\arraystretch}{1.05}",
    "\\begin{adjustbox}{max width=\\textwidth,center}",
    build_par_null_block(t_subset, title),
    "\\end{adjustbox}",
    "\\end{table}",
    "}"
  )
  out_file <- file.path(RES_DIR, paste0("par_null_table_", tag, "_blue.tex"))
  writeLines(lines, out_file)
  out_file
}
mean_out_12 <- write_mean_null_table(c(1, 2), "t12", "table-mean-null-t12-blue")
mean_out_510 <- write_mean_null_table(c(5, 10), "t510", "table-mean-null-t510-blue")
par_out_12 <- write_par_null_table(c(1, 2), "t12", "table-null-par1-t12-blue")
par_out_510 <- write_par_null_table(c(5, 10), "t510", "table-null-par1-t510-blue")

msg("Saved: ", mean_out_12)
msg("Saved: ", mean_out_510)
msg("Saved: ", par_out_12)
msg("Saved: ", par_out_510)
msg("DONE.")
