# ============================================================
# 41_csms_bandwidth_assets.R
# - Build LaTeX tables for focused CSMS bandwidth-sensitivity appendix module
# ============================================================

rm(list = ls())

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

RES_DIR <- file.path(ROOT, "results")
df <- read.csv(file.path(RES_DIR, "csms_bandwidth_sensitivity.csv"), stringsAsFactors = FALSE)

fmt_pct <- function(x) sprintf("%.1f", 100 * x)
fmt_add <- function(x) sprintf("%.1f", x)

make_table <- function(T0, outfile, label) {
  sub <- df[df$T == T0, ]
  sub <- sub[order(ifelse(sub$method == "RSMS", 999, sub$eta)), ]

  row_to_tex <- function(row) {
    if (row[["method"]] == "RSMS") {
      c1 <- "\\textbf{RSMS benchmark}"
      c2 <- "--"
      c3 <- "--"
    } else {
      c1 <- sprintf("$\\eta = %.2f$", row[["eta"]])
      c2 <- sprintf("%d", as.integer(row[["L"]]))
      c3 <- "CSMS"
    }
    paste0(
      c1, " & ", c2, " & ", c3, " & ",
      fmt_pct(row[["size_05"]]), " & ",
      fmt_pct(row[["size_10"]]), " & ",
      fmt_pct(row[["early_DP_post_adj"]]), " / ", fmt_add(row[["early_CADD_adj"]]), " & ",
      fmt_pct(row[["late_DP_post_adj"]]), " / ", fmt_add(row[["late_CADD_adj"]]), " \\\\"
    )
  }

  lines <- vapply(seq_len(nrow(sub)), function(i) row_to_tex(sub[i, ]), character(1))

  cap <- paste0(
    "{\\color{blue}CSMS bandwidth-sensitivity check for DGP1 with $m=250$, $\\gamma=0.15$, and abrupt changes at horizon $T=",
    T0,
    "$. The lag index is $L_{\\eta}=\\lceil \\eta \\cdot \\lceil 1.06 m^{1/5}\\rceil\\rceil$, with kernel bandwidth $h_{\\eta}=L_{\\eta}+1$. Raw size uses the asymptotic critical values already reported in the paper; the early-break alternative uses $\\Delta=0.50$, and the late-break alternative uses $\\Delta=1.50$. The last row reports the corresponding tuning-free RSMS benchmark.}"
  )

  tex <- c(
    "\\begin{table}[htbp]",
    paste0("\\caption{", cap, "}"),
    paste0("\\label{", label, "}"),
    "\\centering",
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{5pt}",
    "\\renewcommand{\\arraystretch}{1.05}",
    "\\begin{tabular}{lcccccc}",
    "\\toprule",
    "$\\eta$ & $L_{\\eta}$ & Method & Size 5\\% & Size 10\\% & Early $\\mathrm{DP}_{\\mathrm{post}}$/CADD & Late $\\mathrm{DP}_{\\mathrm{post}}$/CADD \\\\",
    "\\midrule",
    lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\par\\medskip",
    "\\begin{minipage}{0.94\\linewidth}",
    "\\footnotesize \\emph{Notes:} Size entries are rejection frequencies under the asymptotic CSMS critical values. The early-break column reports the size-adjusted post-break detection probability (in \\%) and conditional average detection delay (CADD) for $k^*=0.1mT$ with $\\Delta=0.50$; the late-break column reports the same quantities for $k^*=0.8mT$ with $\\Delta=1.50$.",
    "\\end{minipage}",
    "\\end{table}"
  )

  writeLines(tex, file.path(RES_DIR, outfile))
}

make_table(1, "csms_bw_table_t1_blue.tex", "table-csms-bw-t1-blue")
make_table(5, "csms_bw_table_t5_blue.tex", "table-csms-bw-t5-blue")
