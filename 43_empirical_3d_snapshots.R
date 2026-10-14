# ============================================================
# 43_empirical_3d_snapshots.R
# - Static 3D snapshots of daily USD/GBP functional curves
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
source(file.path(ROOT, "50_empirical_utils.R"))
jbes_set_root(ROOT)
msg("Project root: ", ROOT)

ensure_pkgs(c("lubridate", "fda"), auto_install = AUTO_INSTALL)

data_path <- env_string(
  "JBES_GBPUSD_PATH",
  file.path(ROOT, "data", "GBPUSD.csv")
)

code_fig_dir <- file.path(ROOT, "figures")
paper_fig_dir <- file.path(dirname(ROOT), "figures")
if (!dir.exists(code_fig_dir)) dir.create(code_fig_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(paper_fig_dir)) dir.create(paper_fig_dir, recursive = TRUE, showWarnings = FALSE)

raw_df <- load_gbpusd_data(data_path)

render_snapshot <- function(year, event_date, event_label, out_name, days_before = 60L, days_after = 20L) {
  year_df <- prepare_gbpusd_year(raw_df, year)
  fd_res <- build_daily_fd(year_df, nbasis_cap = 21L, min_points = 10L)
  eval_res <- evaluate_fd_matrix(fd_res, n_grid = 121L)

  event_idx <- match(as.Date(event_date), fd_res$dates)
  if (is.na(event_idx)) stop("Event date not found in the FD dates for year ", year)

  lo <- max(1L, event_idx - days_before)
  hi <- min(length(fd_res$dates), event_idx + days_after)
  vals <- eval_res$values[lo:hi, , drop = FALSE]
  dates_use <- fd_res$dates[lo:hi]

  zmat <- t(vals)
  x <- eval_res$grid
  y <- seq_len(nrow(vals))

  draw_plot <- function() {
    op <- par(mar = c(2.5, 2.8, 2.5, 0.5) + 0.1)
    on.exit(par(op))
    persp(
      x = x,
      y = y,
      z = zmat,
      theta = 45,
      phi = 28,
      expand = 0.75,
      shade = 0.25,
      col = "grey82",
      border = "grey40",
      ticktype = "detailed",
      xlab = "Normalized intraday time",
      ylab = "Retained date index",
      zlab = "GBP per 100 USD",
      main = sprintf("USD/GBP %s: %s", year, event_label)
    )
    mtext(
      sprintf("Window: %s to %s", format(min(dates_use), "%Y-%m-%d"), format(max(dates_use), "%Y-%m-%d")),
      side = 1,
      line = 1,
      cex = 0.8
    )
  }

  pdf(file.path(code_fig_dir, out_name), width = 7.0, height = 5.0, paper = "special")
  draw_plot()
  dev.off()
  file.copy(file.path(code_fig_dir, out_name), file.path(paper_fig_dir, out_name), overwrite = TRUE)
}

render_snapshot(2016, "2016-06-23", "Brexit referendum", "GBPUSD_2016_brexit_3d_snapshot.pdf")
render_snapshot(2020, "2020-03-17", "COVID-19 shock", "GBPUSD_2020_covid_3d_snapshot.pdf")

msg("Saved empirical 3D snapshots to: ", code_fig_dir)
msg("Saved empirical 3D snapshots to: ", paper_fig_dir)
msg("DONE.")
