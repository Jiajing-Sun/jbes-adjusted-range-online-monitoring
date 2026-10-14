# ============================================================
# 31_empirical_plots.R
# - Rebuild daily USD/GBP plots for 2016 Brexit and 2020 COVID
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

ensure_pkgs(c("ggplot2", "dplyr", "readr", "lubridate"), auto_install = AUTO_INSTALL)

data_path <- env_string(
  "JBES_GBPUSD_PATH",
  file.path(ROOT, "data", "GBPUSD.csv")
)
code_fig_dir <- file.path(ROOT, "figures")
paper_fig_dir <- file.path(dirname(ROOT), "figures")
if (!dir.exists(code_fig_dir)) dir.create(code_fig_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(paper_fig_dir)) dir.create(paper_fig_dir, recursive = TRUE, showWarnings = FALSE)

raw_df <- load_gbpusd_data(data_path)

make_event_plot <- function(year, event_date, event_label, out_name) {
  year_df <- prepare_gbpusd_year(raw_df, year)
  daily_df <- daily_mean_series(year_df)

  p <- ggplot2::ggplot(daily_df, ggplot2::aes(x = Date, y = daily_mean)) +
    ggplot2::geom_line(color = "black", linewidth = 0.55) +
    ggplot2::geom_vline(xintercept = as.Date(event_date), linetype = "dashed", color = "grey30", linewidth = 0.5) +
    ggplot2::labs(
      title = paste0("USD/GBP in ", year, ": ", event_label),
      x = NULL,
      y = "Daily mean of 100 / Close"
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 13)
    )

  ggplot2::ggsave(file.path(code_fig_dir, out_name), p, width = 8, height = 4.2, dpi = 300)
  ggplot2::ggsave(file.path(paper_fig_dir, out_name), p, width = 8, height = 4.2, dpi = 300)
}

make_event_plot(2016, "2016-06-23", "Brexit referendum", "GBPUSD_2016_brexit_plot.png")
make_event_plot(2020, "2020-03-17", "COVID-19 shock", "GBPUSD_2020_covid_plot.png")

msg("Saved empirical figures to: ", code_fig_dir)
msg("Saved empirical figures to: ", paper_fig_dir)
msg("DONE.")
