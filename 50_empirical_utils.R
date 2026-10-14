# ============================================================
# 50_empirical_utils.R
# - Helpers for the USD/GBP empirical application
# ============================================================

load_gbpusd_data <- function(data_path) {
  ensure_pkgs(c("readr", "dplyr", "lubridate"))
  raw <- readr::read_csv(data_path, show_col_types = FALSE, name_repair = "unique_quiet")
  if (ncol(raw) > 0L && grepl("^\\.\\.\\.[0-9]+$", names(raw)[1])) {
    raw <- raw[, -1, drop = FALSE]
  }
  raw$Date <- lubridate::ymd(as.character(raw$Date))
  raw$Year <- lubridate::year(raw$Date)
  raw$Close_inv <- 100 / raw$Close
  raw
}

prepare_gbpusd_year <- function(raw_df, year) {
  ensure_pkgs("dplyr")
  dplyr::filter(raw_df, Year == year) |>
    dplyr::select(Date, DateTime, Close_inv)
}

daily_mean_series <- function(year_df) {
  ensure_pkgs("dplyr")
  dplyr::group_by(year_df, Date) |>
    dplyr::summarise(daily_mean = mean(Close_inv, na.rm = TRUE), .groups = "drop")
}

build_daily_fd <- function(year_df, nbasis_cap = 21L, min_points = 10L) {
  ensure_pkgs(c("dplyr", "fda"))

  daily_list <- dplyr::group_by(year_df, Date) |>
    dplyr::group_split()

  n_by_day <- sapply(daily_list, nrow)
  min_obs <- min(n_by_day)
  max_obs <- max(n_by_day)
  nbasis <- as.integer(nbasis_cap) # prespecified; no monitoring-period dimension choice
  stopifnot(nbasis >= 5L)

  basis <- fda::create.bspline.basis(rangeval = c(0, 1), norder = 4, nbasis = nbasis)
  coefs <- matrix(NA_real_, nrow = nbasis, ncol = length(daily_list))
  day_names <- character(0)
  valid_days <- 0L

  for (i in seq_along(daily_list)) {
    day_data <- daily_list[[i]]
    if (nrow(day_data) < min_points) next

    valid_days <- valid_days + 1L
    day_names[valid_days] <- as.character(unique(day_data$Date))

    t_index <- seq(0, 1, length.out = nrow(day_data))
    fd_obj <- fda::smooth.basis(argvals = t_index, y = day_data$Close_inv, fdParobj = basis)
    coefs[, valid_days] <- fd_obj$fd$coefs
  }

  coefs <- coefs[, seq_len(valid_days), drop = FALSE]
  dates <- as.Date(day_names)

  fd_obj <- fda::fd(
    coef = coefs,
    basisobj = basis,
    fdnames = list(
      time = paste0("t", seq(0, 100, length.out = 101)),
      replicates = as.character(dates),
      values = "GBPUSD_Close_inv"
    )
  )

  list(
    fd = fd_obj,
    dates = dates,
    nbasis = nbasis,
    min_obs = min_obs,
    max_obs = max_obs,
    valid_days = valid_days
  )
}

evaluate_fd_matrix <- function(fd_obj, n_grid = 301L) {
  grid <- seq(0, 1, length.out = n_grid)
  values <- t(fda::eval.fd(grid, fd_obj$fd))
  list(grid = grid, values = values)
}

choose_fpca_rank <- function(sdev, fve_threshold = 0.95, max_dim = 8L) {
  ve <- sdev^2
  keep <- which(cumsum(ve) / sum(ve) >= fve_threshold)[1]
  keep <- max(1L, keep)
  min(keep, max_dim)
}

compute_empirical_score_matrix <- function(value_matrix, dates, m, T, fve_threshold = 0.95, max_dim = 8L) {
  total_len <- as.integer(m + m * T)
  if (total_len > nrow(value_matrix)) {
    stop("Not enough daily curves for m=", m, ", T=", T, ".")
  }

  curves <- value_matrix[seq_len(total_len), , drop = FALSE]
  dates_use <- dates[seq_len(total_len)]

  train_curves <- curves[seq_len(m), , drop = FALSE]
  monitor_curves <- curves[(m + 1L):total_len, , drop = FALSE]

  pca_fit <- stats::prcomp(train_curves, center = TRUE, scale. = FALSE)
  k <- choose_fpca_rank(pca_fit$sdev, fve_threshold = fve_threshold, max_dim = max_dim)

  rotation <- pca_fit$rotation[, seq_len(k), drop = FALSE]
  center_vec <- pca_fit$center
  scores_train <- pca_fit$x[, seq_len(k), drop = FALSE]
  centered_monitor <- sweep(monitor_curves, 2, center_vec, "-")
  scores_monitor <- centered_monitor %*% rotation
  scores_all <- rbind(scores_train, scores_monitor)

  list(
    scores = scores_all,
    dates = dates_use,
    d_eff = ncol(scores_all),
    total_len = total_len
  )
}

feasible_empirical_configs <- function(dates, event_date, m_values, T_values, gamma_values) {
  event_date <- as.Date(event_date)
  event_idx <- match(event_date, dates)
  if (is.na(event_idx)) stop("Event date ", event_date, " not found in the daily curve dates.")

  grid <- expand.grid(
    m = m_values,
    T = T_values,
    gamma = gamma_values,
    stringsAsFactors = FALSE
  )
  grid$total_len <- grid$m + grid$m * grid$T
  subset(grid, total_len <= length(dates) & event_idx > m & event_idx <= total_len)
}

run_empirical_monitoring <- function(score_mat, dates, m, T, gamma, alpha = 0.05) {
  mT <- as.integer(round(m * T))
  d_eff <- ncol(score_mat)
  if (d_eff > 8L) {
    score_mat <- score_mat[, seq_len(8L), drop = FALSE]
    d_eff <- 8L
  }

  methods <- c("RSMS", "SSMS", "CSMS")
  rows <- lapply(methods, function(method) {
    Mvec <- compute_path(method, score_mat, m, T, gamma)
    cval <- switch(method,
      RSMS = get_critical_value_rsms(d = d_eff, T = T, gamma = gamma, alpha_level = alpha),
      SSMS = get_critical_value_ssms(d = d_eff, T = T, gamma = gamma, alpha_level = alpha),
      CSMS = get_critical_value_csms(d = d_eff, T = T, gamma = gamma, alpha_level = alpha)
    )
    stop_k <- first_exceed(Mvec, cval)
    reject <- stop_k <= mT
    first_date <- if (reject) as.character(dates[m + stop_k]) else NA_character_

    data.frame(
      method = method,
      d_eff = d_eff,
      statistic = max(Mvec),
      critical_value = cval,
      reject = reject,
      stop_k = stop_k,
      first_detection_date = first_date,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out$m <- m
  out$T <- T
  out$gamma <- gamma
  out
}
