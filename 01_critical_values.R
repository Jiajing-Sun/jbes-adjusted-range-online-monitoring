# ============================================================
# 01_critical_values.R  -- asymptotic critical values (paper tables)
# ============================================================

# Helpers for finite/open-end critical-value lookup
normalize_cv_T_key <- function(T) {
  if (length(T) != 1L) stop("T must be scalar")
  if (is.character(T)) {
    key <- tolower(trimws(T))
    if (key %in% c("inf", "infty", "infinity")) return("Inf")
    T_num <- suppressWarnings(as.numeric(T))
    if (!is.na(T_num)) T <- T_num
  }
  if (is.infinite(T)) return("Inf")
  if (T %in% c(1, 2, 5, 10)) return(as.character(as.integer(T)))
  stop("T must be 1, 2, 5, 10, or Inf")
}

format_cv_alpha_key <- function(alpha_level) {
  if (length(alpha_level) != 1L) stop("alpha_level must be scalar")
  sprintf("%.2f", as.numeric(alpha_level))
}

open_end_cv_path <- function() {
  root <- getOption("JBES_ROOT")
  if (is.null(root) || !nzchar(root)) root <- getwd()
  file.path(root, "results", "affine_open_end_critical_values.csv")
}

finite_cv_path <- function() {
  root <- getOption("JBES_ROOT")
  if (is.null(root) || !nzchar(root)) root <- getwd()
  file.path(root, "results", "affine_finite_critical_values.csv")
}

load_open_end_critical_values <- function() {
  path <- open_end_cv_path()
  if (!file.exists(path)) {
    stop(
      "Open-end critical values not found at ", path,
      ". Run 00_open_end_critical_values.R first."
    )
  }
  read.csv(path, stringsAsFactors = FALSE)
}

load_finite_critical_values <- function() {
  path <- finite_cv_path()
  if (!file.exists(path)) {
    stop(
      "Finite-horizon critical values not found at ", path,
      ". Run 00_finite_critical_values.R first."
    )
  }
  read.csv(path, stringsAsFactors = FALSE)
}

lookup_open_end_cv <- function(method, d, gamma, alpha_level) {
  if (!gamma %in% c(0, 0.15)) stop("gamma must be 0 or 0.15")
  if (!alpha_level %in% c(0.05, 0.10)) stop("alpha_level must be 0.05 or 0.10")
  tab <- load_open_end_critical_values()
  hit <- tab$method == toupper(method) &
    tab$d == d &
    abs(tab$gamma - gamma) < 1e-12 &
    abs(tab$alpha - alpha_level) < 1e-12
  if (!any(hit)) {
    stop(
      "Missing open-end critical value for method=", method,
      ", d=", d, ", gamma=", gamma, ", alpha=", alpha_level, "."
    )
  }
  as.numeric(tab$critical_value[which(hit)[1]])
}

lookup_finite_cv <- function(method, d, T, gamma, alpha_level) {
  if (!gamma %in% c(0, 0.15)) stop("gamma must be 0 or 0.15")
  if (!alpha_level %in% c(0.05, 0.10)) stop("alpha_level must be 0.05 or 0.10")
  if (!T %in% c(1, 2, 5, 10)) stop("T must be 1, 2, 5, or 10")
  tab <- load_finite_critical_values()
  hit <- tab$method == toupper(method) &
    tab$d == d &
    abs(tab$T - T) < 1e-12 &
    abs(tab$gamma - gamma) < 1e-12 &
    abs(tab$alpha - alpha_level) < 1e-12
  if (!any(hit)) {
    stop(
      "Missing finite critical value for method=", method,
      ", d=", d, ", T=", T, ", gamma=", gamma, ", alpha=", alpha_level, "."
    )
  }
  as.numeric(tab$critical_value[which(hit)[1]])
}

# RSMS (range-based SN) critical values
get_critical_value_rsms <- function(d, T, gamma, alpha_level) {
  # alpha_level should be 0.05 or 0.10
  T_key <- normalize_cv_T_key(T)
  if (!gamma %in% c(0, 0.15)) stop("gamma must be 0 or 0.15")
  if (!alpha_level %in% c(0.05, 0.10)) stop("alpha_level must be 0.05 or 0.10")
  if (d < 1) stop("d must be a positive integer")
  if (T_key == "Inf") return(lookup_open_end_cv("RSMS", d, gamma, alpha_level))
  lookup_finite_cv("RSMS", d, as.numeric(T_key), gamma, alpha_level)
}

# CSMS (classical monitoring) critical values
get_critical_value_csms <- function(d, T, gamma, alpha_level) {
  T_key <- normalize_cv_T_key(T)
  if (!gamma %in% c(0, 0.15)) stop("gamma must be 0 or 0.15")
  if (!alpha_level %in% c(0.05, 0.10)) stop("alpha_level must be 0.05 or 0.10")
  if (d < 1) stop("d must be a positive integer")
  if (T_key == "Inf") return(lookup_open_end_cv("CSMS", d, gamma, alpha_level))
  lookup_finite_cv("CSMS", d, as.numeric(T_key), gamma, alpha_level)
}

# SSMS (Shao-type SN monitoring) critical values
get_critical_value_ssms <- function(d, T, gamma, alpha_level) {
  T_key <- normalize_cv_T_key(T)
  if (!gamma %in% c(0, 0.15)) stop("gamma must be 0 or 0.15")
  if (!alpha_level %in% c(0.05, 0.10)) stop("alpha_level must be 0.05 or 0.10")
  if (d < 1) stop("d must be a positive integer")
  if (T_key == "Inf") return(lookup_open_end_cv("SSMS", d, gamma, alpha_level))
  lookup_finite_cv("SSMS", d, as.numeric(T_key), gamma, alpha_level)
}
