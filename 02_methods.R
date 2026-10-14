# ============================================================
# 02_methods.R  -- RSMS / SSMS / CSMS statistics for simulations
# ============================================================

# Requires: 00_utils.R, 01_critical_values.R

# Affine-equivariant adjusted-range normalization.  No covariance ridge,
# coordinatewise aggregation, or partial prewhitening is used by RSMS.
# Singular paths are rejected explicitly rather than silently regularized.
# Reusable gauge LP with numerical primal/dual certification. A failed warm
# start is re-solved on the SAME path and objective, never discarded/redrawn.
ae_lp_engine <- function(G) {
  G <- as.matrix(G);n <- nrow(G);q <- ncol(G);A <- cbind(G, -1)
  model <- highs::highs_model(
    L = rep(0, q + 1L), lower = rep(-Inf, q + 1L),
    upper = rep(Inf, q + 1L), A = A, lhs = rep(0, n), rhs = rep(1, n), maximum = TRUE
  )
  options <- list(
    threads = 1L, presolve = 'off', solver = 'simplex',
    simplex_strategy = 4L, output_flag = FALSE
  )
  solver <- highs::hi_new_solver(model)
  highs::hi_solver_set_options(solver, options)
  # With output_flag=FALSE, bypass only the package's per-call capture.output
  # sink on the verified package version; the same native solver is called.
  run_solver <- if (as.character(utils::packageVersion('highs')) == '1.14.0.2')
    getFromNamespace('solver_run', 'highs') else highs::hi_solver_run
  function(z) {
    if (max(abs(z)) == 0) return(0)
    attempts <- list()
    for (attempt in 1:3) {
      if (attempt == 2L) {
        highs::hi_solver_clear_solver(solver)
        highs::hi_solver_set_options(solver, list(presolve = 'on', simplex_strategy = 1L))
      }
      if (attempt == 3L) {
        solver <<- highs::hi_new_solver(model)
        highs::hi_solver_set_options(solver, c(
          options[c('threads', 'solver', 'output_flag')],
          list(
            presolve = 'on', simplex_strategy = 1L,
            primal_feasibility_tolerance = 1e-9, dual_feasibility_tolerance = 1e-9
          )
        ))
      }
      highs::hi_solver_set_objective(solver, seq_len(q) - 1L, z)
      run <- run_solver(solver);status <- highs::hi_solver_status(solver)
      sol <- highs::hi_solver_get_solution(solver)
      primal <- sum(sol$col_value[seq_len(q)] * z)
      dual <- sum(pmax(sol$row_dual, 0))
      row_error <- max(c(0, -sol$row_value, sol$row_value - 1))
      dual_error <- max(abs(sol$col_dual))
      gap <- abs(primal - dual)
      certificate <- isTRUE(sol$value_valid) && isTRUE(sol$dual_valid) &&
        all(is.finite(c(primal, dual, row_error, dual_error))) &&
        row_error <= 1e-7 && dual_error <= 1e-7 * (1 + max(abs(z))) &&
        gap <= 1e-7 * (1 + abs(primal) + abs(dual))
      attempts[[attempt]] <- list(
        run = run, status = status, primal = primal, dual = dual,
        row_error = row_error, dual_error = dual_error, gap = gap, certificate = certificate
      )
      if (status == 7L && certificate) {
        if (attempt > 1L || run != 0L) {
          directory <- Sys.getenv('JBES_LP_DIAGNOSTICS', file.path('results', 'lp_diagnostics'))
          dir.create(directory, recursive = TRUE, showWarnings = FALSE)
          saveRDS(list(G = G, z = z, attempts = attempts), tempfile('recovered_', directory, '.rds'))
          highs::hi_solver_set_options(solver, options)
        }
        return(max(0, primal))
      }
    }
    directory <- Sys.getenv('JBES_LP_DIAGNOSTICS', file.path('results', 'lp_diagnostics'))
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    file <- tempfile('failure_', directory, '.rds')
    saveRDS(list(G = G, z = z, attempts = attempts), file)
    stop('Gauge LP failed numerical certification; input retained in ', file)
  }
}

normalize_hac_kernel <- function(kernel = "bartlett") {
  k <- tolower(trimws(kernel))
  if (k %in% c("bartlett", "newey-west")) return("bartlett")
  if (k %in% c("parzen")) return("parzen")
  if (k %in% c("qs", "quadratic spectral", "quadraticspectral")) return("qs")
  stop("Unknown HAC kernel: ", kernel)
}

hac_kernel_weight <- function(x, kernel = "bartlett") {
  k <- normalize_hac_kernel(kernel)
  ax <- abs(x)

  if (k == "bartlett") {
    return(ifelse(ax <= 1, 1 - ax, 0))
  }

  if (k == "parzen") {
    return(ifelse(
      ax <= 0.5,
      1 - 6 * ax^2 + 6 * ax^3,
      ifelse(ax <= 1, 2 * (1 - ax)^3, 0)
    ))
  }

  # Quadratic spectral kernel
  out <- numeric(length(ax))
  out[ax < 1e-12] <- 1
  nz <- ax >= 1e-12
  z <- ax[nz]
  out[nz] <- (25 / (12 * pi^2 * z^2)) *
    (sin(6 * pi * z / 5) / (6 * pi * z / 5) - cos(6 * pi * z / 5))
  out
}

# ------------------------------------------------------------
# CSMS path: HAC LRV computed from training window
# Bandwidth default: "Silverman-type" rule: L = ceil(1.06 * m^(1/5))
# ------------------------------------------------------------
hac_lrv_bartlett <- function(psi_train, L = NULL) {
  hac_lrv(psi_train, L = L, kernel = "bartlett")
}

csms_path <- function(psi_mat, m, T, gamma, L = NULL, kernel = "bartlett", ridge = 1e-8) {
  q <- ncol(psi_mat)
  mT <- as.integer(round(m * T))
  stopifnot(nrow(psi_mat) >= m + mT)

  psi_train <- psi_mat[1:m, , drop = FALSE]
  psi_mon <- psi_mat[(m + 1):(m + mT), , drop = FALSE]

  Gamma <- hac_lrv(psi_train, L = L, kernel = kernel)
  Gamma_inv <- inv_pd(Gamma, ridge = ridge)

  S_mon <- apply(psi_mon, 2, cumsum)
  quad <- rowSums((S_mon %*% Gamma_inv) * S_mon)

  k <- 1:mT
  scaling <- (m^(-1)) * (1 + k / m)^(-2) * ((k / m) / (1 + k / m))^(-2 * gamma)
  as.numeric(quad * scaling)
}

# Convenience wrapper: compute M path by method name
compute_path <- function(method, psi_mat, m, T, gamma, ...) {
  method <- toupper(method)
  if (method == "RSMS") return(rsms_path(psi_mat, m, T, gamma, ...))
  if (method == "SSMS") return(ssms_path(psi_mat, m, T, gamma, ...))
  if (method == "CSMS") return(csms_path(psi_mat, m, T, gamma, ...))
  stop("Unknown method: ", method)
}

# Convenience wrapper: stop time given a critical value
compute_stop <- function(method, psi_mat, m, T, gamma, cval, ...) {
  Mvec <- compute_path(method, psi_mat, m, T, gamma, ...)
  stop_k <- first_exceed(Mvec, cval)
  list(stop_k = stop_k, max_stat = max(Mvec), Mvec = Mvec)
}

# Central package adapter; online wrappers above retain stopping definitions.
source(file.path(Sys.getenv("JBES_PROJECT_ROOT", unset = getwd()), "07_aersn_backend.R"))
