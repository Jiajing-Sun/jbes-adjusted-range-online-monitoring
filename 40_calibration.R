# ============================================================
# 40_calibration.R  -- Monte Carlo calibration for size-adjusted power
# ============================================================

# Calibrate empirical critical values c_hat for each method in METHODS
# given a function gen_psi(rep_id) that returns psi_mat of size (m+mT) x q.
#
# Returns a data.frame with columns:
#   method, c_hat, alpha, m, T, gamma, tag
calibrate_critical_values <- function(METHODS, gen_psi,
                                      m, T, gamma,
                                      alpha = 0.05,
                                      N_CAL = 2000,
                                      n_cores = 1L,
                                      seed = 1L,
                                      tag = "") {
  stopifnot(alpha > 0 && alpha < 1)

  # one replicate -> named vector of max stats by method
  one_rep <- function(i) {
    psi <- gen_psi(i)
    out <- sapply(METHODS, function(mm) {
      if (mm == "RSMS" && ncol(psi) >= 3L && exists("rsms_maximum",mode="function")) {
        rsms_maximum(psi,m,T,gamma)
      } else max(compute_path(mm, psi, m, T, gamma))
    })
    out
  }

  res_list <- par_lapply_safe(1:N_CAL, one_rep, n_cores = n_cores, seed = seed)
  maxmat <- do.call(rbind, res_list) # N_CAL x length(METHODS)
  colnames(maxmat) <- METHODS

  qhat <- apply(maxmat, 2, function(z) as.numeric(quantile(z, probs = 1 - alpha, type = 8)))

  data.frame(method = METHODS,
             c_hat = as.numeric(qhat),
             alpha = alpha,
             m = m,
             T = T,
             gamma = gamma,
             N_CAL = N_CAL,
             tag = tag,
             stringsAsFactors = FALSE)
}
