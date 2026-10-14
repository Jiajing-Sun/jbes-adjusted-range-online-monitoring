# Reporting only: this script reads completed theory-check output.
# It performs no data selection or new Monte Carlo simulation. Every scalar
# configuration appears in the tables. R/S/C means RSMS/SSMS/CSMS; displayed
# probabilities are percentages, while source CSV values are probabilities.
# CADD gaps are divided by sqrt(m), and their standard errors are rescaled
# by the same factor. The plot retains negative as well as positive shifts.

# Tables and figures derived from saved theory-check output; no new simulation.
root <- normalizePath(Sys.getenv('JBES_PROJECT_ROOT', getwd()))
out <- file.path(root, 'results', 'theory_assets')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, 'figures'), showWarnings = FALSE)
methods <- c('RSMS', 'SSMS', 'CSMS')
a <- read.csv(file.path(root, 'results', 'scalar_theory', 'scalar_summary.csv'))
p <- read.csv(file.path(root, 'results', 'scalar_theory', 'scalar_paired.csv'))
sel <- read.csv(file.path(root, 'results', 'scalar_theory', 'scalar_selection.csv'))
stopifnot(file.exists(file.path(root, 'results', 'scalar_theory', 'session.txt')))
configs <- unique(a[c('noise', 'pi', 'm')])
short_noise <- function(s) if (s == 'iid_gaussian') 'Gaussian' else r'(AR(1)-$t_5$)'
percent <- function(x) paste(sprintf('%.1f', 100 * x), collapse = '/')
lookup <- function(z, scenario, calibration, ms = methods) {
  x <- a[a$noise == z$noise & a$pi == z$pi & a$m == z$m & a$scenario == scenario & a$calibration == calibration, ]
  x[match(ms, x$method), ]
}
table_tex <- function(filename, caption, label, columns, header, rows, notes) {
  s <- c(
    r'(\begin{table}[!htbp])', paste0(r'(\caption{)', caption, '}'),
    paste0(r'(\label{)', label, '}'), r'(\centering\small)',
    r'(\setlength{\tabcolsep}{4pt})', r'(\renewcommand{\arraystretch}{1.1})',
    paste0(r'(\begin{tabular}{)', columns, '}'), r'(\toprule)', header, r'(\midrule)', rows,
    r'(\bottomrule\end{tabular})', r'(\par\medskip\begin{minipage}{\linewidth}\footnotesize)',
    paste0(r'(\emph{Notes:} )', notes), r'(\end{minipage}\end{table})'
  )
  writeLines(s, file.path(out, filename))
}
size_rows <- power_rows <- delay_rows <- character()
for (i in seq_len(nrow(configs))) {
  z <- configs[i, ];prefix <- paste(short_noise(z$noise), sprintf('%.2f', z$pi), z$m, sep = ' & ')
  x <- sel[sel$noise == z$noise & sel$pi == z$pi & sel$m == z$m & sel$phase == 2, ]
  size_rows <- c(size_rows, paste(prefix, percent(lookup(z, 'null', 'clean')$alarm),
    percent(lookup(z, 'null', 'oracle')$alarm), percent(lookup(z, 'null', 'fitted', methods[1:2])$alarm),
    sprintf('%.2f', 100 * x$selection),
    sep = ' & '
  ) |> paste0(r'( \\)'))
  gap <- p[p$noise == z$noise & p$pi == z$pi & p$m == z$m & p$scenario == 'local_h2' & p$calibration == 'oracle', ]
  power_rows <- c(power_rows, paste(prefix, percent(lookup(z, 'local_h2', 'oracle')$post),
    sprintf('%.2f (%.2f)', 100 * gap$post_gap_R_minus_S, 100 * gap$paired_se),
    percent(lookup(z, 'local_h2', 'fitted', methods[1:2])$post),
    sep = ' & '
  ) |> paste0(r'( \\)'))
  d <- lookup(z, 'fixed_0.75', 'oracle', methods[1:2]);gap <- p[p$noise == z$noise & p$pi == z$pi & p$m == z$m & p$scenario == 'fixed_0.75' & p$calibration == 'oracle', ]
  delay_rows <- c(delay_rows, paste(prefix, paste(d$post_count, collapse = '/'),
    paste(sprintf('%.1f', d$cadd), collapse = '/'), sprintf('%.3f (%.3f)', gap$scaled_gap, gap$cadd_gap_se / sqrt(z$m)),
    sprintf('%.3f', gap$theory_scaled_gap),
    sep = ' & '
  ) |> paste0(r'( \\)'))
}
table_tex(
  'scalar_size.tex', 'No-additional-change alarms under persistent scalar training contamination.', 'tab:scalar-new-size', 'llrcccc',
  r'(Noise & $\pi$ & $m$ & Clean R/S/C & Oracle R/S/C & Fitted R/S & Split (\%) \\)', size_rows,
  r'(Entries are percentages from 5,000 independent null paths. R/S/C denote RSMS/SSMS/CSMS. Each oracle threshold uses 10,000 separate contaminated-null paths. Clean thresholds are the operational clean-training references; their alarm rates under contamination are not clean-null size. Fitted thresholds and the common split decision use training observations only. At an alarm probability of 5\%, the marginal evaluation standard error is 0.31 percentage points. Calibration uncertainty is separate.)'
)
table_tex(
  'scalar_power.tex', 'Local post-change detection under persistent scalar training contamination, $h=2$.', 'tab:scalar-new-power', 'llrccc',
  r'(Noise & $\pi$ & $m$ & Oracle R/S/C (\%) & R--S (SE), pp & Fitted R/S (\%) \\)', power_rows,
  r'(Each row uses 5,000 alternative paths, independent of calibration and null evaluation. The R--S difference uses paired paths; its standard error conditions on the estimated oracle thresholds. Fitted-calibration power must be read together with its null alarm rate in Table~\ref{tab:scalar-new-size}. CSMS uses the same contaminated-null information for its oracle calibration. No oracle path in these scalar experiments crossed before the additional change.)'
)
table_tex(
  'scalar_delay.tex', 'Conditional delay under a fixed additional scalar change, $d=0.75b$.', 'tab:scalar-new-delay', 'llrcccc',
  r'(Noise & $\pi$ & $m$ & Detections R/S & CADD R/S & Scaled gap (SE) & Limit \\)', delay_rows,
  r'(The denominator is 5,000 independent paths per row; oracle calibration uses 10,000 separate paths. CADD conditions on surviving to the additional change and detecting by the common horizon. There were no pre-change alarms. The scaled gap is $(\mathrm{CADD}_S-\mathrm{CADD}_R)/\sqrt m$; its paired ratio-estimator standard error conditions on the estimated thresholds. The final column is the limit in Proposition~\ref{prop:oracle-delay-comparison}; it is not a finite-sample equality.)'
)
pdf(file.path(root, 'figures', 'scalar_oracle_local_power.pdf'), width = 7.4, height = 7.1, useDingbats = FALSE)
layout(matrix(c(1, 2, 3, 4, 5, 5), 3, 2, byrow = TRUE), heights = c(1, 1, .18))
par(mar = c(4.2, 3.9, 2.5, .6), oma = c(0, 0, 0, 0), cex = .87)
colors <- c('#0072B2', '#D55E00', '#009E73');symbols <- c(16, 17, 15)
for (noise in c('iid_gaussian', 'ar_t5')) for (pi0 in c(.25, .5)) {
  x <- a[a$noise == noise & a$pi == pi0 & a$m == 4000 & a$calibration == 'oracle' & grepl('local_h', a$scenario), ]
  x$h <- as.numeric(sub('local_h', '', x$scenario));plot(NA, xlim = c(-4, 4), ylim = c(0, .8), xlab = 'Local change h', ylab = 'Post-change detection probability', main = paste(if (noise == 'iid_gaussian') 'Gaussian' else 'AR(1)-t5', ', pi =', pi0))
  abline(h = .05, col = 'grey80', lty = 3)
  for (j in 1:3) {
    v <- x[x$method == methods[j], ];v <- v[order(v$h), ];lines(v$h, v$post, col = colors[j], pch = symbols[j], type = 'b', lwd = 1.4)
    if (j <= 2) {
      hh <- seq(-4, 4, length.out = 201);vv <- if (j == 1) 2 + 4 / pi0 else 2 + 4 * (3 * pi0^2 + 2 * pi0 + 1) / (5 * pi0^2)
      lines(hh, pnorm(qnorm(.95) - hh * 1.6 / sqrt(vv), lower.tail = FALSE), col = colors[j], lty = 2, lwd = 1.3)
    }
  }
}
par(mar = c(0, 0, 0, 0));plot.new();legend('center', legend = c(methods, 'R/S theoretical limits'), col = c(colors, 'black'), pch = c(symbols, NA), lty = c(1, 1, 1, 2), horiz = TRUE, bty = 'n', cex = .9)
dev.off()
cat('Scalar tables and figure written.\n')

# The bivariate report is generated only after every configured design finishes.
# Delta=0 tables compare increasing m; a separate table retains the transverse
# shifts at the larger m. The complete CSV retains all eta/delta combinations.
bivariate_dir <- file.path(root, 'results', 'bivariate_theory')
if (file.exists(file.path(bivariate_dir, 'session.txt'))) {
  b <- read.csv(file.path(bivariate_dir, 'bivariate_summary.csv'))
  bp <- read.csv(file.path(bivariate_dir, 'bivariate_paired.csv'))
  designs <- unique(b[c('theta', 'rho', 'm')])
  central_rows <- transverse_rows <- character()
  for (i in seq_len(nrow(designs))) {
    design <- designs[i, ]
    z <- b[b$theta == design$theta & b$rho == design$rho & b$m == design$m, ]
    null <- z[z$phase == 2, ];null <- null[match(methods, null$method), ]
    local <- z[z$phase == 3 & z$eta == 2 & z$delta == 0, ]
    local <- local[match(methods, local$method), ]
    pair <- bp[bp$theta == design$theta & bp$rho == design$rho & bp$m == design$m &
      bp$phase == 3 & bp$eta == 2 & bp$delta == 0 & bp$first == 'RSMS' & bp$second == 'SSMS', ]
    central_rows <- c(central_rows, paste(sprintf('%.2f', design$theta),
      sprintf('%.1f', design$rho), design$m, percent(null$post), percent(local$post),
      sprintf('%.2f (%.2f)', 100 * pair$post_difference, 100 * pair$paired_se),
      percent(local$theory_power[1:2]),
      sep = ' & '
    ) |> paste0(r'( \\)'))
    if (design$m == max(designs$m)) for (delta in c(-2, 0, 2)) {
      local <- z[z$phase == 3 & z$eta == 2 & z$delta == delta, ]
      local <- local[match(methods, local$method), ]
      transverse_rows <- c(transverse_rows, paste(sprintf('%.2f', design$theta),
        sprintf('%.1f', design$rho), delta, percent(local$post),
        percent(local$theory_power[1:2]),
        sep = ' & '
      ) |> paste0(r'( \\)'))
    }
  }
  table_tex(
    'bivariate_central.tex', 'Full-rank bivariate calibration and local detection, $\\eta=2$, $\\delta=0$.',
    'tab:bivariate-new-central', 'rrrcccc',
    r'($\vartheta$ & $\rho$ & $m$ & Null R/S/C (\%) & Detection R/S/C (\%) & R--S (SE), pp & Limit R/S (\%) \\)',
    central_rows, r'(Each oracle boundary uses 10,000 contaminated-null paths. Null and alternative evaluation use separate sets of 5,000 paths. R/S/C denote RSMS/SSMS/CSMS. Standard errors use paired alternative paths and condition on the estimated thresholds. The two theoretical limits are from Proposition~\ref{prop:bivariate-oracle-comparison}; no CSMS limit is asserted here.)'
  )
  table_tex(
    'bivariate_transverse.tex', 'Transverse components of the full-rank bivariate local shift, $m=2000$, $\\eta=2$.',
    'tab:bivariate-new-transverse', 'rrrcc',
    r'($\vartheta$ & $\rho$ & $\delta$ & Detection R/S/C (\%) & Limit R/S (\%) \\)',
    transverse_rows, r'(The additional mean change is $m^{-1/2}(-\eta,\delta)^\top$. The RSMS and SSMS first-order limits are independent of $\delta$ in this specified profile. Finite-sample differences are retained. There are 5,000 alternative paths per setting, shared across the transverse components and methods, with independent contaminated-null calibration.)'
  )
  pdf(file.path(root, 'figures', 'bivariate_oracle_local_power.pdf'), width = 7.4, height = 7.1, useDingbats = FALSE)
  layout(matrix(c(1, 2, 3, 4, 5, 5), 3, 2, byrow = TRUE), heights = c(1, 1, .18))
  par(mar = c(4.2, 3.9, 2.5, .6), oma = c(0, 0, 0, 0), cex = .87)
  for (theta in c(.25, .75)) for (rho in c(0, .6)) {
    z <- b[b$theta == theta & b$rho == rho & b$m == 2000 & b$phase == 3 & b$delta == 0, ]
    plot(NA,
      xlim = c(-4, 4), ylim = c(0, .8), xlab = expression(eta), ylab = 'Post-change detection probability',
      main = bquote(vartheta == .(theta) ~ ',' ~ rho == .(rho))
    )
    abline(h = .05, col = 'grey80', lty = 3)
    for (j in 1:3) {
      points <- z[z$method == methods[j], ];points <- points[order(points$eta), ]
      lines(points$eta, points$post, col = colors[j], pch = symbols[j], type = 'b', lwd = 1.4)
      if (j <= 2) {
        eta <- seq(-4, 4, length.out = 201)
        variance <- if (j == 1) 14 else 6 + 64 * 3 * (23 * theta^2 + 59 * theta + 46) / (160 * (theta^2 + theta + 2)^2)
        lines(eta, pnorm(qnorm(.95) - eta * 1.6 / sqrt(variance), lower.tail = FALSE), col = colors[j], lty = 2, lwd = 1.3)
      }
    }
  }
  par(mar = c(0, 0, 0, 0));plot.new()
  legend('center',
    legend = c(methods, 'R/S theoretical limits'), col = c(colors, 'black'),
    pch = c(symbols, NA), lty = c(1, 1, 1, 2), horiz = TRUE, bty = 'n', cex = .9
  )
  dev.off()
  cat('Bivariate tables and figure written.\n')
}
