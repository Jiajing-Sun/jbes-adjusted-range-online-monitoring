# Necessary threshold check for uniform training-hull containment.
#
# The criterion c_S/c_R >= 2(q-1) is necessary, not sufficient. A ratio below
# this bound rules out the uniform containment event at these thresholds;
# a larger ratio does not prove that the event is attainable or likely.
# These calculations use recorded, unrounded reference values. They do not
# simulate training paths or replace the empirical geometry diagnostics.
root <- normalizePath(Sys.getenv('JBES_PROJECT_ROOT', getwd()))
results <- file.path(root, 'results')
finite <- read.csv(file.path(results, 'affine_finite_critical_values.csv'))
open <- read.csv(file.path(results, 'affine_open_end_critical_values.csv'))
open$T <- Inf
keys <- c('d', 'T', 'gamma', 'alpha')
columns <- c(keys, 'method', 'critical_value', 'n_rep')
reference <- rbind(finite[columns], open[columns])
r <- subset(reference, method == 'RSMS')
s <- subset(reference, method == 'SSMS')
joined <- merge(r, s, by = keys, suffixes = c('_R', '_S'), sort = TRUE)
stopifnot(
  nrow(joined) == 400L, !anyDuplicated(joined[keys]),
  all(joined$n_rep_R == 10000), all(joined$n_rep_S == 10000)
)
joined$threshold_ratio <- joined$critical_value_S / joined$critical_value_R
joined$necessary_lower_bound <- 2 * (joined$d - 1)
joined$uniform_containment_ruled_out <-
  joined$threshold_ratio < joined$necessary_lower_bound
write.csv(joined, file.path(results, 'containment_reference_diagnostics.csv'), row.names = FALSE)
summary <- do.call(rbind, lapply(sort(unique(joined$d)), function(q) {
  z <- subset(joined, d == q)
  data.frame(q,
    lower_bound = 2 * (q - 1), minimum = min(z$threshold_ratio),
    maximum = max(z$threshold_ratio), ruled_out = sum(z$uniform_containment_ruled_out),
    configurations = nrow(z)
  )
}))
write.csv(summary, file.path(results, 'containment_reference_summary.csv'), row.names = FALSE)
shown <- subset(summary, q %in% c(1, 2, 5, 8, 12, 20))
lines <- c(
  '\\begin{table}[!htbp]', '\\centering',
  '\\caption{Necessary containment-bound checks for the dimensions used in the numerical experiments.}',
  '\\label{tab:containment-reference-check}',
  '\\begin{tabular}{rrrrr}', '\\toprule',
  '$q$ & $2(q-1)$ & Minimum $c_S/c_R$ & Maximum $c_S/c_R$ & Ruled out / checked \\\\',
  '\\midrule'
)
for (i in seq_len(nrow(shown))) {
  z <- shown[i, ]
  lines <- c(lines, sprintf(
    '%d & %d & %.2f & %.2f & %d / %d \\\\',
    z$q, z$lower_bound, z$minimum, z$maximum, z$ruled_out, z$configurations
  ))
}
lines <- c(
  lines, '\\bottomrule', '\\end{tabular}',
  '\\par\\medskip\\begin{minipage}{.95\\textwidth}\\small',
  '\\emph{Notes:} Each dimension has 20 threshold combinations: $T\\in\\{1,2,5,10,\\infty\\}$, $\\gamma\\in\\{0,0.15\\}$ and $\\alpha\\in\\{0.05,0.10\\}$. Calculations use unrounded reference values. A failed necessary bound rules out uniform containment; passing it is not an attainability or probability statement. The replication output retains all 400 checks for $q=1,\\ldots,20$.',
  '\\end{minipage}', '\\end{table}'
)
writeLines(lines, file.path(results, 'containment_reference_table.tex'))
cat('Saved all 400 necessary-bound checks and the manuscript summary table.\n')
