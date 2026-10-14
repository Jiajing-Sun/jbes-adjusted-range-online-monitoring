# Online Monitoring of Structural Change with Adjusted-Range Self-Normalization

R code accompanying the manuscript by Meiting Zhu, Yongmiao Hong, Jiajing Sun and Oliver Linton. This revision uses **aersn 0.2.3** for adjusted-range geometry, quadratic self-normalization and HAC normalization.

All R files are directly in this directory. There are no datasets, numerical results, figures, PDFs, archives or compiled binaries. Scripts create local `results/`, `figures/` and `logs/` directories; the Git allowlist excludes these outputs. No variance-profile experiment is included.

## Install and verify

Use R 4.1 or later. Run from the repository root:

```sh
Rscript install_dependencies.R
Rscript verify_backend.R
Rscript smoke_pipeline.R
```

The installer requests aersn 0.2.3 and HiGHS 1.14.0.2 plus public CRAN dependencies. Checks cover public package gauges, independent matrix formulas, affine transformations, exact maxima, first crossings, all-lag QS weights, nested PSOCK workers and explicit failure on singular paths. The smoke run uses synthetic observations in a temporary workspace. These checks do not replace full Monte Carlo evaluation.

## Run the study

List stages and master seeds without starting computation:

```sh
Rscript reproduce_reported.R
```

Run stages separately, or run the complete sequence:

```sh
Rscript reproduce_reported.R --stage finite_q1q2
Rscript reproduce_reported.R --stage open_q1q2
Rscript reproduce_reported.R --stage mean_null
Rscript reproduce_reported.R --stage mean_alternative
Rscript reproduce_reported.R --run
```

The full study is computationally intensive: it includes 10,000-draw references through dimension 20, calibration, evaluation, sensitivity analyses and empirical stages. Generate references before dependent experiments. Moderate-dimensional experiments require merged references. Empirical stages require the separate licensed input and fail if it is absent. Precomputed results are not embedded in the R code.

The driver fixes production replication and worker counts to preserve recorded streams. For a small development run, call an individual script directly using its `JBES_*` environment variables. Do not describe reduced-run output as a full reproduction. Reference checkpoints include implementation hashes; use a new checkpoint directory after source changes.

## Source map

| Files | Purpose |
| --- | --- |
| `00_utils.R`, `01_critical_values.R` | Utilities, parallel streams and reference lookup |
| `02_methods.R`, `07_aersn_backend.R` | Online statistics and package normalizers |
| `03_affine_extrema.R`, `05_affine_stopping.R` | Certified full-grid maxima and first crossings |
| `04_reference_replicates.R`, `06_training_geometry.R` | Checkpoints and training-containment diagnostics |
| `10_dgps_mean.R`, `20_dgps_par1.R`, `30_construct_psi.R`, `40_calibration.R` | Generators, fitted scores and independent calibration |
| `00_finite_critical_values.R`, `00_open_end_critical_values.R`, `01_matched_grid_critical_values.R`, `02_merge_reported_reference.R` | Online reference distributions |
| `01_mean_null_calibrate_and_size.R`, `02_mean_alternatives_power_delay.R` | Mean-monitoring simulations |
| `11_par_null_calibrate_and_size.R`, `12_par_alternatives_power_delay.R` | Poisson autoregression simulations |
| `21_moderate_dim_sparse.R`, `22_dependence_sensitivity.R`, `23_training_contamination.R`, `28_persistent_training.R` | Dimension, dependence and contamination designs |
| `24_csms_bandwidth_sensitivity.R`, `25_csms_kernel_sensitivity.R`, `29_csms_automatic_bandwidth.R` | Prespecified CSMS sensitivity comparisons |
| `27_matched_grid_size_comparison.R` | Reference-grid size assessment |
| `50_empirical_utils.R`, `31_empirical_plots.R`, `32_empirical_monitoring_table.R`, `34_empirical_bootstrap_check.R` | USD/GBP preprocessing, monitoring and bootstrap |
| `30_scalar_theory_checks.R`, `31_bivariate_theory_checks.R` | Independent same-contaminated-null calibration and checks of the new scalar and directional propositions |
| `08_mc_audit.R` | First-crossing counts and paired Monte Carlo uncertainty |
| Other numbered R files | Local tables, plots and empirical surface coordinates |

## Statistical conventions

- Training estimates are frozen before monitoring. Online maxima and stopping rules are implemented here. Package offline p-values never replace sequential critical values.
- Multivariate RSMS uses the full increment hull, without a componentwise replacement or a hull ridge.
- Scalar gauges are analytical. Batched two-dimensional gauges use vertices from public aersn functions. Above dimension two, a reusable HiGHS dual LP evaluates the aersn hull. Its first objective for each hull is compared with `aersn_gauge()`, and each accelerated solve checks primal/dual feasibility and objective gap. Set `JBES_AERSN_GAUGE=reference` when invoking a component script to use the public package LP throughout.
- Screening certifies intervals on the original grid; it does not thin the grid. This source-only repository uses portable R loops; no private compiled accelerator is required.
- Default Bartlett lag is `L = ceiling(1.06 * m^(1/5))`, with package bandwidth `h = L + 1`. SSMS and CSMS retain the documented `1e-8` ridge. QS uses all lags: regenerate its sensitivity results rather than reusing older truncated output.
- Independent calibration, null evaluation and alternative evaluation serve distinct purposes. Size-adjusted power is a model-specific benchmark, not an operational claim that the true data-generating process is known. Retain CSMS advantages and ranking reversals; do not choose tuning rules to lower competitor power.
- Numerical failures abort a run; paths are not silently discarded or replaced. An aborted run is not a completed Monte Carlo result.

The theory-check scripts use separate calibration, null-evaluation and alternative streams. Defaults are 10,000 calibration and 5,000 evaluation replications. Both signs of local changes and the bivariate transverse directions are retained. The scalar fixed-change delay experiment is separate from its local-power experiment. Paired evaluation errors condition on estimated thresholds; calibration-quantile brackets are reported separately. These scripts implement the specified models and do not establish a universal ranking.

## Empirical input

Supply a CSV with `Date`, `DateTime` and positive `Close` columns, with records ordered chronologically within each date. The historical USD/GBP application uses 2016 and 2020, transforms prices to `100 / Close`, excludes dates with fewer than 10 observations, uses a prespecified 21-dimensional cubic B-spline basis, and selects FPCA using only training observations (95% explained variance, cap 8).

```sh
export JBES_GBPUSD_PATH="/absolute/path/to/GBPUSD.csv"
Rscript reproduce_reported.R --stage empirical_monitoring
Rscript reproduce_reported.R --stage empirical_bootstrap
```

Obtain the input under the provider's license. No raw data are included. Fixed-setting bootstrap agreement does not establish validity after selection of training windows. Event dates are not known change-point labels; earlier empirical alarms do not by themselves establish shorter true detection delays.

## References

- Hong, Y., Lin, Z., Linton, O., Newey, W. K. and Sun, J. (2026). [Affine-Equivariant Adjusted-Range Self-Normalization](https://www.janeway.econ.cam.ac.uk/publication/affine-equivariant-adjusted-range-self-normalization). Cambridge Working Papers in Economics 2678; Janeway Institute Working Paper 2637.
- Hong, Y., Lin, Z., Linton, O., Newey, W. K. and Sun, J. (2026). [aersn](https://rdrr.io/cran/aersn/), version 0.2.3. Use `citation("aersn")` for the package and foundational references. Comparator methods have their own references in the package documentation and manuscript.
- Zhu, M., Hong, Y., Sun, J. and Linton, O. *Online Monitoring of Structural Change with Adjusted-Range Self-Normalization*, revised manuscript.

The manuscript and full-study outputs are maintained separately. Public source availability does not imply that every new full-study run has finished. No public reuse license has been selected.
