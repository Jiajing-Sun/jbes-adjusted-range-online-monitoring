# Online Monitoring of Structural Change with Adjusted-Range Self-Normalization

R code accompanying the manuscript by Meiting Zhu, Yongmiao Hong, Jiajing Sun and Oliver Linton. This revision uses **aersn 0.2.3** for adjusted-range geometry, quadratic self-normalization and HAC normalization.

All R files are directly in this directory. There are no datasets, numerical results, figures, PDFs, archives or compiled binaries. Scripts create local `results/`, `figures/` and `logs/` directories; the Git allowlist excludes these outputs. No variance-profile experiment is included.

## Install and verify

Use R 4.1 or later. Run from the repository root:

```sh
Rscript install_dependencies.R
Rscript verify_backend.R
Rscript verify_par_indexing.R
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
| `50_empirical_utils.R`, `31_empirical_plots.R`, `32_empirical_monitoring_table.R`, `34_empirical_bootstrap_check.R`, `35_empirical_calibration_details.R` | USD/GBP preprocessing, monitoring and bootstrap |
| `30_scalar_theory_checks.R`, `31_bivariate_theory_checks.R` | Independent same-contaminated-null calibration and checks of the new scalar and directional propositions |
| `33_dgp_exploration.R` | Complete DGP screening grid and independent confirmation with automatic HAC controls |
| `50_focused_dgp_search.R` | Follow-up model, weight and direction screening with a frozen independent confirmation |
| `51_containment_reference_audit.R` | All 400 necessary threshold-bound checks for uniform hull containment |
| `52_focused_dgp_report.R` | Complete frozen-candidate results and paired differences, including automatic HAC controls |
| `34_dgp_comparison_report.R` | Confirmation report for scenarios frozen by screening, retaining automatic-HAC comparisons and null checks |
| `32_theory_assets.R` | Tables and plots from saved theory-check outputs |
| `08_mc_audit.R` | First-crossing counts and paired Monte Carlo uncertainty |
| Other numbered R files | Local tables, plots and empirical surface coordinates |

## Statistical conventions

- Training estimates are frozen before monitoring. Online maxima and stopping rules are implemented here. Package offline p-values never replace sequential critical values.
- Multivariate RSMS uses the full increment hull, without a componentwise replacement or a hull ridge.
- Scalar gauges are analytical. Batched two-dimensional gauges use vertices from public aersn functions. Above dimension two, a reusable HiGHS dual LP evaluates the aersn hull. Its first objective for each hull is compared with `aersn_gauge()`, and each accelerated solve checks primal/dual feasibility and objective gap. Set `JBES_AERSN_GAUGE=reference` when invoking a component script to use the public package LP throughout.
- Screening certifies intervals on the original grid; it does not thin the grid. This source-only repository uses portable R loops; no private compiled accelerator is required.
- Default Bartlett lag is `L = ceiling(1.06 * m^(1/5))`, with package bandwidth `h = L + 1`. SSMS and CSMS retain the documented `1e-8` ridge. QS uses all lags: regenerate its sensitivity results rather than reusing older truncated output.
- PAR training length counts transitions. The returned observations start with X_0, and a change begins at monitored transition k_star+1. Earlier alternatives began one transition too early; regenerate their power and delay results. A failed conditional MLE stops the experiment rather than dropping the path.
- If a rounding-sized edge makes a bivariate difference-polygon facet unusable, the adapter evaluates the same hull through the support widths of its training-path polygon. It records the recovery and checks the public package gauge; it does not drop or replace the path.
- Correlated innovations use `t(chol(Sigma))` for column-vector draws, giving the declared covariance `Sigma`. Earlier code used the upper factor directly. Regenerate DGP2, DGP4 and sensitivity experiments that use correlated innovations; earlier output from that implementation is not interchangeable with these results.
- Independent calibration, null evaluation and alternative evaluation serve distinct purposes. Size-adjusted power is a model-specific benchmark, not an operational claim that the true data-generating process is known. Retain CSMS advantages and ranking reversals; do not choose tuning rules to lower competitor power.
- Numerical failures abort a run; paths are not silently discarded or replaced. An aborted run is not a completed Monte Carlo result.

The theory-check scripts use separate calibration, null-evaluation and alternative streams. Defaults are 10,000 calibration and 5,000 evaluation replications. Both signs of local changes and the bivariate transverse directions are retained. The scalar fixed-change delay experiment is separate from its local-power experiment. Paired evaluation errors condition on estimated thresholds; calibration-quantile brackets are reported separately. These scripts implement the specified models and do not establish a universal ranking.

## DGP screening and confirmation

`dgp_screen` evaluates a fixed grid of 26 designs covering positive/negative serial correlation, seasonal dependence, heavy tails, GARCH, cross-lag VAR dynamics and persistent training changes. Stationary noise has long-run covariance I_2, so local shifts use a common long-run signal scale. `dgp_screen_q1` runs their first-coordinate marginals, each with long-run variance one. The scalar batch has separate random streams and output files. Calibration and both evaluation phases are independent.

`dgp_confirmation` uses new streams, two training lengths and automatic Andrews-QS and Newey-West-Bartlett controls for the three highest-ranked screening designs. All screening settings, null checks and reversals remain in the generated output. Screening scores are exploratory, not multiplicity-adjusted significance tests; confirmation results must be assessed separately. Model-specific calibration is an oracle benchmark and does not assume the true DGP is known in an application.

```sh
Rscript reproduce_reported.R --stage dgp_screen
Rscript reproduce_reported.R --stage dgp_confirmation
Rscript reproduce_reported.R --stage dgp_screen_q1
Rscript reproduce_reported.R --stage dgp_confirmation_q1
```

A follow-up batch, `dgp_screen_profile`, evaluates ten asymmetric perturbations of the full-rank training-mean profile in the paper. This batch was defined after the stationary screen. It uses iid noise, five unit-length change directions, separate seeds and a separate confirmation stage; it does not vary an information clock. All profiles and transverse signs remain in the output.

```sh
Rscript reproduce_reported.R --stage dgp_screen_profile
Rscript reproduce_reported.R --stage dgp_confirmation_profile
Rscript reproduce_reported.R --stage dgp_report
```

The additional `dgp_focused_screen` stage was designed after inspecting the earlier results. It crosses seven models with two monitoring weights, including an iid reference, two stable triangular VAR models and four full-rank training-mean profiles. Each setting evaluates nine directions, three local effect magnitudes and two change locations. All five methods, including both automatic HAC controls, participate in screening. The top three distinct model/weight settings and their exact alternatives are frozen for `dgp_focused_confirmation`, using fresh calibration and evaluation draws at both training lengths. The full follow-up grid remains available; it must not be relabeled as a prespecified confirmation of the earlier search.

```sh
Rscript reproduce_reported.R --stage dgp_focused_screen
Rscript reproduce_reported.R --stage dgp_focused_confirmation
Rscript 52_focused_dgp_report.R
```

The optional `JBES_DGP_REP_OFFSET` in `50_focused_dgp_search.R` shifts each phase to disjoint replication seed indices for an independent repeat. Use a new `JBES_FOCUS_OUTPUT` prefix and copy the frozen `selected_dgps.csv` into its `_screen` directory; never reuse an existing manifest with changed counts or offsets. The fixed profile follow-up uses ID `profile_L16_g0.00`, `m=500`, offset 30000, 20,000 calibration paths and 10,000 paths per evaluation phase. Its alternative remains the one frozen in the original selection file.

The focused search uses 1,000 calibration, 1,000 independent null and 1,000 alternative paths per screening setting. Confirmation uses 10,000 calibration, 5,000 null and 5,000 alternative paths. Model-specific thresholds are oracle simulation comparisons. Every method uses the same evaluation paths; first crossings distinguish pre-change alarms, post-change detections and no crossing. Paired errors condition on the estimated thresholds. Phase checkpoints include source hashes, design, replication counts and package versions and refuse inconsistent reuse.

PAR null and alternative, moderate-dimensional and persistent-training checkpoints include design and implementation manifests. The PAR alternative cache uses a versioned directory so earlier files without manifests remain intact and are not silently reused. A mismatch stops the run instead of reusing stale paths. Mean innovations use compiled scalar recursions in an orthogonal eigenbasis; direct-recursion checks confirm the same random draws and agreement to numerical precision.

The optional `JBES_RESUME_VERIFIED=1` switch in the mean and bivariate scripts is for locally saved checkpoints only after the design, seeds, implementation and thresholds have been verified unchanged. It is off by default. Do not reuse arbitrary earlier output.

## Empirical input

Supply a CSV with `Date`, `DateTime` and positive `Close` columns, with records ordered chronologically within each date. The historical USD/GBP application uses 2016 and 2020, transforms prices to `100 / Close`, excludes dates with fewer than 10 observations, uses a fixed 21-dimensional cubic B-spline basis, and selects FPCA using only training observations (95% explained variance, cap 8).

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

Run `Rscript 51_containment_reference_audit.R` after generating the finite and open-end references to reproduce the per-dimension threshold checks. A passed necessary bound is not evidence of attainable containment.

Run `Rscript 53_hac_record_audit.R` after the mean, bandwidth, kernel and automatic-bandwidth stages to reconcile their saved first-crossing records and threshold quantiles. The fixed-tuning CSMS samples are independent of the baseline RSMS samples, so this report uses independent-sample comparison errors. It preserves each method's own post-change detection event for CADD and conditions on estimated thresholds.
The same audit reports paired errors for the automatic selectors, which share evaluation paths. Fixed-tuning and automatic-selector comparison errors are therefore kept separate.
