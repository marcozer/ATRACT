# ATRACT Public Statistical Analysis

This repository contains the public analysis code and reproducible outputs for the ATRACT colorectal ESD cohort. The release is intentionally limited to the analyses retained for the manuscript and supplement.

## Repository scope

- Python package for ingestion, anonymization, cohort definition, primary causal analyses, tables, figures, and release checks
- Public analysis code and derived results
- Reproducible outputs in `results/`
- Metadata describing the anonymization and DAG-driven covariate strategy

## Related resources

- ATRACT device website: [the-atract-device.fr](https://www.atract-device.fr/the-atract-device)
- Repository citation metadata: `CITATION.cff`

## Reproducibility model

The public repository distributes the full analysis code and the derived manuscript outputs. Reproducing the analysis from row-level data currently requires authorized access to the deidentified analytic dataset or the source workbook.

```bash
make analysis
make test
```

The deidentified analytic dataset is not distributed in this public repository. At this stage, access is available upon reasonable request to the study team and subject to local governance constraints. The manuscript analysis is locked to the authorized analytic dataset with SHA-256 checksum `2da50256c21f9198c963fa3042a8818faf762c574054b87358f1487e65f3089c`; `make check` fails if a different dataset is supplied.

If you have authorized access to the source workbook, you can regenerate the public dataset locally:

```bash
ATRACT_RAW_XLSX=/absolute/path/to/base24.xlsx make public-data
```

Then rerun the analysis:

```bash
make all
```

## Layout

- `src/atract_analysis/`: analysis code
- `data/public/`: local location for the deidentified analytic dataset when access has been granted
- `metadata/`: data dictionary, anonymization specification, DAGitty model
- `results/tables/`: manuscript-facing tables in CSV
- `results/figures/`: manuscript-facing figures in PNG, PDF, and SVG
- `results/model_summaries/`: machine-readable model outputs
- `tests/`: regression and release checks

## Analytic dataset

Each row of the analytic dataset represents one colorectal ESD procedure from the institutional registry after public-data cleaning. The dataset itself is not public in the current repository and is available upon request only.

The main groups are:

- `atract = 1`: procedure performed with the adaptive ATRACT traction device
- `atract = 0`: procedure performed with a conventional traction method

The public dataset keeps only the variables needed to understand case complexity, treatment allocation, and outcomes. Examples include:

- lesion size and location
- lesion morphology and optical classification
- inflammatory bowel disease and recurrence/scar history
- pseudonymized operator identifier
- procedure duration and derived dissection speed
- R0 resection, perforation, and delayed bleeding

The full list and source mapping are in `metadata/public_data_dictionary.csv`.

## Public dataset notes

- Exact dates are replaced with `study_year_index`
- Operators are pseudonymized as `operator_XX`
- Direct identifiers and free text are removed
- Age is top-coded at `90`
- Rare operators are collapsed into `operator_other`
- JNET III lesions are excluded from the analytic cohort

## Statistical rationale

The study compares adaptive traction (`ATRACT`) with conventional traction in an observational cohort. Treatment allocation was therefore not randomized and depended on operator, calendar period, and lesion complexity. For this reason, the repository does not rely on crude group comparisons or on a single adjusted regression model.

The analysis is structured as a DAG-informed, outcome-blind causal inference workflow defined before final model estimation. Covariates were selected from a directed acyclic graph and restricted to variables plausibly known before or at the time of treatment choice. The core adjustment set includes:

- operator
- calendar time
- lesion size
- lesion location
- lesion morphology
- JNET classification
- inflammatory bowel disease history
- recurrence or scar history

These variables are documented in `metadata/dagitty_model.txt`.

Two complementary estimators are retained in the public workflow.

### 1. Primary analysis: matched comparable-case analysis (`PS-NN`)

The primary analysis for dissection speed uses `1:1` nearest-neighbor propensity score matching without replacement.

- The propensity score models the probability of receiving ATRACT conditional on the DAG-informed pre-treatment covariates.
- Matching is exact on pseudonymized operator stratum; calendar time is included in the propensity score model as an ordered study-year term.
- The final caliper is selected from a locked grid using covariate balance and retained sample size, not significance testing.
- The retained matched sample is the largest one satisfying the locked balance criterion used in the manuscript.
- Effect estimation in the matched cohort uses a parsimonious post-match model adjusted for calendar year and lesion size, with matched-set clustered standard errors.
- Lesion-size subgroup estimates and the interaction estimate are derived from one interaction model.

This estimator was chosen as the primary analysis because it targets the comparison of clinically comparable treated and control procedures while keeping the design transparent and auditable.

### 2. Robustness analysis: overlap-weighted augmented inverse-probability analysis (`OW + DR`)

The robustness analysis for dissection speed, and the primary analyses for binary secondary outcomes, use overlap-weighted augmented inverse-probability estimators.

- Overlap weighting emphasizes procedures with empirical support in both treatment groups and downweights observations at the tails of the propensity score distribution.
- The augmented specification combines the treatment-assignment model with an outcome-regression augmentation.
- For dissection speed, the augmentation model uses the DAG-informed pre-treatment prognostic covariates rather than a post-hoc variable screen.
- Supplementary speed robustness outputs include overlap-weighted estimates by lesion-size category and a restricted analysis among lesions measuring at least 50 mm.
- For binary endpoints, the outcome model additionally adjusts for age, sex, ASA class, anticoagulant use, and antiplatelet use.

This complementary estimator preserves a broader supported cohort and provides a robustness check under a different causal design.

## Analysis overview

The committed workflow implements:

1. Locked broad-cohort construction from the deidentified analytic dataset
2. Baseline descriptive table for ATRACT versus non-ATRACT
3. Minimal broad-cohort cleanup with exclusion of uninterpretable locations and raw location `10`
4. Support restriction to operators with empirical data in both treatment groups
5. Primary `1:1` nearest-neighbor propensity-score matching analysis for `speed_mm2_min`, exact on pseudonymized operator stratum with calendar time handled in the PS
6. Overlap-weighted augmented inverse-probability robustness analysis for `speed_mm2_min`
7. Primary overlap-weighted augmented analyses for `r0`, `perforation`, and `delayed_bleeding`, reported as risk ratios and absolute risk differences
8. Supplementary diagnostics: DAG, balance tables, overlap plot, matching grid, temporal sensitivity, matched-set bootstrap, continuous-size diagnostic, operator-year distributions, matched/unmatched cohort profiles, missingness, and result manifest

## How to read the outputs

- `results/tables/table_1_baseline.csv`:
  descriptive characteristics of the overall cohort and the two treatment groups
- `results/tables/table_2_main_results.csv`:
  main manuscript results, including the primary matched speed analysis and the primary weighted binary analyses
- `results/tables/table_s1_balance.csv`:
  covariate balance before and after matching
- `results/tables/table_s2_speed_robustness.csv`:
  broad-cohort robustness analysis for dissection speed
- `results/tables/table_s3_secondary_comparison.csv`:
  supplementary comparison table for binary outcomes
- `results/tables/table_s4_missingness.csv` and `table_s5_population_accounting.csv`:
  missingness and analysis-population accounting
- `results/tables/table_s10_operator_year_distribution.csv` and `table_s11_operator_adoption.csv`:
  ATRACT/conventional distributions and crude speed summaries by operator and study year
- `results/tables/table_s12_matched_characteristics.csv` and `table_s13_unmatched_profile.csv`:
  characteristics of the primary matched cohort and the speed-support procedures not retained by matching
- `results/tables/table_s14_missingness_by_group.csv`:
  missingness by treatment group
- `results/tables/table_s15_speed_temporal_relaxation_grid.csv`:
  matched speed estimates under progressively relaxed calendar-gap constraints
- `results/model_summaries/manuscript_results_manifest.csv`:
  machine-readable map from reported estimates to generated outputs and dataset checksum

`results/model_summaries/patient_multiplicity_private_summary.csv` is an aggregate, non-identifying diagnostic generated from the private source workbook with `make private-diagnostics`. It is not generated from the public analytic CSV because patient identifiers are intentionally excluded from public data.

The main result to look at is the matched analysis for `speed_mm2_min`. The weighted analyses for `R0`, `perforation`, and `delayed_bleeding` are secondary.

## Data and code licensing

- Code: `MIT`, see `LICENSE`
- Deidentified analytic dataset: not publicly distributed at this stage; available upon request subject to institutional and governance review
