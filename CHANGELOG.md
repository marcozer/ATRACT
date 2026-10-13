# Changelog

## v1.1.0 - 2026-10-13

- Added revision analyses in `revisions/v1.1.0/` without changing any of the 61 archived v1.0.0 output files.
- Retained all fibrosis specifications and exported comparisons with and without fibrosis on identical complete pairs, plus a full-covariate model in the original matched cohort.
- Added the documented high-volume versus lower-volume trained operator comparison, with interaction inference and within-group balance diagnostics.
- Reproduced the archived continuous-size curve and added its four-degree-of-freedom joint Wald interaction test.
- Exported directly estimated procedure-duration differences and confidence intervals.
- Exported effects and confidence intervals for calipers 0.40 and 0.60, alongside their failed primary balance criterion.
- Added one consolidated estimates CSV, publication-ready figure formats, methods, legends and reproducibility tests.
- Updated release citation metadata; v1.0.0 and its DOI remain unchanged.

## v1.0.0 - 2026-06-30

- Removed internal review-facing documentation from the public repository.
- Kept supplementary diagnostic tables and figures as manuscript-facing analysis outputs.
- Updated citation metadata for the cleaned Zenodo archive.

## Supplementary diagnostics - 2026-06-26

- Added operator-year ATRACT/conventional distributions and crude speed summaries in `table_s10_operator_year_distribution.csv`.
- Added operator adoption-period summaries in `table_s11_operator_adoption.csv`, including explicit identification of `operator_other` as a grouped heterogeneous operator stratum.
- Added temporal diagnostic outputs for the post-adoption period and existing/rematched calendar-proximity analyses.
- Added `table_s15_speed_temporal_relaxation_grid.csv` and `speed_temporal_relaxation_grid.csv` to document matched speed estimates under progressively relaxed calendar-gap constraints.
- Added matched-cohort characteristics by treatment group in `table_s12_matched_characteristics.csv`.
- Added a profile of speed-support procedures not retained by matching in `table_s13_unmatched_profile.csv`.
- Added missingness by treatment group in `table_s14_missingness_by_group.csv`.
- Added aggregate private patient-multiplicity diagnostic in `patient_multiplicity_private_summary.csv` without exposing patient identifiers.
- Added operator-year treatment and crude speed figures in PNG, PDF, and SVG.
- Updated README and reproducibility documentation to explain the supplementary diagnostic outputs.

## Manuscript analysis hardening - 2026-06-22

- Locked the authorized analytic dataset checksum.
- Recomputed lesion-size subgroup estimates from one matched interaction model.
- Added matched-pair calendar diagnostics and temporal sensitivity analyses.
- Added matched-set bootstrap confidence intervals for speed estimates.
- Added continuous lesion-size diagnostic output and figure.
- Added binary outcome event counts, absolute risks, risk differences, and effective sample size.
- Added missingness and population-accounting tables.
- Added manuscript result manifest and enriched reproducibility metadata.
- Removed significance-based regression testing.
