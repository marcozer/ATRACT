# Changelog

## 1.0.1-reviewer-diagnostics - 2026-06-26

- Added `REVIEWER_DIAGNOSTICS.md` to map the latest reviewer/statistical comments to generated outputs and manuscript interpretation.
- Added operator-year ATRACT/conventional distributions and crude speed summaries in `table_s10_operator_year_distribution.csv`.
- Added operator adoption-period summaries in `table_s11_operator_adoption.csv`, including explicit identification of `operator_other` as a grouped heterogeneous operator stratum.
- Added temporal diagnostic outputs for the post-adoption period and existing/rematched calendar-proximity analyses.
- Added `table_s15_speed_temporal_relaxation_grid.csv` and `speed_temporal_relaxation_grid.csv` to document matched speed estimates under progressively relaxed calendar-gap constraints.
- Added matched-cohort characteristics by treatment group in `table_s12_matched_characteristics.csv`.
- Added a profile of speed-support procedures not retained by matching in `table_s13_unmatched_profile.csv`.
- Added missingness by treatment group in `table_s14_missingness_by_group.csv`.
- Added aggregate private patient-multiplicity diagnostic in `patient_multiplicity_private_summary.csv` without exposing patient identifiers.
- Added operator-year treatment and crude speed figures in PNG, PDF, and SVG.
- Updated README and reproducibility documentation to explain which outputs answer the reviewer diagnostics.

## 1.0.0-manuscript-submission - 2026-06-22

- Locked the authorized analytic dataset checksum.
- Recomputed lesion-size subgroup estimates from one matched interaction model.
- Added matched-pair calendar diagnostics and temporal sensitivity analyses.
- Added matched-set bootstrap confidence intervals for speed estimates.
- Added continuous lesion-size diagnostic output and figure.
- Added binary outcome event counts, absolute risks, risk differences, and effective sample size.
- Added missingness and population-accounting tables.
- Added manuscript result manifest and enriched reproducibility metadata.
- Removed significance-based regression testing.
