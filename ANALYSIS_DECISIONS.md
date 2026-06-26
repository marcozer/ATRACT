# Analysis Decisions

This repository contains the statistical workflow retained for manuscript submission.

## Locked analytic contract

- Analytic dataset: authorized deidentified CSV, SHA-256 `2da50256c21f9198c963fa3042a8818faf762c574054b87358f1487e65f3089c`.
- Primary endpoint: dissection speed in mm2/min.
- Primary speed estimator: 1:1 nearest-neighbor propensity-score matching without replacement.
- Matching constraint: exact matching within pseudonymized operator stratum.
- Calendar time handling: ordered study-year term in the propensity score; temporal proximity evaluated in supplementary diagnostics.
- Balance rule: largest matched sample satisfying maximum absolute SMD <0.10 and no single operator contributing >60% of matched ATRACT procedures.
- Lesion-size heterogeneity: <50 mm, >=50 mm, and interaction effects estimated from one matched interaction model.

## Covariate strategy

The treatment model uses DAG-informed variables plausibly available before or at traction choice:

- pseudonymized operator stratum
- study year
- major lesion diameter with spline
- lesion location
- lesion morphology
- JNET classification
- inflammatory bowel disease history
- recurrence or scar history

Fibrosis is not included in the propensity score because the final degree is not reliably known before traction choice.

## Secondary outcomes

R0 resection, intraprocedural perforation, and delayed bleeding are analyzed with overlap-weighted augmented inverse-probability estimators and reported as risk ratios plus absolute risk differences.

## Supplementary diagnostics

The workflow exports matching balance, propensity overlap, matched-pair calendar gaps, temporal sensitivity, matched-set bootstrap intervals, continuous lesion-size diagnostics, operator-year treatment distributions, crude operator-year speed summaries, matched/unmatched cohort profiles, missingness, and population accounting.

## Patient-level clustering

The public analytic dataset intentionally excludes patient identifiers. A private aggregate diagnostic based on the source workbook found 2002 unique patient identifiers among 2003 cleaned procedures, with one patient contributing two procedures. The primary matched speed cohort contained no repeated patient identifier. This supports the retained matched-set clustered inference for the public analysis; patient identifiers remain unavailable in the public repository.
