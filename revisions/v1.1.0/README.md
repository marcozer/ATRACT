# Revision analyses, v1.1.0

Additional analyses for the manuscript revision. The primary results archived
in v1.0.0 remain in the root `results/` directory and are checked against 61
frozen file checksums before and after each revision run.

From the repository root:

```sh
make revision-analysis
make revision-test
```

The authorized analytic CSV is required for estimation and remains available
upon request. No row-level dataset is included in this release. To use an
authorized CSV at another location:

```sh
python revisions/v1.1.0/analysis.py --input /path/to/atract_analysis_public.csv
```

`results/tables/estimates_summary.csv` collects the point estimates, confidence
intervals, p-values, model populations and wider-caliper balance diagnostics.
`RESULTS_fr.md` describes the requested comparisons. `METHODS.md` contains the
additional methods and figure legends. Every figure is exported as SVG, PNG
and PDF using the Nord palette.

## Specifications

- The original 391-pair matching allocation and the treatment model are reproduced.
- All fibrosis specifications are retained: parsimonious and fully adjusted,
  with and without fibrosis on the same 372 complete pairs. The full baseline
  adjustment is also estimated in the original 391 pairs to distinguish
  adjustment from complete-case restriction.
- Operator-volume groups are defined by the study team's clinical attestation
  of 13 October 2026, not by outcome-based clustering or a reconstructed
  learning curve. `operator_01` and `operator_02` perform approximately 300 ESDs
  annually and have more than 10 years of ESD experience. The remaining
  operators are trained, perform approximately 50–100 ESDs annually, and have
  approximately 5 years of experience. This is a fixed grouping, not an annual
  measurement of experience at each procedure.
- Operator-volume effects are estimated from one interaction model adjusted
  for all pretreatment covariates, including operator fixed effects. A second
  specification excludes the pooled operator stratum. Balance is reported
  within the volume groups; this analysis does not address beginners.
- The four-degree-of-freedom continuous-size interaction uses the v1.0.0
  spline model. The curve is checked against archived Table S9. Its joint
  Wald test uses matched-set clustered covariance. Curve intervals are
  pointwise, rather than simultaneous.
- Procedure duration is estimated directly, adjusting for year and specimen
  area. Minutes are not inferred by inverting a speed difference.
- Wider calipers are fixed at 0.40 and 0.60 SD of the logit propensity score.
  Both specifications are reported regardless of their results. They do not
  satisfy the primary maximum absolute SMD <0.10 criterion and therefore
  remain exploratory.

The archived label `PS-NN` identifies the original analysis. Its allocation
uses a global minimum-distance assignment within pseudonymized operator
strata under a propensity-logit caliper, rather than sequential greedy
nearest-neighbour selection. These analyses preserve that implementation.

## Reading the summary

`n` is the number of procedures in the reported stratum. `n_model` and
`n_pairs_model` describe the full interaction-model population. Size-specific
effects can involve pairs whose members fall on opposite sides of 50 mm.
An interaction estimate is a difference between stratum effects. The
continuous-size test row reports a Wald statistic, so a confidence interval
is not applicable to that row. Binary rows reproduce archived risk ratios;
their risk differences and event counts remain in `secondary_outcomes.csv`.

Inference uses matched-set clustered sandwich covariance for matched models
and the archived empirical influence variance for overlap augmentation. The
propensity scores are held fixed in the complete-case overlap sensitivities.
All additional analyses are exploratory; no specification is selected by its
p-value.
