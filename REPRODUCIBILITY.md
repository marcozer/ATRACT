# Reproducibility

## Data access

The row-level deidentified analytic dataset is not distributed in this public repository. It is available upon reasonable request to the study team, subject to institutional governance review.

For the manuscript analysis, place the authorized file at:

```bash
data/public/atract_analysis_public.csv
```

The locked SHA-256 checksum is:

```text
2da50256c21f9198c963fa3042a8818faf762c574054b87358f1487e65f3089c
```

`make check` fails if the local dataset does not match this checksum.

## Commands

Use the `make` targets below for reproduction. They invoke the interpreter and runtime configuration used for the frozen manuscript outputs.

Run the release checks:

```bash
make check
```

Regenerate all manuscript outputs:

```bash
make analysis
```

Run regression tests:

```bash
make test
```

If authorized access to the private source workbook has been granted, regenerate the deidentified analytic dataset first:

```bash
ATRACT_RAW_XLSX=/absolute/path/to/base24.xlsx make public-data
make analysis
```

Generate aggregate patient-multiplicity diagnostics from the private workbook:

```bash
ATRACT_RAW_XLSX=/absolute/path/to/base24.xlsx make private-diagnostics
```

This target exports only non-identifying aggregate counts. Patient identifiers are not written to the public dataset.

## Outputs

Generated outputs are written to:

- `results/tables/`
- `results/figures/`
- `results/model_summaries/`

`results/model_summaries/analysis_metadata.json` records dataset checksum, software versions, propensity-score formula, outcome formulas, selected caliper, matching diagnostics, and bootstrap settings.

`results/model_summaries/manuscript_results_manifest.csv` links reported estimates to generated output files.

`results/model_summaries/patient_multiplicity_private_summary.csv` is generated only when the private workbook is available.
