# Anonymization specification

## Goal

Produce a public analytic dataset that reproduces the statistical analysis while removing direct identifiers, exact temporal information, and free text that could support re-identification.

## Source

- Private source workbook: local `base24.xlsx`
- Source sheet: automatically detected from the workbook
- Public cohort: procedures with `TRACTION == 1`, `ATRACT in {0, 1}`, interpretable anatomic location, and no JNET III classification

## Transformations

- Remove direct identifiers:
  `ID`, `NOM`, `PRENOM`, `Initiales Prenom`
- Remove exact dates:
  `Date examen` is converted to ordered `study_year_index`
- Remove free text:
  `Comments`, `Raison souci technique ATRACT`, raw device comments, and any narrative fields
- Pseudonymize operators:
  `PIOCHE`, `RIVORY`, `LUPU`, and `ROSTAIN` are mapped to `operator_01` to `operator_04`
- Collapse rare operators:
  lower-volume operators and missing operator labels are mapped to `operator_other`
- Top-code age:
  any age above `90` becomes `90`
- Standardize categorical values:
  location, lesion morphology, fibrosis, JNET, CONECCT, sex, inflammatory bowel disease history, and recurrence history are normalized into stable public categories

## Disclosure controls

- No exact dates are distributed
- No direct patient or physician names are distributed
- No free-text values are distributed
- Rare operators are grouped
- The public dataset contains only variables needed for the published analyses

## Intended data-use terms

The deidentified analytic dataset is not publicly distributed at this stage. Access is available upon reasonable request to the study team and remains subject to institutional and governance review. Raw clinical source data are not distributed in this repository.
