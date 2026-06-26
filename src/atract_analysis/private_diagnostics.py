from __future__ import annotations

from pathlib import Path

import pandas as pd

from .anonymize import (
    _resolve_location_series,
    _resolve_raw_location_series,
    _resolve_surface_series,
    _resolve_speed_series,
    _to_numeric,
)
from .config import (
    EXCLUDED_JNET_GROUPS,
    EXCLUDED_RAW_LOCATION_CODES,
    JNET_MAP,
    PHYSICIAN_PUBLIC_MAP,
    RAW_COLUMNS,
)
from .ingest import load_raw_workbook, normalize_missing_values


PRIVATE_PATIENT_ID_COLUMN = "ID"


def _private_patient_frame(raw_path: str | Path) -> pd.DataFrame:
    raw = normalize_missing_values(load_raw_workbook(raw_path))
    raw_location = _resolve_raw_location_series(raw)
    dataframe = raw.loc[
        _to_numeric(raw[RAW_COLUMNS["traction"]]).eq(1)
        & _to_numeric(raw[RAW_COLUMNS["atract"]]).isin([0, 1])
        & ~raw_location.isin(EXCLUDED_RAW_LOCATION_CODES)
    ].copy()
    location = _resolve_location_series(dataframe)
    dataframe = dataframe.loc[location.notna()].copy()

    # The public dataset deliberately drops the patient identifier. This private
    # helper keeps it only long enough to produce aggregate clustering checks.
    exam_dates = pd.to_datetime(dataframe[RAW_COLUMNS["date_exam"]], errors="coerce")
    year_order = {
        year: index + 1
        for index, year in enumerate(sorted(int(year) for year in exam_dates.dropna().dt.year.unique()))
    }
    surface = _resolve_surface_series(dataframe)
    private = pd.DataFrame(
        {
            "_patient_id_private": dataframe[PRIVATE_PATIENT_ID_COLUMN].astype("string"),
            "atract": _to_numeric(dataframe[RAW_COLUMNS["atract"]]).astype("Int64"),
            "study_year_index": exam_dates.dt.year.map(year_order).astype("Int64"),
            "operator_id_public": dataframe[RAW_COLUMNS["physician"]]
            .map(PHYSICIAN_PUBLIC_MAP)
            .fillna("operator_other")
            .astype("string"),
            "major_diameter_mm": _to_numeric(dataframe[RAW_COLUMNS["major_diameter"]]),
            "surface_mm2": surface,
            "speed_mm2_min": _resolve_speed_series(dataframe, surface),
            "jnet_group": dataframe[RAW_COLUMNS["jnet"]].map(JNET_MAP).astype("string"),
        }
    )
    private = private.loc[~private["jnet_group"].isin(EXCLUDED_JNET_GROUPS)].copy()
    return private.sort_values(
        by=["study_year_index", "operator_id_public", "atract", "major_diameter_mm", "surface_mm2"],
        kind="mergesort",
    ).reset_index(drop=True)


def _multiplicity_row(label: str, patient_ids: pd.Series) -> dict[str, int | str]:
    observed = patient_ids.dropna().astype("string")
    counts = observed.value_counts()
    repeated = counts.loc[counts.gt(1)]
    return {
        "population": label,
        "rows_n": int(len(patient_ids)),
        "patient_id_missing_n": int(patient_ids.isna().sum()),
        "unique_patient_ids_n": int(counts.size),
        "repeated_patient_ids_n": int(repeated.size),
        "rows_in_repeated_patient_ids_n": int(repeated.sum()) if not repeated.empty else 0,
        "max_rows_per_patient_id": int(counts.max()) if not counts.empty else 0,
    }


def build_patient_multiplicity_summary(
    raw_path: str | Path,
    *,
    public_dataframe: pd.DataFrame,
    primary_speed: dict[str, object],
) -> pd.DataFrame:
    private = _private_patient_frame(raw_path)
    if len(private) != len(public_dataframe):
        raise ValueError(
            f"Private patient frame has {len(private)} rows but public dataset has {len(public_dataframe)} rows."
        )

    scored = primary_speed["scored_frame"].copy()
    matched = primary_speed["matched_frame"].copy()

    scored_ids = private.loc[scored["_analysis_row_id"].astype(int), "_patient_id_private"]
    matched_ids = private.loc[matched["_analysis_row_id"].astype(int), "_patient_id_private"]

    return pd.DataFrame(
        [
            _multiplicity_row("cleaned_public_cohort", private["_patient_id_private"]),
            _multiplicity_row("speed_supported_cohort", scored_ids),
            _multiplicity_row("primary_matched_speed_cohort", matched_ids),
        ]
    )


def write_private_diagnostics(
    raw_path: str | Path,
    output_path: str | Path,
    *,
    public_dataframe: pd.DataFrame,
    primary_speed: dict[str, object],
) -> pd.DataFrame:
    summary = build_patient_multiplicity_summary(
        raw_path,
        public_dataframe=public_dataframe,
        primary_speed=primary_speed,
    )
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(output_path, index=False)
    return summary
