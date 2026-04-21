from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .config import (
    CONECCT_MAP,
    DATA_DICTIONARY_ROWS,
    EXPERIENCE_MAP,
    EXCLUDED_RAW_LOCATION_CODES,
    FIBROSIS_MAP,
    JNET_MAP,
    LESION_TYPE_MAP,
    LOCATION_MAP,
    PHYSICIAN_PUBLIC_MAP,
    PUBLIC_COLUMNS,
    RAW_LOCATION_TO_SIMPLIFIED_MAP,
    RAW_COLUMNS,
    SEX_MAP,
)
from .ingest import load_raw_workbook, normalize_missing_values


def _to_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def _to_binary(series: pd.Series) -> pd.Series:
    numeric = _to_numeric(series)
    return numeric.astype("Int64")


def _resolve_location_series(dataframe: pd.DataFrame) -> pd.Series:
    if RAW_COLUMNS["location_simple"] in dataframe.columns:
        return _to_numeric(dataframe[RAW_COLUMNS["location_simple"]])
    if RAW_COLUMNS["location_raw"] in dataframe.columns:
        raw_location = _to_numeric(dataframe[RAW_COLUMNS["location_raw"]])
        simplified = raw_location.map(RAW_LOCATION_TO_SIMPLIFIED_MAP)
        return simplified
    raise KeyError("No supported location column found in the raw workbook.")


def _resolve_raw_location_series(dataframe: pd.DataFrame) -> pd.Series:
    if RAW_COLUMNS["location_raw"] in dataframe.columns:
        return _to_numeric(dataframe[RAW_COLUMNS["location_raw"]])
    if RAW_COLUMNS["location_simple"] in dataframe.columns:
        return _to_numeric(dataframe[RAW_COLUMNS["location_simple"]])
    raise KeyError("No supported raw location column found in the raw workbook.")


def _resolve_operator_experience_series(dataframe: pd.DataFrame) -> pd.Series:
    if RAW_COLUMNS["operator_simple"] in dataframe.columns:
        return dataframe[RAW_COLUMNS["operator_simple"]].map(EXPERIENCE_MAP).astype("string")
    return pd.Series(pd.NA, index=dataframe.index, dtype="string")


def _resolve_surface_series(dataframe: pd.DataFrame) -> pd.Series:
    if RAW_COLUMNS["surface"] in dataframe.columns:
        return _to_numeric(dataframe[RAW_COLUMNS["surface"]])
    major = _to_numeric(dataframe[RAW_COLUMNS["major_diameter"]])
    minor = _to_numeric(dataframe[RAW_COLUMNS["minor_diameter"]])
    return np.pi * major * minor / 4


def _resolve_speed_series(dataframe: pd.DataFrame, surface_series: pd.Series) -> pd.Series:
    if RAW_COLUMNS["speed"] in dataframe.columns:
        return _to_numeric(dataframe[RAW_COLUMNS["speed"]])
    duration = _to_numeric(dataframe[RAW_COLUMNS["duration"]])
    speed = surface_series / duration
    speed = speed.where(duration.ne(0))
    return speed


def _normalize_recurrence(history_text: pd.Series, history_binary: pd.Series) -> pd.Series:
    text = (
        history_text.fillna("")
        .astype(str)
        .str.lower()
        .str.normalize("NFKD")
        .str.encode("ascii", errors="ignore")
        .str.decode("ascii")
    )
    binary = history_binary.fillna("").astype(str).str.lower()
    recurrence = (
        text.str.contains("recid", na=False)
        | binary.str.contains("1-", na=False)
        | binary.str.contains("oui", na=False)
    )
    return recurrence.astype("Int64")


def build_public_dataset(raw_dataframe: pd.DataFrame) -> pd.DataFrame:
    dataframe = normalize_missing_values(raw_dataframe)
    raw_location = _resolve_raw_location_series(dataframe)
    dataframe = dataframe.loc[
        _to_numeric(dataframe[RAW_COLUMNS["traction"]]).eq(1)
        & _to_numeric(dataframe[RAW_COLUMNS["atract"]]).isin([0, 1])
        & ~raw_location.isin(EXCLUDED_RAW_LOCATION_CODES)
    ].copy()
    location_series = _resolve_location_series(dataframe)
    interpretable_location_mask = location_series.notna()
    dataframe = dataframe.loc[interpretable_location_mask].copy()
    location_series = location_series.loc[interpretable_location_mask]

    exam_dates = pd.to_datetime(dataframe[RAW_COLUMNS["date_exam"]], errors="coerce")
    year_order = {
        year: index + 1
        for index, year in enumerate(sorted(int(year) for year in exam_dates.dropna().dt.year.unique()))
    }

    operator_public = (
        dataframe[RAW_COLUMNS["physician"]]
        .map(PHYSICIAN_PUBLIC_MAP)
        .fillna("operator_other")
    )
    operator_experience = _resolve_operator_experience_series(dataframe)
    surface_series = _resolve_surface_series(dataframe)
    speed_series = _resolve_speed_series(dataframe, surface_series)

    public_dataframe = pd.DataFrame(
        {
            "atract": _to_binary(dataframe[RAW_COLUMNS["atract"]]),
            "study_year_index": exam_dates.dt.year.map(year_order).astype("Int64"),
            "operator_id_public": operator_public.astype("string"),
            "operator_experience": operator_experience,
            "age_years_topcoded": _to_numeric(dataframe[RAW_COLUMNS["age"]]).clip(upper=90).round().astype("Int64"),
            "sex": dataframe[RAW_COLUMNS["sex"]].map(SEX_MAP).astype("string"),
            "asa": _to_numeric(dataframe[RAW_COLUMNS["asa"]]).round().astype("Int64"),
            "anticoagulants": _to_binary(dataframe[RAW_COLUMNS["anticoagulants"]]),
            "antiplatelets": _to_binary(dataframe[RAW_COLUMNS["antiplatelets"]]),
            "location_group": location_series.map(LOCATION_MAP).astype("string"),
            "major_diameter_mm": _to_numeric(dataframe[RAW_COLUMNS["major_diameter"]]),
            "minor_diameter_mm": _to_numeric(dataframe[RAW_COLUMNS["minor_diameter"]]),
            "surface_mm2": surface_series,
            "lesion_type": dataframe[RAW_COLUMNS["lesion_type"]].map(LESION_TYPE_MAP).fillna("protruding_or_other").astype("string"),
            "macronodule": _to_binary(dataframe[RAW_COLUMNS["macronodule"]]),
            "jnet_group": dataframe[RAW_COLUMNS["jnet"]].map(JNET_MAP).astype("string"),
            "conecct_group": dataframe[RAW_COLUMNS["conecct"]].map(CONECCT_MAP).astype("string"),
            "fibrosis": dataframe[RAW_COLUMNS["fibrosis"]].map(FIBROSIS_MAP).astype("string"),
            "recurrence_history": _normalize_recurrence(
                dataframe[RAW_COLUMNS["history_text"]],
                dataframe[RAW_COLUMNS["history_binary"]],
            ),
            "procedure_duration_min": _to_numeric(dataframe[RAW_COLUMNS["duration"]]),
            "speed_mm2_min": speed_series,
            "r0": _to_binary(dataframe[RAW_COLUMNS["r0"]]),
            "perforation": _to_binary(dataframe[RAW_COLUMNS["perforation"]]),
            "delayed_bleeding": _to_binary(dataframe[RAW_COLUMNS["bleeding"]]),
            "curative_resection": _to_binary(dataframe[RAW_COLUMNS["curative_resection"]]),
        }
    )[PUBLIC_COLUMNS]

    public_dataframe = public_dataframe.sort_values(
        by=["study_year_index", "operator_id_public", "atract", "major_diameter_mm", "surface_mm2"],
        kind="mergesort",
    ).reset_index(drop=True)

    return public_dataframe


def build_data_dictionary() -> pd.DataFrame:
    dictionary = pd.DataFrame(
        DATA_DICTIONARY_ROWS,
        columns=["column_name", "dtype", "description", "source_field"],
    )
    return dictionary


def write_public_release(raw_path: str | Path, output_path: str | Path, dictionary_path: str | Path) -> pd.DataFrame:
    raw_dataframe = load_raw_workbook(raw_path)
    public_dataframe = build_public_dataset(raw_dataframe)

    output_path = Path(output_path)
    dictionary_path = Path(dictionary_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    dictionary_path.parent.mkdir(parents=True, exist_ok=True)

    public_dataframe.to_csv(output_path, index=False)
    build_data_dictionary().to_csv(dictionary_path, index=False)
    return public_dataframe
