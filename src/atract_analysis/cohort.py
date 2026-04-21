from __future__ import annotations

import pandas as pd

from .config import LARGE_LESION_CUTOFF_MM


def load_public_dataset(path: str) -> pd.DataFrame:
    return pd.read_csv(path)


def add_large_lesion_flag(dataframe: pd.DataFrame) -> pd.DataFrame:
    enriched = dataframe.copy()
    enriched["large_lesion"] = (enriched["major_diameter_mm"] >= LARGE_LESION_CUTOFF_MM).astype(int)
    return enriched


def complete_case(dataframe: pd.DataFrame, required_columns: list[str], outcome: str | None = None) -> pd.DataFrame:
    columns = ["atract", *required_columns]
    if outcome is not None:
        columns.append(outcome)
    return dataframe.dropna(subset=columns).copy()

