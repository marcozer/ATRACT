from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from .config import BANNED_STRING_TOKENS, PUBLIC_COLUMNS


DATE_PATTERN = re.compile(r"\b(?:19|20)\d{2}-\d{2}-\d{2}\b")


def validate_public_dataset(dataframe: pd.DataFrame) -> list[str]:
    errors: list[str] = []

    if list(dataframe.columns) != PUBLIC_COLUMNS:
        errors.append("Public dataset columns do not match the locked public schema.")

    if dataframe.empty:
        errors.append("Public dataset is empty.")

    if "study_year_index" in dataframe.columns:
        year_values = dataframe["study_year_index"].dropna()
        if (year_values < 1).any() or not year_values.mod(1).eq(0).all():
            errors.append("study_year_index contains invalid values.")

    for column in dataframe.select_dtypes(include=["object", "string"]).columns:
        series = dataframe[column].dropna().astype(str)
        if series.str.contains(DATE_PATTERN).any():
            errors.append(f"Column {column} appears to contain an exact date.")
        for token in BANNED_STRING_TOKENS:
            if series.str.contains(token, case=False, regex=False).any():
                errors.append(f"Column {column} contains banned identifier token: {token}.")
        if series.str.len().gt(40).any():
            errors.append(f"Column {column} contains unusually long free text values.")

    return errors


def run_release_checks(path: str | Path) -> None:
    dataframe = pd.read_csv(path)
    errors = validate_public_dataset(dataframe)
    if errors:
        raise ValueError("\n".join(errors))
