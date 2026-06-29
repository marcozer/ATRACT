from __future__ import annotations

from pathlib import Path

import pandas as pd

from .cohort import complete_case
from .config import (
    BINARY_OUTCOMES,
    CAUSAL_SPECS,
    PRIMARY_CAUSAL_SPEC,
    PRIMARY_OPERATOR_YEAR_MIN_PER_ARM,
    PRIMARY_SUPPORT_SCOPE,
    PUBLIC_COLUMNS,
)


VARIABLE_LABELS = {
    "age_years_topcoded": "Age, years",
    "major_diameter_mm": "Major diameter, mm",
    "surface_mm2": "Specimen area, mm²",
    "procedure_duration_min": "Procedure duration, min",
    "speed_mm2_min": "Dissection speed, mm²/min",
    "sex": "Sex",
    "operator_id_public": "Operator",
    "location_group": "Location",
    "lesion_type": "Lesion type",
    "lesion_morphology": "Lesion morphology",
    "jnet_group": "JNET",
    "conecct_group": "CONECCT",
    "fibrosis": "Fibrosis",
    "mici_history": "IBD history",
    "recurrence_history": "Recurrence or scar history",
}


LEVEL_LABELS = {
    "male": "Male",
    "female": "Female",
    "rectum": "Rectum",
    "sigmoid_left_colon": "Sigmoid/left colon",
    "flexures": "Flexures",
    "transverse_colon": "Transverse colon",
    "right_colon_cecum": "Right colon/cecum",
    "ileocecal_valve": "Ileocecal valve",
    "lst_granular": "LST granular",
    "lst_granular_no_macronodule": "LST granular without macronodule",
    "lst_granular_macronodule": "LST granular with macronodule",
    "lst_nongranular": "LST nongranular",
    "protruding_or_other": "Protruding or other",
    "jnet_i": "JNET I",
    "jnet_iia": "JNET IIA",
    "jnet_iib": "JNET IIB",
    "jnet_iii": "JNET III",
    "is": "Is",
    "iia": "IIA",
    "iic": "IIC",
    "iic_plus": "IIC+",
    "F0": "F0",
    "F1": "F1",
    "F2": "F2",
    0: "No",
    1: "Yes",
}


OUTCOME_LABELS = {
    "r0": "R0 resection",
    "perforation": "Intraprocedural perforation",
    "delayed_bleeding": "Delayed bleeding",
}


def _format_mean_sd(series: pd.Series) -> str:
    return f"{series.mean():.2f} ({series.std(ddof=1):.2f})"


def _format_median_iqr(series: pd.Series) -> str:
    return f"{series.median():.2f} ({series.quantile(0.25):.2f}-{series.quantile(0.75):.2f})"


def _format_count_pct(series: pd.Series, value) -> str:
    count = int(series.eq(value).sum())
    total = int(series.notna().sum())
    pct = 100 * count / total if total else 0
    return f"{count} ({pct:.1f}%)"


def _format_missing(series: pd.Series) -> str:
    missing = int(series.isna().sum())
    pct = 100 * missing / len(series) if len(series) else 0
    return f"{missing} ({pct:.1f}%)"


def _summarize_speed(series: pd.Series, prefix: str) -> dict[str, float | int]:
    observed = series.dropna()
    if observed.empty:
        return {
            f"{prefix}_n": 0,
            f"{prefix}_mean": pd.NA,
            f"{prefix}_sd": pd.NA,
            f"{prefix}_median": pd.NA,
            f"{prefix}_iqr_q1": pd.NA,
            f"{prefix}_iqr_q3": pd.NA,
        }
    return {
        f"{prefix}_n": int(observed.size),
        f"{prefix}_mean": float(observed.mean()),
        f"{prefix}_sd": float(observed.std(ddof=1)) if observed.size > 1 else pd.NA,
        f"{prefix}_median": float(observed.median()),
        f"{prefix}_iqr_q1": float(observed.quantile(0.25)),
        f"{prefix}_iqr_q3": float(observed.quantile(0.75)),
    }


def build_table_one(dataframe: pd.DataFrame) -> pd.DataFrame:
    groups = {
        "overall": dataframe,
        "non_atract": dataframe.loc[dataframe["atract"].eq(0)],
        "atract": dataframe.loc[dataframe["atract"].eq(1)],
    }
    rows: list[dict[str, str]] = [
        {
            "variable": "n",
            **{group_name: str(len(group_df)) for group_name, group_df in groups.items()},
        }
    ]

    continuous = ["age_years_topcoded", "major_diameter_mm", "surface_mm2", "procedure_duration_min", "speed_mm2_min"]
    for column in continuous:
        rows.append(
            {
                "variable": f"{VARIABLE_LABELS[column]}, mean (SD)",
                **{group_name: _format_mean_sd(group_df[column].dropna()) for group_name, group_df in groups.items()},
            }
        )
        rows.append(
            {
                "variable": f"{VARIABLE_LABELS[column]}, median (IQR)",
                **{group_name: _format_median_iqr(group_df[column].dropna()) for group_name, group_df in groups.items()},
            }
        )

    categorical_levels = {
        "sex": ["male", "female"],
        "operator_id_public": sorted(dataframe["operator_id_public"].dropna().unique()),
        "location_group": sorted(dataframe["location_group"].dropna().unique()),
        "lesion_morphology": sorted(dataframe["lesion_morphology"].dropna().unique()),
        "jnet_group": sorted(dataframe["jnet_group"].dropna().unique()),
        "conecct_group": sorted(dataframe["conecct_group"].dropna().unique()),
        "fibrosis": ["F0", "F1", "F2"],
        "mici_history": [0, 1],
        "recurrence_history": [0, 1],
    }
    for column, levels in categorical_levels.items():
        if dataframe[column].notna().sum() == 0:
            continue
        for level in levels:
            rows.append(
                {
                    "variable": f"{VARIABLE_LABELS[column]}: {LEVEL_LABELS.get(level, level)}",
                    **{group_name: _format_count_pct(group_df[column], level) for group_name, group_df in groups.items()},
                }
            )
    return pd.DataFrame(rows)


def build_matched_characteristics_table(matched: pd.DataFrame) -> pd.DataFrame:
    table = build_table_one(matched)
    table.insert(0, "population", "primary_matched_speed_cohort")
    return table


def build_unmatched_profile_table(scored: pd.DataFrame, matched: pd.DataFrame) -> pd.DataFrame:
    matched_ids = set(matched["_analysis_row_id"].dropna().astype(int))
    unmatched = scored.loc[~scored["_analysis_row_id"].astype(int).isin(matched_ids)].copy()
    table = build_table_one(unmatched)
    table.insert(0, "population", "speed_support_cohort_not_matched")
    return table


def build_operator_year_distribution_table(dataframe: pd.DataFrame) -> pd.DataFrame:
    rows = []
    grouped = dataframe.groupby(["study_year_index", "operator_id_public"], dropna=False, sort=True)
    for (year, operator_id), group in grouped:
        control = group.loc[group["atract"].eq(0)]
        treated = group.loc[group["atract"].eq(1)]
        row = {
            "study_year_index": int(year) if pd.notna(year) else pd.NA,
            "operator_id_public": operator_id,
            "total_n": int(len(group)),
            "conventional_n": int(len(control)),
            "atract_n": int(len(treated)),
            "atract_pct": float(100 * len(treated) / len(group)) if len(group) else 0.0,
        }
        row.update(_summarize_speed(group["speed_mm2_min"], "speed_overall"))
        row.update(_summarize_speed(control["speed_mm2_min"], "speed_conventional"))
        row.update(_summarize_speed(treated["speed_mm2_min"], "speed_atract"))
        rows.append(row)
    return pd.DataFrame(rows)


def build_operator_adoption_table(dataframe: pd.DataFrame) -> pd.DataFrame:
    rows = []
    first_atract_year = dataframe.loc[dataframe["atract"].eq(1)].groupby("operator_id_public")["study_year_index"].min()
    for operator_id, group in dataframe.groupby("operator_id_public", sort=True):
        first_year = first_atract_year.get(operator_id, pd.NA)
        if pd.isna(first_year):
            before = group.copy()
            after = group.iloc[0:0].copy()
        else:
            before = group.loc[group["study_year_index"].lt(first_year)]
            after = group.loc[group["study_year_index"].ge(first_year)]
        rows.append(
            {
                "operator_id_public": operator_id,
                "first_atract_study_year_index": int(first_year) if pd.notna(first_year) else pd.NA,
                "total_n": int(len(group)),
                "conventional_before_adoption_n": int(before["atract"].eq(0).sum()),
                "atract_before_adoption_n": int(before["atract"].eq(1).sum()),
                "conventional_from_adoption_n": int(after["atract"].eq(0).sum()),
                "atract_from_adoption_n": int(after["atract"].eq(1).sum()),
                "total_from_adoption_n": int(len(after)),
                "operator_category_note": (
                    "Grouped heterogeneous operator stratum"
                    if operator_id == "operator_other"
                    else "Individual pseudonymized operator stratum"
                ),
            }
        )
    return pd.DataFrame(rows)


def build_main_results_table(
    primary_speed: dict[str, object],
    speed_robustness: dict[str, object],
    primary_binary: dict[str, object],
) -> pd.DataFrame:
    speed_effects = primary_speed["effects"].copy()
    speed_effects["section"] = "Speed"
    speed_effects["method_label"] = "PS-NN"
    speed_effects["outcome_label"] = speed_effects["analysis"].map(
        {
            "overall": "Overall speed",
            "small_lesion_<50mm": "Speed in lesions <50 mm",
            "large_lesion_>=50mm": "Speed in lesions ≥50 mm",
            "interaction_large_lesion": "Interaction with lesions ≥50 mm",
        }
    )

    speed_robustness_row = speed_robustness["effects"].loc[
        speed_robustness["effects"]["analysis"].eq("overall")
    ].copy()
    speed_robustness_row["section"] = "Speed"
    speed_robustness_row["method_label"] = "OW + DR"
    speed_robustness_row["outcome_label"] = "Overall speed"

    binary_effects = primary_binary["effects"].copy()
    binary_effects["section"] = "Secondary outcomes"
    binary_effects["method_label"] = "OW + DR"
    binary_effects["outcome_label"] = binary_effects["outcome"].map(OUTCOME_LABELS)

    combined = pd.concat([speed_effects, speed_robustness_row, binary_effects], ignore_index=True)
    output_columns = [
        "section",
        "method_label",
        "outcome_label",
        "n",
        "estimate",
        "ci_lower",
        "ci_upper",
        "p_value",
        "scale",
        "treated_n",
        "treated_events",
        "treated_risk_observed",
        "control_n",
        "control_events",
        "control_risk_observed",
        "treated_risk",
        "control_risk",
        "risk_difference",
        "risk_difference_ci_lower",
        "risk_difference_ci_upper",
    ]
    return combined.reindex(columns=output_columns)


def _support_restrict(dataframe: pd.DataFrame) -> pd.DataFrame:
    if PRIMARY_SUPPORT_SCOPE != "operator_only":
        raise ValueError(f"Unsupported support scope for reporting: {PRIMARY_SUPPORT_SCOPE}")
    counts = dataframe.groupby(["operator_id_public", "atract"]).size().unstack(fill_value=0)
    eligible = counts.index[counts.min(axis=1) >= PRIMARY_OPERATOR_YEAR_MIN_PER_ARM]
    return dataframe.loc[dataframe["operator_id_public"].isin(eligible)].copy()


def build_missingness_table(dataframe: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for column in PUBLIC_COLUMNS:
        missing = int(dataframe[column].isna().sum())
        rows.append(
            {
                "variable": column,
                "missing_n": missing,
                "missing_pct": 100 * missing / len(dataframe),
            }
        )
    return pd.DataFrame(rows)


def build_missingness_by_group_table(dataframe: pd.DataFrame) -> pd.DataFrame:
    groups = {
        "overall": dataframe,
        "conventional": dataframe.loc[dataframe["atract"].eq(0)],
        "atract": dataframe.loc[dataframe["atract"].eq(1)],
    }
    rows = []
    for column in PUBLIC_COLUMNS:
        rows.append(
            {
                "variable": column,
                **{group_name: _format_missing(group_df[column]) for group_name, group_df in groups.items()},
            }
        )
    return pd.DataFrame(rows)


def build_population_accounting_table(dataframe: pd.DataFrame) -> pd.DataFrame:
    spec = CAUSAL_SPECS[PRIMARY_CAUSAL_SPEC]
    definitions = [
        ("speed", "speed_mm2_min", list(spec["speed_required_columns"])),
        *[(outcome, outcome, list(spec["binary_required_columns"])) for outcome in BINARY_OUTCOMES],
    ]
    rows = []
    for analysis, outcome, required_columns in definitions:
        complete = complete_case(dataframe, required_columns, outcome=outcome)
        supported = _support_restrict(complete)
        rows.append(
            {
                "analysis": analysis,
                "public_cohort_n": int(len(dataframe)),
                "complete_case_n": int(len(complete)),
                "support_restricted_n": int(len(supported)),
                "excluded_for_missingness_n": int(len(dataframe) - len(complete)),
                "excluded_for_support_n": int(len(complete) - len(supported)),
            }
        )
    return pd.DataFrame(rows)


def write_tables(
    tables_dir: Path,
    public_dataframe: pd.DataFrame,
    primary_speed: dict[str, object],
    speed_robustness: dict[str, object],
    primary_binary: dict[str, object],
    binary_comparison: dict[str, object],
) -> None:
    tables_dir.mkdir(parents=True, exist_ok=True)
    stale_files = [
        "table_s2_speed_effects.csv",
        "table_s3_binary_effects.csv",
        "table_s4_speed_sensitivity.csv",
        "table_s5_binary_sensitivity.csv",
        "table_s4_speed_iptw_frontier.csv",
        "table_s5_speed_nn_frontier.csv",
        "table_s4_missingness.csv",
        "table_s5_population_accounting.csv",
        "table_s6_speed_temporal_sensitivity.csv",
        "table_s7_speed_bootstrap.csv",
        "table_s8_binary_balance.csv",
        "table_s9_speed_continuous_size.csv",
        "table_s10_operator_year_distribution.csv",
        "table_s11_operator_adoption.csv",
        "table_s12_matched_characteristics.csv",
        "table_s13_unmatched_profile.csv",
        "table_s14_missingness_by_group.csv",
        "table_s15_speed_temporal_relaxation_grid.csv",
    ]
    for filename in stale_files:
        stale_path = tables_dir / filename
        if stale_path.exists():
            stale_path.unlink()

    build_table_one(public_dataframe).to_csv(tables_dir / "table_1_baseline.csv", index=False)
    build_main_results_table(primary_speed, speed_robustness, primary_binary).to_csv(
        tables_dir / "table_2_main_results.csv",
        index=False,
    )
    primary_speed["balance"].to_csv(tables_dir / "table_s1_balance.csv", index=False)
    speed_robustness["summary"].to_csv(tables_dir / "table_s2_speed_robustness.csv", index=False)
    binary_comparison["summary"].to_csv(tables_dir / "table_s3_secondary_comparison.csv", index=False)
    build_missingness_table(public_dataframe).to_csv(tables_dir / "table_s4_missingness.csv", index=False)
    build_population_accounting_table(public_dataframe).to_csv(tables_dir / "table_s5_population_accounting.csv", index=False)
    primary_speed["temporal_sensitivity"]["summary"].to_csv(
        tables_dir / "table_s6_speed_temporal_sensitivity.csv",
        index=False,
    )
    primary_speed["bootstrap"].to_csv(tables_dir / "table_s7_speed_bootstrap.csv", index=False)
    primary_binary["balance"].to_csv(tables_dir / "table_s8_binary_balance.csv", index=False)
    primary_speed["continuous_size"].to_csv(tables_dir / "table_s9_speed_continuous_size.csv", index=False)
    build_operator_year_distribution_table(public_dataframe).to_csv(
        tables_dir / "table_s10_operator_year_distribution.csv",
        index=False,
    )
    build_operator_adoption_table(public_dataframe).to_csv(
        tables_dir / "table_s11_operator_adoption.csv",
        index=False,
    )
    build_matched_characteristics_table(primary_speed["matched_frame"]).to_csv(
        tables_dir / "table_s12_matched_characteristics.csv",
        index=False,
    )
    build_unmatched_profile_table(primary_speed["scored_frame"], primary_speed["matched_frame"]).to_csv(
        tables_dir / "table_s13_unmatched_profile.csv",
        index=False,
    )
    build_missingness_by_group_table(public_dataframe).to_csv(
        tables_dir / "table_s14_missingness_by_group.csv",
        index=False,
    )
    primary_speed["temporal_relaxation_grid"].to_csv(
        tables_dir / "table_s15_speed_temporal_relaxation_grid.csv",
        index=False,
    )
