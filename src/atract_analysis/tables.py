from __future__ import annotations

from pathlib import Path

import pandas as pd


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
    "jnet_group": "JNET",
    "conecct_group": "CONECCT",
    "fibrosis": "Fibrosis",
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


def _format_count_pct(series: pd.Series, value) -> str:
    count = int(series.eq(value).sum())
    total = int(series.notna().sum())
    pct = 100 * count / total if total else 0
    return f"{count} ({pct:.1f}%)"


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

    categorical_levels = {
        "sex": ["male", "female"],
        "operator_id_public": sorted(dataframe["operator_id_public"].dropna().unique()),
        "location_group": sorted(dataframe["location_group"].dropna().unique()),
        "lesion_type": sorted(dataframe["lesion_type"].dropna().unique()),
        "jnet_group": sorted(dataframe["jnet_group"].dropna().unique()),
        "conecct_group": sorted(dataframe["conecct_group"].dropna().unique()),
        "fibrosis": ["F0", "F1", "F2"],
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
            "large_lesion_>=50mm": "Speed in lesions ≥50 mm",
            "interaction_large_lesion": "Interaction with lesions ≥50 mm",
        }
    )

    speed_robustness_row = speed_robustness["effects"].copy()
    speed_robustness_row["section"] = "Speed"
    speed_robustness_row["method_label"] = "OW + DR"
    speed_robustness_row["outcome_label"] = "Overall speed"

    binary_effects = primary_binary["effects"].copy()
    binary_effects["section"] = "Secondary outcomes"
    binary_effects["method_label"] = "OW + DR"
    binary_effects["outcome_label"] = binary_effects["outcome"].map(OUTCOME_LABELS)

    combined = pd.concat([speed_effects, speed_robustness_row, binary_effects], ignore_index=True)
    return combined[
        [
            "section",
            "method_label",
            "outcome_label",
            "n",
            "estimate",
            "ci_lower",
            "ci_upper",
            "p_value",
            "scale",
        ]
    ]


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
