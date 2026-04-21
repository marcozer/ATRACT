from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.optimize import linear_sum_assignment

from .cohort import add_large_lesion_flag, complete_case
from .config import (
    BINARY_OUTCOMES,
    CAUSAL_SPECS,
    LARGE_LESION_CUTOFF_MM,
    NN_MATCH_RATIO,
    NN_MATCH_SCOPE,
    PRIMARY_CAUSAL_SPEC,
    PRIMARY_MATCH_CALIPER_GRID,
    PRIMARY_MATCH_MAX_OPERATOR_SHARE,
    PRIMARY_MATCH_MAX_SMD,
    PRIMARY_OPERATOR_YEAR_MIN_PER_ARM,
    PRIMARY_PS_STRUCTURE,
    PRIMARY_SUPPORT_SCOPE,
)


def _weighted_mean(values: pd.Series, weights: pd.Series) -> float:
    return float(np.average(values, weights=weights))


def _weighted_var(values: pd.Series, weights: pd.Series) -> float:
    mean = _weighted_mean(values, weights)
    return float(np.average((values - mean) ** 2, weights=weights))


def _smd_for_vector(values: pd.Series, treatment: pd.Series, weights: pd.Series | None = None) -> float:
    treated_mask = treatment.eq(1)
    if weights is None:
        treated_values = values.loc[treated_mask]
        control_values = values.loc[~treated_mask]
        treated_mean = treated_values.mean()
        control_mean = control_values.mean()
        treated_var = treated_values.var(ddof=0)
        control_var = control_values.var(ddof=0)
    else:
        treated_values = values.loc[treated_mask]
        control_values = values.loc[~treated_mask]
        treated_weights = weights.loc[treated_mask]
        control_weights = weights.loc[~treated_mask]
        treated_mean = _weighted_mean(treated_values, treated_weights)
        control_mean = _weighted_mean(control_values, control_weights)
        treated_var = _weighted_var(treated_values, treated_weights)
        control_var = _weighted_var(control_values, control_weights)
    pooled_sd = np.sqrt((treated_var + control_var) / 2)
    if pooled_sd == 0:
        return 0.0
    return float((treated_mean - control_mean) / pooled_sd)


def build_balance_table(dataframe: pd.DataFrame, covariates: list[str], weights_column: str | None = None) -> pd.DataFrame:
    treatment = dataframe["atract"].astype(int)
    weights = dataframe[weights_column] if weights_column else None
    design = pd.get_dummies(dataframe[covariates], drop_first=False, dtype=float)
    rows: list[dict[str, float | str]] = []
    for column in design.columns:
        rows.append(
            {
                "term": column,
                "raw_smd": _smd_for_vector(design[column], treatment),
                "adjusted_smd": _smd_for_vector(design[column], treatment, weights) if weights is not None else np.nan,
            }
        )
    return pd.DataFrame(rows).sort_values("term").reset_index(drop=True)


def build_matching_balance_table(reference: pd.DataFrame, matched: pd.DataFrame, covariates: list[str]) -> pd.DataFrame:
    reference_design = pd.get_dummies(reference[covariates], drop_first=False, dtype=float)
    matched_design = pd.get_dummies(matched[covariates], drop_first=False, dtype=float)
    terms = sorted(set(reference_design.columns) | set(matched_design.columns))
    reference_design = reference_design.reindex(columns=terms, fill_value=0.0)
    matched_design = matched_design.reindex(columns=terms, fill_value=0.0)
    reference_treatment = reference["atract"].astype(int)
    matched_treatment = matched["atract"].astype(int)

    rows: list[dict[str, float | str]] = []
    for term in terms:
        rows.append(
            {
                "term": term,
                "raw_smd": _smd_for_vector(reference_design[term], reference_treatment),
                "adjusted_smd": _smd_for_vector(matched_design[term], matched_treatment),
            }
        )
    return pd.DataFrame(rows).sort_values("term").reset_index(drop=True)


def _max_operator_treated_share(dataframe: pd.DataFrame) -> float:
    treated = dataframe.loc[dataframe["atract"].eq(1), "operator_id_public"]
    if treated.empty:
        return 1.0
    return float(treated.value_counts(normalize=True).iloc[0])


def _restrict_support(dataframe: pd.DataFrame, scope: str, min_per_arm: int) -> pd.DataFrame:
    if scope == "none":
        return dataframe.copy()
    if scope == "operator_only":
        counts = dataframe.groupby(["operator_id_public", "atract"]).size().unstack(fill_value=0)
        eligible = counts.index[counts.min(axis=1) >= min_per_arm]
        return dataframe.loc[dataframe["operator_id_public"].isin(eligible)].copy()
    if scope == "operator_year":
        counts = dataframe.groupby(["operator_id_public", "study_year_index", "atract"]).size().unstack(fill_value=0)
        eligible = counts.index[counts.min(axis=1) >= min_per_arm]
        eligible_index = pd.MultiIndex.from_frame(dataframe[["operator_id_public", "study_year_index"]])
        return dataframe.loc[eligible_index.isin(eligible)].copy()
    raise ValueError(f"Unsupported support scope: {scope}")


def _analysis_cohort(dataframe: pd.DataFrame, *, required_columns: list[str], outcome: str) -> pd.DataFrame:
    complete = complete_case(dataframe, required_columns, outcome=outcome)
    return _restrict_support(
        complete,
        scope=PRIMARY_SUPPORT_SCOPE,
        min_per_arm=PRIMARY_OPERATOR_YEAR_MIN_PER_ARM,
    )


def _prepare_propensity_formula(dataframe: pd.DataFrame, rhs: str, ps_structure: str) -> tuple[pd.DataFrame, str]:
    prepared = dataframe.copy()
    if ps_structure == "operator_plus_year":
        structure = "C(operator_id_public) + C(study_year_index)"
    elif ps_structure == "operator_only_plus_year_numeric":
        structure = "C(operator_id_public) + study_year_index"
    elif ps_structure == "operator_year_strata":
        prepared["ps_stratum"] = (
            prepared["operator_id_public"].astype("string")
            + "__y"
            + prepared["study_year_index"].astype("Int64").astype("string")
        )
        structure = "C(ps_stratum)"
    else:
        raise ValueError(f"Unsupported PS structure: {ps_structure}")
    return prepared, f"atract ~ {structure} + {rhs}"


def _compute_weights(dataframe: pd.DataFrame, method: str) -> np.ndarray:
    scores = dataframe["propensity_score"].to_numpy()
    treated = dataframe["atract"].to_numpy()
    if method == "overlap":
        return np.where(treated == 1, 1 - scores, scores)
    raise ValueError(f"Unsupported weighting method: {method}")


def fit_propensity_model(
    dataframe: pd.DataFrame,
    *,
    ps_rhs: str,
    ps_structure: str,
    weighting_method: str,
) -> tuple[pd.DataFrame, object]:
    modeled, formula = _prepare_propensity_formula(dataframe, rhs=ps_rhs, ps_structure=ps_structure)
    try:
        propensity_model = smf.logit(formula, data=modeled).fit(disp=False, maxiter=500, method="bfgs")
    except Exception:
        propensity_model = smf.logit(formula, data=modeled).fit_regularized(alpha=1e-6, maxiter=1000, disp=False)
    scores = np.asarray(propensity_model.predict(modeled)).clip(1e-4, 1 - 1e-4)
    modeled["propensity_score"] = scores
    modeled["analysis_weight"] = _compute_weights(modeled, weighting_method)
    return modeled, propensity_model


def _extract_linear_effect(result, term: str, scale: str) -> dict[str, float | str]:
    conf_int = result.conf_int().loc[term]
    return {
        "term": term,
        "scale": scale,
        "estimate": float(result.params[term]),
        "ci_lower": float(conf_int.iloc[0]),
        "ci_upper": float(conf_int.iloc[1]),
        "p_value": float(result.pvalues[term]),
    }


def _extract_exponentiated_effect(result, term: str, scale: str) -> dict[str, float | str]:
    effect = _extract_linear_effect(result, term, scale)
    return {
        **effect,
        "estimate": float(np.exp(effect["estimate"])),
        "ci_lower": float(np.exp(effect["ci_lower"])),
        "ci_upper": float(np.exp(effect["ci_upper"])),
    }


def _nearest_neighbor_match(
    dataframe: pd.DataFrame,
    *,
    ratio: int,
    caliper_multiplier: float,
    match_scope: str,
) -> pd.DataFrame:
    if match_scope != "operator_only":
        raise ValueError(f"Unsupported match scope: {match_scope}")

    matched_rows: list[pd.DataFrame] = []
    logit_ps = np.log(dataframe["propensity_score"] / (1 - dataframe["propensity_score"]))
    caliper = float(caliper_multiplier * np.std(logit_ps, ddof=1))

    treated = dataframe.loc[dataframe["atract"].eq(1)].copy()
    controls = dataframe.loc[dataframe["atract"].eq(0)].copy()
    treated_groups = treated.groupby("operator_id_public").groups
    control_groups = controls.groupby("operator_id_public").groups

    for operator_id in sorted(set(treated_groups) & set(control_groups)):
        treated_group = dataframe.loc[list(treated_groups[operator_id])].copy()
        control_group = dataframe.loc[list(control_groups[operator_id])].copy()
        repeated_treated = np.repeat(np.arange(len(treated_group)), ratio)
        cost = np.abs(
            logit_ps.loc[treated_group.index].to_numpy()[repeated_treated][:, None]
            - logit_ps.loc[control_group.index].to_numpy()[None, :]
        )
        cost[cost > caliper] = 1e6
        row_ind, col_ind = linear_sum_assignment(cost)
        keep = cost[row_ind, col_ind] < 1e6
        if not np.any(keep):
            continue
        matched_treated = repeated_treated[row_ind[keep]]
        matched_controls = col_ind[keep]
        complete_treated = set(
            pd.Series(matched_treated).value_counts().loc[lambda counts: counts.eq(ratio)].index
        )
        for treated_position in sorted(complete_treated):
            match_group = f"{operator_id}__{treated_group.index[treated_position]}"
            treated_row = treated_group.iloc[[treated_position]].copy()
            treated_row["_match_group"] = match_group
            matched_rows.append(treated_row)
            control_positions = np.where(matched_treated == treated_position)[0]
            control_rows = control_group.iloc[matched_controls[control_positions]].copy()
            control_rows["_match_group"] = match_group
            matched_rows.append(control_rows)

    if not matched_rows:
        return pd.DataFrame(columns=dataframe.columns)
    return pd.concat(matched_rows, ignore_index=True)


def _fit_primary_speed_models(matched: pd.DataFrame, augment_rhs: str) -> tuple[pd.DataFrame, dict[str, float | int | str]]:
    matched = add_large_lesion_flag(matched)
    cluster_codes = matched["_match_group"].astype("category").cat.codes
    overall_model = smf.ols(
        f"speed_mm2_min ~ atract + {augment_rhs}",
        data=matched,
    ).fit(cov_type="cluster", cov_kwds={"groups": cluster_codes})

    large_lesions = matched.loc[matched["large_lesion"].eq(1)].copy()
    large_cluster_codes = large_lesions["_match_group"].astype("category").cat.codes
    large_model = smf.ols(
        f"speed_mm2_min ~ atract + {augment_rhs}",
        data=large_lesions,
    ).fit(cov_type="cluster", cov_kwds={"groups": large_cluster_codes})

    interaction_model = smf.ols(
        f"speed_mm2_min ~ atract * large_lesion + {augment_rhs}",
        data=matched,
    ).fit(cov_type="cluster", cov_kwds={"groups": cluster_codes})

    effects = pd.DataFrame(
        [
            {
                "method": "ps_nn",
                "analysis": "overall",
                "outcome": "speed_mm2_min",
                "n": int(len(matched)),
                **_extract_linear_effect(overall_model, "atract", "mean_difference"),
            },
            {
                "method": "ps_nn",
                "analysis": f"large_lesion_>={LARGE_LESION_CUTOFF_MM}mm",
                "outcome": "speed_mm2_min",
                "n": int(len(large_lesions)),
                **_extract_linear_effect(large_model, "atract", "mean_difference"),
            },
            {
                "method": "ps_nn",
                "analysis": "interaction_large_lesion",
                "outcome": "speed_mm2_min",
                "n": int(len(matched)),
                **_extract_linear_effect(interaction_model, "atract:large_lesion", "interaction_difference"),
            },
        ]
    )
    metadata = {
        "analysis": "primary_speed_ps_nn",
        "n_complete_case": int(len(matched)),
        "n_large_lesion": int(len(large_lesions)),
        "n_match_groups": int(matched["_match_group"].nunique()),
        "top_operator_treated_share": _max_operator_treated_share(matched),
    }
    return effects, metadata


def run_primary_speed_analysis(dataframe: pd.DataFrame) -> dict[str, object]:
    spec = CAUSAL_SPECS[PRIMARY_CAUSAL_SPEC]
    analysis_cohort = _analysis_cohort(
        dataframe,
        required_columns=list(spec["speed_required_columns"]),
        outcome="speed_mm2_min",
    )
    scored, propensity_model = fit_propensity_model(
        analysis_cohort,
        ps_rhs=str(spec["ps_rhs"]),
        ps_structure=PRIMARY_PS_STRUCTURE,
        weighting_method="overlap",
    )

    candidate_rows: list[dict[str, float | int | bool]] = []
    payloads: dict[float, tuple[pd.DataFrame, pd.DataFrame]] = {}
    for caliper_multiplier in PRIMARY_MATCH_CALIPER_GRID:
        matched = _nearest_neighbor_match(
            scored,
            ratio=NN_MATCH_RATIO,
            caliper_multiplier=float(caliper_multiplier),
            match_scope=NN_MATCH_SCOPE,
        )
        if matched.empty:
            continue
        balance = build_matching_balance_table(scored, matched, list(spec["balance_covariates"]))
        max_abs_smd = float(balance["adjusted_smd"].abs().max())
        top_operator_share = _max_operator_treated_share(matched)
        payloads[float(caliper_multiplier)] = (matched, balance)
        candidate_rows.append(
            {
                "caliper_multiplier": float(caliper_multiplier),
                "n": int(len(matched)),
                "n_treated": int(matched["atract"].eq(1).sum()),
                "n_control": int(matched["atract"].eq(0).sum()),
                "n_match_groups": int(matched["_match_group"].nunique()),
                "max_abs_smd": max_abs_smd,
                "top_operator_treated_share": top_operator_share,
                "valid": max_abs_smd <= PRIMARY_MATCH_MAX_SMD and top_operator_share <= PRIMARY_MATCH_MAX_OPERATOR_SHARE,
            }
        )

    matching_grid = pd.DataFrame(candidate_rows).sort_values("caliper_multiplier").reset_index(drop=True)
    if matching_grid.empty:
        raise RuntimeError("No valid PS matching candidates were generated.")

    valid_candidates = matching_grid.loc[matching_grid["valid"]].copy()
    if valid_candidates.empty:
        selected = matching_grid.sort_values(
            by=["max_abs_smd", "top_operator_treated_share", "n", "caliper_multiplier"],
            ascending=[True, True, False, True],
        ).iloc[0]
    else:
        selected = valid_candidates.sort_values(
            by=["n", "max_abs_smd", "top_operator_treated_share", "caliper_multiplier"],
            ascending=[False, True, True, True],
        ).iloc[0]

    selected_caliper = float(selected["caliper_multiplier"])
    matched_frame, balance_table = payloads[selected_caliper]
    effects, metadata = _fit_primary_speed_models(matched_frame, str(spec["speed_aug_rhs"]))
    metadata.update(
        {
            "spec_key": PRIMARY_CAUSAL_SPEC,
            "support_scope": PRIMARY_SUPPORT_SCOPE,
            "min_per_arm": PRIMARY_OPERATOR_YEAR_MIN_PER_ARM,
            "ps_structure": PRIMARY_PS_STRUCTURE,
            "match_scope": NN_MATCH_SCOPE,
            "ratio": NN_MATCH_RATIO,
            "caliper_multiplier": selected_caliper,
            "max_abs_smd": float(balance_table["adjusted_smd"].abs().max()),
            "n_supported_cohort": int(len(scored)),
        }
    )
    return {
        "scored_frame": scored,
        "matched_frame": matched_frame,
        "propensity_model": propensity_model,
        "balance": balance_table,
        "effects": effects,
        "matching_grid": matching_grid,
        "metadata": metadata,
    }


def _fit_weighted_speed_effects(weighted: pd.DataFrame, augment_rhs: str) -> tuple[pd.DataFrame, dict[str, float | int | str]]:
    model = smf.wls(
        f"speed_mm2_min ~ atract + {augment_rhs}",
        data=weighted,
        weights=weighted["analysis_weight"],
    ).fit(cov_type="HC3")
    effects = pd.DataFrame(
        [
            {
                "method": "ow_dr",
                "analysis": "overall",
                "outcome": "speed_mm2_min",
                "n": int(len(weighted)),
                **_extract_linear_effect(model, "atract", "mean_difference"),
            }
        ]
    )
    metadata = {
        "analysis": "speed_overlap_weighting_robustness",
        "n_complete_case": int(len(weighted)),
        "weight_min": float(weighted["analysis_weight"].min()),
        "weight_max": float(weighted["analysis_weight"].max()),
    }
    return effects, metadata


def run_speed_robustness_analysis(dataframe: pd.DataFrame, primary_speed: dict[str, object]) -> dict[str, object]:
    spec = CAUSAL_SPECS[PRIMARY_CAUSAL_SPEC]
    analysis_cohort = _analysis_cohort(
        dataframe,
        required_columns=list(spec["speed_required_columns"]),
        outcome="speed_mm2_min",
    )
    weighted, propensity_model = fit_propensity_model(
        analysis_cohort,
        ps_rhs=str(spec["ps_rhs"]),
        ps_structure=PRIMARY_PS_STRUCTURE,
        weighting_method="overlap",
    )
    balance = build_balance_table(weighted, list(spec["balance_covariates"]), weights_column="analysis_weight")
    effects, metadata = _fit_weighted_speed_effects(weighted, str(spec["speed_aug_rhs"]))
    metadata.update(
        {
            "spec_key": PRIMARY_CAUSAL_SPEC,
            "support_scope": PRIMARY_SUPPORT_SCOPE,
            "min_per_arm": PRIMARY_OPERATOR_YEAR_MIN_PER_ARM,
            "ps_structure": PRIMARY_PS_STRUCTURE,
            "max_abs_weighted_smd": float(balance["adjusted_smd"].abs().max()),
        }
    )

    primary_overall = primary_speed["effects"].loc[primary_speed["effects"]["analysis"].eq("overall")].iloc[0].to_dict()
    primary_overall.update(
        {
            "method_label": "PS-NN",
            "max_abs_smd": float(primary_speed["metadata"]["max_abs_smd"]),
            "weight_max": np.nan,
            "n_match_groups": int(primary_speed["metadata"]["n_match_groups"]),
        }
    )
    robustness_overall = effects.iloc[0].to_dict()
    robustness_overall.update(
        {
            "method_label": "OW + DR",
            "max_abs_smd": float(balance["adjusted_smd"].abs().max()),
            "weight_max": float(metadata["weight_max"]),
            "n_match_groups": np.nan,
        }
    )

    return {
        "weighted_frame": weighted,
        "propensity_model": propensity_model,
        "balance": balance,
        "effects": effects,
        "metadata": metadata,
        "summary": pd.DataFrame([primary_overall, robustness_overall]),
    }


def _fit_weighted_binary_effect(
    weighted: pd.DataFrame,
    *,
    outcome: str,
    augment_rhs: str,
) -> dict[str, float | str]:
    model = smf.glm(
        f"{outcome} ~ atract + {augment_rhs}",
        data=weighted,
        family=sm.families.Poisson(),
        freq_weights=weighted["analysis_weight"],
    ).fit(cov_type="HC0")
    return _extract_exponentiated_effect(model, "atract", "risk_ratio")


def run_primary_binary_analyses(dataframe: pd.DataFrame) -> dict[str, object]:
    spec = CAUSAL_SPECS[PRIMARY_CAUSAL_SPEC]
    rows: list[dict[str, float | str]] = []
    diagnostics: dict[str, dict[str, float | int | str]] = {}

    for outcome in BINARY_OUTCOMES:
        analysis_cohort = _analysis_cohort(
            dataframe,
            required_columns=list(spec["binary_required_columns"]),
            outcome=outcome,
        )
        weighted, _ = fit_propensity_model(
            analysis_cohort,
            ps_rhs=str(spec["ps_rhs"]),
            ps_structure=PRIMARY_PS_STRUCTURE,
            weighting_method="overlap",
        )
        balance = build_balance_table(weighted, list(spec["balance_covariates"]), weights_column="analysis_weight")
        rows.append(
            {
                "method": "ow_dr",
                "analysis": "overall",
                "outcome": outcome,
                "n": int(len(weighted)),
                **_fit_weighted_binary_effect(weighted, outcome=outcome, augment_rhs=str(spec["binary_aug_rhs"])),
            }
        )
        diagnostics[outcome] = {
            "n_complete_case": int(len(weighted)),
            "weight_min": float(weighted["analysis_weight"].min()),
            "weight_max": float(weighted["analysis_weight"].max()),
            "max_abs_weighted_smd": float(balance["adjusted_smd"].abs().max()),
            "spec_key": PRIMARY_CAUSAL_SPEC,
            "support_scope": PRIMARY_SUPPORT_SCOPE,
            "min_per_arm": PRIMARY_OPERATOR_YEAR_MIN_PER_ARM,
            "ps_structure": PRIMARY_PS_STRUCTURE,
        }

    return {
        "effects": pd.DataFrame(rows),
        "diagnostics": diagnostics,
    }


def run_binary_comparison_analyses(
    dataframe: pd.DataFrame,
    *,
    caliper_multiplier: float,
) -> dict[str, object]:
    spec = CAUSAL_SPECS[PRIMARY_CAUSAL_SPEC]
    rows: list[dict[str, float | str]] = []

    for outcome in BINARY_OUTCOMES:
        analysis_cohort = _analysis_cohort(
            dataframe,
            required_columns=list(spec["binary_required_columns"]),
            outcome=outcome,
        )
        weighted, _ = fit_propensity_model(
            analysis_cohort,
            ps_rhs=str(spec["ps_rhs"]),
            ps_structure=PRIMARY_PS_STRUCTURE,
            weighting_method="overlap",
        )
        ow_balance = build_balance_table(weighted, list(spec["balance_covariates"]), weights_column="analysis_weight")
        rows.append(
            {
                "method_label": "OW + DR",
                "outcome": outcome,
                "n": int(len(weighted)),
                "max_abs_smd": float(ow_balance["adjusted_smd"].abs().max()),
                **_fit_weighted_binary_effect(weighted, outcome=outcome, augment_rhs=str(spec["binary_aug_rhs"])),
            }
        )

        matched = _nearest_neighbor_match(
            weighted,
            ratio=NN_MATCH_RATIO,
            caliper_multiplier=caliper_multiplier,
            match_scope=NN_MATCH_SCOPE,
        )
        matched_balance = build_matching_balance_table(weighted, matched, list(spec["balance_covariates"]))
        cluster_codes = matched["_match_group"].astype("category").cat.codes
        matched_model = smf.glm(
            f"{outcome} ~ atract + {spec['binary_aug_rhs']}",
            data=matched,
            family=sm.families.Poisson(),
        ).fit(cov_type="cluster", cov_kwds={"groups": cluster_codes})
        rows.append(
            {
                "method_label": "PS-NN",
                "outcome": outcome,
                "n": int(len(matched)),
                "n_match_groups": int(matched["_match_group"].nunique()),
                "max_abs_smd": float(matched_balance["adjusted_smd"].abs().max()),
                **_extract_exponentiated_effect(matched_model, "atract", "risk_ratio"),
            }
        )

    return {"summary": pd.DataFrame(rows)}


def write_model_outputs(
    results_dir: Path,
    primary_speed: dict[str, object],
    speed_robustness: dict[str, object],
    primary_binary: dict[str, object],
    binary_comparison: dict[str, object],
) -> None:
    results_dir.mkdir(parents=True, exist_ok=True)
    stale_paths = [
        "binary_effects.csv",
        "binary_sensitivity_summary.csv",
        "speed_sensitivity_summary.csv",
        "speed_effects.csv",
        "speed_iptw_frontier.csv",
        "speed_nn_frontier.csv",
    ]
    for filename in stale_paths:
        stale_path = results_dir / filename
        if stale_path.exists():
            stale_path.unlink()

    pd.DataFrame(
        {
            "term": primary_speed["propensity_model"].params.index,
            "coefficient": primary_speed["propensity_model"].params.values,
        }
    ).to_csv(results_dir / "propensity_coefficients.csv", index=False)
    primary_speed["effects"].to_csv(results_dir / "primary_speed_results.csv", index=False)
    speed_robustness["summary"].to_csv(results_dir / "speed_robustness_results.csv", index=False)
    primary_binary["effects"].to_csv(results_dir / "primary_binary_results.csv", index=False)
    binary_comparison["summary"].to_csv(results_dir / "binary_comparison_results.csv", index=False)
    primary_speed["matching_grid"].to_csv(results_dir / "speed_matching_grid.csv", index=False)

    metadata = {
        "primary_speed": primary_speed["metadata"],
        "speed_robustness": speed_robustness["metadata"],
        "primary_binary": primary_binary["diagnostics"],
        "matching": {
            "ratio": NN_MATCH_RATIO,
            "match_scope": NN_MATCH_SCOPE,
            "caliper_grid": PRIMARY_MATCH_CALIPER_GRID,
            "selected_caliper_multiplier": float(primary_speed["metadata"]["caliper_multiplier"]),
        },
    }
    with (results_dir / "analysis_metadata.json").open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)
