from __future__ import annotations

import json
import platform
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from patsy import build_design_matrices
from scipy.optimize import linear_sum_assignment
from scipy.stats import norm

from .cohort import add_large_lesion_flag, complete_case
from .config import (
    BINARY_OUTCOMES,
    CAUSAL_SPECS,
    EXPECTED_PUBLIC_DATASET_SHA256,
    LARGE_LESION_CUTOFF_MM,
    NN_MATCH_RATIO,
    NN_MATCH_SCOPE,
    PRIMARY_BOOTSTRAP_ITERATIONS,
    PRIMARY_BOOTSTRAP_SEED,
    PRIMARY_CAUSAL_SPEC,
    PRIMARY_MATCH_CALIPER_GRID,
    PRIMARY_MATCH_MAX_OPERATOR_SHARE,
    PRIMARY_MATCH_MAX_SMD,
    PRIMARY_OPERATOR_YEAR_MIN_PER_ARM,
    PRIMARY_PS_STRUCTURE,
    PRIMARY_SUPPORT_SCOPE,
    PRIMARY_TEMPORAL_MAX_YEAR_GAP,
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
    supported = _restrict_support(
        complete,
        scope=PRIMARY_SUPPORT_SCOPE,
        min_per_arm=PRIMARY_OPERATOR_YEAR_MIN_PER_ARM,
    )
    supported = supported.copy()
    supported["_analysis_row_id"] = supported.index.astype(int)
    return supported


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


def _effect_from_influence(estimate: float, influence: np.ndarray, scale: str, term: str) -> dict[str, float | str]:
    standard_error = float(np.std(influence, ddof=1) / np.sqrt(len(influence)))
    z_value = estimate / standard_error if standard_error else np.nan
    p_value = float(2 * (1 - norm.cdf(abs(z_value)))) if np.isfinite(z_value) else np.nan
    return {
        "term": term,
        "scale": scale,
        "estimate": float(estimate),
        "ci_lower": float(estimate - 1.96 * standard_error),
        "ci_upper": float(estimate + 1.96 * standard_error),
        "p_value": p_value,
    }


def _extract_linear_combination(result, terms: list[str], *, label: str, scale: str) -> dict[str, float | str]:
    coefficients = result.params.reindex(terms)
    covariance = result.cov_params().reindex(index=terms, columns=terms)
    estimate = float(coefficients.sum())
    variance = float(np.ones(len(terms)) @ covariance.to_numpy() @ np.ones(len(terms)))
    standard_error = float(np.sqrt(max(variance, 0.0)))
    z_value = estimate / standard_error if standard_error else np.nan
    p_value = float(2 * (1 - norm.cdf(abs(z_value)))) if np.isfinite(z_value) else np.nan
    return {
        "term": label,
        "scale": scale,
        "estimate": estimate,
        "ci_lower": estimate - 1.96 * standard_error,
        "ci_upper": estimate + 1.96 * standard_error,
        "p_value": p_value,
    }


def _aipw_components(
    dataframe: pd.DataFrame,
    *,
    outcome: str,
    outcome_rhs: str,
    family=None,
    include_large_interaction: bool = False,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    modeled = add_large_lesion_flag(dataframe)
    treatment_term = "atract * large_lesion" if include_large_interaction else "atract"
    formula = f"{outcome} ~ {treatment_term} + {outcome_rhs}"
    if family is None:
        outcome_model = smf.ols(formula, data=modeled).fit()
    else:
        outcome_model = smf.glm(formula, data=modeled, family=family).fit()

    treated_frame = modeled.copy()
    control_frame = modeled.copy()
    treated_frame["atract"] = 1
    control_frame["atract"] = 0
    predicted_treated = np.asarray(outcome_model.predict(treated_frame))
    predicted_control = np.asarray(outcome_model.predict(control_frame))
    if family is not None:
        predicted_treated = predicted_treated.clip(1e-6, 1 - 1e-6)
        predicted_control = predicted_control.clip(1e-6, 1 - 1e-6)
    return modeled, predicted_treated, predicted_control


def _aipw_overlap_mean_difference(
    dataframe: pd.DataFrame,
    *,
    outcome: str,
    outcome_rhs: str,
    analysis: str,
    mask: pd.Series | None = None,
    include_large_interaction: bool = False,
) -> tuple[dict[str, float | str], np.ndarray]:
    modeled, predicted_treated, predicted_control = _aipw_components(
        dataframe,
        outcome=outcome,
        outcome_rhs=outcome_rhs,
        include_large_interaction=include_large_interaction,
    )
    treatment = modeled["atract"].to_numpy(dtype=int)
    observed = modeled[outcome].to_numpy(dtype=float)
    propensity = modeled["propensity_score"].to_numpy(dtype=float)
    overlap = propensity * (1 - propensity)
    subgroup = np.ones(len(modeled), dtype=float) if mask is None else mask.to_numpy(dtype=float)
    pseudo_outcome = (
        overlap * (predicted_treated - predicted_control)
        + treatment * overlap / propensity * (observed - predicted_treated)
        - (1 - treatment) * overlap / (1 - propensity) * (observed - predicted_control)
    )
    denominator_terms = subgroup * overlap
    denominator = float(np.mean(denominator_terms))
    estimate = float(np.mean(subgroup * pseudo_outcome) / denominator)
    influence = (subgroup * pseudo_outcome - estimate * denominator_terms) / denominator
    return _effect_from_influence(estimate, influence, "mean_difference", analysis), influence


def _aipw_overlap_mean_decomposition(
    dataframe: pd.DataFrame,
    *,
    outcome: str,
    outcome_rhs: str,
) -> dict[str, float | int]:
    modeled, predicted_treated, predicted_control = _aipw_components(
        dataframe,
        outcome=outcome,
        outcome_rhs=outcome_rhs,
    )
    treatment = modeled["atract"].to_numpy(dtype=int)
    observed = modeled[outcome].to_numpy(dtype=float)
    propensity = modeled["propensity_score"].to_numpy(dtype=float)
    overlap = propensity * (1 - propensity)
    denominator = float(np.mean(overlap))
    model_contrast = float(np.mean(overlap * (predicted_treated - predicted_control)) / denominator)
    treated_residual_correction = float(
        np.mean(treatment * overlap / propensity * (observed - predicted_treated)) / denominator
    )
    control_residual_correction = float(
        np.mean(-(1 - treatment) * overlap / (1 - propensity) * (observed - predicted_control)) / denominator
    )
    return {
        "n": int(len(modeled)),
        "model_contrast": model_contrast,
        "treated_residual_correction": treated_residual_correction,
        "control_residual_correction": control_residual_correction,
        "aipw_estimate": model_contrast + treated_residual_correction + control_residual_correction,
    }


def _aipw_overlap_risk_ratio(
    dataframe: pd.DataFrame,
    *,
    outcome: str,
    outcome_rhs: str,
) -> dict[str, float | str]:
    modeled, predicted_treated, predicted_control = _aipw_components(
        dataframe,
        outcome=outcome,
        outcome_rhs=outcome_rhs,
        family=sm.families.Binomial(),
    )
    treatment = modeled["atract"].to_numpy(dtype=int)
    observed = modeled[outcome].to_numpy(dtype=float)
    propensity = modeled["propensity_score"].to_numpy(dtype=float)
    overlap = propensity * (1 - propensity)
    denominator = float(np.mean(overlap))
    treated_terms = overlap * predicted_treated + treatment * overlap / propensity * (observed - predicted_treated)
    control_terms = overlap * predicted_control + (1 - treatment) * overlap / (1 - propensity) * (observed - predicted_control)
    treated_mean = float(np.mean(treated_terms) / denominator)
    control_mean = float(np.mean(control_terms) / denominator)
    risk_ratio = treated_mean / control_mean
    treated_influence = (treated_terms - treated_mean * overlap) / denominator
    control_influence = (control_terms - control_mean * overlap) / denominator
    log_influence = treated_influence / treated_mean - control_influence / control_mean
    standard_error = float(np.std(log_influence, ddof=1) / np.sqrt(len(modeled)))
    log_ratio = float(np.log(risk_ratio))
    z_value = log_ratio / standard_error if standard_error else np.nan
    p_value = float(2 * (1 - norm.cdf(abs(z_value)))) if np.isfinite(z_value) else np.nan
    risk_difference = treated_mean - control_mean
    risk_difference_influence = treated_influence - control_influence
    risk_difference_se = float(np.std(risk_difference_influence, ddof=1) / np.sqrt(len(modeled)))
    return {
        "term": "aipw_overlap",
        "scale": "risk_ratio",
        "estimate": float(risk_ratio),
        "ci_lower": float(np.exp(log_ratio - 1.96 * standard_error)),
        "ci_upper": float(np.exp(log_ratio + 1.96 * standard_error)),
        "p_value": p_value,
        "treated_risk": treated_mean,
        "control_risk": control_mean,
        "risk_difference": risk_difference,
        "risk_difference_ci_lower": float(risk_difference - 1.96 * risk_difference_se),
        "risk_difference_ci_upper": float(risk_difference + 1.96 * risk_difference_se),
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
    max_year_gap: int | None = None,
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
        if max_year_gap is not None:
            year_gap = np.abs(
                treated_group["study_year_index"].to_numpy()[repeated_treated][:, None]
                - control_group["study_year_index"].to_numpy()[None, :]
            )
            cost[year_gap > max_year_gap] = 1e6
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

    interaction_model = smf.ols(
        f"speed_mm2_min ~ atract * large_lesion + {augment_rhs}",
        data=matched,
    ).fit(cov_type="cluster", cov_kwds={"groups": cluster_codes})
    small_lesions = matched.loc[matched["large_lesion"].eq(0)].copy()
    large_lesions = matched.loc[matched["large_lesion"].eq(1)].copy()

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
                "analysis": f"small_lesion_<{LARGE_LESION_CUTOFF_MM}mm",
                "outcome": "speed_mm2_min",
                "n": int(len(small_lesions)),
                **_extract_linear_combination(
                    interaction_model,
                    ["atract"],
                    label="atract",
                    scale="mean_difference",
                ),
            },
            {
                "method": "ps_nn",
                "analysis": f"large_lesion_>={LARGE_LESION_CUTOFF_MM}mm",
                "outcome": "speed_mm2_min",
                "n": int(len(large_lesions)),
                **_extract_linear_combination(
                    interaction_model,
                    ["atract", "atract:large_lesion"],
                    label="atract_plus_interaction",
                    scale="mean_difference",
                ),
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
        "n_small_lesion": int(len(small_lesions)),
        "n_large_lesion": int(len(large_lesions)),
        "n_match_groups": int(matched["_match_group"].nunique()),
        "top_operator_treated_share": _max_operator_treated_share(matched),
    }
    return effects, metadata


def _select_matched_design(
    scored: pd.DataFrame,
    *,
    balance_covariates: list[str],
    max_year_gap: int | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, float]:
    candidate_rows: list[dict[str, float | int | bool]] = []
    payloads: dict[float, tuple[pd.DataFrame, pd.DataFrame]] = {}
    for caliper_multiplier in PRIMARY_MATCH_CALIPER_GRID:
        matched = _nearest_neighbor_match(
            scored,
            ratio=NN_MATCH_RATIO,
            caliper_multiplier=float(caliper_multiplier),
            match_scope=NN_MATCH_SCOPE,
            max_year_gap=max_year_gap,
        )
        if matched.empty:
            continue
        balance = build_matching_balance_table(scored, matched, balance_covariates)
        max_abs_smd = float(balance["adjusted_smd"].abs().max())
        top_operator_share = _max_operator_treated_share(matched)
        payloads[float(caliper_multiplier)] = (matched, balance)
        candidate_rows.append(
            {
                "caliper_multiplier": float(caliper_multiplier),
                "max_year_gap": np.nan if max_year_gap is None else int(max_year_gap),
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
    return matched_frame, balance_table, matching_grid, selected_caliper


def _build_matched_pair_diagnostics(matched: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, float | int]]:
    matched = add_large_lesion_flag(matched)
    rows: list[dict[str, float | int | str | bool]] = []
    for match_group, group in matched.groupby("_match_group", sort=True):
        treated = group.loc[group["atract"].eq(1)].iloc[0]
        control = group.loc[group["atract"].eq(0)].iloc[0]
        treated_large = int(treated["large_lesion"])
        control_large = int(control["large_lesion"])
        if treated_large == control_large == 0:
            size_pattern = "intact_<50mm"
        elif treated_large == control_large == 1:
            size_pattern = "intact_>=50mm"
        else:
            size_pattern = "split_size"
        year_gap = int(abs(int(treated["study_year_index"]) - int(control["study_year_index"])))
        rows.append(
            {
                "match_group": match_group,
                "operator_id_public": treated["operator_id_public"],
                "treated_year_index": int(treated["study_year_index"]),
                "control_year_index": int(control["study_year_index"]),
                "absolute_year_gap": year_gap,
                "same_year": year_gap == 0,
                "treated_large_lesion": treated_large,
                "control_large_lesion": control_large,
                "size_pair_pattern": size_pattern,
            }
        )
    pair_frame = pd.DataFrame(rows)
    treatment_by_size = pd.crosstab(matched["atract"], matched["large_lesion"])
    summary = {
        "n_pairs": int(len(pair_frame)),
        "intact_small_pairs": int(pair_frame["size_pair_pattern"].eq("intact_<50mm").sum()),
        "intact_large_pairs": int(pair_frame["size_pair_pattern"].eq("intact_>=50mm").sum()),
        "split_size_pairs": int(pair_frame["size_pair_pattern"].eq("split_size").sum()),
        "same_year_pairs": int(pair_frame["same_year"].sum()),
        "year_gap_le_1_pairs": int(pair_frame["absolute_year_gap"].le(PRIMARY_TEMPORAL_MAX_YEAR_GAP).sum()),
        "max_year_gap": int(pair_frame["absolute_year_gap"].max()),
        "median_year_gap": float(pair_frame["absolute_year_gap"].median()),
        "n_atract_small": int(treatment_by_size.loc[1, 0]) if 1 in treatment_by_size.index and 0 in treatment_by_size.columns else 0,
        "n_atract_large": int(treatment_by_size.loc[1, 1]) if 1 in treatment_by_size.index and 1 in treatment_by_size.columns else 0,
        "n_conventional_small": int(treatment_by_size.loc[0, 0]) if 0 in treatment_by_size.index and 0 in treatment_by_size.columns else 0,
        "n_conventional_large": int(treatment_by_size.loc[0, 1]) if 0 in treatment_by_size.index and 1 in treatment_by_size.columns else 0,
    }
    return pair_frame, summary


def _speed_effect_point_estimates(matched: pd.DataFrame, augment_rhs: str) -> dict[str, float]:
    matched = add_large_lesion_flag(matched)
    overall_model = smf.ols(
        f"speed_mm2_min ~ atract + {augment_rhs}",
        data=matched,
    ).fit()
    interaction_model = smf.ols(
        f"speed_mm2_min ~ atract * large_lesion + {augment_rhs}",
        data=matched,
    ).fit()
    small = float(interaction_model.params["atract"])
    interaction = float(interaction_model.params["atract:large_lesion"])
    return {
        "overall": float(overall_model.params["atract"]),
        f"small_lesion_<{LARGE_LESION_CUTOFF_MM}mm": small,
        f"large_lesion_>={LARGE_LESION_CUTOFF_MM}mm": small + interaction,
        "interaction_large_lesion": interaction,
    }


def _bootstrap_matched_speed_effects(
    matched: pd.DataFrame,
    augment_rhs: str,
    *,
    iterations: int,
    seed: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    groups = np.array(sorted(matched["_match_group"].unique()))
    observed = _speed_effect_point_estimates(matched, augment_rhs)
    draws: dict[str, list[float]] = {key: [] for key in observed}
    grouped = {group: frame for group, frame in matched.groupby("_match_group", sort=False)}

    for _ in range(iterations):
        sampled_groups = rng.choice(groups, size=len(groups), replace=True)
        sampled_frames = []
        for position, group in enumerate(sampled_groups):
            frame = grouped[group].copy()
            frame["_match_group"] = f"bootstrap_{position}"
            sampled_frames.append(frame)
        sampled = pd.concat(sampled_frames, ignore_index=True)
        try:
            estimates = _speed_effect_point_estimates(sampled, augment_rhs)
        except Exception:
            continue
        for key, value in estimates.items():
            draws[key].append(value)

    rows: list[dict[str, float | int | str]] = []
    for analysis, estimate in observed.items():
        values = np.asarray(draws[analysis], dtype=float)
        rows.append(
            {
                "method": "ps_nn_matched_set_bootstrap",
                "analysis": analysis,
                "outcome": "speed_mm2_min",
                "estimate": estimate,
                "ci_lower": float(np.percentile(values, 2.5)),
                "ci_upper": float(np.percentile(values, 97.5)),
                "n_bootstrap": int(len(values)),
                "seed": int(seed),
            }
        )
    return pd.DataFrame(rows)


def _continuous_size_effects(matched: pd.DataFrame) -> pd.DataFrame:
    matched = add_large_lesion_flag(matched)
    cluster_codes = matched["_match_group"].astype("category").cat.codes
    model = smf.ols(
        "speed_mm2_min ~ atract * bs(major_diameter_mm, df=4, include_intercept=False) + C(study_year_index)",
        data=matched,
    ).fit(cov_type="cluster", cov_kwds={"groups": cluster_codes})
    lower = float(matched["major_diameter_mm"].quantile(0.05))
    upper = float(matched["major_diameter_mm"].quantile(0.95))
    grid = np.linspace(lower, upper, 25)
    design_info = model.model.data.design_info
    covariance = model.cov_params().to_numpy()
    coefficients = model.params.to_numpy()
    rows: list[dict[str, float | str]] = []
    for size in grid:
        treated = matched.copy()
        control = matched.copy()
        treated["major_diameter_mm"] = size
        control["major_diameter_mm"] = size
        treated["atract"] = 1
        control["atract"] = 0
        treated_design = np.asarray(build_design_matrices([design_info], treated)[0])
        control_design = np.asarray(build_design_matrices([design_info], control)[0])
        contrast = (treated_design - control_design).mean(axis=0)
        estimate = float(contrast @ coefficients)
        standard_error = float(np.sqrt(max(contrast @ covariance @ contrast, 0.0)))
        rows.append(
            {
                "method": "ps_nn_continuous_size_interaction",
                "outcome": "speed_mm2_min",
                "major_diameter_mm": float(size),
                "estimate": estimate,
                "ci_lower": estimate - 1.96 * standard_error,
                "ci_upper": estimate + 1.96 * standard_error,
            }
        )
    return pd.DataFrame(rows)


def _run_temporal_speed_sensitivity(
    scored: pd.DataFrame,
    matched_frame: pd.DataFrame,
    augment_rhs: str,
    balance_covariates: list[str],
) -> dict[str, pd.DataFrame]:
    pair_diagnostics, _ = _build_matched_pair_diagnostics(matched_frame)
    retained_groups = set(
        pair_diagnostics.loc[
            pair_diagnostics["absolute_year_gap"].le(PRIMARY_TEMPORAL_MAX_YEAR_GAP),
            "match_group",
        ]
    )
    restricted = matched_frame.loc[matched_frame["_match_group"].isin(retained_groups)].copy()
    restricted_effects, _ = _fit_primary_speed_models(restricted, augment_rhs)
    restricted_effects["sensitivity"] = f"existing_pairs_year_gap_le_{PRIMARY_TEMPORAL_MAX_YEAR_GAP}"
    restricted_effects["caliper_multiplier"] = np.nan
    restricted_effects["max_abs_smd"] = float(
        build_matching_balance_table(scored, restricted, balance_covariates)["adjusted_smd"].abs().max()
    )
    restricted_effects["n_match_groups"] = int(restricted["_match_group"].nunique())

    rematched, rematched_balance, rematching_grid, selected_caliper = _select_matched_design(
        scored,
        balance_covariates=balance_covariates,
        max_year_gap=PRIMARY_TEMPORAL_MAX_YEAR_GAP,
    )
    rematched_effects, _ = _fit_primary_speed_models(rematched, augment_rhs)
    rematched_effects["sensitivity"] = f"rematched_year_gap_le_{PRIMARY_TEMPORAL_MAX_YEAR_GAP}"
    rematched_effects["caliper_multiplier"] = float(selected_caliper)
    rematched_effects["max_abs_smd"] = float(rematched_balance["adjusted_smd"].abs().max())
    rematched_effects["n_match_groups"] = int(rematched["_match_group"].nunique())

    first_atract_year = scored.loc[scored["atract"].eq(1)].groupby("operator_id_public")["study_year_index"].min()
    contemporary = scored.loc[
        scored["study_year_index"].ge(scored["operator_id_public"].map(first_atract_year))
    ].copy()
    contemporary_matched, contemporary_balance, contemporary_grid, contemporary_caliper = _select_matched_design(
        contemporary,
        balance_covariates=balance_covariates,
    )
    contemporary_effects, _ = _fit_primary_speed_models(contemporary_matched, augment_rhs)
    contemporary_effects["sensitivity"] = "rematched_post_adoption_period"
    contemporary_effects["caliper_multiplier"] = float(contemporary_caliper)
    contemporary_effects["max_abs_smd"] = float(contemporary_balance["adjusted_smd"].abs().max())
    contemporary_effects["n_match_groups"] = int(contemporary_matched["_match_group"].nunique())
    contemporary_effects["contemporary_supported_n"] = int(len(contemporary))

    summary = pd.concat([restricted_effects, rematched_effects, contemporary_effects], ignore_index=True)
    return {
        "summary": summary,
        "rematching_grid": rematching_grid,
        "contemporary_rematching_grid": contemporary_grid,
    }


def run_primary_speed_analysis(
    dataframe: pd.DataFrame,
    *,
    bootstrap_iterations: int = PRIMARY_BOOTSTRAP_ITERATIONS,
) -> dict[str, object]:
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

    matched_frame, balance_table, matching_grid, selected_caliper = _select_matched_design(
        scored,
        balance_covariates=list(spec["balance_covariates"]),
    )
    effects, metadata = _fit_primary_speed_models(matched_frame, str(spec["speed_matched_rhs"]))
    pair_diagnostics, pair_summary = _build_matched_pair_diagnostics(matched_frame)
    temporal_sensitivity = _run_temporal_speed_sensitivity(
        scored,
        matched_frame,
        str(spec["speed_matched_rhs"]),
        list(spec["balance_covariates"]),
    )
    bootstrap_results = _bootstrap_matched_speed_effects(
        matched_frame,
        str(spec["speed_matched_rhs"]),
        iterations=bootstrap_iterations,
        seed=PRIMARY_BOOTSTRAP_SEED,
    )
    continuous_size = _continuous_size_effects(matched_frame)
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
            "pair_diagnostics": pair_summary,
            "bootstrap_iterations": bootstrap_iterations,
            "bootstrap_seed": PRIMARY_BOOTSTRAP_SEED,
        }
    )
    return {
        "scored_frame": scored,
        "matched_frame": matched_frame,
        "propensity_model": propensity_model,
        "balance": balance_table,
        "effects": effects,
        "matching_grid": matching_grid,
        "pair_diagnostics": pair_diagnostics,
        "temporal_sensitivity": temporal_sensitivity,
        "bootstrap": bootstrap_results,
        "continuous_size": continuous_size,
        "metadata": metadata,
    }


def _fit_weighted_speed_effects(weighted: pd.DataFrame, augment_rhs: str) -> tuple[pd.DataFrame, dict[str, float | int | str]]:
    weighted = add_large_lesion_flag(weighted)
    overall_effect, _ = _aipw_overlap_mean_difference(
        weighted,
        outcome="speed_mm2_min",
        outcome_rhs=augment_rhs,
        analysis="overall",
    )
    small_effect, small_influence = _aipw_overlap_mean_difference(
        weighted,
        outcome="speed_mm2_min",
        outcome_rhs=augment_rhs,
        analysis=f"small_lesion_<{LARGE_LESION_CUTOFF_MM}mm_interaction_model",
        mask=weighted["large_lesion"].eq(0),
        include_large_interaction=True,
    )
    large_effect, large_influence = _aipw_overlap_mean_difference(
        weighted,
        outcome="speed_mm2_min",
        outcome_rhs=augment_rhs,
        analysis=f"large_lesion_>={LARGE_LESION_CUTOFF_MM}mm_interaction_model",
        mask=weighted["large_lesion"].eq(1),
        include_large_interaction=True,
    )
    interaction_estimate = float(large_effect["estimate"] - small_effect["estimate"])
    interaction_effect = _effect_from_influence(
        interaction_estimate,
        large_influence - small_influence,
        "interaction_difference",
        "interaction_large_lesion",
    )
    large_weighted = weighted.loc[weighted["large_lesion"].eq(1)].copy()
    large_restricted_effect, _ = _aipw_overlap_mean_difference(
        large_weighted,
        outcome="speed_mm2_min",
        outcome_rhs=augment_rhs,
        analysis=f"large_lesion_>={LARGE_LESION_CUTOFF_MM}mm_restricted",
    )

    effects = pd.DataFrame(
        [
            {
                "method": "ow_dr",
                "analysis": "overall",
                "outcome": "speed_mm2_min",
                "n": int(len(weighted)),
                **overall_effect,
            },
            {
                "method": "ow_dr",
                "analysis": f"small_lesion_<{LARGE_LESION_CUTOFF_MM}mm_interaction_model",
                "outcome": "speed_mm2_min",
                "n": int(weighted["large_lesion"].eq(0).sum()),
                **small_effect,
            },
            {
                "method": "ow_dr",
                "analysis": f"large_lesion_>={LARGE_LESION_CUTOFF_MM}mm_interaction_model",
                "outcome": "speed_mm2_min",
                "n": int(weighted["large_lesion"].eq(1).sum()),
                **large_effect,
            },
            {
                "method": "ow_dr",
                "analysis": "interaction_large_lesion",
                "outcome": "speed_mm2_min",
                "n": int(len(weighted)),
                **interaction_effect,
            },
            {
                "method": "ow_dr",
                "analysis": f"large_lesion_>={LARGE_LESION_CUTOFF_MM}mm_restricted",
                "outcome": "speed_mm2_min",
                "n": int(len(large_weighted)),
                **large_restricted_effect,
            },
        ]
    )
    metadata = {
        "analysis": "speed_overlap_weighting_robustness",
        "n_complete_case": int(len(weighted)),
        "n_small_lesion": int(weighted["large_lesion"].eq(0).sum()),
        "n_large_lesion": int(weighted["large_lesion"].eq(1).sum()),
        "weight_min": float(weighted["analysis_weight"].min()),
        "weight_max": float(weighted["analysis_weight"].max()),
        "aipw_decomposition_overall": _aipw_overlap_mean_decomposition(
            weighted,
            outcome="speed_mm2_min",
            outcome_rhs=augment_rhs,
        ),
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
    effects, metadata = _fit_weighted_speed_effects(weighted, str(spec["speed_dr_rhs"]))
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
    robustness_overall = effects.loc[effects["analysis"].eq("overall")].iloc[0].to_dict()
    robustness_overall.update(
        {
            "method_label": "OW + DR",
            "max_abs_smd": float(balance["adjusted_smd"].abs().max()),
            "weight_max": float(metadata["weight_max"]),
            "n_match_groups": np.nan,
        }
    )
    robustness_supplement = effects.loc[~effects["analysis"].eq("overall")].copy()
    robustness_supplement["method_label"] = "OW + DR"
    robustness_supplement["max_abs_smd"] = float(balance["adjusted_smd"].abs().max())
    robustness_supplement["weight_max"] = float(metadata["weight_max"])
    robustness_supplement["n_match_groups"] = np.nan

    return {
        "weighted_frame": weighted,
        "propensity_model": propensity_model,
        "balance": balance,
        "effects": effects,
        "metadata": metadata,
        "summary": pd.concat([pd.DataFrame([primary_overall, robustness_overall]), robustness_supplement], ignore_index=True),
    }


def _fit_weighted_binary_effect(
    weighted: pd.DataFrame,
    *,
    outcome: str,
    augment_rhs: str,
) -> dict[str, float | str]:
    return _aipw_overlap_risk_ratio(weighted, outcome=outcome, outcome_rhs=augment_rhs)


def _effective_sample_size(weights: pd.Series) -> float:
    values = weights.to_numpy(dtype=float)
    return float(values.sum() ** 2 / np.square(values).sum())


def _binary_event_summary(dataframe: pd.DataFrame, outcome: str) -> dict[str, float | int]:
    treated = dataframe.loc[dataframe["atract"].eq(1), outcome]
    control = dataframe.loc[dataframe["atract"].eq(0), outcome]
    treated_n = int(treated.notna().sum())
    control_n = int(control.notna().sum())
    treated_events = int(treated.fillna(0).sum())
    control_events = int(control.fillna(0).sum())
    treated_risk_observed = treated_events / treated_n if treated_n else np.nan
    control_risk_observed = control_events / control_n if control_n else np.nan
    return {
        "treated_n": treated_n,
        "treated_events": treated_events,
        "treated_risk_observed": float(treated_risk_observed),
        "control_n": control_n,
        "control_events": control_events,
        "control_risk_observed": float(control_risk_observed),
        "observed_risk_difference": float(treated_risk_observed - control_risk_observed),
    }


def run_primary_binary_analyses(dataframe: pd.DataFrame) -> dict[str, object]:
    spec = CAUSAL_SPECS[PRIMARY_CAUSAL_SPEC]
    rows: list[dict[str, float | str]] = []
    diagnostics: dict[str, dict[str, float | int | str]] = {}
    balance_rows: list[pd.DataFrame] = []

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
        balance_with_outcome = balance.copy()
        balance_with_outcome.insert(0, "outcome", outcome)
        balance_rows.append(balance_with_outcome)
        rows.append(
            {
                "method": "ow_dr",
                "analysis": "overall",
                "outcome": outcome,
                "n": int(len(weighted)),
                **_binary_event_summary(weighted, outcome),
                **_fit_weighted_binary_effect(weighted, outcome=outcome, augment_rhs=str(spec["binary_aug_rhs"])),
            }
        )
        treated_weights = weighted.loc[weighted["atract"].eq(1), "analysis_weight"]
        control_weights = weighted.loc[weighted["atract"].eq(0), "analysis_weight"]
        diagnostics[outcome] = {
            "n_complete_case": int(len(weighted)),
            "weight_min": float(weighted["analysis_weight"].min()),
            "weight_max": float(weighted["analysis_weight"].max()),
            "effective_sample_size": _effective_sample_size(weighted["analysis_weight"]),
            "treated_effective_sample_size": _effective_sample_size(treated_weights),
            "control_effective_sample_size": _effective_sample_size(control_weights),
            "max_abs_weighted_smd": float(balance["adjusted_smd"].abs().max()),
            "spec_key": PRIMARY_CAUSAL_SPEC,
            "support_scope": PRIMARY_SUPPORT_SCOPE,
            "min_per_arm": PRIMARY_OPERATOR_YEAR_MIN_PER_ARM,
            "ps_structure": PRIMARY_PS_STRUCTURE,
        }

    return {
        "effects": pd.DataFrame(rows),
        "diagnostics": diagnostics,
        "balance": pd.concat(balance_rows, ignore_index=True),
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
                "effective_sample_size": _effective_sample_size(weighted["analysis_weight"]),
                **_binary_event_summary(weighted, outcome),
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
                "effective_sample_size": np.nan,
                **_binary_event_summary(matched, outcome),
                **_extract_exponentiated_effect(matched_model, "atract", "risk_ratio"),
            }
        )

    return {"summary": pd.DataFrame(rows)}


def _package_versions() -> dict[str, str]:
    packages = ["matplotlib", "numpy", "openpyxl", "pandas", "patsy", "scipy", "statsmodels"]
    resolved: dict[str, str] = {"python": platform.python_version()}
    for package in packages:
        try:
            resolved[package] = version(package)
        except PackageNotFoundError:
            resolved[package] = "not-installed"
    return resolved


def _build_results_manifest(
    *,
    dataset_checksum: str,
    primary_speed: dict[str, object],
    speed_robustness: dict[str, object],
    primary_binary: dict[str, object],
) -> pd.DataFrame:
    rows: list[dict[str, str | int | float]] = []
    manifest_sources = [
        (
            "table_2_main_results",
            "results/tables/table_2_main_results.csv",
            primary_speed["effects"].assign(method_label="PS-NN"),
        ),
        (
            "table_s2_speed_robustness",
            "results/tables/table_s2_speed_robustness.csv",
            speed_robustness["summary"],
        ),
        (
            "primary_binary_results",
            "results/model_summaries/primary_binary_results.csv",
            primary_binary["effects"].assign(method_label="OW + DR"),
        ),
    ]
    for table_name, output_file, frame in manifest_sources:
        for _, row in frame.iterrows():
            rows.append(
                {
                    "dataset_sha256": dataset_checksum,
                    "table_or_source": table_name,
                    "output_file": output_file,
                    "population": row.get("analysis", "overall"),
                    "method": row.get("method_label", row.get("method", "")),
                    "outcome": row.get("outcome", ""),
                    "n": int(row.get("n", 0)),
                    "estimate": float(row.get("estimate", np.nan)),
                    "ci_lower": float(row.get("ci_lower", np.nan)),
                    "ci_upper": float(row.get("ci_upper", np.nan)),
                    "p_value": float(row.get("p_value", np.nan)),
                    "scale": row.get("scale", ""),
                    "script_entrypoint": "python -m atract_analysis analysis",
                }
            )
    return pd.DataFrame(rows)


def write_model_outputs(
    results_dir: Path,
    primary_speed: dict[str, object],
    speed_robustness: dict[str, object],
    primary_binary: dict[str, object],
    binary_comparison: dict[str, object],
    *,
    dataset_checksum: str = EXPECTED_PUBLIC_DATASET_SHA256,
) -> None:
    results_dir.mkdir(parents=True, exist_ok=True)
    stale_paths = [
        "binary_effects.csv",
        "binary_sensitivity_summary.csv",
        "speed_sensitivity_summary.csv",
        "speed_effects.csv",
        "speed_iptw_frontier.csv",
        "speed_nn_frontier.csv",
        "matched_pair_diagnostics.csv",
        "speed_temporal_sensitivity.csv",
        "speed_temporal_rematching_grid.csv",
        "speed_contemporary_rematching_grid.csv",
        "speed_bootstrap_results.csv",
        "speed_continuous_size_effects.csv",
        "binary_balance_diagnostics.csv",
        "manuscript_results_manifest.csv",
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
    primary_speed["pair_diagnostics"].to_csv(results_dir / "matched_pair_diagnostics.csv", index=False)
    primary_speed["temporal_sensitivity"]["summary"].to_csv(results_dir / "speed_temporal_sensitivity.csv", index=False)
    primary_speed["temporal_sensitivity"]["rematching_grid"].to_csv(
        results_dir / "speed_temporal_rematching_grid.csv",
        index=False,
    )
    primary_speed["temporal_sensitivity"]["contemporary_rematching_grid"].to_csv(
        results_dir / "speed_contemporary_rematching_grid.csv",
        index=False,
    )
    primary_speed["bootstrap"].to_csv(results_dir / "speed_bootstrap_results.csv", index=False)
    primary_speed["continuous_size"].to_csv(results_dir / "speed_continuous_size_effects.csv", index=False)
    primary_binary["balance"].to_csv(results_dir / "binary_balance_diagnostics.csv", index=False)
    _build_results_manifest(
        dataset_checksum=dataset_checksum,
        primary_speed=primary_speed,
        speed_robustness=speed_robustness,
        primary_binary=primary_binary,
    ).to_csv(results_dir / "manuscript_results_manifest.csv", index=False)

    metadata = {
        "dataset": {
            "sha256": dataset_checksum,
            "expected_sha256": EXPECTED_PUBLIC_DATASET_SHA256,
        },
        "software": _package_versions(),
        "formulas": {
            "propensity_score": primary_speed["propensity_model"].model.formula,
            "speed_matched_outcome": f"speed_mm2_min ~ atract + {CAUSAL_SPECS[PRIMARY_CAUSAL_SPEC]['speed_matched_rhs']}",
            "speed_matched_interaction": f"speed_mm2_min ~ atract * large_lesion + {CAUSAL_SPECS[PRIMARY_CAUSAL_SPEC]['speed_matched_rhs']}",
            "speed_overlap_augmented_outcome": f"speed_mm2_min ~ atract + {CAUSAL_SPECS[PRIMARY_CAUSAL_SPEC]['speed_dr_rhs']}",
            "binary_augmented_outcome": str(CAUSAL_SPECS[PRIMARY_CAUSAL_SPEC]["binary_aug_rhs"]),
        },
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
