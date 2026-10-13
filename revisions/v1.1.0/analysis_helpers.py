from __future__ import annotations

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from patsy import build_design_matrices

from atract_analysis.config import CAUSAL_SPECS, PRIMARY_CAUSAL_SPEC
from atract_analysis.models import (
    _aipw_overlap_mean_difference,
    _continuous_size_effects,
    _effect_from_influence,
    _extract_linear_combination,
    _extract_linear_effect,
    _nearest_neighbor_match,
    build_balance_table,
    build_matching_balance_table,
)

SPEC = CAUSAL_SPECS[PRIMARY_CAUSAL_SPEC]
SPEED_RHS = SPEC["speed_matched_rhs"]
CORE_RHS = SPEC["speed_dr_rhs"]
DURATION_RHS = "C(study_year_index) + bs(surface_mm2, df=4, include_intercept=False)"
STRATUM_LABELS = {
    "overall": "overall",
    "small_lesion_<50mm": "<50 mm",
    "large_lesion_>=50mm": ">=50 mm",
    "interaction_large_lesion": "size interaction",
}


def fit_matched_model(formula: str, frame: pd.DataFrame):
    return smf.ols(formula, data=frame, missing="raise").fit(
        cov_type="cluster", cov_kwds={"groups": frame["_match_group"]}
    )


def matched_effects(frame: pd.DataFrame, outcome: str, rhs: str, analysis: str) -> pd.DataFrame:
    overall_formula = f"{outcome} ~ atract + {rhs}"
    interaction_formula = f"{outcome} ~ atract * large_lesion + {rhs}"
    overall = fit_matched_model(overall_formula, frame)
    interaction = fit_matched_model(interaction_formula, frame)
    rows = []
    for stratum, subset, effect in [
        ("overall", frame, _extract_linear_effect(overall, "atract", "mean_difference")),
        ("<50 mm", frame.loc[frame.large_lesion.eq(0)], _extract_linear_combination(interaction, ["atract"], label="atract", scale="mean_difference")),
        (">=50 mm", frame.loc[frame.large_lesion.eq(1)], _extract_linear_combination(interaction, ["atract", "atract:large_lesion"], label="atract_plus_interaction", scale="mean_difference")),
        ("size interaction", frame, _extract_linear_effect(interaction, "atract:large_lesion", "interaction_difference")),
    ]:
        rows.append({
            "analysis": analysis, "stratum": stratum, "outcome": outcome,
            "n": len(subset), "n_model": len(frame), "n_pairs_model": frame._match_group.nunique(),
            "treated_n": int(subset.atract.eq(1).sum()), "control_n": int(subset.atract.eq(0).sum()),
            "formula": overall_formula if stratum == "overall" else interaction_formula,
            "unit": "min" if outcome == "procedure_duration_min" else "mm2/min",
            **effect,
        })
    return pd.DataFrame(rows)


def fibrosis_models(matched: pd.DataFrame):
    complete = matched.groupby("_match_group").fibrosis.apply(lambda values: values.notna().all())
    retained = matched.loc[matched._match_group.isin(complete.index[complete])].copy()
    specifications = [
        ("matched_primary", matched, SPEED_RHS),
        ("matched_all_pairs_full_covariates_no_fibrosis", matched, CORE_RHS),
        ("fibrosis_complete_pairs_year_size", retained, SPEED_RHS),
        ("fibrosis_complete_pairs_year_size_fibrosis", retained, SPEED_RHS + " + C(fibrosis)"),
        ("fibrosis_complete_pairs_full_covariates", retained, CORE_RHS),
        ("fibrosis_complete_pairs_full_covariates_fibrosis", retained, CORE_RHS + " + C(fibrosis)"),
    ]
    effects = pd.concat([
        matched_effects(frame, "speed_mm2_min", rhs, name) for name, frame, rhs in specifications
    ], ignore_index=True)
    return effects, retained


def overlap_fibrosis_models(scored: pd.DataFrame) -> pd.DataFrame:
    rows = []
    complete = scored.loc[scored.fibrosis.notna()].copy()
    for name, frame, rhs in [
        ("overlap_primary", scored, CORE_RHS),
        ("overlap_fibrosis_complete_cases", complete, CORE_RHS),
        ("overlap_fibrosis_complete_cases_adjusted", complete, CORE_RHS + " + C(fibrosis)"),
    ]:
        effects = {}
        for stratum, mask in [("overall", None), ("<50 mm", frame.large_lesion.eq(0)), (">=50 mm", frame.large_lesion.eq(1))]:
            effect, influence = _aipw_overlap_mean_difference(
                frame, outcome="speed_mm2_min", outcome_rhs=rhs,
                analysis=stratum, mask=mask, include_large_interaction=mask is not None,
            )
            effects[stratum] = (effect, influence)
            subset = frame if mask is None else frame.loc[mask]
            weights = subset.analysis_weight.to_numpy()
            rows.append({
                "analysis": name, "stratum": stratum, "outcome": "speed_mm2_min",
                "n": len(subset), "n_model": len(frame),
                "treated_n": int(subset.atract.eq(1).sum()), "control_n": int(subset.atract.eq(0).sum()),
                "effective_sample_size": weights.sum() ** 2 / np.square(weights).sum(),
                "formula": f"speed_mm2_min ~ atract{' * large_lesion' if mask is not None else ''} + {rhs}",
                "unit": "mm2/min", **effect,
            })
        small, large = effects["<50 mm"], effects[">=50 mm"]
        interaction = _effect_from_influence(
            large[0]["estimate"] - small[0]["estimate"], large[1] - small[1],
            "interaction_difference", "interaction_large_lesion",
        )
        rows.append({"analysis": name, "stratum": "size interaction", "outcome": "speed_mm2_min", "n": len(frame), "n_model": len(frame), "unit": "mm2/min", **interaction})
    return pd.DataFrame(rows)


def operator_volume_models(matched: pd.DataFrame, volume_groups: pd.Series):
    rows, balance_rows = [], []
    for name, frame in [
        ("operator_volume_full_covariates", matched.copy()),
        ("operator_volume_exclude_pooled_stratum", matched.loc[matched.operator_id_public.ne("operator_other")].copy()),
    ]:
        frame["high_volume"] = frame.operator_id_public.map(volume_groups).eq("high_volume").astype(int)
        # Operator fixed effects already contain the volume-group main effect.
        formula = f"speed_mm2_min ~ atract + atract:high_volume + {CORE_RHS}"
        model = fit_matched_model(formula, frame)
        lower = frame.loc[frame.high_volume.eq(0)]
        high = frame.loc[frame.high_volume.eq(1)]
        for stratum, subset, effect in [
            ("lower-volume trained operators", lower, _extract_linear_effect(model, "atract", "mean_difference")),
            ("high-volume operators", high, _extract_linear_combination(model, ["atract", "atract:high_volume"], label="high_volume_effect", scale="mean_difference")),
            ("operator-volume interaction", frame, _extract_linear_effect(model, "atract:high_volume", "interaction_difference")),
        ]:
            rows.append({
                "analysis": name, "stratum": stratum, "outcome": "speed_mm2_min",
                "n": len(subset), "n_model": len(frame), "n_pairs_model": frame._match_group.nunique(),
                "n_pairs_stratum": subset._match_group.nunique(),
                "treated_n": int(subset.atract.eq(1).sum()), "control_n": int(subset.atract.eq(0).sum()),
                "unit": "mm2/min", "formula": formula, **effect,
            })
        for label, subset in [("high-volume operators", high), ("lower-volume trained operators", lower)]:
            balance = build_balance_table(subset, SPEC["balance_covariates"])
            balance["analysis"], balance["stratum"] = name, label
            balance_rows.append(balance)
    return pd.DataFrame(rows), pd.concat(balance_rows, ignore_index=True)


def caliper_models(scored: pd.DataFrame):
    effects, balance_rows = [], []
    for caliper in [0.40, 0.60]:
        matched = _nearest_neighbor_match(scored, ratio=1, caliper_multiplier=caliper, match_scope="operator_only")
        balance = build_matching_balance_table(scored, matched, SPEC["balance_covariates"])
        maximum_smd = float(balance.adjusted_smd.abs().max())
        result = matched_effects(matched, "speed_mm2_min", SPEED_RHS, f"matching_caliper_{caliper:.2f}")
        result["caliper"] = caliper
        result["max_abs_smd"] = maximum_smd
        result["balance_below_0_10"] = maximum_smd < 0.10
        effects.append(result)
        balance["caliper"] = caliper
        balance_rows.append(balance)
    return pd.concat(effects, ignore_index=True), pd.concat(balance_rows, ignore_index=True)


def continuous_size_test(matched: pd.DataFrame):
    formula = "speed_mm2_min ~ atract * bs(major_diameter_mm, df=4, include_intercept=False) + C(study_year_index)"
    model = fit_matched_model(formula, matched)
    terms = [term for term in model.params.index if term.startswith("atract:bs(")]
    contrasts = np.zeros((len(terms), len(model.params)))
    for row, term in enumerate(terms):
        contrasts[row, model.params.index.get_loc(term)] = 1
    test = model.wald_test(contrasts, use_f=False, scalar=True)
    return _continuous_size_effects(matched), pd.DataFrame([{
        "analysis": "continuous_size_joint_interaction", "stratum": "continuous size interaction",
        "outcome": "speed_mm2_min", "n": len(matched), "n_model": len(matched),
        "estimate": float(test.statistic), "unit": "Wald chi-square statistic", "df": len(terms),
        "p_value": float(test.pvalue), "formula": formula, "test": "joint robust Wald chi-square",
        "null_hypothesis": "all four treatment-by-size-spline coefficients equal zero",
    }])


def operator_identity_models(matched: pd.DataFrame) -> pd.DataFrame:
    formula = "speed_mm2_min ~ atract * C(operator_id_public) + " + SPEED_RHS
    model = fit_matched_model(formula, matched)
    design = model.model.data.design_info
    rows = []
    for operator, frame in matched.groupby("operator_id_public"):
        treated, control = frame.copy(), frame.copy()
        treated["atract"], control["atract"] = 1, 0
        vector = (
            np.asarray(build_design_matrices([design], treated)[0])
            - np.asarray(build_design_matrices([design], control)[0])
        ).mean(axis=0)
        effect = model.t_test(vector)
        rows.append({
            "analysis": "individual_operator_exploratory", "stratum": operator,
            "outcome": "speed_mm2_min", "n": len(frame), "n_model": len(matched),
            "n_pairs_model": matched._match_group.nunique(), "n_pairs_stratum": frame._match_group.nunique(),
            "treated_n": int(frame.atract.eq(1).sum()), "control_n": int(frame.atract.eq(0).sum()),
            "estimate": float(effect.effect[0]), "ci_lower": float(effect.conf_int()[0][0]),
            "ci_upper": float(effect.conf_int()[0][1]), "p_value": float(effect.pvalue),
            "formula": formula, "unit": "mm2/min",
        })
    return pd.DataFrame(rows)
