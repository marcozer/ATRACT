# coding: utf-8
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
from importlib.metadata import version

REVISION_ROOT = Path(__file__).resolve().parent
REPO_ROOT = REVISION_ROOT.parents[1]
sys.dont_write_bytecode = True
sys.path.insert(0, str(REPO_ROOT / "src"))
os.environ.setdefault("MPLCONFIGDIR", "/tmp/atract-matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from atract_analysis.config import EXPECTED_PUBLIC_DATASET_SHA256, PRIMARY_PS_STRUCTURE
from atract_analysis.cohort import add_large_lesion_flag
from atract_analysis.models import (
    _analysis_cohort, _nearest_neighbor_match,
    build_balance_table, build_matching_balance_table, fit_propensity_model,
)
from analysis_helpers import (
    CORE_RHS, DURATION_RHS, SPEC, SPEED_RHS, STRATUM_LABELS,
    caliper_models, continuous_size_test, fibrosis_models, matched_effects,
    operator_identity_models, operator_volume_models, overlap_fibrosis_models,
)

BASELINE_COMMIT = "8a74c400ef2a26a3db93e395f634a8f1a4d30b86"
BLUE, ORANGE, DARK = "#5E81AC", "#D08770", "#2E3440"


def checksum(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_baseline() -> dict[str, str]:
    manifest = json.loads((REVISION_ROOT / "baseline_output_checksums.json").read_text())
    for relative, expected in manifest.items():
        if checksum(REPO_ROOT / relative) != expected:
            raise ValueError(f"Archived v1.0.0 output changed: {relative}")
    return manifest


def verify_primary(effects: pd.DataFrame) -> None:
    archived = pd.read_csv(REPO_ROOT / "results/model_summaries/primary_speed_results.csv")
    reference = effects.loc[effects.analysis.eq("matched_primary")].set_index("stratum")
    for column in ["n", "estimate", "ci_lower", "ci_upper", "p_value"]:
        actual = reference.loc[archived.analysis.map(STRATUM_LABELS), column]
        np.testing.assert_allclose(actual, archived[column], atol=1e-10, rtol=1e-10)


def save_figure(fig, output: Path, stem: str) -> None:
    for extension in ["svg", "png", "pdf"]:
        fig.savefig(output / f"{stem}.{extension}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def figures(curve: pd.DataFrame, joint: pd.DataFrame, volume: pd.DataFrame, output: Path) -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10, "text.color": DARK,
        "axes.spines.top": False, "axes.spines.right": False,
        "pdf.fonttype": 42, "svg.fonttype": "none", "svg.hashsalt": "atract-v1.1.0",
    })
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(curve.major_diameter_mm, curve.estimate, color=BLUE, linewidth=2)
    ax.fill_between(curve.major_diameter_mm, curve.ci_lower, curve.ci_upper, color=BLUE, alpha=.18, label="Pointwise 95% CI")
    ax.axhline(0, color=DARK, linewidth=.8)
    ax.axvline(50, color=ORANGE, linestyle="--", linewidth=1)
    ax.set(xlabel="Major diameter (mm)", ylabel="ATRACT - conventional speed (mm²/min)", title=f"Continuous size interaction: Wald test, 4 df, p={joint.p_value.iloc[0]:.3f}")
    ax.legend(frameon=False)
    fig.tight_layout()
    save_figure(fig, output, "continuous_size_interaction")
    frame = volume.loc[volume.analysis.eq("operator_volume_full_covariates") & volume.stratum.ne("operator-volume interaction")]
    fig, ax = plt.subplots(figsize=(7, 3.6))
    positions = np.arange(len(frame))
    ax.errorbar(frame.estimate, positions, xerr=[frame.estimate-frame.ci_lower, frame.ci_upper-frame.estimate], fmt="o", color=BLUE, capsize=3)
    ax.set_yticks(positions, [f"{row.stratum}\n({row.n_pairs_stratum} pairs)" for row in frame.itertuples()])
    ax.axvline(0, color=DARK, linewidth=.8)
    ax.set_xlabel("Speed difference (mm²/min), 95% CI")
    fig.tight_layout()
    save_figure(fig, output, "operator_volume_effects")


def effect_text(row) -> str:
    return f"{row.estimate:+.2f} (IC95 % {row.ci_lower:.2f} à {row.ci_upper:.2f} ; p = {row.p_value:.3f})"


def write_response(tables: dict[str, pd.DataFrame], output: Path) -> None:
    fibrosis = tables["fibrosis_effects"]
    def get(name, stratum="overall"):
        return fibrosis.loc[fibrosis.analysis.eq(name) & fibrosis.stratum.eq(stratum)].iloc[0]
    joint = tables["continuous_size_test"].iloc[0]
    volume = tables["operator_volume_effects"].query("analysis == 'operator_volume_full_covariates'")
    duration = tables["duration_effects"].set_index("stratum")
    calipers = tables["caliper_effects"]
    paragraphs = [
        "# Chiffres complémentaires pour la révision",
        "",
        "Les résultats de v1.0.0 sont conservés. Les analyses complémentaires utilisent le même jeu analytique et des erreurs standards clusterisées par paire pour les modèles appariés.",
        "",
        "## Fibrose",
        "",
        f"Sur les mêmes 372 paires, le modèle complet sans fibrose donne {effect_text(get('fibrosis_complete_pairs_full_covariates'))} mm²/min, contre {effect_text(get('fibrosis_complete_pairs_full_covariates_fibrosis'))} avec fibrose. Dans les grandes lésions, le modèle année + taille + fibrose donne {effect_text(get('fibrosis_complete_pairs_year_size_fibrosis', '>=50 mm'))} mm²/min.",
        "",
        f"La restriction seule, avec le modèle année + taille, donne {effect_text(get('fibrosis_complete_pairs_year_size'))} mm²/min. Avec toutes les covariables et sans fibrose, les 391 paires initiales donnent {effect_text(get('matched_all_pairs_full_covariates_no_fibrosis'))}. L'atténuation n'est donc pas expliquée par l'ajout de la fibrose. Elle est associée à l'ajustement supplémentaire ; la seule comparaison des modèles ne permet pas d'attribuer cette différence à la localisation ou à la restriction exclusivement.",
        "",
        "## Taille continue",
        "",
        f"Il s'agit d'un test de Wald conjoint, avec covariance robuste clusterisée par paire, sur les quatre coefficients d'interaction entre ATRACT et les quatre fonctions de base du B-spline de taille (χ² = {joint.estimate:.3f}, 4 degrés de liberté ; p = {joint.p_value:.6f}). La courbe est numériquement identique à table_s9 de v1.0.0. Les intervalles affichés sont ponctuels à 95 %, non simultanés.",
        "",
        "## Volume des opérateurs",
        "",
        "La classification est issue de la définition clinique fournie le 13 octobre 2026 : operator_01 et operator_02, environ 300 ESD/an et plus de 10 ans d'expérience ; les autres, tous formés, environ 50–100 ESD/an et environ 5 ans d'expérience. Cette classification fixe n'est pas une mesure annuelle de l'expérience au moment de chaque geste. Le modèle d'interaction ajuste sur les covariables préopératoires, y compris les effets fixes des opérateurs, sans fibrose.",
        "",
    ]
    for row in volume.itertuples():
        paragraphs.extend([f"{row.stratum} : n = {row.n}, {effect_text(row)} mm²/min.", ""])
    paragraphs.extend(["## Durée directement mesurée", ""])
    for name in ["overall", ">=50 mm"]:
        row = duration.loc[name]
        paragraphs.extend([f"{name} : différence de durée ATRACT moins traction conventionnelle = {effect_text(row)} minutes. Une valeur négative correspond à une durée plus courte.", ""])
    paragraphs.extend(["Le modèle ajuste sur l'année et un B-spline à quatre degrés de liberté de la surface du spécimen. Les effets par taille proviennent d'un seul modèle d'interaction. Cette estimation de durée totale est distincte d'une conversion de vitesse en minutes.", "", "## Calipers élargis", ""])
    for row in calipers.loc[calipers.stratum.isin(["overall", ">=50 mm", "size interaction"])].itertuples():
        paragraphs.extend([f"Caliper {row.caliper:.2f}, {row.stratum} : {row.n_pairs_model} paires, n = {row.n}, {effect_text(row)} mm²/min ; SMD maximal {row.max_abs_smd:.3f}.", ""])
    paragraphs.extend([
        "Les calipers 0,40 et 0,60 ne satisfont pas le critère de balance SMD <0,10. Ils restent des analyses exploratoires de relâchement, sans remplacer l'analyse principale ni permettre de conclure à une meilleure estimation causale.", "",
        "## Sources", "",
        "Toutes les estimations sont dans results/tables/estimates_summary.csv. Les effectifs n correspondent aux procédures de la strate ; n_model et n_pairs_model donnent la population utilisée pour estimer le modèle. Les termes d'interaction correspondent à une différence entre effets, non à un effet dans une strate supplémentaire. Aucun test n'est sélectionné selon sa significativité.",
    ])
    output.write_text("\n".join(paragraphs) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=REPO_ROOT / "data/public/atract_analysis_public.csv")
    parser.add_argument("--output", type=Path, default=REVISION_ROOT / "results")
    args = parser.parse_args()
    frozen = verify_baseline()
    dataset_hash = checksum(args.input)
    if dataset_hash != EXPECTED_PUBLIC_DATASET_SHA256:
        raise ValueError("The analytic CSV does not match the locked dataset checksum.")
    public = pd.read_csv(args.input)
    definition = pd.read_csv(REVISION_ROOT / "operator_volume_definition.csv")
    groups = definition.set_index("operator_id_public").volume_group
    if not set(public.operator_id_public.dropna()).issubset(groups.index):
        raise ValueError("An operator stratum has no attested volume classification.")
    supported = _analysis_cohort(public, required_columns=SPEC["speed_required_columns"], outcome="speed_mm2_min")
    scored, propensity = fit_propensity_model(supported, ps_rhs=SPEC["ps_rhs"], ps_structure=PRIMARY_PS_STRUCTURE, weighting_method="overlap")
    scored = add_large_lesion_flag(scored)
    matched = _nearest_neighbor_match(scored, ratio=1, caliper_multiplier=.20, match_scope="operator_only")
    fibrosis, complete_pairs = fibrosis_models(matched)
    verify_primary(fibrosis)
    overlap = overlap_fibrosis_models(scored)
    volume, volume_balance = operator_volume_models(matched, groups)
    calipers, caliper_balance = caliper_models(scored)
    curve, continuous_test = continuous_size_test(matched)
    archived_curve = pd.read_csv(REPO_ROOT / "results/tables/table_s9_speed_continuous_size.csv")
    pd.testing.assert_frame_equal(curve, archived_curve, check_exact=False, atol=1e-10, rtol=1e-10)
    duration = matched_effects(matched, "procedure_duration_min", DURATION_RHS, "matched_duration_year_area")
    secondary = pd.read_csv(REPO_ROOT / "results/model_summaries/primary_binary_results.csv")
    secondary["analysis"], secondary["stratum"], secondary["unit"] = "archived_overlap_secondary", "overall", "risk ratio"
    exclude_pooled = matched_effects(matched.loc[matched.operator_id_public.ne("operator_other")], "speed_mm2_min", SPEED_RHS, "matched_exclude_pooled_operator")
    fibrosis_balance = build_matching_balance_table(scored, complete_pairs, SPEC["balance_covariates"] + ["fibrosis"])
    ow_complete = scored.loc[scored.fibrosis.notna()]
    overlap_balance = build_balance_table(ow_complete, SPEC["balance_covariates"] + ["fibrosis"], "analysis_weight")
    tables = {
        "fibrosis_effects": fibrosis, "overlap_fibrosis_effects": overlap,
        "operator_volume_effects": volume, "operator_volume_balance": volume_balance,
        "caliper_effects": calipers, "caliper_balance": caliper_balance,
        "continuous_size_curve": curve, "continuous_size_test": continuous_test,
        "duration_effects": duration, "secondary_outcomes": secondary,
        "individual_operator_effects": operator_identity_models(matched),
        "exclude_pooled_operator_effects": exclude_pooled,
        "fibrosis_pair_balance": fibrosis_balance, "fibrosis_overlap_balance": overlap_balance,
    }
    summary = pd.concat([tables[name] for name in [
        "fibrosis_effects", "overlap_fibrosis_effects", "operator_volume_effects", "caliper_effects",
        "continuous_size_test", "duration_effects", "secondary_outcomes", "individual_operator_effects", "exclude_pooled_operator_effects",
    ]], ignore_index=True)
    columns = ["analysis", "stratum", "outcome", "n", "n_model", "n_pairs_model", "treated_n", "control_n", "estimate", "ci_lower", "ci_upper", "p_value", "unit", "scale", "caliper", "max_abs_smd", "balance_below_0_10", "effective_sample_size", "df", "formula"]
    tables["estimates_summary"] = summary.reindex(columns=columns)
    metadata = {
        "version": "v1.1.0", "analysis_date": "2026-10-13", "baseline_commit": BASELINE_COMMIT,
        "analytic_dataset_sha256": dataset_hash, "operator_definition_sha256": checksum(REVISION_ROOT / "operator_volume_definition.csv"),
        "script_sha256": checksum(Path(__file__)), "helpers_sha256": checksum(REVISION_ROOT / "analysis_helpers.py"),
        "software": {"python": platform.python_version(), **{name: version(name) for name in ["numpy", "pandas", "scipy", "statsmodels", "patsy", "matplotlib"]}},
        "propensity_formula": propensity.model.formula, "primary_pairs": matched._match_group.nunique(),
        "fibrosis_complete_pairs": complete_pairs._match_group.nunique(),
        "continuous_curve_matches_v1_0_0": True,
        "continuous_curve_sha256_v1_0_0": checksum(REPO_ROOT / "results/tables/table_s9_speed_continuous_size.csv"),
        "continuous_interaction": continuous_test.iloc[0].to_dict(),
        "all_v1_0_0_outputs_verified": True, "frozen_output_count": len(frozen),
        "operator_classification_source": "clinical study-team attestation dated 2026-10-13; fixed group, not time-varying experience",
        "inference": "matched-set clustered sandwich for matched models; archived empirical influence variance for overlap augmentation",
        "matching": "global minimum-distance assignment within pseudonymized operator strata; logit calipers in SD units",
        "wider_calipers_are_exploratory": True,
    }
    for name in ["tables", "figures", "metadata"]:
        (args.output / name).mkdir(parents=True, exist_ok=True)
    for name, frame in tables.items():
        frame.to_csv(args.output / "tables" / f"{name}.csv", index=False)
    (args.output / "metadata/analysis_metadata.json").write_text(json.dumps(metadata, indent=2, default=lambda x: x.item() if hasattr(x, "item") else str(x)) + "\n")
    figures(curve, continuous_test, volume, args.output / "figures")
    write_response(tables, args.output.parent / "RESULTS_fr.md")
    verify_baseline()
    if checksum(args.input) != dataset_hash:
        raise ValueError("The analytic input changed during analysis.")
    print(f"v1.1.0: {len(summary)} estimates exported; {len(frozen)} archived outputs unchanged.")


if __name__ == "__main__":
    main()
