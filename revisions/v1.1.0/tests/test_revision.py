from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
TABLES = ROOT / "results/tables"
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(ROOT))


class RevisionTests(unittest.TestCase):
    def test_archived_results_are_bitwise_unchanged(self):
        manifest = json.loads((ROOT / "baseline_output_checksums.json").read_text())
        self.assertEqual(len(manifest), 61)
        for relative, expected in manifest.items():
            actual = hashlib.sha256((REPO / relative).read_bytes()).hexdigest()
            self.assertEqual(actual, expected, relative)

    def test_primary_estimates_reproduce_v1_0_0(self):
        from analysis_helpers import STRATUM_LABELS
        archived = pd.read_csv(REPO / "results/model_summaries/primary_speed_results.csv")
        table = pd.read_csv(TABLES / "fibrosis_effects.csv").query("analysis == 'matched_primary'").set_index("stratum")
        for column in ["n", "estimate", "ci_lower", "ci_upper", "p_value"]:
            np.testing.assert_allclose(table.loc[archived.analysis.map(STRATUM_LABELS), column], archived[column], rtol=1e-10, atol=1e-10)

    def test_fibrosis_comparisons_use_the_same_complete_pairs(self):
        table = pd.read_csv(TABLES / "fibrosis_effects.csv").set_index(["analysis", "stratum"])
        reference = table.loc["fibrosis_complete_pairs_year_size"]
        self.assertEqual(reference.loc["overall", "n_pairs_model"], 372)
        for name in ["fibrosis_complete_pairs_year_size_fibrosis", "fibrosis_complete_pairs_full_covariates", "fibrosis_complete_pairs_full_covariates_fibrosis"]:
            for column in ["n", "n_model", "n_pairs_model", "treated_n", "control_n"]:
                np.testing.assert_array_equal(reference[column], table.loc[name, column])

    def test_interaction_contrasts_are_consistent(self):
        for filename in ["fibrosis_effects", "overlap_fibrosis_effects", "caliper_effects", "duration_effects"]:
            table = pd.read_csv(TABLES / f"{filename}.csv")
            for _, frame in table.groupby("analysis"):
                effects = frame.set_index("stratum").estimate
                self.assertAlmostEqual(effects[">=50 mm"], effects["<50 mm"] + effects["size interaction"], places=10)
        volume = pd.read_csv(TABLES / "operator_volume_effects.csv")
        for _, frame in volume.groupby("analysis"):
            effects = frame.set_index("stratum").estimate
            self.assertAlmostEqual(effects["high-volume operators"], effects["lower-volume trained operators"] + effects["operator-volume interaction"], places=10)

    def test_continuous_curve_and_wald_degrees_of_freedom(self):
        archived = pd.read_csv(REPO / "results/tables/table_s9_speed_continuous_size.csv")
        curve = pd.read_csv(TABLES / "continuous_size_curve.csv")
        pd.testing.assert_frame_equal(curve, archived, atol=1e-10, rtol=1e-10)
        joint = pd.read_csv(TABLES / "continuous_size_test.csv").iloc[0]
        self.assertEqual(joint.df, 4)
        self.assertEqual(joint.n, 782)
        self.assertEqual(joint.test, "joint robust Wald chi-square")

    def test_operator_volume_definition_and_counts(self):
        definition = pd.read_csv(ROOT / "operator_volume_definition.csv")
        high = set(definition.loc[definition.volume_group.eq("high_volume"), "operator_id_public"])
        self.assertEqual(high, {"operator_01", "operator_02"})
        table = pd.read_csv(TABLES / "operator_volume_effects.csv").query("analysis == 'operator_volume_full_covariates'").set_index("stratum")
        self.assertEqual(table.loc["high-volume operators", "n_pairs_stratum"], 293)
        self.assertEqual(table.loc["lower-volume trained operators", "n_pairs_stratum"], 98)
        self.assertTrue(table.formula.str.contains("C(operator_id_public)", regex=False).all())

    def test_wider_calipers_report_failed_balance(self):
        grid = pd.read_csv(REPO / "results/model_summaries/speed_matching_grid.csv").set_index("caliper_multiplier")
        table = pd.read_csv(TABLES / "caliper_effects.csv").query("stratum == 'overall'").set_index("caliper")
        for caliper in [.4, .6]:
            self.assertEqual(table.loc[caliper, "n_pairs_model"], grid.loc[caliper, "n_match_groups"])
            self.assertAlmostEqual(table.loc[caliper, "max_abs_smd"], grid.loc[caliper, "max_abs_smd"], places=10)
            self.assertFalse(table.loc[caliper, "balance_below_0_10"])

    def test_summary_intervals_and_private_column_exclusion(self):
        table = pd.read_csv(TABLES / "estimates_summary.csv")
        effects = table.loc[table.unit.ne("Wald chi-square statistic")]
        self.assertTrue(effects[["estimate", "ci_lower", "ci_upper", "p_value"]].notna().all().all())
        self.assertTrue(effects.ci_lower.le(effects.estimate).all())
        self.assertTrue(effects.estimate.le(effects.ci_upper).all())
        self.assertTrue(table.p_value.between(0, 1).all())
        banned = {"ID", "NOM", "PRENOM", "Date examen", "Physician", "_analysis_row_id", "_match_group"}
        for path in TABLES.glob("*.csv"):
            self.assertFalse(set(pd.read_csv(path, nrows=0).columns) & banned, path.name)

    def test_figures_have_all_three_formats(self):
        for stem in ["continuous_size_interaction", "operator_volume_effects"]:
            for extension, signature in [("svg", b"<?xml"), ("png", b"\x89PNG"), ("pdf", b"%PDF")]:
                self.assertTrue((ROOT / "results/figures" / f"{stem}.{extension}").read_bytes().startswith(signature))

    @unittest.skipUnless((REPO / "data/public/atract_analysis_public.csv").exists(), "Authorized row-level data are not publicly distributed.")
    def test_matching_reproduces_from_authorized_data(self):
        from atract_analysis.cohort import add_large_lesion_flag
        from atract_analysis.config import PRIMARY_PS_STRUCTURE
        from atract_analysis.models import _analysis_cohort, _nearest_neighbor_match, fit_propensity_model
        from analysis_helpers import SPEC, fibrosis_models
        public = pd.read_csv(REPO / "data/public/atract_analysis_public.csv")
        supported = _analysis_cohort(public, required_columns=SPEC["speed_required_columns"], outcome="speed_mm2_min")
        scored, _ = fit_propensity_model(supported, ps_rhs=SPEC["ps_rhs"], ps_structure=PRIMARY_PS_STRUCTURE, weighting_method="overlap")
        matched = _nearest_neighbor_match(add_large_lesion_flag(scored), ratio=1, caliper_multiplier=.20, match_scope="operator_only")
        recomputed, _ = fibrosis_models(matched)
        saved = pd.read_csv(TABLES / "fibrosis_effects.csv")
        for column in ["n", "estimate", "ci_lower", "ci_upper", "p_value"]:
            np.testing.assert_allclose(recomputed[column], saved[column], rtol=1e-10, atol=1e-10)


if __name__ == "__main__":
    unittest.main()
