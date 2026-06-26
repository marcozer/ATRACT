from __future__ import annotations

import unittest

import pandas as pd

from atract_analysis.anonymize import build_data_dictionary
from atract_analysis.checks import validate_public_dataset
from atract_analysis.config import CAUSAL_SPECS, PRIMARY_CAUSAL_SPEC, PUBLIC_COLUMNS, REPO_ROOT
from atract_analysis.models import (
    run_binary_comparison_analyses,
    run_primary_binary_analyses,
    run_primary_speed_analysis,
    run_speed_robustness_analysis,
)
from atract_analysis.tables import (
    build_missingness_by_group_table,
    build_operator_adoption_table,
    build_operator_year_distribution_table,
    build_unmatched_profile_table,
)


PUBLIC_DATA_PATH = REPO_ROOT / "data" / "public" / "atract_analysis_public.csv"


class PublicReleaseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.public_df = pd.read_csv(PUBLIC_DATA_PATH) if PUBLIC_DATA_PATH.exists() else None

    def require_public_data(self) -> pd.DataFrame:
        if self.public_df is None:
            self.skipTest("Authorized analytic dataset is not present.")
        return self.public_df

    def test_public_schema_is_locked(self) -> None:
        public_df = self.require_public_data()
        self.assertEqual(list(public_df.columns), PUBLIC_COLUMNS)

    def test_public_counts_are_locked(self) -> None:
        public_df = self.require_public_data()
        self.assertEqual(len(public_df), 2003)
        counts = public_df["atract"].value_counts().to_dict()
        self.assertEqual(counts.get(0), 1503)
        self.assertEqual(counts.get(1), 500)

    def test_jnet_iii_is_excluded(self) -> None:
        public_df = self.require_public_data()
        self.assertNotIn("jnet_iii", set(public_df["jnet_group"].dropna()))

    def test_macronodule_is_applicable_only_to_lst_granular(self) -> None:
        public_df = self.require_public_data()
        non_granular = public_df.loc[~public_df["lesion_type"].eq("lst_granular")]
        self.assertTrue(non_granular["macronodule"].isna().all())

    def test_release_checks_pass(self) -> None:
        public_df = self.require_public_data()
        errors = validate_public_dataset(public_df)
        self.assertEqual(errors, [])

    def test_dictionary_matches_schema(self) -> None:
        dictionary = build_data_dictionary()
        self.assertEqual(dictionary["column_name"].tolist(), PUBLIC_COLUMNS)

    def test_primary_spec_uses_preoperative_fibrosis_risk_variables(self) -> None:
        spec = CAUSAL_SPECS[PRIMARY_CAUSAL_SPEC]
        formula_text = " ".join(str(value) for value in spec.values())
        self.assertIn("jnet_group", formula_text)
        self.assertIn("lesion_morphology", formula_text)
        self.assertIn("mici_history", formula_text)
        self.assertIn("recurrence_history", formula_text)
        self.assertNotIn("conecct_group", formula_text)
        self.assertNotIn("C(fibrosis)", formula_text)

    def test_analysis_pipeline_generates_expected_frames(self) -> None:
        public_df = self.require_public_data()
        primary_speed = run_primary_speed_analysis(public_df, bootstrap_iterations=5)
        speed_robustness = run_speed_robustness_analysis(public_df, primary_speed)
        primary_binary = run_primary_binary_analyses(public_df)
        binary_comparison = run_binary_comparison_analyses(
            public_df,
            caliper_multiplier=float(primary_speed["metadata"]["caliper_multiplier"]),
        )
        self.assertIn("matched_frame", primary_speed)
        self.assertIn("scored_frame", primary_speed)
        self.assertIn("_analysis_row_id", primary_speed["scored_frame"].columns)
        self.assertIn("_analysis_row_id", primary_speed["matched_frame"].columns)
        self.assertGreater(len(primary_speed["matched_frame"]), 0)
        self.assertGreater(len(primary_speed["matching_grid"]), 0)
        self.assertEqual(primary_speed["metadata"]["support_scope"], "operator_only")
        self.assertEqual(primary_speed["metadata"]["ps_structure"], "operator_only_plus_year_numeric")
        self.assertEqual(primary_speed["metadata"]["n_match_groups"], 391)
        self.assertEqual(primary_speed["metadata"]["n_supported_cohort"], 1916)
        primary_overall = primary_speed["effects"].loc[primary_speed["effects"]["analysis"].eq("overall")].iloc[0]
        self.assertAlmostEqual(primary_overall["estimate"], 2.9986546485814127)
        self.assertIn("summary", speed_robustness)
        self.assertEqual(set(speed_robustness["summary"]["method_label"]), {"PS-NN", "OW + DR"})
        self.assertIn("large_lesion_>=50mm_interaction_model", set(speed_robustness["summary"]["analysis"]))
        self.assertIn("large_lesion_>=50mm_restricted", set(speed_robustness["summary"]["analysis"]))
        self.assertIn("aipw_decomposition_overall", speed_robustness["metadata"])
        speed_effects = primary_speed["effects"].set_index("analysis")
        large_identity = (
            speed_effects.loc["small_lesion_<50mm", "estimate"]
            + speed_effects.loc["interaction_large_lesion", "estimate"]
        )
        self.assertAlmostEqual(speed_effects.loc["large_lesion_>=50mm", "estimate"], large_identity)
        self.assertEqual(primary_speed["metadata"]["pair_diagnostics"]["split_size_pairs"], 168)
        self.assertEqual(primary_speed["metadata"]["pair_diagnostics"]["year_gap_le_1_pairs"], 332)
        self.assertEqual(set(primary_speed["temporal_sensitivity"]["summary"]["sensitivity"]), {
            "existing_pairs_year_gap_le_1",
            "rematched_post_adoption_period",
            "rematched_year_gap_le_1",
        })
        self.assertEqual(len(primary_speed["bootstrap"]), 4)
        self.assertGreater(len(primary_speed["continuous_size"]), 0)
        self.assertEqual(set(primary_binary["effects"]["outcome"]), {"r0", "perforation", "delayed_bleeding"})
        self.assertEqual(set(primary_binary["effects"]["term"]), {"aipw_overlap"})
        self.assertIn("risk_difference", primary_binary["effects"].columns)
        self.assertIn("effective_sample_size", primary_binary["diagnostics"]["r0"])
        self.assertGreater(len(primary_binary["balance"]), 0)
        self.assertEqual(set(binary_comparison["summary"]["outcome"]), {"r0", "perforation", "delayed_bleeding"})
        self.assertEqual(set(binary_comparison["summary"]["method_label"]), {"OW + DR", "PS-NN"})
        self.assertLess(primary_speed["metadata"]["max_abs_smd"], 0.10)
        self.assertIn("small_lesion_<50mm", set(primary_speed["effects"]["analysis"]))
        self.assertIn("large_lesion_>=50mm", set(primary_speed["effects"]["analysis"]))

        unmatched_profile = build_unmatched_profile_table(primary_speed["scored_frame"], primary_speed["matched_frame"])
        unmatched_n = unmatched_profile.loc[unmatched_profile["variable"].eq("n"), "overall"].iloc[0]
        self.assertEqual(unmatched_n, "1134")

        operator_year = build_operator_year_distribution_table(public_df)
        self.assertIn("atract_pct", operator_year.columns)
        self.assertIn("speed_atract_median", operator_year.columns)

        adoption = build_operator_adoption_table(public_df)
        self.assertIn("operator_category_note", adoption.columns)
        self.assertIn("operator_other", set(adoption["operator_id_public"]))

        missingness_by_group = build_missingness_by_group_table(public_df)
        self.assertEqual(set(["overall", "conventional", "atract"]).issubset(missingness_by_group.columns), True)


if __name__ == "__main__":
    unittest.main()
