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


PUBLIC_DATA_PATH = REPO_ROOT / "data" / "public" / "atract_analysis_public.csv"


class PublicReleaseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.public_df = pd.read_csv(PUBLIC_DATA_PATH)

    def test_public_schema_is_locked(self) -> None:
        self.assertEqual(list(self.public_df.columns), PUBLIC_COLUMNS)

    def test_public_counts_are_locked(self) -> None:
        self.assertEqual(len(self.public_df), 2003)
        counts = self.public_df["atract"].value_counts().to_dict()
        self.assertEqual(counts.get(0), 1503)
        self.assertEqual(counts.get(1), 500)

    def test_jnet_iii_is_excluded(self) -> None:
        self.assertNotIn("jnet_iii", set(self.public_df["jnet_group"].dropna()))

    def test_macronodule_is_applicable_only_to_lst_granular(self) -> None:
        non_granular = self.public_df.loc[~self.public_df["lesion_type"].eq("lst_granular")]
        self.assertTrue(non_granular["macronodule"].isna().all())

    def test_release_checks_pass(self) -> None:
        errors = validate_public_dataset(self.public_df)
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
        primary_speed = run_primary_speed_analysis(self.public_df)
        speed_robustness = run_speed_robustness_analysis(self.public_df, primary_speed)
        primary_binary = run_primary_binary_analyses(self.public_df)
        binary_comparison = run_binary_comparison_analyses(
            self.public_df,
            caliper_multiplier=float(primary_speed["metadata"]["caliper_multiplier"]),
        )
        self.assertIn("matched_frame", primary_speed)
        self.assertIn("scored_frame", primary_speed)
        self.assertGreater(len(primary_speed["matched_frame"]), 0)
        self.assertGreater(len(primary_speed["matching_grid"]), 0)
        self.assertEqual(primary_speed["metadata"]["support_scope"], "operator_only")
        self.assertEqual(primary_speed["metadata"]["ps_structure"], "operator_only_plus_year_numeric")
        self.assertEqual(primary_speed["metadata"]["n_match_groups"], 391)
        primary_overall = primary_speed["effects"].loc[primary_speed["effects"]["analysis"].eq("overall")].iloc[0]
        self.assertLess(primary_overall["p_value"], 0.05)
        self.assertIn("summary", speed_robustness)
        self.assertEqual(set(speed_robustness["summary"]["method_label"]), {"PS-NN", "OW + DR"})
        self.assertEqual(set(primary_binary["effects"]["outcome"]), {"r0", "perforation", "delayed_bleeding"})
        self.assertEqual(set(binary_comparison["summary"]["outcome"]), {"r0", "perforation", "delayed_bleeding"})
        self.assertEqual(set(binary_comparison["summary"]["method_label"]), {"OW + DR", "PS-NN"})
        self.assertLess(primary_speed["metadata"]["max_abs_smd"], 0.10)
        self.assertIn("small_lesion_<50mm", set(primary_speed["effects"]["analysis"]))
        self.assertIn("large_lesion_>=50mm", set(primary_speed["effects"]["analysis"]))


if __name__ == "__main__":
    unittest.main()
