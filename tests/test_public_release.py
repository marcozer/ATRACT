from __future__ import annotations

import unittest

import pandas as pd

from atract_analysis.anonymize import build_data_dictionary
from atract_analysis.checks import validate_public_dataset
from atract_analysis.config import PUBLIC_COLUMNS, REPO_ROOT
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
        self.assertEqual(len(self.public_df), 2189)
        counts = self.public_df["atract"].value_counts().to_dict()
        self.assertEqual(counts.get(0), 1617)
        self.assertEqual(counts.get(1), 572)

    def test_release_checks_pass(self) -> None:
        errors = validate_public_dataset(self.public_df)
        self.assertEqual(errors, [])

    def test_dictionary_matches_schema(self) -> None:
        dictionary = build_data_dictionary()
        self.assertEqual(dictionary["column_name"].tolist(), PUBLIC_COLUMNS)

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
        self.assertIn("summary", speed_robustness)
        self.assertEqual(set(speed_robustness["summary"]["method_label"]), {"PS-NN", "OW + DR"})
        self.assertEqual(set(primary_binary["effects"]["outcome"]), {"r0", "perforation", "delayed_bleeding"})
        self.assertEqual(set(binary_comparison["summary"]["outcome"]), {"r0", "perforation", "delayed_bleeding"})
        self.assertEqual(set(binary_comparison["summary"]["method_label"]), {"OW + DR", "PS-NN"})
        self.assertLess(primary_speed["metadata"]["max_abs_smd"], 0.07)


if __name__ == "__main__":
    unittest.main()
