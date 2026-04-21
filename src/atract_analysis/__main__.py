from __future__ import annotations

import argparse
from pathlib import Path

from .anonymize import write_public_release
from .checks import run_release_checks
from .cohort import load_public_dataset
from .config import REPO_ROOT
from .models import (
    run_binary_comparison_analyses,
    run_primary_binary_analyses,
    run_primary_speed_analysis,
    run_speed_robustness_analysis,
    write_model_outputs,
)
from .tables import write_tables


def _default_public_path() -> Path:
    return REPO_ROOT / "data" / "public" / "atract_analysis_public.csv"


def run_analysis(public_path: Path) -> None:
    from .figures import write_figures

    run_release_checks(public_path)
    public_dataframe = load_public_dataset(str(public_path))
    primary_speed = run_primary_speed_analysis(public_dataframe)
    speed_robustness = run_speed_robustness_analysis(public_dataframe, primary_speed)
    primary_binary = run_primary_binary_analyses(public_dataframe)
    binary_comparison = run_binary_comparison_analyses(
        public_dataframe,
        caliper_multiplier=float(primary_speed["metadata"]["caliper_multiplier"]),
    )

    write_tables(
        REPO_ROOT / "results" / "tables",
        public_dataframe=public_dataframe,
        primary_speed=primary_speed,
        speed_robustness=speed_robustness,
        primary_binary=primary_binary,
        binary_comparison=binary_comparison,
    )
    write_figures(
        REPO_ROOT / "results" / "figures",
        public_dataframe=public_dataframe,
        primary_speed=primary_speed,
        speed_robustness=speed_robustness,
        primary_binary=primary_binary,
    )
    write_model_outputs(
        REPO_ROOT / "results" / "model_summaries",
        primary_speed,
        speed_robustness,
        primary_binary,
        binary_comparison,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="ATRACT public statistical analysis")
    subparsers = parser.add_subparsers(dest="command", required=True)

    public_data = subparsers.add_parser("public-data", help="Generate the public deidentified dataset from the raw workbook")
    public_data.add_argument("--raw", required=True, type=Path)
    public_data.add_argument("--output", default=_default_public_path(), type=Path)
    public_data.add_argument("--dictionary", default=REPO_ROOT / "metadata" / "public_data_dictionary.csv", type=Path)

    analysis = subparsers.add_parser("analysis", help="Run the public analysis from the committed public dataset")
    analysis.add_argument("--input", default=_default_public_path(), type=Path)

    check = subparsers.add_parser("check", help="Run release checks against the public dataset")
    check.add_argument("--input", default=_default_public_path(), type=Path)

    all_cmd = subparsers.add_parser("all", help="Regenerate the public dataset if requested, then run analysis")
    all_cmd.add_argument("--raw", type=Path, default=None)
    all_cmd.add_argument("--input", default=_default_public_path(), type=Path)

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "public-data":
        write_public_release(args.raw, args.output, args.dictionary)
        return

    if args.command == "analysis":
        run_analysis(args.input)
        return

    if args.command == "check":
        run_release_checks(args.input)
        return

    if args.command == "all":
        if args.raw is not None:
            write_public_release(args.raw, args.input, REPO_ROOT / "metadata" / "public_data_dictionary.csv")
        run_analysis(args.input)
        return
