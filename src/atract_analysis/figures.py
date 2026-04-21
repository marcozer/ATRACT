from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

NORD = {
    "bg": "#ECEFF4",
    "fg": "#2E3440",
    "muted": "#4C566A",
    "grid": "#D8DEE9",
    "frost_1": "#8FBCBB",
    "frost_2": "#88C0D0",
    "frost_3": "#81A1C1",
    "frost_4": "#5E81AC",
    "red": "#BF616A",
    "orange": "#D08770",
    "yellow": "#EBCB8B",
    "green": "#A3BE8C",
}

plt.rcParams.update(
    {
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.edgecolor": NORD["muted"],
        "axes.labelcolor": NORD["fg"],
        "axes.titlecolor": NORD["fg"],
        "xtick.color": NORD["muted"],
        "ytick.color": NORD["muted"],
        "grid.color": NORD["grid"],
        "grid.linewidth": 0.8,
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "legend.fontsize": 9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)


def _finalize_axis(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(NORD["muted"])
    ax.spines["bottom"].set_color(NORD["muted"])
    ax.tick_params(length=3.5, width=0.8)


def _save_figure(fig, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".png", ".pdf", ".svg"):
        fig.savefig(output_path.with_suffix(suffix), dpi=300, bbox_inches="tight")


def plot_cohort_flow(flow_counts: dict[str, int], output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.axis("off")

    boxes = [
        ("Broad cleaned traction cohort", flow_counts["public_cohort"], (0.50, 0.88)),
        ("Speed analysis cohort after complete-case and support restriction", flow_counts["speed_supported"], (0.50, 0.62)),
        ("Primary PS-NN matched speed cohort", flow_counts["speed_matched"], (0.50, 0.36)),
        ("Binary OW+DR cohorts\nR0 / perforation / delayed bleeding", None, (0.50, 0.10)),
    ]

    for label, count, (x_pos, y_pos) in boxes:
        text = label if count is None else f"{label}\n n={count}"
        bbox = dict(boxstyle="round,pad=0.4", facecolor=NORD["bg"], edgecolor=NORD["frost_4"], linewidth=1.2)
        ax.text(x_pos, y_pos, text, ha="center", va="center", bbox=bbox, transform=ax.transAxes)

    ax.text(
        0.50,
        0.04,
        (
            f"R0 n={flow_counts['r0_supported']}   "
            f"Perforation n={flow_counts['perforation_supported']}   "
            f"Delayed bleeding n={flow_counts['bleeding_supported']}"
        ),
        ha="center",
        va="center",
        transform=ax.transAxes,
    )

    arrow_style = dict(arrowstyle="->", color=NORD["muted"], linewidth=1.3)
    ax.annotate("", xy=(0.50, 0.70), xytext=(0.50, 0.77), xycoords=ax.transAxes, textcoords=ax.transAxes, arrowprops=arrow_style)
    ax.annotate("", xy=(0.50, 0.44), xytext=(0.50, 0.53), xycoords=ax.transAxes, textcoords=ax.transAxes, arrowprops=arrow_style)
    ax.annotate("", xy=(0.50, 0.18), xytext=(0.50, 0.27), xycoords=ax.transAxes, textcoords=ax.transAxes, arrowprops=arrow_style)

    fig.tight_layout()
    _save_figure(fig, output_path)
    plt.close(fig)


def plot_speed_effects(primary_speed: pd.DataFrame, speed_robustness: pd.DataFrame, output_path: Path) -> None:
    primary_plot = primary_speed.loc[primary_speed["analysis"].isin(["overall", "large_lesion_>=50mm"])].copy()
    primary_plot["label"] = primary_plot["analysis"].map(
        {
            "overall": "PS-NN overall",
            "large_lesion_>=50mm": "PS-NN lesions ≥50 mm",
        }
    )
    robustness_plot = speed_robustness.copy()
    robustness_plot["label"] = "OW + DR overall"

    plot_data = pd.concat([primary_plot, robustness_plot], ignore_index=True)
    fig, ax = plt.subplots(figsize=(7, 4))
    y_positions = np.arange(len(plot_data))
    ax.errorbar(
        plot_data["estimate"],
        y_positions,
        xerr=[
            plot_data["estimate"] - plot_data["ci_lower"],
            plot_data["ci_upper"] - plot_data["estimate"],
        ],
        fmt="o",
        color=NORD["frost_4"],
        ecolor=NORD["frost_3"],
        elinewidth=1.4,
        capsize=3,
        markersize=6,
    )
    ax.axvline(0, color=NORD["muted"], linestyle="--", linewidth=1)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(plot_data["label"])
    ax.set_xlabel("Adjusted mean difference in speed (mm²/min)")
    ax.set_title("Primary and robustness estimates for dissection speed")
    ax.grid(axis="x", alpha=0.6)
    _finalize_axis(ax)
    fig.tight_layout()
    _save_figure(fig, output_path)
    plt.close(fig)


def plot_love(balance_table: pd.DataFrame, output_path: Path) -> None:
    plot_data = balance_table.copy()
    plot_data["adjusted_abs"] = plot_data["adjusted_smd"].abs()
    plot_data = plot_data.sort_values("adjusted_abs", ascending=True).tail(15)

    fig, ax = plt.subplots(figsize=(8, 6))
    y_positions = np.arange(len(plot_data))
    ax.scatter(plot_data["raw_smd"], y_positions, label="Raw", color=NORD["red"], s=34)
    ax.scatter(plot_data["adjusted_smd"], y_positions, label="Adjusted", color=NORD["frost_4"], s=34)
    ax.axvline(0.1, color=NORD["muted"], linestyle="--", linewidth=1)
    ax.axvline(-0.1, color=NORD["muted"], linestyle="--", linewidth=1)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(plot_data["term"])
    ax.set_xlabel("Standardized mean difference")
    ax.set_title("Covariate balance before and after matching")
    ax.legend()
    ax.grid(axis="x", alpha=0.6)
    _finalize_axis(ax)
    fig.tight_layout()
    _save_figure(fig, output_path)
    plt.close(fig)


def plot_propensity_overlap(scored_frame: pd.DataFrame, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(
        scored_frame.loc[scored_frame["atract"].eq(0), "propensity_score"],
        bins=20,
        alpha=0.6,
        label="Non-ATRACT",
        color=NORD["frost_2"],
        edgecolor="white",
    )
    ax.hist(
        scored_frame.loc[scored_frame["atract"].eq(1), "propensity_score"],
        bins=20,
        alpha=0.6,
        label="ATRACT",
        color=NORD["frost_4"],
        edgecolor="white",
    )
    ax.set_xlabel("Propensity score")
    ax.set_ylabel("Count")
    ax.set_title("Propensity score overlap in the supported speed cohort")
    ax.legend()
    ax.grid(axis="y", alpha=0.6)
    _finalize_axis(ax)
    fig.tight_layout()
    _save_figure(fig, output_path)
    plt.close(fig)


def plot_dag(output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.axis("off")

    positions = {
        "Patient factors": (0.08, 0.82),
        "Lesion factors": (0.32, 0.82),
        "Calendar time": (0.55, 0.82),
        "Operator": (0.78, 0.82),
        "ATRACT": (0.45, 0.50),
        "Dissection speed": (0.22, 0.18),
        "R0 resection": (0.45, 0.18),
        "Perforation": (0.66, 0.18),
        "Delayed bleeding": (0.86, 0.18),
    }

    for label, (x_pos, y_pos) in positions.items():
        bbox = dict(boxstyle="round,pad=0.3", facecolor=NORD["bg"], edgecolor=NORD["frost_4"], linewidth=1.1)
        ax.text(x_pos, y_pos, label, ha="center", va="center", bbox=bbox, transform=ax.transAxes)

    edges = [
        ("Lesion factors", "ATRACT"),
        ("Calendar time", "ATRACT"),
        ("Operator", "ATRACT"),
        ("Patient factors", "Dissection speed"),
        ("Lesion factors", "Dissection speed"),
        ("Operator", "Dissection speed"),
        ("ATRACT", "Dissection speed"),
        ("Lesion factors", "R0 resection"),
        ("Operator", "R0 resection"),
        ("ATRACT", "R0 resection"),
        ("Patient factors", "Perforation"),
        ("Lesion factors", "Perforation"),
        ("Operator", "Perforation"),
        ("ATRACT", "Perforation"),
        ("Patient factors", "Delayed bleeding"),
        ("Lesion factors", "Delayed bleeding"),
        ("Operator", "Delayed bleeding"),
        ("ATRACT", "Delayed bleeding"),
    ]

    for start, end in edges:
        x1, y1 = positions[start]
        x2, y2 = positions[end]
        ax.annotate(
            "",
            xy=(x2, y2),
            xytext=(x1, y1),
            xycoords=ax.transAxes,
            textcoords=ax.transAxes,
            arrowprops=dict(arrowstyle="->", color=NORD["muted"], linewidth=1.2),
        )

    fig.tight_layout()
    _save_figure(fig, output_path)
    plt.close(fig)


def write_figures(
    figures_dir: Path,
    public_dataframe: pd.DataFrame,
    primary_speed: dict[str, object],
    speed_robustness: dict[str, object],
    primary_binary: dict[str, object],
) -> None:
    figures_dir.mkdir(parents=True, exist_ok=True)
    stale_files = [
        "dag_structure.png",
        "binary_effects.png",
        "love_plot.png",
        "propensity_overlap.png",
        "speed_effects.png",
        "speed_boxplot.png",
    ]
    for filename in stale_files:
        stale_path = figures_dir / filename
        if stale_path.exists():
            stale_path.unlink()

    flow_counts = {
        "public_cohort": int(len(public_dataframe)),
        "speed_supported": int(len(primary_speed["scored_frame"])),
        "speed_matched": int(len(primary_speed["matched_frame"])),
        "r0_supported": int(primary_binary["diagnostics"]["r0"]["n_complete_case"]),
        "perforation_supported": int(primary_binary["diagnostics"]["perforation"]["n_complete_case"]),
        "bleeding_supported": int(primary_binary["diagnostics"]["delayed_bleeding"]["n_complete_case"]),
    }

    plot_cohort_flow(flow_counts, figures_dir / "figure_1_cohort_flow.png")
    plot_speed_effects(primary_speed["effects"], speed_robustness["effects"], figures_dir / "figure_2_speed_effects.png")
    plot_love(primary_speed["balance"], figures_dir / "figure_s1_love_plot.png")
    plot_propensity_overlap(primary_speed["scored_frame"], figures_dir / "figure_s2_propensity_overlap.png")
    plot_dag(figures_dir / "figure_s3_dag.png")
