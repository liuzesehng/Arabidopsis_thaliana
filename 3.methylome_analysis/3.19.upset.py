#!/usr/bin/env python3
"""Build a publication-ready composite UpSet figure for AT2G39730 features."""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Sequence

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rca-upset-matplotlib-cache")
)

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib import colors as mpl_colors
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Circle, Rectangle
from matplotlib.transforms import ScaledTranslation
import numpy as np
import pandas as pd


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[2]
DEFAULT_RCA_DIR = PROJECT_ROOT / "list" / "RCA"
DEFAULT_INPUT_DIR = DEFAULT_RCA_DIR / "AT2G39730_shap_location_outputs"
DEFAULT_OUTPUT_DIR = DEFAULT_RCA_DIR / "upset"
DEFAULT_BASENAME = "AT2G39730_top10_methylation_upset_composite"

DATASETS = ("A_TPM", "A_per_TPM", "B_TPM", "B_per_TPM", "total_TPM")
DATASET_LABELS = {
    "A_TPM": "α",
    "A_per_TPM": "α%",
    "B_TPM": "β",
    "B_per_TPM": "β%",
    "total_TPM": "total",
}
CONTEXTS = ("CHH", "CHG", "CG")
ALL_CONTEXTS = (*CONTEXTS, "snp")
REGIONS = ("upstream", "gene_body", "downstream")
REGION_LABELS = {
    "upstream": "Upstream",
    "gene_body": "Gene body",
    "downstream": "Downstream",
}
TOP_PERCENT = 0.10

PANEL_SPECS = (
    ("a", "upstream", "CHH", "Upstream CHH"),
    ("b", "upstream", "CHG", "Upstream CHG"),
    ("c", "upstream", "CG", "Upstream CG"),
    ("d", "downstream", "CG", "Downstream CG"),
)

REQUIRED_COLUMNS = (
    "dataset",
    "context",
    "feature",
    "chrom",
    "start",
    "end",
    "plot_pos",
    "shap_value",
    "region_class",
)

CONTEXT_COLORS = {
    "CHH": "#C65A5A",
    "CHG": "#4E79A7",
    "CG": "#59A14F",
}
INK = "#202124"
MUTED_INK = "#5F6368"
INACTIVE = "#D5D8DC"
ROW_SHADE = "#F4F5F6"
GRID_COLOR = "#D9DDE1"

FIGURE_WIDTH_MM = 183.0
FIGURE_HEIGHT_MM = 170.0
FIGURE_WIDTH_IN = FIGURE_WIDTH_MM / 25.4
FIGURE_HEIGHT_IN = FIGURE_HEIGHT_MM / 25.4
DEFAULT_DPI = 600


# Optional external QA: keep each plotting script usable on its own.
_OPTIONAL_QA_SKIPPED = []


def _skip_optional_qa(reason, output=None):
    """Record an unavailable check without claiming that it passed."""
    if reason not in _OPTIONAL_QA_SKIPPED:
        _OPTIONAL_QA_SKIPPED.append(reason)
        print(f"[QA SKIPPED] {reason}", file=sys.stderr)
    report = {"verdict": "SKIPPED", "status": "SKIPPED", "auditable": False,
              "reason": reason}
    if output is not None:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return report


def _skip_alignment(*args, json_out=None, **kwargs):
    return _skip_optional_qa("audit_panel_alignment.py is unavailable", json_out)


def configure_style() -> None:
    """Apply compact publication typography with editable vector text."""
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
            "font.size": 6.2,
            "axes.titlesize": 7.2,
            "axes.labelsize": 6.4,
            "axes.linewidth": 0.65,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.labelsize": 5.8,
            "ytick.labelsize": 5.8,
            "xtick.major.size": 2.2,
            "ytick.major.size": 2.2,
            "xtick.major.width": 0.55,
            "ytick.major.width": 0.55,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
        }
    )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot a four-panel composite showing UpSet "
            "membership of top-decile AT2G39730 methylation features."
        )
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help=f"Directory containing five feature-detail CSV files (default: {DEFAULT_INPUT_DIR})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR})",
    )
    parser.add_argument(
        "--basename",
        default=DEFAULT_BASENAME,
        help=f"Output filename stem (default: {DEFAULT_BASENAME})",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=DEFAULT_DPI,
        help="Raster export resolution; must be at least 300 (default: 600)",
    )
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> None:
    if not args.input_dir.is_dir():
        raise FileNotFoundError(f"input directory not found: {args.input_dir}")
    if args.dpi < 300:
        raise ValueError("--dpi must be at least 300")
    if not args.basename or Path(args.basename).name != args.basename:
        raise ValueError("--basename must be a non-empty filename stem, not a path")


def load_dataset_table(input_dir: Path, dataset: str) -> pd.DataFrame:
    """Read one pre-annotated feature table and enforce its data contract."""
    path = input_dir / f"{dataset}_feature_details.csv"
    if not path.is_file():
        raise FileNotFoundError(f"feature-detail CSV not found: {path}")

    table = pd.read_csv(path, encoding="utf-8")
    missing = [column for column in REQUIRED_COLUMNS if column not in table.columns]
    if missing:
        raise ValueError(f"{path.name} is missing required columns: {', '.join(missing)}")
    table = table.loc[:, REQUIRED_COLUMNS].copy()

    if table.empty:
        raise ValueError(f"{path.name} contains no feature rows")
    if table["feature"].isna().any() or (table["feature"].astype(str).str.strip() == "").any():
        raise ValueError(f"{path.name} contains a missing or empty feature identifier")
    table["feature"] = table["feature"].astype(str)
    duplicated = table.loc[table["feature"].duplicated(), "feature"].tolist()
    if duplicated:
        raise ValueError(f"{path.name} contains duplicate features: {', '.join(duplicated[:5])}")

    if table["dataset"].isna().any():
        raise ValueError(f"{path.name} contains a missing dataset value")
    observed_datasets = set(table["dataset"].astype(str))
    if observed_datasets != {dataset}:
        raise ValueError(
            f"{path.name} dataset column must contain only {dataset}; observed "
            f"{sorted(observed_datasets)}"
        )

    table["context"] = table["context"].astype(str)
    invalid_contexts = sorted(set(table["context"]) - set(ALL_CONTEXTS))
    if invalid_contexts:
        raise ValueError(f"{path.name} contains invalid contexts: {invalid_contexts}")

    table["region_class"] = table["region_class"].astype(str)
    invalid_regions = sorted(set(table["region_class"]) - set(REGIONS))
    if invalid_regions:
        raise ValueError(f"{path.name} contains invalid regions: {invalid_regions}")

    for column in ("start", "end", "plot_pos", "shap_value"):
        table[column] = pd.to_numeric(table[column], errors="raise")
    numeric = table.loc[:, ["start", "end", "plot_pos", "shap_value"]].to_numpy(float)
    if not np.isfinite(numeric).all():
        raise ValueError(f"{path.name} contains missing or non-finite numeric values")
    if (table["start"] < 1).any() or (table["end"] < table["start"]).any():
        raise ValueError(f"{path.name} contains invalid genomic coordinates")
    if (table["shap_value"] < 0).any():
        raise ValueError(f"{path.name} contains negative mean absolute SHAP values")
    if not (table["shap_value"] > 0).any():
        raise ValueError(f"{path.name} contains no positive SHAP features")

    table["chrom"] = table["chrom"].astype(str)
    if (table["chrom"].str.strip() == "").any():
        raise ValueError(f"{path.name} contains an empty chromosome value")
    table["start"] = table["start"].astype(int)
    table["end"] = table["end"].astype(int)
    return table


def validate_cross_dataset_metadata(tables: dict[str, pd.DataFrame]) -> None:
    """Require shared feature identifiers to describe the same genomic site."""
    combined = pd.concat(
        [table.assign(_source_dataset=dataset) for dataset, table in tables.items()],
        ignore_index=True,
    )
    metadata_columns = ("chrom", "start", "end", "context", "region_class")
    for column in metadata_columns:
        inconsistent = combined.groupby("feature", sort=False)[column].nunique(dropna=False)
        bad = inconsistent[inconsistent > 1]
        if not bad.empty:
            examples = ", ".join(bad.index.astype(str).tolist()[:5])
            raise ValueError(
                f"shared features have inconsistent {column} values across datasets: {examples}"
            )


def select_top_features(
    tables: dict[str, pd.DataFrame],
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame]:
    """Select exactly the top 10% positive features within each response."""
    selections: dict[str, pd.DataFrame] = {}
    summary_rows: list[dict[str, object]] = []

    for dataset in DATASETS:
        positive = tables[dataset].loc[tables[dataset]["shap_value"] > 0].copy()
        positive = positive.sort_values(
            ["shap_value", "feature"], ascending=[False, True], kind="mergesort"
        ).reset_index(drop=True)
        top_n = max(1, math.floor(len(positive) * TOP_PERCENT))
        selected = positive.iloc[:top_n].copy()
        cutoff = float(selected["shap_value"].min())
        boundary_tie_count = int(np.isclose(positive["shap_value"], cutoff, rtol=0, atol=0).sum())
        selected["selection_rank"] = np.arange(1, top_n + 1, dtype=int)
        selections[dataset] = selected
        summary_rows.append(
            {
                "dataset": dataset,
                "display_label": DATASET_LABELS[dataset],
                "positive_feature_count": int(len(positive)),
                "top_percent": TOP_PERCENT,
                "selected_feature_count": top_n,
                "cutoff_shap_value": cutoff,
                "boundary_tie_count": boundary_tie_count,
                "selected_methylation_count": int(selected["context"].isin(CONTEXTS).sum()),
                "selected_snp_count": int((selected["context"] == "snp").sum()),
            }
        )

    return selections, pd.DataFrame(summary_rows)


def build_membership_table(selections: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Create one row per feature in the union of all response selections."""
    selected_union = sorted(
        set().union(*(set(table["feature"]) for table in selections.values()))
    )
    lookups = {dataset: table.set_index("feature", drop=False) for dataset, table in selections.items()}
    records: list[dict[str, object]] = []

    for feature in selected_union:
        source_rows = [lookups[dataset].loc[feature] for dataset in DATASETS if feature in lookups[dataset].index]
        source = source_rows[0]
        record: dict[str, object] = {
            "feature": feature,
            "context": str(source["context"]),
            "chrom": str(source["chrom"]),
            "start": int(source["start"]),
            "end": int(source["end"]),
            "plot_pos": float(source["plot_pos"]),
            "region_class": str(source["region_class"]),
            "plotted_methylation_feature": int(str(source["context"]) in CONTEXTS),
        }
        for dataset in DATASETS:
            present = feature in lookups[dataset].index
            record[f"in_{dataset}"] = int(present)
            record[f"{dataset}_shap_value"] = (
                float(lookups[dataset].loc[feature, "shap_value"]) if present else np.nan
            )
            record[f"{dataset}_selection_rank"] = (
                int(lookups[dataset].loc[feature, "selection_rank"]) if present else pd.NA
            )
        record["dataset_membership_count"] = sum(
            int(record[f"in_{dataset}"]) for dataset in DATASETS
        )
        record["membership_mask"] = "".join(
            str(record[f"in_{dataset}"]) for dataset in DATASETS
        )
        records.append(record)

    memberships = pd.DataFrame(records)
    region_rank = {region: index for index, region in enumerate(REGIONS)}
    context_rank = {context: index for index, context in enumerate(ALL_CONTEXTS)}
    memberships["_region_rank"] = memberships["region_class"].map(region_rank)
    memberships["_context_rank"] = memberships["context"].map(context_rank)
    memberships = memberships.sort_values(
        ["plotted_methylation_feature", "_region_rank", "_context_rank", "start", "feature"],
        ascending=[False, True, True, True, True],
        kind="mergesort",
    ).drop(columns=["_region_rank", "_context_rank"])
    return memberships.reset_index(drop=True)


def build_overview_table(memberships: pd.DataFrame) -> pd.DataFrame:
    plotted = memberships.loc[memberships["plotted_methylation_feature"] == 1]
    counts = plotted.groupby(["region_class", "context"], observed=False).size()
    rows = []
    for region in REGIONS:
        for context in CONTEXTS:
            rows.append(
                {
                    "region_class": region,
                    "region_label": REGION_LABELS[region],
                    "context": context,
                    "unique_feature_count": int(counts.get((region, context), 0)),
                }
            )
    return pd.DataFrame(rows)


def build_intersection_table(memberships: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    membership_columns = [f"in_{dataset}" for dataset in DATASETS]

    for panel, region, context, title in PANEL_SPECS:
        subset = memberships.loc[
            (memberships["plotted_methylation_feature"] == 1)
            & (memberships["region_class"] == region)
            & (memberships["context"] == context)
        ].copy()
        if subset.empty:
            raise ValueError(f"planned panel {panel} ({region} {context}) has no features")
        grouped = (
            subset.groupby(membership_columns, dropna=False, sort=False)
            .agg(
                feature_count=("feature", "size"),
                features=("feature", lambda values: ",".join(sorted(values))),
            )
            .reset_index()
        )
        grouped["membership_mask"] = grouped.apply(
            lambda row: "".join(str(int(row[column])) for column in membership_columns),
            axis=1,
        )
        grouped = grouped.sort_values(
            ["feature_count", "membership_mask"],
            ascending=[False, False],
            kind="mergesort",
        ).reset_index(drop=True)
        for rank, row in grouped.iterrows():
            record: dict[str, object] = {
                "panel": panel,
                "panel_title": title,
                "region_class": region,
                "context": context,
                "intersection_rank": rank + 1,
                "membership_mask": row["membership_mask"],
                "feature_count": int(row["feature_count"]),
                "features": row["features"],
            }
            for dataset, column in zip(DATASETS, membership_columns, strict=True):
                record[dataset] = int(row[column])
            rows.append(record)

    result = pd.DataFrame(rows)
    plotted = memberships.loc[memberships["plotted_methylation_feature"] == 1]
    actual_nonempty = set(
        plotted.groupby(["region_class", "context"], observed=False).size().loc[lambda x: x > 0].index
    )
    planned_nonempty = {(region, context) for _, region, context, _ in PANEL_SPECS}
    if actual_nonempty != planned_nonempty:
        missing = sorted(actual_nonempty - planned_nonempty)
        unexpected_empty = sorted(planned_nonempty - actual_nonempty)
        raise ValueError(
            "non-empty region/context cells no longer match the reviewed composite layout; "
            f"unplotted={missing}, empty_planned={unexpected_empty}"
        )
    return result


def add_panel_heading(ax: Axes, panel: str, title: str) -> None:
    """Draw a parenthesized label while keeping the letter auditable."""
    base = ax.transAxes
    point_transform = ax.figure.dpi_scale_trans
    y_offset = 7.0 / 72.0
    pieces = (
        ("(", -4.0, "bold"),
        (panel, 0.0, "bold"),
        (")", 4.8, "bold"),
    )
    for text, x_offset, weight in pieces:
        ax.text(
            0,
            1,
            text,
            transform=base + ScaledTranslation(x_offset / 72.0, y_offset, point_transform),
            ha="left",
            va="bottom",
            fontsize=7.2,
            fontweight=weight,
            color=INK,
            clip_on=False,
        )
    ax.text(
        0.5, 1, title,
        transform=base + ScaledTranslation(0, y_offset, point_transform),
        ha="center", va="bottom", fontsize=7.2,
        fontweight="semibold", color=INK, clip_on=False,
    )


def plot_overview(ax: Axes, overview: pd.DataFrame) -> None:
    matrix = np.zeros((len(REGIONS), len(CONTEXTS)), dtype=int)
    for row in overview.itertuples(index=False):
        matrix[REGIONS.index(row.region_class), CONTEXTS.index(row.context)] = int(
            row.unique_feature_count
        )

    color_map = mpl_colors.LinearSegmentedColormap.from_list(
        "muted_blue", ["#F5F7FA", "#B8CEE5", "#2F5D8A"]
    )
    normalizer = mpl_colors.Normalize(vmin=0, vmax=max(1, int(matrix.max())))
    for row_index in range(matrix.shape[0]):
        for column_index in range(matrix.shape[1]):
            value = int(matrix[row_index, column_index])
            facecolor = color_map(normalizer(value))
            ax.add_patch(
                Rectangle(
                    (column_index - 0.5, row_index - 0.5),
                    1,
                    1,
                    facecolor=facecolor,
                    edgecolor="white",
                    linewidth=1.2,
                )
            )
            ax.text(
                column_index,
                row_index,
                str(value),
                ha="center",
                va="center",
                fontsize=7.2,
                fontweight="semibold",
                color="white" if normalizer(value) > 0.52 else INK,
            )

    ax.set_xlim(-0.5, len(CONTEXTS) - 0.5)
    ax.set_ylim(len(REGIONS) - 0.5, -0.5)
    ax.set_xticks(range(len(CONTEXTS)), CONTEXTS)
    ax.set_yticks(range(len(REGIONS)), [REGION_LABELS[region] for region in REGIONS])
    ax.tick_params(axis="both", length=0, pad=2.5)
    for label, context in zip(ax.get_xticklabels(), CONTEXTS, strict=True):
        label.set_color(CONTEXT_COLORS[context])
        label.set_fontweight("semibold")
    ax.set_xlabel("Methylation context", labelpad=4)
    for spine in ax.spines.values():
        spine.set_visible(False)
    add_panel_heading(ax, "a", "Regional overview")


def plot_upset_panel(ax: Axes, panel_rows: pd.DataFrame, panel: str, title: str) -> None:
    """Draw a compact UpSet chart in one axes for reliable panel alignment."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    add_panel_heading(ax, panel, title)

    context = str(panel_rows["context"].iloc[0])
    signal_color = CONTEXT_COLORS[context]
    member_y = np.linspace(0.43, 0.09, len(DATASETS))
    matrix_left = 0.37
    matrix_right = 0.985
    matrix_width = matrix_right - matrix_left
    bar_base = 0.56
    bar_height = 0.27

    counts = panel_rows["feature_count"].astype(int).to_numpy()
    number_of_intersections = len(panel_rows)
    slot_width = matrix_width / number_of_intersections
    x_positions = matrix_left + slot_width * (np.arange(number_of_intersections) + 0.5)
    maximum_count = int(counts.max())

    ax.plot(
        [matrix_left, matrix_right],
        [bar_base, bar_base],
        transform=ax.transAxes,
        color=GRID_COLOR,
        linewidth=0.65,
        clip_on=False,
    )
    for x_value, count in zip(x_positions, counts, strict=True):
        height = bar_height * count / maximum_count
        width = min(slot_width * 0.68, 0.055)
        ax.add_patch(
            Rectangle(
                (x_value - width / 2, bar_base),
                width,
                height,
                transform=ax.transAxes,
                facecolor=signal_color,
                edgecolor="none",
            )
        )
        ax.text(
            x_value,
            bar_base + height + 0.014,
            str(count),
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=5.3,
            color=INK,
        )

    set_sizes = {
        dataset: int(
            (
                panel_rows[dataset].astype(int)
                * panel_rows["feature_count"].astype(int)
            ).sum()
        )
        for dataset in DATASETS
    }
    maximum_set_size = max(set_sizes.values())
    ax.text(
        0.02,
        0.50,
        "Set size",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=5.4,
        color=MUTED_INK,
    )

    for row_index, (dataset, y_value) in enumerate(zip(DATASETS, member_y, strict=True)):
        if row_index % 2 == 0:
            ax.add_patch(
                Rectangle(
                    (0.0, y_value - 0.034),
                    1.0,
                    0.068,
                    transform=ax.transAxes,
                    facecolor=ROW_SHADE,
                    edgecolor="none",
                    zorder=0,
                )
            )
        size = set_sizes[dataset]
        set_width = 0.14 * size / maximum_set_size if maximum_set_size else 0.0
        if size > 0:
            ax.add_patch(
                Rectangle(
                    (0.17 - set_width, y_value - 0.014),
                    set_width,
                    0.028,
                    transform=ax.transAxes,
                    facecolor=signal_color,
                    edgecolor="none",
                    zorder=2,
                )
            )
        ax.text(
            0.185,
            y_value,
            str(size),
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=5.2,
            color=MUTED_INK,
        )
        ax.text(
            0.32,
            y_value,
            DATASET_LABELS[dataset],
            transform=ax.transAxes,
            ha="right",
            va="center",
            fontsize=5.8,
            color=INK,
        )

    for x_value, row in zip(x_positions, panel_rows.itertuples(index=False), strict=True):
        active_y = [
            y_value
            for dataset, y_value in zip(DATASETS, member_y, strict=True)
            if int(getattr(row, dataset)) == 1
        ]
        if len(active_y) > 1:
            ax.plot(
                [x_value, x_value],
                [min(active_y), max(active_y)],
                transform=ax.transAxes,
                color=INK,
                linewidth=0.75,
                zorder=3,
            )
        for dataset, y_value in zip(DATASETS, member_y, strict=True):
            active = int(getattr(row, dataset)) == 1
            ax.add_patch(
                Circle(
                    (x_value, y_value),
                    radius=0.0105 if active else 0.0082,
                    transform=ax.transAxes,
                    facecolor=INK if active else INACTIVE,
                    edgecolor="none",
                    zorder=4,
                )
            )


def build_figure(overview: pd.DataFrame, intersections: pd.DataFrame) -> tuple[Figure, list[Axes]]:
    configure_style()
    fig = plt.figure(figsize=(FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN), facecolor="white")
    grid = fig.add_gridspec(
        2,
        2,
        left=0.075,
        right=0.985,
        bottom=0.105,
        top=0.945,
        wspace=0.22,
        hspace=0.27,
        height_ratios=[1.0, 1.0],
    )
    axes = [
        fig.add_subplot(grid[0, 0]),
        fig.add_subplot(grid[0, 1]),
        fig.add_subplot(grid[1, 0]),
        fig.add_subplot(grid[1, 1]),
    ]

    for axis, (panel, _region, _context, title) in zip(axes, PANEL_SPECS, strict=True):
        plot_upset_panel(
            axis,
            intersections.loc[intersections["panel"] == panel].sort_values("intersection_rank"),
            panel,
            title,
        )

    fig.text(
        0.5,
        0.035,
        "Top 10%: highest positive mean(|SHAP|) features within each response, ranked across "
        "CHH, CHG, CG and SNP; SNPs are excluded after selection.",
        ha="center",
        va="center",
        fontsize=5.5,
        color=MUTED_INK,
    )
    fig.text(
        0.5,
        0.017,
        "Counts represent unique genomic features; matrix rows denote membership in α, α%, β, β% and total response sets.",
        ha="center",
        va="center",
        fontsize=5.5,
        color=MUTED_INK,
    )
    return fig, axes


def qa_scripts_directory() -> Path:
    return Path(
        os.environ.get(
            "NATURE_FIGURE_QA_DIR",
            Path.home() / ".codex" / "skills" / "nature-figure" / "scripts",
        )
    )


def load_alignment_helper():
    directory = qa_scripts_directory()
    helper_path = directory / "audit_panel_alignment.py"
    if not helper_path.is_file():
        _skip_optional_qa("Optional panel-alignment helper not found")
        return _skip_alignment
    sys.path.insert(0, str(directory))
    try:
        from audit_panel_alignment import require_matplotlib_panel_alignment
    except ImportError as exc:
        _skip_optional_qa(str(exc))
        return _skip_alignment
    return require_matplotlib_panel_alignment


def run_pdf_text_audit(pdf_path: Path, report_path: Path, final_pdf_path: Path) -> None:
    script = qa_scripts_directory() / "audit_pdf_text.py"
    if not script.is_file():
        _skip_optional_qa(f"Optional QA tool not found: {script}", report_path)
        return
    completed = subprocess.run(
        [sys.executable, str(script), str(pdf_path), "--min-pt", "5", "--json"],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        details = (completed.stdout + "\n" + completed.stderr).strip()
        raise RuntimeError(f"PDF text audit failed: {details}")
    report = json.loads(completed.stdout)
    report["pdf"] = str(final_pdf_path)
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def run_collision_audit(pdf_path: Path, report_path: Path, final_pdf_path: Path) -> None:
    script = qa_scripts_directory() / "audit_figure_collisions.py"
    if not script.is_file():
        _skip_optional_qa(f"Optional QA tool not found: {script}", report_path)
        return
    completed = subprocess.run(
        [
            sys.executable,
            str(script),
            str(pdf_path),
            "--json-out",
            str(report_path),
            "--strict",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        details = (completed.stdout + "\n" + completed.stderr).strip()
        raise RuntimeError(f"strict rendered-PDF collision audit failed: {details}")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if report.get("verdict") != "PASS":
        raise RuntimeError(f"collision audit did not pass: {report.get('verdict', 'unknown')}")
    report["pdf"] = str(final_pdf_path)
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def write_table(table: pd.DataFrame, path: Path) -> None:
    table.to_csv(path, sep="\t", index=False, encoding="utf-8", float_format="%.12g")


def export_outputs(
    fig: Figure,
    axes: Sequence[Axes],
    selection_summary: pd.DataFrame,
    memberships: pd.DataFrame,
    overview: pd.DataFrame,
    intersections: pd.DataFrame,
    *,
    output_dir: Path,
    basename: str,
    dpi: int,
) -> list[Path]:
    """Stage, validate and atomically publish the complete figure bundle."""
    output_dir.mkdir(parents=True, exist_ok=True)
    suffixes = (
        ".pdf",
        ".svg",
        ".png",
        ".tiff",
        ".selection_summary.tsv",
        ".memberships.tsv",
        ".overview.tsv",
        ".intersections.tsv",
        ".alignment.json",
        ".alignment.svg",
        ".pdf-text-audit.json",
        ".collision-audit.json",
    )
    final_paths = [output_dir / f"{basename}{suffix}" for suffix in suffixes]
    final_pdf = final_paths[0]

    with tempfile.TemporaryDirectory(prefix=f".{basename}.", dir=output_dir) as temp_name:
        temp_dir = Path(temp_name)
        staged = {path.name: temp_dir / path.name for path in final_paths}
        temporary_paths = [staged[path.name] for path in final_paths]
        (
            pdf_path,
            svg_path,
            png_path,
            tiff_path,
            selection_path,
            memberships_path,
            overview_path,
            intersections_path,
            alignment_path,
            alignment_overlay_path,
            text_audit_path,
            collision_path,
        ) = temporary_paths

        write_table(selection_summary, selection_path)
        write_table(memberships, memberships_path)
        write_table(overview, overview_path)
        write_table(intersections, intersections_path)

        require_matplotlib_panel_alignment = load_alignment_helper()
        alignment_report = require_matplotlib_panel_alignment(
            fig,
            axes=list(axes),
            panel_ids=["a", "b", "c", "d"],
            row_groups=[["a", "b"], ["c", "d"]],
            json_out=alignment_path,
            overlay_svg=alignment_overlay_path,
            tolerance_pt=1.5,
            gutter_tolerance_pt=1.5,
            require_panel_labels=True,
            strict=True,
        )
        if alignment_report.get("verdict") == "SKIPPED":
            # An unavailable auditor does not produce a diagnostic SVG overlay.
            temporary_paths.remove(alignment_overlay_path)
            final_paths.remove(output_dir / alignment_overlay_path.name)

        fig.savefig(pdf_path, format="pdf")
        fig.savefig(svg_path, format="svg")
        fig.savefig(png_path, format="png", dpi=dpi)
        fig.savefig(
            tiff_path,
            format="tiff",
            dpi=dpi,
            pil_kwargs={"compression": "tiff_lzw"},
        )
        run_pdf_text_audit(pdf_path, text_audit_path, final_pdf)
        run_collision_audit(pdf_path, collision_path, final_pdf)

        for source, destination in zip(temporary_paths, final_paths, strict=True):
            os.replace(source, destination)
    return final_paths


def validate_expected_current_data(
    selection_summary: pd.DataFrame,
    memberships: pd.DataFrame,
    overview: pd.DataFrame,
) -> None:
    """Guard the reviewed figure against silent changes in the active RCA inputs."""
    observed_selected = selection_summary.set_index("dataset")["selected_feature_count"].to_dict()
    if observed_selected != {dataset: 29 for dataset in DATASETS}:
        raise ValueError(f"expected 29 selected features per dataset; observed {observed_selected}")
    if len(memberships) != 52:
        raise ValueError(f"expected a 52-feature top-decile union; observed {len(memberships)}")
    plotted_count = int(memberships["plotted_methylation_feature"].sum())
    if plotted_count != 48:
        raise ValueError(f"expected 48 plotted methylation features; observed {plotted_count}")
    observed_overview = {
        (row.region_class, row.context): int(row.unique_feature_count)
        for row in overview.itertuples(index=False)
    }
    expected_overview = {
        ("upstream", "CHH"): 36,
        ("upstream", "CHG"): 8,
        ("upstream", "CG"): 2,
        ("gene_body", "CHH"): 0,
        ("gene_body", "CHG"): 0,
        ("gene_body", "CG"): 0,
        ("downstream", "CHH"): 0,
        ("downstream", "CHG"): 0,
        ("downstream", "CG"): 2,
    }
    if observed_overview != expected_overview:
        raise ValueError(
            "regional/context feature counts changed from the reviewed layout; "
            f"observed {observed_overview}"
        )


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        validate_args(args)
        tables = {dataset: load_dataset_table(args.input_dir, dataset) for dataset in DATASETS}
        validate_cross_dataset_metadata(tables)
        selections, selection_summary = select_top_features(tables)
        memberships = build_membership_table(selections)
        overview = build_overview_table(memberships)
        intersections = build_intersection_table(memberships)
        validate_expected_current_data(selection_summary, memberships, overview)
        fig, axes = build_figure(overview, intersections)
        try:
            outputs = export_outputs(
                fig,
                axes,
                selection_summary,
                memberships,
                overview,
                intersections,
                output_dir=args.output_dir,
                basename=args.basename,
                dpi=args.dpi,
            )
        finally:
            plt.close(fig)
    except (FileNotFoundError, ImportError, OSError, RuntimeError, TypeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(
        f"Selected {int(selection_summary['selected_feature_count'].sum())} dataset-feature "
        f"memberships; union={len(memberships)}, plotted methylation union="
        f"{int(memberships['plotted_methylation_feature'].sum())}."
    )
    for output in outputs:
        print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
