#!/usr/bin/env python3
"""Plot RCA abundance and alpha/beta composition across temperatures."""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
from pathlib import Path
from typing import Sequence

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rca-matplotlib-cache")
)

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.transforms import ScaledTranslation
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[2]
DEFAULT_INPUT = PROJECT_ROOT / "list" / "RCA" / "RCA.rnadata.tsv"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "list" / "RCA" / "experiment"
DEFAULT_BASENAME = "RCA_temperature_response"

EXPECTED_TEMPERATURES = np.array([4, 8, 12, 16, 20, 24, 28, 32], dtype=float)
EXPECTED_REPLICATES = 3
REQUIRED_COLUMNS = ("Temperature", "total", "α", "β", "α/%", "β/%")
EXPRESSION_COLUMNS = ("total", "α", "β")
PROPORTION_COLUMNS = ("α/%", "β/%")
PROPORTION_SUM_TOLERANCE = 0.02

COLORS = {
    "total": "#4D4D4D",
    "alpha": "#2F6DAE",
    "beta": "#D9822B",
    "control": "#E8E8E8",
    "reference": "#B8B8B8",
    "grid": "#E3E3E3",
}

FIGURE_WIDTH_IN = 183.0 / 25.4
FIGURE_HEIGHT_IN = 95.0 / 25.4
POINT_OFFSETS = np.linspace(-0.65, 0.65, EXPECTED_REPLICATES)
GENE_LABEL = r"$\mathbf{RCA}$"
GENE_LABEL_ITALIC = r"$\mathit{RCA}$"


def configure_style() -> None:
    """Apply a compact, editable-text publication style."""
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "mathtext.fontset": "custom",
            "mathtext.rm": "DejaVu Sans",
            "mathtext.it": "DejaVu Sans:italic",
            "mathtext.bf": "DejaVu Sans:italic:bold",
            "mathtext.cal": "DejaVu Sans",
            "mathtext.sf": "DejaVu Sans",
            "mathtext.tt": "DejaVu Sans",
            "font.size": 6.5,
            "axes.titlesize": 7.3,
            "axes.titleweight": "semibold",
            "axes.labelsize": 7.0,
            "axes.linewidth": 0.7,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.labelsize": 6.2,
            "ytick.labelsize": 6.2,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "legend.frameon": False,
            "savefig.facecolor": "white",
        }
    )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot total RCA, RCA α/β abundance, and α/β proportions across "
            "eight treatment temperatures."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"Input TSV (default: {DEFAULT_INPUT})",
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
        help=f"Output basename without extension (default: {DEFAULT_BASENAME})",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=600,
        help="PNG resolution in dots per inch (default: 600)",
    )
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> None:
    if args.dpi < 300:
        raise ValueError("--dpi must be at least 300 for the high-resolution PNG")
    if not args.basename or Path(args.basename).name != args.basename:
        raise ValueError("--basename must be a non-empty filename stem, not a path")


def load_and_validate_data(path: Path) -> pd.DataFrame:
    """Read the source table and enforce the declared experimental design."""
    if not path.is_file():
        raise FileNotFoundError(f"input TSV not found: {path}")

    data = pd.read_csv(path, sep="\t", encoding="utf-8")
    missing = [column for column in REQUIRED_COLUMNS if column not in data.columns]
    if missing:
        raise ValueError(f"input TSV is missing required columns: {', '.join(missing)}")

    data = data.loc[:, REQUIRED_COLUMNS].copy()
    for column in REQUIRED_COLUMNS:
        data[column] = pd.to_numeric(data[column], errors="raise")

    numeric_values = data.loc[:, REQUIRED_COLUMNS].to_numpy(dtype=float)
    if not np.isfinite(numeric_values).all():
        raise ValueError("input TSV contains missing, infinite, or non-numeric values")

    observed_temperatures = np.sort(data["Temperature"].unique().astype(float))
    if not np.array_equal(observed_temperatures, EXPECTED_TEMPERATURES):
        raise ValueError(
            "temperatures must be exactly: "
            + ", ".join(f"{temperature:g}" for temperature in EXPECTED_TEMPERATURES)
        )

    replicate_counts = (
        data.groupby("Temperature", sort=True, observed=False).size().reindex(EXPECTED_TEMPERATURES)
    )
    if replicate_counts.isna().any() or not (replicate_counts == EXPECTED_REPLICATES).all():
        observed = ", ".join(
            f"{temperature:g} °C={int(count)}"
            for temperature, count in replicate_counts.fillna(0).items()
        )
        raise ValueError(
            f"each temperature must contain {EXPECTED_REPLICATES} rows; observed {observed}"
        )

    if (data.loc[:, EXPRESSION_COLUMNS] < 0).any().any():
        raise ValueError("TPM-derived expression columns must be non-negative")
    if ((data.loc[:, PROPORTION_COLUMNS] < 0) | (data.loc[:, PROPORTION_COLUMNS] > 100)).any().any():
        raise ValueError("proportion columns must lie within 0–100%")

    if not np.allclose(
        data["total"].to_numpy(dtype=float),
        (data["α"] + data["β"]).to_numpy(dtype=float),
        rtol=1e-6,
        atol=0.02,
    ):
        max_difference = float(np.max(np.abs(data["total"] - data["α"] - data["β"])))
        raise ValueError(
            "total is inconsistent with α + β beyond tolerance "
            f"(maximum absolute difference: {max_difference:.6g})"
        )

    proportion_sum = data["α/%"] + data["β/%"]
    proportion_difference = np.abs(proportion_sum.to_numpy(dtype=float) - 100.0)
    # Include the declared ±0.02 boundary. One ULP at 100 guards against
    # binary floating-point representation turning 0.02 into 0.02000000000001.
    comparison_tolerance = PROPORTION_SUM_TOLERANCE + np.spacing(100.0)
    if np.any(proportion_difference > comparison_tolerance):
        max_difference = float(np.max(proportion_difference))
        raise ValueError(
            "α/% + β/% is inconsistent with 100% beyond tolerance "
            f"(maximum absolute difference: {max_difference:.6g})"
        )

    data["_input_order"] = np.arange(len(data))
    data = data.sort_values(["Temperature", "_input_order"], kind="stable").reset_index(drop=True)
    for column in EXPRESSION_COLUMNS:
        data[f"{column}_ln"] = np.log1p(data[column].to_numpy(dtype=float))

    return data


def summarize(data: pd.DataFrame, value_column: str) -> pd.DataFrame:
    summary = (
        data.groupby("Temperature", sort=True, observed=False)[value_column]
        .agg(mean="mean", sd="std", n="size")
        .reindex(EXPECTED_TEMPERATURES)
    )
    if summary.isna().any().any() or not (summary["n"] == EXPECTED_REPLICATES).all():
        raise RuntimeError(f"failed to summarize all temperatures for {value_column}")
    return summary


def add_panel_label(ax: Axes, label: str) -> None:
    offset = ScaledTranslation(-8.0 / 72.0, 3.0 / 72.0, ax.figure.dpi_scale_trans)
    artist = ax.text(
        0,
        1,
        f"({label})",
        transform=ax.transAxes + offset,
        fontsize=8.0,
        fontweight="bold",
        ha="right",
        va="bottom",
        clip_on=False,
    )
    artist.set_gid("parenthesized-panel-label")


def validate_parenthesized_panel_labels(
    fig: Figure,
    axes: Sequence[Axes],
    *,
    panel_ids: Sequence[str],
    tolerance_pt: float = 1.5,
) -> dict[str, object]:
    """Validate '(a)' labels that the generic single-letter auditor cannot parse."""
    if len(axes) != len(panel_ids):
        raise RuntimeError("panel-label validation received mismatched axes and panel ids")

    fig.canvas.draw()
    figure_width_pt, _ = fig.get_size_inches() * 72.0
    anchors: dict[str, list[float]] = {}
    x_offsets: dict[str, float] = {}
    for ax, panel_id in zip(axes, panel_ids):
        candidates = [
            artist
            for artist in ax.texts
            if artist.get_gid() == "parenthesized-panel-label"
        ]
        if len(candidates) != 1:
            raise RuntimeError(
                f"panel {panel_id} must contain exactly one parenthesized panel label"
            )
        artist = candidates[0]
        expected = f"({panel_id})"
        if artist.get_text() != expected or str(artist.get_fontweight()).lower() != "bold":
            raise RuntimeError(f"panel {panel_id} label must be bold text {expected}")

        display_xy = artist.get_transform().transform(artist.get_position())
        inches_xy = fig.dpi_scale_trans.inverted().transform(display_xy)
        anchor_pt = [float(inches_xy[0] * 72.0), float(inches_xy[1] * 72.0)]
        axes_left_pt = float(ax.get_position().x0 * figure_width_pt)
        anchors[panel_id] = anchor_pt
        x_offsets[panel_id] = anchor_pt[0] - axes_left_pt

    row_groups = [list(panel_ids)]
    row_y_spreads: dict[str, float] = {}
    for index, row in enumerate(row_groups, 1):
        y_values = [anchors[panel_id][1] for panel_id in row]
        spread = float(max(y_values) - min(y_values))
        row_y_spreads[f"row-{index}"] = spread
        if spread > tolerance_pt:
            raise RuntimeError(
                f"parenthesized panel-label y anchors differ by {spread:.3f} pt "
                f"in row {index}; tolerance is {tolerance_pt:.3f} pt"
            )

    x_offset_spread = float(max(x_offsets.values()) - min(x_offsets.values()))
    if x_offset_spread > tolerance_pt:
        raise RuntimeError(
            "parenthesized panel-label x offsets are inconsistent by "
            f"{x_offset_spread:.3f} pt; tolerance is {tolerance_pt:.3f} pt"
        )

    return {
        "verdict": "PASS",
        "format": f"({panel_ids[0]})-({panel_ids[-1]})",
        "fontweight": "bold",
        "tolerance_pt": tolerance_pt,
        "anchors_pt": anchors,
        "row_y_spreads_pt": row_y_spreads,
        "x_offset_spread_pt": x_offset_spread,
    }


def style_axis(ax: Axes, *, proportion: bool) -> None:
    ax.axvspan(22, 26, color=COLORS["control"], alpha=0.75, linewidth=0, zorder=0)
    if proportion:
        ax.axhline(
            50,
            color=COLORS["reference"],
            linestyle=(0, (3, 2)),
            linewidth=0.7,
            zorder=1,
        )
    ax.set_xlim(2, 34)
    ax.set_xticks(EXPECTED_TEMPERATURES)
    ax.set_xticklabels([f"{temperature:g}" for temperature in EXPECTED_TEMPERATURES])
    ax.grid(axis="y", color=COLORS["grid"], linewidth=0.45, zorder=0)
    ax.set_axisbelow(True)


def plot_series(
    ax: Axes,
    data: pd.DataFrame,
    *,
    value_column: str,
    color: str,
    marker: str,
) -> pd.DataFrame:
    summary = summarize(data, value_column)

    for temperature in EXPECTED_TEMPERATURES:
        values = data.loc[data["Temperature"] == temperature, value_column].to_numpy(dtype=float)
        if len(values) != EXPECTED_REPLICATES:
            raise RuntimeError(
                f"expected {EXPECTED_REPLICATES} samples at {temperature:g} °C in {value_column}"
            )
        ax.scatter(
            temperature + POINT_OFFSETS,
            values,
            s=15,
            facecolor=color,
            edgecolor="white",
            marker=marker,
            linewidth=0.45,
            alpha=0.68,
            zorder=3,
        )

    means = summary["mean"].to_numpy(dtype=float)
    standard_deviations = summary["sd"].to_numpy(dtype=float)
    ax.plot(
        EXPECTED_TEMPERATURES,
        means,
        color=color,
        linewidth=1.05,
        alpha=0.92,
        zorder=4,
    )
    ax.errorbar(
        EXPECTED_TEMPERATURES,
        means,
        yerr=standard_deviations,
        fmt=marker,
        markersize=3.4,
        markerfacecolor=color,
        markeredgecolor="white",
        markeredgewidth=0.45,
        color=color,
        ecolor=color,
        elinewidth=0.85,
        capsize=2.3,
        capthick=0.85,
        zorder=5,
    )
    return summary


def add_series_legend(
    ax: Axes,
    specs: Sequence[tuple[str, str, str, str]],
    *,
    location: str = "lower center",
    anchor: tuple[float, float] | None = None,
) -> None:
    """Add a compact legend for series overlaid in one panel."""
    handles = [
        Line2D(
            [0],
            [0],
            color=color,
            marker=marker,
            markersize=3.8,
            markeredgecolor="white",
            markeredgewidth=0.45,
            linewidth=1.05,
            label=legend_label,
        )
        for _, legend_label, color, marker in specs
    ]
    ax.legend(
        handles=handles,
        loc=location,
        bbox_to_anchor=anchor,
        ncol=len(handles),
        fontsize=6.2,
        handlelength=1.8,
        columnspacing=1.2,
        borderaxespad=0.45,
    )


def format_p_value(p_value: float) -> str:
    """Format a two-sided P value compactly without undersized math glyphs."""
    if p_value < 0.001:
        mantissa, exponent = f"{p_value:.2e}".split("e")
        return f"{mantissa}e{int(exponent):+d}".replace("-", "−")
    return f"{p_value:.3f}"


def add_spearman_annotation(
    ax: Axes,
    data: pd.DataFrame,
    specs: Sequence[tuple[str, str, str, str]],
    *,
    x: float,
    y_top: float,
    horizontal_alignment: str,
    line_step: float = 0.065,
) -> dict[str, tuple[float, float]]:
    """Annotate two-sided Spearman correlations using temperature-group means."""
    statistics: dict[str, tuple[float, float]] = {}

    for index, (value_column, label, color, _) in enumerate(specs):
        temperature_summary = summarize(data, value_column)
        if len(temperature_summary) != len(EXPECTED_TEMPERATURES):
            raise RuntimeError(
                f"Spearman correlation for {value_column} must use exactly "
                f"{len(EXPECTED_TEMPERATURES)} temperature means"
            )
        temperature_means = temperature_summary["mean"].to_numpy(dtype=float)
        result = spearmanr(
            EXPECTED_TEMPERATURES,
            temperature_means,
            alternative="two-sided",
        )
        rho = float(result.statistic)
        p_value = float(result.pvalue)
        if not np.isfinite(rho) or not np.isfinite(p_value):
            raise RuntimeError(
                f"Spearman correlation is undefined for Temperature and {value_column}"
            )
        statistics[label] = (rho, p_value)
        line_y = y_top - line_step * index
        rho_text = f"{rho:.3f}".replace("-", "−")
        ax.text(
            x,
            line_y,
            f"{label}: ρ = {rho_text}, P = {format_p_value(p_value)}",
            transform=ax.transAxes,
            fontsize=6.2,
            color=color,
            ha=horizontal_alignment,
            va="top",
            clip_on=False,
        )
    return statistics


def expression_limits(
    data: pd.DataFrame, summaries: Sequence[pd.DataFrame]
) -> tuple[float, float]:
    raw_values = np.concatenate(
        [data[f"{column}_ln"].to_numpy(dtype=float) for column in EXPRESSION_COLUMNS]
    )
    lower_extents = np.concatenate(
        [(summary["mean"] - summary["sd"]).to_numpy(dtype=float) for summary in summaries]
    )
    upper_extents = np.concatenate(
        [(summary["mean"] + summary["sd"]).to_numpy(dtype=float) for summary in summaries]
    )
    lower = float(min(np.min(raw_values), np.min(lower_extents)))
    upper = float(max(np.max(raw_values), np.max(upper_extents)))
    padding = max(0.15, 0.07 * (upper - lower))
    return lower - padding, upper + padding


def proportion_limits(
    data: pd.DataFrame, summaries: Sequence[pd.DataFrame]
) -> tuple[tuple[float, float], np.ndarray]:
    """Return shared percentage limits that include all points and mean ± SD."""
    raw_values = data.loc[:, PROPORTION_COLUMNS].to_numpy(dtype=float).ravel()
    lower_extents = np.concatenate(
        [(summary["mean"] - summary["sd"]).to_numpy(dtype=float) for summary in summaries]
    )
    upper_extents = np.concatenate(
        [(summary["mean"] + summary["sd"]).to_numpy(dtype=float) for summary in summaries]
    )
    lower = float(min(np.min(raw_values), np.min(lower_extents)))
    upper = float(max(np.max(raw_values), np.max(upper_extents)))
    padding = max(0.75, 0.08 * (upper - lower))
    lower_limit = max(0.0, float(np.floor(lower - padding)))
    upper_limit = min(100.0, float(np.ceil(upper + padding)))
    if upper_limit <= lower_limit:
        raise RuntimeError("could not derive a valid proportion-axis range")

    span = upper_limit - lower_limit
    tick_step = 2.0 if span <= 16 else 5.0 if span <= 35 else 10.0
    tick_start = np.ceil(lower_limit / tick_step) * tick_step
    ticks = np.arange(tick_start, upper_limit + 0.5 * tick_step, tick_step)
    return (lower_limit, upper_limit), ticks


def build_figure(data: pd.DataFrame) -> tuple[Figure, list[Axes]]:
    configure_style()
    fig = plt.figure(figsize=(FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN))
    grid = fig.add_gridspec(
        1,
        2,
        left=0.078,
        right=0.985,
        bottom=0.16,
        top=0.70,
        wspace=0.24,
    )

    axes = [
        fig.add_subplot(grid[0, 0]),
        fig.add_subplot(grid[0, 1]),
    ]

    expression_specs = [
        ("total_ln", "Total", COLORS["total"], "D"),
        ("α_ln", "α", COLORS["alpha"], "o"),
        ("β_ln", "β", COLORS["beta"], "s"),
    ]
    expression_ax = axes[0]
    style_axis(expression_ax, proportion=False)
    expression_summaries = [
        plot_series(
            expression_ax,
            data,
            value_column=value_column,
            color=color,
            marker=marker,
        )
        for value_column, _, color, marker in expression_specs
    ]
    expression_ax.set_title(f"{GENE_LABEL} expression", loc="center", pad=5.0)
    expression_ax.set_ylabel("ln(TPM + 1)")
    add_panel_label(expression_ax, "a")
    add_series_legend(
        expression_ax,
        expression_specs,
        location="lower right",
    )

    shared_expression_limits = expression_limits(data, expression_summaries)
    expression_span = shared_expression_limits[1] - shared_expression_limits[0]
    expression_ax.set_ylim(
        shared_expression_limits[0],
        shared_expression_limits[1] + 0.25 * expression_span,
    )
    expression_tick_start = np.ceil(shared_expression_limits[0] * 2.0) / 2.0
    expression_tick_stop = np.floor(shared_expression_limits[1] * 2.0) / 2.0
    expression_ax.set_yticks(
        np.arange(expression_tick_start, expression_tick_stop + 0.25, 0.5)
    )
    add_spearman_annotation(
        expression_ax,
        data,
        expression_specs,
        x=0.02,
        y_top=0.97,
        horizontal_alignment="left",
    )

    proportion_specs = [
        ("α/%", "α", COLORS["alpha"], "o"),
        ("β/%", "β", COLORS["beta"], "s"),
    ]
    proportion_ax = axes[1]
    style_axis(proportion_ax, proportion=True)
    proportion_summaries = [
        plot_series(
            proportion_ax,
            data,
            value_column="α/%",
            color=COLORS["alpha"],
            marker="o",
        ),
        plot_series(
            proportion_ax,
            data,
            value_column="β/%",
            color=COLORS["beta"],
            marker="s",
        ),
    ]
    proportion_ax.set_title(f"{GENE_LABEL} α/β proportion", loc="center", pad=5.0)
    proportion_ax.set_ylabel(f"Proportion of {GENE_LABEL_ITALIC} (%)")
    add_panel_label(proportion_ax, "b")
    add_series_legend(
        proportion_ax,
        proportion_specs,
        location="lower left",
    )
    shared_proportion_limits, proportion_ticks = proportion_limits(
        data, proportion_summaries
    )
    proportion_ax.set_ylim(*shared_proportion_limits)
    proportion_ax.set_yticks(proportion_ticks)
    add_spearman_annotation(
        proportion_ax,
        data,
        proportion_specs,
        x=0.02,
        y_top=0.97,
        horizontal_alignment="left",
    )

    fig.suptitle(
        f"{GENE_LABEL} expression and α/β composition across treatment temperatures\n"
        "3 d treatment; n = 3 independent samples per temperature; "
        "light band = 24 °C control\n"
        "translucent points = independent samples; lines and error bars = mean ± SD\n"
        f"Spearman correlations: two-sided; "
        f"n = {len(EXPECTED_TEMPERATURES)} temperature means",
        fontsize=7.2,
        fontweight="normal",
        linespacing=1.35,
        y=0.975,
    )
    fig.supxlabel("Temperature (°C)", fontsize=7.0, y=0.045)
    return fig, axes


def load_alignment_helper():
    qa_directory = Path(
        os.environ.get(
            "NATURE_FIGURE_QA_DIR",
            Path.home() / ".codex" / "skills" / "nature-figure" / "scripts",
        )
    )
    helper_path = qa_directory / "audit_panel_alignment.py"
    if not helper_path.is_file():
        raise RuntimeError(
            "panel-alignment helper not found; set NATURE_FIGURE_QA_DIR to the "
            "directory containing audit_panel_alignment.py"
        )
    sys.path.insert(0, str(qa_directory))
    try:
        from audit_panel_alignment import require_matplotlib_panel_alignment
    except ImportError as exc:
        raise RuntimeError(f"could not import panel-alignment helper: {exc}") from exc
    return require_matplotlib_panel_alignment


def export_figure(
    fig: Figure,
    axes: Sequence[Axes],
    *,
    output_dir: Path,
    basename: str,
    dpi: int,
) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    output_base = output_dir / basename
    panel_ids = ["a", "b"]
    validate_parenthesized_panel_labels(
        fig,
        axes,
        panel_ids=panel_ids,
        tolerance_pt=1.5,
    )
    require_matplotlib_panel_alignment = load_alignment_helper()
    require_matplotlib_panel_alignment(
        fig,
        axes=list(axes),
        panel_ids=panel_ids,
        row_groups=[panel_ids],
        tolerance_pt=1.5,
        gutter_tolerance_pt=1.5,
        require_panel_labels=False,
        strict=True,
    )

    outputs = [
        output_base.with_suffix(".pdf"),
        output_base.with_suffix(".png"),
    ]
    fig.savefig(outputs[0], format="pdf")
    fig.savefig(outputs[1], format="png", dpi=dpi)
    return outputs


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        validate_args(args)
        data = load_and_validate_data(args.input)
        fig, axes = build_figure(data)
        try:
            outputs = export_figure(
                fig,
                axes,
                output_dir=args.output_dir,
                basename=args.basename,
                dpi=args.dpi,
            )
        finally:
            plt.close(fig)
    except (FileNotFoundError, OSError, RuntimeError, TypeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(
        f"Validated {len(data)} samples across {len(EXPECTED_TEMPERATURES)} "
        f"temperatures ({EXPECTED_REPLICATES} replicates per temperature)."
    )
    for output in outputs:
        print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
