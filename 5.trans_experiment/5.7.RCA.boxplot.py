#!/usr/bin/env python3
"""Plot RCA temperature-response bars with ANOVA, Dunnett, and Spearman tests."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Sequence

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rca-boxplot-matplotlib-cache")
)
os.environ.setdefault(
    "XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "rca-boxplot-xdg-cache")
)

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import PercentFormatter
from matplotlib.transforms import ScaledTranslation
import numpy as np
import pandas as pd
from scipy import stats


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[2]
DEFAULT_INPUT = PROJECT_ROOT / "list" / "RCA" / "RCA.rnadata.tsv"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "list" / "RCA" / "experiment"
DEFAULT_BASENAME = "RCA_temperature_boxplot"

EXPECTED_TEMPERATURES = np.array([4, 8, 12, 16, 20, 24, 28, 32], dtype=float)
CONTROL_TEMPERATURE = 24.0
CONTROL_INDEX = int(np.flatnonzero(EXPECTED_TEMPERATURES == CONTROL_TEMPERATURE)[0])
TREATMENT_TEMPERATURES = tuple(
    float(value) for value in EXPECTED_TEMPERATURES if value != CONTROL_TEMPERATURE
)
EXPECTED_REPLICATES = 3
REQUIRED_COLUMNS = ("Temperature", "total", "α", "β", "α/%", "β/%")
EXPRESSION_COLUMNS = ("total", "α", "β")
PROPORTION_COLUMNS = ("α/%", "β/%")
PROPORTION_SUM_TOLERANCE = 0.02

COLORS = {
    "total": "#4D4D4D",
    "alpha": "#2F6DAE",
    "beta": "#D9822B",
    "control": "#E9E9E9",
    "grid": "#E1E1E1",
    "reference": "#AFAFAF",
}

FIGURE_WIDTH_IN = 183.0 / 25.4
FIGURE_HEIGHT_IN = 125.0 / 25.4
BOX_POSITIONS = np.arange(len(EXPECTED_TEMPERATURES), dtype=float)
POINT_OFFSETS = np.linspace(-0.13, 0.13, EXPECTED_REPLICATES)
PANEL_IDS = ("a", "b", "c", "d", "e")


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


@dataclass(frozen=True)
class MetricSpec:
    """Configuration for one plotted and tested RCA response."""

    key: str
    source_column: str
    analysis_column: str
    display_label: str
    panel_title: str
    color: str
    analysis_scale: str
    is_proportion: bool


@dataclass(frozen=True)
class MetricAnalysis:
    """ANOVA, Dunnett, and Spearman results for one response metric."""

    spec: MetricSpec
    groups: dict[float, np.ndarray]
    anova_f: float
    anova_p: float
    spearman_rho: float
    spearman_p: float
    spearman_n: int
    dunnett_statistics: dict[float, float]
    dunnett_pvalues: dict[float, float]
    dunnett_seed: int


METRIC_SPECS = (
    MetricSpec(
        key="total",
        source_column="total",
        analysis_column="total_ln",
        display_label="RCA total expression",
        panel_title=r"$\mathbf{RCA}$ total expression",
        color=COLORS["total"],
        analysis_scale="ln(TPM + 1)",
        is_proportion=False,
    ),
    MetricSpec(
        key="alpha",
        source_column="α",
        analysis_column="α_ln",
        display_label="RCA α expression",
        panel_title=r"$\mathbf{RCA}$ α expression",
        color=COLORS["alpha"],
        analysis_scale="ln(TPM + 1)",
        is_proportion=False,
    ),
    MetricSpec(
        key="beta",
        source_column="β",
        analysis_column="β_ln",
        display_label="RCA β expression",
        panel_title=r"$\mathbf{RCA}$ β expression",
        color=COLORS["beta"],
        analysis_scale="ln(TPM + 1)",
        is_proportion=False,
    ),
    MetricSpec(
        key="alpha_percent",
        source_column="α/%",
        analysis_column="α/%",
        display_label="RCA α proportion",
        panel_title=r"$\mathbf{RCA}$ α proportion",
        color=COLORS["alpha"],
        analysis_scale="percentage",
        is_proportion=True,
    ),
    MetricSpec(
        key="beta_percent",
        source_column="β/%",
        analysis_column="β/%",
        display_label="RCA β proportion",
        panel_title=r"$\mathbf{RCA}$ β proportion",
        color=COLORS["beta"],
        analysis_scale="percentage",
        is_proportion=True,
    ),
)


def configure_style() -> None:
    """Apply a compact publication style with editable PDF text."""
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
            "font.size": 6.4,
            "axes.titlesize": 6.8,
            "axes.titleweight": "semibold",
            "axes.labelsize": 6.8,
            "axes.linewidth": 0.7,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.labelsize": 5.8,
            "ytick.labelsize": 5.9,
            "xtick.major.size": 2.4,
            "ytick.major.size": 2.4,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "savefig.facecolor": "white",
        }
    )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot five RCA temperature-response bar charts and perform one-way "
            "ANOVA followed by two-sided Dunnett comparisons against 24 °C."
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
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help=(
            "Family-wise significance level recorded in the statistics table "
            "(default: 0.05; figure stars use the conventional fixed cutoffs)"
        ),
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260903,
        help="Base seed for reproducible Dunnett integration (default: 20260903)",
    )
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> None:
    if args.dpi < 300:
        raise ValueError("--dpi must be at least 300 for the high-resolution PNG")
    if not (0.0 < args.alpha < 1.0):
        raise ValueError("--alpha must be strictly between 0 and 1")
    if args.seed < 0:
        raise ValueError("--seed must be a non-negative integer")
    if not args.basename or Path(args.basename).name != args.basename:
        raise ValueError("--basename must be a non-empty filename stem, not a path")
    if Path(args.basename).suffix:
        raise ValueError("--basename must not contain a filename extension")


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
        try:
            data[column] = pd.to_numeric(data[column], errors="raise")
        except (TypeError, ValueError) as exc:
            raise ValueError(f"column {column} contains non-numeric values") from exc

    values = data.loc[:, REQUIRED_COLUMNS].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("input TSV contains missing or infinite numeric values")

    observed_temperatures = np.sort(data["Temperature"].unique().astype(float))
    if not np.array_equal(observed_temperatures, EXPECTED_TEMPERATURES):
        expected = ", ".join(f"{value:g}" for value in EXPECTED_TEMPERATURES)
        raise ValueError(f"temperatures must be exactly: {expected}")

    counts = (
        data.groupby("Temperature", sort=True, observed=False)
        .size()
        .reindex(EXPECTED_TEMPERATURES)
    )
    if counts.isna().any() or not (counts == EXPECTED_REPLICATES).all():
        observed = ", ".join(
            f"{temperature:g} °C={int(count)}"
            for temperature, count in counts.fillna(0).items()
        )
        raise ValueError(
            f"each temperature must contain {EXPECTED_REPLICATES} rows; "
            f"observed {observed}"
        )

    if (data.loc[:, EXPRESSION_COLUMNS] < 0).any().any():
        raise ValueError("TPM-derived expression columns must be non-negative")
    proportions = data.loc[:, PROPORTION_COLUMNS]
    if ((proportions < 0) | (proportions > 100)).any().any():
        raise ValueError("proportion columns must lie within 0–100%")

    component_sum = data["α"] + data["β"]
    if not np.allclose(
        data["total"].to_numpy(dtype=float),
        component_sum.to_numpy(dtype=float),
        rtol=1e-6,
        atol=0.02,
    ):
        maximum = float(np.max(np.abs(data["total"] - component_sum)))
        raise ValueError(
            "total is inconsistent with α + β beyond tolerance "
            f"(maximum absolute difference: {maximum:.6g})"
        )

    proportion_difference = np.abs(
        (data["α/%"] + data["β/%"]).to_numpy(dtype=float) - 100.0
    )
    tolerance = PROPORTION_SUM_TOLERANCE + np.spacing(100.0)
    if np.any(proportion_difference > tolerance):
        maximum = float(np.max(proportion_difference))
        raise ValueError(
            "α/% + β/% is inconsistent with 100% beyond tolerance "
            f"(maximum absolute difference: {maximum:.6g})"
        )

    data["_input_order"] = np.arange(len(data))
    data = data.sort_values(
        ["Temperature", "_input_order"], kind="stable"
    ).reset_index(drop=True)
    for column in EXPRESSION_COLUMNS:
        data[f"{column}_ln"] = np.log1p(data[column].to_numpy(dtype=float))
    return data


def p_value_display(p_value: float) -> str:
    """Format a P value without ever presenting a numerical zero."""
    if p_value < 0.0001:
        return "<0.0001"
    if p_value < 0.001:
        return f"{p_value:.4f}"
    return f"{p_value:.3f}"


def p_value_statement(p_value: float) -> str:
    """Format a complete P-value statement for a panel title."""
    display = p_value_display(p_value)
    return f"P {display}" if display.startswith("<") else f"P = {display}"


def significance_label(p_value: float) -> str:
    """Return the conventional fixed-threshold significance label."""
    if p_value < 0.0001:
        return "****"
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return "ns"


def analyze_metrics(
    data: pd.DataFrame, *, alpha: float, base_seed: int
) -> tuple[list[MetricAnalysis], pd.DataFrame]:
    """Run ANOVA, Dunnett control comparisons, and Spearman tests per metric."""
    analyses: list[MetricAnalysis] = []
    rows: list[dict[str, object]] = []
    temperature_values = data["Temperature"].to_numpy(dtype=float)

    for metric_index, spec in enumerate(METRIC_SPECS):
        groups = {
            float(temperature): data.loc[
                data["Temperature"] == temperature, spec.analysis_column
            ].to_numpy(dtype=float)
            for temperature in EXPECTED_TEMPERATURES
        }
        ordered_groups = [groups[float(value)] for value in EXPECTED_TEMPERATURES]
        anova = stats.f_oneway(*ordered_groups)
        if not np.isfinite(anova.statistic) or not np.isfinite(anova.pvalue):
            raise ValueError(
                f"ANOVA is undefined for {spec.display_label}; check within-group variation"
            )

        spearman = stats.spearmanr(
            temperature_values,
            data[spec.analysis_column].to_numpy(dtype=float),
            alternative="two-sided",
        )
        spearman_rho = float(spearman.statistic)
        spearman_p = float(spearman.pvalue)
        spearman_n = len(data)
        if (
            not np.isfinite(spearman_rho)
            or not np.isfinite(spearman_p)
            or not -1.0 <= spearman_rho <= 1.0
            or not 0.0 <= spearman_p <= 1.0
        ):
            raise ValueError(
                f"Spearman correlation is undefined for {spec.display_label}"
            )

        dunnett_seed = base_seed + metric_index
        treatment_groups = [groups[value] for value in TREATMENT_TEMPERATURES]
        dunnett = stats.dunnett(
            *treatment_groups,
            control=groups[CONTROL_TEMPERATURE],
            alternative="two-sided",
            random_state=dunnett_seed,
        )
        dunnett_statistics = {
            temperature: float(statistic)
            for temperature, statistic in zip(
                TREATMENT_TEMPERATURES, dunnett.statistic, strict=True
            )
        }
        dunnett_pvalues = {
            temperature: float(p_value)
            for temperature, p_value in zip(
                TREATMENT_TEMPERATURES, dunnett.pvalue, strict=True
            )
        }
        pvalues = np.fromiter(dunnett_pvalues.values(), dtype=float)
        if not np.isfinite(pvalues).all() or np.any((pvalues < 0) | (pvalues > 1)):
            raise RuntimeError(f"Dunnett returned invalid P values for {spec.display_label}")

        analysis = MetricAnalysis(
            spec=spec,
            groups=groups,
            anova_f=float(anova.statistic),
            anova_p=float(anova.pvalue),
            spearman_rho=spearman_rho,
            spearman_p=spearman_p,
            spearman_n=spearman_n,
            dunnett_statistics=dunnett_statistics,
            dunnett_pvalues=dunnett_pvalues,
            dunnett_seed=dunnett_seed,
        )
        analyses.append(analysis)

        control = groups[CONTROL_TEMPERATURE]
        for temperature in TREATMENT_TEMPERATURES:
            treatment = groups[temperature]
            adjusted_p = dunnett_pvalues[temperature]
            rows.append(
                {
                    "metric": spec.key,
                    "display_label": spec.display_label,
                    "source_column": spec.source_column,
                    "analysis_scale": spec.analysis_scale,
                    "treatment_temperature_c": temperature,
                    "control_temperature_c": CONTROL_TEMPERATURE,
                    "n_treatment": len(treatment),
                    "n_control": len(control),
                    "treatment_mean": float(np.mean(treatment)),
                    "control_mean": float(np.mean(control)),
                    "mean_difference": float(np.mean(treatment) - np.mean(control)),
                    "anova_df_between": len(EXPECTED_TEMPERATURES) - 1,
                    "anova_df_within": len(data) - len(EXPECTED_TEMPERATURES),
                    "anova_f": float(anova.statistic),
                    "anova_p": float(anova.pvalue),
                    "correlation_test": "two-sided Spearman",
                    "spearman_n": spearman_n,
                    "spearman_rho": spearman_rho,
                    "spearman_p": spearman_p,
                    "spearman_p_display": p_value_display(spearman_p),
                    "comparison_test": "two-sided Dunnett vs 24 C",
                    "multiplicity_scope": "7 comparisons within metric",
                    "dunnett_statistic": dunnett_statistics[temperature],
                    "adjusted_p": adjusted_p,
                    "adjusted_p_display": p_value_display(adjusted_p),
                    "significance": significance_label(adjusted_p),
                    "alpha": alpha,
                    "reject_at_alpha": bool(adjusted_p < alpha),
                    "dunnett_seed": dunnett_seed,
                }
            )

    statistics_table = pd.DataFrame.from_records(rows)
    expected_rows = len(METRIC_SPECS) * len(TREATMENT_TEMPERATURES)
    if len(statistics_table) != expected_rows:
        raise RuntimeError(
            f"statistics table must contain {expected_rows} rows; "
            f"observed {len(statistics_table)}"
        )
    return analyses, statistics_table


def shared_axis_limits(
    analyses: Sequence[MetricAnalysis], *, proportion: bool
) -> tuple[float, float]:
    """Derive shared limits containing raw points, mean ± SD, and labels."""
    raw_values = [
        group for analysis in analyses for group in analysis.groups.values()
    ]
    uncertainty_extents = [
        extent
        for analysis in analyses
        for group in analysis.groups.values()
        for extent in (
            np.mean(group) - np.std(group, ddof=1),
            np.mean(group) + np.std(group, ddof=1),
        )
    ]
    values = np.concatenate(
        [*raw_values, np.asarray(uncertainty_extents, dtype=float)]
    )
    lower = float(np.min(values))
    upper = float(np.max(values))
    span = upper - lower
    if span <= 0:
        raise ValueError("cannot derive plot limits from a constant response")

    if proportion:
        step = 2.0 if span <= 16.0 else 5.0
        lower_limit = max(0.0, float(np.floor((lower - 0.08 * span) / step) * step))
        upper_limit = min(
            100.0, float(np.ceil((upper + 0.35 * span) / step) * step)
        )
    else:
        step = 0.2
        lower_limit = float(np.floor((lower - 0.07 * span) / step) * step)
        upper_limit = float(np.ceil((upper + 0.35 * span) / step) * step)

    if upper_limit <= lower_limit:
        raise RuntimeError("could not derive a valid shared y-axis range")
    return lower_limit, upper_limit


def style_axis(ax: Axes, *, proportion: bool) -> None:
    """Apply common categorical-temperature axis styling."""
    ax.axvspan(
        CONTROL_INDEX - 0.48,
        CONTROL_INDEX + 0.48,
        color=COLORS["control"],
        alpha=0.85,
        linewidth=0,
        zorder=0,
    )
    if proportion:
        ax.axhline(
            50.0,
            color=COLORS["reference"],
            linestyle=(0, (3, 2)),
            linewidth=0.65,
            zorder=1,
        )
    ax.set_xlim(-0.6, len(EXPECTED_TEMPERATURES) - 0.4)
    ax.set_xticks(BOX_POSITIONS)
    labels = [f"{value:g}" for value in EXPECTED_TEMPERATURES]
    ax.set_xticklabels(labels)
    ax.set_axisbelow(True)


def plot_metric(
    ax: Axes,
    analysis: MetricAnalysis,
    *,
    panel_id: str,
    y_limits: tuple[float, float],
) -> None:
    """Draw one metric's mean ± SD bars, replicates, and statistical labels."""
    spec = analysis.spec
    style_axis(ax, proportion=spec.is_proportion)
    ordered_groups = [
        analysis.groups[float(temperature)] for temperature in EXPECTED_TEMPERATURES
    ]
    means = np.asarray([np.mean(group) for group in ordered_groups], dtype=float)
    standard_deviations = np.asarray(
        [np.std(group, ddof=1) for group in ordered_groups], dtype=float
    )
    ax.bar(
        BOX_POSITIONS,
        means - y_limits[0],
        bottom=y_limits[0],
        width=0.58,
        color=spec.color,
        alpha=0.68,
        edgecolor=spec.color,
        linewidth=0.8,
        yerr=standard_deviations,
        capsize=2.4,
        error_kw={
            "ecolor": spec.color,
            "elinewidth": 0.85,
            "capthick": 0.85,
            "zorder": 3,
        },
        zorder=2,
    )

    for position, values in zip(BOX_POSITIONS, ordered_groups, strict=True):
        if len(values) != EXPECTED_REPLICATES:
            raise RuntimeError(
                f"expected {EXPECTED_REPLICATES} values at position {position:g}"
            )
        ax.scatter(
            position + POINT_OFFSETS,
            values,
            s=15,
            facecolor=spec.color,
            edgecolor="white",
            linewidth=0.45,
            alpha=0.92,
            zorder=4,
        )

    y_span = y_limits[1] - y_limits[0]
    for temperature in TREATMENT_TEMPERATURES:
        position = int(np.flatnonzero(EXPECTED_TEMPERATURES == temperature)[0])
        annotation_fraction = 0.82 if position % 2 == 0 else 0.90
        ax.text(
            position,
            y_limits[0] + annotation_fraction * y_span,
            significance_label(analysis.dunnett_pvalues[temperature]),
            ha="center",
            va="bottom",
            fontsize=5.8,
            fontweight="semibold",
            color="#242424",
            clip_on=False,
            zorder=5,
        )

    ax.set_title(spec.panel_title, loc="center", pad=4.5)
    offset = ScaledTranslation(-8.0 / 72.0, 3.0 / 72.0, ax.figure.dpi_scale_trans)
    label_artist = ax.text(
        0,
        1,
        f"({panel_id})",
        transform=ax.transAxes + offset,
        ha="right",
        va="bottom",
        fontsize=8.0,
        fontweight="bold",
        clip_on=False,
    )
    label_artist.set_gid("parenthesized-panel-label")


def validate_parenthesized_panel_labels(
    fig: Figure,
    axes: Sequence[Axes],
    *,
    panel_ids: Sequence[str],
    row_groups: Sequence[Sequence[str]],
    tolerance_pt: float = 1.5,
) -> None:
    """Validate fixed-point anchors for bold parenthesized panel labels."""
    if len(axes) != len(panel_ids):
        raise RuntimeError("panel-label validation received mismatched axes and ids")

    fig.canvas.draw()
    figure_width_pt = float(fig.get_size_inches()[0] * 72.0)
    anchors: dict[str, tuple[float, float]] = {}
    x_offsets: dict[str, float] = {}
    for ax, panel_id in zip(axes, panel_ids, strict=True):
        candidates = [
            artist
            for artist in ax.texts
            if artist.get_gid() == "parenthesized-panel-label"
        ]
        if len(candidates) != 1:
            raise RuntimeError(
                f"panel {panel_id} must contain exactly one parenthesized label"
            )
        artist = candidates[0]
        if artist.get_text() != f"({panel_id})":
            raise RuntimeError(f"panel {panel_id} label must be ({panel_id})")
        if str(artist.get_fontweight()).lower() != "bold":
            raise RuntimeError(f"panel {panel_id} label must be bold")

        display_xy = artist.get_transform().transform(artist.get_position())
        inches_xy = fig.dpi_scale_trans.inverted().transform(display_xy)
        anchor = (float(inches_xy[0] * 72.0), float(inches_xy[1] * 72.0))
        axes_left_pt = float(ax.get_position().x0 * figure_width_pt)
        anchors[panel_id] = anchor
        x_offsets[panel_id] = anchor[0] - axes_left_pt

    for row in row_groups:
        y_values = [anchors[panel_id][1] for panel_id in row]
        spread = float(max(y_values) - min(y_values))
        if spread > tolerance_pt:
            raise RuntimeError(
                "parenthesized panel-label y anchors differ by "
                f"{spread:.3f} pt; tolerance is {tolerance_pt:.3f} pt"
            )

    x_spread = float(max(x_offsets.values()) - min(x_offsets.values()))
    if x_spread > tolerance_pt:
        raise RuntimeError(
            "parenthesized panel-label x offsets differ by "
            f"{x_spread:.3f} pt; tolerance is {tolerance_pt:.3f} pt"
        )


def anova_figure_note(analyses: Sequence[MetricAnalysis]) -> str:
    """Return a compact two-line figure-level ANOVA annotation."""
    entries = [
        f"({panel_id}) F(7,16)={analysis.anova_f:.2f}, "
        f"{p_value_statement(analysis.anova_p)}"
        for panel_id, analysis in zip(PANEL_IDS, analyses, strict=True)
    ]
    return "One-way ANOVA: " + "; ".join(entries[:3]) + "\n" + "; ".join(entries[3:])


def spearman_figure_note(analyses: Sequence[MetricAnalysis]) -> str:
    """Return a compact figure-level temperature-correlation annotation."""
    sample_sizes = {analysis.spearman_n for analysis in analyses}
    if len(sample_sizes) != 1:
        raise RuntimeError("Spearman sample sizes must match across all five metrics")
    sample_size = sample_sizes.pop()
    entries = [
        f"({panel_id}) ρ={analysis.spearman_rho:.2f}, "
        f"{p_value_statement(analysis.spearman_p)}"
        for panel_id, analysis in zip(PANEL_IDS, analyses, strict=True)
    ]
    return (
        f"Spearman correlation with temperature (n = {sample_size}): "
        + "; ".join(entries)
    )


def build_figure(analyses: Sequence[MetricAnalysis]) -> tuple[Figure, list[Axes]]:
    """Build the 3+2 quantitative grid."""
    if len(analyses) != len(METRIC_SPECS):
        raise ValueError("five metric analyses are required to build the figure")

    configure_style()
    fig = plt.figure(figsize=(FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN))
    grid = fig.add_gridspec(
        2,
        6,
        left=0.075,
        right=0.99,
        bottom=0.115,
        top=0.70,
        wspace=0.92,
        hspace=0.78,
    )
    axes = [
        fig.add_subplot(grid[0, 0:2]),
        fig.add_subplot(grid[0, 2:4]),
        fig.add_subplot(grid[0, 4:6]),
        fig.add_subplot(grid[1, 1:3]),
        fig.add_subplot(grid[1, 3:5]),
    ]

    expression_limits = shared_axis_limits(analyses[:3], proportion=False)
    proportion_limits = shared_axis_limits(analyses[3:], proportion=True)
    for ax, analysis, panel_id in zip(axes, analyses, PANEL_IDS, strict=True):
        limits = proportion_limits if analysis.spec.is_proportion else expression_limits
        plot_metric(ax, analysis, panel_id=panel_id, y_limits=limits)
        if analysis.spec.is_proportion:
            ax.set_ylim(*proportion_limits)
            step = 2.0 if proportion_limits[1] - proportion_limits[0] <= 16 else 5.0
            ax.set_yticks(
                np.arange(proportion_limits[0], proportion_limits[1] + 0.5 * step, step)
            )
            ax.yaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))
        else:
            ax.set_ylim(*expression_limits)

    axes[0].set_ylabel("ln(TPM + 1)")
    axes[3].set_ylabel(r"Proportion of $\mathit{RCA}$ (%)")
    for ax in (axes[1], axes[2], axes[4]):
        ax.tick_params(axis="y", labelleft=False)

    fig.suptitle(
        r"Temperature-dependent $\mathit{RCA}$ expression and α/β composition",
        x=0.5,
        y=0.985,
        fontsize=7.2,
        fontweight="normal",
    )
    fig.text(
        0.5,
        0.93,
        "Bars = mean; error bars = SD; points = independent biological "
        "replicates (n = 3); stars = two-sided Dunnett vs 24 °C",
        ha="center",
        va="center",
        fontsize=6.2,
    )
    fig.text(
        0.5,
        0.865,
        anova_figure_note(analyses),
        ha="center",
        va="center",
        fontsize=5.8,
        linespacing=1.4,
    )
    fig.text(
        0.5,
        0.785,
        spearman_figure_note(analyses),
        ha="center",
        va="center",
        fontsize=5.8,
    )
    fig.supxlabel("Treatment temperature (°C)", fontsize=6.8, y=0.025)
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


def run_collision_audit(pdf_path: Path, json_path: Path) -> None:
    """Require a strict rendered-PDF collision audit before publishing outputs."""
    audit_script = qa_scripts_directory() / "audit_figure_collisions.py"
    if not audit_script.is_file():
        _skip_optional_qa(f"Optional QA tool not found: {audit_script}", json_path)
        return
    completed = subprocess.run(
        [
            sys.executable,
            str(audit_script),
            str(pdf_path),
            "--json-out",
            str(json_path),
            "--strict",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        details = (completed.stdout + "\n" + completed.stderr).strip()
        raise RuntimeError(
            "strict rendered-PDF collision audit failed"
            + (f": {details}" if details else "")
        )
    try:
        report = json.loads(json_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"could not read collision-audit JSON: {exc}") from exc
    if report.get("verdict") != "PASS":
        raise RuntimeError(
            f"collision audit did not pass (verdict: {report.get('verdict', 'unknown')})"
        )


def write_statistics(table: pd.DataFrame, path: Path) -> None:
    table.to_csv(path, sep="\t", index=False, encoding="utf-8", float_format="%.12g")


def export_outputs(
    fig: Figure,
    axes: Sequence[Axes],
    statistics_table: pd.DataFrame,
    *,
    output_dir: Path,
    basename: str,
    dpi: int,
) -> list[Path]:
    """Stage, audit, and atomically publish all requested output files."""
    output_dir.mkdir(parents=True, exist_ok=True)
    suffixes = (".pdf", ".png", ".statistics.tsv", ".alignment.json", ".collision-audit.json")
    final_paths = [output_dir / f"{basename}{suffix}" for suffix in suffixes]

    with tempfile.TemporaryDirectory(prefix=f".{basename}.", dir=output_dir) as temp_name:
        temp_dir = Path(temp_name)
        temp_paths = [temp_dir / path.name for path in final_paths]
        pdf_path, png_path, statistics_path, alignment_path, collision_path = temp_paths

        row_groups = [list(PANEL_IDS[:3]), list(PANEL_IDS[3:])]
        validate_parenthesized_panel_labels(
            fig,
            axes,
            panel_ids=PANEL_IDS,
            row_groups=row_groups,
            tolerance_pt=1.5,
        )
        require_matplotlib_panel_alignment = load_alignment_helper()
        require_matplotlib_panel_alignment(
            fig,
            axes=list(axes),
            panel_ids=list(PANEL_IDS),
            row_groups=row_groups,
            json_out=alignment_path,
            tolerance_pt=1.5,
            gutter_tolerance_pt=1.5,
            require_panel_labels=False,
            strict=True,
        )
        write_statistics(statistics_table, statistics_path)
        fig.savefig(pdf_path, format="pdf")
        fig.savefig(png_path, format="png", dpi=dpi)
        run_collision_audit(pdf_path, collision_path)

        collision_report = json.loads(collision_path.read_text(encoding="utf-8"))
        collision_report["pdf"] = str(final_paths[0])
        collision_path.write_text(
            json.dumps(collision_report, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )

        for staged, final in zip(temp_paths, final_paths, strict=True):
            os.replace(staged, final)
    return final_paths


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        validate_args(args)
        data = load_and_validate_data(args.input)
        analyses, statistics_table = analyze_metrics(
            data, alpha=args.alpha, base_seed=args.seed
        )
        fig, axes = build_figure(analyses)
        try:
            outputs = export_outputs(
                fig,
                axes,
                statistics_table,
                output_dir=args.output_dir,
                basename=args.basename,
                dpi=args.dpi,
            )
        finally:
            plt.close(fig)
    except (
        FileNotFoundError,
        OSError,
        RuntimeError,
        subprocess.SubprocessError,
        TypeError,
        ValueError,
    ) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    significant = int(statistics_table["reject_at_alpha"].sum())
    print(
        f"Validated {len(data)} samples across {len(EXPECTED_TEMPERATURES)} "
        f"temperatures ({EXPECTED_REPLICATES} replicates per temperature)."
    )
    print(
        f"Completed {len(analyses)} one-way ANOVAs and {len(statistics_table)} "
        f"Dunnett comparisons; {significant} comparisons have adjusted P < {args.alpha:g}."
    )
    print(
        f"Completed {len(analyses)} two-sided Spearman correlations using all "
        f"{len(data)} samples."
    )
    for output in outputs:
        print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
