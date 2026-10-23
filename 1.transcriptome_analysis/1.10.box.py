#!/usr/bin/env python3
"""Plot five RCA climate-temperature boxplots with rank-based statistics."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from typing import Sequence


# Matplotlib/fontconfig cannot write to the normal user cache in some batch jobs.
os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rca-climate-boxplot-mpl")
)
os.environ.setdefault(
    "XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "rca-climate-boxplot-cache")
)

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.colors import to_rgba
from matplotlib.figure import Figure
from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, TextArea
from matplotlib.ticker import MaxNLocator, PercentFormatter
from matplotlib.transforms import ScaledTranslation
import numpy as np
import pandas as pd
from scipy import stats


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[2]
DEFAULT_WORKDIR = PROJECT_ROOT / "list" / "RCA"
DEFAULT_INPUT = DEFAULT_WORKDIR / "RCA.climate.tsv"
DEFAULT_OUTPUT_DIR = DEFAULT_WORKDIR
DEFAULT_BASENAME = "RCA_climate_temperature_boxplot"

EXPECTED_TEMPERATURES = (10.0, 16.0, 22.0)
PAIRWISE_COMPARISONS = ((10.0, 16.0), (16.0, 22.0), (10.0, 22.0))
REQUIRED_COLUMNS = ("Temperature", "total", "α", "β", "α/%", "β/%")
EXPRESSION_COLUMNS = ("total", "α", "β")
PROPORTION_COLUMNS = ("α/%", "β/%")
PANEL_IDS = ("a", "b", "c", "d", "e")

FIGURE_WIDTH_IN = 183.0 / 25.4
FIGURE_HEIGHT_IN = 125.0 / 25.4
TOTAL_COMPONENT_ATOL = 0.10
TOTAL_COMPONENT_RTOL = 5e-5
PROPORTION_SUM_ATOL = 0.03
BOX_LINEWIDTH = 0.9

COLORS = {
    "total": "#4D4D4D",
    "alpha": "#2F6DAE",
    "beta": "#D9822B",
}


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
    """Configuration for one response panel."""

    key: str
    source_column: str
    analysis_column: str
    title: str
    ylabel: str
    scale: str
    color: str
    proportion: bool


@dataclass(frozen=True)
class PairwiseResult:
    """One two-sided Dunn comparison."""

    group1: float
    group2: float
    z_statistic: float
    p_raw: float
    p_adjusted: float


@dataclass(frozen=True)
class MetricAnalysis:
    """Rank-based analysis results and plotted groups for one metric."""

    spec: MetricSpec
    groups: dict[float, np.ndarray]
    kruskal_h: float
    kruskal_p: float
    pairwise: dict[tuple[float, float], PairwiseResult]


METRIC_SPECS = (
    MetricSpec(
        key="total",
        source_column="total",
        analysis_column="total_ln",
        title="Total RCA expression",
        ylabel="ln(TPM + 1)",
        scale="ln(TPM + 1)",
        color=COLORS["total"],
        proportion=False,
    ),
    MetricSpec(
        key="alpha",
        source_column="α",
        analysis_column="alpha_ln",
        title="RCA α expression",
        ylabel="ln(TPM + 1)",
        scale="ln(TPM + 1)",
        color=COLORS["alpha"],
        proportion=False,
    ),
    MetricSpec(
        key="beta",
        source_column="β",
        analysis_column="beta_ln",
        title="RCA β expression",
        ylabel="ln(TPM + 1)",
        scale="ln(TPM + 1)",
        color=COLORS["beta"],
        proportion=False,
    ),
    MetricSpec(
        key="alpha_percent",
        source_column="α/%",
        analysis_column="α/%",
        title="RCA α proportion",
        ylabel="Proportion (%)",
        scale="percentage",
        color=COLORS["alpha"],
        proportion=True,
    ),
    MetricSpec(
        key="beta_percent",
        source_column="β/%",
        analysis_column="β/%",
        title="RCA β proportion",
        ylabel="Proportion (%)",
        scale="percentage",
        color=COLORS["beta"],
        proportion=True,
    ),
)


def configure_style() -> None:
    """Apply a compact publication style with editable vector text."""
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "mathtext.fontset": "dejavusans",
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
            "Plot a 3+2 RCA temperature boxplot figure and run two-sided "
            "Kruskal-Wallis/Dunn tests without normality or variance tests."
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
        help="Significance threshold recorded in the statistics TSV (default: 0.05)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260905,
        help="Retained for command-line compatibility; no random marks are drawn",
    )
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> None:
    if args.dpi < 300:
        raise ValueError("--dpi must be at least 300")
    if not 0.0 < args.alpha < 1.0:
        raise ValueError("--alpha must be strictly between 0 and 1")
    if args.seed < 0:
        raise ValueError("--seed must be a non-negative integer")
    if not args.basename or Path(args.basename).name != args.basename:
        raise ValueError("--basename must be a filename stem, not a path")
    if Path(args.basename).suffix:
        raise ValueError("--basename must not include an extension")


def _validate_approximate_identity(
    observed: np.ndarray,
    expected: np.ndarray,
    *,
    label: str,
    atol: float,
    rtol: float,
) -> None:
    matches = np.isclose(observed, expected, atol=atol, rtol=rtol)
    if bool(np.all(matches)):
        return
    differences = np.abs(observed - expected)
    worst = int(np.argmax(differences))
    raise ValueError(
        f"{label} failed for row {worst + 2}: observed={observed[worst]:.12g}, "
        f"expected={expected[worst]:.12g}, absolute difference={differences[worst]:.12g}"
    )


def load_and_validate_data(path: Path) -> pd.DataFrame:
    """Load only the six analysis fields and validate their data contract."""
    if not path.is_file():
        raise FileNotFoundError(f"input TSV not found: {path}")

    raw = pd.read_csv(path, sep="\t", encoding="utf-8")
    missing = [column for column in REQUIRED_COLUMNS if column not in raw.columns]
    if missing:
        raise ValueError(f"input TSV is missing required columns: {', '.join(missing)}")

    data = raw.loc[:, REQUIRED_COLUMNS].copy()
    for column in REQUIRED_COLUMNS:
        try:
            data[column] = pd.to_numeric(data[column], errors="raise")
        except (TypeError, ValueError) as exc:
            raise ValueError(f"column {column} contains non-numeric values") from exc

    numeric = data.loc[:, REQUIRED_COLUMNS].to_numpy(dtype=float)
    if not np.isfinite(numeric).all():
        raise ValueError("required analysis columns contain missing or infinite values")

    observed_temperatures = tuple(sorted(data["Temperature"].unique().astype(float)))
    if observed_temperatures != EXPECTED_TEMPERATURES:
        raise ValueError(
            "Temperature must contain exactly 10, 16, and 22 °C; observed "
            + ", ".join(f"{value:g}" for value in observed_temperatures)
        )
    group_counts = data.groupby("Temperature", observed=False).size()
    if any(int(group_counts.get(temperature, 0)) == 0 for temperature in EXPECTED_TEMPERATURES):
        raise ValueError("each expected temperature group must contain observations")

    if (data.loc[:, EXPRESSION_COLUMNS] < 0).any().any():
        raise ValueError("total, α, and β expression values must be non-negative")
    proportions = data.loc[:, PROPORTION_COLUMNS]
    if ((proportions < 0) | (proportions > 100)).any().any():
        raise ValueError("α/% and β/% values must lie between 0 and 100")

    _validate_approximate_identity(
        data["total"].to_numpy(dtype=float),
        (data["α"] + data["β"]).to_numpy(dtype=float),
        label="total ≈ α + β",
        atol=TOTAL_COMPONENT_ATOL,
        rtol=TOTAL_COMPONENT_RTOL,
    )
    _validate_approximate_identity(
        (data["α/%"] + data["β/%"]).to_numpy(dtype=float),
        np.full(len(data), 100.0),
        label="α/% + β/% ≈ 100",
        atol=PROPORTION_SUM_ATOL,
        rtol=0.0,
    )

    data["total_ln"] = np.log1p(data["total"])
    data["alpha_ln"] = np.log1p(data["α"])
    data["beta_ln"] = np.log1p(data["β"])
    return data


def dunn_pairwise(
    groups: dict[float, np.ndarray],
) -> dict[tuple[float, float], PairwiseResult]:
    """Calculate two-sided Dunn tests with within-metric Bonferroni correction."""
    ordered_values = [np.asarray(groups[t], dtype=float) for t in EXPECTED_TEMPERATURES]
    pooled = np.concatenate(ordered_values)
    ranks = stats.rankdata(pooled, method="average")
    tie_correction = float(stats.tiecorrect(ranks))
    if tie_correction <= 0.0:
        raise ValueError("Dunn test is undefined because all observations are tied")

    rank_groups: dict[float, np.ndarray] = {}
    start = 0
    for temperature, values in zip(EXPECTED_TEMPERATURES, ordered_values, strict=True):
        stop = start + len(values)
        rank_groups[temperature] = ranks[start:stop]
        start = stop

    n_total = len(pooled)
    base_variance = tie_correction * n_total * (n_total + 1.0) / 12.0
    number_of_comparisons = len(PAIRWISE_COMPARISONS)
    results: dict[tuple[float, float], PairwiseResult] = {}
    for group1, group2 in PAIRWISE_COMPARISONS:
        ranks1 = rank_groups[group1]
        ranks2 = rank_groups[group2]
        standard_error = np.sqrt(
            base_variance * (1.0 / len(ranks1) + 1.0 / len(ranks2))
        )
        z_statistic = float((np.mean(ranks1) - np.mean(ranks2)) / standard_error)
        p_raw = float(2.0 * stats.norm.sf(abs(z_statistic)))
        p_adjusted = float(min(1.0, p_raw * number_of_comparisons))
        results[(group1, group2)] = PairwiseResult(
            group1=group1,
            group2=group2,
            z_statistic=z_statistic,
            p_raw=p_raw,
            p_adjusted=p_adjusted,
        )
    return results


def analyze_metrics(
    data: pd.DataFrame, *, alpha: float
) -> tuple[list[MetricAnalysis], pd.DataFrame, pd.DataFrame]:
    """Create plot-ready groups, descriptive summaries, and rank-test tables."""
    analyses: list[MetricAnalysis] = []
    summary_rows: list[dict[str, object]] = []
    statistics_rows: list[dict[str, object]] = []

    for spec in METRIC_SPECS:
        groups = {
            temperature: data.loc[
                data["Temperature"] == temperature, spec.analysis_column
            ].to_numpy(dtype=float)
            for temperature in EXPECTED_TEMPERATURES
        }
        pooled = np.concatenate(list(groups.values()))
        if float(np.ptp(pooled)) == 0.0:
            raise ValueError(f"metric {spec.source_column} is constant; rank tests are undefined")

        kruskal = stats.kruskal(
            *(groups[temperature] for temperature in EXPECTED_TEMPERATURES),
            nan_policy="raise",
        )
        pairwise = dunn_pairwise(groups)
        analysis = MetricAnalysis(
            spec=spec,
            groups=groups,
            kruskal_h=float(kruskal.statistic),
            kruskal_p=float(kruskal.pvalue),
            pairwise=pairwise,
        )
        analyses.append(analysis)

        statistics_rows.append(
            {
                "metric": spec.key,
                "source_column": spec.source_column,
                "analysis_scale": spec.scale,
                "test": "Kruskal-Wallis",
                "group1": "all",
                "group2": "all",
                "n_group1": len(pooled),
                "n_group2": np.nan,
                "statistic": analysis.kruskal_h,
                "p_raw": analysis.kruskal_p,
                "p_adjusted": np.nan,
                "adjustment": "none",
                "alpha": alpha,
                "reject_at_alpha": analysis.kruskal_p < alpha,
                "significance": significance_label(analysis.kruskal_p),
            }
        )
        for comparison in PAIRWISE_COMPARISONS:
            result = pairwise[comparison]
            statistics_rows.append(
                {
                    "metric": spec.key,
                    "source_column": spec.source_column,
                    "analysis_scale": spec.scale,
                    "test": "Dunn two-sided",
                    "group1": f"{result.group1:g}",
                    "group2": f"{result.group2:g}",
                    "n_group1": len(groups[result.group1]),
                    "n_group2": len(groups[result.group2]),
                    "statistic": result.z_statistic,
                    "p_raw": result.p_raw,
                    "p_adjusted": result.p_adjusted,
                    "adjustment": "Bonferroni within metric (3 comparisons)",
                    "alpha": alpha,
                    "reject_at_alpha": result.p_adjusted < alpha,
                    "significance": significance_label(result.p_adjusted),
                }
            )

        for temperature in EXPECTED_TEMPERATURES:
            values = groups[temperature]
            q1, median, q3 = np.quantile(values, [0.25, 0.50, 0.75])
            summary_rows.append(
                {
                    "metric": spec.key,
                    "source_column": spec.source_column,
                    "analysis_scale": spec.scale,
                    "temperature_c": temperature,
                    "n": len(values),
                    "mean": float(np.mean(values)),
                    "sd": float(np.std(values, ddof=1)),
                    "median": float(median),
                    "q1": float(q1),
                    "q3": float(q3),
                    "min": float(np.min(values)),
                    "max": float(np.max(values)),
                }
            )

    return analyses, pd.DataFrame(summary_rows), pd.DataFrame(statistics_rows)


def significance_label(p_value: float) -> str:
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return "ns"


def _box_whisker_extent(values: np.ndarray) -> tuple[float, float]:
    """Return the observed endpoints used by a 1.5-IQR boxplot whisker."""
    q1, q3 = np.quantile(values, [0.25, 0.75])
    iqr = q3 - q1
    lower_candidates = values[values >= q1 - 1.5 * iqr]
    upper_candidates = values[values <= q3 + 1.5 * iqr]
    return float(np.min(lower_candidates)), float(np.max(upper_candidates))


def _analysis_limits(
    analysis: MetricAnalysis,
) -> tuple[float, float, tuple[float, float, float], float]:
    """Scale one panel to its box-and-whisker data plus annotation headroom."""
    extents = [
        _box_whisker_extent(analysis.groups[temperature])
        for temperature in EXPECTED_TEMPERATURES
    ]
    whisker_lower = min(extent[0] for extent in extents)
    whisker_upper = max(extent[1] for extent in extents)
    span = max(whisker_upper - whisker_lower, 1e-9)
    lower = max(0.0, whisker_lower - 0.12 * span)
    if analysis.spec.proportion:
        step = max(1.5, 0.12 * span)
    else:
        step = max(0.25, 0.12 * span)
    # The separated levels keep nested comparison labels unambiguous at the
    # final 183 x 125 mm output size.
    bracket_levels = tuple(
        whisker_upper + step * index for index in (1.0, 2.5, 4.0)
    )
    upper = bracket_levels[-1] + step * 1.25
    return lower, upper, bracket_levels, whisker_upper


def add_significance_bracket(
    ax: Axes,
    x1: float,
    x2: float,
    y: float,
    label: str,
    *,
    step: float,
) -> None:
    height = 0.12 * step
    midpoint = (x1 + x2) / 2.0
    half_gap = 0.17
    ax.plot(
        [x1, x1, midpoint - half_gap],
        [y, y + height, y + height],
        color="#333333",
        linewidth=0.65,
        clip_on=False,
        zorder=5,
    )
    ax.plot(
        [midpoint + half_gap, x2, x2],
        [y + height, y + height, y],
        color="#333333",
        linewidth=0.65,
        clip_on=False,
        zorder=5,
    )
    ax.text(
        midpoint,
        y + height,
        label,
        ha="center",
        va="center",
        fontsize=6.0,
        color="#222222",
        clip_on=False,
        zorder=6,
    )


def add_panel_label(ax: Axes, panel_id: str) -> None:
    offset = ScaledTranslation(-8.0 / 72.0, 3.0 / 72.0, ax.figure.dpi_scale_trans)
    artist = ax.text(
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
    artist.set_gid("parenthesized-panel-label")


def rca_title_box(text: str, *, fontsize: float) -> HPacker:
    """Build a title in which only RCA is bold italic."""
    if text.count("RCA") != 1:
        raise ValueError("an RCA title must contain exactly one RCA token")
    prefix, suffix = text.split("RCA")
    children: list[TextArea] = []
    common = {"fontsize": fontsize, "fontfamily": "sans-serif", "color": "#111111"}
    if prefix.strip():
        children.append(
            TextArea(
                prefix.strip(),
                textprops={**common, "fontweight": "semibold"},
            )
        )
    children.append(
        TextArea(
            "RCA",
            textprops={
                **common,
                "fontweight": "bold",
                "fontstyle": "italic",
            },
        )
    )
    if suffix.strip():
        children.append(
            TextArea(
                suffix.strip(),
                textprops={**common, "fontweight": "semibold"},
            )
        )
    return HPacker(children=children, align="baseline", pad=0, sep=2.0)


def add_axes_rca_title(ax: Axes, text: str) -> None:
    title = AnchoredOffsetbox(
        loc="lower center",
        child=rca_title_box(text, fontsize=6.8),
        pad=0,
        borderpad=0,
        frameon=False,
        bbox_to_anchor=(0.5, 1.045),
        bbox_transform=ax.transAxes,
    )
    ax.add_artist(title)


def add_figure_rca_title(fig: Figure, text: str) -> None:
    title = AnchoredOffsetbox(
        loc="upper center",
        child=rca_title_box(text, fontsize=7.6),
        pad=0,
        borderpad=0,
        frameon=False,
        bbox_to_anchor=(0.5, 0.982),
        bbox_transform=fig.transFigure,
    )
    fig.add_artist(title)


def plot_metric(
    ax: Axes,
    analysis: MetricAnalysis,
    *,
    panel_id: str,
    limits: tuple[float, float, tuple[float, float, float], float],
) -> None:
    lower, upper, bracket_levels, _data_upper = limits
    positions = np.arange(len(EXPECTED_TEMPERATURES), dtype=float)
    values_by_group = [analysis.groups[t] for t in EXPECTED_TEMPERATURES]

    ax.boxplot(
        values_by_group,
        positions=positions,
        widths=0.50,
        whis=1.5,
        showfliers=False,
        patch_artist=True,
        boxprops={
            "facecolor": to_rgba(analysis.spec.color, 0.24),
            "edgecolor": analysis.spec.color,
            "linewidth": BOX_LINEWIDTH,
        },
        whiskerprops={"color": analysis.spec.color, "linewidth": BOX_LINEWIDTH},
        capprops={"color": analysis.spec.color, "linewidth": BOX_LINEWIDTH},
        medianprops={"color": "#1F1F1F", "linewidth": BOX_LINEWIDTH},
        zorder=2,
    )

    position_for = {
        temperature: positions[index]
        for index, temperature in enumerate(EXPECTED_TEMPERATURES)
    }
    step = bracket_levels[1] - bracket_levels[0]
    for comparison, y in zip(PAIRWISE_COMPARISONS, bracket_levels, strict=True):
        result = analysis.pairwise[comparison]
        add_significance_bracket(
            ax,
            position_for[comparison[0]],
            position_for[comparison[1]],
            y,
            significance_label(result.p_adjusted),
            step=step,
        )

    ax.set_xticks(positions)
    ax.set_xticklabels([f"{temperature:g}" for temperature in EXPECTED_TEMPERATURES])
    ax.set_xlim(-0.55, len(EXPECTED_TEMPERATURES) - 0.45)
    ax.set_ylim(lower, upper)
    add_axes_rca_title(ax, analysis.spec.title)
    ax.tick_params(direction="out")
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    if analysis.spec.proportion:
        ax.yaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))
    add_panel_label(ax, panel_id)


def build_figure(analyses: Sequence[MetricAnalysis]) -> tuple[Figure, list[Axes]]:
    if len(analyses) != 5:
        raise ValueError("five metric analyses are required")

    configure_style()
    fig = plt.figure(figsize=(FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN))
    grid = fig.add_gridspec(
        2,
        6,
        left=0.085,
        right=0.985,
        bottom=0.135,
        top=0.825,
        wspace=0.86,
        hspace=0.86,
    )
    axes = [
        fig.add_subplot(grid[0, 0:2]),
        fig.add_subplot(grid[0, 2:4]),
        fig.add_subplot(grid[0, 4:6]),
        fig.add_subplot(grid[1, 1:3]),
        fig.add_subplot(grid[1, 3:5]),
    ]

    for ax, analysis, panel_id in zip(axes, analyses, PANEL_IDS, strict=True):
        plot_metric(
            ax,
            analysis,
            panel_id=panel_id,
            limits=_analysis_limits(analysis),
        )

    axes[0].set_ylabel("ln(TPM + 1)")
    axes[3].set_ylabel("Proportion (%)")

    counts = [len(analyses[0].groups[t]) for t in EXPECTED_TEMPERATURES]
    sample_size_note = "; ".join(
        f"{temperature:g} °C, n = {count}"
        for temperature, count in zip(EXPECTED_TEMPERATURES, counts, strict=True)
    )

    add_figure_rca_title(
        fig,
        "RCA expression and α/β composition across temperature groups",
    )
    fig.text(
        0.5,
        0.925,
        "Boxes show median and IQR; whiskers extend to 1.5×IQR; outlier symbols are omitted. "
        "Brackets show two-sided Dunn tests with within-metric Bonferroni correction.\n"
        f"Sample sizes: {sample_size_note}.",
        ha="center",
        va="center",
        fontsize=5.7,
        linespacing=1.35,
    )
    fig.supxlabel("Temperature (°C)", fontsize=6.8, y=0.035)
    return fig, axes


def validate_parenthesized_panel_labels(
    fig: Figure,
    axes: Sequence[Axes],
    *,
    panel_ids: Sequence[str],
    row_groups: Sequence[Sequence[str]],
    tolerance_pt: float = 1.5,
) -> None:
    """Check bold parenthesized labels and their fixed-point alignment."""
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
            raise RuntimeError(f"panel {panel_id} must contain one parenthesized label")
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
                f"parenthesized panel-label y anchors differ by {spread:.3f} pt"
            )
    x_spread = float(max(x_offsets.values()) - min(x_offsets.values()))
    if x_spread > tolerance_pt:
        raise RuntimeError(
            f"parenthesized panel-label x offsets differ by {x_spread:.3f} pt"
        )


def qa_scripts_directory() -> Path:
    return Path(
        os.environ.get(
            "NATURE_FIGURE_QA_DIR",
            Path.home() / ".codex" / "skills" / "nature-figure" / "scripts",
        )
    )


def load_alignment_helper():
    directory = qa_scripts_directory()
    helper = directory / "audit_panel_alignment.py"
    if not helper.is_file():
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
    report = json.loads(json_path.read_text(encoding="utf-8"))
    if report.get("verdict") != "PASS":
        raise RuntimeError(
            f"collision audit did not pass: {report.get('verdict', 'unknown')}"
        )


def write_tsv(table: pd.DataFrame, path: Path) -> None:
    table.to_csv(path, sep="\t", index=False, encoding="utf-8", float_format="%.12g")


def export_outputs(
    fig: Figure,
    axes: Sequence[Axes],
    summary_table: pd.DataFrame,
    statistics_table: pd.DataFrame,
    *,
    output_dir: Path,
    basename: str,
    dpi: int,
) -> list[Path]:
    """Audit staged files and publish the complete output bundle atomically."""
    output_dir.mkdir(parents=True, exist_ok=True)
    suffixes = (
        ".svg",
        ".pdf",
        ".png",
        ".summary.tsv",
        ".statistics.tsv",
        ".alignment.json",
        ".collision-audit.json",
    )
    final_paths = [output_dir / f"{basename}{suffix}" for suffix in suffixes]

    with tempfile.TemporaryDirectory(prefix=f".{basename}.", dir=output_dir) as name:
        temp_dir = Path(name)
        staged = [temp_dir / path.name for path in final_paths]
        (
            svg_path,
            pdf_path,
            png_path,
            summary_path,
            statistics_path,
            alignment_path,
            collision_path,
        ) = staged

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

        write_tsv(summary_table, summary_path)
        write_tsv(statistics_table, statistics_path)
        fig.savefig(svg_path, format="svg")
        fig.savefig(pdf_path, format="pdf")
        fig.savefig(png_path, format="png", dpi=dpi)
        run_collision_audit(pdf_path, collision_path)

        report = json.loads(collision_path.read_text(encoding="utf-8"))
        report["pdf"] = str(final_paths[1])
        collision_path.write_text(
            json.dumps(report, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        for source, destination in zip(staged, final_paths, strict=True):
            os.replace(source, destination)
    return final_paths


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        validate_args(args)
        data = load_and_validate_data(args.input)
        analyses, summary_table, statistics_table = analyze_metrics(
            data, alpha=args.alpha
        )
        fig, axes = build_figure(analyses)
        try:
            outputs = export_outputs(
                fig,
                axes,
                summary_table,
                statistics_table,
                output_dir=args.output_dir,
                basename=args.basename,
                dpi=args.dpi,
            )
        finally:
            plt.close(fig)
    except (
        FileNotFoundError,
        json.JSONDecodeError,
        OSError,
        RuntimeError,
        subprocess.SubprocessError,
        TypeError,
        ValueError,
    ) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    counts = data.groupby("Temperature", observed=False).size()
    count_text = ", ".join(
        f"{temperature:g} °C: n={int(counts[temperature])}"
        for temperature in EXPECTED_TEMPERATURES
    )
    dunn_rows = statistics_table[statistics_table["test"] == "Dunn two-sided"]
    significant = int(dunn_rows["reject_at_alpha"].sum())
    print(f"Validated {len(data)} observations ({count_text}).")
    print(
        f"Completed 5 Kruskal-Wallis tests and {len(dunn_rows)} two-sided Dunn "
        f"comparisons; {significant} adjusted comparisons have P < {args.alpha:g}."
    )
    for output in outputs:
        print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
