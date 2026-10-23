#!/usr/bin/env python3
"""Plot methylation features across expression tertiles.

The scatter panels retain every finite observation and show mean +/- sample SD.
Lines join group means across SHAP-ranked features only to describe visual trends;
they do not represent a continuous-feature model or an inferential test.
"""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

PLOT_CACHE_DIR = os.path.join(tempfile.gettempdir(), "methylation-scatter-box-cache")
os.makedirs(PLOT_CACHE_DIR, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", PLOT_CACHE_DIR)
os.environ.setdefault("XDG_CACHE_HOME", PLOT_CACHE_DIR)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
from scipy import stats


FEATURE_GROUPS = ("CHH", "CHG", "CG")
SCATTER_TOP_N = 10
EXPRESSION_ORDER = ("Low", "Medium", "High")
STAT_GROUP_ORDER = ("High", "Medium", "Low")
PAIRWISE_COMPARISONS = (("High", "Medium"), ("High", "Low"), ("Medium", "Low"))
PAIRWISE_LEVELS = {"High_vs_Medium": 1, "Medium_vs_Low": 2, "High_vs_Low": 3}
GROUP_COLORS = {"Low": "#4E6F9E", "Medium": "#D3A35D", "High": "#79A9B8"}
METHYLATION_AXIS_LABEL = "Methylation level (%)"
GROUP_MARKERS = {"Low": "o", "Medium": "s", "High": "D"}
SCATTER_OFFSETS = {"Low": -0.22, "Medium": 0.0, "High": 0.22}
RANDOM_SEED = 20260905
ALIGNMENT_COLUMNS = (
    "Temperature",
    "total",
    "α",
    "β",
    "β1",
    "β2",
    "α/%",
    "β/%",
    "β1/%",
    "β2/%",
)
TARGET_CONFIGS = (
    {"result_file": "total_TPM_import_results.txt", "expr_col": "total", "output_dir": "total_TPM", "title": "total_TPM"},
    {"result_file": "A_TPM_import_results.txt", "expr_col": "α", "output_dir": "A_TPM", "title": "A_TPM"},
    {"result_file": "B_TPM_import_results.txt", "expr_col": "β", "output_dir": "B_TPM", "title": "B_TPM"},
    {"result_file": "A_per_TPM_import_results.txt", "expr_col": "α/%", "output_dir": "A_per_TPM", "title": "A_per_TPM"},
    {"result_file": "B_per_TPM_import_results.txt", "expr_col": "β/%", "output_dir": "B_per_TPM", "title": "B_per_TPM"},
)
FEATURE_PATTERN = re.compile(r"^\s*((?:snp|CG|CHG|CHH)_[^:]+):\s*([-+0-9.eE]+)\s*$")
REPO_ROOT_MARKER = "Arabidopsis_thaliana"
PAIRWISE_COLUMNS = (
    "comparison",
    "group1",
    "group2",
    "test_method",
    "raw_p_value",
    "adj_p_value",
    "star",
    "feature_group",
    "feature_name",
)


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
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "axes.unicode_minus": False,
        }
    )


def resolve_project_root() -> Path:
    script_path = Path(__file__).resolve()
    for parent in script_path.parents:
        if parent.name == REPO_ROOT_MARKER:
            return parent
    raise RuntimeError(f"Unable to locate {REPO_ROOT_MARKER!r} above {script_path}")


def first_existing(paths: Iterable[Path]) -> Path | None:
    return next((path for path in paths if path.exists()), None)


def resolve_io_paths() -> Tuple[Path, Path]:
    project_root = resolve_project_root()
    rca_dir = first_existing(
        (
            project_root / "list" / "RCA",
            Path("/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana/list/RCA"),
            Path("/datapool/life-gongl/zesheng/Arabidopsis_thaliana/list/RCA"),
        )
    )
    result_dir = first_existing(
        (
            project_root / "list" / "xgboot" / "TPM_4.5_all",
            Path("/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana/list/xgboot/TPM_4.5_all"),
            Path("/datapool/life-gongl/zesheng/Arabidopsis_thaliana/list/xgboot/TPM_4.5_all"),
        )
    )
    if rca_dir is None or result_dir is None:
        raise FileNotFoundError("Could not resolve the RCA or import_results directory.")
    return rca_dir, result_dir


def require_columns(data_frame: pd.DataFrame, columns: Sequence[str], source: Path | str) -> None:
    missing = [column for column in columns if column not in data_frame.columns]
    if missing:
        raise KeyError(f"{source} is missing required columns: {', '.join(missing)}")


def load_numeric_table(tsv_path: Path) -> pd.DataFrame:
    data_frame = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    if data_frame.empty:
        raise ValueError(f"Input table is empty: {tsv_path}")
    return data_frame.apply(pd.to_numeric, errors="coerce")


def align_coverage_rows(
    methylation_df: pd.DataFrame,
    coverage_df: pd.DataFrame,
    methylation_path: Path | str = "methylation table",
    coverage_path: Path | str = "coverage table",
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, int]]:
    require_columns(methylation_df, ALIGNMENT_COLUMNS, methylation_path)
    require_columns(coverage_df, ALIGNMENT_COLUMNS, coverage_path)

    key_columns = list(ALIGNMENT_COLUMNS)
    if methylation_df[key_columns].isna().any().any():
        raise ValueError(f"{methylation_path} contains missing alignment-key values")
    if coverage_df[key_columns].isna().any().any():
        raise ValueError(f"{coverage_path} contains missing alignment-key values")

    methylation_duplicates = int(methylation_df.duplicated(key_columns, keep=False).sum())
    coverage_duplicates = int(coverage_df.duplicated(key_columns, keep=False).sum())
    if methylation_duplicates:
        raise ValueError(f"{methylation_path} has {methylation_duplicates} rows with duplicate alignment keys")
    if coverage_duplicates:
        raise ValueError(f"{coverage_path} has {coverage_duplicates} rows with duplicate alignment keys")

    left = methylation_df.reset_index(drop=True).reset_index(names="methylation_index")
    right_keys = coverage_df.reset_index(drop=True).reset_index(names="coverage_index")
    try:
        matched = left[["methylation_index", *key_columns]].merge(
            right_keys[["coverage_index", *key_columns]],
            how="left",
            on=key_columns,
            validate="one_to_one",
        )
    except pd.errors.MergeError as exc:
        raise ValueError(f"Coverage alignment is ambiguous: {exc}") from exc

    if len(matched) != len(methylation_df):
        raise ValueError("Coverage alignment changed the methylation row count")
    unmatched = int(matched["coverage_index"].isna().sum())
    if unmatched:
        raise ValueError(f"Coverage alignment failed for {unmatched} methylation rows")

    coverage_indices = matched.sort_values("methylation_index")["coverage_index"].astype(int).to_numpy()
    aligned_methylation = methylation_df.reset_index(drop=True).copy()
    aligned_coverage = coverage_df.iloc[coverage_indices].reset_index(drop=True).copy()
    if not np.allclose(
        aligned_methylation[key_columns].to_numpy(dtype=float),
        aligned_coverage[key_columns].to_numpy(dtype=float),
        rtol=0.0,
        atol=1e-10,
        equal_nan=False,
    ):
        raise ValueError("Coverage alignment keys do not match after reordering")

    report = {
        "methylation_rows": int(len(aligned_methylation)),
        "coverage_rows": int(len(coverage_df)),
        "matched_rows": int(len(coverage_indices)),
        "extra_coverage_rows": int(len(coverage_df) - len(np.unique(coverage_indices))),
    }
    return aligned_methylation, aligned_coverage, report


def parse_shap_features(result_path: Path) -> Dict[str, List[Tuple[str, float]]]:
    grouped: Dict[str, List[Tuple[str, float]]] = {group: [] for group in FEATURE_GROUPS}
    found_header = False
    in_feature_block = False
    seen: set[str] = set()

    with result_path.open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if "所有特征及其平均|SHAP|值" in line:
                found_header = True
                in_feature_block = True
                continue
            if not in_feature_block:
                continue
            if not line.strip():
                break
            match = FEATURE_PATTERN.match(line)
            if match is None:
                raise ValueError(f"Malformed SHAP feature line at {result_path}:{line_number}: {line.strip()!r}")
            feature_name, shap_text = match.groups()
            shap_value = float(shap_text)
            if not np.isfinite(shap_value):
                raise ValueError(f"Non-finite SHAP value for {feature_name} in {result_path}")
            if feature_name in seen:
                raise ValueError(f"Duplicate SHAP feature {feature_name!r} in {result_path}")
            seen.add(feature_name)
            feature_group = feature_name.split("_", 1)[0]
            if shap_value > 0 and feature_group in grouped:
                grouped[feature_group].append((feature_name, shap_value))

    if not found_header:
        raise ValueError(f"SHAP feature section was not found in {result_path}")
    missing_groups = [group for group in FEATURE_GROUPS if not grouped[group]]
    if missing_groups:
        raise ValueError(f"{result_path.name} has no positive-SHAP features for: {', '.join(missing_groups)}")
    for group in FEATURE_GROUPS:
        grouped[group].sort(key=lambda item: item[1], reverse=True)
    return grouped


def validate_feature_columns(
    grouped_features: Dict[str, List[Tuple[str, float]]],
    methylation_df: pd.DataFrame,
    coverage_df: pd.DataFrame,
    result_file: str,
) -> None:
    features = [feature for group in FEATURE_GROUPS for feature, _ in grouped_features[group]]
    missing_methylation = [feature for feature in features if feature not in methylation_df.columns]
    missing_coverage = [feature for feature in features if feature not in coverage_df.columns]
    if missing_methylation or missing_coverage:
        messages = []
        if missing_methylation:
            messages.append(f"methylation table: {', '.join(missing_methylation[:10])}")
        if missing_coverage:
            messages.append(f"coverage table: {', '.join(missing_coverage[:10])}")
        raise KeyError(f"Features from {result_file} are missing in " + "; ".join(messages))


def add_expression_groups(data_frame: pd.DataFrame, expr_col: str) -> Tuple[pd.DataFrame, float, float]:
    require_columns(data_frame, [expr_col], "methylation table")
    expression = pd.to_numeric(data_frame[expr_col], errors="coerce")
    eligible_mask = expression.notna() & np.isfinite(expression) & expression.gt(0)
    eligible = expression.loc[eligible_mask]
    if eligible.empty:
        raise ValueError(f"No finite positive samples were found for {expr_col}")

    transformed = eligible if "%" in expr_col else np.log1p(eligible)
    if not np.isfinite(transformed).all():
        raise ValueError(f"The grouping transform produced non-finite values for {expr_col}")
    lower_threshold = float(transformed.quantile(1 / 3))
    upper_threshold = float(transformed.quantile(2 / 3))
    if not lower_threshold < upper_threshold:
        raise ValueError(f"Expression tertile thresholds are not distinct for {expr_col}")

    grouped = data_frame.loc[eligible_mask].copy()
    grouped["expression_group"] = np.where(
        transformed < lower_threshold,
        "Low",
        np.where(transformed > upper_threshold, "High", "Medium"),
    )
    counts = grouped["expression_group"].value_counts()
    empty_groups = [group for group in EXPRESSION_ORDER if int(counts.get(group, 0)) == 0]
    if empty_groups:
        raise ValueError(f"Expression grouping produced empty groups for {expr_col}: {', '.join(empty_groups)}")
    return grouped, lower_threshold, upper_threshold


def bh_adjust(p_values: Sequence[float]) -> List[float]:
    if not p_values:
        return []
    pvals = np.asarray(p_values, dtype=float)
    if np.any(~np.isfinite(pvals)) or np.any((pvals < 0) | (pvals > 1)):
        raise ValueError("P values supplied to BH adjustment must be finite and between 0 and 1")
    order = np.argsort(pvals)
    ranked = pvals[order]
    adjusted_ranked = np.empty(len(pvals), dtype=float)
    cumulative = 1.0
    for index in range(len(pvals) - 1, -1, -1):
        cumulative = min(cumulative, ranked[index] * len(pvals) / (index + 1))
        adjusted_ranked[index] = cumulative
    adjusted = np.empty(len(pvals), dtype=float)
    adjusted[order] = np.clip(adjusted_ranked, 0.0, 1.0)
    return adjusted.tolist()


def significance_stars(p_value: float | None) -> str:
    if p_value is None or not np.isfinite(p_value):
        return ""
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return ""


def dunn_pairwise(feature_values: Dict[str, pd.Series]) -> List[Dict[str, object]]:
    arrays = [feature_values[group].dropna().to_numpy(dtype=float) for group in STAT_GROUP_ORDER]
    if any(len(values) == 0 for values in arrays):
        return []
    all_values = np.concatenate(arrays)
    if np.all(all_values == all_values[0]):
        raw_p_values = [1.0] * len(PAIRWISE_COMPARISONS)
    else:
        ranks = stats.rankdata(all_values)
        group_sizes = [len(values) for values in arrays]
        offsets = np.cumsum([0, *group_sizes])
        mean_ranks = [float(np.mean(ranks[offsets[index] : offsets[index + 1]])) for index in range(len(arrays))]
        n_total = len(all_values)
        tie_counts = np.unique(all_values, return_counts=True)[1]
        tie_term = float(np.sum(tie_counts**3 - tie_counts) / (12 * (n_total - 1))) if n_total > 1 else 0.0
        rank_variance = n_total * (n_total + 1) / 12 - tie_term
        if rank_variance <= 0:
            raw_p_values = [1.0] * len(PAIRWISE_COMPARISONS)
        else:
            raw_p_values = []
            for first, second in PAIRWISE_COMPARISONS:
                first_index = STAT_GROUP_ORDER.index(first)
                second_index = STAT_GROUP_ORDER.index(second)
                denominator = np.sqrt(
                    rank_variance * (1 / group_sizes[first_index] + 1 / group_sizes[second_index])
                )
                z_score = abs(mean_ranks[first_index] - mean_ranks[second_index]) / denominator
                raw_p_values.append(float(2 * stats.norm.sf(z_score)))

    adjusted = bh_adjust(raw_p_values)
    records: List[Dict[str, object]] = []
    for (first, second), raw_p, adjusted_p in zip(PAIRWISE_COMPARISONS, raw_p_values, adjusted):
        records.append(
            {
                "comparison": f"{first}_vs_{second}",
                "group1": first,
                "group2": second,
                "test_method": "Dunn",
                "raw_p_value": raw_p,
                "adj_p_value": adjusted_p,
                "star": significance_stars(adjusted_p),
            }
        )
    return records


def compute_coverage_weighted_scores(
    grouped_data: pd.DataFrame,
    grouped_coverage: pd.DataFrame,
    grouped_features: Dict[str, List[Tuple[str, float]]],
) -> pd.DataFrame:
    if len(grouped_data) != len(grouped_coverage):
        raise ValueError("Grouped methylation and coverage tables have different row counts")
    weighted = grouped_data[["expression_group"]].reset_index(drop=True).copy()
    for feature_group in FEATURE_GROUPS:
        feature_names = [feature for feature, _ in grouped_features[feature_group]]
        methylation_values = grouped_data[feature_names].to_numpy(dtype=float)
        coverage_values = grouped_coverage[feature_names].to_numpy(dtype=float)
        valid = np.isfinite(methylation_values) & np.isfinite(coverage_values) & (coverage_values > 0)
        weight_sums = np.where(valid, coverage_values, 0.0).sum(axis=1)
        numerators = np.where(valid, methylation_values * coverage_values, 0.0).sum(axis=1)
        scores = np.full(len(grouped_data), np.nan, dtype=float)
        usable = weight_sums > 0
        scores[usable] = numerators[usable] / weight_sums[usable]
        weighted[feature_group] = scores
    return weighted


def summarize_weighted_group(feature_group: str, feature_values: Dict[str, pd.Series]) -> Dict[str, object]:
    row: Dict[str, object] = {
        "feature_group": feature_group,
        "feature_name": feature_group,
        "overall_test": "Kruskal_Wallis",
    }
    for group in STAT_GROUP_ORDER:
        values = feature_values[group].dropna()
        row[f"n_{group.lower()}"] = int(len(values))
        row[f"mean_{group.lower()}"] = float(values.mean()) if not values.empty else np.nan
        row[f"median_{group.lower()}"] = float(values.median()) if not values.empty else np.nan
    return row


def run_group_level_statistics(weighted_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    feature_rows: List[Dict[str, object]] = []
    pending_pairwise: Dict[str, List[Dict[str, object]]] = {}
    tested_groups: List[str] = []
    overall_p_values: List[float] = []

    for feature_group in FEATURE_GROUPS:
        values_by_group = {
            group: weighted_df.loc[weighted_df["expression_group"] == group, feature_group].dropna()
            for group in STAT_GROUP_ORDER
        }
        row = summarize_weighted_group(feature_group, values_by_group)
        arrays = [values_by_group[group].to_numpy(dtype=float) for group in STAT_GROUP_ORDER]
        if any(len(values) == 0 for values in arrays):
            row.update({"overall_p_value": np.nan, "overall_adj_p_value": np.nan, "overall_significant": False})
            feature_rows.append(row)
            continue
        combined = np.concatenate(arrays)
        overall_p = 1.0 if np.all(combined == combined[0]) else float(stats.kruskal(*arrays).pvalue)
        row.update({"overall_p_value": overall_p, "overall_adj_p_value": np.nan, "overall_significant": False})
        feature_rows.append(row)
        tested_groups.append(feature_group)
        overall_p_values.append(overall_p)
        pending_pairwise[feature_group] = dunn_pairwise(values_by_group)

    feature_df = pd.DataFrame(feature_rows)
    for feature_group, adjusted_p in zip(tested_groups, bh_adjust(overall_p_values)):
        mask = feature_df["feature_name"].eq(feature_group)
        feature_df.loc[mask, "overall_adj_p_value"] = adjusted_p
        feature_df.loc[mask, "overall_significant"] = adjusted_p < 0.05

    pairwise_rows: List[Dict[str, object]] = []
    significant_groups = set(feature_df.loc[feature_df["overall_significant"].eq(True), "feature_name"])
    for feature_group in FEATURE_GROUPS:
        if feature_group not in significant_groups:
            continue
        for record in pending_pairwise.get(feature_group, []):
            record["feature_group"] = feature_group
            record["feature_name"] = feature_group
            pairwise_rows.append(record)
    pairwise_df = pd.DataFrame(pairwise_rows, columns=PAIRWISE_COLUMNS)
    return feature_df, pairwise_df


def dynamic_limits(values: Sequence[float], lower_anchor: float | None = None, padding_fraction: float = 0.06) -> Tuple[float, float]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        raise ValueError("Cannot determine axis limits because no finite values are available")
    lower = float(np.min(finite))
    upper = float(np.max(finite))
    if lower_anchor is not None:
        lower = min(lower, lower_anchor)
        upper = max(upper, lower_anchor)
    span = upper - lower
    if span <= 0:
        span = max(abs(upper), 1.0)
    padding = span * padding_fraction
    plotted_lower = lower - padding if lower < 0 or lower_anchor is None else lower
    return plotted_lower, upper + padding


def add_significance_bracket(
    axis: plt.Axes,
    x1: float,
    x2: float,
    y: float,
    height: float,
    label: str,
) -> None:
    axis.plot([x1, x1, x2, x2], [y, y + height, y + height, y], color="#343434", linewidth=0.8, clip_on=False)
    axis.annotate(
        label,
        xy=((x1 + x2) / 2, y + height),
        xytext=(0, 7),
        textcoords="offset points",
        ha="center",
        va="bottom",
        fontsize=8,
        fontweight="bold",
        color="#202020",
        annotation_clip=False,
    )


def plot_scatter_figure(
    grouped_features: Dict[str, List[Tuple[str, float]]],
    grouped_data: pd.DataFrame,
) -> Tuple[plt.Figure, List[plt.Axes]]:
    figure, axes_array = plt.subplots(3, 1, figsize=(7.2, 8.8), constrained_layout=False)
    axes = list(np.atleast_1d(axes_array))
    figure.subplots_adjust(left=0.09, right=0.985, top=0.95, bottom=0.07, hspace=0.40)
    rng = np.random.default_rng(RANDOM_SEED)

    for axis, feature_group in zip(axes, FEATURE_GROUPS):
        feature_items = grouped_features[feature_group][:SCATTER_TOP_N]
        if not feature_items:
            raise ValueError(f"No scatter features are available for {feature_group}")
        axis_extent_values: List[float] = []

        for expression_group in EXPRESSION_ORDER:
            subset = grouped_data.loc[grouped_data["expression_group"].eq(expression_group)]
            mean_x: List[float] = []
            means: List[float] = []
            standard_deviations: List[float] = []
            for feature_index, (feature_name, _shap_value) in enumerate(feature_items):
                values = pd.to_numeric(subset[feature_name], errors="coerce").to_numpy(dtype=float)
                values = values[np.isfinite(values)]
                if values.size == 0:
                    continue
                center_x = feature_index + SCATTER_OFFSETS[expression_group]
                jitter = rng.uniform(-0.055, 0.055, size=values.size)
                axis.scatter(
                    np.full(values.size, center_x) + jitter,
                    values,
                    s=8,
                    marker=GROUP_MARKERS[expression_group],
                    color=GROUP_COLORS[expression_group],
                    alpha=0.18,
                    edgecolors="none",
                    rasterized=True,
                    zorder=1,
                )
                mean_value = float(np.mean(values))
                sd_value = float(np.std(values, ddof=1)) if values.size > 1 else np.nan
                mean_x.append(center_x)
                means.append(mean_value)
                standard_deviations.append(sd_value)
                axis_extent_values.extend(values.tolist())
                axis_extent_values.append(mean_value)
                if np.isfinite(sd_value):
                    axis_extent_values.extend([mean_value - sd_value, mean_value + sd_value])

            if not means:
                raise ValueError(f"No finite {feature_group} values are available for {expression_group}")
            axis.plot(
                mean_x,
                means,
                color=GROUP_COLORS[expression_group],
                linewidth=1.2,
                alpha=0.92,
                zorder=3,
            )
            axis.errorbar(
                mean_x,
                means,
                yerr=standard_deviations,
                fmt=GROUP_MARKERS[expression_group],
                markersize=4.8,
                color=GROUP_COLORS[expression_group],
                markeredgecolor="white",
                markeredgewidth=0.55,
                elinewidth=1.0,
                capsize=2.4,
                capthick=1.0,
                zorder=4,
            )

        y_lower, y_upper = dynamic_limits(axis_extent_values, lower_anchor=0.0, padding_fraction=0.055)
        axis.set_ylim(y_lower, y_upper)
        axis.yaxis.set_major_locator(MaxNLocator(nbins=6))
        axis.set_xlim(-0.55, len(feature_items) - 0.45)
        axis.set_xticks(np.arange(len(feature_items)))
        axis.set_xticklabels(
            [feature for feature, _ in feature_items],
            rotation=0,
            rotation_mode="anchor",
            ha="center",
            fontsize=6.5,
        )
        axis.set_ylabel(METHYLATION_AXIS_LABEL)
        axis.set_title(f"{feature_group} features ({len(feature_items)})", loc="left", fontweight="bold")
        axis.grid(axis="y", color="#D9D9D9", linestyle="--", linewidth=0.6, alpha=0.75)
        axis.tick_params(axis="both", length=3, width=0.7, color="#555555")
        axis.set_axisbelow(True)

    legend_handles = [
        Line2D(
            [0],
            [0],
            color=GROUP_COLORS[group],
            marker=GROUP_MARKERS[group],
            markersize=5,
            linewidth=1.2,
            markeredgecolor="white",
            markeredgewidth=0.5,
            label=group,
        )
        for group in EXPRESSION_ORDER
    ]
    for axis in axes:
        axis.legend(
            handles=legend_handles, loc="lower right", bbox_to_anchor=(1.0, 1.01),
            ncol=3, fontsize=8, handlelength=1.5, columnspacing=1.2, borderaxespad=0,
        )
    return figure, axes


def plot_weighted_boxplot(
    weighted_df: pd.DataFrame,
    feature_stats: pd.DataFrame,
    pairwise_stats: pd.DataFrame,
    feature_group: str,
) -> Tuple[plt.Figure, plt.Axes]:
    figure, axis = plt.subplots(figsize=(7.2, 5.2), constrained_layout=False)
    figure.subplots_adjust(left=0.14, right=0.97, top=0.86, bottom=0.16)
    positions = np.arange(len(EXPRESSION_ORDER), dtype=float)
    values_by_group: List[np.ndarray] = []
    all_values: List[float] = []

    for position, expression_group in zip(positions, EXPRESSION_ORDER):
        values = weighted_df.loc[weighted_df["expression_group"].eq(expression_group), feature_group].dropna().to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        if values.size == 0:
            raise ValueError(f"No finite weighted {feature_group} values are available for {expression_group}")
        values_by_group.append(values)
        all_values.extend(values.tolist())

    boxplot = axis.boxplot(
        values_by_group,
        positions=positions,
        widths=0.46,
        patch_artist=True,
        showfliers=False,
        manage_ticks=False,
        medianprops={"color": "#1C1C1C", "linewidth": 1.3},
        whiskerprops={"color": "#555555", "linewidth": 1.0},
        capprops={"color": "#555555", "linewidth": 1.0},
        boxprops={"color": "#555555", "linewidth": 1.0},
        zorder=2,
    )
    for patch, expression_group in zip(boxplot["boxes"], EXPRESSION_ORDER):
        patch.set_facecolor(GROUP_COLORS[expression_group])
        patch.set_alpha(1.0)

    # Hidden outliers remain in the statistics. Scale the display to visible
    # Tukey whiskers and annotations now that sample points are not drawn.
    visible_values = np.concatenate([line.get_ydata() for line in boxplot["whiskers"]])
    data_min = float(np.min(visible_values))
    data_max = float(np.max(visible_values))
    data_span = max(data_max - min(0.0, data_min), 1.0)
    bracket_step = max(data_span * 0.13, 1.2)
    bracket_height = max(data_span * 0.018, 0.22)
    annotation_max = data_max
    overall_row = feature_stats.loc[feature_stats["feature_name"].eq(feature_group)]
    if overall_row.empty:
        raise ValueError(f"Overall statistics are missing for {feature_group}")
    overall_significant = bool(overall_row.iloc[0]["overall_significant"])

    if overall_significant:
        comparisons = pairwise_stats.loc[pairwise_stats["feature_name"].eq(feature_group)].copy()
        if len(comparisons) != len(PAIRWISE_COMPARISONS):
            raise ValueError(f"Expected three Dunn comparisons for significant group {feature_group}")
        comparisons["level"] = comparisons["comparison"].map(PAIRWISE_LEVELS)
        comparisons = comparisons.sort_values("level")
        position_by_group = {group: float(index) for index, group in enumerate(EXPRESSION_ORDER)}
        for row in comparisons.itertuples(index=False):
            y = data_max + bracket_step * row.level
            label = row.star if row.star else "ns"
            add_significance_bracket(axis, position_by_group[row.group1], position_by_group[row.group2], y, bracket_height, label)
            annotation_max = max(annotation_max, y + bracket_height + bracket_step * 0.18)
    else:
        annotation_y = data_max + bracket_step
        axis.text(
            1.0,
            annotation_y,
            "Overall: ns",
            ha="center",
            va="bottom",
            fontsize=8,
            color="#333333",
        )
        annotation_max = annotation_y + bracket_step * 0.35

    y_lower, y_upper = dynamic_limits([*visible_values, annotation_max], lower_anchor=0.0, padding_fraction=0.055)
    axis.set_ylim(y_lower, y_upper)
    axis.yaxis.set_major_locator(MaxNLocator(nbins=6))
    axis.set_xlim(-0.55, 2.55)
    axis.set_xticks(positions)
    axis.set_xticklabels(
        list(EXPRESSION_ORDER),
        fontsize=9,
    )
    sample_counts = "; ".join(
        f"{group}: n = {len(values)}" for group, values in zip(EXPRESSION_ORDER, values_by_group)
    )
    figure.text(0.14, 0.065, sample_counts, ha="left", va="bottom", fontsize=8)
    figure.text(
        0.14, 0.025,
        "n = samples with a valid coverage-weighted methylation value.",
        ha="left", va="bottom", fontsize=8,
    )
    axis.set_ylabel(METHYLATION_AXIS_LABEL)
    axis.set_title(feature_group, loc="left", fontweight="bold", pad=18)
    axis.grid(False)
    axis.tick_params(axis="both", length=3, width=0.7, color="#555555")
    axis.set_axisbelow(True)
    legend_handles = [
        Patch(facecolor=GROUP_COLORS[group], edgecolor="#555555", alpha=1.0, label=group)
        for group in EXPRESSION_ORDER
    ]
    axis.legend(handles=legend_handles, loc="lower right", bbox_to_anchor=(1.0, 1.075),
                ncol=3, fontsize=8, borderaxespad=0)
    return figure, axis


def write_statistics(
    output_dir: Path,
    prefix: str,
    feature_stats: pd.DataFrame,
    pairwise_stats: pd.DataFrame,
    group_counts: Dict[str, int],
    thresholds: Tuple[float, float],
) -> List[Path]:
    feature_path = output_dir / f"{prefix}_weighted_group_stats.csv"
    pairwise_path = output_dir / f"{prefix}_weighted_group_pairwise_stats.csv"
    summary_path = output_dir / f"{prefix}_weighted_group_stats_summary.txt"
    feature_stats.to_csv(feature_path, index=False)
    pairwise_stats.to_csv(pairwise_path, index=False)

    lines = [
        f"=== {prefix} coverage-weighted group statistics ===",
        "",
        "grouping\tpositive finite values; tertiles at 1/3 and 2/3",
        f"lower_threshold\t{thresholds[0]:.12g}",
        f"upper_threshold\t{thresholds[1]:.12g}",
        f"n_low\t{group_counts['Low']}",
        f"n_medium\t{group_counts['Medium']}",
        f"n_high\t{group_counts['High']}",
        "overall_test\tKruskal-Wallis",
        "overall_correction\tBenjamini-Hochberg across CHH, CHG, CG",
        "posthoc_test\tDunn after significant adjusted overall test",
        "posthoc_correction\tBenjamini-Hochberg within feature group",
        "",
        f"tested_groups\t{len(feature_stats)}",
        f"overall_significant_groups\t{int(feature_stats['overall_significant'].eq(True).sum())}",
        f"pairwise_significant_comparisons\t{int(pairwise_stats['star'].ne('').sum()) if not pairwise_stats.empty else 0}",
        "",
        "feature_group\toverall_adj_p_value\toverall_significant",
    ]
    for row in feature_stats.itertuples(index=False):
        lines.append(f"{row.feature_name}\t{row.overall_adj_p_value}\t{row.overall_significant}")
    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return [feature_path, pairwise_path, summary_path]


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
        details = "\n".join(part.strip() for part in (completed.stdout, completed.stderr) if part.strip())
        raise RuntimeError("Strict rendered-PDF collision audit failed" + (f":\n{details}" if details else ""))
    report = json.loads(json_path.read_text(encoding="utf-8"))
    if report.get("verdict") != "PASS":
        raise RuntimeError(f"Collision audit did not pass for {pdf_path.name}: {report.get('verdict', 'unknown')}")


def export_figure(
    figure: plt.Figure,
    axes: Sequence[plt.Axes],
    output_dir: Path,
    stem: str,
    require_matplotlib_panel_alignment,
) -> List[Path]:
    png_path = output_dir / f"{stem}.png"
    pdf_path = output_dir / f"{stem}.pdf"
    svg_path = output_dir / f"{stem}.svg"
    tiff_path = output_dir / f"{stem}.tiff"
    alignment_path = output_dir / f"{stem}.alignment.json"
    collision_path = output_dir / f"{stem}.collision-audit.json"
    text_path = output_dir / f"{stem}.text-audit.json"
    panel_ids = [chr(ord("a") + index) for index in range(len(axes))]
    require_matplotlib_panel_alignment(
        figure,
        axes=list(axes),
        panel_ids=panel_ids,
        json_out=alignment_path,
        tolerance_pt=1.5,
        gutter_tolerance_pt=1.5,
        require_panel_labels=False,
        strict=True,
    )
    figure.savefig(svg_path, format="svg", bbox_inches="tight", facecolor="white")
    figure.savefig(pdf_path, format="pdf", bbox_inches="tight", facecolor="white")
    figure.savefig(png_path, format="png", dpi=600, bbox_inches="tight", facecolor="white")
    figure.savefig(
        tiff_path,
        format="tiff",
        dpi=600,
        bbox_inches="tight",
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(figure)
    run_collision_audit(pdf_path, collision_path)
    text_script = qa_scripts_directory() / "audit_pdf_text.py"
    if not text_script.is_file():
        _skip_optional_qa(f"Optional QA tool not found: {text_script}", text_path)
    else:
        text_check = subprocess.run(
            [sys.executable, str(qa_scripts_directory() / "audit_pdf_text.py"),
             str(pdf_path), "--min-pt", "5", "--json"],
            capture_output=True, text=True, check=True,
        )
        text_path.write_text(text_check.stdout, encoding="utf-8")
    return [png_path, pdf_path, svg_path, tiff_path, alignment_path, collision_path, text_path]


def add_summary_statistics(
    summary: pd.DataFrame,
    feature_group: str,
    feature_stats: pd.DataFrame,
    pairwise_stats: pd.DataFrame,
) -> pd.DataFrame:
    """Repeat context-level tests on the three descriptive group rows.

    Significance uses BH-adjusted P: * < .05, ** < .01, *** < .001.
    NA denotes an unavailable/not-performed test, never non-significance.
    """
    result = summary.copy()
    overall = feature_stats.loc[feature_stats["feature_name"].eq(feature_group)].iloc[0]
    result["overall_test"] = overall["overall_test"]
    result["overall_p_value"] = overall["overall_p_value"]
    result["overall_adj_p_value"] = overall["overall_adj_p_value"]
    result["overall_correction"] = "BH across CHH, CHG, CG within target"
    p = float(overall["overall_adj_p_value"])
    result["overall_significance"] = (significance_stars(p) or "n.s.") if np.isfinite(p) else "NA"
    result["pairwise_test"] = "Dunn"
    result["pairwise_correction"] = "BH across three comparisons within feature group"
    pairs = pairwise_stats.loc[pairwise_stats["feature_name"].eq(feature_group)]
    for first, second in PAIRWISE_COMPARISONS:
        comparison = f"{first}_vs_{second}"
        selected = pairs.loc[pairs["comparison"].eq(comparison)]
        raw_p = float(selected.iloc[0]["raw_p_value"]) if len(selected) else np.nan
        adj_p = float(selected.iloc[0]["adj_p_value"]) if len(selected) else np.nan
        result[f"{comparison}_p_value"] = raw_p
        result[f"{comparison}_adj_p_value"] = adj_p
        result[f"{comparison}_significance"] = (
            (significance_stars(adj_p) or "n.s.") if np.isfinite(adj_p) else "NA"
        )
        result[f"{comparison}_status"] = "tested" if np.isfinite(adj_p) else "not_tested"
    return result


def write_plot_data(
    output_dir: Path,
    grouped_data: pd.DataFrame,
    weighted_df: pd.DataFrame,
    grouped_features: Dict[str, List[Tuple[str, float]]],
    config: Dict[str, str],
) -> List[Path]:
    """Export actual observations and plot summaries; missing values stay explicit.

    sample_row is the one-based data-row number in the filtered input TSV
    (excluding its header), not a biological sample identifier. All methylation
    measurements remain on their original percentage scale.
    """
    from matplotlib.cbook import boxplot_stats

    outputs: List[Path] = []
    feature_stats, pairwise_stats = run_group_level_statistics(weighted_df)

    def save(frame: pd.DataFrame, name: str) -> None:
        path = output_dir / name
        frame.to_csv(path, sep="\t", index=False, na_rep="NA")
        outputs.append(path)

    metadata = pd.DataFrame({
        "target": config["title"],
        "sample_row": grouped_data.index.to_numpy() + 1,
        "temperature": grouped_data["Temperature"].to_numpy(),
        "expression_column": config["expr_col"],
        "expression_value": grouped_data[config["expr_col"]].to_numpy(),
        "expression_group": grouped_data["expression_group"].to_numpy(),
    })
    scatter_frames = []
    scatter_summary = []
    for feature_group in FEATURE_GROUPS:
        for rank, (feature, shap) in enumerate(grouped_features[feature_group][:SCATTER_TOP_N], 1):
            frame = metadata.copy()
            frame["feature_group"] = feature_group
            frame["feature_name"] = feature
            frame["shap_rank"] = rank
            frame["mean_abs_shap"] = shap
            frame["methylation_percent"] = grouped_data[feature].to_numpy(dtype=float)
            frame["included_in_plot"] = np.isfinite(frame["methylation_percent"])
            scatter_frames.append(frame)
            for group in EXPRESSION_ORDER:
                subset = frame.loc[frame["expression_group"].eq(group)]
                values = subset.loc[subset["included_in_plot"], "methylation_percent"].to_numpy()
                mean = float(np.mean(values)) if len(values) else np.nan
                sd = float(np.std(values, ddof=1)) if len(values) > 1 else np.nan
                scatter_summary.append({
                    "target": config["title"], "feature_group": feature_group,
                    "feature_name": feature, "shap_rank": rank, "mean_abs_shap": shap,
                    "expression_group": group, "n_total": len(subset), "n": len(values),
                    "n_missing": len(subset) - len(values), "mean": mean, "sd": sd,
                    "mean_minus_sd": mean - sd, "mean_plus_sd": mean + sd,
                })
    save(pd.concat(scatter_frames, ignore_index=True), "me_contribution_scatter.source_data.tsv")
    save(pd.DataFrame(scatter_summary), "me_contribution_scatter.summary.tsv")

    for feature_group in FEATURE_GROUPS:
        stem = f"me_contribution_weighted_boxplot_{feature_group}"
        frame = metadata.copy()
        frame["feature_group"] = feature_group
        frame["weighted_methylation_percent"] = weighted_df[feature_group].to_numpy(dtype=float)
        frame["included_in_boxplot"] = np.isfinite(frame["weighted_methylation_percent"])
        save(frame, f"{stem}.source_data.tsv")
        summary = []
        for group in EXPRESSION_ORDER:
            subset = frame.loc[frame["expression_group"].eq(group)]
            values = subset.loc[subset["included_in_boxplot"], "weighted_methylation_percent"].to_numpy()
            box = boxplot_stats(values, whis=1.5)[0]
            summary.append({
                "target": config["title"], "feature_group": feature_group,
                "expression_group": group, "n_total": len(subset), "n": len(values),
                "n_missing": len(subset) - len(values), "mean": box["mean"],
                "sd": float(np.std(values, ddof=1)) if len(values) > 1 else np.nan,
                "q1": box["q1"], "median": box["med"], "q3": box["q3"],
                "iqr": box["iqr"], "whisker_low": box["whislo"],
                "whisker_high": box["whishi"], "n_hidden_outliers": len(box["fliers"]),
            })
        save(add_summary_statistics(pd.DataFrame(summary), feature_group, feature_stats, pairwise_stats),
             f"{stem}.summary.tsv")
    return outputs


def process_target(
    grouped_features: Dict[str, List[Tuple[str, float]]],
    methylation_df: pd.DataFrame,
    coverage_df: pd.DataFrame,
    config: Dict[str, str],
    staging_dir: Path,
    require_matplotlib_panel_alignment,
) -> List[Path]:
    grouped_data, lower_threshold, upper_threshold = add_expression_groups(methylation_df, config["expr_col"])
    grouped_coverage = coverage_df.loc[grouped_data.index].copy()
    group_counts = {
        group: int(grouped_data["expression_group"].eq(group).sum()) for group in EXPRESSION_ORDER
    }
    weighted_df = compute_coverage_weighted_scores(grouped_data, grouped_coverage, grouped_features)
    feature_stats, pairwise_stats = run_group_level_statistics(weighted_df)

    outputs = write_statistics(
        staging_dir,
        config["title"],
        feature_stats,
        pairwise_stats,
        group_counts,
        (lower_threshold, upper_threshold),
    )
    outputs.extend(write_plot_data(staging_dir, grouped_data, weighted_df, grouped_features, config))
    scatter_figure, scatter_axes = plot_scatter_figure(grouped_features, grouped_data)
    outputs.extend(
        export_figure(
            scatter_figure,
            scatter_axes,
            staging_dir,
            "me_contribution_scatter",
            require_matplotlib_panel_alignment,
        )
    )
    for feature_group in FEATURE_GROUPS:
        box_figure, box_axis = plot_weighted_boxplot(weighted_df, feature_stats, pairwise_stats, feature_group)
        outputs.extend(
            export_figure(
                box_figure,
                [box_axis],
                staging_dir,
                f"me_contribution_weighted_boxplot_{feature_group}",
                require_matplotlib_panel_alignment,
            )
        )
    print(
        f"[{config['title']}] groups="
        + ", ".join(f"{group}:{group_counts[group]}" for group in EXPRESSION_ORDER)
    )
    return outputs


def main() -> None:
    configure_style()
    rca_dir, result_dir = resolve_io_paths()
    methylation_path = rca_dir / "Alt.all.snp_meth.filtered.tsv"
    coverage_path = rca_dir / "Alt.weighted.snp_meth.tsv"
    if not methylation_path.is_file():
        raise FileNotFoundError(f"Missing methylation matrix: {methylation_path}")
    if not coverage_path.is_file():
        raise FileNotFoundError(f"Missing coverage matrix: {coverage_path}")

    methylation_raw = load_numeric_table(methylation_path)
    coverage_raw = load_numeric_table(coverage_path)
    methylation_df, coverage_df, alignment_report = align_coverage_rows(
        methylation_raw,
        coverage_raw,
        methylation_path,
        coverage_path,
    )
    print(
        "[coverage alignment] "
        + ", ".join(f"{key}={value}" for key, value in alignment_report.items())
    )

    parsed_targets = []
    for config in TARGET_CONFIGS:
        result_path = result_dir / config["result_file"]
        if not result_path.is_file():
            raise FileNotFoundError(f"Missing import_results file: {result_path}")
        grouped_features = parse_shap_features(result_path)
        validate_feature_columns(grouped_features, methylation_df, coverage_df, config["result_file"])
        parsed_targets.append((config, grouped_features))

    require_matplotlib_panel_alignment = load_alignment_helper()
    staged_to_final: List[Tuple[Path, Path]] = []
    with tempfile.TemporaryDirectory(prefix=".me_contribution_scatter_box.", dir=rca_dir) as temp_name:
        temp_root = Path(temp_name)
        for config, grouped_features in parsed_targets:
            staging_dir = temp_root / config["output_dir"]
            staging_dir.mkdir(parents=True)
            staged_paths = process_target(
                grouped_features,
                methylation_df,
                coverage_df,
                config,
                staging_dir,
                require_matplotlib_panel_alignment,
            )
            final_dir = rca_dir / config["output_dir"]
            for staged_path in staged_paths:
                staged_to_final.append((staged_path, final_dir / staged_path.name))

        for staged_path, final_path in staged_to_final:
            final_path.parent.mkdir(parents=True, exist_ok=True)
            if staged_path.name.endswith((".collision-audit.json", ".text-audit.json")):
                report = json.loads(staged_path.read_text(encoding="utf-8"))
                report["pdf"] = str(final_path.with_name(
                    final_path.name.replace(".collision-audit.json", ".pdf").replace(".text-audit.json", ".pdf")
                ))
                staged_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
            os.replace(staged_path, final_path)

    print(f"Published {len(staged_to_final)} files under {rca_dir}")


if __name__ == "__main__":
    main()
