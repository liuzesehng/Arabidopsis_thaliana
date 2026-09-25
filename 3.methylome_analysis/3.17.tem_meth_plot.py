#!/usr/bin/env python3
"""Plot temperature-associated methylation for SHAP-ranked RCA features.

Each expression target is exported as one six-panel figure.  The upper row
shows feature-wise means ± sample SD for 10, 16 and 22 degrees C on one
shared methylation axis.  Lines only connect SHAP-ranked feature means as a
descriptive visual guide.  The lower row shows coverage-weighted Top-10
methylation distributions calculated from every usable sample, and non-parametric
temperature comparisons.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Mapping, Sequence


PLOT_CACHE_DIR = Path(tempfile.gettempdir()) / "tem-meth-matplotlib-cache"
PLOT_CACHE_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(PLOT_CACHE_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(PLOT_CACHE_DIR))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
from matplotlib.transforms import ScaledTranslation
from scipy import stats


PROJECT_MARKER = "Arabidopsis_thaliana"
METHYLATION_FILENAME = "Alt.all.snp_meth.filtered.tsv"
COVERAGE_FILENAME = "Alt.weighted.snp_meth.tsv"
IMPORT_RESULTS_DIRNAME = "TPM_4.5_all"
OUTPUT_DIRNAME = "tem_meth"

FEATURE_GROUPS = ("CHH", "CHG", "CG")
TEMPERATURES = (10, 16, 22)
PAIRWISE_COMPARISONS = ((10, 16), (16, 22), (10, 22))
TEMPERATURE_COLORS = {10: "#4E6F9E", 16: "#D3A35D", 22: "#79A9B8"}
TEMPERATURE_MARKERS = {10: "o", 16: "s", 22: "D"}
TEMPERATURE_OFFSETS = {10: -0.20, 16: 0.0, 22: 0.20}

TARGET_CONFIGS = (
    {"name": "A_per_TPM", "result": "A_per_TPM_import_results.txt", "expression": "α/%"},
    {"name": "A_TPM", "result": "A_TPM_import_results.txt", "expression": "α"},
    {"name": "B_per_TPM", "result": "B_per_TPM_import_results.txt", "expression": "β/%"},
    {"name": "B_TPM", "result": "B_TPM_import_results.txt", "expression": "β"},
    {"name": "total_TPM", "result": "total_TPM_import_results.txt", "expression": "total"},
)
ABSOLUTE_EXPRESSION_COLUMNS = {"total", "α", "β"}
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
FEATURE_PATTERN = re.compile(
    r"^\s*((?:CHH|CHG|CG|snp)_[^:]+):\s*([-+0-9.eE]+)\s*$"
)

FIGURE_WIDTH_IN = 183.0 / 25.4
FIGURE_HEIGHT_IN = 180.0 / 25.4
MIN_FONT_PT = 5.0


def configure_style() -> None:
    """Apply a restrained publication style with editable vector text."""
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 6.6,
            "axes.titlesize": 7.2,
            "axes.titleweight": "semibold",
            "axes.labelsize": 6.8,
            "axes.linewidth": 0.65,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.labelsize": 5.7,
            "ytick.labelsize": 5.9,
            "xtick.major.size": 2.4,
            "ytick.major.size": 2.4,
            "xtick.major.width": 0.55,
            "ytick.major.width": 0.55,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "axes.unicode_minus": False,
            "savefig.facecolor": "white",
        }
    )


def resolve_project_root() -> Path:
    script_path = Path(__file__).resolve()
    for parent in script_path.parents:
        if parent.name == PROJECT_MARKER:
            return parent
    raise RuntimeError(f"Unable to locate {PROJECT_MARKER!r} above {script_path}")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    project_root = resolve_project_root()
    default_rca_dir = project_root / "list" / "RCA"
    parser = argparse.ArgumentParser(
        description="Create combined temperature methylation scatter/boxplot figures."
    )
    parser.add_argument("--rca-dir", type=Path, default=default_rca_dir)
    parser.add_argument(
        "--result-dir",
        type=Path,
        default=project_root / "list" / "xgboot" / IMPORT_RESULTS_DIRNAME,
    )
    parser.add_argument(
        "--output-dir", type=Path, default=default_rca_dir / OUTPUT_DIRNAME
    )
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument(
        "--qa-tools-dir",
        type=Path,
        default=None,
        help=(
            "Directory containing the figure QA scripts. If omitted, use "
            "NATURE_FIGURE_QA_DIR or the standard local skill installation."
        ),
    )
    args = parser.parse_args(argv)
    if args.top_n < 2:
        raise ValueError("--top-n must be at least 2")
    if not 0 < args.alpha < 1:
        raise ValueError("--alpha must lie strictly between 0 and 1")
    if args.dpi < 300:
        raise ValueError("--dpi must be at least 300")
    return args


def read_table(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Input table not found: {path}")
    data = pd.read_csv(path, sep="\t", low_memory=False)
    if data.empty:
        raise ValueError(f"Input table is empty: {path}")
    return data


def require_columns(data: pd.DataFrame, columns: Sequence[str], source: Path | str) -> None:
    missing = [column for column in columns if column not in data.columns]
    if missing:
        raise KeyError(f"{source} is missing required columns: {', '.join(missing)}")


def numeric_alignment_keys(data: pd.DataFrame, source: Path | str) -> pd.DataFrame:
    require_columns(data, ALIGNMENT_COLUMNS, source)
    converted = data.copy()
    for column in ALIGNMENT_COLUMNS:
        try:
            converted[column] = pd.to_numeric(converted[column], errors="raise")
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{source} has a non-numeric alignment key in {column}") from exc
    if converted.loc[:, ALIGNMENT_COLUMNS].isna().any().any():
        raise ValueError(f"{source} contains missing alignment-key values")
    return converted


def align_coverage_rows(
    methylation_raw: pd.DataFrame,
    coverage_raw: pd.DataFrame,
    methylation_path: Path,
    coverage_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    """Align coverage to methylation rows using the reference script's keys."""
    methylation = numeric_alignment_keys(methylation_raw, methylation_path)
    coverage = numeric_alignment_keys(coverage_raw, coverage_path)
    keys = list(ALIGNMENT_COLUMNS)

    methylation_duplicates = int(methylation.duplicated(keys, keep=False).sum())
    coverage_duplicates = int(coverage.duplicated(keys, keep=False).sum())
    if methylation_duplicates:
        raise ValueError(
            f"Methylation table has {methylation_duplicates} rows with duplicate alignment keys"
        )
    if coverage_duplicates:
        raise ValueError(
            f"Coverage table has {coverage_duplicates} rows with duplicate alignment keys"
        )

    left_keys = methylation.reset_index(drop=True).reset_index(names="methylation_index")
    right_keys = coverage.reset_index(drop=True).reset_index(names="coverage_index")
    try:
        matched = left_keys[["methylation_index", *keys]].merge(
            right_keys[["coverage_index", *keys]],
            on=keys,
            how="left",
            validate="one_to_one",
        )
    except pd.errors.MergeError as exc:
        raise ValueError(f"Coverage alignment is ambiguous: {exc}") from exc

    if len(matched) != len(methylation):
        raise ValueError("Coverage alignment changed the methylation row count")
    unmatched = int(matched["coverage_index"].isna().sum())
    if unmatched:
        raise ValueError(f"Coverage alignment failed for {unmatched} methylation rows")

    coverage_indices = (
        matched.sort_values("methylation_index")["coverage_index"].astype(int).to_numpy()
    )
    aligned_methylation = methylation.reset_index(drop=True).copy()
    aligned_coverage = coverage.iloc[coverage_indices].reset_index(drop=True).copy()
    aligned_methylation["_source_row"] = np.arange(2, len(aligned_methylation) + 2)

    if not np.allclose(
        aligned_methylation.loc[:, keys].to_numpy(dtype=float),
        aligned_coverage.loc[:, keys].to_numpy(dtype=float),
        rtol=0.0,
        atol=1e-10,
    ):
        raise ValueError("Coverage alignment keys do not match after reordering")

    report = {
        "methylation_rows": int(len(aligned_methylation)),
        "coverage_rows": int(len(coverage)),
        "matched_rows": int(len(coverage_indices)),
        "extra_coverage_rows": int(len(coverage) - len(np.unique(coverage_indices))),
    }
    return aligned_methylation, aligned_coverage, report


def validate_temperature_design(data: pd.DataFrame) -> dict[int, int]:
    require_columns(data, ["Temperature"], "methylation table")
    temperatures = pd.to_numeric(data["Temperature"], errors="coerce")
    if temperatures.isna().any() or not np.isfinite(temperatures).all():
        raise ValueError("Temperature contains missing or non-finite values")
    observed = tuple(sorted(int(value) for value in temperatures.unique()))
    if observed != TEMPERATURES:
        raise ValueError(f"Expected temperatures {TEMPERATURES}, observed {observed}")
    counts = temperatures.astype(int).value_counts().reindex(TEMPERATURES, fill_value=0)
    if (counts < 2).any():
        raise ValueError("Each temperature must contain at least two rows for sample SD")
    return {int(temp): int(counts.loc[temp]) for temp in TEMPERATURES}


def parse_shap_features(result_path: Path) -> dict[str, list[tuple[str, float]]]:
    if not result_path.is_file():
        raise FileNotFoundError(f"SHAP result file not found: {result_path}")
    grouped: dict[str, list[tuple[str, float]]] = {
        feature_group: [] for feature_group in FEATURE_GROUPS
    }
    seen: set[str] = set()
    found_header = False
    in_feature_block = False

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
                raise ValueError(
                    f"Malformed SHAP line at {result_path}:{line_number}: {line.strip()!r}"
                )
            feature_name, shap_text = match.groups()
            if feature_name in seen:
                raise ValueError(f"Duplicate SHAP feature {feature_name!r} in {result_path}")
            seen.add(feature_name)
            shap_value = float(shap_text)
            if not np.isfinite(shap_value):
                raise ValueError(f"Non-finite SHAP value for {feature_name}")
            feature_group = feature_name.split("_", 1)[0]
            if feature_group in grouped and shap_value > 0:
                grouped[feature_group].append((feature_name, shap_value))

    if not found_header:
        raise ValueError(f"SHAP feature block was not found in {result_path}")
    for feature_group in FEATURE_GROUPS:
        grouped[feature_group].sort(key=lambda item: item[1], reverse=True)
        if not grouped[feature_group]:
            raise ValueError(f"No positive-SHAP {feature_group} features in {result_path}")
    return grouped


def select_top_features(
    grouped: Mapping[str, Sequence[tuple[str, float]]], top_n: int
) -> dict[str, list[tuple[str, float]]]:
    selected: dict[str, list[tuple[str, float]]] = {}
    for feature_group in FEATURE_GROUPS:
        items = list(grouped[feature_group][:top_n])
        if len(items) != top_n:
            raise ValueError(
                f"{feature_group} contains {len(items)} positive-SHAP features; "
                f"{top_n} are required"
            )
        selected[feature_group] = items
    return selected


def validate_selected_features(
    selected: Mapping[str, Sequence[tuple[str, float]]],
    methylation: pd.DataFrame,
    coverage: pd.DataFrame,
) -> None:
    names = [name for group in FEATURE_GROUPS for name, _ in selected[group]]
    require_columns(methylation, names, "methylation table")
    require_columns(coverage, names, "coverage table")

    methylation_values = methylation.loc[:, names].apply(
        pd.to_numeric, errors="coerce"
    ).to_numpy(dtype=float)
    coverage_values = coverage.loc[:, names].apply(
        pd.to_numeric, errors="coerce"
    ).to_numpy(dtype=float)
    finite_methylation = methylation_values[np.isfinite(methylation_values)]
    finite_coverage = coverage_values[np.isfinite(coverage_values)]
    if finite_methylation.size == 0 or finite_coverage.size == 0:
        raise ValueError("Selected feature columns contain no usable numeric values")
    if np.any((finite_methylation < 0) | (finite_methylation > 100)):
        raise ValueError("Selected methylation values must lie within 0-100 percent")
    if np.any(finite_coverage < 0):
        raise ValueError("Selected coverage weights must be non-negative")


def subset_for_target(
    methylation: pd.DataFrame,
    coverage: pd.DataFrame,
    expression_column: str,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[int, int]]:
    require_columns(methylation, [expression_column], "methylation table")
    expression = pd.to_numeric(methylation[expression_column], errors="coerce")
    eligible = expression.notna() & np.isfinite(expression)
    if expression_column in ABSOLUTE_EXPRESSION_COLUMNS:
        eligible &= expression > 0
    if not eligible.any():
        raise ValueError(f"No eligible samples for expression column {expression_column}")

    selected_methylation = methylation.loc[eligible].reset_index(drop=True).copy()
    selected_coverage = coverage.loc[eligible].reset_index(drop=True).copy()
    counts = validate_temperature_design(selected_methylation)
    return selected_methylation, selected_coverage, counts


def finite_feature_values(data: pd.DataFrame, feature: str, temperature: int) -> np.ndarray:
    values = pd.to_numeric(
        data.loc[data["Temperature"].eq(temperature), feature], errors="coerce"
    ).to_numpy(dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 2:
        raise ValueError(
            f"{feature} at {temperature} degrees C has fewer than two finite observations"
        )
    return values


def build_source_data(
    dataset_name: str,
    selected: Mapping[str, Sequence[tuple[str, float]]],
    methylation: pd.DataFrame,
    coverage: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_rows: list[dict[str, object]] = []
    weighted_rows: list[dict[str, object]] = []

    for feature_group in FEATURE_GROUPS:
        feature_items = list(selected[feature_group])
        feature_names = [feature for feature, _ in feature_items]
        for rank, (feature, shap_value) in enumerate(feature_items, start=1):
            for temperature in TEMPERATURES:
                values = finite_feature_values(methylation, feature, temperature)
                summary_rows.append(
                    {
                        "dataset": dataset_name,
                        "feature_group": feature_group,
                        "shap_rank": rank,
                        "feature": feature,
                        "mean_abs_shap": shap_value,
                        "temperature_c": temperature,
                        "n": int(values.size),
                        "mean_methylation_percent": float(np.mean(values)),
                        "sample_sd_methylation_percent": float(np.std(values, ddof=1)),
                        "minimum_methylation_percent": float(np.min(values)),
                        "maximum_methylation_percent": float(np.max(values)),
                        "median_methylation_percent": float(np.median(values)),
                    }
                )

        methylation_values = methylation.loc[:, feature_names].apply(
            pd.to_numeric, errors="coerce"
        ).to_numpy(dtype=float)
        coverage_values = coverage.loc[:, feature_names].apply(
            pd.to_numeric, errors="coerce"
        ).to_numpy(dtype=float)
        valid = (
            np.isfinite(methylation_values)
            & np.isfinite(coverage_values)
            & (coverage_values > 0)
        )
        weight_sums = np.where(valid, coverage_values, 0.0).sum(axis=1)
        weighted_sums = np.where(
            valid, methylation_values * coverage_values, 0.0
        ).sum(axis=1)
        scores = np.divide(
            weighted_sums,
            weight_sums,
            out=np.full(len(methylation), np.nan, dtype=float),
            where=weight_sums > 0,
        )
        valid_feature_counts = valid.sum(axis=1)

        for row_index in range(len(methylation)):
            weighted_rows.append(
                {
                    "dataset": dataset_name,
                    "feature_group": feature_group,
                    "source_row": int(methylation.loc[row_index, "_source_row"]),
                    "temperature_c": int(methylation.loc[row_index, "Temperature"]),
                    "weighted_methylation_percent": float(scores[row_index]),
                    "valid_feature_count": int(valid_feature_counts[row_index]),
                    "coverage_weight_sum": float(weight_sums[row_index]),
                    "selected_feature_count": len(feature_names),
                }
            )

    summary = pd.DataFrame.from_records(summary_rows)
    weighted = pd.DataFrame.from_records(weighted_rows)
    return summary, weighted


def add_statistics_medians(statistics: pd.DataFrame, weighted: pd.DataFrame) -> pd.DataFrame:
    medians = weighted.groupby(["feature_group", "temperature_c"])[
        "weighted_methylation_percent"
    ].median()
    result = statistics.copy()
    result["median_methylation_percent"] = [
        "; ".join(
            f"{temperature}°C: {medians.loc[(row.feature_group, temperature)]:.15g}"
            for temperature in (
                TEMPERATURES if row.test_scope == "overall"
                else (int(row.temperature_1_c), int(row.temperature_2_c))
            )
        )
        for row in statistics.itertuples(index=False)
    ]
    return result


def bh_adjust(p_values: Sequence[float]) -> list[float]:
    pvals = np.asarray(p_values, dtype=float)
    if pvals.size == 0:
        return []
    if np.any(~np.isfinite(pvals)) or np.any((pvals < 0) | (pvals > 1)):
        raise ValueError("P values for BH correction must be finite and within 0-1")
    order = np.argsort(pvals)
    ranked = pvals[order]
    adjusted_ranked = np.empty(len(pvals), dtype=float)
    running = 1.0
    for index in range(len(pvals) - 1, -1, -1):
        running = min(running, ranked[index] * len(pvals) / (index + 1))
        adjusted_ranked[index] = running
    adjusted = np.empty(len(pvals), dtype=float)
    adjusted[order] = np.clip(adjusted_ranked, 0.0, 1.0)
    return adjusted.tolist()


def significance_label(p_value: float, alpha: float) -> str:
    if p_value >= alpha:
        return "ns"
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    return "*"


def dunn_pairwise_raw(groups: Mapping[int, np.ndarray]) -> list[dict[str, float | int]]:
    arrays = [np.asarray(groups[temp], dtype=float) for temp in TEMPERATURES]
    if any(array.size == 0 for array in arrays):
        raise ValueError("Dunn test requires a non-empty group for every temperature")
    all_values = np.concatenate(arrays)
    group_sizes = [len(array) for array in arrays]

    if np.all(all_values == all_values[0]):
        z_scores = [0.0] * len(PAIRWISE_COMPARISONS)
        p_values = [1.0] * len(PAIRWISE_COMPARISONS)
    else:
        ranks = stats.rankdata(all_values)
        offsets = np.cumsum([0, *group_sizes])
        mean_ranks = [
            float(np.mean(ranks[offsets[index] : offsets[index + 1]]))
            for index in range(len(arrays))
        ]
        n_total = len(all_values)
        tie_counts = np.unique(all_values, return_counts=True)[1]
        tie_term = (
            float(np.sum(tie_counts**3 - tie_counts) / (12 * (n_total - 1)))
            if n_total > 1
            else 0.0
        )
        rank_variance = n_total * (n_total + 1) / 12 - tie_term
        if rank_variance <= 0:
            z_scores = [0.0] * len(PAIRWISE_COMPARISONS)
            p_values = [1.0] * len(PAIRWISE_COMPARISONS)
        else:
            z_scores = []
            p_values = []
            for first, second in PAIRWISE_COMPARISONS:
                first_index = TEMPERATURES.index(first)
                second_index = TEMPERATURES.index(second)
                denominator = np.sqrt(
                    rank_variance
                    * (1 / group_sizes[first_index] + 1 / group_sizes[second_index])
                )
                z_score = abs(mean_ranks[first_index] - mean_ranks[second_index]) / denominator
                z_scores.append(float(z_score))
                p_values.append(float(2 * stats.norm.sf(z_score)))

    records: list[dict[str, float | int]] = []
    for (first, second), z_score, p_value in zip(
        PAIRWISE_COMPARISONS, z_scores, p_values
    ):
        records.append(
            {
                "temperature_1_c": first,
                "temperature_2_c": second,
                "n_1": int(len(groups[first])),
                "n_2": int(len(groups[second])),
                "statistic": z_score,
                "raw_p": p_value,
            }
        )
    return records


def run_statistics(
    dataset_name: str, weighted: pd.DataFrame, alpha: float
) -> pd.DataFrame:
    overall_records: list[dict[str, object]] = []
    pairwise_records: list[dict[str, object]] = []

    for feature_group in FEATURE_GROUPS:
        subset = weighted.loc[weighted["feature_group"].eq(feature_group)]
        groups = {
            temperature: subset.loc[
                subset["temperature_c"].eq(temperature),
                "weighted_methylation_percent",
            ]
            .to_numpy(dtype=float)
            for temperature in TEMPERATURES
        }
        groups = {
            temperature: values[np.isfinite(values)]
            for temperature, values in groups.items()
        }
        if any(values.size == 0 for values in groups.values()):
            raise ValueError(f"No weighted {feature_group} values for one temperature")
        all_values = np.concatenate(list(groups.values()))
        if np.all(all_values == all_values[0]):
            statistic, raw_p = 0.0, 1.0
        else:
            result = stats.kruskal(*(groups[temp] for temp in TEMPERATURES))
            statistic, raw_p = float(result.statistic), float(result.pvalue)
        overall_records.append(
            {
                "dataset": dataset_name,
                "feature_group": feature_group,
                "test_scope": "overall",
                "comparison": "10_vs_16_vs_22",
                "temperature_1_c": np.nan,
                "temperature_2_c": np.nan,
                "n_1": int(len(groups[10])),
                "n_2": int(len(groups[16])),
                "n_3": int(len(groups[22])),
                "test_method": "Kruskal-Wallis",
                "statistic": statistic,
                "raw_p": raw_p,
            }
        )
        for record in dunn_pairwise_raw(groups):
            pairwise_records.append(
                {
                    "dataset": dataset_name,
                    "feature_group": feature_group,
                    "test_scope": "pairwise",
                    "comparison": (
                        f"{record['temperature_1_c']}_vs_{record['temperature_2_c']}"
                    ),
                    "temperature_1_c": record["temperature_1_c"],
                    "temperature_2_c": record["temperature_2_c"],
                    "n_1": record["n_1"],
                    "n_2": record["n_2"],
                    "n_3": np.nan,
                    "test_method": "Dunn",
                    "statistic": record["statistic"],
                    "raw_p": record["raw_p"],
                }
            )

    overall_adjusted = bh_adjust([float(row["raw_p"]) for row in overall_records])
    for row, adjusted_p in zip(overall_records, overall_adjusted):
        row.update(
            {
                "adjusted_p": adjusted_p,
                "correction": "Benjamini-Hochberg across 3 context tests within dataset",
                "correction_family_size": 3,
                "alpha": alpha,
                "significance": significance_label(adjusted_p, alpha),
            }
        )

    pairwise_adjusted = bh_adjust([float(row["raw_p"]) for row in pairwise_records])
    for row, adjusted_p in zip(pairwise_records, pairwise_adjusted):
        row.update(
            {
                "adjusted_p": adjusted_p,
                "correction": "Benjamini-Hochberg across 9 pairwise tests within dataset",
                "correction_family_size": 9,
                "alpha": alpha,
                "significance": significance_label(adjusted_p, alpha),
            }
        )

    return add_statistics_medians(
        pd.DataFrame.from_records([*overall_records, *pairwise_records]), weighted
    )


def format_p_value(value: float) -> str:
    if value < 0.001:
        return "< 0.001"
    return f"= {value:.3f}"


def shared_scatter_limits(summary: pd.DataFrame) -> tuple[float, float]:
    centers = summary["mean_methylation_percent"].to_numpy(dtype=float)
    spread = summary["sample_sd_methylation_percent"].to_numpy(dtype=float)
    extents = np.concatenate([centers - spread, centers + spread])
    if not np.isfinite(extents).all():
        raise ValueError("Scatter summary contains non-finite uncertainty extents")
    lower = float(np.min(extents))
    upper = float(np.max(extents))
    span = upper - lower
    if span <= 0:
        raise ValueError("Cannot derive scatter limits from constant summaries")
    padding = max(2.0, 0.055 * span)
    step = 5.0 if span <= 60 else 10.0
    return (
        float(np.floor((lower - padding) / step) * step),
        float(np.ceil((upper + padding) / step) * step),
    )


def add_panel_label(axis: Axes, label: str) -> None:
    offset = ScaledTranslation(-9.0 / 72.0, 4.0 / 72.0, axis.figure.dpi_scale_trans)
    axis.text(
        0,
        1,
        f"({label})",
        transform=axis.transAxes + offset,
        fontsize=8.0,
        fontweight="bold",
        ha="right",
        va="bottom",
        clip_on=False,
    )


def draw_scatter_panel(
    axis: Axes,
    feature_group: str,
    selected: Mapping[str, Sequence[tuple[str, float]]],
    summary: pd.DataFrame,
    y_limits: tuple[float, float],
) -> None:
    feature_items = list(selected[feature_group])
    feature_names = [name for name, _ in feature_items]
    for temperature in TEMPERATURES:
        rows = (
            summary.loc[
                summary["feature_group"].eq(feature_group)
                & summary["temperature_c"].eq(temperature)
            ]
            .sort_values("shap_rank")
            .reset_index(drop=True)
        )
        if len(rows) != len(feature_items):
            raise RuntimeError(f"Incomplete scatter summary for {feature_group}, {temperature}")
        x_values = np.arange(len(feature_items), dtype=float) + TEMPERATURE_OFFSETS[temperature]
        means = rows["mean_methylation_percent"].to_numpy(dtype=float)
        deviations = rows["sample_sd_methylation_percent"].to_numpy(dtype=float)
        axis.plot(
            x_values,
            means,
            color=TEMPERATURE_COLORS[temperature],
            linewidth=0.9,
            alpha=0.88,
            zorder=2,
        )
        axis.errorbar(
            x_values,
            means,
            yerr=deviations,
            fmt=TEMPERATURE_MARKERS[temperature],
            markersize=3.6,
            color=TEMPERATURE_COLORS[temperature],
            markeredgecolor="white",
            markeredgewidth=0.45,
            elinewidth=0.75,
            capsize=1.8,
            capthick=0.75,
            zorder=3,
        )

    axis.set_xlim(-0.55, len(feature_items) - 0.45)
    axis.set_ylim(*y_limits)
    axis.yaxis.set_major_locator(MaxNLocator(nbins=6))
    axis.set_xticks(np.arange(len(feature_items)))
    axis.set_xticklabels(
        feature_names,
        rotation=82,
        rotation_mode="anchor",
        ha="right",
        va="top",
        fontsize=5.2,
    )
    axis.set_title(f"{feature_group} features ({len(feature_items)})", loc="left")
    axis.grid(False)
    axis.set_axisbelow(True)
    axis.tick_params(axis="both", color="#555555")


def draw_box_panel(
    axis: Axes,
    feature_group: str,
    weighted: pd.DataFrame,
    statistics: pd.DataFrame,
    dataset_index: int,
) -> None:
    subset = weighted.loc[weighted["feature_group"].eq(feature_group)]
    arrays: list[np.ndarray] = []

    for position, temperature in enumerate(TEMPERATURES):
        values = subset.loc[
            subset["temperature_c"].eq(temperature), "weighted_methylation_percent"
        ].to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        if values.size == 0:
            raise ValueError(f"No weighted {feature_group} values at {temperature} degrees C")
        arrays.append(values)

    boxes = axis.boxplot(
        arrays,
        positions=np.arange(len(TEMPERATURES), dtype=float),
        widths=0.48,
        patch_artist=True,
        showfliers=False,
        manage_ticks=False,
        medianprops={"color": "#1C1C1C", "linewidth": 0.75},
        whiskerprops={"color": "#555555", "linewidth": 0.75},
        capprops={"color": "#555555", "linewidth": 0.75},
        boxprops={"color": "#555555", "linewidth": 0.75},
        zorder=2,
    )
    for patch, temperature in zip(boxes["boxes"], TEMPERATURES):
        patch.set_facecolor(TEMPERATURE_COLORS[temperature])
        patch.set_alpha(1.0)

    whisker_values = np.concatenate([line.get_ydata() for line in boxes["whiskers"]])
    data_min = float(np.min(whisker_values))
    data_max = float(np.max(whisker_values))
    data_span = max(data_max - data_min, 1.0)
    lower = max(0.0, data_min - 0.06 * data_span)
    panel_height_pt = axis.get_position().height * axis.figure.get_figheight() * 72
    data_height_pt = panel_height_pt - 64.0
    if data_height_pt <= 20:
        raise ValueError("Boxplot panel is too short for significance annotations")
    units_per_pt = max(data_max - lower, 1.0) / data_height_pt
    upper = lower + panel_height_pt * units_per_pt
    axis.set_ylim(lower, upper)
    pairwise = statistics.loc[
        statistics["feature_group"].eq(feature_group)
        & statistics["test_scope"].eq("pairwise")
    ]
    if len(pairwise) != len(PAIRWISE_COMPARISONS):
        raise RuntimeError(f"Expected three pairwise rows for {feature_group}")
    bracket_order = ((16, 22), (10, 16), (10, 22))
    for level, (first, second) in enumerate(bracket_order):
        row = pairwise.loc[
            pairwise["temperature_1_c"].eq(first)
            & pairwise["temperature_2_c"].eq(second)
        ]
        if len(row) != 1:
            raise RuntimeError(f"Missing pairwise statistics for {first} vs {second}")
        label = str(row.iloc[0]["significance"])
        first_position = TEMPERATURES.index(first)
        second_position = TEMPERATURES.index(second)
        bracket_y = data_max + (8 + 14 * level) * units_per_pt
        cap_height = 2 * units_per_pt
        axis.plot(
            [first_position, first_position, second_position, second_position],
            [bracket_y - cap_height, bracket_y, bracket_y, bracket_y - cap_height],
            color="#555555", linewidth=0.65,
        )
        axis.annotate(
            label,
            xy=((first_position + second_position) / 2, bracket_y),
            xytext=(0, 3), textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=5.7,
            fontweight="bold" if label != "ns" else "normal",
            linespacing=1.05,
            color="#303030",
            zorder=6,
        )

    if upper <= lower:
        raise RuntimeError(f"Invalid boxplot limits for {feature_group}")
    axis.set_ylim(lower, upper)
    axis.yaxis.set_major_locator(MaxNLocator(nbins=5))
    axis.set_xlim(-0.55, len(TEMPERATURES) - 0.45)
    axis.set_xticks(np.arange(len(TEMPERATURES), dtype=float))
    axis.set_xticklabels(
        [str(temp) for temp in TEMPERATURES]
    )
    overall = statistics.loc[
        statistics["feature_group"].eq(feature_group)
        & statistics["test_scope"].eq("overall")
    ]
    if len(overall) != 1:
        raise RuntimeError(f"Expected one overall statistic for {feature_group}")
    axis.set_title(feature_group, loc="left")
    axis.legend(
        handles=[Patch(facecolor=TEMPERATURE_COLORS[temp], edgecolor="#555555",
                       linewidth=0.75, label=f"{temp} °C") for temp in TEMPERATURES],
        loc="upper right", ncol=3, fontsize=5.2, handlelength=1.2,
        columnspacing=0.7, handletextpad=0.4, borderaxespad=0.5,
    )
    axis.set_xlabel("Temperature (°C)")
    # Keep the annotation band free of horizontal strokes so that significance
    # text cannot be mistaken for, or visually crossed by, a grid line.
    axis.grid(False)
    axis.set_axisbelow(True)
    axis.tick_params(axis="both", color="#555555")


def build_figure(
    dataset_name: str,
    selected: Mapping[str, Sequence[tuple[str, float]]],
    summary: pd.DataFrame,
    weighted: pd.DataFrame,
    statistics: pd.DataFrame,
    dataset_index: int,
) -> tuple[Figure, list[Axes]]:
    figure = plt.figure(
        figsize=(FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN), constrained_layout=False
    )
    grid = figure.add_gridspec(2, 3, height_ratios=(1.08, 1.0))
    top_axes: list[Axes] = []
    for column in range(3):
        axis = figure.add_subplot(
            grid[0, column], sharey=top_axes[0] if top_axes else None
        )
        top_axes.append(axis)
    bottom_axes = [figure.add_subplot(grid[1, column]) for column in range(3)]
    axes = [*top_axes, *bottom_axes]
    figure.subplots_adjust(
        left=0.088,
        right=0.985,
        top=0.915,
        bottom=0.13,
        wspace=0.30,
        hspace=0.63,
    )

    scatter_limits = shared_scatter_limits(summary)
    scatter_height_pt = top_axes[0].get_position().height * FIGURE_HEIGHT_IN * 72
    scatter_limits = (
        scatter_limits[0],
        scatter_limits[0] + (scatter_limits[1] - scatter_limits[0])
        / (1 - 21 / scatter_height_pt),
    )
    for index, (axis, feature_group) in enumerate(zip(top_axes, FEATURE_GROUPS)):
        draw_scatter_panel(axis, feature_group, selected, summary, scatter_limits)
        if index == 0:
            axis.set_ylabel("Methylation level (%)\n(mean ± sample SD)")
        else:
            axis.tick_params(labelleft=False)

    for axis, feature_group in zip(bottom_axes, FEATURE_GROUPS):
        draw_box_panel(axis, feature_group, weighted, statistics, dataset_index)
        axis.set_ylabel("Weighted methylation (%)")

    for label, axis in zip("abcdef", axes):
        add_panel_label(axis, label)

    temperature_handles = [
        Line2D(
            [0],
            [0],
            color=TEMPERATURE_COLORS[temp],
            marker=TEMPERATURE_MARKERS[temp],
            markersize=4.2,
            markeredgecolor="white",
            markeredgewidth=0.4,
            linewidth=0.9,
            label=f"{temp} °C",
        )
        for temp in TEMPERATURES
    ]
    figure.suptitle(
        f"{dataset_name}: temperature-associated methylation of SHAP-ranked features",
        x=0.5,
        y=0.982,
        fontsize=9.2,
        fontweight="semibold",
    )
    for axis in top_axes:
        axis.legend(handles=temperature_handles, loc="upper right", ncol=3,
                    fontsize=5.2, handlelength=1.2, columnspacing=0.7,
                    handletextpad=0.4, borderaxespad=0.5)
    figure.text(
        0.5,
        0.035,
        (
            "Top: mean ± sample SD; connecting lines are descriptive only.\n"
            "Bottom: median, IQR and 1.5×IQR whiskers; sample points hidden.\n"
            + "; ".join(
                group + " n=" + "/".join(str(int(weighted.loc[
                    weighted["feature_group"].eq(group)
                    & weighted["temperature_c"].eq(temp),
                    "weighted_methylation_percent"].notna().sum())) for temp in TEMPERATURES)
                for group in FEATURE_GROUPS
            ) + " (10/16/22 °C). n = valid weighted samples."
        ),
        ha="center",
        va="center",
        fontsize=5.5,
        color="#4A4A4A",
    )
    return figure, axes


def qa_scripts_directory(argument: Path | None) -> Path:
    if argument is not None:
        directory = argument
    elif os.environ.get("NATURE_FIGURE_QA_DIR"):
        directory = Path(os.environ["NATURE_FIGURE_QA_DIR"])
    else:
        directory = Path.home() / ".codex" / "skills" / "nature-figure" / "scripts"
    if not directory.is_dir():
        raise RuntimeError(
            "Figure QA directory not found; provide --qa-tools-dir or set "
            "NATURE_FIGURE_QA_DIR"
        )
    return directory


def load_alignment_helper(qa_directory: Path):
    helper = qa_directory / "audit_panel_alignment.py"
    if not helper.is_file():
        raise RuntimeError(f"Panel-alignment helper not found in {qa_directory}")
    sys.path.insert(0, str(qa_directory))
    try:
        from audit_panel_alignment import require_matplotlib_panel_alignment
    except ImportError as exc:
        raise RuntimeError(f"Could not import panel-alignment helper: {exc}") from exc
    return require_matplotlib_panel_alignment


def subprocess_failure(prefix: str, completed: subprocess.CompletedProcess[str]) -> RuntimeError:
    details = "\n".join(
        part.strip()
        for part in (completed.stdout, completed.stderr)
        if part and part.strip()
    )
    return RuntimeError(prefix + (f":\n{details}" if details else ""))


def run_rendered_qa(pdf_path: Path, qa_directory: Path, stem: str) -> list[Path]:
    collision_path = pdf_path.with_name(f"{stem}.collision-audit.json")
    collision_script = qa_directory / "audit_figure_collisions.py"
    if not collision_script.is_file():
        raise RuntimeError(f"Collision-audit helper not found in {qa_directory}")
    collision = subprocess.run(
        [
            sys.executable,
            str(collision_script),
            str(pdf_path),
            "--json-out",
            str(collision_path),
            "--strict",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    if collision.returncode != 0:
        raise subprocess_failure("Strict rendered-PDF collision audit failed", collision)
    collision_report = json.loads(collision_path.read_text(encoding="utf-8"))
    if collision_report.get("verdict") != "PASS":
        raise RuntimeError(
            f"Collision audit did not pass: {collision_report.get('verdict', 'unknown')}"
        )

    text_path = pdf_path.with_name(f"{stem}.pdf-text-audit.json")
    text_script = qa_directory / "audit_pdf_text.py"
    if not text_script.is_file():
        raise RuntimeError(f"PDF-text audit helper not found in {qa_directory}")
    text_audit = subprocess.run(
        [
            sys.executable,
            str(text_script),
            str(pdf_path),
            "--min-pt",
            str(MIN_FONT_PT),
            "--json",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    if text_audit.returncode != 0:
        raise subprocess_failure("Rendered PDF text audit failed", text_audit)
    text_report = json.loads(text_audit.stdout)
    if not text_report.get("auditable") or text_report.get("below_minimum_count") != 0:
        raise RuntimeError("Rendered PDF contains text below the 5 pt floor")
    text_path.write_text(
        json.dumps(text_report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return [collision_path, text_path]


def validate_panel_label_anchors(figure: Figure, axes: Sequence[Axes]) -> dict:
    figure.canvas.draw()
    anchors = []
    for axis, label in zip(axes, "abcdef"):
        matches = [artist for artist in axis.texts if artist.get_text() == f"({label})"]
        if len(matches) != 1 or matches[0].get_fontweight() != "bold":
            raise ValueError(f"Expected one bold ({label}) panel label")
        artist = matches[0]
        display_position = artist.get_transform().transform(artist.get_position())
        anchors.append(figure.dpi_scale_trans.inverted().transform(display_position) * 72)
    if len(anchors) != 6:
        raise ValueError("Expected six panel-label anchors")
    coordinates = np.asarray(anchors)
    deviations = [np.ptp(coordinates[row, 1]) for row in (slice(0, 3), slice(3, 6))]
    deviations.extend(np.ptp(coordinates[[column, column + 3], 0]) for column in range(3))
    if max(deviations) > 1.5:
        raise ValueError("Parenthesized panel labels exceed the 1.5 pt alignment tolerance")
    return {"verdict": "PASS", "anchors_pt": coordinates.tolist(),
            "labels": [f"({label})" for label in "abcdef"], "tolerance_pt": 1.5}


def export_figure(
    figure: Figure,
    axes: Sequence[Axes],
    staging_dir: Path,
    stem: str,
    dpi: int,
    qa_directory: Path,
    require_matplotlib_panel_alignment,
) -> list[Path]:
    paths = {
        "pdf": staging_dir / f"{stem}.pdf",
        "svg": staging_dir / f"{stem}.svg",
        "png": staging_dir / f"{stem}.png",
        "tiff": staging_dir / f"{stem}.tiff",
        "alignment": staging_dir / f"{stem}.alignment.json",
    }
    label_report = validate_panel_label_anchors(figure, axes)
    alignment_report = require_matplotlib_panel_alignment(
        figure,
        axes=list(axes),
        panel_ids=list("abcdef"),
        json_out=paths["alignment"],
        tolerance_pt=1.5,
        gutter_tolerance_pt=1.5,
        require_panel_labels=False,
        strict=True,
    )
    alignment_report["parenthesized_panel_labels"] = label_report
    paths["alignment"].write_text(json.dumps(alignment_report, indent=2) + "\n", encoding="utf-8")
    figure.savefig(paths["pdf"], format="pdf", facecolor="white")
    figure.savefig(paths["svg"], format="svg", facecolor="white")
    figure.savefig(paths["png"], format="png", dpi=dpi, facecolor="white")
    figure.savefig(
        paths["tiff"],
        format="tiff",
        dpi=dpi,
        facecolor="white",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(figure)
    qa_paths = run_rendered_qa(paths["pdf"], qa_directory, stem)
    return [*paths.values(), *qa_paths]


def write_tsv(data: pd.DataFrame, path: Path) -> None:
    data.to_csv(path, sep="\t", index=False, na_rep="NA")


def rewrite_qa_paths(path: Path, final_pdf: Path) -> None:
    if not path.name.endswith((".collision-audit.json", ".pdf-text-audit.json")):
        return
    report = json.loads(path.read_text(encoding="utf-8"))
    if "pdf" in report:
        report["pdf"] = str(final_pdf)
    path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    configure_style()
    qa_directory = qa_scripts_directory(args.qa_tools_dir)
    require_matplotlib_panel_alignment = load_alignment_helper(qa_directory)

    methylation_path = args.rca_dir / METHYLATION_FILENAME
    coverage_path = args.rca_dir / COVERAGE_FILENAME
    methylation_raw = read_table(methylation_path)
    coverage_raw = read_table(coverage_path)
    methylation, coverage, alignment_report = align_coverage_rows(
        methylation_raw, coverage_raw, methylation_path, coverage_path
    )
    temperature_counts = validate_temperature_design(methylation)

    parsed_targets: list[tuple[dict[str, str], dict[str, list[tuple[str, float]]]]] = []
    for config in TARGET_CONFIGS:
        result_path = args.result_dir / config["result"]
        grouped = parse_shap_features(result_path)
        selected = select_top_features(grouped, args.top_n)
        validate_selected_features(selected, methylation, coverage)
        parsed_targets.append((config, selected))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    staged_to_final: list[tuple[Path, Path]] = []
    run_targets: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory(prefix=".tem_meth.", dir=args.output_dir) as temp_name:
        staging_root = Path(temp_name)
        for dataset_index, (config, selected) in enumerate(parsed_targets):
            dataset_name = config["name"]
            target_methylation, target_coverage, target_counts = subset_for_target(
                methylation, coverage, config["expression"]
            )
            summary, weighted = build_source_data(
                dataset_name, selected, target_methylation, target_coverage
            )
            statistics = run_statistics(dataset_name, weighted, args.alpha)
            figure, axes = build_figure(
                dataset_name,
                selected,
                summary,
                weighted,
                statistics,
                dataset_index,
            )
            stem = f"{dataset_name}_temperature_methylation_top{args.top_n}"
            generated = export_figure(
                figure,
                axes,
                staging_root,
                stem,
                args.dpi,
                qa_directory,
                require_matplotlib_panel_alignment,
            )
            table_paths = [
                staging_root / f"{stem}.scatter-summary.tsv",
                staging_root / f"{stem}.weighted-samples.tsv",
                staging_root / f"{stem}.statistics.tsv",
            ]
            write_tsv(summary, table_paths[0])
            write_tsv(weighted, table_paths[1])
            write_tsv(statistics, table_paths[2])
            generated.extend(table_paths)
            caption_path = staging_root / f"{stem}.caption.txt"
            caption_lines = [
                f"{dataset_name}: temperature-associated methylation.",
                "a-c: CHH, CHG and CG features ranked by mean absolute SHAP value; mean ± sample SD (ddof=1).",
                "Lines connect feature means for description only, not a fitted temperature response.",
                "d-f: coverage-weighted methylation across the same selected features; median, IQR and 1.5×IQR whiskers.",
                "All valid samples contribute to the boxes and tests; individual points and outlier symbols are not displayed.",
                "n denotes usable input rows, not an independently verified biological-replicate count. Missing values make n feature-specific.",
                "Kruskal-Wallis overall tests: BH correction across 3 contexts per target; Dunn two-sided comparisons: BH across 9 pairs per target.",
                f"Significance: *** q<0.001; ** q<0.01; * q<{args.alpha:g}; ns otherwise. Overall tests are retained in the statistics TSV.",
                "Mean ± SD may extend outside 0-100%; these error bars are not clipped. Box axes follow visible whiskers, not hidden outliers.",
                "Boxplot n at 10/16/22 °C:",
            ]
            for feature_group in FEATURE_GROUPS:
                counts = [int(weighted.loc[
                    weighted["feature_group"].eq(feature_group)
                    & weighted["temperature_c"].eq(temp),
                    "weighted_methylation_percent"].notna().sum()) for temp in TEMPERATURES]
                caption_lines.append(f"{feature_group}: " + "/".join(map(str, counts)))
            caption_lines.append("Scatter n per feature at 10/16/22 °C:")
            for feature_group in FEATURE_GROUPS:
                for feature, _ in selected[feature_group]:
                    counts = summary.loc[summary["feature"].eq(feature)].sort_values("temperature_c")["n"]
                    caption_lines.append(feature + ": " + "/".join(str(int(count)) for count in counts))
            caption_path.write_text("\n".join(caption_lines) + "\n", encoding="utf-8")
            generated.append(caption_path)

            final_pdf = args.output_dir / f"{stem}.pdf"
            for staged_path in generated:
                rewrite_qa_paths(staged_path, final_pdf)
                staged_to_final.append((staged_path, args.output_dir / staged_path.name))
            run_targets.append(
                {
                    "dataset": dataset_name,
                    "expression_eligibility_column": config["expression"],
                    "eligible_temperature_counts": target_counts,
                    "selected_features_per_context": args.top_n,
                    "overall_significant_contexts": int(
                        (
                            statistics.loc[statistics["test_scope"].eq("overall"), "adjusted_p"]
                            < args.alpha
                        ).sum()
                    ),
                    "significant_pairwise_comparisons": int(
                        (
                            statistics.loc[statistics["test_scope"].eq("pairwise"), "adjusted_p"]
                            < args.alpha
                        ).sum()
                    ),
                }
            )

        run_summary_path = staging_root / "temperature_methylation_run_summary.json"
        run_summary = {
            "methylation_input": str(methylation_path),
            "coverage_input": str(coverage_path),
            "import_results_directory": str(args.result_dir),
            "output_directory": str(args.output_dir),
            "temperature_counts_before_target_filter": temperature_counts,
            "coverage_alignment": alignment_report,
            "scatter_center": "mean",
            "scatter_spread": "sample SD (ddof=1)",
            "scatter_lines": "descriptive only",
            "boxplot_definition": "median, IQR, 1.5xIQR whiskers from all usable samples; sample and outlier points hidden",
            "weighted_summary": "coverage-weighted mean across the same SHAP Top-N features",
            "overall_test": "Kruskal-Wallis",
            "overall_correction": "BH across 3 context tests within each dataset",
            "pairwise_test": "Dunn",
            "pairwise_correction": "BH across 9 pairwise tests within each dataset",
            "alpha": args.alpha,
            "targets": run_targets,
        }
        run_summary_path.write_text(
            json.dumps(run_summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )
        staged_to_final.append(
            (run_summary_path, args.output_dir / run_summary_path.name)
        )

        for staged_path, final_path in staged_to_final:
            os.replace(staged_path, final_path)

    print(
        "Coverage alignment: "
        + ", ".join(f"{key}={value}" for key, value in alignment_report.items())
    )
    for target in run_targets:
        print(
            f"[{target['dataset']}] counts={target['eligible_temperature_counts']}; "
            f"overall significant={target['overall_significant_contexts']}/3; "
            f"pairwise significant={target['significant_pairwise_comparisons']}/9"
        )
    print(f"Published {len(staged_to_final)} files to {args.output_dir}")


if __name__ == "__main__":
    main()
