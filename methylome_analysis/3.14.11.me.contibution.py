#!/usr/bin/env python3
"""Generate SHAP-positive feature scatter plots and SHAP-weighted boxplots."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams
from matplotlib.patches import Patch
from scipy import stats


FEATURE_GROUPS = ("CHH", "CHG", "CG", "snp")
SCATTER_FEATURE_GROUPS = ("CHH", "CHG", "CG")
BOXPLOT_FEATURE_GROUPS = ("CHH", "CHG", "CG")
DISPLAY_GROUPS = {"CHH": "CHH", "CHG": "CHG", "CG": "CG", "snp": "SNP"}
SCATTER_TOP_N = 10
EXPRESSION_ORDER = ("High", "Medium", "Low")
EXPRESSION_LABELS = {"High": "High", "Medium": "Medium", "Low": "Low"}
PAIRWISE_COMPARISONS = (("High", "Medium"), ("High", "Low"), ("Medium", "Low"))
PAIRWISE_LEVELS = {"High_vs_Medium": 1, "Medium_vs_Low": 2, "High_vs_Low": 3}
GROUP_COLORS = {"High": "#d1495b", "Medium": "#2a9d8f", "Low": "#4c6ef5"}
BOX_OFFSETS = {"High": -0.24, "Medium": 0.0, "Low": 0.24}
TARGET_CONFIGS = (
    {"result_file": "total_TPM_import_results.txt", "expr_col": "total", "output_dir": "total_TPM", "title": "total_TPM"},
    {"result_file": "A_TPM_import_results.txt", "expr_col": "α", "output_dir": "A_TPM", "title": "A_TPM"},
    {"result_file": "B_TPM_import_results.txt", "expr_col": "β", "output_dir": "B_TPM", "title": "B_TPM"},
    {"result_file": "A_per_TPM_import_results.txt", "expr_col": "α/%", "output_dir": "A_per_TPM", "title": "A_per_TPM"},
    {"result_file": "B_per_TPM_import_results.txt", "expr_col": "β/%", "output_dir": "B_per_TPM", "title": "B_per_TPM"},
)
FEATURE_PATTERN = re.compile(r"^\s*((?:snp|CG|CHG|CHH)_[^:]+):\s*([-+0-9.eE]+)\s*$")
REPO_ROOT_MARKER = "Arabidopsis_thaliana"


def configure_fonts() -> None:
    rcParams["font.family"] = "DejaVu Sans"
    rcParams["axes.unicode_minus"] = False


def resolve_base_dir() -> Path:
    script_path = Path(__file__).resolve()
    for parent in script_path.parents:
        if parent.name == REPO_ROOT_MARKER:
            return parent
    raise RuntimeError(f"Unable to locate project root above {script_path}")


def candidate_paths(*paths: str) -> Iterable[Path]:
    for raw_path in paths:
        path = Path(raw_path)
        if path.exists():
            yield path


def resolve_io_paths() -> Tuple[Path, Path]:
    base_dir = resolve_base_dir()
    home_rca = base_dir / "list" / "RCA"
    home_results = base_dir / "list" / "xgboot" / "TPM_4.5_all"
    external_rca = Path("/datapool/life-gongl/zesheng/Arabidopsis_thaliana/list/RCA")
    external_results = Path("/datapool/life-gongl/zesheng/Arabidopsis_thaliana/list/xgboot/TPM_4.5_all")
    rca_dir = next(candidate_paths(str(external_rca), str(home_rca)), None)
    result_dir = next(candidate_paths(str(external_results), str(home_results)), None)
    if rca_dir is None or result_dir is None:
        raise FileNotFoundError("Could not resolve RCA or import_results directories.")
    return rca_dir, result_dir


def parse_shap_features(result_path: Path) -> Dict[str, List[Tuple[str, float]]]:
    grouped_features: Dict[str, List[Tuple[str, float]]] = {group: [] for group in FEATURE_GROUPS}
    in_feature_block = False

    with result_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if "所有特征及其平均|SHAP|值" in line:
                in_feature_block = True
                continue
            if not in_feature_block:
                continue
            match = FEATURE_PATTERN.match(line)
            if not match:
                if line.strip():
                    continue
                break
            feature_name, shap_text = match.groups()
            shap_value = float(shap_text)
            if shap_value <= 0:
                continue
            prefix = feature_name.split("_", 1)[0]
            grouped_features[prefix].append((feature_name, shap_value))

    for group_name in FEATURE_GROUPS:
        grouped_features[group_name] = sorted(grouped_features[group_name], key=lambda item: item[1], reverse=True)
    return grouped_features


def load_methylation_data(tsv_path: Path) -> pd.DataFrame:
    data_frame = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    return data_frame.apply(pd.to_numeric, errors="coerce")


def add_expression_groups(data_frame: pd.DataFrame, expr_col: str) -> Tuple[pd.DataFrame, float, float]:
    expr_series = pd.to_numeric(data_frame[expr_col], errors="coerce")
    eligible = expr_series.dropna()
    if eligible.empty:
        raise ValueError(f"No valid samples found for {expr_col}")
    q1 = eligible.quantile(0.20)
    q3 = eligible.quantile(0.80)
    grouped = data_frame.loc[eligible.index].copy()

    def classify(value: float) -> str:
        if value < q1:
            return "Low"
        if value > q3:
            return "High"
        return "Medium"

    grouped["expression_group"] = eligible.map(classify)
    return grouped, float(q1), float(q3)


def filter_scatter_features(grouped_features: Dict[str, List[Tuple[str, float]]]) -> Dict[str, List[Tuple[str, float]]]:
    return {group_name: grouped_features[group_name][:SCATTER_TOP_N] for group_name in SCATTER_FEATURE_GROUPS}


def filter_boxplot_features(grouped_features: Dict[str, List[Tuple[str, float]]]) -> Dict[str, List[Tuple[str, float]]]:
    return {group_name: list(grouped_features[group_name]) for group_name in BOXPLOT_FEATURE_GROUPS}


def bh_adjust(p_values: List[float]) -> List[float]:
    if not p_values:
        return []
    pvals = np.asarray(p_values, dtype=float)
    n_tests = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    adjusted = np.empty(n_tests, dtype=float)
    cumulative = 1.0
    for idx in range(n_tests - 1, -1, -1):
        rank = idx + 1
        cumulative = min(cumulative, ranked[idx] * n_tests / rank)
        adjusted[idx] = cumulative
    restored = np.empty(n_tests, dtype=float)
    restored[order] = np.clip(adjusted, 0.0, 1.0)
    return restored.tolist()


def significance_stars(p_value: float | None) -> str:
    if p_value is None or pd.isna(p_value):
        return ""
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return ""


def dunn_pairwise(feature_values: Dict[str, pd.Series]) -> List[Dict[str, object]]:
    arrays = [feature_values[group].dropna().to_numpy(dtype=float) for group in EXPRESSION_ORDER]
    if any(len(arr) == 0 for arr in arrays):
        return []

    all_values = np.concatenate(arrays)
    ranks = stats.rankdata(all_values)
    group_sizes = [len(arr) for arr in arrays]
    offsets = np.cumsum([0] + group_sizes)
    mean_ranks = [np.mean(ranks[offsets[idx] : offsets[idx + 1]]) for idx in range(len(arrays))]

    n_total = len(all_values)
    unique_counts = np.unique(all_values, return_counts=True)[1]
    tie_term = float(np.sum(unique_counts**3 - unique_counts) / (12 * (n_total - 1))) if n_total > 1 else 0.0
    variance = n_total * (n_total + 1) / 12 - tie_term
    if variance <= 0:
        return []

    records: List[Dict[str, object]] = []
    raw_pvals: List[float] = []
    for first, second in PAIRWISE_COMPARISONS:
        i = EXPRESSION_ORDER.index(first)
        j = EXPRESSION_ORDER.index(second)
        z_score = abs(mean_ranks[i] - mean_ranks[j]) / np.sqrt(variance * (1 / group_sizes[i] + 1 / group_sizes[j]))
        p_value = float(2 * stats.norm.sf(abs(z_score)))
        raw_pvals.append(p_value)
        records.append(
            {
                "comparison": f"{first}_vs_{second}",
                "group1": first,
                "group2": second,
                "test_method": "Dunn",
                "raw_p_value": p_value,
            }
        )

    for record, adjusted in zip(records, bh_adjust(raw_pvals)):
        record["adj_p_value"] = adjusted
        record["star"] = significance_stars(adjusted)
    return records


def compute_shap_weighted_scores(
    grouped_data: pd.DataFrame,
    grouped_features: Dict[str, List[Tuple[str, float]]],
) -> pd.DataFrame:
    weighted_df = grouped_data[["expression_group"]].copy()
    for group_name in BOXPLOT_FEATURE_GROUPS:
        feature_items = [(feature, shap) for feature, shap in grouped_features[group_name] if feature in grouped_data.columns]
        if not feature_items:
            weighted_df[group_name] = np.nan
            continue

        feature_names = [feature for feature, _ in feature_items]
        weights = np.asarray([shap for _, shap in feature_items], dtype=float)
        values = grouped_data[feature_names].to_numpy(dtype=float)
        valid_mask = ~np.isnan(values)
        normalized_weights = np.broadcast_to(weights, values.shape)
        weighted_values = np.full(values.shape[0], np.nan, dtype=float)
        weight_sums = np.where(valid_mask, normalized_weights, 0.0).sum(axis=1)
        valid_rows = weight_sums > 0
        if np.any(valid_rows):
            numerators = np.where(valid_mask[valid_rows], values[valid_rows] * normalized_weights[valid_rows], 0.0).sum(axis=1)
            weighted_values[valid_rows] = numerators / weight_sums[valid_rows]
        weighted_df[group_name] = weighted_values
    return weighted_df


def summarize_boxplot_group(feature_name: str, feature_values: Dict[str, pd.Series]) -> Dict[str, object]:
    row: Dict[str, object] = {"feature_group": feature_name, "feature_name": feature_name, "overall_test": "Kruskal_Wallis"}
    for group in EXPRESSION_ORDER:
        values = feature_values[group].dropna()
        row[f"n_{group.lower()}"] = int(len(values))
        row[f"mean_{group.lower()}"] = float(values.mean()) if not values.empty else np.nan
        row[f"median_{group.lower()}"] = float(values.median()) if not values.empty else np.nan
    return row


def run_group_level_statistics(weighted_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    feature_rows: List[Dict[str, object]] = []
    pairwise_rows: List[Dict[str, object]] = []
    overall_pvals: List[float] = []
    tested_groups: List[str] = []

    for group_name in BOXPLOT_FEATURE_GROUPS:
        feature_values = {
            expr_group: weighted_df.loc[weighted_df["expression_group"] == expr_group, group_name].dropna()
            for expr_group in EXPRESSION_ORDER
        }
        row = summarize_boxplot_group(group_name, feature_values)
        arrays = [feature_values[group].to_numpy(dtype=float) for group in EXPRESSION_ORDER]
        if any(len(arr) == 0 for arr in arrays):
            row["overall_p_value"] = np.nan
            row["overall_adj_p_value"] = np.nan
            row["overall_significant"] = False
            feature_rows.append(row)
            continue

        overall = stats.kruskal(*arrays)
        row["overall_p_value"] = float(overall.pvalue)
        row["overall_adj_p_value"] = np.nan
        row["overall_significant"] = False
        feature_rows.append(row)
        overall_pvals.append(float(overall.pvalue))
        tested_groups.append(group_name)

        for pair_row in dunn_pairwise(feature_values):
            pair_row["feature_group"] = group_name
            pair_row["feature_name"] = group_name
            pairwise_rows.append(pair_row)

    feature_df = pd.DataFrame(feature_rows)
    adjusted = bh_adjust(overall_pvals)
    for group_name, adj_p in zip(tested_groups, adjusted):
        mask = feature_df["feature_name"] == group_name
        feature_df.loc[mask, "overall_adj_p_value"] = adj_p
        feature_df.loc[mask, "overall_significant"] = adj_p < 0.05

    pairwise_df = pd.DataFrame(pairwise_rows)
    if not pairwise_df.empty:
        significant_groups = set(feature_df.loc[feature_df["overall_significant"], "feature_name"])
        pairwise_df = pairwise_df[pairwise_df["feature_name"].isin(significant_groups)].copy()

    return feature_df, pairwise_df


def add_significance_bracket(axis: plt.Axes, x1: float, x2: float, y: float, height: float, label: str) -> None:
    axis.plot([x1, x1, x2, x2], [y, y + height, y + height, y], color="#333333", linewidth=1.0)
    axis.text((x1 + x2) / 2, y + height, label, ha="center", va="bottom", fontsize=15, fontweight="bold")


def plot_scatter_figure(grouped_features: Dict[str, List[Tuple[str, float]]], grouped_data: pd.DataFrame, output_dir: Path) -> None:
    scatter_features = filter_scatter_features(grouped_features)
    max_features = max((len(items) for items in scatter_features.values()), default=1)
    fig_width = max(14.0, min(24.0, 6.0 + max_features * 1.35))
    fig_height = max(15.0, len(SCATTER_FEATURE_GROUPS) * 5.6)
    fig, axes = plt.subplots(len(SCATTER_FEATURE_GROUPS), 1, figsize=(fig_width, fig_height), constrained_layout=False)
    for axis, group_name in zip(axes, SCATTER_FEATURE_GROUPS):
        feature_items = scatter_features[group_name]
        if not feature_items:
            axis.axis("off")
            continue

        band_height = 100.0
        band_gap = 15.0
        band_bases = {"Low": 0.0, "Medium": band_height + band_gap, "High": 2 * (band_height + band_gap)}
        median_lines: Dict[str, List[float]] = {expr_group: [] for expr_group in EXPRESSION_ORDER}
        for expr_group in EXPRESSION_ORDER:
            subset = grouped_data[grouped_data["expression_group"] == expr_group]
            x_values: List[float] = []
            y_values: List[float] = []
            for idx, (feature_name, _shap) in enumerate(feature_items, start=1):
                values = pd.to_numeric(subset[feature_name], errors="coerce").dropna().to_numpy(dtype=float)
                group_base = band_bases[expr_group]
                jitter = np.linspace(-0.12, 0.12, num=len(values)) if len(values) > 1 else np.zeros(len(values))
                x_values.extend((idx + jitter).tolist())
                y_values.extend((group_base + values).tolist())
                median_lines[expr_group].append(float(np.nanmedian(values)) if len(values) else np.nan)
            axis.scatter(x_values, y_values, s=18, alpha=0.45, color=GROUP_COLORS[expr_group])
            valid_points = [
                (idx, group_base + value)
                for idx, value in enumerate(median_lines[expr_group], start=1)
                if not np.isnan(value)
            ]
            if valid_points:
                median_x, median_y = zip(*valid_points)
                axis.plot(median_x, median_y, color=GROUP_COLORS[expr_group], linewidth=2.4, alpha=0.95)

        tick_positions = list(range(1, len(feature_items) + 1))
        tick_labels = [feature for feature, _ in feature_items]
        axis.set_xticks(tick_positions)
        axis.set_xticklabels(tick_labels, rotation=0, ha="center", fontsize=15)
        axis.set_title(f"{DISPLAY_GROUPS[group_name]} features ({len(feature_items)})", fontsize=24, fontweight="bold", loc="left")
        axis.set_ylabel("Methylation level", fontsize=24, labelpad=4)
        axis.tick_params(axis="y", labelsize=20)
        y_ticks = [0.0, 50.0, 100.0, band_bases["Medium"], band_bases["Medium"] + 50.0, band_bases["Medium"] + 100.0, band_bases["High"], band_bases["High"] + 50.0, band_bases["High"] + 100.0]
        axis.set_yticks(y_ticks)
        axis.set_yticklabels(["0", "50", "100", "0", "50", "100", "0", "50", "100"], fontsize=20)
        axis.axhline(band_bases["Medium"], color="#bbbbbb", linestyle="--", linewidth=0.8, alpha=0.6)
        axis.axhline(band_bases["High"], color="#bbbbbb", linestyle="--", linewidth=0.8, alpha=0.6)
        axis.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.35)
        axis.set_ylim(-2.0, band_bases["High"] + band_height + 2.0)
        axis.set_xlim(0.65, len(feature_items) + 0.35)
        axis.text(1.008, 0.84, "High", transform=axis.transAxes, color=GROUP_COLORS["High"], fontsize=20, va="center")
        axis.text(1.008, 0.50, "Medium", transform=axis.transAxes, color=GROUP_COLORS["Medium"], fontsize=20, va="center")
        axis.text(1.008, 0.16, "Low", transform=axis.transAxes, color=GROUP_COLORS["Low"], fontsize=20, va="center")
    fig.subplots_adjust(left=0.04, right=0.985, top=0.985, bottom=0.06, hspace=0.30)
    fig.savefig(output_dir / "me_contribution_scatter.png", dpi=300, bbox_inches="tight")
    fig.savefig(output_dir / "me_contribution_scatter.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_single_boxplot(weighted_df: pd.DataFrame, pairwise_df: pd.DataFrame, output_dir: Path, feature_name: str) -> None:
    fig, axis = plt.subplots(figsize=(9, 6), constrained_layout=False)
    box_width = 0.12
    feature_max = 0.0

    for expr_group in EXPRESSION_ORDER:
        values = weighted_df.loc[weighted_df["expression_group"] == expr_group, feature_name].dropna()
        if values.empty:
            continue
        feature_max = max(feature_max, float(values.max()))
        axis.boxplot(
            values.to_numpy(dtype=float),
            positions=[1 + BOX_OFFSETS[expr_group]],
            widths=box_width,
            patch_artist=True,
            manage_ticks=False,
            boxprops={"facecolor": GROUP_COLORS[expr_group], "edgecolor": GROUP_COLORS[expr_group], "alpha": 0.35, "linewidth": 1.2},
            medianprops={"color": "#111111", "linewidth": 1.3},
            whiskerprops={"color": GROUP_COLORS[expr_group], "linewidth": 1.0},
            capprops={"color": GROUP_COLORS[expr_group], "linewidth": 1.0},
            flierprops={"marker": "o", "markersize": 2.0, "markerfacecolor": GROUP_COLORS[expr_group], "markeredgecolor": GROUP_COLORS[expr_group], "alpha": 0.35},
        )

    y_top = max(feature_max * 1.20, 5.0)
    bracket_base = max(y_top * 0.06, 2.0)
    bracket_height = max(y_top * 0.03, 1.0)
    subset = pairwise_df[pairwise_df["feature_name"] == feature_name].copy()
    significant_rows = subset[subset["star"] != ""].copy()
    significant_rows["level"] = significant_rows["comparison"].map(PAIRWISE_LEVELS).fillna(1).astype(int)
    significant_rows = significant_rows.sort_values("level")
    for row in significant_rows.itertuples(index=False):
        x1 = 1 + BOX_OFFSETS[row.group1]
        x2 = 1 + BOX_OFFSETS[row.group2]
        y = min(feature_max, y_top) + bracket_base * row.level
        add_significance_bracket(axis, x1, x2, y, bracket_height, row.star)

    axis.set_xlim(0.55, 1.45)
    axis.set_ylim(0, y_top + bracket_base * max(1, len(significant_rows) + 1))
    axis.set_xticks([1])
    axis.set_xticklabels([feature_name], fontsize=18)
    axis.set_ylabel("Methylation level", fontsize=22)
    axis.tick_params(axis="y", labelsize=16)
    axis.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.35)
    axis.set_title(feature_name, fontsize=22, fontweight="bold", loc="left")

    legend_handles = [
        Patch(facecolor=GROUP_COLORS[group], edgecolor=GROUP_COLORS[group], alpha=0.35, label=EXPRESSION_LABELS[group])
        for group in EXPRESSION_ORDER
    ]
    axis.legend(handles=legend_handles, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=16)
    fig.subplots_adjust(left=0.12, right=0.82, top=0.92, bottom=0.15)
    fig.savefig(output_dir / f"me_contribution_weighted_boxplot_{feature_name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(output_dir / f"me_contribution_weighted_boxplot_{feature_name}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_boxplot_figures(weighted_df: pd.DataFrame, pairwise_df: pd.DataFrame, output_dir: Path) -> None:
    for feature_name in BOXPLOT_FEATURE_GROUPS:
        plot_single_boxplot(weighted_df, pairwise_df, output_dir, feature_name)
    for suffix in ("png", "pdf"):
        for legacy_name in ("me_contribution_weighted_boxplot", "me_contribution_weighted_boxplot_snp"):
            legacy_path = output_dir / f"{legacy_name}.{suffix}"
            if legacy_path.exists():
                legacy_path.unlink()


def write_weighted_stats_outputs(output_dir: Path, prefix: str, feature_stats: pd.DataFrame, pairwise_stats: pd.DataFrame) -> None:
    feature_stats.to_csv(output_dir / f"{prefix}_weighted_group_stats.csv", index=False)
    pairwise_stats.to_csv(output_dir / f"{prefix}_weighted_group_pairwise_stats.csv", index=False)

    lines = [
        f"=== {prefix} weighted group statistics summary ===",
        "",
        f"tested_groups\t{len(feature_stats)}",
        f"overall_significant_groups\t{int(feature_stats['overall_significant'].eq(True).sum())}",
        f"pairwise_significant_comparisons\t{int((pairwise_stats['star'] != '').sum()) if not pairwise_stats.empty else 0}",
        "",
        "feature_group\toverall_adj_p_value\toverall_significant",
    ]
    for row in feature_stats.itertuples(index=False):
        lines.append(f"{row.feature_name}\t{row.overall_adj_p_value}\t{row.overall_significant}")

    (output_dir / f"{prefix}_weighted_group_stats_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def process_target(
    grouped_features: Dict[str, List[Tuple[str, float]]],
    methylation_df: pd.DataFrame,
    config: Dict[str, str],
    output_dir: Path,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    grouped_df, _q1, _q3 = add_expression_groups(methylation_df, config["expr_col"])
    boxplot_features = filter_boxplot_features(grouped_features)
    weighted_df = compute_shap_weighted_scores(grouped_df, boxplot_features)
    weighted_stats, weighted_pairwise = run_group_level_statistics(weighted_df)
    write_weighted_stats_outputs(output_dir, config["title"], weighted_stats, weighted_pairwise)
    plot_scatter_figure(grouped_features, grouped_df, output_dir)
    plot_boxplot_figures(weighted_df, weighted_pairwise, output_dir)


def validate_features(grouped_features: Dict[str, List[Tuple[str, float]]], result_file: str) -> None:
    missing_groups = [group for group in BOXPLOT_FEATURE_GROUPS if not grouped_features[group]]
    if missing_groups:
        raise ValueError(f"{result_file} has no positive-SHAP features for groups: {', '.join(missing_groups)}")


def main() -> None:
    configure_fonts()
    rca_dir, result_dir = resolve_io_paths()
    tsv_path = rca_dir / "Alt.all.snp_meth.filtered.tsv"
    if not tsv_path.exists():
        raise FileNotFoundError(f"Missing methylation matrix: {tsv_path}")

    methylation_df = load_methylation_data(tsv_path)
    for config in TARGET_CONFIGS:
        result_path = result_dir / config["result_file"]
        if not result_path.exists():
            raise FileNotFoundError(f"Missing import_results file: {result_path}")
        grouped_features = parse_shap_features(result_path)
        validate_features(grouped_features, config["result_file"])
        process_target(grouped_features, methylation_df, config, rca_dir / config["output_dir"])
        print(f"Saved scatter plots and SHAP-weighted boxplots for {config['title']} -> {rca_dir / config['output_dir']}")


if __name__ == "__main__":
    main()
