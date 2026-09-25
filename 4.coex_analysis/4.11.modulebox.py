#!/usr/bin/env python3
"""Recreate module expression boxplots without points.

Run from list/RCA after ``conda activate XGBoot`` using this script's full path.
Default outputs: list/RCA/modulebox. Override with --output-dir.
The observation unit is a gene's mean TPM within a temperature group.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cbook import boxplot_stats
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy.stats import kruskal
import scikit_posthocs as sp

BASE_DIR = Path(__file__).resolve().parents[2]
OUTPUT_DIR = BASE_DIR / "list/RCA/modulebox"
MODULE_FILES = {
    name: BASE_DIR / f"list/coex2/total/{name}_module_genes.txt"
    for name in ("yellow", "negative_grey")
}
EXPRESSION_FILES = {
    temp: BASE_DIR / f"list/coex2/RCA.{temp}.tran.txt"
    for temp in ("10C", "16C", "22C")
}
TEMPERATURE_ORDER = ["10C", "16C", "22C"]
PAIR_ORDER: List[Tuple[str, str]] = [("10C", "16C"), ("10C", "22C"), ("16C", "22C")]
PLOT_COLORS = {"10C": "#4C78A8", "16C": "#E69F55", "22C": "#409B8A"}
EXPECTED_MODULE_SIZES = {"yellow": 902, "negative_grey": 558}
EXPECTED_SAMPLE_SIZES = {"10C": 160, "16C": 163, "22C": 728}


def check_gene_ids(ids: pd.Series, context: str) -> None:
    if ids.isna().any() or ids.astype(str).str.strip().eq("").any():
        raise ValueError(f"{context}: missing gene ID")
    if ids.duplicated().any():
        raise ValueError(f"{context}: duplicate gene ID")


def read_module_genes(module_file: Path) -> List[str]:
    df = pd.read_csv(module_file, sep="\t", usecols=[0], header=0, dtype=str)
    check_gene_ids(df.iloc[:, 0], str(module_file))
    return df.iloc[:, 0].tolist()


def read_expression_matrix(expression_file: Path) -> pd.DataFrame:
    with expression_file.open() as handle:
        header = next(csv.reader(handle, delimiter="\t"))
    if len(header) < 2 or len(header) != len(set(header)) or any(not h.strip() for h in header):
        raise ValueError(f"{expression_file}: missing or duplicate sample column")
    df = pd.read_csv(expression_file, sep="\t", dtype={header[0]: str})
    return df.rename(columns={df.columns[0]: "gene_id"})


def validate_inputs(module_genes, expression_tables) -> None:
    for name, expected in EXPECTED_MODULE_SIZES.items():
        genes = pd.Series(module_genes[name], dtype="object")
        check_gene_ids(genes, name)
        if len(genes) != expected:
            raise ValueError(f"{name}: unexpected gene count {len(genes)} != {expected}")
    for temp, expected in EXPECTED_SAMPLE_SIZES.items():
        df = expression_tables[temp]
        if df.columns.duplicated().any() or "gene_id" not in df:
            raise ValueError(f"{temp}: invalid expression columns")
        check_gene_ids(df["gene_id"], temp)
        if df.shape[1] - 1 != expected:
            raise ValueError(f"{temp}: unexpected sample count {df.shape[1] - 1} != {expected}")
        values = df.drop(columns="gene_id")
        if not all(pd.api.types.is_numeric_dtype(dtype) for dtype in values.dtypes):
            raise ValueError(f"{temp}: nonnumeric expression value")
        array = values.to_numpy(dtype=float)
        if not np.isfinite(array).all():
            raise ValueError(f"{temp}: missing or nonfinite expression value")
        if (array < 0).any():
            raise ValueError(f"{temp}: negative expression value")
        available = set(df["gene_id"])
        for name, genes in module_genes.items():
            missing = set(genes) - available
            if missing:
                raise ValueError(f"{temp}/{name}: {len(missing)} missing module genes")


def build_module_long_table(module_name, module_genes, expression_tables):
    records = []
    for temp in TEMPERATURE_ORDER:
        subset = expression_tables[temp].loc[
            expression_tables[temp]["gene_id"].isin(module_genes)
        ].copy()
        means = subset.drop(columns="gene_id").mean(axis=1)
        if not np.isfinite(means).all() or (means < 0).any():
            raise ValueError(f"{module_name}/{temp}: invalid mean expression")
        records.append(pd.DataFrame({
            "module": module_name, "temperature": temp,
            "gene_id": subset["gene_id"], "mean_expression": means,
            "log1p_expression": np.log1p(means),
        }))
    result = pd.concat(records, ignore_index=True)
    result["temperature"] = pd.Categorical(
        result["temperature"], categories=TEMPERATURE_ORDER, ordered=True
    )
    return result


def p_to_stars(p_value: float) -> str:
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return "ns"


def run_statistics(long_df: pd.DataFrame, module_name: str) -> pd.DataFrame:
    groups = [
        long_df.loc[long_df["temperature"] == temperature, "log1p_expression"].to_numpy()
        for temperature in TEMPERATURE_ORDER
    ]
    if all(np.all(group == groups[0][0]) for group in groups):
        kw_stat, kw_p = 0.0, 1.0
    else:
        kw_stat, kw_p = kruskal(*groups)

    stat_rows = []

    for temperature in TEMPERATURE_ORDER:
        temp_df = long_df.loc[long_df["temperature"] == temperature]
        stat_rows.append(
            {
                "module": module_name,
                "record_type": "summary",
                "temperature": temperature,
                "test": "",
                "group1": "",
                "group2": "",
                "group3": "",
                "statistic": np.nan,
                "p_value": np.nan,
                "p_adjusted": np.nan,
                "significance": "",
                "mean_expression_mean": temp_df["mean_expression"].mean(),
                "mean_expression_median": temp_df["mean_expression"].median(),
                "log1p_expression_mean": temp_df["log1p_expression"].mean(),
                "log1p_expression_median": temp_df["log1p_expression"].median(),
            }
        )

    stat_rows.append(
        {
            "module": module_name,
            "record_type": "test",
            "temperature": "",
            "test": "Kruskal-Wallis",
            "group1": "10C",
            "group2": "16C",
            "group3": "22C",
            "statistic": kw_stat,
            "p_value": kw_p,
            "p_adjusted": kw_p,
            "significance": p_to_stars(kw_p) if kw_p < 0.05 else "ns",
            "mean_expression_mean": np.nan,
            "mean_expression_median": np.nan,
            "log1p_expression_mean": np.nan,
            "log1p_expression_median": np.nan,
        }
    )

    if kw_p < 0.05:
        dunn_df = sp.posthoc_dunn(
            long_df,
            val_col="log1p_expression",
            group_col="temperature",
            p_adjust="fdr_bh",
        )
        for group1, group2 in PAIR_ORDER:
            p_adj = float(dunn_df.loc[group1, group2])
            stat_rows.append(
                {
                    "module": module_name,
                    "record_type": "test",
                    "temperature": "",
                    "test": "Dunn_BH",
                    "group1": group1,
                    "group2": group2,
                    "group3": "",
                    "statistic": np.nan,
                    "p_value": np.nan,
                    "p_adjusted": p_adj,
                    "significance": p_to_stars(p_adj) if p_adj < 0.05 else "ns",
                    "mean_expression_mean": np.nan,
                    "mean_expression_median": np.nan,
                    "log1p_expression_mean": np.nan,
                    "log1p_expression_median": np.nan,
                }
            )

    return pd.DataFrame(stat_rows)


def plot_module(long_df, stats_df, module_name):
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": 7, "axes.labelsize": 8, "axes.titlesize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7,
        "pdf.fonttype": 42, "svg.fonttype": "none",
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.7, "savefig.facecolor": "white",
    })
    fig, ax = plt.subplots(figsize=(89 / 25.4, 108 / 25.4))
    fig.subplots_adjust(left=0.19, right=0.96, bottom=0.27, top=0.80)
    values = [long_df.loc[long_df["temperature"] == t, "log1p_expression"].to_numpy()
              for t in TEMPERATURE_ORDER]
    boxes = boxplot_stats(values, whis=1.5)
    artists = ax.bxp(
        boxes, positions=[0, 1, 2], widths=0.52, showfliers=False,
        patch_artist=True, manage_ticks=False,
        boxprops={"linewidth": 0.8, "edgecolor": "#333333"},
        medianprops={"linewidth": 1.2, "color": "#222222"},
        whiskerprops={"linewidth": 0.8, "color": "#333333"},
        capprops={"linewidth": 0.8, "color": "#333333"},
    )
    for patch, temp in zip(artists["boxes"], TEMPERATURE_ORDER):
        patch.set_facecolor(PLOT_COLORS[temp])
    ax.set_xticks([0, 1, 2], ["10 °C", "16 °C", "22 °C"])
    ax.set_xlim(-0.65, 2.65)
    ax.set_xlabel("Temperature", labelpad=5)
    ax.set_ylabel("ln(mean TPM + 1)", labelpad=5)
    title = {"yellow": "Yellow module", "negative_grey": "Negative grey module"}[module_name]
    fig.text(0.5, 0.95, title, ha="center", va="top", fontsize=9)
    fig.legend(
        handles=[Patch(facecolor=PLOT_COLORS[t], edgecolor="#333333",
                       linewidth=0.6, label=f"{t[:-1]} °C") for t in TEMPERATURE_ORDER],
        loc="upper right", bbox_to_anchor=(0.96, 0.89), ncol=3,
        frameon=False, fontsize=6, handlelength=1.1, columnspacing=0.8,
    )
    fig.text(
        0.5, 0.095,
        f"n = {long_df['gene_id'].nunique()} genes per temperature.\n"
        "Boxes: median and IQR; whiskers: 1.5 × IQR.\n"
        "ns: BH-adjusted P ≥ 0.05 (pairwise comparisons).",
        ha="center", va="top", fontsize=6, linespacing=1.5,
    )
    ax.tick_params(direction="out", width=0.7, length=3)
    lower = min(b["whislo"] for b in boxes)
    upper = max(b["whishi"] for b in boxes)
    span = max(upper - lower, 0.5)
    pairwise = stats_df.loc[stats_df["test"] == "Dunn_BH"]
    comparisons = sorted(pairwise.to_dict("records"), key=lambda r: (
        TEMPERATURE_ORDER.index(r["group2"]) - TEMPERATURE_ORDER.index(r["group1"]),
        TEMPERATURE_ORDER.index(r["group1"]),
    ))
    top = upper + 0.10 * span
    for row in comparisons:
        x1, x2 = [TEMPERATURE_ORDER.index(row[k]) for k in ("group1", "group2")]
        height = 0.035 * span
        ax.plot([x1, x1, x2, x2], [top, top + height, top + height, top],
                color="#222222", linewidth=0.7)
        ax.text((x1 + x2) / 2, top + height + 0.025 * span,
                row["significance"], ha="center", va="bottom", fontsize=8)
        top += 0.18 * span
    if not comparisons:
        overall = stats_df.loc[stats_df["test"] == "Kruskal-Wallis"]
        if not overall.empty and float(overall.iloc[0]["p_value"]) >= 0.05:
            ax.text(1, top, "Overall: ns", ha="center", va="bottom", fontsize=8)
            top += 0.18 * span
    # Keep zero-valued lower caps separate from the bottom axis spine.
    ax.set_ylim(lower - 0.05 * span, top + 0.08 * span)
    fig.canvas.draw()
    if ax.collections or artists["fliers"]:
        raise AssertionError("Unexpected scatter or outlier markers")
    return fig, {
        "module": module_name, "gene_count": long_df["gene_id"].nunique(),
        "record_count": len(long_df), "figure_size_mm": [89, 108],
        "scatter_count": len(ax.collections), "flier_count": len(artists["fliers"]),
        "annotated_comparisons": len(comparisons),
        "significant_comparisons": int((pairwise["p_adjusted"] < 0.05).sum()),
        "nonsignificant_comparisons": int((pairwise["p_adjusted"] >= 0.05).sum()),
        "alignment": "NOT APPLICABLE: single plot area",
        "box_statistics": [{k: float(b[k]) for k in ("q1", "med", "q3", "whislo", "whishi")}
                           for b in boxes],
        "rendered_pdf_audits": "Run external text and collision audits on final PDF",
    }


def save_outputs(module_name, long_df, stats_df, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = output_dir / f"{module_name}_module_expression_boxplot"
    long_df.sort_values(["temperature", "gene_id"]).to_csv(
        output_dir / f"{module_name}_module_gene_mean_log1p.tsv", sep="\t", index=False)
    stats_df.to_csv(output_dir / f"{module_name}_module_stats.tsv", sep="\t", index=False)
    fig, qa = plot_module(long_df, stats_df, module_name)
    try:
        # Fixed canvas preserves the declared 89 x 108 mm dimensions.
        fig.savefig(stem.with_suffix(".pdf"))
        fig.savefig(stem.with_suffix(".svg"))
        fig.savefig(stem.with_suffix(".png"), dpi=600)
        fig.savefig(stem.with_suffix(".tiff"), dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    finally:
        plt.close(fig)
    stem.with_suffix(".qa.json").write_text(json.dumps(qa, indent=2) + "\n")
    caption = (
        f"{module_name}: module gene expression across 10, 16 and 22 °C. "
        f"n = {qa['gene_count']} genes per temperature; sample counts are 160, 163 and 728. "
        "Each observation is ln(mean TPM + 1), with mean TPM computed across all "
        "samples within the temperature before transformation. Boxes show the "
        "25th-75th percentiles and median; whiskers reach the most extreme values "
        "within 1.5 times the IQR. Individual points and outlier markers are hidden; "
        "all observations remain in box summaries and statistical tests. "
        "Kruskal-Wallis is followed, only when P < 0.05, by Dunn pairwise tests "
        "with Benjamini-Hochberg adjustment across three comparisons within each module. "
        "All tested pairs are annotated: ns, adjusted P >= 0.05; * adjusted P < 0.05, ** < 0.01, *** < 0.001. "
        "If the overall test is nonsignificant, Overall: ns is shown and pairwise tests are not performed. "
        "Exact P values are provided in the statistics TSV. The analysis unit is a gene, "
        "not a biological replicate. The inherited tests do not model matching of genes "
        "across temperatures or dependence among coexpressed genes.\n"
    )
    stem.with_suffix(".caption.txt").write_text(caption)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()
    module_genes = {name: read_module_genes(path) for name, path in MODULE_FILES.items()}
    tables = {temp: read_expression_matrix(path) for temp, path in EXPRESSION_FILES.items()}
    validate_inputs(module_genes, tables)
    for name, genes in module_genes.items():
        long_df = build_module_long_table(name, genes, tables)
        stats_df = run_statistics(long_df, name)
        save_outputs(name, long_df, stats_df, args.output_dir)
        print(f"[OK] {name}: {len(genes)} genes, {len(long_df)} records -> {args.output_dir}")


if __name__ == "__main__":
    main()
