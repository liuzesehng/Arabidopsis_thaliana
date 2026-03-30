#!/usr/bin/env python3
"""
Plot module-level gene expression boxplots across RCA temperature groups.

Run in the target environment:
    ssh c61
    conda activate XGBoot
    python /datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana/pipline/coex_analysis/4.11.module.py
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal
import scikit_posthocs as sp


BASE_DIR = Path(
    "/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana"
)
OUTPUT_DIR = BASE_DIR / "list/RCA/total_TPM"

MODULE_FILES: Dict[str, Path] = {
    "yellow": BASE_DIR / "list/coex2/total/yellow_module_genes.txt",
    "negative_grey": BASE_DIR / "list/coex2/total/negative_grey_module_genes.txt",
}

EXPRESSION_FILES: Dict[str, Path] = {
    "10C": BASE_DIR / "list/coex2/RCA.10C.tran.txt",
    "16C": BASE_DIR / "list/coex2/RCA.16C.tran.txt",
    "22C": BASE_DIR / "list/coex2/RCA.22C.tran.txt",
}

TEMPERATURE_ORDER = ["10C", "16C", "22C"]
TEMPERATURE_LABELS = {"10C": "10℃", "16C": "16℃", "22C": "22℃"}
PAIR_ORDER: List[Tuple[str, str]] = [("10C", "16C"), ("10C", "22C"), ("16C", "22C")]
PLOT_COLORS = {"10C": "#4C72B0", "16C": "#DD8452", "22C": "#55A868"}


def read_module_genes(module_file: Path) -> List[str]:
    module_df = pd.read_csv(module_file, sep="\t", usecols=[0], header=0)
    genes = module_df.iloc[:, 0].dropna().astype(str).tolist()
    return genes


def read_expression_matrix(expression_file: Path) -> pd.DataFrame:
    expr_df = pd.read_csv(expression_file, sep="\t")
    expr_df = expr_df.rename(columns={expr_df.columns[0]: "gene_id"})
    expr_df["gene_id"] = expr_df["gene_id"].astype(str)
    return expr_df


def build_module_long_table(
    module_name: str,
    module_genes: Sequence[str],
    expression_tables: Dict[str, pd.DataFrame],
) -> pd.DataFrame:
    records: List[pd.DataFrame] = []

    for temperature in TEMPERATURE_ORDER:
        expr_df = expression_tables[temperature]
        subset = expr_df[expr_df["gene_id"].isin(module_genes)].copy()
        subset["mean_expression"] = subset.iloc[:, 1:].mean(axis=1)
        subset["log1p_expression"] = np.log1p(subset["mean_expression"])
        subset["module"] = module_name
        subset["temperature"] = temperature
        records.append(
            subset[["module", "temperature", "gene_id", "mean_expression", "log1p_expression"]]
        )

    long_df = pd.concat(records, ignore_index=True)
    long_df["temperature"] = pd.Categorical(
        long_df["temperature"], categories=TEMPERATURE_ORDER, ordered=True
    )
    return long_df


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


def add_significance_annotations(
    ax: plt.Axes,
    long_df: pd.DataFrame,
    stats_df: pd.DataFrame,
) -> None:
    dunn_rows = stats_df[
        (stats_df["test"] == "Dunn_BH") & (stats_df["p_adjusted"] < 0.05)
    ].copy()
    if dunn_rows.empty:
        return

    ymax = long_df["log1p_expression"].max()
    ymin = long_df["log1p_expression"].min()
    data_span = max(ymax - ymin, 0.2)
    base_height = ymax + data_span * 0.08
    step = data_span * 0.12
    x_positions = {temp: idx for idx, temp in enumerate(TEMPERATURE_ORDER)}

    for idx, (_, row) in enumerate(dunn_rows.iterrows()):
        x1 = x_positions[row["group1"]]
        x2 = x_positions[row["group2"]]
        y = base_height + idx * step
        h = step * 0.35
        ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.4, c="black")
        ax.text(
            (x1 + x2) / 2,
            y + h + step * 0.05,
            row["significance"],
            ha="center",
            va="bottom",
            fontsize=12,
        )

    ax.set_ylim(top=base_height + len(dunn_rows) * step + step * 0.8)


def plot_module(long_df: pd.DataFrame, stats_df: pd.DataFrame, output_path: Path) -> None:
    sns.set_theme(style="whitegrid")
    fig, ax = plt.subplots(figsize=(7, 6), dpi=300)

    sns.boxplot(
        data=long_df,
        x="temperature",
        y="log1p_expression",
        hue="temperature",
        order=TEMPERATURE_ORDER,
        palette=PLOT_COLORS,
        width=0.55,
        fliersize=0,
        linewidth=1.2,
        ax=ax,
        legend=False,
    )
    sns.stripplot(
        data=long_df,
        x="temperature",
        y="log1p_expression",
        order=TEMPERATURE_ORDER,
        color="black",
        alpha=0.45,
        jitter=0.24,
        size=3,
        ax=ax,
    )

    ax.set_xlabel("Temperature", fontsize=12)
    ax.set_ylabel("log1p(mean expression)", fontsize=12)
    ax.set_xticklabels([TEMPERATURE_LABELS[temp] for temp in TEMPERATURE_ORDER])

    add_significance_annotations(ax, long_df, stats_df)

    plt.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def validate_inputs(
    module_genes: Dict[str, Sequence[str]],
    expression_tables: Dict[str, pd.DataFrame],
) -> None:
    expected_module_sizes = {"yellow": 902, "negative_grey": 558}
    expected_sample_sizes = {"10C": 160, "16C": 163, "22C": 728}

    for module_name, expected_size in expected_module_sizes.items():
        observed_size = len(module_genes[module_name])
        if observed_size != expected_size:
            raise ValueError(
                f"Unexpected gene count for {module_name}: {observed_size} != {expected_size}"
            )

    for temperature, expected_size in expected_sample_sizes.items():
        observed_size = expression_tables[temperature].shape[1] - 1
        if observed_size != expected_size:
            raise ValueError(
                f"Unexpected sample count for {temperature}: {observed_size} != {expected_size}"
            )

    reference_genes = set(expression_tables["10C"]["gene_id"])
    for module_name, genes in module_genes.items():
        missing = sorted(set(genes) - reference_genes)
        if missing:
            raise ValueError(f"{module_name} has {len(missing)} genes missing from expression matrix.")


def save_outputs(module_name: str, long_df: pd.DataFrame, stats_df: pd.DataFrame) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    long_path = OUTPUT_DIR / f"{module_name}_module_gene_mean_log1p.tsv"
    stats_path = OUTPUT_DIR / f"{module_name}_module_stats.tsv"
    plot_path = OUTPUT_DIR / f"{module_name}_module_expression_boxplot.png"

    long_df.sort_values(["temperature", "gene_id"]).to_csv(long_path, sep="\t", index=False)
    stats_df.to_csv(stats_path, sep="\t", index=False)
    plot_module(long_df, stats_df, plot_path)


def main() -> None:
    module_genes = {
        module_name: read_module_genes(path) for module_name, path in MODULE_FILES.items()
    }
    expression_tables = {
        temperature: read_expression_matrix(path)
        for temperature, path in EXPRESSION_FILES.items()
    }

    validate_inputs(module_genes, expression_tables)

    for module_name, genes in module_genes.items():
        long_df = build_module_long_table(module_name, genes, expression_tables)
        if long_df["log1p_expression"].isna().any():
            raise ValueError(f"Missing log1p expression values detected for {module_name}.")
        if (long_df["log1p_expression"] < 0).any():
            raise ValueError(f"Negative log1p expression values detected for {module_name}.")

        stats_df = run_statistics(long_df, module_name)
        save_outputs(module_name, long_df, stats_df)

        print(
            f"[OK] {module_name}: {len(genes)} genes, "
            f"{len(long_df)} records, outputs written to {OUTPUT_DIR}"
        )


if __name__ == "__main__":
    main()
