#!/usr/bin/env python3
"""
Plot TF-level TPM boxplots across RCA temperature groups.

Run in the target environment:
    ssh c61
    conda activate XGBoot
    python /datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana/pipline/coex_analysis/4.12.module.contibution.py
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scikit_posthocs as sp
import seaborn as sns
from matplotlib.legend_handler import HandlerPatch
from matplotlib.patches import Rectangle
from scipy.stats import kruskal


BASE_DIR = Path(
    "/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana"
)
OUTPUT_DIR = BASE_DIR / "list/RCA/total_TPM"

TF_SET_CONFIGS: Dict[str, Dict[str, object]] = {
    "coexpression": {
        "tf_ids_file": BASE_DIR / "list/coex2/total/gene.TF_ids.txt",
        "annotation_file": BASE_DIR / "list/coex2/total/total.gene.TF_ids.annotated.txt",
        "expected_gene_count": 5,
        "plot_title": "Coexpression TF TPM Across RCA Temperatures",
        "long_filename": "coexpression_tf_expression_long.tsv",
        "stats_filename": "coexpression_tf_statistics.tsv",
        "plot_filename": "coexpression_tf_boxplot.png",
    },
    "antagonistic": {
        "tf_ids_file": BASE_DIR / "list/coex2/total/gene.TF_negative_ids.txt",
        "annotation_file": BASE_DIR / "list/coex2/total/total.negative.gene.TF_ids.annotated.txt",
        "expected_gene_count": 1,
        "plot_title": "Antagonistic TF TPM Across RCA Temperatures",
        "long_filename": "antagonistic_tf_expression_long.tsv",
        "stats_filename": "antagonistic_tf_statistics.tsv",
        "plot_filename": "antagonistic_tf_boxplot.png",
    },
}

EXPRESSION_FILES: Dict[str, Path] = {
    "10C": BASE_DIR / "list/coex2/RCA.10C.tran.txt",
    "16C": BASE_DIR / "list/coex2/RCA.16C.tran.txt",
    "22C": BASE_DIR / "list/coex2/RCA.22C.tran.txt",
}

TEMPERATURE_ORDER = ["10C", "16C", "22C"]
TEMPERATURE_LABELS = {"10C": "10℃", "16C": "16℃", "22C": "22℃"}
TEMPERATURE_DISPLAY_ORDER = ["10℃", "16℃", "22℃"]
PAIR_ORDER: List[Tuple[str, str]] = [("10C", "16C"), ("10C", "22C"), ("16C", "22C")]
PLOT_COLORS = {"10℃": "#4C72B0", "16℃": "#DD8452", "22℃": "#55A868"}
EXPECTED_SAMPLE_COUNTS = {"10C": 160, "16C": 163, "22C": 728}


class VerticalPatchHandler(HandlerPatch):
    def create_artists(
        self,
        legend,
        orig_handle,
        xdescent,
        ydescent,
        width,
        height,
        fontsize,
        trans,
    ):
        rect_width = width * 0.42
        rect_height = height * 0.95
        x = xdescent + (width - rect_width) / 2
        y = ydescent + (height - rect_height) / 2
        patch = Rectangle(
            (x, y),
            rect_width,
            rect_height,
            facecolor=orig_handle.get_facecolor(),
            edgecolor=orig_handle.get_edgecolor(),
            linewidth=orig_handle.get_linewidth(),
            transform=trans,
        )
        return [patch]


def read_tf_gene_ids(tf_ids_file: Path) -> List[str]:
    tf_df = pd.read_csv(tf_ids_file, sep="\t")
    genes = tf_df["motif_id"].dropna().astype(str).drop_duplicates().tolist()
    return genes


def read_annotation_map(annotation_file: Path, gene_ids: Sequence[str]) -> pd.DataFrame:
    annotation_df = pd.read_csv(annotation_file, sep="\t")
    annotation_df["Gene_ID"] = annotation_df["Gene_ID"].astype(str)
    annotation_df["Symbol"] = annotation_df["Symbol"].fillna("").astype(str).str.strip()
    annotation_df = annotation_df.drop_duplicates(subset=["Gene_ID"]).copy()

    missing_genes = sorted(set(gene_ids) - set(annotation_df["Gene_ID"]))
    if missing_genes:
        raise ValueError(
            f"{annotation_file.name} is missing {len(missing_genes)} TF annotations: "
            f"{', '.join(missing_genes)}"
        )

    subset = annotation_df.loc[annotation_df["Gene_ID"].isin(gene_ids), ["Gene_ID", "Symbol"]].copy()
    missing_symbol_mask = subset["Symbol"].eq("")
    if missing_symbol_mask.any():
        for gene_id in subset.loc[missing_symbol_mask, "Gene_ID"]:
            print(f"[WARN] Missing Symbol for {gene_id}; fallback to Gene_ID.")
        subset.loc[missing_symbol_mask, "Symbol"] = subset.loc[missing_symbol_mask, "Gene_ID"]

    return subset.rename(columns={"Gene_ID": "gene_id", "Symbol": "symbol"})


def read_expression_matrix(expression_file: Path) -> pd.DataFrame:
    expr_df = pd.read_csv(expression_file, sep="\t")
    expr_df = expr_df.rename(columns={expr_df.columns[0]: "gene_id"})
    expr_df["gene_id"] = expr_df["gene_id"].astype(str)
    return expr_df


def build_tf_long_table(
    set_name: str,
    gene_ids: Sequence[str],
    annotation_df: pd.DataFrame,
    expression_tables: Dict[str, pd.DataFrame],
) -> pd.DataFrame:
    records: List[pd.DataFrame] = []
    gene_order_map = {gene_id: idx for idx, gene_id in enumerate(gene_ids)}

    for temperature in TEMPERATURE_ORDER:
        expr_df = expression_tables[temperature]
        subset = expr_df.loc[expr_df["gene_id"].isin(gene_ids)].copy()
        subset["gene_order"] = subset["gene_id"].map(gene_order_map)
        subset = subset.sort_values("gene_order").drop(columns="gene_order")
        subset = subset.merge(annotation_df, on="gene_id", how="left", validate="many_to_one")

        long_df = subset.melt(
            id_vars=["gene_id", "symbol"],
            var_name="sample_id",
            value_name="TPM",
        )
        long_df["set_name"] = set_name
        long_df["temperature"] = temperature
        long_df["temperature_label"] = TEMPERATURE_LABELS[temperature]
        long_df["TPM"] = pd.to_numeric(long_df["TPM"], errors="coerce")
        long_df["log1p_TPM"] = np.log1p(long_df["TPM"])
        records.append(
            long_df[
                [
                    "set_name",
                    "gene_id",
                    "symbol",
                    "temperature",
                    "temperature_label",
                    "sample_id",
                    "TPM",
                    "log1p_TPM",
                ]
            ]
        )

    long_df = pd.concat(records, ignore_index=True)
    symbol_order = annotation_df.set_index("gene_id").loc[list(gene_ids), "symbol"].tolist()
    long_df["temperature"] = pd.Categorical(
        long_df["temperature"], categories=TEMPERATURE_ORDER, ordered=True
    )
    long_df["temperature_label"] = pd.Categorical(
        long_df["temperature_label"], categories=TEMPERATURE_DISPLAY_ORDER, ordered=True
    )
    long_df["symbol"] = pd.Categorical(long_df["symbol"], categories=symbol_order, ordered=True)
    return long_df


def p_to_stars(p_value: float) -> str:
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return "ns"


def run_symbol_statistics(long_df: pd.DataFrame, annotation_df: pd.DataFrame, set_name: str) -> pd.DataFrame:
    stat_rows: List[Dict[str, object]] = []

    for row in annotation_df.itertuples(index=False):
        symbol_df = long_df.loc[long_df["gene_id"] == row.gene_id].copy()
        for temperature in TEMPERATURE_ORDER:
            temp_df = symbol_df.loc[symbol_df["temperature"] == temperature]
            stat_rows.append(
                {
                    "set_name": set_name,
                    "symbol": row.symbol,
                    "gene_id": row.gene_id,
                    "record_type": "summary",
                    "temperature": temperature,
                    "test": "",
                    "group1": "",
                    "group2": "",
                    "group3": "",
                    "sample_count": int(len(temp_df)),
                    "statistic": np.nan,
                    "p_value": np.nan,
                    "p_adjusted": np.nan,
                    "significance": "",
                    "TPM_mean": temp_df["TPM"].mean(),
                    "TPM_median": temp_df["TPM"].median(),
                    "log1p_TPM_mean": temp_df["log1p_TPM"].mean(),
                    "log1p_TPM_median": temp_df["log1p_TPM"].median(),
                }
            )

        groups = [
            symbol_df.loc[symbol_df["temperature"] == temperature, "log1p_TPM"].to_numpy()
            for temperature in TEMPERATURE_ORDER
        ]
        kw_stat, kw_p = kruskal(*groups)

        stat_rows.append(
            {
                "set_name": set_name,
                "symbol": row.symbol,
                "gene_id": row.gene_id,
                "record_type": "test",
                "temperature": "",
                "test": "Kruskal-Wallis",
                "group1": "10C",
                "group2": "16C",
                "group3": "22C",
                "sample_count": np.nan,
                "statistic": kw_stat,
                "p_value": kw_p,
                "p_adjusted": kw_p,
                "significance": p_to_stars(kw_p) if kw_p < 0.05 else "ns",
                "TPM_mean": np.nan,
                "TPM_median": np.nan,
                "log1p_TPM_mean": np.nan,
                "log1p_TPM_median": np.nan,
            }
        )

        if kw_p < 0.05:
            dunn_df = sp.posthoc_dunn(
                symbol_df,
                val_col="log1p_TPM",
                group_col="temperature",
                p_adjust="fdr_bh",
            )
            for group1, group2 in PAIR_ORDER:
                p_adj = float(dunn_df.loc[group1, group2])
                stat_rows.append(
                    {
                        "set_name": set_name,
                        "symbol": row.symbol,
                        "gene_id": row.gene_id,
                        "record_type": "test",
                        "temperature": "",
                        "test": "Dunn_BH",
                        "group1": group1,
                        "group2": group2,
                        "group3": "",
                        "sample_count": np.nan,
                        "statistic": np.nan,
                        "p_value": np.nan,
                        "p_adjusted": p_adj,
                        "significance": p_to_stars(p_adj) if p_adj < 0.05 else "ns",
                        "TPM_mean": np.nan,
                        "TPM_median": np.nan,
                        "log1p_TPM_mean": np.nan,
                        "log1p_TPM_median": np.nan,
                    }
                )

    stats_df = pd.DataFrame(stat_rows)
    stats_df["symbol"] = pd.Categorical(
        stats_df["symbol"], categories=annotation_df["symbol"].tolist(), ordered=True
    )
    stats_df["record_type"] = pd.Categorical(
        stats_df["record_type"], categories=["summary", "test"], ordered=True
    )
    stats_df["temperature"] = pd.Categorical(
        stats_df["temperature"], categories=TEMPERATURE_ORDER, ordered=True
    )
    return stats_df.sort_values(
        ["symbol", "record_type", "temperature", "test", "group1", "group2"]
    ).reset_index(drop=True)


def add_symbol_significance_annotations(
    ax: plt.Axes,
    long_df: pd.DataFrame,
    stats_df: pd.DataFrame,
    symbol_order: Sequence[str],
    hue_offsets: Dict[str, float] | None = None,
) -> None:
    dunn_rows = stats_df[
        (stats_df["test"] == "Dunn_BH") & (stats_df["p_adjusted"] < 0.05)
    ].copy()
    if dunn_rows.empty:
        return

    if hue_offsets is None:
        hue_offsets = {"10℃": -0.266, "16℃": 0.0, "22℃": 0.266}
    base_ymax = float(long_df["log1p_TPM"].max())
    base_ymin = float(long_df["log1p_TPM"].min())
    data_span = max(base_ymax - base_ymin, 0.2)
    group_max = long_df.groupby("symbol", observed=False)["log1p_TPM"].max().to_dict()
    row_counts = dunn_rows.groupby("symbol", observed=False).size().to_dict()
    per_symbol_index = {symbol: 0 for symbol in symbol_order}

    top_limit = base_ymax
    for symbol, count in row_counts.items():
        symbol_max = float(group_max.get(symbol, base_ymax))
        top_limit = max(top_limit, symbol_max + data_span * (0.18 + count * 0.16))

    for _, row in dunn_rows.iterrows():
        symbol = row["symbol"]
        symbol_idx = symbol_order.index(symbol)
        level = per_symbol_index[symbol]
        per_symbol_index[symbol] += 1
        x1 = symbol_idx + hue_offsets[TEMPERATURE_LABELS[row["group1"]]]
        x2 = symbol_idx + hue_offsets[TEMPERATURE_LABELS[row["group2"]]]
        symbol_max = float(group_max.get(symbol, base_ymax))
        y = symbol_max + data_span * (0.08 + level * 0.14)
        h = data_span * 0.04
        ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.2, c="black")
        ax.text((x1 + x2) / 2, y + h + data_span * 0.01, row["significance"], ha="center", va="bottom", fontsize=11)

    ax.set_ylim(top=top_limit)


def add_temperature_legend(ax: plt.Axes) -> None:
    legend_handles = [
        Rectangle((0, 0), 1, 1, facecolor=PLOT_COLORS[temp], edgecolor="#444444", linewidth=1.0)
        for temp in TEMPERATURE_DISPLAY_ORDER
    ]
    ax.legend(
        legend_handles,
        TEMPERATURE_DISPLAY_ORDER,
        title="Temperature",
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        handler_map={Rectangle: VerticalPatchHandler()},
        handlelength=1.2,
        handleheight=1.8,
        borderaxespad=0.0,
    )


def plot_tf_set(
    long_df: pd.DataFrame,
    stats_df: pd.DataFrame,
    symbol_order: Sequence[str],
    title: str,
    output_path: Path,
) -> None:
    sns.set_theme(style="whitegrid")
    single_symbol = len(symbol_order) == 1
    fig_width = max(6.5, 1.4 * len(symbol_order) + 2.2)
    if single_symbol:
        fig_width = 6.6
    fig, ax = plt.subplots(figsize=(fig_width, 6), dpi=300)

    if single_symbol:
        center = 0.0
        hue_offsets = {"10℃": -0.42, "16℃": 0.0, "22℃": 0.42}
        positions = [center + hue_offsets[temp] for temp in TEMPERATURE_DISPLAY_ORDER]
        grouped_values = [
            long_df.loc[long_df["temperature_label"] == temp, "log1p_TPM"].to_numpy()
            for temp in TEMPERATURE_DISPLAY_ORDER
        ]

        box = ax.boxplot(
            grouped_values,
            positions=positions,
            widths=0.24,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "#444444", "linewidth": 1.2},
            whiskerprops={"color": "#555555", "linewidth": 1.1},
            capprops={"color": "#555555", "linewidth": 1.1},
            boxprops={"edgecolor": "#555555", "linewidth": 1.1},
        )
        for patch, temp in zip(box["boxes"], TEMPERATURE_DISPLAY_ORDER):
            patch.set_facecolor(PLOT_COLORS[temp])
            patch.set_alpha(0.9)

        rng = np.random.default_rng(20260327)
        for pos, temp in zip(positions, TEMPERATURE_DISPLAY_ORDER):
            temp_df = long_df.loc[long_df["temperature_label"] == temp]
            jitter = rng.uniform(-0.08, 0.08, size=len(temp_df))
            ax.scatter(
                np.full(len(temp_df), pos) + jitter,
                temp_df["log1p_TPM"].to_numpy(),
                color="black",
                alpha=0.32,
                s=7,
                linewidths=0,
                zorder=3,
            )

        ax.set_xticks([center])
        ax.set_xticklabels(symbol_order)
        ax.set_xlim(-0.72, 0.72)
        add_symbol_significance_annotations(
            ax,
            long_df,
            stats_df,
            list(symbol_order),
            hue_offsets=hue_offsets,
        )
    else:
        sns.boxplot(
            data=long_df,
            x="symbol",
            y="log1p_TPM",
            hue="temperature_label",
            order=symbol_order,
            hue_order=TEMPERATURE_DISPLAY_ORDER,
            palette=PLOT_COLORS,
            width=0.72,
            fliersize=0,
            linewidth=1.1,
            ax=ax,
        )
        sns.stripplot(
            data=long_df,
            x="symbol",
            y="log1p_TPM",
            hue="temperature_label",
            order=symbol_order,
            hue_order=TEMPERATURE_DISPLAY_ORDER,
            dodge=True,
            jitter=0.2,
            alpha=0.38,
            size=2.3,
            palette=["black", "black", "black"],
            ax=ax,
        )
        ax.margins(x=0.03)
        add_symbol_significance_annotations(ax, long_df, stats_df, list(symbol_order))

    add_temperature_legend(ax)
    ax.set_xlabel("Symbol", fontsize=12)
    ax.set_ylabel("ln(TPM+1)", fontsize=12)
    ax.tick_params(axis="x", rotation=0)

    fig.subplots_adjust(right=0.8)
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def validate_inputs(
    set_name: str,
    gene_ids: Sequence[str],
    annotation_df: pd.DataFrame,
    expression_tables: Dict[str, pd.DataFrame],
) -> None:
    expected_gene_count = int(TF_SET_CONFIGS[set_name]["expected_gene_count"])
    if len(gene_ids) != expected_gene_count:
        raise ValueError(f"Unexpected TF count for {set_name}: {len(gene_ids)} != {expected_gene_count}")

    if annotation_df["gene_id"].nunique() != expected_gene_count:
        raise ValueError(
            f"Unexpected annotation count for {set_name}: "
            f"{annotation_df['gene_id'].nunique()} != {expected_gene_count}"
        )

    for temperature, expected_sample_count in EXPECTED_SAMPLE_COUNTS.items():
        observed_sample_count = expression_tables[temperature].shape[1] - 1
        if observed_sample_count != expected_sample_count:
            raise ValueError(
                f"Unexpected sample count for {temperature}: "
                f"{observed_sample_count} != {expected_sample_count}"
            )

    for temperature, expr_df in expression_tables.items():
        missing_genes = sorted(set(gene_ids) - set(expr_df["gene_id"]))
        if missing_genes:
            raise ValueError(
                f"{set_name} has {len(missing_genes)} TFs missing from {temperature}: "
                f"{', '.join(missing_genes)}"
            )


def save_outputs(
    set_name: str,
    long_df: pd.DataFrame,
    stats_df: pd.DataFrame,
    symbol_order: Sequence[str],
) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    config = TF_SET_CONFIGS[set_name]

    long_path = OUTPUT_DIR / str(config["long_filename"])
    stats_path = OUTPUT_DIR / str(config["stats_filename"])
    plot_path = OUTPUT_DIR / str(config["plot_filename"])

    long_df.sort_values(["symbol", "temperature", "sample_id"]).to_csv(long_path, sep="\t", index=False)
    stats_df.to_csv(stats_path, sep="\t", index=False)
    plot_tf_set(long_df, stats_df, symbol_order, str(config["plot_title"]), plot_path)


def main() -> None:
    expression_tables = {
        temperature: read_expression_matrix(path)
        for temperature, path in EXPRESSION_FILES.items()
    }

    for set_name, config in TF_SET_CONFIGS.items():
        gene_ids = read_tf_gene_ids(Path(config["tf_ids_file"]))
        annotation_df = read_annotation_map(Path(config["annotation_file"]), gene_ids)
        validate_inputs(set_name, gene_ids, annotation_df, expression_tables)

        long_df = build_tf_long_table(set_name, gene_ids, annotation_df, expression_tables)
        if long_df["TPM"].isna().any():
            raise ValueError(f"Missing TPM values detected for {set_name}.")
        if (long_df["TPM"] < 0).any():
            raise ValueError(f"Negative TPM values detected for {set_name}.")
        if long_df["log1p_TPM"].isna().any():
            raise ValueError(f"Missing log1p TPM values detected for {set_name}.")

        stats_df = run_symbol_statistics(long_df, annotation_df, set_name)
        symbol_order = annotation_df.set_index("gene_id").loc[gene_ids, "symbol"].tolist()
        save_outputs(set_name, long_df, stats_df, symbol_order)

        print(
            f"[OK] {set_name}: {len(gene_ids)} TFs, {len(long_df)} records, "
            f"outputs written to {OUTPUT_DIR}"
        )


if __name__ == "__main__":
    main()
