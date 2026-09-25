#!/usr/bin/env python3
"""Recreate TF temperature boxplots without points; run with conda activate XGBoot.

Each sample column contributes one observation. No values are excluded from
box statistics or tests. Sample independence follows the reference analysis;
these tables do not establish biological versus technical replicate identity.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cbook import boxplot_stats
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import scikit_posthocs as sp
from scipy.stats import kruskal

BASE_DIR = Path(__file__).resolve().parents[2]
OUTPUT_DIR = BASE_DIR / "list/RCA/tf_boxplot"
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
EXPECTED_SAMPLE_COUNTS = {"10C": 160, "16C": 163, "22C": 728}
COLORS = ["#4477AA", "#EE9944", "#228877"]
EXPECTED_GENES = {
    "coexpression": ["AT1G01060", "AT2G22540", "AT3G09600", "AT5G02840", "AT5G17300"],
    "antagonistic": ["AT5G10140"],
}

def read_tf_gene_ids(tf_ids_file: Path) -> List[str]:
    tf_df = pd.read_csv(tf_ids_file, sep="\t")
    if tf_df["motif_id"].isna().any() or tf_df["motif_id"].astype(str).str.strip().eq("").any():
        raise ValueError(f"{tf_ids_file.name}: missing TF identifier")
    genes = tf_df["motif_id"].astype(str).drop_duplicates().tolist()
    return genes


def read_annotation_map(annotation_file: Path, gene_ids: Sequence[str]) -> pd.DataFrame:
    annotation_df = pd.read_csv(annotation_file, sep="\t")
    annotation_df["Gene_ID"] = annotation_df["Gene_ID"].astype(str)
    annotation_df["Symbol"] = annotation_df["Symbol"].fillna("").astype(str).str.strip()
    conflicts = annotation_df.groupby("Gene_ID")["Symbol"].nunique()
    if (conflicts.loc[conflicts.index.isin(gene_ids)] > 1).any():
        raise ValueError(f"{annotation_file.name}: conflicting TF symbols")
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
    with expression_file.open() as handle:
        header = next(csv.reader(handle, delimiter="\t"))
    if len(header) < 2 or any(not x.strip() for x in header):
        raise ValueError(f"{expression_file.name}: empty gene/sample column name")
    if len(set(header)) != len(header) or "gene_id" in header[1:]:
        raise ValueError(f"{expression_file.name}: duplicate or reserved sample identifier")
    expr_df = pd.read_csv(expression_file, sep="\t")
    expr_df = expr_df.rename(columns={expr_df.columns[0]: "gene_id"})
    if expr_df["gene_id"].isna().any() or expr_df["gene_id"].duplicated().any():
        raise ValueError(f"{expression_file.name}: missing or duplicate gene identifier")
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
        kw_stat, kw_p = (0.0, 1.0) if np.ptp(np.concatenate(groups)) == 0 else kruskal(*groups)

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


def prepare_set(set_name, matrices):
    config = TF_SET_CONFIGS[set_name]
    genes = read_tf_gene_ids(Path(config['tf_ids_file']))
    if genes != EXPECTED_GENES[set_name]:
        raise ValueError(f'{set_name}: unexpected TF identities or order')
    annotations = read_annotation_map(Path(config['annotation_file']), genes)
    if annotations.symbol.duplicated().any():
        raise ValueError(f'{set_name}: duplicate TF symbols')
    for temp, matrix in matrices.items():
        if matrix.shape[1] - 1 != EXPECTED_SAMPLE_COUNTS[temp]:
            raise ValueError(f'{temp}: unexpected sample count')
        selected = matrix.loc[matrix.gene_id.isin(genes)]
        if len(selected) != len(genes) or set(selected.gene_id) != set(genes):
            raise ValueError(f'{temp}: missing or duplicate TF genes')
        values = selected.drop(columns='gene_id').apply(pd.to_numeric, errors='coerce').to_numpy()
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError(f'{set_name}/{temp}: TPM must be numeric, finite and nonnegative')
    long_df = build_tf_long_table(set_name, genes, annotations, matrices)
    if long_df.duplicated(['gene_id', 'temperature', 'sample_id']).any():
        raise ValueError(f'{set_name}: duplicate observations')
    stats = run_symbol_statistics(long_df, annotations, set_name)
    symbols = annotations.set_index('gene_id').loc[genes, 'symbol'].tolist()
    return long_df, stats, symbols


def plot_tf_set(long_df, stats_df, symbols, set_name):
    """No point markers; full observations enter quartiles and rank tests."""
    plt.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['DejaVu Sans'],
        'font.size': 7, 'axes.labelsize': 8, 'axes.titlesize': 9,
        'pdf.fonttype': 42, 'svg.fonttype': 'none',
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
        'ytick.major.width': 0.6, 'legend.frameon': False,
    })
    single = len(symbols) == 1
    width_mm = 89 if single else 183
    fig, ax = plt.subplots(figsize=(width_mm / 25.4, 108 / 25.4))
    fig.subplots_adjust(left=0.17 if single else 0.085, right=0.97, bottom=0.27, top=0.80)
    offsets = np.array([-0.24, 0, 0.24])
    tops = {}
    for i, symbol in enumerate(symbols):
        values = [long_df.loc[(long_df.symbol == symbol) & (long_df.temperature == t),
                             'log1p_TPM'].to_numpy() for t in TEMPERATURE_ORDER]
        boxes = boxplot_stats(values, whis=1.5)
        tops[symbol] = max(b['whishi'] for b in boxes)
        artists = ax.bxp(boxes, positions=i + offsets, widths=0.18,
                         showfliers=False, showmeans=False, patch_artist=True,
                         manage_ticks=False,
                         boxprops={'edgecolor': '#333333', 'linewidth': 0.65},
                         medianprops={'color': '#222222', 'linewidth': 1.0},
                         whiskerprops={'color': '#444444', 'linewidth': 0.65},
                         capprops={'color': '#444444', 'linewidth': 0.65})
        for patch, color in zip(artists['boxes'], COLORS):
            patch.set_facecolor(color)
            patch.set_alpha(0.82)
    span = max(max(tops.values()), 1.0)
    ylim = span * 1.08
    for i, symbol in enumerate(symbols):
        comparisons = stats_df.loc[(stats_df.symbol == symbol) &
                                   (stats_df.test == 'Dunn_BH')]
        level = 0
        for t1, t2 in [('10C', '16C'), ('16C', '22C'), ('10C', '22C')]:
            row = comparisons.loc[(comparisons.group1 == t1) & (comparisons.group2 == t2)]
            if row.empty:
                continue
            x1, x2 = [i + offsets[TEMPERATURE_ORDER.index(t)] for t in (t1, t2)]
            y = tops[symbol] + span * (0.09 + 0.14 * level)
            h = span * 0.025
            ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], color='#333333', lw=0.6)
            ax.text((x1+x2)/2, y+h+span*0.012, row.iloc[0].significance,
                    ha='center', va='bottom', fontsize=7)
            ylim = max(ylim, y + span * 0.14)
            level += 1
    ax.set_ylim(0, ylim)
    ax.set_xlim(-0.55, len(symbols)-0.45)
    ax.set_xticks(range(len(symbols)), symbols)
    ax.set_ylabel('ln(TPM + 1)')
    ax.set_xlabel('Transcription factor', labelpad=8)
    ax.tick_params(axis='x', length=0, pad=7)
    fig.text(0.5, 0.95, 'Antagonistic TF' if single else 'Coexpression TFs',
             ha='center', va='top', fontsize=9)
    fig.legend(handles=[Patch(facecolor=c, edgecolor='#333333', linewidth=0.6,
                              label=f'{t[:-1]} °C')
                        for t, c in zip(TEMPERATURE_ORDER, COLORS)],
               loc='upper right', bbox_to_anchor=(0.97, 0.89), ncol=3,
               handlelength=1.1, columnspacing=1.2 if not single else 0.8,
               fontsize=7 if not single else 6)
    fig.text(
        0.5, 0.095,
        'n = 160 (10 °C), 163 (16 °C), 728 (22 °C) per TF.\n'
        'n denotes sample columns. Boxes: median and IQR;\n'
        'whiskers: 1.5 × IQR. Points are hidden.',
        ha='center', va='top', fontsize=6, linespacing=1.5,
    )
    fig.canvas.draw()
    return fig, ax


def save_outputs(set_name, long_df, stats_df, symbols, output_dir):
    prefix = output_dir / f'{set_name}_tf_boxplot'
    long_df.sort_values(['symbol', 'temperature', 'sample_id']).to_csv(
        output_dir / f'{set_name}_tf_expression_long.tsv', sep='\t', index=False)
    stats_df.to_csv(output_dir / f'{set_name}_tf_statistics.tsv', sep='\t', index=False)
    fig, ax = plot_tf_set(long_df, stats_df, symbols, set_name)
    # Single plotting axes: multi-panel alignment is not applicable.
    prefix.with_suffix('.alignment.json').write_text(json.dumps({
        'verdict': 'NOT APPLICABLE', 'reason': 'Single plotting axes',
        'axes_count': len(fig.axes), 'axes_bounds': list(ax.get_position().bounds),
    }, indent=2) + '\n')
    fig.savefig(prefix.with_suffix('.pdf'))
    fig.savefig(prefix.with_suffix('.svg'))
    fig.savefig(prefix.with_suffix('.png'), dpi=600)
    fig.savefig(prefix.with_suffix('.tiff'), dpi=600, pil_kwargs={'compression': 'tiff_lzw'})
    plt.close(fig)
    prefix.with_suffix('.caption.txt').write_text(
        f'{set_name.capitalize()} TF expression across RCA temperature groups. '
        f'TF order: {", ".join(symbols)}. Expression is ln(TPM + 1). '
        'Boxes show Q1–Q3 and the median; whiskers end at the most extreme observations '
        'within 1.5×IQR. Raw points and outlier markers are hidden; all values enter '
        'box statistics and tests. Sample sizes (n, sample columns) are 160, 163 and 728 '
        'at 10, 16 and 22 °C, respectively. Biological/technical replicate identity and sample '
        'independence are not established by these matrices; independence follows '
        'the reference analysis. Per-gene Kruskal–Wallis tests are followed, only '
        'when P < 0.05, by two-sided Dunn tests with Benjamini–Hochberg correction '
        'over the three temperature comparisons within each gene. All tested pairwise '
        'comparisons are annotated: ns, adjusted P >= 0.05; * adjusted P < 0.05, '
        '** < 0.01, *** < 0.001. '
        'Exact P values and all tested comparisons are in the statistics TSV.\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()
    try:
        matrices = {t: read_expression_matrix(p) for t, p in EXPRESSION_FILES.items()}
        prepared = {name: prepare_set(name, matrices) for name in TF_SET_CONFIGS}
        args.output_dir.mkdir(parents=True, exist_ok=True)
        for name, (long_df, stats, symbols) in prepared.items():
            save_outputs(name, long_df, stats, symbols, args.output_dir)
            print(f'[OK] {name}: {len(long_df)} observations; {args.output_dir}')
    except (ValueError, KeyError, OSError) as exc:
        parser.exit(1, f'[ERROR] {exc}\n')


if __name__ == '__main__':
    main()
