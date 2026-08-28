#!/usr/bin/env python3
"""Plot CHH motif/window match counts by methylation feature."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


DEFAULT_INPUT = Path(
    "/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/"
    "Arabidopsis_thaliana/list/RCA/AT2G39730_motif_gene_outputs/"
    "total_TPM_AT2G39730_motif_gene_matches.tsv"
)
DEFAULT_OUTPUT_DIR = Path(
    "/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/"
    "Arabidopsis_thaliana/list/RCA/TF_CHH"
)
DEFAULT_OUTPUT_NAME = "AT2G39730_CHH_motif_number.png"
DEFAULT_PDF_NAME = "AT2G39730_CHH_motif_number.pdf"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Count CHH rows where in_motif or in_window is True, then plot "
            "deduplicated methylation_feature + symbol counts."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"Input TSV file. Default: {DEFAULT_INPUT}",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory. Default: {DEFAULT_OUTPUT_DIR}",
    )
    parser.add_argument(
        "--output-name",
        default=DEFAULT_OUTPUT_NAME,
        help=f"Output PNG filename. Default: {DEFAULT_OUTPUT_NAME}",
    )
    parser.add_argument(
        "--pdf-name",
        default=DEFAULT_PDF_NAME,
        help=f"Output PDF filename. Default: {DEFAULT_PDF_NAME}",
    )
    return parser.parse_args()


def parse_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes", "y"})


def load_chh_counts(input_path: Path) -> pd.DataFrame:
    df = pd.read_csv(input_path, sep="\t")
    required_columns = {
        "methylation_feature",
        "symbol",
        "context",
        "in_motif",
        "in_window",
    }
    missing_columns = required_columns.difference(df.columns)
    if missing_columns:
        raise ValueError(
            "Input file is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    is_chh = df["context"].astype(str).str.strip().eq("CHH")
    in_motif = parse_bool(df["in_motif"])
    in_window = parse_bool(df["in_window"])
    hits = df.loc[
        is_chh & (in_motif | in_window), ["methylation_feature", "symbol"]
    ].copy()
    hits = hits.dropna(subset=["methylation_feature", "symbol"])

    if hits.empty:
        return pd.DataFrame(columns=["methylation_feature", "count"])

    hits["methylation_feature"] = hits["methylation_feature"].astype(str).str.strip()
    hits["symbol"] = hits["symbol"].astype(str).str.strip()
    hits = hits.drop_duplicates(subset=["methylation_feature", "symbol"])

    counts = (
        hits.groupby("methylation_feature", as_index=False)
        .size()
        .rename(columns={"size": "count"})
    )
    feature_number = counts["methylation_feature"].str.extract(r"(\d+)$")[0]
    counts["_feature_number"] = pd.to_numeric(feature_number, errors="coerce")
    counts = counts.sort_values(
        ["count", "_feature_number", "methylation_feature"],
        ascending=[False, True, True],
        na_position="last",
    ).drop(columns="_feature_number")
    return counts


def plot_counts(counts: pd.DataFrame, output_paths: list[Path]) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    fig_width = max(10, min(18, 0.55 * max(len(counts), 1) + 5))
    fig, ax = plt.subplots(figsize=(fig_width, 6.5))

    if counts.empty:
        ax.text(
            0.5,
            0.5,
            "No CHH motif/window matches",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=16,
        )
        ax.set_axis_off()
    else:
        palette = sns.color_palette("crest", n_colors=len(counts))
        sns.barplot(
            data=counts,
            x="methylation_feature",
            y="count",
            hue="methylation_feature",
            palette=palette,
            legend=False,
            ax=ax,
        )
        ax.set_xlabel("")
        ax.set_ylabel("Motif number")
        ax.set_ylim(0, counts["count"].max() * 1.15 + 1)
        ax.tick_params(axis="x", labelrotation=45, labelsize=10)
        for label in ax.get_xticklabels():
            label.set_horizontalalignment("right")
        sns.despine(ax=ax)

    fig.tight_layout()
    for output_path in output_paths:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    png_path = output_dir / args.output_name
    pdf_path = output_dir / args.pdf_name

    counts = load_chh_counts(args.input)
    plot_counts(counts, [png_path, pdf_path])

    print(f"Input: {args.input}")
    print(f"PNG output: {png_path}")
    print(f"PDF output: {pdf_path}")
    total_count = int(counts["count"].sum()) if not counts.empty else 0
    print(f"Total CHH methylation_feature + symbol pairs: {total_count}")
    print(f"CHH methylation_feature count: {len(counts)}")


if __name__ == "__main__":
    main()
