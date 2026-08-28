#!/usr/bin/env python3
"""Plot AT2G39730 CHH co-expression points by symbol and region."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd


BASE_DIR = Path(
    "/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/"
    "Arabidopsis_thaliana/list/RCA"
)
INPUT_TSV = (
    BASE_DIR
    / "AT2G39730_motif_gene_outputs"
    / "total_TPM_AT2G39730_motif_gene_matches.tsv"
)
OUTPUT_DIR = BASE_DIR / "TF_CHH"
OUTPUT_PNG = OUTPUT_DIR / "AT2G39730_CHH_symbol_region_coexp.png"
OUTPUT_PDF = OUTPUT_DIR / "AT2G39730_CHH_symbol_region_coexp.pdf"
OUTPUT_POINTS_TSV = OUTPUT_DIR / "AT2G39730_CHH_symbol_region_coexp_points.tsv"

REGION_ORDER = ["in", "flanking"]
REGION_PRIORITY = {"in": 0, "flanking": 1}
REGION_STYLES = {
    "in": {"facecolor": "#2A9D8F", "edgecolor": "#1E6F66"},
    "flanking": {"facecolor": "#E76F51", "edgecolor": "#A64F3A"},
}
SET_TYPE_ORDER = ["negative", "positive"]
SET_TYPE_LABELS = {"negative": "non-co-exp", "positive": "co-exp"}


def to_bool(value: object) -> bool:
    """Handle real booleans and TSV string booleans consistently."""
    if isinstance(value, bool):
        return value
    if pd.isna(value):
        return False
    return str(value).strip().lower() == "true"


def assign_region(row: pd.Series) -> str | None:
    in_motif = to_bool(row["in_motif"])
    in_window = to_bool(row["in_window"])
    if in_motif:
        return "in"
    if in_window and not in_motif:
        return "flanking"
    return None


def load_points(input_tsv: Path) -> pd.DataFrame:
    df = pd.read_csv(input_tsv, sep="\t")

    required_cols = {"symbol", "context", "set_type", "in_window", "in_motif"}
    missing_cols = required_cols.difference(df.columns)
    if missing_cols:
        raise ValueError(f"Missing required columns: {sorted(missing_cols)}")

    chh = df[df["context"].eq("CHH")].copy()
    chh["region"] = chh.apply(assign_region, axis=1)
    chh = chh[
        chh["region"].isin(REGION_ORDER)
        & chh["set_type"].isin(SET_TYPE_ORDER)
        & chh["symbol"].notna()
    ].copy()

    points = chh.drop_duplicates(subset=["symbol", "region", "set_type"]).copy()
    points["_region_priority"] = points["region"].map(REGION_PRIORITY)
    points = (
        points.sort_values("_region_priority", kind="stable")
        .drop_duplicates(subset=["symbol", "set_type"], keep="first")
        .sort_index()
        .drop(columns="_region_priority")
    )
    points["region"] = pd.Categorical(
        points["region"], categories=REGION_ORDER, ordered=True
    )
    points["set_type"] = pd.Categorical(
        points["set_type"], categories=SET_TYPE_ORDER, ordered=True
    )
    return points.reset_index(drop=True)


def plot_points(points: pd.DataFrame, output_png: Path, output_pdf: Path) -> None:
    symbols = list(dict.fromkeys(points["symbol"].astype(str)))
    symbol_to_x = {symbol: idx for idx, symbol in enumerate(symbols)}
    set_type_to_y = {set_type: idx for idx, set_type in enumerate(SET_TYPE_ORDER)}

    fig_width = max(7.0, 0.38 * len(symbols) + 2.0)
    fig, ax = plt.subplots(figsize=(fig_width, 3.6))
    ax.set_facecolor("#FFFFFF")
    fig.patch.set_facecolor("#FFFFFF")

    marker_size = 52
    for region in REGION_ORDER:
        style = REGION_STYLES[region]
        subset = points[points["region"].astype(str).eq(region)]
        ax.scatter(
            [symbol_to_x[str(symbol)] for symbol in subset["symbol"]],
            [set_type_to_y[str(set_type)] for set_type in subset["set_type"]],
            s=marker_size,
            marker="o",
            facecolors=style["facecolor"],
            edgecolors=style["edgecolor"],
            linewidths=1.1,
            label=region,
            zorder=3,
        )

    ax.set_xlim(-0.6, len(symbols) - 0.4)
    ax.set_ylim(-0.45, len(SET_TYPE_ORDER) - 0.55)
    ax.set_xticks(range(len(symbols)))
    ax.set_xticklabels(symbols, rotation=60, ha="right")
    ax.set_yticks(range(len(SET_TYPE_ORDER)))
    ax.set_yticklabels([SET_TYPE_LABELS[set_type] for set_type in SET_TYPE_ORDER])
    ax.invert_yaxis()
    ax.set_xlabel("symbol")
    ax.set_ylabel("")
    ax.grid(axis="both", color="#E2E6E7", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(axis="both", length=0)
    for spine in ax.spines.values():
        spine.set_color("#B8C2C4")
        spine.set_linewidth(0.8)

    ax.legend(
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0,
    )

    fig.tight_layout()
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    fig.savefig(output_pdf, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    points = load_points(INPUT_TSV)
    points.to_csv(OUTPUT_POINTS_TSV, sep="\t", index=False)
    plot_points(points, OUTPUT_PNG, OUTPUT_PDF)
    print(points[["symbol", "region", "set_type"]].to_string(index=False))
    print(f"Saved point table: {OUTPUT_POINTS_TSV}")
    print(f"Saved plot: {OUTPUT_PNG}")
    print(f"Saved plot: {OUTPUT_PDF}")


if __name__ == "__main__":
    main()
