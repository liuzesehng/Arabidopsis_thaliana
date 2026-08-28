#!/usr/bin/env python3
"""Plot CHH counts for AT2G39730 motif hits in motif and flanking regions."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
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
OUTPUT_PREFIX = OUTPUT_DIR / "AT2G39730_CHH_motif_region_counts"

REGION_ORDER = ["in", "flanking"]
SET_TYPE_ORDER = ["negative", "positive"]
LABELS = {"negative": "non-co-exp", "positive": "co-exp"}
REGION_COLORS = {"in": "#2A9D8F", "flanking": "#E76F51"}
EXPECTED_COUNTS = {
    ("in", "negative"): 1,
    ("in", "positive"): 0,
    ("flanking", "negative"): 17,
    ("flanking", "positive"): 1,
}


def assign_region(row: pd.Series) -> str | None:
    """Return plotting region for a methylated motif row."""
    if bool(row["in_motif"]):
        return "in"
    if bool(row["in_window"]) and not bool(row["in_motif"]):
        return "flanking"
    return None


def build_counts(input_tsv: Path) -> pd.DataFrame:
    df = pd.read_csv(input_tsv, sep="\t")

    required_cols = {
        "context",
        "set_type",
        "in_motif",
        "in_window",
        "feature_pos",
        "methylation_feature",
        "shap_value",
    }
    missing_cols = required_cols.difference(df.columns)
    if missing_cols:
        raise ValueError(f"Missing required columns: {sorted(missing_cols)}")

    chh = df[df["context"].eq("CHH")].copy()
    chh = chh[chh["in_window"] | chh["in_motif"]]

    motif_hit_cols = [
        col
        for col in chh.columns
        if col
        not in {
            "feature_pos",
            "methylation_feature",
            "shap_value",
            "in_window",
            "in_motif",
        }
    ]
    chh = (
        chh.groupby(motif_hit_cols, as_index=False, dropna=False)
        .agg({"in_window": "any", "in_motif": "any"})
    )
    chh["region"] = chh.apply(assign_region, axis=1)
    chh = chh[chh["region"].isin(REGION_ORDER)]

    counts = (
        chh.groupby(["region", "set_type"], observed=False)
        .size()
        .rename("count")
        .reindex(
            pd.MultiIndex.from_product(
                [REGION_ORDER, SET_TYPE_ORDER], names=["region", "set_type"]
            ),
            fill_value=0,
        )
        .reset_index()
    )
    counts["label"] = counts["set_type"].map(LABELS)
    return counts[["region", "set_type", "label", "count"]]


def validate_counts(counts: pd.DataFrame) -> None:
    count_map = {
        (row.region, row.set_type): int(row.count)
        for row in counts.itertuples(index=False)
    }
    if count_map != EXPECTED_COUNTS:
        raise ValueError(
            "Unexpected counts after CHH filtering and motif-level deduplication: "
            f"{count_map}; expected {EXPECTED_COUNTS}"
        )


def plot_counts(counts: pd.DataFrame, output_prefix: Path) -> None:
    fig, ax = plt.subplots(figsize=(4.4, 3.4))
    x_positions = list(range(len(SET_TYPE_ORDER)))
    width = 0.28
    bottoms = [0] * len(SET_TYPE_ORDER)

    for region in REGION_ORDER:
        subset = counts[counts["region"].eq(region)].set_index("set_type")
        values = [int(subset.loc[set_type, "count"]) for set_type in SET_TYPE_ORDER]
        ax.bar(
            x_positions,
            values,
            width=width,
            bottom=bottoms,
            label=region,
            facecolor=REGION_COLORS[region],
            edgecolor="#333333",
            linewidth=1.0,
        )
        bottoms = [bottom + value for bottom, value in zip(bottoms, values)]

    ax.set_xticks(x_positions)
    ax.set_xticklabels([LABELS[set_type] for set_type in SET_TYPE_ORDER])
    ax.set_ylim(0, max(max(bottoms) + 2, 3))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", direction="out")

    fig.tight_layout()
    fig.savefig(f"{output_prefix}.png", dpi=300)
    fig.savefig(f"{output_prefix}.pdf")
    plt.close(fig)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    counts = build_counts(INPUT_TSV)
    validate_counts(counts)
    counts.to_csv(f"{OUTPUT_PREFIX}.tsv", sep="\t", index=False)
    plot_counts(counts, OUTPUT_PREFIX)
    print(counts.to_string(index=False))
    print(f"Saved outputs to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
