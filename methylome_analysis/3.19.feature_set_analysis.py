from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import gridspec


BASE_DIR = Path(
    "/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana/list/xgboot/TPM_4.5_all"
)
FEATURE_DATA_DIR = BASE_DIR / "feature_data_extraction"
OUTPUT_DIR = BASE_DIR / "temperature_meth_plots"
GFF_PATH = Path(
    "/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana/refgen/GCF_000001735.4_TAIR10.1_genomic.gff"
)

DATASETS = ["A_TPM", "A_per_TPM", "B_TPM", "B_per_TPM", "total_TPM"]
TYPE_ORDER = ["CHH", "CHG", "CG", "SNP"]
PLOT_TYPES = ["CHH", "CHG", "CG"]
DATASET_ORDER = ["A_TPM", "A_per_TPM", "B_TPM", "B_per_TPM", "total_TPM"]
UPSET_ORDER = DATASET_ORDER + TYPE_ORDER
TOP_PERCENT = 0.10
DATASET_LABELS = {
    "A_TPM": "α",
    "A_per_TPM": "α%",
    "B_TPM": "β",
    "B_per_TPM": "β%",
    "total_TPM": "total",
}
TYPE_LABELS = {"CHH": "CHH", "CHG": "CHG", "CG": "CG", "SNP": "SNP"}
REGION_ORDER = ["upstream", "gene body", "downstream"]
REGION_FILENAME = {
    "upstream": "upstream",
    "gene body": "gene_body",
    "downstream": "downstream",
}
REGION_COLORS = {
    "upstream": "#5B8FF9",
    "gene body": "#5AD8A6",
    "downstream": "#F6BD16",
}
TYPE_COLORS = {
    "CHH": "#C44E52",
    "CHG": "#4E79A7",
    "CG": "#59A14F",
    "SNP": "#B07AA1",
}


@dataclass(frozen=True)
class GeneRegion:
    chrom: str
    strand: str
    gene_start: int
    gene_end: int
    upstream_start: int
    upstream_end: int
    downstream_start: int
    downstream_end: int


def parse_feature_importances(path: Path, top_percent: float = TOP_PERCENT) -> Dict[str, float]:
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    capture = False
    feature_values: Dict[str, float] = {}

    for raw_line in lines:
        line = raw_line.rstrip()
        if "所有特征及其平均|SHAP|值" in line:
            capture = True
            continue
        if not capture:
            continue
        if not line.strip():
            continue
        if line.startswith("SHAP特征重要性分析"):
            break
        match = re.match(
            r"^\s*([A-Za-z]+_\d+|snp_\d+)\s*:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*$",
            line,
        )
        if match:
            feature_values[match.group(1)] = float(match.group(2))

    if not feature_values:
        raise ValueError(f"Failed to parse SHAP feature values from {path}")
    ranked_features = [(feature, value) for feature, value in feature_values.items() if value > 0]
    ranked_features = sorted(ranked_features, key=lambda item: item[1], reverse=True)
    top_n = max(1, int(len(ranked_features) * top_percent))
    ranked_features = ranked_features[:top_n]
    return dict(ranked_features)


def load_dataset_coords(dataset: str) -> pd.DataFrame:
    frames: List[pd.DataFrame] = []

    for feature_type in ["CHH", "CHG", "CG"]:
        path = FEATURE_DATA_DIR / dataset / f"{dataset}_{feature_type}_data.csv"
        df = pd.read_csv(path)
        renamed = df.rename(columns={"Feature": "feature"})
        renamed["type"] = feature_type
        renamed["chr"] = renamed["chr"].astype(str)
        renamed["start"] = renamed["start"].astype(int)
        renamed["end"] = renamed["end"].astype(int)
        renamed["position"] = renamed["start"].astype(int)
        frames.append(renamed[["feature", "type", "chr", "start", "end", "position"]])

    snp_path = FEATURE_DATA_DIR / dataset / f"{dataset}_snp_data.csv"
    snp_df = pd.read_csv(snp_path)
    snp_df = snp_df.rename(
        columns={"Feature": "feature", "#Chromosome": "chr", "Position": "position"}
    )
    snp_df["type"] = "SNP"
    snp_df["position"] = snp_df["position"].astype(int)
    snp_df["start"] = snp_df["position"]
    snp_df["end"] = snp_df["position"]
    frames.append(snp_df[["feature", "type", "chr", "start", "end", "position"]])

    merged = pd.concat(frames, ignore_index=True).drop_duplicates("feature")
    return merged


def load_gene_region(gff_path: Path, gene_id: str = "AT2G39730", flank_size: int = 2000) -> GeneRegion:
    with gff_path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            if f"ID=gene-{gene_id}" not in parts[8] and f"locus_tag={gene_id}" not in parts[8]:
                continue

            chrom = normalize_chrom(parts[0])
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            if strand == "-":
                upstream_start = end + 1
                upstream_end = end + flank_size
                downstream_start = max(1, start - flank_size)
                downstream_end = start - 1
            else:
                upstream_start = max(1, start - flank_size)
                upstream_end = start - 1
                downstream_start = end + 1
                downstream_end = end + flank_size
            return GeneRegion(
                chrom=chrom,
                strand=strand,
                gene_start=start,
                gene_end=end,
                upstream_start=upstream_start,
                upstream_end=upstream_end,
                downstream_start=downstream_start,
                downstream_end=downstream_end,
            )
    raise ValueError(f"Gene {gene_id} not found in {gff_path}")


def normalize_chrom(chrom: str) -> str:
    chrom = str(chrom).strip()
    if chrom.startswith("NC_"):
        mapping = {
            "NC_003070.9": "Chr1",
            "NC_003071.7": "Chr2",
            "NC_003074.8": "Chr3",
            "NC_003075.7": "Chr4",
            "NC_003076.8": "Chr5",
        }
        return mapping.get(chrom, chrom)
    if chrom.lower().startswith("chr"):
        suffix = chrom[3:]
        return f"Chr{suffix}"
    return chrom


def classify_region(chrom: str, position: int, region: GeneRegion) -> str | None:
    if normalize_chrom(chrom) != region.chrom:
        return None
    if region.upstream_start <= position <= region.upstream_end:
        return "upstream"
    if region.gene_start <= position <= region.gene_end:
        return "gene body"
    if region.downstream_start <= position <= region.downstream_end:
        return "downstream"
    return None


def build_membership_table() -> pd.DataFrame:
    region_def = load_gene_region(GFF_PATH)
    coords_by_dataset = {dataset: load_dataset_coords(dataset) for dataset in DATASETS}
    shap_by_dataset = {
        dataset: parse_feature_importances(BASE_DIR / f"{dataset}_import_results.txt")
        for dataset in DATASETS
    }

    feature_records: Dict[str, Dict[str, object]] = {}

    for dataset in DATASETS:
        coords = coords_by_dataset[dataset].set_index("feature", drop=False)
        for feature, shap_value in shap_by_dataset[dataset].items():
            if feature not in coords.index:
                continue
            row = coords.loc[feature]
            if str(row["type"]) not in PLOT_TYPES:
                continue
            region = classify_region(str(row["chr"]), int(row["position"]), region_def)
            if region is None:
                continue

            if feature not in feature_records:
                feature_type = str(row["type"])
                record = {
                    "feature": feature,
                    "type": feature_type,
                    "chr": normalize_chrom(str(row["chr"])),
                    "position": int(row["position"]),
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "region": region,
                }
                for ds in DATASETS:
                    record[f"{ds}_shap_value"] = 0.0
                    record[f"in_{ds}"] = 0
                for feature_type_name in TYPE_ORDER:
                    record[f"is_{feature_type_name}"] = int(feature_type == feature_type_name)
                feature_records[feature] = record

            feature_records[feature][f"{dataset}_shap_value"] = float(shap_value)
            feature_records[feature][f"in_{dataset}"] = 1

    memberships = pd.DataFrame(feature_records.values())
    if memberships.empty:
        raise ValueError("No SHAP-positive features found inside AT2G39730 upstream/body/downstream regions.")

    memberships["dataset_membership_count"] = memberships[[f"in_{ds}" for ds in DATASETS]].sum(axis=1)
    memberships["intersection_label"] = memberships.apply(
        lambda row: build_intersection_label(row, DATASET_ORDER),
        axis=1,
    )
    memberships = memberships.sort_values(
        ["region", "dataset_membership_count", "type", "position", "feature"],
        ascending=[True, False, True, True, True],
    ).reset_index(drop=True)
    return memberships


def build_intersection_label(row: pd.Series, order: Iterable[str]) -> str:
    labels = []
    for name in order:
        if name in DATASETS and int(row[f"in_{name}"]) == 1:
            labels.append(DATASET_LABELS[name])
        elif name in TYPE_ORDER and int(row[f"is_{name}"]) == 1:
            labels.append(TYPE_LABELS[name])
    return " | ".join(labels)


def summarize_region_type(df: pd.DataFrame, region_name: str, feature_type: str) -> pd.DataFrame:
    region_df = df[(df["region"] == region_name) & (df["type"] == feature_type)].copy()
    if region_df.empty:
        columns = ["region", "type", "intersection_label", "feature_count", "features"] + DATASET_ORDER
        return pd.DataFrame(columns=columns)

    group_cols = [f"in_{ds}" for ds in DATASETS]
    summary = region_df.groupby(group_cols, dropna=False).agg(
        feature_count=("feature", "size"),
        features=("feature", lambda values: ",".join(sorted(values))),
    ).reset_index()
    summary["region"] = region_name
    summary["type"] = feature_type
    summary["intersection_label"] = summary.apply(
        lambda row: build_summary_label(row, DATASET_ORDER),
        axis=1,
    )
    for dataset in DATASETS:
        summary[dataset] = summary[f"in_{dataset}"].astype(int)

    summary = summary.sort_values(
        ["feature_count", "intersection_label"],
        ascending=[False, True],
    ).reset_index(drop=True)
    keep_cols = ["region", "type", "intersection_label", "feature_count", "features"] + DATASET_ORDER
    return summary[keep_cols]


def build_summary_label(row: pd.Series, member_order: Iterable[str]) -> str:
    labels = []
    for dataset in member_order:
        value = row[f"in_{dataset}"] if f"in_{dataset}" in row.index else row[dataset]
        if int(value) == 1:
            labels.append(DATASET_LABELS[dataset])
    return " | ".join(labels)


def plot_region_upset(
    summary_df: pd.DataFrame,
    region_name: str,
    feature_type: str,
    output_prefix: Path,
    member_order: List[str],
) -> None:
    plt.rcParams["font.family"] = "DejaVu Sans"

    display_df = summary_df.copy()
    if display_df.empty:
        fig, ax = plt.subplots(figsize=(6, 5), facecolor="white")
        ax.axis("off")
        ax.text(0.5, 0.5, f"No {feature_type} top-{int(TOP_PERCENT * 100)}% features in {region_name}", ha="center", va="center", fontsize=18)
        fig.savefig(output_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
        fig.savefig(output_prefix.with_suffix(".pdf"), bbox_inches="tight")
        plt.close(fig)
        return

    counts = display_df["feature_count"].astype(int).tolist()
    set_sizes = {
        name: int((display_df[name].astype(int) * display_df["feature_count"].astype(int)).sum())
        for name in member_order
    }
    display_row_order = [
        name for name, size in sorted(set_sizes.items(), key=lambda item: (item[1], item[0])) if size > 0
    ]
    y_labels = [
        DATASET_LABELS[name] if name in DATASETS else TYPE_LABELS[name]
        for name in display_row_order
    ]
    y_spacing = 0.38
    y_positions = [idx * y_spacing for idx in range(len(y_labels))]

    base_fig_width = 11.0
    base_left_width = base_fig_width * (1.35 + 0.95) / (1.35 + 0.95 + 5.9)
    base_right_width = base_fig_width * 5.9 / (1.35 + 0.95 + 5.9)
    reference_intersections = 5
    fig_width = base_left_width + base_right_width * max(1, len(display_df)) / reference_intersections
    fig = plt.figure(figsize=(fig_width, 13.4), facecolor="white")
    grid = gridspec.GridSpec(
        nrows=2,
        ncols=3,
        width_ratios=[1.35, 0.95, 5.9],
        height_ratios=[2.0, 3.25],
        wspace=0.0,
        hspace=0.08,
    )

    ax_bar = fig.add_subplot(grid[0, 2])
    ax_matrix = fig.add_subplot(grid[1, 2], sharex=ax_bar)
    ax_setsize = fig.add_subplot(grid[1, 0], sharey=ax_matrix)
    ax_label = fig.add_subplot(grid[1, 1], sharey=ax_matrix)
    ax_blank_left = fig.add_subplot(grid[0, 0])
    ax_blank_mid = fig.add_subplot(grid[0, 1])
    ax_blank_left.axis("off")
    ax_blank_mid.axis("off")

    x_spacing = 1.55
    x = [idx * x_spacing for idx in range(len(display_df))]
    bar_color = "#111111"
    ax_bar.bar(x, counts, color=bar_color, edgecolor=bar_color, linewidth=0.4, width=0.90)
    ax_bar.set_ylabel("Intersection size", fontsize=16)
    ax_bar.grid(axis="y", color="#E5E5E5", linewidth=0.8)
    ax_bar.set_axisbelow(True)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)
    ax_bar.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_bar.tick_params(axis="y", labelsize=18)
    ax_bar.set_ylim(0, max(counts) * 1.12 if counts else 1)

    for x_pos, count in zip(x, counts):
        ax_bar.text(
            x_pos,
            count + max(counts) * 0.035,
            str(count),
            ha="center",
            va="bottom",
            fontsize=18,
            clip_on=False,
        )

    for idx, y in enumerate(y_positions):
        color = "#F6F6F6" if idx % 2 == 0 else "#FFFFFF"
        ax_matrix.axhspan(y - y_spacing * 0.45, y + y_spacing * 0.45, color=color, zorder=0)
        ax_setsize.axhspan(y - y_spacing * 0.45, y + y_spacing * 0.45, color=color, zorder=0)

    for idx, (_, row) in enumerate(display_df.iterrows()):
        active_y: List[int] = []
        x_pos = x[idx]
        for y, member_name in zip(y_positions, display_row_order):
            is_active = int(row[member_name]) == 1
            if is_active:
                active_y.append(y)
                ax_matrix.scatter(x_pos, y, s=300, color="#111111", zorder=3)
            else:
                ax_matrix.scatter(x_pos, y, s=185, color="#D0D0D0", zorder=2)
        if len(active_y) >= 2:
            ax_matrix.plot([x_pos, x_pos], [min(active_y), max(active_y)], color="#111111", linewidth=2.2, zorder=1)

    set_size_values = [set_sizes[name] for name in display_row_order]
    ax_setsize.barh(y_positions, set_size_values, height=0.15, color="#111111", edgecolor="#111111", linewidth=0.3)
    for y, value in zip(y_positions, set_size_values):
        ax_setsize.text(
            value + max(set_size_values) * 0.05,
            y,
            str(value),
            va="center",
            ha="right",
            fontsize=13,
        )

    ax_setsize.set_xlabel("")
    ax_setsize.invert_xaxis()
    ax_setsize.tick_params(axis="x", labelsize=18, pad=2)
    ax_setsize.tick_params(axis="y", left=False, labelleft=False)
    ax_setsize.spines["top"].set_visible(False)
    ax_setsize.spines["right"].set_visible(False)
    ax_setsize.spines["left"].set_visible(False)
    ax_setsize.set_xlim(max(set_size_values) * 1.38 if set_size_values else 1, 0)
    guide_x = 100 if set_size_values and max(set_size_values) >= 100 else (max(set_size_values) * 0.5 if set_size_values else 0)
    if guide_x:
        ax_setsize.set_xticks([guide_x, 0])
        ax_setsize.set_xticklabels([str(int(round(guide_x))), "0"])
    else:
        ax_setsize.set_xticks([0])
        ax_setsize.set_xticklabels(["0"])

    ax_label.set_xlim(0, 1)
    ax_label.set_xticks([])
    ax_label.set_yticks([])
    ax_label.spines["top"].set_visible(False)
    ax_label.spines["right"].set_visible(False)
    ax_label.spines["bottom"].set_visible(False)
    ax_label.spines["left"].set_visible(False)
    for y, label in zip(y_positions, y_labels):
        ax_label.text(0.84, y, label, ha="right", va="center", fontsize=16, color="#111111")

    ax_matrix.set_yticks(y_positions)
    ax_matrix.set_yticklabels([])
    ax_matrix.set_ylim((len(y_labels) - 1) * y_spacing + y_spacing * 0.75, -y_spacing * 0.75)
    ax_matrix.set_xlim(-x_spacing * 0.5, (len(display_df) - 1) * x_spacing + x_spacing * 0.5 if len(display_df) else x_spacing * 0.5)
    ax_matrix.set_xticks(x)
    ax_matrix.set_xticklabels([])
    ax_matrix.tick_params(axis="x", bottom=False)
    ax_matrix.tick_params(axis="y", length=0)
    ax_matrix.spines["top"].set_visible(False)
    ax_matrix.spines["right"].set_visible(False)
    ax_matrix.spines["bottom"].set_visible(False)
    ax_matrix.spines["left"].set_visible(False)

    fig.subplots_adjust(left=0.05, right=0.99, top=0.985, bottom=0.06)
    fig.savefig(output_prefix.with_suffix(".png"), dpi=300)
    fig.savefig(output_prefix.with_suffix(".pdf"))
    plt.close(fig)


def write_outputs(memberships: pd.DataFrame) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    membership_path = OUTPUT_DIR / "AT2G39730_region_feature_memberships.tsv"
    memberships.to_csv(membership_path, sep="\t", index=False)

    for region in REGION_ORDER:
        for feature_type in PLOT_TYPES:
            summary = summarize_region_type(memberships, region, feature_type)
            prefix = OUTPUT_DIR / f"AT2G39730_{REGION_FILENAME[region]}_{feature_type}_upset"
            summary.to_csv(prefix.with_name(prefix.name + "_summary.tsv"), sep="\t", index=False)
            plot_region_upset(summary, region, feature_type, prefix, DATASET_ORDER)


def main() -> None:
    memberships = build_membership_table()
    write_outputs(memberships)


if __name__ == "__main__":
    main()
