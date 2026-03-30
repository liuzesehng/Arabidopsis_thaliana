#!/usr/bin/env python3
"""AT2G39730 seven-set feature analysis and UpSet-style plotting."""

from __future__ import annotations

import argparse
import re
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd


TARGET_GENE = "AT2G39730"
FLANK_BP = 2000
ALPHA = 0.05
PREFIXES = ["A_TPM", "A_per_TPM", "B_TPM", "B_per_TPM", "total_TPM"]

ROOT_DIR = Path("/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng")
ARABIDOPSIS_DIR = ROOT_DIR / "Arabidopsis_thaliana"
TEMPERATURE_PLOT_DIR = ARABIDOPSIS_DIR / "list/xgboot/TPM_4.5_all/temperature_meth_plots"
FEATURE_DATA_DIR = ARABIDOPSIS_DIR / "list/xgboot/TPM_4.5_all/feature_data_extraction"
GFF_PATH = ROOT_DIR / "ref/Arabidopsis_thaliana/refgen/GCF_000001735.4_TAIR10.1_genomic.gff"
IMPORT_RESULT_DIR = ARABIDOPSIS_DIR / "list/xgboot/TPM_4.5_all"

SEQID_TO_CHR = {
    "NC_003070.9": "Chr1",
    "NC_003071.7": "Chr2",
    "NC_003074.8": "Chr3",
    "NC_003075.7": "Chr4",
    "NC_003076.8": "Chr5",
    "NC_037304.1": "ChrC",
    "NC_000932.1": "ChrC",
    "NC_000933.1": "ChrM",
}
TYPE_ORDER = ["CHH", "CHG", "CG", "SNP"]
REGION_ORDER = ["upstream", "gene body", "downstream"]
SET_ORDER = TYPE_ORDER + REGION_ORDER
TYPE_COLORS = {
    "CHH": "#4575b4",
    "CHG": "#53b27f",
    "CG": "#d8842a",
    "SNP": "#9b579c",
}
REGION_COLORS = {
    "upstream": "#d9ead3",
    "gene body": "#fce5cd",
    "downstream": "#cfe2f3",
}


@dataclass(frozen=True)
class GeneRegion:
    gene_id: str
    chromosome: str
    seqid: str
    start: int
    end: int
    strand: str

    def region_bounds(self, flank_bp: int) -> dict[str, tuple[int, int]]:
        if self.strand == "-":
            return {
                "upstream": (self.end + 1, self.end + flank_bp),
                "gene body": (self.start, self.end),
                "downstream": (self.start - flank_bp, self.start - 1),
            }
        return {
            "upstream": (self.start - flank_bp, self.start - 1),
            "gene body": (self.start, self.end),
            "downstream": (self.end + 1, self.end + flank_bp),
        }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot seven-set AT2G39730 significant feature intersections.")
    parser.add_argument(
        "--prefixes",
        nargs="+",
        choices=PREFIXES,
        default=PREFIXES,
        help="Dataset prefixes to process.",
    )
    parser.add_argument(
        "--flank-bp",
        type=int,
        default=FLANK_BP,
        help="Upstream/downstream flank size in base pairs.",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=ALPHA,
        help="BH-adjusted significance threshold.",
    )
    return parser.parse_args()


def normalize_feature_type(raw_value: object) -> str | None:
    if raw_value is None or pd.isna(raw_value):
        return None
    value = str(raw_value).strip().upper()
    if value == "SNP":
        return "SNP"
    if value in {"CHH", "CHG", "CG"}:
        return value
    return None


def read_target_gene(gff_path: Path, target_gene: str) -> GeneRegion:
    matches: list[GeneRegion] = []
    with gff_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = fields[8]
            if f"locus_tag={target_gene}" not in attrs:
                continue
            seqid, _source, _feature_type, start, end, _score, strand, _phase, _attrs = fields
            chromosome = SEQID_TO_CHR.get(seqid)
            if chromosome is None:
                raise KeyError(f"Unknown seqid to chromosome mapping for {seqid}")
            matches.append(
                GeneRegion(
                    gene_id=target_gene,
                    chromosome=chromosome,
                    seqid=seqid,
                    start=int(start),
                    end=int(end),
                    strand=strand,
                )
            )
    if not matches:
        raise ValueError(f"Target gene {target_gene} not found in {gff_path}")
    if len(matches) > 1:
        raise ValueError(f"Target gene {target_gene} matched multiple gene records in {gff_path}")
    return matches[0]


def load_selected_features_from_import_results(import_result_path: Path) -> pd.DataFrame:
    rows = []
    in_shap_section = False
    feature_line_pattern = re.compile(r"^\s*([A-Za-z0-9_]+):\s*([+-]?\d+(?:\.\d+)?)\s*$")

    with import_result_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if "所有特征及其平均|SHAP|值" in line:
                in_shap_section = True
                continue
            if not in_shap_section:
                continue
            if not line.strip():
                continue
            match = feature_line_pattern.match(line)
            if not match:
                continue
            feature = match.group(1)
            shap_value = float(match.group(2))
            feature_type = normalize_feature_type(feature.split("_", 1)[0])
            if feature_type is None or shap_value <= 0:
                continue
            rows.append({"feature": feature, "type": feature_type, "shap_value": shap_value})

    if not rows:
        return pd.DataFrame(columns=["feature", "type", "shap_value"])
    return pd.DataFrame(rows).drop_duplicates(subset=["feature", "type"], keep="first")


def load_all_selected_features(prefixes: list[str]) -> dict[str, pd.DataFrame]:
    selected_by_prefix: dict[str, pd.DataFrame] = {}
    for prefix in prefixes:
        import_result_path = IMPORT_RESULT_DIR / f"{prefix}_import_results.txt"
        selected_by_prefix[prefix] = load_selected_features_from_import_results(import_result_path)
    return selected_by_prefix


def intersect_selected_features(selected_by_prefix: dict[str, pd.DataFrame]) -> set[str]:
    feature_sets = []
    for selected_df in selected_by_prefix.values():
        feature_sets.append(set(selected_df["feature"].tolist()))
    if not feature_sets:
        return set()
    return set.intersection(*feature_sets)


def build_common_feature_shap_table(selected_by_prefix: dict[str, pd.DataFrame], common_features: set[str]) -> pd.DataFrame:
    merged_df: pd.DataFrame | None = None
    for prefix, selected_df in selected_by_prefix.items():
        subset = (
            selected_df[selected_df["feature"].isin(common_features)][["feature", "type", "shap_value"]]
            .rename(columns={"shap_value": f"{prefix}_shap_value"})
            .copy()
        )
        if merged_df is None:
            merged_df = subset
        else:
            merged_df = merged_df.merge(subset, on=["feature", "type"], how="inner", validate="one_to_one")
    if merged_df is None:
        return pd.DataFrame(columns=["feature", "type"])
    return merged_df


def load_feature_coordinates(prefix: str, feature_type: str) -> pd.DataFrame:
    feature_dir = FEATURE_DATA_DIR / prefix
    suffix = "snp" if feature_type == "SNP" else feature_type
    csv_path = feature_dir / f"{prefix}_{suffix}_data.csv"
    df = pd.read_csv(csv_path)

    if feature_type == "SNP":
        renamed = df.rename(
            columns={
                "Feature": "feature",
                "#Chromosome": "chr",
                "Position": "position",
            }
        ).copy()
        bed_start = pd.to_numeric(renamed["position"], errors="coerce")
        renamed["start"] = bed_start + 1
        renamed["end"] = renamed["start"]
    else:
        renamed = df.rename(
            columns={
                "Feature": "feature",
                "chr": "chr",
                "start": "start",
                "end": "end",
            }
        ).copy()
        bed_start = pd.to_numeric(renamed["start"], errors="coerce")
        bed_end = pd.to_numeric(renamed["end"], errors="coerce")
        renamed["start"] = bed_start + 1
        renamed["end"] = bed_end

    renamed["position"] = renamed["start"]
    renamed["type"] = feature_type
    renamed["chr"] = renamed["chr"].astype(str).str.replace("^chr", "Chr", regex=True)
    return renamed[["feature", "type", "chr", "position", "start", "end"]].dropna()


def assign_region(position: int, region_bounds: dict[str, tuple[int, int]]) -> str | None:
    for region_name in REGION_ORDER:
        region_start, region_end = region_bounds[region_name]
        if region_start <= position <= region_end:
            return region_name
    return None


def build_common_feature_table(
    gene_region: GeneRegion,
    flank_bp: int,
    common_feature_shap_df: pd.DataFrame,
) -> pd.DataFrame:
    if common_feature_shap_df.empty:
        return pd.DataFrame(
            columns=[
                "feature",
                "type",
                "chr",
                "position",
                "start",
                "end",
                "region",
                "membership_label",
            ]
        )

    coordinate_prefix = PREFIXES[0]
    coord_frames = [load_feature_coordinates(coordinate_prefix, feature_type) for feature_type in TYPE_ORDER]
    coords_df = pd.concat(coord_frames, ignore_index=True)
    merged = common_feature_shap_df.merge(coords_df, on=["feature", "type"], how="left", validate="one_to_one")
    if merged[["chr", "start", "end"]].isna().any(axis=None):
        missing = merged.loc[merged["chr"].isna(), "feature"].tolist()
        raise KeyError(f"Missing coordinates for shared features: {', '.join(missing[:10])}")

    merged = merged.loc[merged["chr"] == gene_region.chromosome].copy()
    region_bounds = gene_region.region_bounds(flank_bp)
    merged["position"] = pd.to_numeric(merged["position"], errors="coerce").astype(int)
    merged["region"] = merged["position"].map(lambda pos: assign_region(pos, region_bounds))
    merged = merged[merged["region"].notna()].copy()
    merged["membership_label"] = merged["region"] + " | " + merged["type"]
    shap_columns = [column for column in merged.columns if column.endswith("_shap_value")]
    merged["mean_shap_value"] = merged[shap_columns].mean(axis=1)
    merged = merged.sort_values(["region", "type", "mean_shap_value", "position", "feature"], ascending=[True, True, False, True, True]).reset_index(drop=True)
    return merged[
        ["feature", "type", *shap_columns, "mean_shap_value", "chr", "position", "start", "end", "region", "membership_label"]
    ]


def build_membership_summary(detail_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for region_name in REGION_ORDER:
        for feature_type in TYPE_ORDER:
            subset = detail_df[(detail_df["region"] == region_name) & (detail_df["type"] == feature_type)]
            feature_count = int(subset.shape[0])
            row = {
                "region": region_name,
                "type": feature_type,
                "intersection_label": f"{region_name} | {feature_type}",
                "feature_count": feature_count,
            }
            for set_name in SET_ORDER:
                row[set_name] = int(set_name in {region_name, feature_type})
            rows.append(row)

    summary_df = pd.DataFrame(rows)
    summary_df["region_rank"] = summary_df["region"].map({name: idx for idx, name in enumerate(REGION_ORDER)})
    summary_df["type_rank"] = summary_df["type"].map({name: idx for idx, name in enumerate(TYPE_ORDER)})
    summary_df = summary_df.sort_values(
        ["feature_count", "region_rank", "type_rank"],
        ascending=[False, True, True],
    ).reset_index(drop=True)
    return summary_df.drop(columns=["region_rank", "type_rank"])


def configure_matplotlib() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titleweight": "regular",
            "axes.labelsize": 14,
            "xtick.labelsize": 11,
            "ytick.labelsize": 12,
        }
    )


def plot_region_type_bars(summary_df: pd.DataFrame, output_prefix: Path) -> None:
    plot_df = (
        summary_df.pivot(index="region", columns="type", values="feature_count")
        .reindex(REGION_ORDER)
        .reindex(columns=TYPE_ORDER)
        .fillna(0)
        .astype(int)
    )

    fig, ax = plt.subplots(figsize=(9.5, 6.2))
    x_positions = list(range(len(REGION_ORDER)))
    width = 0.18
    offsets = [-1.5, -0.5, 0.5, 1.5]

    for offset, feature_type in zip(offsets, TYPE_ORDER):
        values = plot_df[feature_type].tolist()
        bar_positions = [x + offset * width for x in x_positions]
        bars = ax.bar(
            bar_positions,
            values,
            width=width,
            label=feature_type,
            color=TYPE_COLORS[feature_type],
            edgecolor="#333333",
            linewidth=0.9,
        )
        for bar, value in zip(bars, values):
            if value <= 0:
                continue
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                value + max(plot_df.to_numpy().max() * 0.01, 0.3),
                str(value),
                ha="center",
                va="bottom",
                fontsize=10,
            )

    ax.set_xticks(x_positions)
    ax.set_xticklabels(REGION_ORDER)
    ax.set_ylabel("Feature count")
    ax.grid(axis="y", linestyle="--", linewidth=0.7, alpha=0.35)
    ax.legend(frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.08))

    fig.subplots_adjust(left=0.1, right=0.98, bottom=0.12, top=0.88)
    fig.savefig(output_prefix.with_suffix(".png"), dpi=300)
    fig.savefig(output_prefix.with_suffix(".pdf"))
    plt.close(fig)


def write_outputs(
    detail_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    gene_region: GeneRegion,
    flank_bp: int,
) -> list[Path]:
    output_paths: list[Path] = []
    summary_path = TEMPERATURE_PLOT_DIR / f"{gene_region.gene_id}_shared_shap_region_type_summary.tsv"
    detail_path = TEMPERATURE_PLOT_DIR / f"{gene_region.gene_id}_shared_shap_region_type_sites.tsv"
    figure_prefix = TEMPERATURE_PLOT_DIR / f"{gene_region.gene_id}_shared_shap_region_type_bars"

    summary_df.to_csv(summary_path, sep="\t", index=False)
    detail_df.to_csv(detail_path, sep="\t", index=False)
    plot_region_type_bars(summary_df, figure_prefix)

    output_paths.extend([summary_path, detail_path, figure_prefix.with_suffix(".png"), figure_prefix.with_suffix(".pdf")])
    return output_paths


def validate_summary(summary_df: pd.DataFrame, detail_df: pd.DataFrame) -> None:
    active_counts = summary_df[SET_ORDER].sum(axis=1)
    if not (active_counts == 2).all():
        raise AssertionError("Each summary row must include exactly one feature-type set and one region set.")
    if int(summary_df["feature_count"].sum()) != int(detail_df.shape[0]):
        raise AssertionError("Summary counts do not sum to the number of detail rows.")


def process_common_features(
    gene_region: GeneRegion,
    flank_bp: int,
    alpha: float,
    common_feature_shap_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, list[Path]]:
    _ = alpha
    detail_df = build_common_feature_table(gene_region, flank_bp, common_feature_shap_df)
    if detail_df.empty:
        raise ValueError("No shared SHAP>0 features within target regions.")

    summary_df = build_membership_summary(detail_df)
    validate_summary(summary_df, detail_df)
    output_paths = write_outputs(detail_df, summary_df, gene_region, flank_bp)
    return detail_df, summary_df, output_paths


def main() -> None:
    args = parse_args()
    configure_matplotlib()
    gene_region = read_target_gene(GFF_PATH, TARGET_GENE)
    selected_by_prefix = load_all_selected_features(args.prefixes)
    common_features = intersect_selected_features(selected_by_prefix)
    if not common_features:
        raise ValueError("No shared SHAP>0 features found across the selected import_results datasets.")
    common_feature_shap_df = build_common_feature_shap_table(selected_by_prefix, common_features)
    detail_df, summary_df, generated_paths = process_common_features(
        gene_region,
        args.flank_bp,
        args.alpha,
        common_feature_shap_df,
    )

    print(f"Target gene: {gene_region.gene_id} ({gene_region.chromosome}:{gene_region.start}-{gene_region.end}, {gene_region.strand})")
    print(f"Flank size: {args.flank_bp} bp | Alpha: {args.alpha}")
    print(f"Shared SHAP>0 features across {len(args.prefixes)} datasets: {len(common_features)}")
    print(f"Shared SHAP>0 features within target regions: {len(detail_df)}")
    for path in generated_paths:
        print(f"Generated: {path}")


if __name__ == "__main__":
    main()
