#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch, FancyArrowPatch, Rectangle
import pandas as pd


PROJECT_ROOT = Path(
    "/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana"
)
DEFAULT_WORKDIR = PROJECT_ROOT / "list/RCA"
DEFAULT_OUTPUT_DIR = DEFAULT_WORKDIR / "AT2G39730_motif_gene_outputs"
DEFAULT_GFF = (
    Path("/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/ref/Arabidopsis_thaliana/refgen")
    / "GCF_000001735.4_TAIR10.1_genomic.gff"
)
DEFAULT_TF_POSITIVE = PROJECT_ROOT / "list/coex2/total/gene.TF_ids.txt"
DEFAULT_TF_NEGATIVE = PROJECT_ROOT / "list/coex2/all.fimo.out.filt"
DEFAULT_POSITIVE_ANNOT = PROJECT_ROOT / "list/coex2/total/total.gene.TF_ids.annotated.txt"
DEFAULT_NEGATIVE_ANNOT = PROJECT_ROOT / "list/coex2/total/total.motif.gene.TF_ids.annotated.txt"
DEFAULT_IMPORT_DIR = PROJECT_ROOT / "list/xgboot/TPM_4.5_all"
DEFAULT_FEATURE_DIR = DEFAULT_IMPORT_DIR / "feature_data_extraction"
DEFAULT_PREFIXES = ["total_TPM"]
DEFAULT_GENE_ID = "AT2G39730"
DEFAULT_GENE_WINDOW_BP = 2000
DEFAULT_FLANK_BP = 50
CONTEXT_COLORS = {"CG": "#6A3D9A", "CHG": "#1B9E77", "CHH": "#3182BD"}
BASE_COLORS = {"A": "#3182BD", "T": "#D95F0E", "C": "#31A354", "G": "#756BB1", "N": "#636363"}
GENE_BODY_COLOR = "#587565"
GENE_EXON_COLOR = "#496656"
GENE_CDS_COLOR = "#395548"
ARABIDOPSIS_CHR_MAP = {
    "NC_003070": "Chr1",
    "NC_003071": "Chr2",
    "NC_003074": "Chr3",
    "NC_003075": "Chr4",
    "NC_003076": "Chr5",
    "NC_001284": "ChrC",
    "NC_000932": "ChrM",
}


@dataclass
class GeneModel:
    gene_id: str
    gene_name: str
    chrom: str
    start: int
    end: int
    strand: str
    exons: list[tuple[int, int]]
    cds: list[tuple[int, int]]
    utrs: list[tuple[int, int]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot RCA gene structure with nearby TF motifs and SHAP-selected methylation sites."
    )
    parser.add_argument("--workdir", type=Path, default=DEFAULT_WORKDIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--gff", type=Path, default=DEFAULT_GFF)
    parser.add_argument("--gene-id", default=DEFAULT_GENE_ID)
    parser.add_argument("--prefixes", nargs="+", default=DEFAULT_PREFIXES)
    parser.add_argument("--gene-window-bp", type=int, default=DEFAULT_GENE_WINDOW_BP)
    parser.add_argument("--flank-bp", type=int, default=DEFAULT_FLANK_BP)
    parser.add_argument("--tf-positive", type=Path, default=DEFAULT_TF_POSITIVE)
    parser.add_argument("--tf-negative", type=Path, default=DEFAULT_TF_NEGATIVE)
    parser.add_argument("--positive-annot", type=Path, default=DEFAULT_POSITIVE_ANNOT)
    parser.add_argument("--negative-annot", type=Path, default=DEFAULT_NEGATIVE_ANNOT)
    parser.add_argument("--import-dir", type=Path, default=DEFAULT_IMPORT_DIR)
    parser.add_argument("--feature-dir", type=Path, default=DEFAULT_FEATURE_DIR)
    return parser.parse_args()


def ensure_paths_exist(paths: list[Path]) -> None:
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required files:\n" + "\n".join(missing))


def parse_attributes(attr_text: str) -> dict[str, str]:
    attributes = {}
    for item in attr_text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attributes[key] = value
    return attributes


def normalize_chromosome(raw_value: str) -> str:
    text = str(raw_value).strip()
    if not text:
        return text

    base = text.split(".", 1)[0]
    if base in ARABIDOPSIS_CHR_MAP:
        return ARABIDOPSIS_CHR_MAP[base]

    lowered = text.lower().replace("chromosome", "chr")
    lowered = lowered.replace("_", "").replace("-", "")
    if lowered.startswith("chr"):
        suffix = lowered[3:]
        if suffix in {"1", "2", "3", "4", "5"}:
            return f"Chr{suffix}"
        if suffix in {"c", "chloro", "chloroplast"}:
            return "ChrC"
        if suffix in {"m", "mito", "mitochondria"}:
            return "ChrM"
    return text


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    ordered = sorted((int(start), int(end)) for start, end in intervals)
    merged = [ordered[0]]
    for start, end in ordered[1:]:
        prev_start, prev_end = merged[-1]
        if start <= prev_end + 1:
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))
    return merged


def load_gene_model(gff_path: Path, gene_id: str) -> GeneModel:
    gene_record: dict[str, object] | None = None
    exons: list[tuple[int, int]] = []
    cds: list[tuple[int, int]] = []
    utrs: list[tuple[int, int]] = []

    with gff_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, _source, feature_type, start, end, _score, strand, _phase, attr_text = parts
            attrs = parse_attributes(attr_text)
            if attrs.get("locus_tag") != gene_id:
                continue

            start_i = int(start)
            end_i = int(end)
            if feature_type == "gene":
                gene_record = {
                    "gene_name": attrs.get("gene", gene_id),
                    "chrom": normalize_chromosome(seqid),
                    "start": start_i,
                    "end": end_i,
                    "strand": strand,
                }
            elif feature_type == "exon":
                exons.append((start_i, end_i))
            elif feature_type == "CDS":
                cds.append((start_i, end_i))
            elif "UTR" in feature_type.upper():
                utrs.append((start_i, end_i))

    if gene_record is None:
        raise ValueError(f"Gene {gene_id} not found in {gff_path}")

    return GeneModel(
        gene_id=gene_id,
        gene_name=str(gene_record["gene_name"]),
        chrom=str(gene_record["chrom"]),
        start=int(gene_record["start"]),
        end=int(gene_record["end"]),
        strand=str(gene_record["strand"]),
        exons=merge_intervals(exons),
        cds=merge_intervals(cds),
        utrs=merge_intervals(utrs),
    )


def parse_import_results(path: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    in_shap_section = False
    pattern = re.compile(r"^\s*([A-Za-z0-9_]+):\s*([-+0-9.eE]+)\s*$")

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if "所有特征及其平均|SHAP|值" in line:
                in_shap_section = True
                continue
            if not in_shap_section:
                continue
            if not line.strip():
                continue

            match = pattern.match(line)
            if not match:
                if rows:
                    break
                continue

            feature = match.group(1)
            shap_value = float(match.group(2))
            if shap_value <= 0:
                continue
            context = feature.split("_", 1)[0].upper()
            rows.append({"feature": feature, "context": context, "shap_value": shap_value})

    result = pd.DataFrame(rows).drop_duplicates(subset=["feature"], keep="first")
    if result.empty:
        return result
    result = result.sort_values(["shap_value", "context", "feature"], ascending=[False, True, True])
    top_n = max(1, math.ceil(len(result) * 0.10))
    top5_n = max(1, math.ceil(len(result) * 0.05))
    top_features = set(result.head(top_n)["feature"].astype(str))
    top5_features = set(result.head(top5_n)["feature"].astype(str))
    result = result[result["feature"].astype(str).isin(top_features) & result["context"].isin(CONTEXT_COLORS)].copy()
    result["is_top10_percent"] = result["feature"].astype(str).isin(top_features)
    result["is_top5_percent"] = result["feature"].astype(str).isin(top5_features)
    return result.reset_index(drop=True)


def read_motif_table(path: Path) -> pd.DataFrame:
    expected_columns = [
        "motif_id",
        "motif_alt_id",
        "sequence_name",
        "start",
        "stop",
        "strand",
        "score",
        "p-value",
        "q-value",
        "matched_sequence",
    ]
    with path.open("r", encoding="utf-8") as handle:
        first_line = handle.readline().strip()
    if first_line.startswith("motif_id\t"):
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path, sep="\t", names=expected_columns)


def build_symbol_map(*annot_paths: Path) -> dict[str, str]:
    symbol_map: dict[str, str] = {}
    for annot_path in annot_paths:
        if not annot_path.exists():
            continue
        annot_df = pd.read_csv(annot_path, sep="\t")
        if "Gene_ID" not in annot_df.columns or "Symbol" not in annot_df.columns:
            continue
        for gene_id, symbol in zip(annot_df["Gene_ID"].astype(str), annot_df["Symbol"].fillna("").astype(str)):
            if gene_id not in symbol_map and symbol.strip():
                symbol_map[gene_id] = symbol
    return symbol_map


def motif_occurrence_key(df: pd.DataFrame) -> pd.Series:
    return (
        df["motif_id"].astype(str)
        + "|"
        + df["motif_alt_id"].astype(str)
        + "|"
        + df["sequence_name"].astype(str)
        + "|"
        + df["start"].astype(str)
        + "|"
        + df["stop"].astype(str)
        + "|"
        + df["strand"].astype(str)
        + "|"
        + df["matched_sequence"].astype(str).str.upper()
    )


def merge_adjacent_same_motifs(motif_df: pd.DataFrame) -> pd.DataFrame:
    if motif_df.empty:
        return motif_df.copy()
    sort_cols = ["motif_id", "sequence_name", "start", "stop", "strand"]
    ordered = motif_df.sort_values(sort_cols).reset_index(drop=True)
    merged_rows: list[dict[str, object]] = []
    current: dict[str, object] | None = None

    for row in ordered.to_dict("records"):
        row["matched_sequence"] = str(row["matched_sequence"]).upper()
        row["start"] = int(row["start"])
        row["stop"] = int(row["stop"])
        row["score"] = float(row["score"])
        row["p-value"] = float(row["p-value"])
        row["q-value"] = float(row["q-value"])
        if current is None:
            current = dict(row)
            continue

        same_group = (
            row["motif_id"] == current["motif_id"]
            and row["sequence_name"] == current["sequence_name"]
            and row["start"] <= int(current["stop"]) + 5
        )
        if same_group:
            current["start"] = min(int(current["start"]), row["start"])
            current["stop"] = max(int(current["stop"]), row["stop"])
            if row["score"] > float(current["score"]):
                current["score"] = row["score"]
            current["p-value"] = min(float(current["p-value"]), row["p-value"])
            current["q-value"] = min(float(current["q-value"]), row["q-value"])
            current["strand"] = current["strand"] if row["strand"] == current["strand"] else "."
            if len(row["matched_sequence"]) > len(str(current["matched_sequence"])):
                current["matched_sequence"] = row["matched_sequence"]
        else:
            merged_rows.append(current)
            current = dict(row)

    if current is not None:
        merged_rows.append(current)
    return pd.DataFrame(merged_rows)


def load_feature_positions(prefix: str, feature_dir: Path) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for context in ["CHH", "CHG", "CG"]:
        csv_path = feature_dir / prefix / f"{prefix}_{context}_data.csv"
        df = pd.read_csv(csv_path)
        frame = pd.DataFrame(
            {
                "feature": df["Feature"].astype(str),
                "context": context,
                "chrom": df["chr"].map(normalize_chromosome),
                "bed_start": pd.to_numeric(df["start"], errors="coerce"),
                "bed_end": pd.to_numeric(df["end"], errors="coerce"),
            }
        )
        frame["plot_start"] = frame["bed_start"] + 1
        frame["plot_end"] = frame["bed_end"]
        frame["position"] = ((frame["plot_start"] + frame["plot_end"]) / 2.0).astype(float)
        frames.append(frame)
    feature_df = pd.concat(frames, ignore_index=True)
    return feature_df.dropna(subset=["plot_start", "plot_end", "position"]).reset_index(drop=True)


def load_motif_hits(path: Path, symbol_map: dict[str, str], set_type: str) -> pd.DataFrame:
    motif_df = read_motif_table(path)

    motif_df = motif_df.copy()
    motif_df["chrom"] = motif_df["sequence_name"].map(normalize_chromosome)
    motif_df["motif_start"] = pd.to_numeric(motif_df["start"], errors="coerce")
    motif_df["motif_end"] = pd.to_numeric(motif_df["stop"], errors="coerce")
    motif_df["strand"] = motif_df["strand"].fillna(".").astype(str)
    motif_df["matched_sequence"] = motif_df["matched_sequence"].fillna("").astype(str).str.upper()
    motif_df["symbol"] = motif_df["motif_id"].astype(str).map(symbol_map).fillna("")
    motif_df["symbol"] = motif_df.apply(
        lambda row: row["symbol"] if str(row["symbol"]).strip() else row["motif_id"], axis=1
    )
    motif_df["set_type"] = set_type
    motif_df["motif_center"] = (
        pd.to_numeric(motif_df["motif_start"], errors="coerce") + pd.to_numeric(motif_df["motif_end"], errors="coerce")
    ) / 2.0
    return motif_df.dropna(subset=["motif_start", "motif_end"]).reset_index(drop=True)


def load_rca_motifs(args: argparse.Namespace, gene_model: GeneModel) -> pd.DataFrame:
    window_start = gene_model.start - args.gene_window_bp
    window_end = gene_model.end + args.gene_window_bp

    symbol_map = build_symbol_map(args.positive_annot, args.negative_annot)

    positive_raw = read_motif_table(args.tf_positive)
    positive_df = load_motif_hits(args.tf_positive, symbol_map, "positive")

    negative_raw = read_motif_table(args.tf_negative)
    negative_raw = negative_raw.copy()
    positive_keys = set(motif_occurrence_key(positive_raw))
    negative_raw = negative_raw.loc[~motif_occurrence_key(negative_raw).isin(positive_keys)].copy()
    negative_raw = merge_adjacent_same_motifs(negative_raw)
    negative_tmp_path = args.workdir / "_merged_negative_tmp.tsv"
    negative_raw.to_csv(negative_tmp_path, sep="\t", index=False)
    negative_df = load_motif_hits(negative_tmp_path, symbol_map, "negative")
    if negative_tmp_path.exists():
        negative_tmp_path.unlink()

    motif_df = pd.concat([positive_df, negative_df], ignore_index=True)
    motif_df = motif_df[
        (motif_df["chrom"] == gene_model.chrom)
        & (motif_df["motif_end"] >= window_start)
        & (motif_df["motif_start"] <= window_end)
    ].copy()
    motif_df["motif_start"] = motif_df["motif_start"].astype(int)
    motif_df["motif_end"] = motif_df["motif_end"].astype(int)
    motif_df = motif_df.sort_values(["motif_center", "motif_id"]).reset_index(drop=True)
    return motif_df


def load_all_fimo_motifs(args: argparse.Namespace, gene_model: GeneModel) -> pd.DataFrame:
    window_start = gene_model.start - args.gene_window_bp
    window_end = gene_model.end + args.gene_window_bp
    symbol_map = build_symbol_map(args.positive_annot, args.negative_annot)
    all_fimo_raw = read_motif_table(args.tf_negative).copy()
    all_fimo_raw = merge_adjacent_same_motifs(all_fimo_raw)
    all_fimo_tmp_path = args.workdir / "_all_fimo_merged_tmp.tsv"
    all_fimo_raw.to_csv(all_fimo_tmp_path, sep="\t", index=False)
    motif_df = load_motif_hits(all_fimo_tmp_path, symbol_map, "all_fimo")
    if all_fimo_tmp_path.exists():
        all_fimo_tmp_path.unlink()
    motif_df = motif_df[
        (motif_df["chrom"] == gene_model.chrom)
        & (motif_df["motif_end"] >= window_start)
        & (motif_df["motif_start"] <= window_end)
    ].copy()
    motif_df["motif_start"] = motif_df["motif_start"].astype(int)
    motif_df["motif_end"] = motif_df["motif_end"].astype(int)
    motif_df = motif_df.sort_values(["motif_center", "motif_id"]).reset_index(drop=True)
    return motif_df


def build_motif_feature_matches(
    motif_df: pd.DataFrame,
    feature_df: pd.DataFrame,
    shap_df: pd.DataFrame,
    flank_bp: int,
    dataset: str,
) -> pd.DataFrame:
    merged_features = feature_df.merge(shap_df, on=["feature", "context"], how="inner", validate="one_to_one")
    records: list[dict[str, object]] = []

    for motif in motif_df.itertuples(index=False):
        window_start = int(motif.motif_start) - flank_bp
        window_end = int(motif.motif_end) + flank_bp
        subset = merged_features[
            (merged_features["chrom"] == motif.chrom)
            & (merged_features["position"] >= window_start)
            & (merged_features["position"] <= window_end)
        ].copy()

        if subset.empty:
            records.append(
                {
                    "dataset": dataset,
                    "set_type": motif.set_type,
                    "motif_id": motif.motif_id,
                    "symbol": motif.symbol,
                    "strand": motif.strand,
                    "chrom": motif.chrom,
                    "motif_start": int(motif.motif_start),
                    "motif_end": int(motif.motif_end),
                    "matched_sequence": motif.matched_sequence,
                    "window_start": window_start,
                    "window_end": window_end,
                    "methylation_feature": pd.NA,
                    "context": pd.NA,
                    "feature_pos": pd.NA,
                    "shap_value": pd.NA,
                    "is_top10_percent": pd.NA,
                    "is_top5_percent": pd.NA,
                    "in_window": False,
                    "in_motif": False,
                }
            )
            continue

        for row in subset.itertuples(index=False):
            position = float(row.position)
            records.append(
                {
                    "dataset": dataset,
                    "set_type": motif.set_type,
                    "motif_id": motif.motif_id,
                    "symbol": motif.symbol,
                    "strand": motif.strand,
                    "chrom": motif.chrom,
                    "motif_start": int(motif.motif_start),
                    "motif_end": int(motif.motif_end),
                    "matched_sequence": motif.matched_sequence,
                    "window_start": window_start,
                    "window_end": window_end,
                    "methylation_feature": row.feature,
                    "context": row.context,
                    "feature_pos": position,
                    "shap_value": float(row.shap_value),
                    "is_top10_percent": bool(row.is_top10_percent),
                    "is_top5_percent": bool(row.is_top5_percent),
                    "in_window": True,
                    "in_motif": int(motif.motif_start) <= position <= int(motif.motif_end),
                }
            )

    return pd.DataFrame(records)


def assign_box_layout(motif_df: pd.DataFrame) -> pd.DataFrame:
    layout_df = motif_df.copy()
    layout_df["box_width"] = layout_df["matched_sequence"].str.len().clip(lower=10).apply(
        lambda n: min(0.58, max(0.40, 0.34 + 0.0125 * float(n)))
    )
    layout_df["box_height"] = 0.18

    top_ys = [0.83, 0.64]
    bottom_ys = [0.07, 0.26]
    lane_state = {("top", idx): -1.0 for idx in range(len(top_ys))}
    lane_state.update({("bottom", idx): -1.0 for idx in range(len(bottom_ys))})
    margins = (0.03, 0.97)

    box_lefts: list[float] = []
    box_bottoms: list[float] = []
    sides: list[str] = []
    lanes: list[int] = []

    for order_idx, motif in enumerate(layout_df.sort_values(["motif_center", "motif_id"]).itertuples(index=False)):
        total = max(1, len(layout_df))
        side_idx = order_idx // 2
        edge_steps = [0.10, 0.90, 0.24, 0.76, 0.38, 0.62, 0.50]
        desired_center = edge_steps[min(side_idx, len(edge_steps) - 1)]
        width = float(motif.box_width)
        preferred_side = "top" if order_idx % 2 == 0 else "bottom"
        candidate_specs = []
        for side, y_positions in [("top", top_ys), ("bottom", bottom_ys)]:
            for lane_idx, lane_y in enumerate(y_positions):
                prev_right = lane_state[(side, lane_idx)]
                left = max(desired_center - width / 2.0, prev_right + 0.02)
                left = max(margins[0], min(left, margins[1] - width))
                shift = abs((left + width / 2.0) - desired_center)
                penalty = 0.0 if side == preferred_side else 0.05
                candidate_specs.append((shift + penalty, side, lane_idx, lane_y, left))
        candidate_specs.sort(key=lambda item: item[0])
        _, side, lane_idx, lane_y, left = candidate_specs[0]
        lane_state[(side, lane_idx)] = left + width
        box_lefts.append(left)
        box_bottoms.append(lane_y)
        sides.append(side)
        lanes.append(lane_idx)

    ordered = layout_df.sort_values(["motif_center", "motif_id"]).copy()
    ordered["box_left"] = box_lefts
    ordered["box_bottom"] = box_bottoms
    ordered["box_side"] = sides
    ordered["box_lane"] = lanes
    ordered["box_center"] = ordered["box_left"] + ordered["box_width"] / 2.0
    ordered.loc[(ordered["box_side"] == "top") & (ordered["symbol"] == "LHY"), "box_bottom"] = 0.62
    ordered.loc[ordered["symbol"] == "RVE1", "box_bottom"] = 0.48
    ordered.loc[ordered["symbol"] == "AGL25", "box_bottom"] = 0.56
    ordered.loc[ordered["symbol"] == "AGL25", "box_left"] = 0.02
    if set(ordered["set_type"].astype(str)) == {"positive"}:
        ordered["box_width"] = ordered["box_width"].clip(lower=0.205, upper=0.225)
        manual_positions = {
            "AGL22": (0.60, 0.80),
            "LHY": (0.84, 0.78),
            "RVE1": (0.69, 0.47),
            "LCL5": (0.60, 0.13),
            "LCL1": (0.84, 0.13),
        }
        for symbol, (left, bottom) in manual_positions.items():
            ordered.loc[ordered["symbol"] == symbol, "box_left"] = left
            ordered.loc[ordered["symbol"] == symbol, "box_bottom"] = bottom
    ordered["box_center"] = ordered["box_left"] + ordered["box_width"] / 2.0
    return ordered.sort_index()


def spread_anchor_positions(values: list[float], min_gap: float) -> list[float]:
    if not values:
        return []
    placed = [values[0]]
    for value in values[1:]:
        placed.append(max(value, placed[-1] + min_gap))
    midpoint_shift = (placed[-1] - values[-1]) / 2.0 if len(values) > 1 else 0.0
    return [item - midpoint_shift for item in placed]


def to_fig_coords(fig: plt.Figure, ax: plt.Axes, xy: tuple[float, float], coords: str) -> tuple[float, float]:
    if coords == "data":
        display_xy = ax.transData.transform(xy)
    elif coords == "axes":
        display_xy = ax.transAxes.transform(xy)
    else:
        raise ValueError(f"Unsupported coords: {coords}")
    fig_xy = fig.transFigure.inverted().transform(display_xy)
    return float(fig_xy[0]), float(fig_xy[1])


def add_elbow_connector(
    fig: plt.Figure,
    main_ax: plt.Axes,
    inset_ax: plt.Axes,
    start_x: float,
    start_y: float,
    lane_y: float,
    corridor_x: float,
    end_anchor: tuple[float, float],
    color: str = "#4D4D4D",
    linewidth: float = 1.2,
) -> None:
    start_fig = to_fig_coords(fig, main_ax, (start_x, start_y), "data")
    end_fig = to_fig_coords(fig, inset_ax, end_anchor, "axes")
    lane_fig = to_fig_coords(fig, main_ax, (start_x, lane_y), "data")[1]

    segments = [
        ((start_fig[0], start_fig[1]), (start_fig[0], lane_fig)),
        ((start_fig[0], lane_fig), (corridor_x, lane_fig)),
        ((corridor_x, lane_fig), (corridor_x, end_fig[1])),
    ]
    for (x0, y0), (x1, y1) in segments:
        if abs(x0 - x1) < 1e-6 and abs(y0 - y1) < 1e-6:
            continue
        fig.add_artist(
            Line2D(
                [x0, x1],
                [y0, y1],
                transform=fig.transFigure,
                color=color,
                linewidth=linewidth,
                zorder=2,
            )
        )

    fig.add_artist(
        FancyArrowPatch(
            posA=(corridor_x, end_fig[1]),
            posB=end_fig,
            transform=fig.transFigure,
            arrowstyle="->",
            mutation_scale=12,
            linewidth=linewidth,
            color=color,
            zorder=2,
            shrinkA=0.0,
            shrinkB=0.0,
        )
    )


def mix_with_white(base_color: str, strength: float) -> tuple[float, float, float]:
    base_rgb = to_rgb(base_color)
    strength = max(0.22, min(1.0, strength))
    return tuple(1.0 - (1.0 - channel) * strength for channel in base_rgb)


def compute_display_positions(motif_row: pd.Series, valid_points: pd.DataFrame) -> dict[float, float]:
    motif_start = int(motif_row["motif_start"])
    motif_end = int(motif_row["motif_end"])
    flank_start = motif_start - int(motif_row["flank_bp"])
    flank_end = motif_end + int(motif_row["flank_bp"])
    anchors: set[float] = {float(flank_start), float(motif_start), float(motif_end), float(flank_end)}
    anchors.update(float(pos) for pos in range(motif_start, motif_end + 1))
    anchors.update(float(pos) for pos in valid_points["feature_pos"].dropna().tolist())
    sorted_anchors = sorted(anchors)
    display_map = {sorted_anchors[0]: 0.0}

    for prev, curr in zip(sorted_anchors, sorted_anchors[1:]):
        gap = curr - prev
        prev_in = motif_start <= prev <= motif_end
        curr_in = motif_start <= curr <= motif_end
        if prev_in and curr_in:
            increment = 1.9
        elif gap <= 1:
            increment = 1.45
        elif gap <= 5:
            increment = 1.95
        else:
            increment = 2.4 + min(2.1, math.log1p(gap) * 0.5)
        display_map[curr] = display_map[prev] + increment
    return display_map


def draw_motif_inset(ax: plt.Axes, motif_row: pd.Series, match_df: pd.DataFrame, flank_bp: int) -> None:
    motif_start = int(motif_row["motif_start"])
    motif_end = int(motif_row["motif_end"])
    is_negative = str(motif_row["strand"]) == "-"
    valid_points = match_df[
        (match_df["motif_id"] == motif_row["motif_id"]) &
        (match_df["set_type"] == motif_row["set_type"]) &
        (match_df["in_window"] == True)
    ].copy()
    display_map = compute_display_positions(
        pd.Series(
            {
                "motif_start": motif_start,
                "motif_end": motif_end,
                "flank_bp": flank_bp,
            }
        ),
        valid_points,
    )
    flank_start = float(motif_start - flank_bp)
    flank_end = float(motif_end + flank_bp)
    x_min = display_map[flank_start]
    x_max = display_map[flank_end]

    ax.set_xlim(x_min, x_max)
    if is_negative:
        ax.invert_xaxis()
    ax.set_ylim(-1.25, 1.35)
    ax.axhline(0.0, color="#444444", linewidth=1.0, zorder=1)
    ax.axvspan(
        display_map[float(motif_start)] - 0.3,
        display_map[float(motif_end)] + 0.3,
        color="#F2F0F7",
        alpha=0.9,
        zorder=0,
    )
    ax.axvline(display_map[float(motif_start)], color="#969696", linestyle="--", linewidth=0.8)
    ax.axvline(display_map[float(motif_end)], color="#969696", linestyle="--", linewidth=0.8)

    seq = str(motif_row["matched_sequence"]).upper()
    for offset, base in enumerate(seq):
        genomic_pos = motif_end - offset if is_negative else motif_start + offset
        xpos = display_map[float(genomic_pos)]
        ax.text(
            xpos,
            0.55 + (0.10 if offset % 2 else 0.0),
            base,
            ha="center",
            va="center",
            fontsize=16,
            fontweight="heavy",
            color=BASE_COLORS.get(base, BASE_COLORS["N"]),
            zorder=3,
        )

    if not valid_points.empty:
        placed_labels: list[tuple[float, int]] = []
        label_levels = [-0.10, -0.19, -0.28, -0.37]
        for row in valid_points.sort_values(["feature_pos", "context", "methylation_feature"]).itertuples(index=False):
            xpos = display_map[float(row.feature_pos)]
            context = str(row.context)
            color = mix_with_white(CONTEXT_COLORS[context], 1.0 if bool(row.is_top5_percent) else 0.34)
            size = 220 if bool(row.in_motif) else 170
            ax.scatter(
                [xpos],
                [-0.52],
                marker="*",
                s=size,
                color=color,
                alpha=1.0,
                edgecolors="none",
                zorder=4,
            )

            level_idx = 0
            while level_idx < len(label_levels):
                conflict = False
                for prev_x, prev_level in placed_labels:
                    if prev_level == level_idx and abs(prev_x - xpos) < 2.8:
                        conflict = True
                        break
                if not conflict:
                    break
                level_idx += 1
            if level_idx >= len(label_levels):
                level_idx = len(label_levels) - 1
            placed_labels.append((xpos, level_idx))

            ax.text(
                xpos,
                label_levels[level_idx],
                str(row.methylation_feature).split("_", 1)[-1],
                ha="center",
                va="top",
                fontsize=6.1,
                fontweight="bold",
                color="#333333",
                rotation=0,
                zorder=5,
            )

    ax.text(
        0.02,
        0.96,
        f"{motif_row['symbol']}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        fontweight="bold",
        color="#222222",
    )
    ax.text(
        0.98,
        0.96,
        (
            f"(-) {motif_row['chrom']}:{motif_end}-{motif_start}"
            if is_negative
            else f"(+) {motif_row['chrom']}:{motif_start}-{motif_end}"
        ),
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        fontweight="bold",
        color="#666666",
    )
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_edgecolor("#4D4D4D")
        spine.set_linewidth(0.9)


def draw_gene_model(
    ax: plt.Axes,
    gene_model: GeneModel,
    window_start: int,
    window_end: int,
    coord_transform,
) -> None:
    ax.axhline(0.0, color="black", linewidth=1.2, zorder=1)
    for start, end in gene_model.utrs:
        x0 = coord_transform(start)
        x1 = coord_transform(end)
        ax.add_patch(
            Rectangle((x0, -0.018), x1 - x0, 0.036, facecolor=GENE_BODY_COLOR, edgecolor="#2E463B", linewidth=0.8)
        )
    for start, end in gene_model.exons:
        x0 = coord_transform(start)
        x1 = coord_transform(end)
        ax.add_patch(
            Rectangle((x0, -0.032), x1 - x0, 0.064, facecolor=GENE_EXON_COLOR, edgecolor="#2E463B", linewidth=0.9)
        )
    for start, end in gene_model.cds:
        x0 = coord_transform(start)
        x1 = coord_transform(end)
        ax.add_patch(
            Rectangle((x0, -0.038), x1 - x0, 0.076, facecolor=GENE_CDS_COLOR, edgecolor="#2E463B", linewidth=0.9)
        )

    label_x = coord_transform(gene_model.start) + (coord_transform(gene_model.end) - coord_transform(gene_model.start)) * 0.04
    ax.text(
        window_end,
        -0.09,
        f"Chr2:{gene_model.end}-{gene_model.start}",
        ha="right",
        va="center",
        fontsize=13,
        fontweight="bold",
        color="#222222",
    )


def add_context_legend(ax: plt.Axes) -> None:
    x0 = 0.68
    y0 = 0.95
    labels = [("CG", "CG methylation"), ("CHG", "CHG methylation"), ("CHH", "CHH methylation")]
    for idx, (context, label) in enumerate(labels):
        color = CONTEXT_COLORS[context]
        ax.scatter(
            [x0 + idx * 0.07],
            [y0],
            transform=ax.transAxes,
            marker="*",
            s=220,
            color=color,
            clip_on=False,
            zorder=10,
        )
        ax.text(
            x0 + idx * 0.07 + 0.018,
            y0,
            context,
            transform=ax.transAxes,
            va="center",
            ha="left",
            fontsize=14,
            fontweight="bold",
        )


def sanitize_label(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(text)).strip("_") or "motif"


def plot_gene_overview(
    gene_model: GeneModel,
    motif_df: pd.DataFrame,
    prefix: str,
    set_type: str,
    output_dir: Path,
    gene_window_bp: int,
) -> None:
    if motif_df.empty:
        return
    output_dir.mkdir(parents=True, exist_ok=True)
    left_extra_pad = 400
    right_extra_pad = 400
    window_start = gene_model.start - gene_window_bp - left_extra_pad
    window_end = gene_model.end + gene_window_bp + right_extra_pad

    fig = plt.figure(figsize=(24, 8.0))
    main_ax = fig.add_axes([0.04, 0.16, 0.93, 0.74])
    main_ax.set_xlim(window_start, window_end)
    main_ax.set_ylim(-0.66, 0.66)
    span = window_end - window_start

    gene_mid = (gene_model.start + gene_model.end) / 2.0
    def keep_x(value: float) -> float:
        return gene_mid + (float(value) - gene_mid) * 0.72

    draw_gene_model(main_ax, gene_model, window_start, window_end, keep_x)
    add_context_legend(main_ax)

    main_ax.set_yticks([])
    main_ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    main_ax.spines["left"].set_visible(False)
    main_ax.spines["right"].set_visible(False)
    main_ax.spines["top"].set_visible(False)
    main_ax.spines["bottom"].set_visible(False)
    main_ax.tick_params(axis="x", labelsize=11, pad=1, length=0)

    fig.savefig(output_dir / f"{prefix}_{set_type}_AT2G39730_gene_overview.png", dpi=300, bbox_inches="tight")
    fig.savefig(output_dir / f"{prefix}_{set_type}_AT2G39730_gene_overview.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_single_motif_figure(
    motif_row: pd.Series,
    match_df: pd.DataFrame,
    prefix: str,
    set_type: str,
    output_dir: Path,
    flank_bp: int,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(12, 3.8))
    ax = fig.add_axes([0.06, 0.18, 0.88, 0.70])
    draw_motif_inset(ax, motif_row, match_df, flank_bp)
    stem = sanitize_label(str(motif_row["symbol"]))
    fig.savefig(output_dir / f"{prefix}_{set_type}_{stem}_motif.png", dpi=300, bbox_inches="tight")
    fig.savefig(output_dir / f"{prefix}_{set_type}_{stem}_motif.pdf", bbox_inches="tight")
    plt.close(fig)


def motif_has_methylation(motif_row: pd.Series, match_df: pd.DataFrame) -> bool:
    subset = match_df[
        (match_df["motif_id"] == motif_row["motif_id"])
        & (match_df["motif_start"] == motif_row["motif_start"])
        & (match_df["motif_end"] == motif_row["motif_end"])
        & (match_df["set_type"] == motif_row["set_type"])
        & (match_df["in_window"] == True)
        & (match_df["context"].notna())
    ]
    return not subset.empty


def main() -> None:
    args = parse_args()
    output_root = args.output_dir
    output_root.mkdir(parents=True, exist_ok=True)

    ensure_paths_exist(
        [
            args.gff,
            args.tf_positive,
            args.tf_negative,
            args.positive_annot,
            args.negative_annot,
            args.import_dir,
            args.feature_dir,
        ]
    )

    gene_model = load_gene_model(args.gff, args.gene_id)
    motif_df = load_rca_motifs(args, gene_model)
    all_fimo_motif_df = load_all_fimo_motifs(args, gene_model)
    if motif_df.empty:
        raise ValueError(f"No motifs from TF lists overlap {gene_model.gene_id} within ±{args.gene_window_bp} bp")

    all_fimo_hits_across_models: list[pd.DataFrame] = []
    shap_top10_percent_stats: list[dict[str, object]] = []

    for prefix in args.prefixes:
        import_result_path = args.import_dir / f"{prefix}_import_results.txt"
        ensure_paths_exist([import_result_path])

        shap_df = parse_import_results(import_result_path)
        top10_df = shap_df[shap_df["is_top10_percent"] == True].copy()
        for context in ["CHH", "CHG", "CG"]:
            context_df = top10_df[top10_df["context"] == context].copy()
            shap_top10_percent_stats.append(
                {
                    "model": prefix,
                    "context": context,
                    "top10_percent_feature_count": len(top10_df),
                    "context_feature_count": len(context_df),
                    "features": ",".join(context_df["feature"].astype(str).tolist()),
                }
            )
        feature_df = load_feature_positions(prefix, args.feature_dir)
        detail_df = build_motif_feature_matches(
            motif_df=motif_df,
            feature_df=feature_df,
            shap_df=shap_df,
            flank_bp=args.flank_bp,
            dataset=prefix,
        )
        all_fimo_detail_df = build_motif_feature_matches(
            motif_df=all_fimo_motif_df,
            feature_df=feature_df,
            shap_df=shap_df,
            flank_bp=args.flank_bp,
            dataset=prefix,
        )
        all_fimo_hits_across_models.append(all_fimo_detail_df)
        detail_df.to_csv(output_root / f"{prefix}_AT2G39730_motif_gene_matches.tsv", sep="\t", index=False)
        for set_type in ["positive", "negative"]:
            set_motif_df = motif_df[motif_df["set_type"] == set_type].copy()
            set_detail_df = detail_df[detail_df["set_type"] == set_type].copy()
            set_detail_df.to_csv(
                output_root / f"{prefix}_{set_type}_AT2G39730_motif_gene_matches.tsv",
                sep="\t",
                index=False,
            )
            plot_gene_overview(
                gene_model=gene_model,
                motif_df=set_motif_df,
                prefix=prefix,
                set_type=set_type,
                output_dir=output_root,
                gene_window_bp=args.gene_window_bp,
            )
            for suffix in ("png", "pdf"):
                for old_file in output_root.glob(f"{prefix}_{set_type}_*_motif.{suffix}"):
                    old_file.unlink()
            for motif in set_motif_df.sort_values(["motif_center", "symbol", "motif_id"]).itertuples(index=False):
                motif_row = pd.Series(motif._asdict())
                if not motif_has_methylation(motif_row, set_detail_df):
                    continue
                plot_single_motif_figure(
                    motif_row=motif_row,
                    match_df=set_detail_df,
                    prefix=prefix,
                    set_type=set_type,
                    output_dir=output_root,
                    flank_bp=args.flank_bp,
                )
        print(f"[done] {prefix}: outputs written to {output_root}")

    pd.DataFrame(shap_top10_percent_stats).to_csv(
        output_root / "shap_top10_percent_CHH_CHG_CG_counts.tsv",
        sep="\t",
        index=False,
    )

    if all_fimo_hits_across_models:
        all_fimo_union_df = pd.concat(all_fimo_hits_across_models, ignore_index=True)
        all_fimo_union_df = all_fimo_union_df[
            (all_fimo_union_df["in_window"] == True) & (all_fimo_union_df["context"].notna())
        ].copy()
        if not all_fimo_union_df.empty:
            grouped = (
                all_fimo_union_df.groupby(
                    ["motif_id", "symbol", "chrom", "motif_start", "motif_end", "strand", "matched_sequence"],
                    dropna=False,
                )
                .agg(
                    models=("dataset", lambda s: ",".join(sorted(set(map(str, s))))),
                    model_count=("dataset", lambda s: len(set(map(str, s)))),
                    methylation_site_count=("methylation_feature", lambda s: len(set(map(str, s)))),
                    contexts=("context", lambda s: ",".join(sorted(set(map(str, s))))),
                )
                .reset_index()
                .sort_values(["motif_start", "motif_end", "motif_id"])
            )
        else:
            grouped = pd.DataFrame(
                columns=[
                    "motif_id",
                    "symbol",
                    "chrom",
                    "motif_start",
                    "motif_end",
                    "strand",
                    "matched_sequence",
                    "models",
                    "model_count",
                    "methylation_site_count",
                    "contexts",
                ]
            )
        grouped.to_csv(output_root / "all_fimo_any_model_methylated_motifs.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
