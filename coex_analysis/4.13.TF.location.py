#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


DEFAULT_WORKDIR = Path(
    "/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana/list/RCA"
)
PROJECT_ROOT = Path(
    "/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana"
)
DEFAULT_TF_POSITIVE = PROJECT_ROOT / "list/coex2/total/gene.TF_ids.txt"
DEFAULT_TF_NEGATIVE = PROJECT_ROOT / "list/coex2/total/gene.TF_negative_ids.txt"
DEFAULT_POSITIVE_ANNOT = PROJECT_ROOT / "list/coex2/total/total.gene.TF_ids.annotated.txt"
DEFAULT_NEGATIVE_ANNOT = PROJECT_ROOT / "list/coex2/total/total.negative.gene.TF_ids.annotated.txt"
DEFAULT_IMPORT_DIR = PROJECT_ROOT / "list/xgboot/TPM_4.5_all"
DEFAULT_FEATURE_DIR = DEFAULT_IMPORT_DIR / "feature_data_extraction"
DEFAULT_PREFIXES = ["A_TPM", "B_TPM", "total_TPM", "A_per_TPM", "B_per_TPM"]
CONTEXT_COLORS = {"CHH": "#5e3c99", "CHG": "#1b9e77", "CG": "#d95f02"}
BASE_COLORS = {"A": "#2c7fb8", "T": "#d95f0e", "C": "#238b45", "G": "#7b3294", "N": "#4d4d4d"}
ARABIDOPSIS_CHR_MAP = {
    "NC_003070": "Chr1",
    "NC_003071": "Chr2",
    "NC_003074": "Chr3",
    "NC_003075": "Chr4",
    "NC_003076": "Chr5",
    "NC_001284": "ChrC",
    "NC_000932": "ChrM",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot SHAP-selected methylation sites around TF motif hits."
    )
    parser.add_argument("--workdir", type=Path, default=DEFAULT_WORKDIR)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--prefixes", nargs="+", default=DEFAULT_PREFIXES)
    parser.add_argument("--flank-bp", type=int, default=50)
    parser.add_argument("--tf-positive", type=Path, default=DEFAULT_TF_POSITIVE)
    parser.add_argument("--tf-negative", type=Path, default=DEFAULT_TF_NEGATIVE)
    parser.add_argument("--positive-annot", type=Path, default=DEFAULT_POSITIVE_ANNOT)
    parser.add_argument("--negative-annot", type=Path, default=DEFAULT_NEGATIVE_ANNOT)
    parser.add_argument("--import-dir", type=Path, default=DEFAULT_IMPORT_DIR)
    parser.add_argument("--feature-dir", type=Path, default=DEFAULT_FEATURE_DIR)
    return parser.parse_args()


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


def parse_import_results(path: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    in_shap_section = False
    pattern = re.compile(r"^\s*([A-Za-z0-9_]+):\s*([+-]?\d+(?:\.\d+)?)\s*$")

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
            context = feature.split("_", 1)[0].upper()
            if context not in CONTEXT_COLORS:
                continue
            if shap_value <= 0:
                continue

            rows.append({"feature": feature, "context": context, "shap_value": shap_value})

    if not rows:
        return pd.DataFrame(columns=["feature", "context", "shap_value"])
    return pd.DataFrame(rows).drop_duplicates(subset=["feature"], keep="first")


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
    return feature_df.dropna(subset=["plot_start", "plot_end", "position"])


def load_motif_hits(path: Path, annotated_path: Path) -> pd.DataFrame:
    motif_df = pd.read_csv(path, sep="\t")
    annot_df = pd.read_csv(annotated_path, sep="\t")
    symbol_map = dict(zip(annot_df["Gene_ID"].astype(str), annot_df["Symbol"].fillna("").astype(str)))

    motif_df = motif_df.copy()
    motif_df["chrom"] = motif_df["sequence_name"].map(normalize_chromosome)
    motif_df["motif_start"] = pd.to_numeric(motif_df["start"], errors="coerce").astype("Int64")
    motif_df["motif_end"] = pd.to_numeric(motif_df["stop"], errors="coerce").astype("Int64")
    motif_df["strand"] = motif_df["strand"].fillna(".").astype(str)
    motif_df["matched_sequence"] = motif_df["matched_sequence"].fillna("").astype(str).str.upper()
    motif_df["symbol"] = motif_df["motif_id"].astype(str).map(symbol_map).fillna("")
    motif_df["symbol"] = motif_df.apply(
        lambda row: row["symbol"] if str(row["symbol"]).strip() else row["motif_id"], axis=1
    )
    return motif_df.dropna(subset=["motif_start", "motif_end"]).reset_index(drop=True)


def build_motif_feature_matches(
    motif_df: pd.DataFrame, feature_df: pd.DataFrame, shap_df: pd.DataFrame, flank_bp: int
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
                    "motif_id": motif.motif_id,
                    "symbol": motif.symbol,
                    "matched_sequence": motif.matched_sequence,
                    "chrom": motif.chrom,
                    "strand": motif.strand,
                    "motif_start": int(motif.motif_start),
                    "motif_end": int(motif.motif_end),
                    "window_start": window_start,
                    "window_end": window_end,
                    "methylation_feature": pd.NA,
                    "context": pd.NA,
                    "shap_value": pd.NA,
                    "position": pd.NA,
                    "plot_start": pd.NA,
                    "plot_end": pd.NA,
                    "in_window": False,
                    "in_motif": False,
                }
            )
            continue

        for row in subset.itertuples(index=False):
            position = float(row.position)
            records.append(
                {
                    "motif_id": motif.motif_id,
                    "symbol": motif.symbol,
                    "matched_sequence": motif.matched_sequence,
                    "chrom": motif.chrom,
                    "strand": motif.strand,
                    "motif_start": int(motif.motif_start),
                    "motif_end": int(motif.motif_end),
                    "window_start": window_start,
                    "window_end": window_end,
                    "methylation_feature": row.feature,
                    "context": row.context,
                    "shap_value": float(row.shap_value),
                    "position": position,
                    "plot_start": int(row.plot_start),
                    "plot_end": int(row.plot_end),
                    "in_window": True,
                    "in_motif": int(motif.motif_start) <= position <= int(motif.motif_end),
                }
            )

    match_df = pd.DataFrame(records)
    if match_df.empty:
        return pd.DataFrame(
            columns=[
                "motif_id",
                "symbol",
                "matched_sequence",
                "chrom",
                "strand",
                "motif_start",
                "motif_end",
                "window_start",
                "window_end",
                "methylation_feature",
                "context",
                "shap_value",
                "position",
                "plot_start",
                "plot_end",
                "in_window",
                "in_motif",
            ]
        )
    return match_df


def compute_display_positions(motif_row: dict[str, object], valid_points: pd.DataFrame) -> dict[float, float]:
    window_start = float(motif_row["window_start"])
    window_end = float(motif_row["window_end"])
    motif_start = int(motif_row["motif_start"])
    motif_end = int(motif_row["motif_end"])

    anchors: set[float] = {window_start, float(motif_start), float(motif_end), window_end}
    anchors.update(float(pos) for pos in valid_points["position"].dropna().tolist())
    anchors.update(float(pos) for pos in range(motif_start, motif_end + 1))
    sorted_anchors = sorted(anchors)

    display_map: dict[float, float] = {sorted_anchors[0]: 0.0}
    for prev, curr in zip(sorted_anchors, sorted_anchors[1:]):
        gap = curr - prev
        prev_in_motif = motif_start <= prev <= motif_end
        curr_in_motif = motif_start <= curr <= motif_end
        if prev_in_motif and curr_in_motif:
            increment = 2.7
        elif gap <= 1:
            increment = 2.2
        elif gap <= 3:
            increment = 2.6
        elif gap <= 8:
            increment = 3.0
        else:
            increment = 3.1 + min(2.2, math.log1p(gap) * 0.55)
        display_map[curr] = display_map[prev] + increment
    return display_map


def draw_motif_sequence(ax: plt.Axes, motif_row: pd.Series, display_map: dict[float, float]) -> None:
    seq = str(motif_row["matched_sequence"])
    motif_start = int(motif_row["motif_start"])
    motif_end = int(motif_row["motif_end"])
    is_negative = str(motif_row["strand"]) == "-"

    ax.axvspan(display_map[float(motif_start)] - 1.0, display_map[float(motif_end)] + 1.0, color="#f1eef6", alpha=0.8, zorder=0)
    for offset, base in enumerate(seq):
        genomic_pos = motif_end - offset if is_negative else motif_start + offset
        xpos = display_map[float(genomic_pos)]
        y_offset = 1.12 + (0.16 if offset % 2 else 0.0)
        ax.text(
            xpos,
            y_offset,
            base,
            ha="center",
            va="center",
            fontsize=26,
            fontweight="bold",
            color=BASE_COLORS.get(base, BASE_COLORS["N"]),
            zorder=3,
        )


def plot_motif_panels(match_df: pd.DataFrame, output_prefix: Path, set_name: str, prefix: str, flank_bp: int) -> None:
    motif_order = (
        match_df[["motif_id", "symbol", "matched_sequence", "chrom", "strand", "motif_start", "motif_end", "window_start", "window_end"]]
        .drop_duplicates()
        .sort_values(["motif_start", "motif_id"])
        .reset_index(drop=True)
    )
    if motif_order.empty:
        return

    fig_height = max(2.4 * len(motif_order), 3.0)
    fig, axes = plt.subplots(len(motif_order), 1, figsize=(18, fig_height), squeeze=False)
    axes_flat = axes.flatten()
    y_positions = {"CHH": 2.0, "CHG": 2.45, "CG": 2.9}

    for ax, motif_row in zip(axes_flat, motif_order.to_dict("records")):
        motif_slice = match_df[match_df["motif_id"] == motif_row["motif_id"]].copy()
        valid_points = motif_slice[motif_slice["in_window"] == True].copy()
        display_map = compute_display_positions(motif_row, valid_points)
        x_min = display_map[float(motif_row["window_start"])]
        x_max = display_map[float(motif_row["window_end"])]

        ax.set_xlim(x_min, x_max)
        if str(motif_row["strand"]) == "-":
            ax.invert_xaxis()
        ax.set_ylim(-2.35, 3.35)
        ax.axhline(0, color="black", linewidth=1.2)
        ax.axvline(display_map[float(motif_row["motif_start"])], color="#666666", linestyle="--", linewidth=1.0)
        ax.axvline(display_map[float(motif_row["motif_end"])], color="#666666", linestyle="--", linewidth=1.0)

        draw_motif_sequence(ax, pd.Series(motif_row), display_map)
        loc_text = (
            f"{motif_row['chrom']}:{motif_row['motif_end']}-{motif_row['motif_start']}"
            if str(motif_row["strand"]) == "-"
            else f"{motif_row['chrom']}:{motif_row['motif_start']}-{motif_row['motif_end']}"
        )
        ax.text(
            0.01,
            0.975,
            f"{motif_row['symbol']} ({motif_row['strand']}) {loc_text}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=11,
            fontweight="bold",
            color="#222222",
        )

        if valid_points.empty:
            ax.text(
                (x_min + x_max) / 2,
                -0.55,
                f"no SHAP-selected methylation site in ±{flank_bp} bp",
                ha="center",
                va="center",
                fontsize=12,
                color="#666666",
            )
        else:
            valid_points = valid_points.sort_values(["position", "context", "methylation_feature"]).reset_index(drop=True)
            for idx, row in enumerate(valid_points.itertuples(index=False)):
                ypos = y_positions[str(row.context)]
                xpos = display_map[float(row.position)]
                color = CONTEXT_COLORS[str(row.context)]
                ax.vlines(xpos, 0.1, ypos - 0.1, color=color, linewidth=1.5, alpha=0.9, zorder=2)
                ax.scatter([xpos], [ypos], s=90, color=color, edgecolor="black", linewidth=0.7, zorder=3)
                label_y = -1.0 if idx % 2 == 0 else -1.24
                ax.text(
                    xpos,
                    label_y,
                    str(row.methylation_feature),
                    ha="center",
                    va="center",
                    fontsize=8,
                    color="#222222",
                )
                if bool(row.in_motif):
                    ax.text(
                        xpos,
                        -1.95,
                        "★",
                        ha="center",
                        va="center",
                        fontsize=22,
                        color="#fdb863",
                        zorder=5,
                    )

        ax.set_xticks([])
        ax.set_yticks([])
        ax.tick_params(axis="x", length=0, labelbottom=False)
        ax.tick_params(axis="y", length=0)
        ax.text(0.99, 0.975, f"±{flank_bp} bp", transform=ax.transAxes, fontsize=10, ha="right", va="top", color="#666666")
        for spine in ["left", "right"]:
            ax.spines[spine].set_visible(False)

    fig.tight_layout(rect=(0, 0, 1, 1))
    fig.savefig(output_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(output_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def process_set(
    prefix: str,
    set_name: str,
    motif_path: Path,
    annotated_path: Path,
    shap_df: pd.DataFrame,
    feature_df: pd.DataFrame,
    output_dir: Path,
    flank_bp: int,
) -> pd.DataFrame:
    motif_df = load_motif_hits(motif_path, annotated_path)
    match_df = build_motif_feature_matches(motif_df, feature_df, shap_df, flank_bp)
    if match_df.empty:
        return match_df

    match_df = match_df.copy()
    match_df.insert(0, "prefix", prefix)
    match_df.insert(1, "set_type", set_name)

    figure_prefix = output_dir / f"{prefix}_{set_name.lower()}_tf_motif_methylation"
    plot_motif_panels(match_df, figure_prefix, set_name, prefix, flank_bp)
    return match_df


def ensure_paths_exist(paths: list[Path]) -> None:
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required files:\n" + "\n".join(missing))


def main() -> None:
    args = parse_args()
    output_root = args.output_dir or (args.workdir / "TF_location_outputs")
    output_root.mkdir(parents=True, exist_ok=True)

    ensure_paths_exist(
        [
            args.tf_positive,
            args.tf_negative,
            args.positive_annot,
            args.negative_annot,
            args.import_dir,
            args.feature_dir,
        ]
    )

    for prefix in args.prefixes:
        import_result_path = args.import_dir / f"{prefix}_import_results.txt"
        ensure_paths_exist([import_result_path])

        shap_df = parse_import_results(import_result_path)
        feature_df = load_feature_positions(prefix, args.feature_dir)
        prefix_output_dir = output_root / prefix
        prefix_output_dir.mkdir(parents=True, exist_ok=True)

        all_details: list[pd.DataFrame] = []
        for set_name, motif_path, annot_path in [
            ("positive", args.tf_positive, args.positive_annot),
            ("negative", args.tf_negative, args.negative_annot),
        ]:
            detail_df = process_set(
                prefix=prefix,
                set_name=set_name,
                motif_path=motif_path,
                annotated_path=annot_path,
                shap_df=shap_df,
                feature_df=feature_df,
                output_dir=prefix_output_dir,
                flank_bp=args.flank_bp,
            )
            if not detail_df.empty:
                all_details.append(detail_df)

        if all_details:
            combined_df = pd.concat(all_details, ignore_index=True)
            combined_df.to_csv(prefix_output_dir / f"{prefix}_tf_motif_methylation_details.tsv", sep="\t", index=False)
        else:
            pd.DataFrame().to_csv(prefix_output_dir / f"{prefix}_tf_motif_methylation_details.tsv", sep="\t", index=False)

        print(f"[done] {prefix}: outputs written to {prefix_output_dir}")


if __name__ == "__main__":
    main()
