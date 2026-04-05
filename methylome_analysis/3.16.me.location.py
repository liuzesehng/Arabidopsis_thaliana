from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import ScalarFormatter
import pandas as pd


PROJECT_ROOT = Path("/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana")
SHARED_ROOT = Path("/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng")
GFF_PATH = SHARED_ROOT / "ref/Arabidopsis_thaliana/refgen/GCF_000001735.4_TAIR10.1_genomic.gff"
BASE_DIR = PROJECT_ROOT / "list/xgboot/TPM_4.5_all"
OUTPUT_DIR = PROJECT_ROOT / "pipline/methylome_analysis/AT2G39730_shap_location_outputs"
TARGET_GENE_ID = "AT2G39730"
GENE_CHR = "Chr2"
FLANK_BP = 2000
DATASETS = ["A_TPM", "A_per_TPM", "B_TPM", "B_per_TPM", "total_TPM"]
CONTEXT_ORDER = ["CHH", "CHG", "CG", "snp"]
REGION_ORDER = ["upstream", "gene_body", "downstream"]
REGION_COLORS = {
    "upstream": "#D95F02",
    "gene_body": "#1B9E77",
    "downstream": "#7570B3",
}
@dataclass
class TranscriptModel:
    transcript_id: str
    orig_transcript_id: str
    chrom: str
    start: int
    end: int
    strand: str
    exons: list[tuple[int, int]]
    cds: list[tuple[int, int]]


@dataclass
class GeneModel:
    gene_id: str
    chrom: str
    start: int
    end: int
    strand: str
    transcripts: list[TranscriptModel]

    @property
    def window_start(self) -> int:
        return self.start - FLANK_BP

    @property
    def window_end(self) -> int:
        return self.end + FLANK_BP


def parse_attributes(attr_text: str) -> dict[str, str]:
    attributes = {}
    for item in attr_text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attributes[key] = value
    return attributes


def load_gene_model(gff_path: Path, gene_id: str) -> GeneModel:
    gene_record = None
    transcripts: dict[str, dict[str, str | int]] = {}
    exon_map: dict[str, list[tuple[int, int]]] = {}
    cds_map: dict[str, list[tuple[int, int]]] = {}

    with gff_path.open() as handle:
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
                    "chrom": seqid,
                    "start": start_i,
                    "end": end_i,
                    "strand": strand,
                }
            elif feature_type == "mRNA":
                transcript_key = attrs["ID"]
                transcripts[transcript_key] = {
                    "transcript_id": attrs.get("transcript_id", transcript_key),
                    "orig_transcript_id": attrs.get("orig_transcript_id", transcript_key),
                    "chrom": seqid,
                    "start": start_i,
                    "end": end_i,
                    "strand": strand,
                }
                exon_map[transcript_key] = []
                cds_map[transcript_key] = []
            elif feature_type == "exon":
                exon_map.setdefault(attrs["Parent"], []).append((start_i, end_i))
            elif feature_type == "CDS":
                cds_map.setdefault(attrs["Parent"], []).append((start_i, end_i))

    if gene_record is None:
        raise ValueError(f"Gene {gene_id} not found in {gff_path}")

    transcript_models = []
    for transcript_key, meta in sorted(transcripts.items(), key=lambda item: item[1]["orig_transcript_id"]):
        transcript_models.append(
            TranscriptModel(
                transcript_id=str(meta["transcript_id"]),
                orig_transcript_id=str(meta["orig_transcript_id"]),
                chrom=str(meta["chrom"]),
                start=int(meta["start"]),
                end=int(meta["end"]),
                strand=str(meta["strand"]),
                exons=sorted(exon_map.get(transcript_key, [])),
                cds=sorted(cds_map.get(transcript_key, [])),
            )
        )

    return GeneModel(
        gene_id=gene_id,
        chrom=str(gene_record["chrom"]),
        start=int(gene_record["start"]),
        end=int(gene_record["end"]),
        strand=str(gene_record["strand"]),
        transcripts=transcript_models,
    )


def parse_import_results(import_txt: Path) -> pd.DataFrame:
    pattern = re.compile(r"\s+([A-Za-z0-9_]+):\s+([-+0-9.eE]+)")
    records = []
    in_section = False

    with import_txt.open() as handle:
        for line in handle:
            if "所有特征及其平均|SHAP|值:" in line:
                in_section = True
                continue
            if in_section and not line.startswith("    "):
                break
            if not in_section:
                continue
            match = pattern.match(line)
            if match:
                feature = match.group(1)
                shap_value = float(match.group(2))
                if shap_value > 0:
                    records.append({"feature": feature, "shap_value": shap_value})

    if not records:
        raise ValueError(f"No positive SHAP features found in {import_txt}")

    shap_df = pd.DataFrame(records).sort_values("shap_value", ascending=False).reset_index(drop=True)
    shap_df["context"] = shap_df["feature"].str.split("_").str[0]
    return shap_df


def load_feature_positions(dataset_prefix: str, base_dir: Path) -> pd.DataFrame:
    dataset_dir = base_dir / "feature_data_extraction" / dataset_prefix
    frames = []
    file_specs = {
        "snp": dataset_dir / f"{dataset_prefix}_snp_data.csv",
        "CG": dataset_dir / f"{dataset_prefix}_CG_data.csv",
        "CHG": dataset_dir / f"{dataset_prefix}_CHG_data.csv",
        "CHH": dataset_dir / f"{dataset_prefix}_CHH_data.csv",
    }

    for context, file_path in file_specs.items():
        df = pd.read_csv(file_path)
        if context == "snp":
            frame = pd.DataFrame(
                {
                    "feature": df["Feature"],
                    "chrom": df["#Chromosome"],
                    "start": df["Position"].astype(int),
                    "end": df["Position"].astype(int),
                    "context": context,
                }
            )
        else:
            # Methylation feature tables use BED-style coordinates:
            # start is 0-based and end is 1-based exclusive. Convert to
            # GFF-style 1-based inclusive coordinates before comparison.
            frame = pd.DataFrame(
                {
                    "feature": df["Feature"],
                    "chrom": df["chr"],
                    "start": df["start"].astype(int) + 1,
                    "end": df["end"].astype(int),
                    "context": context,
                }
            )
        frames.append(frame)

    return pd.concat(frames, ignore_index=True)


def annotate_regions(features_df: pd.DataFrame, gene_model: GeneModel, flank_bp: int = FLANK_BP) -> pd.DataFrame:
    df = features_df.copy()
    df["plot_pos"] = ((df["start"] + df["end"]) / 2).astype(float)

    window_start = gene_model.start - flank_bp
    window_end = gene_model.end + flank_bp
    df = df[(df["chrom"] == GENE_CHR) & (df["end"] >= window_start) & (df["start"] <= window_end)].copy()

    def classify_region(row: pd.Series) -> str:
        if row["start"] >= gene_model.start and row["end"] <= gene_model.end:
            return "gene_body"
        if gene_model.strand == "-":
            if row["start"] > gene_model.end and row["end"] <= gene_model.end + flank_bp:
                return "upstream"
            if row["end"] < gene_model.start and row["start"] >= gene_model.start - flank_bp:
                return "downstream"
        else:
            if row["end"] < gene_model.start and row["start"] >= gene_model.start - flank_bp:
                return "upstream"
            if row["start"] > gene_model.end and row["end"] <= gene_model.end + flank_bp:
                return "downstream"
        return "outside_window"

    df["region_class"] = df.apply(classify_region, axis=1)
    df = df[df["region_class"].isin(REGION_ORDER)].copy()
    return df


def build_dataset_feature_table(dataset_prefix: str, base_dir: Path, gene_model: GeneModel) -> pd.DataFrame:
    shap_df = parse_import_results(base_dir / f"{dataset_prefix}_import_results.txt")
    position_df = load_feature_positions(dataset_prefix, base_dir)
    merged = shap_df.merge(position_df, on=["feature", "context"], how="inner", validate="one_to_one")
    annotated = annotate_regions(merged, gene_model, flank_bp=FLANK_BP)
    annotated["dataset"] = dataset_prefix
    annotated["context"] = pd.Categorical(annotated["context"], categories=CONTEXT_ORDER, ordered=True)
    annotated = annotated[
        ["dataset", "context", "feature", "chrom", "start", "end", "plot_pos", "shap_value", "region_class"]
    ].sort_values(["context", "plot_pos", "feature"])
    return annotated.reset_index(drop=True)


def summarize_dataset(dataset_df: pd.DataFrame, dataset_prefix: str) -> pd.DataFrame:
    total_shap = dataset_df["shap_value"].sum()
    total_feature_count = len(dataset_df)
    rows = []
    context_total_counts = dataset_df.groupby("context", observed=False).size().to_dict()

    overall = dataset_df.groupby("region_class", observed=False)["shap_value"].sum().reindex(REGION_ORDER, fill_value=0.0)
    for region in REGION_ORDER:
        shap_sum = float(overall.loc[region])
        feature_count = int((dataset_df["region_class"] == region).sum())
        rows.append(
            {
                "dataset": dataset_prefix,
                "context": "all",
                "region_class": region,
                "feature_count": feature_count,
                "shap_sum": shap_sum,
                "shap_percent": (shap_sum / total_shap * 100.0) if total_shap else 0.0,
                "feature_percent": (feature_count / total_feature_count * 100.0) if total_feature_count else 0.0,
            }
        )

    grouped = (
        dataset_df.groupby(["context", "region_class"], observed=False)["shap_value"]
        .agg(["sum", "count"])
        .reset_index()
        .rename(columns={"sum": "shap_sum", "count": "feature_count"})
    )
    for context in CONTEXT_ORDER:
        context_subset = grouped[grouped["context"] == context].set_index("region_class")
        for region in REGION_ORDER:
            if region in context_subset.index:
                shap_sum = float(context_subset.loc[region, "shap_sum"])
                feature_count = int(context_subset.loc[region, "feature_count"])
            else:
                shap_sum = 0.0
                feature_count = 0
            context_total_count = int(context_total_counts.get(context, 0))
            rows.append(
                {
                    "dataset": dataset_prefix,
                    "context": context,
                    "region_class": region,
                    "feature_count": feature_count,
                    "shap_sum": shap_sum,
                    "shap_percent": (shap_sum / total_shap * 100.0) if total_shap else 0.0,
                    "feature_percent": (feature_count / context_total_count * 100.0) if context_total_count else 0.0,
                }
            )

    return pd.DataFrame(rows)


def add_region_background(ax: plt.Axes, gene_model: GeneModel) -> None:
    region_spans = []
    if gene_model.strand == "-":
        region_spans = [
            ("downstream", gene_model.start - FLANK_BP, gene_model.start),
            ("gene_body", gene_model.start, gene_model.end),
            ("upstream", gene_model.end, gene_model.end + FLANK_BP),
        ]
    else:
        region_spans = [
            ("upstream", gene_model.start - FLANK_BP, gene_model.start),
            ("gene_body", gene_model.start, gene_model.end),
            ("downstream", gene_model.end, gene_model.end + FLANK_BP),
        ]

    for region, start, end in region_spans:
        ax.axvspan(start, end, color=REGION_COLORS[region], alpha=0.08, lw=0)


def draw_transcript(ax: plt.Axes, transcript: TranscriptModel, y_pos: float) -> None:
    ax.hlines(y_pos, transcript.start, transcript.end, color="#303030", linewidth=1.8, zorder=2)

    for start, end in transcript.exons:
        ax.add_patch(
            Rectangle((start, y_pos - 0.18), end - start, 0.36, facecolor="#BFC9CA", edgecolor="#566573", lw=0.8, zorder=3)
        )
    for start, end in transcript.cds:
        ax.add_patch(
            Rectangle((start, y_pos - 0.22), end - start, 0.44, facecolor="#2C3E50", edgecolor="#1B2631", lw=0.8, zorder=4)
        )

    arrow_start = min(transcript.end, transcript.start + 260)
    arrow_end = max(transcript.start, arrow_start - 180)
    ax.annotate(
        "",
        xy=(arrow_end, y_pos),
        xytext=(arrow_start, y_pos),
        arrowprops={"arrowstyle": "-|>", "lw": 1.2, "color": "#1F1F1F"},
        zorder=5,
    )

def plot_dataset(
    dataset_prefix: str,
    dataset_df: pd.DataFrame,
    stats_df: pd.DataFrame,
    gene_model: GeneModel,
    output_png_path: Path,
    output_pdf_path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(14, 6.5))
    ax.set_xlim(gene_model.start - FLANK_BP, gene_model.end + FLANK_BP)
    ax.set_ylim(0.5, 7.8)
    add_region_background(ax, gene_model)

    transcript_y = [7.0, 6.2, 5.4]
    for transcript, y_pos in zip(gene_model.transcripts, transcript_y):
        draw_transcript(ax, transcript, y_pos)

    context_y = {"CHH": 4.2, "CHG": 3.2, "CG": 2.2, "snp": 1.2}
    for context, y_pos in context_y.items():
        ax.hlines(y_pos, gene_model.start - FLANK_BP, gene_model.end + FLANK_BP, color="#D0D0D0", linewidth=0.8, zorder=1)
        ax.text(
            ax.get_xlim()[0] - 220,
            y_pos,
            context.upper() if context == "snp" else context,
            ha="right",
            va="center",
            fontsize=13,
        )
        subset = dataset_df[dataset_df["context"] == context]
        if subset.empty:
            continue
        ax.scatter(
            subset["plot_pos"],
            [y_pos] * len(subset),
            c=subset["region_class"].map(REGION_COLORS),
            s=(subset["shap_value"] / subset["shap_value"].max() * 120).clip(lower=25),
            edgecolors="#222222",
            linewidths=0.4,
            alpha=0.95,
            zorder=6,
        )

    context_stats = {
        context: stats_df[stats_df["context"] == context].set_index("region_class").reindex(REGION_ORDER)
        for context in CONTEXT_ORDER
    }

    legend_handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label=region,
            markerfacecolor=color,
            markeredgecolor="#222222",
            markersize=8,
        )
        for region, color in REGION_COLORS.items()
    ]
    ax.legend(
        handles=legend_handles,
        loc="upper left",
        frameon=False,
        title="Region",
        fontsize=11,
        title_fontsize=12,
    )

    inset_ax = fig.add_axes([0.77, 0.62, 0.21, 0.27])
    x = list(range(len(CONTEXT_ORDER)))
    bar_width = 0.22
    region_offsets = {"upstream": -bar_width, "gene_body": 0.0, "downstream": bar_width}

    for region in REGION_ORDER:
        heights = [float(context_stats[context].loc[region, "feature_percent"]) for context in CONTEXT_ORDER]
        bar_positions = [idx + region_offsets[region] for idx in x]
        bars = inset_ax.bar(
            bar_positions,
            heights,
            width=bar_width,
            color=REGION_COLORS[region],
            edgecolor="#333333",
            linewidth=0.5,
        )
        for bar, height in zip(bars, heights):
            inset_ax.text(
                bar.get_x() + bar.get_width() / 2,
                height + 1.5,
                f"{height:.1f}",
                ha="center",
                va="bottom",
                fontsize=7.5,
                color="#222222",
            )

    inset_ax.set_xticks(x)
    inset_ax.set_xticklabels(["CHH", "CHG", "CG", "SNP"], fontsize=9)
    inset_ax.text(0.0, 1.03, "%", transform=inset_ax.transAxes, ha="center", va="bottom", fontsize=9)
    inset_ax.set_title("Feature Percent", fontsize=10)
    inset_ax.tick_params(axis="y", labelsize=8)
    inset_ax.set_ylim(0, 110)
    inset_ax.grid(axis="y", color="#D8D8D8", linewidth=0.6, alpha=0.8)
    inset_ax.set_axisbelow(True)
    inset_ax.set_facecolor("white")
    for spine in ["top", "right"]:
        inset_ax.spines[spine].set_visible(False)

    formatter = ScalarFormatter(useOffset=False)
    formatter.set_scientific(False)
    ax.xaxis.set_major_formatter(formatter)
    ax.ticklabel_format(style="plain", axis="x", useOffset=False)

    ax.set_xlabel(f"Genomic position on {GENE_CHR}", fontsize=13)
    ax.tick_params(axis="x", labelsize=12)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    fig.subplots_adjust(left=0.08, right=0.98, top=0.96, bottom=0.12)
    fig.savefig(output_png_path, dpi=300, bbox_inches="tight")
    fig.savefig(output_pdf_path, bbox_inches="tight")
    plt.close(fig)


def write_outputs(dataset_prefix: str, dataset_df: pd.DataFrame, stats_df: pd.DataFrame, output_dir: Path) -> None:
    dataset_df.to_csv(output_dir / f"{dataset_prefix}_feature_details.csv", index=False)
    stats_df.to_csv(output_dir / f"{dataset_prefix}_region_stats.csv", index=False)


def validate_gene_model(gene_model: GeneModel) -> None:
    if len(gene_model.transcripts) != 3:
        raise AssertionError(f"Expected 3 transcripts, found {len(gene_model.transcripts)}")

    expected_exon_counts = {
        "gnl|JCVI|mRNA.AT2G39730.1": 7,
        "gnl|JCVI|mRNA.AT2G39730.2": 7,
        "gnl|JCVI|mRNA.AT2G39730.3": 6,
    }
    expected_cds_counts = expected_exon_counts
    for transcript in gene_model.transcripts:
        if len(transcript.exons) != expected_exon_counts[transcript.orig_transcript_id]:
            raise AssertionError(f"Unexpected exon count for {transcript.orig_transcript_id}")
        if len(transcript.cds) != expected_cds_counts[transcript.orig_transcript_id]:
            raise AssertionError(f"Unexpected CDS count for {transcript.orig_transcript_id}")

    transcript_three = next(t for t in gene_model.transcripts if t.orig_transcript_id.endswith(".3"))
    other_ends = [t.end for t in gene_model.transcripts if t.orig_transcript_id != transcript_three.orig_transcript_id]
    if not all(transcript_three.end < end for end in other_ends):
        raise AssertionError("AT2G39730.3 should terminate earlier than the other transcripts")


def validate_dataset(dataset_prefix: str, dataset_df: pd.DataFrame, stats_df: pd.DataFrame) -> None:
    overall = stats_df[stats_df["context"] == "all"].set_index("region_class").reindex(REGION_ORDER)
    percent_sum = overall["feature_percent"].sum()
    if abs(percent_sum - 100.0) > 1e-6:
        raise AssertionError(f"{dataset_prefix} percentages do not sum to 100: {percent_sum}")

    if not dataset_df["region_class"].isin(REGION_ORDER).all():
        raise AssertionError(f"{dataset_prefix} contains invalid region labels")

    if dataset_df["context"].nunique() < 3:
        raise AssertionError(f"{dataset_prefix} lost too many contexts after filtering")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    gene_model = load_gene_model(GFF_PATH, TARGET_GENE_ID)
    validate_gene_model(gene_model)

    all_stats = []
    dataset_summaries = []

    for dataset_prefix in DATASETS:
        dataset_df = build_dataset_feature_table(dataset_prefix, BASE_DIR, gene_model)
        stats_df = summarize_dataset(dataset_df, dataset_prefix)
        validate_dataset(dataset_prefix, dataset_df, stats_df)
        write_outputs(dataset_prefix, dataset_df, stats_df, OUTPUT_DIR)
        plot_dataset(
            dataset_prefix,
            dataset_df,
            stats_df,
            gene_model,
            OUTPUT_DIR / f"{dataset_prefix}_feature_map.png",
            OUTPUT_DIR / f"{dataset_prefix}_feature_map.pdf",
        )
        all_stats.append(stats_df)

        overall = stats_df[stats_df["context"] == "all"].set_index("region_class").reindex(REGION_ORDER)
        dataset_summaries.append(
            {
                "dataset": dataset_prefix,
                "feature_count": int(len(dataset_df)),
                "total_shap_sum": float(dataset_df["shap_value"].sum()),
                "upstream_percent": float(overall.loc["upstream", "feature_percent"]),
                "gene_body_percent": float(overall.loc["gene_body", "feature_percent"]),
                "downstream_percent": float(overall.loc["downstream", "feature_percent"]),
            }
        )

    pd.concat(all_stats, ignore_index=True).to_csv(OUTPUT_DIR / "AT2G39730_region_context_summary.csv", index=False)
    pd.DataFrame(dataset_summaries).to_csv(OUTPUT_DIR / "AT2G39730_dataset_summary.csv", index=False)
    print(f"Finished. Outputs written to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
