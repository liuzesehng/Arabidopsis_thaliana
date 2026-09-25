#!/usr/bin/env python3
"""Reconstruct the historical, unmerged TF–CHH analysis and draw three panels.

Counts describe motif occurrences, TF symbols, or feature–TF pairs, as labelled;
they are not replicate means and have no error bars or significance tests.
The historical tables are validation references, never plotting constants.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "tf-chh-mpl"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "tf-chh-cache"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
from matplotlib.transforms import ScaledTranslation
import numpy as np
import pandas as pd


DEFAULT_BASE = Path(__file__).resolve().parents[2] / "list" / "RCA"
DEFAULT_PREFIX = "AT2G39730_CHH_combined_historical"
DEFAULT_QA_TOOLS = Path.home() / ".codex" / "skills" / "nature-figure" / "scripts"
SET_ORDER = ["negative", "positive"]
REGION_ORDER = ["in", "flanking"]
SET_LABELS = {"negative": "non-co-exp", "positive": "co-exp"}
COLORS = {"in": "#2A9D8F", "flanking": "#E76F51"}
# Only occurrence-level identifiers belong in this key, never SHAP ranks.
MOTIF_KEY = [
    "dataset", "set_type", "motif_id", "symbol", "strand", "chrom",
    "motif_start", "motif_end", "matched_sequence", "window_start", "window_end",
]
REQUIRED = set(MOTIF_KEY + [
    "context", "methylation_feature", "feature_pos", "shap_value",
    "in_motif", "in_window",
])


def require_columns(df: pd.DataFrame, columns: set[str], source: str) -> None:
    missing = columns.difference(df.columns)
    if missing:
        raise ValueError(f"{source}: missing required columns: {sorted(missing)}")


def parse_bool(series: pd.Series, source: str) -> pd.Series:
    values = series.astype("string").str.strip().str.lower()
    mapping = {"true": True, "1": True, "yes": True, "y": True,
               "false": False, "0": False, "no": False, "n": False}
    bad = values.isna() | ~values.isin(mapping)
    if bad.any():
        raise ValueError(f"{source}: invalid Boolean values at rows "
                         f"{(np.flatnonzero(bad.to_numpy()) + 2).tolist()}")
    return values.map(mapping).astype(bool)


def prepare_hits(df: pd.DataFrame, symbol_reference: pd.DataFrame) -> pd.DataFrame:
    require_columns(df, REQUIRED, "Match table")
    require_columns(symbol_reference, {"motif_id", "symbol"}, "Symbol reference")
    df = df.copy()
    for column in ["in_motif", "in_window"]:
        df[column] = parse_bool(df[column], column)
    df = df.loc[df["context"].astype("string").str.strip().eq("CHH")
                & (df["in_motif"] | df["in_window"])].copy()
    if df.empty:
        raise ValueError("No CHH motif/window matches in the supplied historical data")
    for column in ["dataset", "set_type", "motif_id", "strand", "chrom",
                   "matched_sequence", "methylation_feature"]:
        if df[column].isna().any() or df[column].astype(str).str.strip().eq("").any():
            raise ValueError(f"CHH hits: missing identifier in {column}")
    if not df["set_type"].isin(SET_ORDER).all():
        raise ValueError("CHH hits: set_type must be negative or positive")
    for column in ["motif_start", "motif_end", "window_start", "window_end", "feature_pos"]:
        df[column] = pd.to_numeric(df[column], errors="raise")
        if not np.isfinite(df[column]).all():
            raise ValueError(f"CHH hits: non-finite {column}")
    if (df["motif_start"] > df["motif_end"]).any():
        raise ValueError("CHH hits: motif_start exceeds motif_end")
    symbols = symbol_reference[["motif_id", "symbol"]].drop_duplicates()
    if symbols.isna().any().any() or symbols.astype(str).apply(lambda s: s.str.strip().eq("")).any().any():
        raise ValueError("Symbol reference: missing motif_id or symbol")
    if symbols["motif_id"].duplicated().any():
        raise ValueError("Symbol reference: conflicting motif_id to symbol mapping")
    df["symbol"] = df["motif_id"].map(symbols.set_index("motif_id")["symbol"])
    if df["symbol"].isna().any():
        missing = df.loc[df["symbol"].isna(), "motif_id"].unique().tolist()
        raise ValueError(f"Symbol reference: missing mapping for {missing}")
    df["region"] = np.where(df["in_motif"], "in", "flanking")
    return df.reset_index(drop=True)


def build_tables(hits: pd.DataFrame, symbol_order: list[str]) -> tuple[pd.DataFrame, ...]:
    motifs = hits.groupby(MOTIF_KEY, as_index=False, dropna=False).agg(
        in_motif=("in_motif", "any"), in_window=("in_window", "any")
    )
    motifs["region"] = np.where(motifs["in_motif"], "in", "flanking")
    counts = (motifs.groupby(["region", "set_type"]).size()
              .reindex(pd.MultiIndex.from_product([REGION_ORDER, SET_ORDER],
                       names=["region", "set_type"]), fill_value=0)
              .rename("count").reset_index())
    counts["label"] = counts["set_type"].map(SET_LABELS)
    counts = counts[["region", "set_type", "label", "count"]]

    points = hits[["motif_id", "symbol", "set_type", "region"]].copy()
    points["_priority"] = points["region"].map({"in": 0, "flanking": 1})
    points = (points.sort_values("_priority", kind="stable")
              .drop_duplicates(["symbol", "set_type"]).drop(columns="_priority"))
    points["_order"] = points["symbol"].map({s: i for i, s in enumerate(symbol_order)})
    if points["_order"].isna().any():
        raise ValueError("Historical symbol order does not cover every TF")
    points = points.sort_values(["_order", "set_type"]).drop(columns="_order").reset_index(drop=True)

    pairs = hits[["methylation_feature", "symbol"]].drop_duplicates()
    features = pairs.groupby("methylation_feature").size().rename("count").reset_index()
    features["_number"] = pd.to_numeric(features["methylation_feature"].str.extract(r"(\d+)$")[0])
    features = features.sort_values(["count", "_number", "methylation_feature"],
                                    ascending=[False, True, True]).drop(columns="_number")
    return counts, points, features.reset_index(drop=True)


def check_historical_results(hits: pd.DataFrame, tables: tuple[pd.DataFrame, ...],
                             references: tuple[pd.DataFrame, ...]) -> None:
    counts, points, features = tables
    old_counts, old_points, old_features = references
    for current, old, keys, value in [
        (counts, old_counts, ["region", "set_type"], "count"),
        (features, old_features.rename(columns={"motif_number": "count"}),
         ["methylation_feature"], "count"),
    ]:
        require_columns(old, set(keys + [value]), "Historical count reference")
        a = current.set_index(keys)[value].sort_index()
        b = old.set_index(keys)[value].sort_index()
        if b.index.duplicated().any() or not a.equals(b):
            raise ValueError(f"Historical result mismatch for {keys}")
    point_columns = ["symbol", "set_type", "region"]
    require_columns(old_points, set(point_columns), "Historical point reference")
    if (len(points) != len(old_points)
            or set(map(tuple, points[point_columns].values)) != set(map(tuple, old_points[point_columns].values))):
        raise ValueError("Historical point/region assignments do not match")
    require_columns(old_features, {"methylation_feature", "symbols"}, "Historical pair reference")
    expected_pairs = {(row.methylation_feature, symbol)
                      for row in old_features.itertuples(index=False)
                      for symbol in str(row.symbols).split(";")}
    actual_pairs = set(map(tuple, hits[["methylation_feature", "symbol"]].values))
    if actual_pairs != expected_pairs:
        raise ValueError("Historical feature–TF pairs do not match")


def load_historical(base_dir: Path) -> tuple[pd.DataFrame, tuple[pd.DataFrame, ...], dict]:
    match_dir = base_dir / "AT2G39730_motif_gene_outputs"
    old_dir = base_dir / "TF_CHH"
    paths = {
        "negative_unmerged": match_dir / "total_TPM_negative_AT2G39730_motif_gene_matches.in_window_or_in_motif.tsv",
        "positive": match_dir / "total_TPM_positive_AT2G39730_motif_gene_matches.tsv",
        "historical_counts": old_dir / "AT2G39730_CHH_motif_region_counts.tsv",
        "historical_points_and_symbols": old_dir / "AT2G39730_CHH_symbol_region_coexp_points.tsv",
        "historical_features": old_dir / "AT2G39730_CHH_methylation_feature_motifs.tsv",
    }
    sources = {key: pd.read_csv(path, sep="\t") for key, path in paths.items()}
    for source, expected_set in [("negative_unmerged", "negative"), ("positive", "positive")]:
        require_columns(sources[source], REQUIRED, source)
        if not sources[source]["set_type"].eq(expected_set).all():
            raise ValueError(f"{source}: unexpected set_type")
    raw = pd.concat([sources["negative_unmerged"], sources["positive"]], ignore_index=True)
    old_points = sources["historical_points_and_symbols"]
    hits = prepare_hits(raw, old_points)
    tables = build_tables(hits, old_points["symbol"].drop_duplicates().tolist())
    check_historical_results(hits, tables, (sources["historical_counts"], old_points,
                                          sources["historical_features"]))
    provenance = {
        "analysis": "Historical unmerged TF motifs; no AGL/AG/AP3 collapsing",
        "sources": {key: {"sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                           "rows": len(sources[key])} for key, path in paths.items()},
        "input_rows": len(raw), "CHH_window_or_motif_rows": len(hits),
        "excluded_rows": len(raw) - len(hits),
        "filter": "context == CHH and (in_window or in_motif)",
        "symbol_mapping": "motif_id to historical symbol; exact unique mapping required",
        "historical_checks": "PASS: region counts, point assignments, feature counts and all feature–TF pairs",
        "motif_key": MOTIF_KEY,
        "region_priority": "in_motif before flanking at motif and symbol/set levels",
        "feature_count_unit": "unique methylation_feature + symbol pairs",
        "uncertainty": "Not applicable: descriptive counts, not replicate estimates",
        "statistics": "No hypothesis tests or causal/ enrichment claims",
        "versions": {"python": sys.version.split()[0], "pandas": pd.__version__,
                     "matplotlib": matplotlib.__version__, "numpy": np.__version__},
    }
    return hits, tables, provenance


def draw_figure(tables: tuple[pd.DataFrame, ...], prefix: Path, qa_tools: Path) -> None:
    sys.path.insert(0, str(qa_tools))
    from audit_panel_alignment import (
        require_matplotlib_panel_alignment, matplotlib_layout_manifest,
        audit_layout_manifest, write_json_report,
    )

    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
        "font.size": 7, "axes.labelsize": 8, "xtick.labelsize": 7,
        "ytick.labelsize": 7, "legend.fontsize": 7, "axes.linewidth": 0.65,
        "axes.spines.top": False, "axes.spines.right": False,
        "pdf.fonttype": 42, "svg.fonttype": "none", "savefig.dpi": 600,
    })
    counts, points, features = tables
    # 183 x 145 mm expressed in inches, with no export-time cropping.
    fig = plt.figure(figsize=(7.204724409448819, 5.708661417322835), facecolor="white")
    grid = fig.add_gridspec(2, 2, width_ratios=[1, 2.6], height_ratios=[1, 1.02],
                           left=0.10, right=0.975, bottom=0.155, top=0.88,
                           wspace=0.36, hspace=0.88)
    axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1]),
            fig.add_subplot(grid[1, :])]
    label_artists = []
    for ax, letter in zip(axes, "abc"):
        ax.set_gid(letter)
        label_artists.append(ax.text(0, 1, f"({letter})", ha="left", va="bottom",
            fontweight="bold", fontsize=9,
            transform=ax.transAxes + ScaledTranslation(0, 7 / 72, fig.dpi_scale_trans)))
        ax.set_axisbelow(True)
        ax.tick_params(direction="out", width=0.65, length=3)

    ax = axes[0]
    bottoms = np.zeros(2, dtype=int)
    for region in REGION_ORDER:
        values = counts[counts["region"].eq(region)].set_index("set_type").loc[SET_ORDER, "count"].to_numpy()
        ax.bar(np.arange(2), values, bottom=bottoms, width=0.48,
               color=COLORS[region], edgecolor="none")
        bottoms += values
    for x, total in enumerate(bottoms):
        ax.annotate(str(total), (x, total), xytext=(0, 4), textcoords="offset points",
                    ha="center", va="bottom", fontsize=7)
    ax.set_xticks([0, 1], [SET_LABELS[s] for s in SET_ORDER])
    ax.set_xlim(-0.65, 1.65)
    ax.set_ylim(0, max(bottoms) * 1.17)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
    ax.set_ylabel("Motif count")

    ax = axes[1]
    symbols = points["symbol"].drop_duplicates().tolist()
    x_positions = {symbol: i for i, symbol in enumerate(symbols)}
    for region, marker in [("in", "D"), ("flanking", "o")]:
        subset = points.loc[points["region"].eq(region)]
        ax.scatter(subset["symbol"].map(x_positions), subset["set_type"].map({"negative": 0, "positive": 1}),
                   s=20, marker=marker, facecolor=COLORS[region], edgecolor="white", linewidth=0.35, zorder=3)
    ax.set_xticks(range(len(symbols)), symbols, rotation=90, ha="right", va="center", rotation_mode="anchor")
    ax.set_yticks([0, 1], [SET_LABELS[s] for s in SET_ORDER])
    ax.set_ylim(1.45, -0.45)
    ax.set_xlim(-0.7, len(symbols) - 0.3)
    ax.set_xlabel("TF", labelpad=4)
    ax.grid(axis="both", color="#E5E8EA", linewidth=0.45)

    ax = axes[2]
    bars = ax.bar(range(len(features)), features["count"], width=0.65, color="#64869B", edgecolor="none")
    for bar, count in zip(bars, features["count"]):
        ax.annotate(str(count), (bar.get_x() + bar.get_width() / 2, count),
                    xytext=(0, 4), textcoords="offset points", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(range(len(features)), features["methylation_feature"],
                  rotation=70, ha="right", rotation_mode="anchor")
    ax.set_ylabel("TF count")
    ax.set_xlabel("CHH methylation feature")
    ax.set_ylim(0, features["count"].max() * 1.18)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
    ax.grid(axis="y", color="#E5E8EA", linewidth=0.45)
    fig.legend(handles=[Patch(facecolor=COLORS[r], label=label)
                        for r, label in [("in", "In motif"), ("flanking", "Flanking")]],
               loc="upper center", bbox_to_anchor=(0.54, 0.975), ncol=2, frameon=False,
               handlelength=1.4, columnspacing=1.8)

    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig, strict=True, tolerance_pt=1.5,
                                        gutter_tolerance_pt=1.5, panel_ids=list("abc"))
    # The general detector recognizes bare letters; measure parenthesized labels explicitly.
    manifest = matplotlib_layout_manifest(fig, panel_ids=list("abc"))
    for panel, artist in zip(manifest["panels"], label_artists):
        display_xy = artist.get_transform().transform(artist.get_position())
        panel["panel_label"] = panel["id"]
        panel["panel_label_anchor_pt"] = (fig.dpi_scale_trans.inverted().transform(display_xy) * 72).tolist()
    report = audit_layout_manifest(manifest, tolerance_pt=1.5, gutter_tolerance_pt=1.5,
                                   require_panel_labels=True)
    write_json_report(report, f"{prefix}.alignment.json")
    if report["summary"]["fail"] or report["summary"]["warn"]:
        plt.close(fig)
        raise ValueError("Panel alignment failed; see alignment JSON")
    # Retain the exact physical canvas; tight cropping would change its dimensions.
    fig.savefig(f"{prefix}.pdf")
    fig.savefig(f"{prefix}.svg")
    fig.savefig(f"{prefix}.png", dpi=600)
    fig.savefig(f"{prefix}.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def run_pdf_qa(prefix: Path, qa_tools: Path) -> None:
    commands = [
        ([sys.executable, str(qa_tools / "audit_pdf_text.py"), f"{prefix}.pdf",
          "--min-pt", "5", "--json"], f"{prefix}.pdf-text-audit.json"),
        ([sys.executable, str(qa_tools / "audit_figure_collisions.py"), f"{prefix}.pdf",
          "--json-out", f"{prefix}.collision-audit.json", "--strict"], None),
    ]
    for command, capture_path in commands:
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        if capture_path:
            Path(capture_path).write_text(result.stdout, encoding="utf-8")
        if result.returncode:
            raise RuntimeError(f"PDF QA failed or NOT AUDITABLE:\n{result.stdout}\n{result.stderr}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-dir", type=Path, default=DEFAULT_BASE)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--output-prefix", default=DEFAULT_PREFIX)
    parser.add_argument("--qa-tools", type=Path, default=DEFAULT_QA_TOOLS)
    args = parser.parse_args()
    if Path(args.output_prefix).name != args.output_prefix or args.output_prefix in {"", ".", ".."}:
        parser.error("--output-prefix must be a filename stem")
    for script in ["audit_panel_alignment.py", "audit_pdf_text.py", "audit_figure_collisions.py"]:
        if not (args.qa_tools / script).is_file():
            parser.error(f"QA tool missing: {script}; specify --qa-tools")
    hits, tables, provenance = load_historical(args.base_dir)
    output_dir = args.output_dir or args.base_dir / "TF_CHH"
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = output_dir / args.output_prefix
    draw_figure(tables, prefix, args.qa_tools)
    run_pdf_qa(prefix, args.qa_tools)
    for table, name in zip(tables, ["panel-a-counts", "panel-b-points", "panel-c-features"]):
        table.to_csv(f"{prefix}.{name}.tsv", sep="\t", index=False)
    hits.to_csv(f"{prefix}.source-hits.tsv", sep="\t", index=False)
    provenance["automated_qa"] = "PASS: alignment, PDF text >= 5 pt, strict collision audit"
    provenance["visual_review"] = "Separate final-size review required; not inferred from automated checks"
    Path(f"{prefix}.provenance.json").write_text(json.dumps(provenance, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    caption = (
        "Historical unmerged TF–CHH associations at the RCA locus.\n\n"
        "(a) Unique motif occurrences in non-co-exp and co-exp groups, partitioned into in-motif "
        "and flanking associations. Bar labels give total occurrences; a motif with any in-motif "
        "CHH match is assigned to In motif. (b) TF symbols by co-expression group; each symbol "
        "appears once per group, with in-motif matches taking priority. Diamonds mark In motif; "
        "circles mark Flanking. (c) Number of distinct TF symbols associated with each CHH "
        "methylation feature, counting each feature–symbol pair once.\n\n"
        "The reconstruction retains separate AGL/AG/AP3 entries. Negative and positive labels "
        "denote non-co-expression and co-expression, not effect direction. These are descriptive "
        "counts, without replicate-based error bars or hypothesis tests. Differences in counts "
        "alone do not establish enrichment or regulation.\n"
    )
    Path(f"{prefix}.legend.txt").write_text(caption, encoding="utf-8")
    print(tables[0].to_string(index=False))
    print(f"TF points: {len(tables[1])}; CHH features: {len(tables[2])}; feature–TF pairs: {tables[2]['count'].sum()}")
    print(f"Saved: {prefix}")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, RuntimeError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
