#!/usr/bin/env python3
"""RCA total-expression SHAP methylation / predicted TF motif spatial maps.

All displayed genomic positions are 1-based inclusive. Raw minus-strand FIMO
sequences are reversed once to follow increasing genomic position.
Only sites inside motif intervals, including both endpoints, are matched.
These descriptive maps do not establish TF binding or causal regulation.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rca-motif-mpl"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Rectangle

PROJECT = Path(__file__).resolve().parents[2]
DEFAULT_WORKDIR = PROJECT / "list"
DEFAULT_GFF = PROJECT.parent / "ref/Arabidopsis_thaliana/refgen/GCF_000001735.4_TAIR10.1_genomic.gff"
QA_SCRIPTS = Path("/datapool/home/2023102768/.codex/skills/nature-figure/scripts")
DEFAULT_DPI = 600
SELECT_TOP_FRACTION = .50
FLANK_BP = 0
TRANSCRIPT = "rna-NM_129531.3"
BLUE = "#2997FF"
BASE_COLORS = {"A": "#3182BD", "T": "#D95F0E", "C": "#31A354", "G": "#756BB1"}
CONTEXT_COLORS = {"CG": "#6A3D9A", "CHG": "#1B9E77", "CHH": "#3182BD"}
TF_IDS = {"LHY": "AT1G01060", "AGL22": "AT2G22540", "LCL5": "AT3G09600",
          "LCL1": "AT5G02840", "RVE1": "AT5G17300"}


def read_table(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def normalize_chrom(chrom):
    accession = {"NC_003070": "Chr1", "NC_003071": "Chr2", "NC_003074": "Chr3", "NC_003075": "Chr4", "NC_003076": "Chr5"}
    if chrom.split(".")[0] in accession:
        return accession[chrom.split(".")[0]]
    if re.fullmatch(r"(?:chr)?[1-5]", chrom, re.I):
        return "Chr" + chrom[-1]
    raise ValueError(f"Unrecognized chromosome: {chrom}")


def parse_shap(path):
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    headings = [i for i, line in enumerate(lines) if line.strip() == "所有特征及其平均|SHAP|值:"]
    if len(headings) != 1 or not any("total_TPM 导入模型结果" in x for x in lines):
        raise ValueError("Expected exactly one total_TPM full SHAP section")
    rows, seen = [], set()
    for line in lines[headings[0]+1:]:
        if not line.strip():
            if rows:
                break
            continue
        match = re.fullmatch(r"\s+((?:snp|CG|CHG|CHH)_\d+):\s*(\S+)\s*", line)
        if not match:
            if rows and not line.startswith("    "):
                break
            raise ValueError(f"Malformed full SHAP row: {line}")
        feature, raw = match.groups()
        value = float(raw)
        if feature in seen or not math.isfinite(value) or value < 0:
            raise ValueError(f"Duplicate or invalid SHAP value: {feature}")
        seen.add(feature)
        rows.append(dict(feature=feature, context=feature.split("_")[0].upper(), mean_abs_shap=value))
    rows.sort(key=lambda r: (-r["mean_abs_shap"], r["context"], r["feature"]))
    n50 = math.ceil(len(rows)*SELECT_TOP_FRACTION)
    for rank, row in enumerate(rows, 1):
        row.update(rank=rank, ranking_denominator=len(rows),
                   is_top50=rank <= n50, selected=rank <= n50 and row["context"] in CONTEXT_COLORS)
    return rows


def load_sites(workdir, ranking):
    with (workdir / "RCA/Alt.all.snp_meth.filtered.tsv").open(encoding="utf-8-sig") as handle:
        header = next(csv.reader(handle, delimiter="\t"))
    if len(header) != len(set(header)):
        raise ValueError("Duplicate matrix columns")
    selected = [dict(r) for r in ranking if r["selected"]]
    missing = {r["feature"] for r in selected} - set(header)
    if missing:
        raise ValueError(f"Selected methylation columns absent from matrix: {missing}")
    lookup = {}
    for context in CONTEXT_COLORS:
        with (workdir / f"meth/{context}.CV.before.bedGraph").open() as handle:
            for index, line in enumerate(handle, 1):
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 5:
                    raise ValueError(f"Expected five-column headerless bedGraph: {context}, row {index}")
                chrom, start, end = normalize_chrom(fields[0]), int(fields[1]), int(fields[2])
                if start < 0 or end != start+1:
                    raise ValueError(f"Not a single-base BED interval: {context}, row {index}")
                lookup[f"{context}_{index}"] = dict(chrom=chrom, bed_start=start, bed_end=end, position_1based=end)
    for site in selected:
        if site["feature"] not in lookup:
            raise ValueError(f"Numbered site absent: {site['feature']}")
        site.update(lookup[site["feature"]])
    return selected


def normalize_motif(row):
    result = dict(row)
    result.update(chrom=normalize_chrom(row["chrom"]), motif_start=int(row["motif_start"]), motif_end=int(row["motif_end"]),
                  matched_sequence=row["matched_sequence"].upper())
    if result["strand"] not in ("+", "-") or result["motif_start"] < 1:
        raise ValueError(f"Invalid motif coordinate/strand: {row}")
    if len(result["matched_sequence"]) != result["motif_end"]-result["motif_start"]+1 or set(result["matched_sequence"]) - set(BASE_COLORS):
        raise ValueError(f"Invalid motif sequence: {row}")
    return result




def load_motifs(workdir):
    annotations = {r["Gene_ID"]: r for r in read_table(workdir / "coex2/total/total.gene.TF_ids.annotated.txt")}
    groups = {name: [] for name in TF_IDS}
    id_to_name = {gene: name for name, gene in TF_IDS.items()}
    seen = set()
    with (workdir / "coex2/all.fimo.out.filt").open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip().split("\t")
            if len(f) != 10:
                raise ValueError("Expected ten-column headerless FIMO table")
            if f[0] not in id_to_name:
                continue
            name = id_to_name[f[0]]
            raw_sequence = f[9].upper()
            row = normalize_motif(dict(motif_id=f[0], motif_alt_id=f[1], symbol=name,
                chrom=f[2], motif_start=f[3], motif_end=f[4], strand=f[5],
                matched_sequence=raw_sequence if f[5] == "+" else raw_sequence[::-1]))
            key = tuple(row[k] for k in ("motif_id", "motif_alt_id", "chrom", "motif_start", "motif_end", "strand", "matched_sequence"))
            if key in seen:
                continue
            seen.add(key)
            annotation = annotations.get(f[0])
            if annotation is None or annotation["Symbol"] != name:
                raise ValueError(f"Missing or mismatched total-expression TF annotation: {name}")
            row.update(group="coexpressed", family=annotation["Family"], score=float(f[6]),
                       p_value=float(f[7]), q_value=float(f[8]), raw_fimo_sequence=raw_sequence,
                       sequence_orientation="increasing_genomic_coordinates")
            groups[name].append(row)
    if any(len(rows) != 1 for rows in groups.values()):
        raise ValueError(f"Expected one FIMO locus for each requested TF: { {k: len(v) for k, v in groups.items()} }")
    return groups


def load_gene(gff):
    model = dict(exons=[], cds=[])
    with gff.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip().split("\t")
            if len(f) != 9:
                continue
            attrs = dict(v.split("=", 1) for v in f[8].split(";") if "=" in v)
            if attrs.get("locus_tag") != "AT2G39730":
                continue
            if f[2] == "mRNA" and attrs.get("ID") == TRANSCRIPT:
                model.update(chrom=normalize_chrom(f[0]), start=int(f[3]), end=int(f[4]), strand=f[6])
            if attrs.get("Parent") == TRANSCRIPT and f[2] in ("exon", "CDS"):
                model["exons" if f[2] == "exon" else "cds"].append((int(f[3]), int(f[4])))
    model["exons"].sort()
    model["cds"].sort()
    if (model.get("start"), model.get("end"), model.get("strand"), len(model["exons"])) != (16570746, 16573692, "-", 7):
        raise ValueError("RCA transcript differs from the seven-exon NM_129531.3 model")
    return model


def match_sites(record, sites):
    start, end = record["motif_start"], record["motif_end"]
    result = []
    for site in sites:
        pos = site["position_1based"]
        if site["chrom"] != record["chrom"] or not start <= pos <= end:
            continue
        result.append(dict(site, motif_id=record["motif_id"], tf_id=record["motif_id"], symbol=record["symbol"],
                           family=record["family"], group=record["group"], motif_start=start, motif_end=end,
                           strand=record["strand"], distance_bp=0, relation="motif", side="inside"))
    return sorted(result, key=lambda r: (r["position_1based"], r["feature"]))


def load_data(workdir, gff):
    ranking = parse_shap(workdir / "xgboot/TPM_4.5_all/total_TPM_import_results.txt")
    sites = load_sites(workdir, ranking)
    groups = load_motifs(workdir)
    matches = [m for rows in groups.values() for r in rows for m in match_sites(r, sites)]
    plotted = {r["feature"] for r in matches}
    for site in sites:
        site["plotted"] = site["feature"] in plotted
        site["exclusion_reason"] = "" if site["plotted"] else "outside_requested_five_TF_motif_intervals"
    counts = dict(total_ranked_features=len(ranking), top50_rank_count=math.ceil(len(ranking)*SELECT_TOP_FRACTION),
                  selected_methylation_sites=len(sites), plotted_unique_sites=len(plotted),
                  context_counts={c: sum(s["context"] == c for s in sites) for c in CONTEXT_COLORS},
                  unplotted_sites=[s["feature"] for s in sites if not s["plotted"]], motif_site_pairs=len(matches),
                  tf_count=len(groups), tf_order=list(groups),
                  motif_site_counts={name: len(match_sites(rows[0], sites)) for name, rows in groups.items()},
                  flank_bp=FLANK_BP)
    if not sites:
        raise ValueError("Empty SHAP methylation selection")
    return ranking, sites, groups, matches, load_gene(gff), counts




def site_color(site):
    return CONTEXT_COLORS[site["context"]]


def draw_motif(ax, record, sites, compact=False):
    points = match_sites(record, sites)
    start, end = record["motif_start"], record["motif_end"]
    ax.set(xlim=(start-1.2, end+1.2), ylim=(0, 1), xticks=[], yticks=[])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#777777")
        spine.set_linewidth(.6)
    ax.add_patch(Rectangle((start-.5, .50), end-start+1, .25, facecolor="#F1EFF6", edgecolor="none"))
    ax.axhline(.50, color="#777777", lw=.6)
    ax.text(.02, .91, record["symbol"], transform=ax.transAxes, fontsize=7, weight="bold", va="center")
    label = f"({record['strand']}) {record['chrom']}:{start}–{end}"
    ax.text(.98, .81 if compact else .91, label, transform=ax.transAxes, ha="right", va="center",
            fontsize=5.2 if compact else 5.5, color="#555555")
    for position, base in enumerate(record["matched_sequence"], start):
        ax.text(position, .625, base, color=BASE_COLORS[base], fontsize=8, weight="bold", ha="center", va="center")
    left, right = ("5′", "3′") if record["strand"] == "+" else ("3′", "5′")
    ax.text(.02, .625, left, transform=ax.transAxes, fontsize=5.5, ha="center", va="center", color="#555555")
    ax.text(.98, .625, right, transform=ax.transAxes, fontsize=5.5, ha="center", va="center", color="#555555")
    ax.figure.canvas.draw()
    renderer = ax.figure.canvas.get_renderer()
    lanes = [[], [], []]
    levels = [.23, .135, .065]
    for site in points:
        x = site["position_1based"]
        ax.plot(x, .37, marker="*", ms=7, markeredgewidth=0, color=site_color(site), zorder=4)
        text = ax.text(x, levels[0], site["feature"], fontsize=6, ha="center", va="center", color="#333333")
        bbox = text.get_window_extent(renderer)
        px = ax.transData.transform((x, 0))[0]
        half = bbox.width/2
        px_label = min(max(px, ax.bbox.x0+half+3), ax.bbox.x1-half-3)
        label_x = ax.transData.inverted().transform((px_label, 0))[0]
        interval = (px_label-half-2, px_label+half+2)
        lane = next((i for i, used in enumerate(lanes) if all(interval[1] <= a or interval[0] >= b for a, b in used)), None)
        if lane is None:
            raise ValueError(f"Horizontal site IDs do not fit in {record['symbol']}")
        lanes[lane].append(interval)
        y = levels[lane]
        text.set_position((label_x, y))
        if lane or abs(px_label-px) > 1:
            top = y+bbox.height/ax.bbox.height/2+.025
            ax.plot([x, label_x, label_x], [.315, top+.02, top], lw=.4, color="#999999", zorder=1)
    return points


def draw_gene(ax, gene, records):
    right = max(r["motif_end"] for r in records)+150
    ax.set(xlim=(16570580, right), ylim=(0, 1))
    ax.axis("off")
    height_pt = ax.figure.get_figheight()*72*ax.get_position().height
    box_h = 8/height_pt
    y = .62
    ax.hlines(y, gene["end"]+.5, right-50, color=BLUE, lw=.7)
    for start, end in gene["exons"]:
        ax.add_patch(Rectangle((start-.5, y-box_h/2), end-start+1, box_h, facecolor="white", edgecolor=BLUE, lw=.8))
    for start, end in gene["cds"]:
        ax.add_patch(Rectangle((start-.5, y-box_h/2), end-start+1, box_h, facecolor=BLUE, edgecolor="none"))
    for (_, end), (start, _) in zip(gene["exons"], gene["exons"][1:]):
        ax.hlines(y, end+.5, start-.5, color=BLUE, lw=.7)
    ax.text(gene["start"], .08, "RCA", fontsize=8, style="italic", ha="left", va="center", color="#333333")
    ax.text((gene["start"]+gene["end"])/2+170, .08, "Chr2:16573692–16570746  (−)",
            fontsize=6.5, ha="center", va="center")
    anchors = {}
    union_center = (min(r["motif_start"] for r in records)+max(r["motif_end"] for r in records))/2
    plotted_centers = set()
    for record in records:
        x = union_center if len(records) > 1 else (record["motif_start"]+record["motif_end"])/2
        if x not in plotted_centers:
            if len(records) > 1:
                width_pt = ax.figure.get_figwidth()*72*ax.get_position().width
                marker_width = 8*(right-16570580)/width_pt
                marker_height = 8/height_pt
                ax.add_patch(Rectangle((x-marker_width/2, y-marker_height/2), marker_width, marker_height,
                                       facecolor="#FFB000", edgecolor="none", zorder=4))
            else:
                ax.plot(x, y, marker="s", ms=5, color="#FFB000", markeredgewidth=0, zorder=4)
            plotted_centers.add(x)
        anchors[record["symbol"]] = (x, y)
    return anchors


def fig_point(fig, ax, point):
    return fig.transFigure.inverted().transform(ax.transData.transform(point))


def connector(fig, points, arrow=True):
    xs, ys = zip(*points)
    fig.add_artist(Line2D(xs[:-1] if arrow else xs, ys[:-1] if arrow else ys, transform=fig.transFigure,
                         color="#555555", lw=.65))
    if arrow:
        fig.add_artist(FancyArrowPatch(points[-2], points[-1], transform=fig.transFigure, arrowstyle="->",
                                      mutation_scale=6, lw=.65, color="#555555", shrinkA=0, shrinkB=0))




def add_legend(fig):
    handles = [Line2D([], [], marker="*", color=color, markersize=7, markeredgewidth=0, lw=0)
               for color in CONTEXT_COLORS.values()]
    fig.legend(handles, list(CONTEXT_COLORS), ncol=3, loc="upper right", bbox_to_anchor=(.955, .995), frameon=False,
               prop={"size": 7, "weight": "bold"}, handletextpad=.7, columnspacing=2.5)
    multi = fig.get_figheight()*25.4 > 100
    arrow_y = .575 if multi else .82
    arrow_start, arrow_end = (.80, .65) if multi else (.93, .78)
    fig.add_artist(FancyArrowPatch((arrow_start, arrow_y), (arrow_end, arrow_y), transform=fig.transFigure,
                                  arrowstyle="-|>", mutation_scale=9, lw=.9, color="#333333"))
    fig.text(.5, .016, "SHAP rank: top 50% of all model features. Motif overlap only; no flanking sites.",
             ha="center", va="center", fontsize=5.5, color="#555555")




def render(family, records, sites, gene, out, dpi):
    sys.path.insert(0, str(QA_SCRIPTS))
    from audit_panel_alignment import require_matplotlib_panel_alignment
    multi = len(records) > 1
    width_mm = 183
    height = 145 if multi else 80
    fig = plt.figure(figsize=(width_mm/25.4, height/25.4))
    add_legend(fig)
    gene_ax = fig.add_axes([.045, .475, .91, .105] if multi else [.045, .65, .91, .15], label="gene_overview")
    anchors = draw_gene(gene_ax, gene, records)
    axes = []
    for index, record in enumerate(records):
        rect = [.08 if index % 2 == 0 else .54, .64 if index < 2 else .15, .38, .20] if multi else [.12, .17, .76, .32]
        ax = fig.add_axes(rect, label=record["symbol"])
        draw_motif(ax, record, sites, compact=multi)
        axes.append(ax)
    fig.canvas.draw()
    if not multi:
        x, y = fig_point(fig, gene_ax, anchors[records[0]["symbol"]])
        connector(fig, [(x, y-3/(height/25.4*72)), (x, .56), (.5, .56), (.5, .495)])
    else:
        # Independent paths leave distinct points on the same locus glyph;
        # the marker represents the overlapping intervals, not shifted loci.
        x, y = fig_point(fig, gene_ax, anchors[records[0]["symbol"]])
        dx = 2/(fig.get_figwidth()*72)
        dy = 4/(fig.get_figheight()*72)
        routes = {}
        for index, ax in enumerate(axes):
            box = ax.get_position()
            center = (box.x0+box.x1)/2
            source_x = x-dx if index % 2 == 0 else x+dx
            upper = index < 2
            elbow = [.60, .62, .44, .42][index]
            target_y = box.y0-.004 if upper else box.y1+.004
            routes[records[index]["symbol"]] = [(source_x, y+dy if upper else y-dy),
                                                (source_x, elbow), (center, elbow), (center, target_y)]
        audit_routes(routes, out / "RCA_MYB_related_motif_colocalization.connector-audit.json")
        for route in routes.values():
            connector(fig, route)
    prefix = out / f"RCA_{family}_motif_colocalization"
    names = [r["symbol"] for r in records]
    require_matplotlib_panel_alignment(fig, json_out=f"{prefix}.alignment.json", exclude_axes=[gene_ax],
                                       row_groups=[names[:2], names[2:]] if multi else [],
                                       column_groups=[names[::2], names[1::2]] if multi else [],
                                       tolerance_pt=1.5, gutter_tolerance_pt=1.5, strict=True)
    fig.savefig(f"{prefix}.pdf", facecolor="white")
    fig.savefig(f"{prefix}.svg", facecolor="white")
    fig.savefig(f"{prefix}.png", dpi=dpi, facecolor="white")
    fig.savefig(f"{prefix}.tiff", dpi=dpi, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    return prefix


def audit_routes(routes, output):
    findings = []
    names = list(routes)
    for i, name in enumerate(names):
        for other in names[i+1:]:
            for a, b in zip(routes[name], routes[name][1:]):
                for c, d in zip(routes[other], routes[other][1:]):
                    dx = min(max(a[0], b[0]), max(c[0], d[0]))-max(min(a[0], b[0]), min(c[0], d[0]))
                    dy = min(max(a[1], b[1]), max(c[1], d[1]))-max(min(a[1], b[1]), min(c[1], d[1]))
                    if dx >= -1e-9 and dy >= -1e-9:
                        findings.append(dict(first=name, second=other, segments=[[a, b], [c, d]]))
    report = dict(verdict="FAIL" if findings else "PASS", route_count=len(routes),
                  crossing_or_overlap_count=len(findings), findings=findings, routes=routes)
    output.write_text(json.dumps(report, indent=2)+"\n")
    if findings:
        raise ValueError(f"Connector paths intersect: {findings}")


def write_tsv(path, rows, fieldnames=None):
    columns = fieldnames or list(dict.fromkeys(key for row in rows for key in row))
    if not columns:
        raise ValueError(f"Explicit headers required for an empty table: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_colocalization_tables(out, groups, matches):
    fields = ["context", "feature", "chrom", "position_1based", "tf_id", "symbol", "family",
              "motif_start", "motif_end", "strand", "rank", "mean_abs_shap", "relation"]
    write_tsv(out / "RCA_motif_site_matches.tsv", matches, fields)
    summary, tf_summary = [], []
    for context in CONTEXT_COLORS:
        subset = [r for r in matches if r["context"] == context]
        symbols = [name for name in groups if any(r["symbol"] == name for r in subset)]
        ids = [groups[name][0]["motif_id"] for name in symbols]
        site_ids = list(dict.fromkeys(r["feature"] for r in sorted(subset, key=lambda r: r["position_1based"])))
        summary.append(dict(context=context, tf_count=len(symbols), tf_symbols=";".join(symbols), tf_ids=";".join(ids),
                            unique_site_count=len(site_ids), site_ids=";".join(site_ids), motif_site_pair_count=len(subset)))
        write_tsv(out / f"RCA_{context}_motif_colocalization.tsv", subset, fields)
    for name, records in groups.items():
        record = records[0]
        for context in CONTEXT_COLORS:
            subset = [r for r in matches if r["symbol"] == name and r["context"] == context]
            tf_summary.append(dict(tf_id=record["motif_id"], symbol=name, family=record["family"], context=context,
                                   site_count=len({r["feature"] for r in subset}),
                                   site_ids=";".join(r["feature"] for r in subset), chrom=record["chrom"],
                                   motif_start=record["motif_start"], motif_end=record["motif_end"], strand=record["strand"]))
    write_tsv(out / "RCA_context_TF_colocalization_summary.tsv", summary)
    write_tsv(out / "RCA_TF_context_colocalization_counts.tsv", tf_summary)


def run_qa(prefix):
    commands = [
        ("pdf-text", [sys.executable, str(QA_SCRIPTS / "audit_pdf_text.py"), f"{prefix}.pdf", "--min-pt", "5", "--json"]),
        ("collision", [sys.executable, str(QA_SCRIPTS / "audit_figure_collisions.py"), f"{prefix}.pdf", "--strict",
                       "--json-out", f"{prefix}.collision-audit.json", "--overlay-pdf", f"{prefix}.collision-audit.pdf"]),
    ]
    for name, cmd in commands:
        result = subprocess.run(cmd, capture_output=True, text=True)
        Path(f"{prefix}.{name}.txt").write_text(result.stdout+result.stderr, encoding="utf-8")
        if result.returncode:
            raise RuntimeError(f"{prefix.name} {name} QA failed: {result.stdout[-4000:]} {result.stderr[-500:]}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, default=DEFAULT_WORKDIR)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--gff", type=Path, default=DEFAULT_GFF)
    parser.add_argument("--dpi", type=int, default=DEFAULT_DPI)
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    if args.dpi < 300:
        parser.error("Publication raster exports require --dpi >= 300")
    ranking, sites, groups, matches, gene, counts = load_data(args.workdir, args.gff)
    if args.validate_only:
        print(json.dumps(dict(counts=counts, gene=gene), indent=2))
        return
    out = args.output_dir or args.workdir / "meth/motif"
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"], "font.size": 7,
                         "pdf.fonttype": 42, "svg.fonttype": "none", "axes.linewidth": .6})
    write_tsv(out / "RCA_all_feature_SHAP_ranking.tsv", ranking)
    write_tsv(out / "RCA_selected_methylation_sites.tsv", sites)
    write_colocalization_tables(out, groups, matches)
    write_tsv(out / "RCA_plot_motifs.tsv", [r for rows in groups.values() for r in rows])
    write_tsv(out / "RCA_gene_model.tsv", [dict(transcript="NM_129531.3", chrom=gene["chrom"], feature=kind,
              start=s, end=e, strand="-") for kind in ("exons", "cds") for s, e in gene[kind]])
    (out / "selection_counts.json").write_text(json.dumps(counts, indent=2)+"\n")
    family_groups = {"MIKC_MADS": groups["AGL22"],
                     "MYB_related": [groups[name][0] for name in ("LHY", "LCL5", "LCL1", "RVE1")]}
    for family, records in family_groups.items():
        if any(record["family"] != family for record in records):
            raise ValueError(f"TF annotation does not agree with requested family: {family}")
        prefix = render(family, records, sites, gene, out, args.dpi)
        run_qa(prefix)
        print(f"Exported and audited: {prefix}", flush=True)
    (out / "README.md").write_text(
        "# RCA motif-internal methylation colocalization by TF family\n\n"
        "Exactly two current figures are generated: MIKC_MADS (AGL22) and MYB_related (LHY, LCL5, LCL1, RVE1), "
        "with the latter keeping four separate motif panels. TF family annotations are validated against the "
        "total-expression annotation table. Filtered FIMO records supply the motif coordinates and sequences.\n\n"
        "## Selection and colocalization\n\n"
        f"Rank all {len(ranking)} model features by descending mean absolute SHAP, breaking ties by context and feature name. "
        f"Select the first {counts['top50_rank_count']} features (ceil(N × 50%)), then retain the {len(sites)} methylation "
        "features in the filtered matrix. Every selected site remains in the source table, with an exclusion reason "
        "if it does not overlap one of the requested motifs. CG/CHG/CHH have uniform purple/green/blue colors; "
        "no color intensity encodes SHAP rank or methylation magnitude.\n\n"
        "Colocalization requires the same chromosome and motif_start <= position_1based <= motif_end. "
        "Both endpoints count. No upstream/downstream sites are matched or drawn. "
        f"There are {counts['plotted_unique_sites']} distinct overlapping sites and {len(matches)} TF–site pairs. "
        "Shared sites are counted once in a context's unique-site total and once per corresponding TF in pair totals. "
        "These are descriptive spatial overlaps, not proof of TF binding or causal regulation.\n\n"
        "## Tables\n\n"
        "RCA_context_TF_colocalization_summary.tsv has three rows (CG, CHG, CHH), with TF names/IDs, unique-site "
        "counts and IDs, and TF–site pair counts. RCA_TF_context_colocalization_counts.tsv has all 15 TF-by-context "
        "combinations, including explicit zero counts. RCA_CG_motif_colocalization.tsv, "
        "RCA_CHG_motif_colocalization.tsv and RCA_CHH_motif_colocalization.tsv provide the corresponding detailed "
        "pairs; an empty context retains its header. RCA_motif_site_matches.tsv combines all detailed pairs.\n\n"
        "## Coordinates and display\n\n"
        "The Nth physical line of each headerless bedGraph identifies context_N. Single-base BED intervals are "
        "0-based half-open; displayed position = start+1 = end. FIMO/GFF intervals are 1-based inclusive. "
        "NC_003071.7 is normalized to Chr2. Raw minus-strand FIMO sequences are reversed once, not complemented, "
        "to display increasing genomic coordinates; both raw and displayed sequences are retained. "
        "Plus motif panels run 5′–3′; minus panels run 3′–5′.\n\n"
        "Local motif panels now use linear genomic positions without the former flank compression. Only small "
        "blank margins are reserved for strand-end labels. Full methylation IDs are horizontal; stars remain "
        "at their exact coordinate and leaders are used only to keep text within its panel. "
        "The RCA NM_129531.3 seven-exon drawing stays at left, with solid-blue CDS and hollow UTRs. "
        "The unlabelled left-pointing negative-strand arrow is below the methylation legend in the AGL22 figure; "
        "in MYB_related it is in the gap below LHY/LCL5 and above RCA, clear of motif connectors. "
        "No family headings appear in the figures. MYB_related has LHY/LCL5 above "
        "the gene and LCL1/RVE1 below. Four independent nonintersecting connectors depart from distinct points on "
        "a single marker representing the overlapping motif intervals; no shared junctions or displaced loci are implied.\n\n"
        "MIKC_MADS is 183 × 80 mm; MYB_related is 183 × 145 mm. PDF/SVG retain editable text; PNG/TIFF default "
        "to 600 dpi. The four MYB_related panels pass a 1.5 pt alignment gate; all final PDFs undergo glyph-size "
        "and strict collision audits. Earlier outputs are retained separately under archive/.\n",
        encoding="utf-8")



if __name__ == "__main__":
    main()
