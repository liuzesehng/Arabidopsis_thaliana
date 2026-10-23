#!/usr/bin/env python3
"""RCA motif/CHH spatial maps; all displayed genomic coordinates are 1-based.

The sequence on a minus-strand hit is reversed (not complemented) to show
the actual bound strand from left to right in increasing genomic coordinates.
No inference of binding or methylation-dependent regulation is made.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, FancyArrowPatch

PROJECT = Path(__file__).resolve().parents[2]
DEFAULT_WORKDIR = PROJECT / "list"
DEFAULT_GFF = PROJECT.parent / "ref/Arabidopsis_thaliana/refgen/GCF_000001735.4_TAIR10.1_genomic.gff"
QA_SCRIPTS = Path("/datapool/home/2023102768/.agents/skills/nature-figure/scripts")
DEFAULT_DPI = 600
TRANSCRIPT = "rna-NM_129531.3"
BLUE = "#2997FF"
CHH_BLUE = "#3182BD"
BASE_COLORS = {"A": "#3182BD", "T": "#D95F0E", "C": "#31A354", "G": "#756BB1"}
WINDOW = (16574872, 16574890)
ORDER = ["LHY", "EPR1", "RVE7L", "RVE5", "LCL1", "RVE1", "RVE6"]


# Optional external QA: keep each plotting script usable on its own.
_OPTIONAL_QA_SKIPPED = []


def _skip_optional_qa(reason, output=None):
    """Record an unavailable check without claiming that it passed."""
    if reason not in _OPTIONAL_QA_SKIPPED:
        _OPTIONAL_QA_SKIPPED.append(reason)
        print(f"[QA SKIPPED] {reason}", file=sys.stderr)
    report = {"verdict": "SKIPPED", "status": "SKIPPED", "auditable": False,
              "reason": reason}
    if output is not None:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return report


def _skip_alignment(*args, json_out=None, **kwargs):
    return _skip_optional_qa("audit_panel_alignment.py is unavailable", json_out)


def _skip_missing_qa_command(command, output=None):
    """Only missing optional scripts are skipped; executed checks stay strict."""
    script = Path(command[1])
    if script.is_file():
        return False
    reason = f"Optional QA tool not found: {script}"
    _skip_optional_qa(reason, output)
    if "--json-out" in command:
        _skip_optional_qa(reason, command[command.index("--json-out") + 1])
    return True


def read_table(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_data(workdir):
    source = workdir / "meth_experiment/CHH_TF_loci"
    symbols = {}
    for row in read_table(source / "TF.non_coexpressed.annotated.txt"):
        key, symbol = row["Gene_ID"].strip(), row["Symbol"].strip()
        if symbol:
            if key in symbols and symbols[key] != symbol:
                raise ValueError(f"Conflicting symbols for {key}")
            symbols[key] = symbol
    sites = {}
    with (workdir / "meth_experiment/CHH.CV.before.numbered.tsv").open() as handle:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 6:
                raise ValueError("Expected six columns in numbered CHH table")
            site, chrom, start, end, _, _ = fields
            start, end = int(start), int(end)
            if end != start + 1 or site in sites:
                raise ValueError(f"Invalid or duplicate single-base site: {site}")
            sites[site] = (chrom, start, end)
    hits = []
    with (workdir / "coex2/all.fimo.out.filt").open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip().split("\t")
            if len(f) != 10:
                raise ValueError("Expected ten columns in FIMO table")
            hits.append(dict(tf_id=f[0], motif_id=f[1], chrom=f[2],
                             start=int(f[3]), end=int(f[4]), strand=f[5],
                             score=f[6], p_value=f[7], q_value=f[8], sequence=f[9].upper()))
    groups, counts = {}, {}
    for group in ("coexpressed", "non_coexpressed"):
        rows = read_table(source / f"CHH_TF_loci.{group}.tsv")
        selected = [r for r in rows if int(r["motif"]) == 1]
        records, seen = [], set()
        for row in selected:
            chrom, bed_start, position = sites[row["site_id"]]
            if int(row["location"]) != position:
                raise ValueError(f"Coordinate mismatch: {row['site_id']}")
            symbol = row["symble"].strip() if group == "coexpressed" else symbols[row["tf_id"]]
            matched = [h for h in hits if h["tf_id"] == row["tf_id"]
                       and h["chrom"] == chrom and h["start"] <= position <= h["end"]]
            if not matched:
                raise ValueError(f"No overlapping FIMO hit for {row['tf_id']}")
            for hit in matched:
                sequence = hit["sequence"]
                if len(sequence) != hit["end"] - hit["start"] + 1 or set(sequence) - set(BASE_COLORS):
                    raise ValueError(f"Invalid motif sequence: {hit}")
                if hit["strand"] not in ("+", "-"):
                    raise ValueError("Invalid FIMO strand")
                key = (row["site_id"], *[hit[k] for k in ("tf_id", "motif_id", "chrom", "start", "end", "strand", "sequence")])
                if key in seen:
                    continue
                seen.add(key)
                records.append(dict(group=group, symbol=symbol, site_id=row["site_id"],
                                    bed_start=bed_start, position=position, **hit,
                                    display_sequence=sequence if hit["strand"] == "+" else sequence[::-1]))
        groups[group] = sorted(records, key=lambda r: ORDER.index(r["symbol"]) if group == "non_coexpressed" else r["symbol"])
        counts[group] = dict(input_rows=len(rows), motif_1_rows=len(selected), excluded_motif_0=len(rows)-len(selected),
                             unique_tfs=len({r["tf_id"] for r in records}), plotted_hits=len(records))
    if len(groups["coexpressed"]) != 1 or len(groups["non_coexpressed"]) != 7:
        raise ValueError("This figure contract expects one coexpressed and seven non-coexpressed motifs")
    if groups["coexpressed"][0]["symbol"] != "LCL5":
        raise ValueError("Expected coexpressed LCL5")
    if {r["site_id"] for rows in groups.values() for r in rows} != {"CHH_1651"}:
        raise ValueError("Expected only CHH_1651")
    return groups, counts


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
                model.update(chrom=f[0], start=int(f[3]), end=int(f[4]), strand=f[6])
            if attrs.get("Parent") == TRANSCRIPT and f[2] in ("exon", "CDS"):
                model["exons" if f[2] == "exon" else "cds"].append((int(f[3]), int(f[4])))
    model["exons"].sort()
    model["cds"].sort()
    if (model.get("start"), model.get("end"), model.get("strand"), len(model["exons"])) != (16570746, 16573692, "-", 7):
        raise ValueError("RCA transcript does not match the approved seven-exon model")
    return model


def draw_gene(ax, gene, records):
    ax.set(xlim=(16570580, 16575100), ylim=(0, 1))
    ax.axis("off")
    # Draw only introns and flanking DNA: no baseline inside hollow UTRs.
    ax.hlines(.62, gene["end"]+.5, 16575020, color=BLUE, lw=.7)
    height_pt = ax.figure.get_figheight()*72*ax.get_position().height
    box_height = 8 / height_pt
    box_bottom = .62-box_height/2
    for start, end in gene["exons"]:
        ax.add_patch(Rectangle((start-.5, box_bottom), end-start+1, box_height, facecolor="white", edgecolor=BLUE, lw=.8))
    for start, end in gene["cds"]:
        ax.add_patch(Rectangle((start-.5, box_bottom), end-start+1, box_height, facecolor=BLUE, edgecolor="none"))
    for (_, end), (start, _) in zip(gene["exons"], gene["exons"][1:]):
        ax.hlines(.62, end+.5, start-.5, color=BLUE, lw=.7)
    lo, hi = min(r["start"] for r in records), max(r["end"] for r in records)
    # A compact 5 pt square marks the locus; local coordinates remain exact.
    width_pt = ax.figure.get_figwidth()*72*ax.get_position().width
    marker_width = max(hi-lo+1, 5*(ax.get_xlim()[1]-ax.get_xlim()[0])/width_pt)
    marker_height = 5 / height_pt
    ax.add_patch(Rectangle(((lo+hi)/2-marker_width/2, .62-marker_height/2),
                           marker_width, marker_height, color="#FFB000", lw=0, zorder=4))
    ax.text(gene["start"], .12, "RCA", fontsize=8, style="italic", va="center", color="#333333")
    ax.text((gene["start"]+gene["end"])/2+220, .12,
            "Chr2:16570746–16573692  (−)", ha="center", va="center", fontsize=6.5)
    return (lo+hi)/2, .62-marker_height/2


def draw_motif(ax, record):
    ax.set(xlim=WINDOW, ylim=(0, 1), xticks=[], yticks=[])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#777777")
        spine.set_linewidth(.6)
    lo, hi = record["start"]-.5, record["end"]+.5
    ax.add_patch(Rectangle((lo, .43), hi-lo, .35, facecolor="#F1EFF6", edgecolor="none"))
    ax.axhline(.43, color="#777777", lw=.6)
    ax.text(.02, .91, record["symbol"], transform=ax.transAxes, ha="left", va="center", weight="bold", fontsize=7)
    label = f"({record['strand']}) Chr2:{record['start']}–{record['end']}"
    ax.text(.98, .91, label, transform=ax.transAxes, ha="right", va="center", fontsize=5.5, color="#555555")
    for position, base in enumerate(record["display_sequence"], record["start"]):
        ax.text(position, .605, base, color=BASE_COLORS[base], fontsize=8, weight="bold", ha="center", va="center")
    left, right = ("5′", "3′") if record["strand"] == "+" else ("3′", "5′")
    ax.text(WINDOW[0]+.3, .605, left, fontsize=5.5, ha="left", va="center", color="#555555")
    ax.text(WINDOW[1]-.3, .605, right, fontsize=5.5, ha="right", va="center", color="#555555")
    ax.plot(record["position"], .265, marker="*", markersize=8, markeredgewidth=0, color=CHH_BLUE)
    ax.text(record["position"], .10, record["site_id"], ha="center", va="center", fontsize=6)


def connect(fig, gene_ax, anchor, axes, multi):
    """A shared connector reflects the overlapping loci without position jitter."""
    point = fig.transFigure.inverted().transform(gene_ax.transData.transform(anchor))
    if not multi:
        box = axes[0].get_position()
        center = (box.x0+box.x1)/2
        yy = box.y1+.045
        xs, ys = [point[0], point[0], center, center], [point[1], yy, yy, box.y1+.007]
        fig.add_artist(Line2D(xs, ys, transform=fig.transFigure, color="#444444", lw=.7))
        fig.add_artist(FancyArrowPatch((center, box.y1+.025), (center, box.y1+.002), transform=fig.transFigure,
                                      arrowstyle="->", mutation_scale=6, lw=.7, color="#444444"))
        return
    junction_y = .823
    fig.add_artist(Line2D([point[0], point[0], .5], [point[1], junction_y, junction_y],
                         transform=fig.transFigure, lw=.7, color="#555555"))
    bottom = min(a.get_position().y1+.018 for a in axes)
    fig.add_artist(Line2D([.5, .5], [junction_y, bottom], transform=fig.transFigure, lw=.7, color="#555555"))
    for ax in axes:
        box = ax.get_position()
        center, yy = (box.x0+box.x1)/2, box.y1+.018
        fig.add_artist(Line2D([.5, center, center], [yy, yy, box.y1+.006],
                             transform=fig.transFigure, lw=.7, color="#555555"))
        fig.add_artist(FancyArrowPatch((center, box.y1+.016), (center, box.y1+.002),
                                      transform=fig.transFigure, arrowstyle="->", mutation_scale=5, lw=.7, color="#555555"))


def merge_strand(records, strand):
    """Merge overlapping hits by genomic position, rejecting base conflicts."""
    selected = [r for r in records if r["strand"] == strand]
    bases = {}
    for record in selected:
        for position, base in enumerate(record["display_sequence"], record["start"]):
            if position in bases and bases[position] != base:
                raise ValueError(f"Inconsistent {strand}-strand base at {position}")
            bases[position] = base
    start, end = min(bases), max(bases)
    if len(bases) != end-start+1:
        raise ValueError("Merged motifs contain an unobserved sequence gap")
    return dict(selected[0], symbol=f"strand: {strand}", start=start, end=end,
                display_sequence="".join(bases[p] for p in range(start, end+1)))


def connect_strands(fig, gene_ax, anchor, axes, records):
    for index, (strand, ax) in enumerate(zip(("+", "-"), axes)):
        top = index == 0
        # Use the upper/lower marker surface, never a line through its interior.
        surface = (anchor[0], 1.24-anchor[1] if top else anchor[1])
        point = fig.transFigure.inverted().transform(gene_ax.transData.transform(surface))
        edge, elbow = (.63, .59) if top else (.40, .435)
        fig.add_artist(Line2D([point[0], point[0], .5], [point[1], elbow, elbow],
                             transform=fig.transFigure, color="#444444", lw=.8))
        fig.add_artist(FancyArrowPatch((.5, elbow), (.5, edge + (-.003 if top else .003)),
                                      transform=fig.transFigure, arrowstyle="->", mutation_scale=7,
                                      color="#444444", lw=.8, shrinkA=0, shrinkB=0))
        names = [r["symbol"] for r in records if r["strand"] == strand]
        width, gap = .12, .025
        left = .5-(len(names)*width+(len(names)-1)*gap)/2
        bottom = .885 if top else .065
        fig.text(left, .95 if top else .14, f"strand: {strand}", fontsize=7, weight="bold")
        for i, name in enumerate(names):
            x = left+i*(width+gap)
            fig.add_artist(Rectangle((x, bottom), width, .045, transform=fig.transFigure,
                                     facecolor="white", edgecolor="#333333", lw=.9))
            fig.text(x+width/2, bottom+.0225, name, ha="center", va="center", fontsize=7, weight="bold")
        start_y, end_y = (.833, .882) if top else (.197, .113)
        fig.add_artist(FancyArrowPatch((.5, start_y), (.5, end_y), transform=fig.transFigure,
                                      arrowstyle="->", mutation_scale=7, color="#444444", lw=.8))


def render(group, records, gene, out, dpi):
    sys.path.insert(0, str(QA_SCRIPTS))
    try:
        if not (QA_SCRIPTS / "audit_panel_alignment.py").is_file():
            raise ImportError(f"Optional QA tool not found: {QA_SCRIPTS / 'audit_panel_alignment.py'}")
        from audit_panel_alignment import require_matplotlib_panel_alignment
    except ImportError as exc:
        _skip_optional_qa(f"Optional alignment dependency unavailable: {exc}")
        require_matplotlib_panel_alignment = _skip_alignment
    multi = group == "non_coexpressed"
    width_inches, height_inches = 183/25.4, (145 if multi else 80)/25.4
    fig = plt.figure(figsize=(width_inches, height_inches))
    fig.legend([Line2D([], [], marker="*", color=CHH_BLUE, markersize=7, lw=0)], ["CHH"],
               loc="center right", bbox_to_anchor=(.96, .96), fontsize=7, frameon=False)
    arrow_y = .555 if multi else .82
    arrow_start, arrow_end = (.80, .65) if multi else (.93, .78)
    fig.add_artist(FancyArrowPatch((arrow_start, arrow_y), (arrow_end, arrow_y),
                                  transform=fig.transFigure, arrowstyle="-|>",
                                  mutation_scale=9, lw=.9, color="#333333"))
    fig.text((arrow_start+arrow_end)/2, arrow_y+3/(height_inches*72),
             "Transcription direction", ha="center", va="bottom", fontsize=6.5,
             color="#333333")
    gene_ax = fig.add_axes([.045, .46, .91, .10] if multi else [.045, .65, .91, .16])
    anchor = draw_gene(gene_ax, gene, records)
    axes = []
    plotted = [merge_strand(records, strand) for strand in ("+", "-")] if multi else records
    for index, record in enumerate(plotted):
        if multi:
            rect = [.12, .63 if index == 0 else .20, .76, .20]
        else:
            rect = [.12, .12, .76, .37]
        ax = fig.add_axes(rect, label=record["symbol"])
        draw_motif(ax, record)
        axes.append(ax)
    fig.canvas.draw()
    if multi:
        connect_strands(fig, gene_ax, anchor, axes, records)
    else:
        connect(fig, gene_ax, anchor, axes, False)
    prefix = out / f"CHH_TF_{group}"
    names = [r["symbol"] for r in plotted]
    row_groups = []
    column_groups = [names] if multi else []
    require_matplotlib_panel_alignment(fig, json_out=f"{prefix}.alignment.json", exclude_axes=[gene_ax],
                                       row_groups=row_groups, column_groups=column_groups,
                                       tolerance_pt=1.5, gutter_tolerance_pt=1.5, strict=True)
    fig.savefig(f"{prefix}.pdf", facecolor="white")
    fig.savefig(f"{prefix}.svg", facecolor="white")
    fig.savefig(f"{prefix}.png", dpi=dpi, facecolor="white")
    fig.savefig(f"{prefix}.tiff", dpi=dpi, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    return prefix


def write_tsv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def run_qa(prefix):
    pdf = f"{prefix}.pdf"
    commands = [
        ("pdf-text", [sys.executable, str(QA_SCRIPTS / "audit_pdf_text.py"), pdf, "--min-pt", "5", "--json"]),
        ("collision", [sys.executable, str(QA_SCRIPTS / "audit_figure_collisions.py"), pdf, "--strict", "--json-out", f"{prefix}.collision-audit.json", "--overlay-pdf", f"{prefix}.collision-audit.pdf"]),
    ]
    for name, command in commands:
        if _skip_missing_qa_command(command, Path(f"{prefix}.{name}.txt")):
            continue
        result = subprocess.run(command, capture_output=True, text=True)
        Path(f"{prefix}.{name}.txt").write_text(result.stdout + result.stderr)
        if result.returncode:
            raise RuntimeError(f"{name} QA failed for {prefix.name}: {result.stdout[-3000:]}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, default=DEFAULT_WORKDIR)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--gff", type=Path, default=DEFAULT_GFF)
    parser.add_argument("--dpi", type=int, default=DEFAULT_DPI)
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    groups, counts = load_data(args.workdir)
    gene = load_gene(args.gff)
    if args.validate_only:
        print(json.dumps(dict(counts=counts, gene=gene), indent=2))
        return
    out = args.output_dir or args.workdir / "meth_experiment/CHH_TF_coexloci"
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"], "font.size": 7, "pdf.fonttype": 42,
                         "svg.fonttype": "none", "axes.linewidth": .6})
    rows = [r for records in groups.values() for r in records]
    write_tsv(out / "CHH_TF_plot_data.tsv", rows)
    merged_rows = []
    for strand in ("+", "-"):
        merged = merge_strand(groups["non_coexpressed"], strand)
        merged_rows.append(dict(strand=strand, start=merged["start"], end=merged["end"],
                                display_sequence=merged["display_sequence"],
                                symbols=",".join(r["symbol"] for r in groups["non_coexpressed"] if r["strand"] == strand),
                                site_id="CHH_1651", position=16574877))
    write_tsv(out / "CHH_TF_merged_strands.tsv", merged_rows)
    write_tsv(out / "RCA_gene_model.tsv", [dict(transcript="NM_129531.3", feature=kind, start=s, end=e, strand="-")
              for kind in ("exons", "cds") for s, e in gene[kind]])
    (out / "selection_counts.json").write_text(json.dumps(counts, indent=2)+"\n")
    for group, records in groups.items():
        prefix = render(group, records, gene, out, args.dpi)
        run_qa(prefix)
        print(f"Exported: {prefix}; QA {'INCOMPLETE' if _OPTIONAL_QA_SKIPPED else 'PASS'}")
    (out / "README.md").write_text(
        "# RCA motif and CHH spatial maps\n\n"
        "These figures show predicted motif overlap with CHH_1651, not experimental TF binding or causal regulation. "
        "Coexpression classifications come from the supplied tables. Only motif=1 rows are included; "
        "other CHH sites and flanking-only matches are excluded as requested. No methylation magnitude or statistics are encoded.\n\n"
        "All displayed positions are 1-based inclusive. Numbered CHH input is single-base BED: position = start + 1 = end. "
        "FIMO intervals and GFF coordinates are 1-based inclusive. NC_003071.7 is displayed as Chr2. "
        "CHH_1651 is at 16574877. Every local panel uses Chr2:16574872–16574890. "
        "Half-base drawing boundaries enclose nucleotide centers; they do not alter genomic coordinates.\n\n"
        "FIMO matched_sequence is recorded in motif-strand 5′–3′ orientation. Plus-strand strings are shown unchanged; "
        "minus-strand strings are reversed, not complemented, so increasing x retains increasing genomic position "
        "and the displayed minus strand runs 3′–5′. A G opposite a methylated C on the complementary strand is expected.\n\n"
        "RCA is AT2G39730.1 / NM_129531.3, Chr2:16570746–16573692, negative strand. "
        "Seven exons are taken from the reference GFF; CDS is solid blue and untranslated exon segments are white with blue borders. "
        "Dashed lines and group headings are omitted. Amber marks the motif interval union center, "
        "with a compact 5 pt square overview marker for legibility (not to genomic scale); "
        "local motif panels retain exact coordinates. Exon boxes have a consistent 8 pt height; "
        "blue connecting lines occur only in introns and flanking DNA, never inside hollow UTRs, and have no arrowheads. "
        "A separate dark-gray left-pointing arrow indicates RCA negative-strand transcription: below the legend "
        "in the coexpressed figure and above the gene overview in the non-coexpressed figure. "
        "Black connectors meet the top or bottom marker surface without crossing its interior; elbow segments join without shrink gaps. "
        "the shared connector links overlapping hits without shifting genomic coordinates.\n\n"
        "PDF/SVG retain editable text; PNG/TIFF are 600 dpi by default. "
        "Alignment tolerance is 1.5 pt. Available external PDF text and collision audits run on export; missing tools are recorded as SKIPPED. "
        "Coexpressed motif alignment is not applicable because it has one local panel. "
        "The coexpressed figure is 183 x 80 mm. The non-coexpressed figure is 183 x 145 mm with a merged plus-strand panel above RCA and "
        "a merged minus-strand panel below RCA. Overlapping bases must agree before merging; "
        "no majority consensus or inferred sequence is used. TF names appear in separate boxes "
        "above and below the corresponding merged panels. Individual hit boundaries remain in the source TSV.\n",
        encoding="utf-8")


if __name__ == "__main__":
    main()
