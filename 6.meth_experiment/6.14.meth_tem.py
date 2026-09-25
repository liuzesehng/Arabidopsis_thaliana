#!/usr/bin/env python3
"""RCA negative-strand methylation profiles: per-site temperature means.

No binning, smoothing, imputation, coverage weighting or significance tests.
Coordinates are BED (0-based, half-open); sample rows are replicate units.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rca-profile-mpl"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parents[2]
DEFAULT_DIRECTORY = PROJECT / "list/meth_experiment"
DEFAULT_INPUT = DEFAULT_DIRECTORY / "Alt.all.snp_meth.filtered.tsv"
DEFAULT_QA_TOOLS = Path("/datapool/home/2023102768/.codex/skills/nature-figure/scripts")
DEFAULT_DPI = 600
FIG_WIDTH_MM = 183
FIG_HEIGHT_MM = 200
PREFIX = "RCA_temperature_methylation_profile"
CONTEXTS = ("CG", "CHG", "CHH")
TEMPERATURES = tuple(range(4, 33, 4))
CHROMOSOME = "NC_003071.7"
GENE_START, GENE_END, FLANK = 16570745, 16573692, 2000
COLORS = ("#3B4CC0", "#2878B5", "#249EA0", "#629B55",
          "#B49125", "#E17C30", "#CC4245", "#823B85")


def read_rows(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return [row for row in csv.reader(handle, delimiter="\t")
                if row and not row[0].startswith("#")]


def locate(start, end):
    """Return region, interval midpoint and monotonic transcript display x.

    Display intervals are [-1, 0], [0, 3], [3, 4], hence 1:3:1 widths.
    """
    if not (GENE_START - FLANK <= start < end <= GENE_END + FLANK):
        raise ValueError(f"interval outside RCA +/-2 kb: {start}, {end}")
    for boundary in (GENE_START, GENE_END):
        if start < boundary < end:
            raise ValueError(f"site interval straddles a gene boundary: {start}, {end}")
    midpoint = (start + end) / 2
    if start >= GENE_END:
        return "upstream", midpoint, -(midpoint - GENE_END) / FLANK
    if end <= GENE_START:
        return "downstream", midpoint, 3 + (GENE_START - midpoint) / FLANK
    return "gene_body", midpoint, 3 * (GENE_END - midpoint) / (GENE_END - GENE_START)


def load_data(input_path, position_paths):
    rows = read_rows(input_path)
    if not rows or len(rows) < 2:
        raise ValueError("methylation table is empty")
    header = rows[0]
    if any(not h.strip() for h in header) or len(set(header)) != len(header):
        raise ValueError("blank or duplicate table column names")
    if any(len(r) != len(header) for r in rows[1:]):
        raise ValueError("inconsistent TSV field counts")
    if "Temperature" not in header:
        raise ValueError("missing Temperature column")
    data = pd.DataFrame(rows[1:], columns=header)
    sample = data.Temperature.str.extract(r"^T(\d+)-(\d+)$")
    if sample.isna().any().any() or data.Temperature.duplicated().any():
        raise ValueError("sample names must be unique T<temperature>-<replicate> IDs")
    data["temperature"] = sample[0].astype(int)
    data["replicate"] = sample[1].astype(int)
    if (data.replicate < 1).any() or data.duplicated(["temperature", "replicate"]).any():
        raise ValueError("replicate numbers must be positive and unique within temperature")
    if set(data.temperature) != set(TEMPERATURES):
        raise ValueError("expected temperature groups (4, 8, 12, 16, 20, 24, 28, 32)")
    features = [c for c in header if re.fullmatch(r"(?:CG|CHG|CHH)_\d+", c)]
    invalid = [c for c in header if c.startswith(("CG_", "CHG_", "CHH_")) and c not in features]
    if invalid:
        raise ValueError(f"invalid site identifiers: {invalid}")
    for feature in features:
        raw = data[feature].str.strip()
        missing = raw.isin(["", "NA"])
        values = pd.to_numeric(raw.mask(missing), errors="coerce")
        if ((~missing) & values.isna()).any() or np.isinf(values).any():
            raise ValueError(f"invalid methylation value in {feature}")
        if ((values < 0) | (values > 100)).any():
            raise ValueError(f"{feature}: methylation must be 0-100 percent")
        data[feature] = values
    sites, counts = [], {}
    for context in CONTEXTS:
        wanted = {f for f in features if f.startswith(context + "_")}
        if not wanted:
            raise ValueError(f"no {context} columns")
        records, seen, coordinates = {}, set(), set()
        for row in read_rows(position_paths[context]):
            if len(row) != 6:
                raise ValueError(f"{context}: expected six columns in numbered position table")
            identifier, chrom, begin, finish = row[:4]
            if not re.fullmatch(context + r"_\d+", identifier) or identifier in seen:
                raise ValueError(f"{context}: invalid or duplicate ID {identifier}")
            seen.add(identifier)
            start, end = int(begin), int(finish)
            if start < 0 or end <= start:
                raise ValueError(f"{identifier}: invalid BED interval")
            key = chrom, start, end
            if key in coordinates:
                raise ValueError(f"{context}: duplicate coordinate {key}")
            coordinates.add(key)
            if identifier not in wanted:
                continue
            if chrom != CHROMOSOME:
                raise ValueError(f"{identifier}: wrong chromosome {chrom}")
            region, midpoint, x = locate(start, end)
            records[identifier] = dict(feature=identifier, context=context, chromosome=chrom,
                                       start=start, end=end, midpoint=midpoint, strand="-",
                                       region=region, x_plot=x)
        if wanted - records.keys():
            raise ValueError(f"{context}: unmapped sites {sorted(wanted - records.keys())}")
        sites.extend(records.values())
        counts[context] = dict(position_rows=len(seen), selected_sites=len(wanted),
                               unselected_position_rows=len(seen) - len(wanted))
    sites = pd.DataFrame(sites)
    sites["context"] = pd.Categorical(sites.context, CONTEXTS, ordered=True)
    sites = sites.sort_values(["context", "x_plot"]).reset_index(drop=True)
    raw = data.rename(columns={"Temperature": "sample"}).melt(
        id_vars=["sample", "temperature", "replicate"], value_vars=features,
        var_name="feature", value_name="methylation_pct").merge(sites, on="feature", validate="many_to_one")
    raw["missing"] = raw.methylation_pct.isna()
    raw = raw.sort_values(["context", "x_plot", "temperature", "replicate"]).reset_index(drop=True)
    return data, sites, raw, counts, [c for c in header if c not in features and c != "Temperature"]


def summarize(raw, sites):
    summary = raw.groupby(["feature", "temperature"], sort=False).methylation_pct.agg(
        mean_pct="mean", sd_pct="std", n_valid="count", n_total="size").reset_index()
    summary["n_missing"] = summary.n_total - summary.n_valid
    return summary.merge(sites, on="feature", validate="many_to_one").sort_values(
        ["context", "temperature", "x_plot"]).reset_index(drop=True)


def load_auditor(directory):
    spec = importlib.util.spec_from_file_location("rca_profile_alignment", directory / "audit_panel_alignment.py")
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load alignment auditor")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def make_figure(summary, auditor, qa_dir):
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"], "font.size": 7,
                         "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.6, "pdf.fonttype": 42, "svg.fonttype": "none",
                         "legend.frameon": False, "savefig.facecolor": "white"})
    fig, axes = plt.subplots(3, 1, figsize=(FIG_WIDTH_MM / 25.4, FIG_HEIGHT_MM / 25.4))
    # Reserve 20 mm below the original layout without compressing the panels.
    fig.subplots_adjust(left=0.105, right=0.975, top=177.5 / FIG_HEIGHT_MM,
                        bottom=38.9 / FIG_HEIGHT_MM, hspace=0.66)
    records = []
    panel_labels = []
    for index, (ax, context) in enumerate(zip(axes, CONTEXTS)):
        ax.set_gid(chr(97 + index))
        block = summary[summary.context == context]
        for temperature, color in zip(TEMPERATURES, COLORS):
            series = block[block.temperature == temperature].sort_values("x_plot")
            x, y = series.x_plot.to_numpy(), series.mean_pct.to_numpy()
            if not np.all(np.diff(x) > 0):
                raise ValueError(f"{context}: non-increasing site coordinates")
            line, = ax.plot(x, y, color=color, linewidth=0.65, marker="o", markersize=1.5,
                            markeredgewidth=0, label=f"{temperature} °C", clip_on=False)
            # Numerical plot/table equivalence is checked on the actual artists.
            np.testing.assert_allclose(line.get_xdata(), x)
            np.testing.assert_allclose(line.get_ydata(), y, equal_nan=True)
        for boundary in (0, 3):
            ax.axvline(boundary, color="#A8A8A8", linewidth=0.55, linestyle=(0, (3, 3)), zorder=0)
        ax.set(xlim=(-1, 4), ylim=(0, 100), ylabel="DNA methylation level (%)")
        ax.set_yticks([0, 20, 40, 60, 80, 100])
        ax.set_xticks([-1, 0, 0.6, 1.2, 1.8, 2.4, 3, 4],
                      ["−2 kb", "TSS\n0%", "20%", "40%", "60%", "80%", "TES\n100%", "+2 kb"])
        ax.tick_params(length=2.5, width=0.6, pad=2)
        panel_labels.append(ax.text(0, 1.08, f"({chr(97 + index)})", transform=ax.transAxes,
                                    fontsize=9, fontweight="bold"))
        ax.text(0.055, 1.08, context, transform=ax.transAxes, fontsize=8, fontweight="bold")
        for left, width, color, label in [(-1, 1, "#B4B4B4", "Upstream"),
                                          (0, 3, "#5367A0", "Gene body"),
                                          (3, 1, "#B4B4B4", "Downstream")]:
            ax.add_patch(Rectangle((left, -0.29), width, 0.035, transform=ax.get_xaxis_transform(),
                                   facecolor=color, edgecolor="none", clip_on=False))
            ax.text(left + width / 2, -0.36, label, transform=ax.get_xaxis_transform(),
                    ha="center", va="center", fontsize=6.5)
        records.append(dict(panel=chr(97 + index), context=context, curves=8,
                            sites=int(block.feature.nunique()), finite_means=int(block.mean_pct.notna().sum()),
                            summary="per-site arithmetic mean", uncertainty="none; sample SD in source table",
                            replicate_unit="input sample row", y_limits=[0, 100]))
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=4, loc="upper right", bbox_to_anchor=(0.975, 197.3 / FIG_HEIGHT_MM),
               fontsize=7, handlelength=2, handletextpad=0.5, columnspacing=1.5)
    sample_counts = summary.groupby("temperature").n_total.first()
    count_text = "; ".join(f"{int(t)} °C: {int(n)}" for t, n in sample_counts.items())
    figure_caption = [
        "Points show per-site temperature means for filtered CG, CHG and CHH sites in RCA and its 2 kb flanks.",
        "Lines connect retained sites only; no binning or smoothing. They do not represent unmeasured positions.",
        "Sample rows per temperature (n): " + count_text + ".",
        "Zeros are retained; missing values are omitted without imputation. Sample SD and valid n are in the source table.",
        "Negative-strand RCA is shown in transcription direction; the gene body is scaled to 0–100% (flank:body:flank = 1:3:1).",
    ]
    fig.text(0.105, 22 / FIG_HEIGHT_MM, "\n".join(figure_caption), ha="left", va="top",
             fontsize=6.5, linespacing=1.4)
    fig.canvas.draw()
    auditor.require_matplotlib_panel_alignment(fig,
        tolerance_pt=1.5, gutter_tolerance_pt=1.5, strict=True)
    # Supply actual anchors because the auditor does not infer parenthesized labels.
    manifest = auditor.matplotlib_layout_manifest(fig)
    for panel, label in zip(manifest["panels"], panel_labels):
        anchor = label.get_transform().transform(label.get_position())
        anchor = fig.dpi_scale_trans.inverted().transform(anchor) * 72
        panel["panel_label"] = label.get_text().strip("()")
        panel["panel_label_anchor_pt"] = list(map(float, anchor))
    alignment = auditor.audit_layout_manifest(manifest, tolerance_pt=1.5,
        gutter_tolerance_pt=1.5, require_panel_labels=True)
    auditor.write_json_report(alignment, qa_dir / "alignment.json")
    if auditor.exit_code(alignment, strict=True):
        raise RuntimeError(auditor.render_text(alignment, strict=True))
    return fig, records


def run_audit(command, output):
    result = subprocess.run(command, text=True, capture_output=True)
    output.write_text(result.stdout + result.stderr, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(f"QA failed ({result.returncode}): {output}\n{result.stdout}{result.stderr}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--positions-dir", type=Path, default=DEFAULT_DIRECTORY)
    for context in CONTEXTS:
        parser.add_argument(f"--{context.lower()}-positions", type=Path)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_DIRECTORY)
    parser.add_argument("--qa-tools", type=Path, default=DEFAULT_QA_TOOLS)
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args(argv)
    paths = {c: getattr(args, f"{c.lower()}_positions") or args.positions_dir / f"{c}.CV.before.numbered.tsv"
             for c in CONTEXTS}
    data, sites, raw, counts, ignored = load_data(args.input, paths)
    summary = summarize(raw, sites)
    sample_counts = {int(t): int(n) for t, n in data.groupby("temperature").size().items()}
    if args.validate_only:
        print(json.dumps(dict(samples=sample_counts, sites=counts, raw_rows=len(raw), means=len(summary)), indent=2))
        return 0
    auditor = load_auditor(args.qa_tools)
    if importlib.util.find_spec("fitz") is None:
        raise RuntimeError("PDF collision QA requires PyMuPDF (fitz)")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    qa_dir = args.output_dir / f"{PREFIX}_QA"
    qa_dir.mkdir(exist_ok=True)
    stem = args.output_dir / PREFIX
    run_audit([sys.executable, str(args.qa_tools / "validate_figure.py"), __file__], qa_dir / "source-preflight.txt")
    fig, records = make_figure(summary, auditor, qa_dir)
    try:
        fig.savefig(stem.with_suffix(".pdf"))
        run_audit([sys.executable, str(args.qa_tools / "audit_pdf_text.py"), str(stem.with_suffix(".pdf")),
                   "--min-pt", "5", "--json"], qa_dir / "pdf-text.json")
        run_audit([sys.executable, str(args.qa_tools / "audit_figure_collisions.py"), str(stem.with_suffix(".pdf")),
                   "--strict", "--json-out", str(qa_dir / "collision-audit.json")], qa_dir / "collision-audit.txt")
        fig.savefig(stem.with_suffix(".svg"))
        fig.savefig(stem.with_suffix(".png"), dpi=DEFAULT_DPI)
        fig.savefig(stem.with_suffix(".tiff"), dpi=DEFAULT_DPI, pil_kwargs={"compression": "tiff_lzw"})
    finally:
        plt.close(fig)
    for name, table in [("sites", sites), ("replicates", raw), ("summary", summary)]:
        table.to_csv(stem.with_suffix(f".{name}.tsv"), sep="\t", index=False, na_rep="NA", float_format="%.12g")
    caption = """# RCA 温度分组甲基化分布 / Temperature-grouped RCA methylation profiles

Panels a–c show CG, CHG and CHH, respectively. Points are per-site arithmetic
means across the available sample rows at each temperature; lines connect
adjacent retained sites and do not represent measurements at unobserved positions.
No binning, smoothing, extrapolation, coverage weighting or significance tests
are applied. No uncertainty band is shown; sample SD (ddof=1), valid n and missing
n are provided in the summary table. Zero is a measured value, not missing;
blank/NA observations are omitted without imputation. A wholly missing site
breaks the curve; SD is NA when fewer than two replicates are available.

RCA (AT2G39730), NC_003071.7:[16570745,16573692), is on the negative strand.
BED coordinates are 0-based and half-open. Interval midpoints locate each site.
The left-to-right axis follows transcription, from high to low genomic coordinates:
upstream [16573692,16575692), gene body [16570745,16573692), downstream
[16568745,16570745). The body is scaled to 0–100%; the flanks each span 2 kb.
Display widths are 1:3:1 and thus are not a uniform genomic-distance scale.

本图仅描述过滤后保留位点，不代表 RCA 区域全部胞嘧啶的无偏平均水平。
曲线仅连接保留位点的温度均值，不代表未测位置的甲基化水平。
样本行作为重复单位；不从样本名称推断其为生物学重复或技术重复。
"""
    caption += "\nSample rows per temperature: " + ", ".join(f"{t} °C: n={n}" for t, n in sample_counts.items()) + ".\n"
    stem.with_suffix(".caption.md").write_text(caption, encoding="utf-8")
    report = dict(status="automated_QA_PASS", visual_review="pending", panels=records,
                  samples=sample_counts, sample_ids=data.Temperature.tolist(), site_counts=counts,
                  region_counts=sites.groupby(["context", "region"], observed=True).size().rename("n").reset_index().to_dict("records"),
                  source_rows=len(data), raw_rows=len(raw), mean_rows=len(summary),
                  missing_measurements=int(raw.missing.sum()), ignored_non_methylation_columns=ignored,
                  exclusions="No methylation columns or input sample rows excluded; unused coordinate rows lack filtered-table columns.",
                  coordinate_system="BED 0-based half-open", chromosome=CHROMOSOME,
                  gene_start=GENE_START, gene_end=GENE_END, strand="-", flank_bp=FLANK,
                  aggregation="equal-weight sample mean at each retained site; sample SD ddof=1",
                  inputs={name: dict(path=str(path.resolve()), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
                          for name, path in {"methylation": args.input, **paths}.items()},
                  figure_mm=[FIG_WIDTH_MM, FIG_HEIGHT_MM], raster_dpi=DEFAULT_DPI,
                  versions=dict(python=sys.version.split()[0], numpy=np.__version__, pandas=pd.__version__, matplotlib=matplotlib.__version__))
    (qa_dir / "summary.json").write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(f"PASS: {len(data)} samples, {len(sites)} sites, {len(summary)} means; PDF/SVG/PNG/TIFF: {stem}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, RuntimeError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
