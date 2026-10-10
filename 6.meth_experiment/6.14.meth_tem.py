#!/usr/bin/env python3
"""RCA negative-strand methylation profiles: cgmaptools bins of all input sites.

Sample-level bin means followed by temperature means; PCHIP display interpolation.
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
from fractions import Fraction

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rca-profile-mpl"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.transforms import ScaledTranslation
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

PROJECT = Path(__file__).resolve().parents[2]
DEFAULT_DIRECTORY = PROJECT / "list/meth_experiment"
DEFAULT_INPUT = DEFAULT_DIRECTORY / "Alt.all.snp_meth.tsv"
DEFAULT_QA_TOOLS = Path("/datapool/home/2023102768/.codex/skills/nature-figure/scripts")
DEFAULT_DPI = 600
DEFAULT_BIN_SIZE = 400
DEFAULT_CGMAPTOOLS = Path("/datapool/life-gongl/zesheng/program/cgmaptools-0.1.3/cgmaptools")
FIG_WIDTH_MM = 183
FIG_HEIGHT_MM = 240
BODY_DISPLAY_WIDTH = 2.0
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

    Display intervals are [-1, 0], [0, 2], [2, 3], hence 1:2:1 widths.
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
        return "downstream", midpoint, BODY_DISPLAY_WIDTH + (GENE_START - midpoint) / FLANK
    return "gene_body", midpoint, BODY_DISPLAY_WIDTH * (GENE_END - midpoint) / (GENE_END - GENE_START)


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


def bin_with_cgmaptools(raw, sites, executable, bin_size, work_dir):
    """Use mbin on lossless ratio encodings, never invent biological coverage.

    Region-local 1-based coordinates run in transcription direction. MC/NC
    encode the input percentage exactly as a rational number; they are NOT
    read counts. All denominators pass -c 1; mbin averages sites equally.
    """
    if bin_size <= 0:
        raise ValueError("bin size must be positive")
    if not executable.is_file():
        raise ValueError(f"cgmaptools not found: {executable}")
    if not (sites.end - sites.start == 1).all():
        raise ValueError("cgmaptools adapter requires single-base site intervals")
    regions = {"upstream": (GENE_END, GENE_END + FLANK),
               "gene_body": (GENE_START, GENE_END),
               "downstream": (GENE_START - FLANK, GENE_START)}
    bins = []
    for region, (lower, upper) in regions.items():
        for offset in range(0, upper - lower, bin_size):
            end = upper - offset
            start = max(lower, end - bin_size)
            _, midpoint, x = locate(start, end)
            bins.append(dict(bin_id=f"{region}_{offset // bin_size + 1:03d}", region=region,
                             chromosome=CHROMOSOME, start=start, end=end, width_bp=end-start,
                             midpoint=midpoint, x_plot=x, local_start=offset+1,
                             local_end=min(offset+bin_size, upper-lower)))
    bins = pd.DataFrame(bins)
    assigned = sites.copy()
    assigned["local_position"] = [regions[r.region][1] - int(r.end) + 1 for r in assigned.itertuples()]
    assigned["bin_id"] = [f"{r.region}_{(r.local_position-1)//bin_size+1:03d}" for r in assigned.itertuples()]
    sample_sites = raw.merge(assigned[["feature", "bin_id", "local_position"]], on="feature", validate="many_to_one")
    work_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / "README.txt").write_text(
        "These CGmap files are computational adapters of input-table percentages.\n"
        "MC/NC are exact rational encodings, NOT methylated reads or coverage.\n"
        "Chromosome labels are region names and positions are region-local, 1-based,\n"
        "in transcription direction. Use exported sites/bins tables for genomic coordinates.\n"
        "No coverage selection or weighting is performed; mbin -c 1 includes every encoded value.\n"
        "cgmaptools mbin rounds fractions to four decimals (0.01 percentage-point resolution).\n")
    records, commands = [], []
    for (sample, context), block in sample_sites.groupby(["sample", "context"], observed=True, sort=True):
        adapter = work_dir / f"{sample}.{context}.ratio_encoded.CGmap"
        lines = []
        for region in regions:
            for row in block[block.region == region].sort_values("local_position").itertuples():
                if pd.isna(row.methylation_pct):
                    continue
                ratio = Fraction(str(row.methylation_pct)) / 100
                lines.append(f"{region}\tC\t{row.local_position}\t{context}\tNA\t{float(ratio):.12g}\t{ratio.numerator}\t{ratio.denominator}\n")
        adapter.write_text("".join(lines))
        command = [sys.executable, str(executable), "mbin", "-i", str(adapter),
                   "-B", str(bin_size), "-c", "1", "-C", str(context)]
        result = subprocess.run(command, text=True, capture_output=True)
        (work_dir / f"{sample}.{context}.mbin.tsv").write_text(result.stdout)
        (work_dir / f"{sample}.{context}.stderr.txt").write_text(result.stderr)
        commands.append(dict(command=command, returncode=result.returncode))
        if result.returncode or result.stderr.strip():
            raise RuntimeError(f"cgmaptools failed for {sample}/{context}: {result.stderr}")
        values = {}
        for line in result.stdout.splitlines():
            fields = line.split("\t")
            if len(fields) != 4:
                raise RuntimeError(f"unexpected cgmaptools output: {line}")
            region, left, right, value = fields
            key = region, int(left)
            if key in values or region not in regions or int(right)-int(left)+1 != bin_size:
                raise RuntimeError(f"invalid or duplicate cgmaptools bin: {line}")
            values[key] = np.nan if value.lower() == "na" else float(value) * 100
        for b in bins.itertuples(index=False):
            subset = block[block.bin_id == b.bin_id].methylation_pct
            valid = int(subset.count())
            value = values.get((b.region, b.local_start), np.nan)
            expected = float(subset.mean()) if valid else np.nan
            if valid:
                if not np.isfinite(value) or abs(value-expected) > 0.005001:
                    raise RuntimeError(f"cgmaptools/Python mean mismatch: {sample}/{context}/{b.bin_id}")
            elif not np.isnan(value):
                raise RuntimeError("empty bin must remain missing")
            records.append(dict(sample=sample, temperature=int(block.temperature.iloc[0]),
                                replicate=int(block.replicate.iloc[0]), context=str(context),
                                bin_id=b.bin_id, methylation_pct=value,
                                n_sites_total=len(subset), n_sites_valid=valid,
                                n_sites_missing=len(subset)-valid))
    (work_dir / "commands.json").write_text(json.dumps(commands, indent=2)+"\n")
    per_sample = pd.DataFrame(records)
    summary = per_sample.groupby(["context", "bin_id", "temperature"]).methylation_pct.agg(
        mean_pct="mean", sd_pct="std", n_valid="count", n_total="size").reset_index()
    summary["n_missing"] = summary.n_total-summary.n_valid
    summary = summary.merge(bins, on="bin_id", validate="many_to_one").sort_values(["context", "temperature", "x_plot"])
    return assigned, bins, per_sample.merge(bins, on="bin_id", validate="many_to_one"), summary


def smooth_segments(x, y):
    """Interpolate consecutive finite bins, preserving anchors and missing gaps."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if not np.all(np.diff(x) > 0):
        raise ValueError("interpolation coordinates must be strictly increasing")
    valid = np.flatnonzero(np.isfinite(y))
    runs = np.split(valid, np.flatnonzero(np.diff(valid) > 1) + 1)
    xs, ys = [], []
    for run in runs:
        if not len(run):
            continue
        sx, sy = x[run], y[run]
        if len(run) >= 2:
            interpolator = PchipInterpolator(sx, sy, extrapolate=False)
            np.testing.assert_allclose(interpolator(sx), sy, atol=1e-12)
            dense_x = np.unique(np.concatenate([sx, np.linspace(sx[0], sx[-1], max(160, 80 * len(run)))]))
            dense_y = interpolator(dense_x)
            if dense_y.min() < sy.min()-1e-10 or dense_y.max() > sy.max()+1e-10:
                raise ValueError("interpolated curve overshoots observed values")
        else:
            dense_x, dense_y = sx, sy
        xs.extend([*dense_x, np.nan])
        ys.extend([*dense_y, np.nan])
    return np.asarray(xs), np.asarray(ys)


def make_figure(summary, auditor, qa_dir, bin_size):
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"], "font.size": 7,
                         "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.6, "pdf.fonttype": 42, "svg.fonttype": "none",
                         "legend.frameon": False, "savefig.facecolor": "white"})
    fig = plt.figure(figsize=(FIG_WIDTH_MM / 25.4, FIG_HEIGHT_MM / 25.4))
    grid = fig.add_gridspec(5, 1, height_ratios=[33, 22, 33, 22, 70],
                           left=0.105, right=0.975, top=217.5 / FIG_HEIGHT_MM,
                           bottom=37.5 / FIG_HEIGHT_MM, hspace=0)
    axes = [fig.add_subplot(grid[i, 0]) for i in (0, 2, 4)]
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
            curve_x, curve_y = smooth_segments(x, y)
            line, = ax.plot(curve_x, curve_y, color=color, linewidth=0.95, marker=None,
                            solid_capstyle="round", solid_joinstyle="round",
                            label=f"{temperature} °C", clip_on=False)
            np.testing.assert_allclose(line.get_xdata(), curve_x, equal_nan=True)
            np.testing.assert_allclose(line.get_ydata(), curve_y, equal_nan=True)
        for boundary in (0, BODY_DISPLAY_WIDTH):
            ax.axvline(boundary, color="#A8A8A8", linewidth=0.55, linestyle=(0, (3, 3)), zorder=0)
        finite_means = block.mean_pct.dropna()
        if finite_means.empty:
            raise ValueError(f"{context}: no finite means for y-axis scaling")
        peak_value = float(finite_means.max())
        y_max = {"CG": 100, "CHG": 60, "CHH": 30}[context]
        tick_step = {"CG": 20, "CHG": 20, "CHH": 10}[context]
        if peak_value > y_max:
            raise ValueError(f"{context}: peak {peak_value} exceeds requested y-axis upper limit {y_max}")
        ticks = list(range(0, y_max + 1, tick_step))
        tick_labels = [f"{t:.2f}".rstrip("0").rstrip(".") for t in ticks]
        ax.set(xlim=(-1, BODY_DISPLAY_WIDTH + 1), ylim=(0, y_max), ylabel="DNA methylation level (%)")
        ax.set_yticks(ticks, tick_labels)
        ax.set_xticks([-1, 0, *[BODY_DISPLAY_WIDTH * p for p in (0.2, 0.4, 0.6, 0.8)],
                       BODY_DISPLAY_WIDTH, BODY_DISPLAY_WIDTH + 1],
                      ["−2 kb", "TSS\n0%", "20%", "40%", "60%", "80%", "TES\n100%", "+2 kb"])
        ax.tick_params(length=2.5, width=0.6, pad=2)
        label_transform = ax.transAxes + ScaledTranslation(0, 3 / 25.4, fig.dpi_scale_trans)
        panel_labels.append(ax.text(0, 1, f"({chr(97 + index)})", transform=label_transform,
                                    fontsize=9, fontweight="bold"))
        ax.text(0.055, 1, context, transform=label_transform, fontsize=8, fontweight="bold")
        axes_height_mm = ax.get_position().height * FIG_HEIGHT_MM
        for left, width, color, label in [(-1, 1, "#B4B4B4", "Upstream"),
                                          (0, BODY_DISPLAY_WIDTH, "#5367A0", "Gene body"),
                                          (BODY_DISPLAY_WIDTH, 1, "#B4B4B4", "Downstream")]:
            ax.add_patch(Rectangle((left, -9.5 / axes_height_mm), width, 1 / axes_height_mm, transform=ax.get_xaxis_transform(),
                                   facecolor=color, edgecolor="none", clip_on=False))
            ax.text(left + width / 2, -12.5 / axes_height_mm, label, transform=ax.get_xaxis_transform(),
                    ha="center", va="center", fontsize=6.5)
        if context == "CHH":
            peak = block[block.region == "upstream"].groupby("bin_id").mean_pct.mean().idxmax()
            peak_row = block[block.bin_id == peak].iloc[0]
            peak_start, peak_end = int(peak_row.start), int(peak_row.end)
        records.append(dict(panel=chr(97 + index), context=context, curves=8,
                            bins=int(block.bin_id.nunique()), finite_means=int(block.mean_pct.notna().sum()),
                            summary="temperature mean of sample-level cgmaptools bin means", uncertainty="none; sample SD in source table",
                            replicate_unit="input sample row", y_limits=[0, y_max], y_ticks=ticks,
                            peak_mean_pct=peak_value, height_mm=axes_height_mm))
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=4, loc="upper right", bbox_to_anchor=(0.975, 237.3 / FIG_HEIGHT_MM),
               fontsize=7, handlelength=2, handletextpad=0.5, columnspacing=1.5)
    valid_counts = summary.loc[summary.n_valid.gt(0)].groupby("temperature").n_valid.agg(["min", "max"])
    count_text = "; ".join(
        f"{int(t)} °C: {int(row['min'])}" + (f"–{int(row['max'])}" if row['min'] != row['max'] else "")
        for t, row in valid_counts.iterrows())
    figure_caption = [
        f"cgmaptools: {bin_size} bp bins of all input CG, CHG and CHH sites; sample bin means are averaged within temperature.",
        "Curves use shape-preserving PCHIP interpolation through bin means; empty bins remain gaps. No coverage weighting.",
        "Valid samples per nonempty bin: " + count_text + ".",
        "Zeros are retained; missing values are omitted without imputation. Sample SD and valid n are in the source table.",
        "Negative-strand RCA follows transcription direction; display widths = 1:2:1. Y-axis ranges differ between panels.",
        f"Peak bin (BED): {CHROMOSOME}:[{peak_start}, {peak_end}); "
        f"TSS upstream {(peak_start-GENE_END)/1000:.1f}–{(peak_end-GENE_END)/1000:.1f} kb.",
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
    parser.add_argument("--cgmaptools", type=Path, default=DEFAULT_CGMAPTOOLS)
    parser.add_argument("--bin-size", type=int, default=DEFAULT_BIN_SIZE)
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args(argv)
    paths = {c: getattr(args, f"{c.lower()}_positions") or args.positions_dir / f"{c}.CV.before.numbered.tsv"
             for c in CONTEXTS}
    data, sites, raw, counts, ignored = load_data(args.input, paths)
    site_summary = summarize(raw, sites)
    sample_counts = {int(t): int(n) for t, n in data.groupby("temperature").size().items()}
    if args.validate_only:
        print(json.dumps(dict(samples=sample_counts, sites=counts, raw_rows=len(raw), site_means=len(site_summary)), indent=2))
        return 0
    auditor = load_auditor(args.qa_tools)
    if importlib.util.find_spec("fitz") is None:
        raise RuntimeError("PDF collision QA requires PyMuPDF (fitz)")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    qa_dir = args.output_dir / f"{PREFIX}_QA"
    qa_dir.mkdir(exist_ok=True)
    stem = args.output_dir / PREFIX
    sites, bins, bin_replicates, summary = bin_with_cgmaptools(
        raw, sites, args.cgmaptools, args.bin_size, qa_dir / "cgmaptools")
    run_audit([sys.executable, str(args.qa_tools / "validate_figure.py"), __file__], qa_dir / "source-preflight.txt")
    fig, records = make_figure(summary, auditor, qa_dir, args.bin_size)
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
    for name, table in [("sites", sites), ("replicates", raw), ("site_summary", site_summary),
                        ("bins", bins), ("bin_replicates", bin_replicates), ("summary", summary)]:
        table.to_csv(stem.with_suffix(f".{name}.tsv"), sep="\t", index=False, na_rep="NA", float_format="%.12g")
    caption = """# RCA 温度分组甲基化分布 / Temperature-grouped RCA methylation profiles

Panels (a)–(c) show CG, CHG and CHH. cgmaptools mbin calculates the arithmetic
mean of all input site percentages within each sample and bin. These sample-level
bin means are then averaged with equal sample weights within each temperature.
Curves use shape-preserving PCHIP interpolation through bin-center means for display;
no point markers are shown. This does not alter the bin estimates or add observations.
No extrapolation, coverage weighting or significance tests are applied. No uncertainty band is shown;
sample SD (ddof=1), valid sample n, and missing n are in the summary table.
Zero is measured, not missing; blank/NA observations are omitted without imputation.
Empty bins break the curves; SD is NA for fewer than two valid sample bin means.
Y-axis limits and ticks are CG 0–100% (step 20), CHG 0–60% (step 20),
and CHH 0–30% (step 10).

RCA (AT2G39730), NC_003071.7:[16570745,16573692), is on the negative strand.
BED coordinates are 0-based and half-open. Bin midpoints anchor the interpolated curves.
The left-to-right axis follows transcription, from high to low genomic coordinates:
upstream [16573692,16575692), gene body [16570745,16573692), downstream
[16568745,16570745). The body is scaled to 0–100%; the flanks each span 2 kb.
Display widths are 1:2:1 and thus are not a uniform genomic-distance scale.
The CHH plot area is 70 mm high; CG and CHG are each 33 mm high.
The caption reports the upstream bin with the highest CHH mean across temperatures.

本图使用输入表中的全部甲基化位点，不额外进行位点过滤；不代表输入表之外的胞嘧啶。
分箱分别从三个区域的转录方向起点开始，不跨越 TSS/TES；每区末尾不足整箱时保留短箱。
曲线以 PCHIP 保形插值连接 bin 均值，仅改善显示，不改变统计量或制造波动。
空箱不补零、不跨空箱插值、不外推；单个孤立有效箱无法形成无标记曲线，仍完整保留在数据表中。
样本行作为重复单位；不从样本名称推断其为生物学重复或技术重复。

cgmaptools adapter: the input table contains percentages but no read counts.
Each percentage/100 is encoded losslessly as an integer numerator/denominator
for mbin input. These computational ratios are NOT sequencing read counts or
coverage. No coverage inference or selection is made; -c 1 retains every encoded
value, including zero. Region-local coordinates follow transcription direction.
mbin rounds ratios to four decimals (0.01 percentage-point resolution).
Every sample/bin result is independently checked against the input percentages
within that rounding tolerance. Adapter inputs, commands and raw mbin outputs
are retained in the QA folder; use sites/bins TSVs for original genomic positions.
"""
    caption += f"\nBin width: {args.bin_size} bp; bins per context: {len(bins)}.\n"
    chh_upstream = summary[(summary.context == "CHH") & (summary.region == "upstream")]
    peak_id = chh_upstream.groupby("bin_id").mean_pct.mean().idxmax()
    peak_table = chh_upstream[chh_upstream.bin_id == peak_id]
    peak_table.to_csv(stem.with_suffix(".CHH_peak_region.tsv"), sep="\t", index=False, na_rep="NA", float_format="%.12g")
    peak_start, peak_end = int(peak_table.iloc[0].start), int(peak_table.iloc[0].end)
    caption += (f"\nCHH peak bin: {CHROMOSOME}:[{peak_start},{peak_end}) (BED); "
                f"1-based inclusive: {CHROMOSOME}:{peak_start+1}–{peak_end}; "
                f"TSS upstream {peak_start-GENE_END}–{peak_end-GENE_END} bp. "
                "This is the peak-containing bin, not a precisely delimited biological peak.\n")
    caption += "\nSample rows per temperature: " + ", ".join(f"{t} °C: n={n}" for t, n in sample_counts.items()) + ".\n"
    missing_samples = raw.groupby("sample").missing.all()
    missing_sample_ids = missing_samples.index[missing_samples].tolist()
    caption += f"\nInput methylation table: {args.input.name}; {len(sites)} sites; {len(data)} sample rows.\n"
    caption += ("\nSamples with all methylation measurements missing: " +
                (", ".join(missing_sample_ids) or "none") +
                ". Retained in source tables; omitted from means without imputation. "
                "The figure reports valid sample counts in nonempty bins.\n")
    stem.with_suffix(".caption.md").write_text(caption, encoding="utf-8")
    report = dict(status="automated_QA_PASS", visual_review="pending", panels=records,
                  samples=sample_counts, sample_ids=data.Temperature.tolist(), site_counts=counts,
                  region_counts=sites.groupby(["context", "region"], observed=True).size().rename("n").reset_index().to_dict("records"),
                  source_rows=len(data), raw_rows=len(raw), mean_rows=len(summary),
                  missing_measurements=int(raw.missing.sum()), ignored_non_methylation_columns=ignored,
                  wholly_missing_samples=missing_sample_ids,
                  exclusions="No methylation columns or input sample rows excluded; unused coordinate rows lack input-table columns. Missing values are omitted only in numerical aggregation.",
                  coordinate_system="BED 0-based half-open", chromosome=CHROMOSOME,
                  gene_start=GENE_START, gene_end=GENE_END, strand="-", flank_bp=FLANK,
                  aggregation="cgmaptools equal-weight site mean per sample/bin, then equal-weight sample mean per temperature; SD ddof=1",
                  curve_display="PCHIP through consecutive finite bin centers; no markers, extrapolation or bridging empty bins",
                  bin_size_bp=args.bin_size, bins_per_context=len(bins),
                  bin_sample_rows=len(bin_replicates), empty_sample_bins=int(bin_replicates.methylation_pct.isna().sum()),
                  cgmaptools=dict(path=str(args.cgmaptools), subcommand="mbin", min_encoded_denominator=1,
                                  adapter="exact rational encoding of percentages, not read counts", rounding_tolerance_pct=0.005001),
                  inputs={name: dict(path=str(path.resolve()), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
                          for name, path in {"methylation": args.input, **paths}.items()},
                  figure_mm=[FIG_WIDTH_MM, FIG_HEIGHT_MM], raster_dpi=DEFAULT_DPI,
                  display_width_ratio=[1, BODY_DISPLAY_WIDTH, 1],
                  chh_peak_bin=dict(chromosome=CHROMOSOME, start_bed=peak_start, end_bed=peak_end,
                                    upstream_distance_bp=[peak_start-GENE_END, peak_end-GENE_END]),
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
