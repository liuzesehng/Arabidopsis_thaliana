#!/usr/bin/env python3
"""Regional CHH methylation versus temperature, with exact Spearman inference.

BED site coordinates are zero-based, half-open; the requested region is
one-based and inclusive. For a single-base site, position_1based = start + 1
= end. Average sites within each sample first, then summarize sample means
by temperature. Empty cells are missing,
not zero. Lines connecting means describe trends only.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import itertools
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "chh-region-mpl"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.transforms import ScaledTranslation
import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata

PROJECT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = PROJECT / "list/meth_experiment"
DEFAULT_INPUT = DEFAULT_OUTPUT / "Alt.all.snp_meth.tsv"
DEFAULT_BED = DEFAULT_OUTPUT / "CHH.CV.before.numbered.tsv"
DEFAULT_QA_TOOLS = Path("/datapool/home/2023102768/.codex/skills/nature-figure/scripts")
# Requested 1-based closed interval: [16574892 - 150, 16575292 + 150].
CHROM, START, END = "NC_003071.7", 16574742, 16575442
TEMPERATURES = tuple(range(4, 33, 4))
PREFIX = f"CHH_region_{START}_{END}_temperature"
BLUE = "#2F6DAE"


def load_sites(path: Path) -> pd.DataFrame:
    """Read ID/chrom/BED-start/BED-end; select single-base sites in 1-based space."""
    records, seen_ids, seen_positions = [], set(), set()
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for line, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if len(row) < 4:
                raise ValueError(f"coordinate row {line}: expected ID, chromosome, start, end")
            feature, chrom = row[:2]
            try:
                start, end = int(row[2]), int(row[3])
            except ValueError as exc:
                raise ValueError(f"coordinate row {line}: invalid BED coordinates") from exc
            if not feature.startswith("CHH_") or not chrom or start < 0 or end <= start:
                raise ValueError(f"coordinate row {line}: invalid CHH ID or BED interval")
            if end != start + 1:
                raise ValueError(f"coordinate row {line}: expected a single-base BED site (end = start + 1)")
            if feature in seen_ids or (chrom, start, end) in seen_positions:
                raise ValueError(f"coordinate row {line}: duplicate site ID or interval")
            seen_ids.add(feature)
            seen_positions.add((chrom, start, end))
            position_1based = start + 1
            if chrom == CHROM and START <= position_1based <= END:
                records.append(dict(feature=feature, chromosome=chrom, start=start, end=end,
                                    position_1based=position_1based))
    if not records:
        raise ValueError("no CHH sites fall within the requested 1-based inclusive region")
    return pd.DataFrame(records).sort_values(["start", "end", "feature"]).reset_index(drop=True)


def load_data(path: Path, sites: pd.DataFrame) -> pd.DataFrame:
    features = sites.feature.tolist()
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows or not rows[0] or any(not h.strip() for h in rows[0]):
        raise ValueError("missing or blank column header")
    if len(set(rows[0])) != len(rows[0]):
        raise ValueError("duplicate column names")
    if any(len(row) != len(rows[0]) for row in rows[1:]):
        raise ValueError("inconsistent number of TSV fields")
    absent = set(["Temperature", *features]) - set(rows[0])
    if absent:
        raise ValueError(f"missing required columns: {', '.join(sorted(absent))}")
    data = pd.DataFrame(rows[1:], columns=rows[0])[["Temperature", *features]].copy()
    data = data.rename(columns={"Temperature": "sample"})
    parsed = data["sample"].str.extract(r"^T(\d+)-(\d+)$")
    if parsed.isna().any().any() or data["sample"].duplicated().any():
        raise ValueError("sample IDs must be unique and match T<temperature>-<replicate>")
    data["temperature"] = parsed[0].astype(int)
    data["replicate"] = parsed[1].astype(int)
    if (data.replicate < 1).any() or data.duplicated(["temperature", "replicate"]).any():
        raise ValueError("replicate numbers must be positive and unique per temperature")
    counts = data.groupby("temperature").size().to_dict()
    if counts != {t: 3 for t in TEMPERATURES}:
        raise ValueError(f"expected three sample rows per temperature at 4–32 °C; got {counts}")
    for feature in features:
        raw = data[feature].str.strip()
        numeric = pd.to_numeric(raw.replace("", np.nan), errors="coerce").astype(float)
        if ((raw != "") & numeric.isna()).any() or np.isinf(numeric).any():
            raise ValueError(f"{feature}: invalid numeric value; only empty fields denote missingness")
        if ((numeric < 0) | (numeric > 100)).any():
            raise ValueError(f"{feature}: methylation must be between 0 and 100 percent")
        data[feature] = numeric
    return data.sort_values(["temperature", "replicate"], kind="stable").reset_index(drop=True)


def summarize(data: pd.DataFrame, sites: pd.DataFrame):
    values = data[sites.feature.tolist()]
    points = data[["sample", "temperature", "replicate"]].copy()
    points["n_sites_total"] = len(sites)
    points["n_sites_valid"] = values.count(axis=1)
    points["n_sites_missing"] = len(sites) - points.n_sites_valid
    points["methylation"] = values.mean(axis=1)
    points["missing"] = points.n_sites_valid.eq(0)
    points["partial_missing"] = points.n_sites_missing.gt(0) & ~points.missing
    points["unit"] = "%"
    points["x_plot"] = points.temperature.astype(float)
    groups = []
    for temperature, block in points.groupby("temperature", sort=True):
        valid = block.methylation.count()
        if valid < 2:
            raise ValueError(f"{temperature} °C: at least two valid regional replicates required for sample SD")
        # Deterministic display offsets never enter summaries or correlation.
        points.loc[block.index, "x_plot"] = temperature + np.linspace(-0.55, 0.55, len(block))
        groups.append(dict(temperature=int(temperature), n_valid=int(valid),
                           n_missing=int(block.missing.sum()), n_sites=len(sites),
                           mean=float(block.methylation.mean()),
                           sd=float(block.methylation.std(ddof=1)), unit="%"))
    return points, pd.DataFrame(groups)


def exact_spearman(means):
    """Two-sided exhaustive pairing permutation test, including tied ranks."""
    means = np.asarray(means, dtype=float)
    if means.shape != (8,) or not np.isfinite(means).all():
        raise ValueError("Spearman requires eight finite temperature means")
    x, y = rankdata(TEMPERATURES) - 4.5, rankdata(means) - 4.5
    denominator = np.linalg.norm(x) * np.linalg.norm(y)
    result = dict(n_groups=8, total_permutations=40320,
                  test="two-sided exhaustive pairing permutation Spearman",
                  correlation_unit="temperature means", adjustment="none")
    if denominator == 0:
        return dict(result, rho=np.nan, p=np.nan, extreme_permutations=None, status="constant_variable")
    permutations = np.array(list(itertools.permutations(range(8))), dtype=np.int16)
    observed = float(x @ y)
    null = np.einsum("ij,j->i", x[permutations], y)
    extreme = int(np.count_nonzero(np.abs(null) >= abs(observed) - 1e-12))
    return dict(result, rho=float(np.clip(observed / denominator, -1, 1)),
                p=extreme / len(permutations), extreme_permutations=extreme, status="ok")


def load_auditor(directory):
    spec = importlib.util.spec_from_file_location("chh_region_alignment", directory / "audit_panel_alignment.py")
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load panel alignment auditor")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def make_figure(points, groups, statistics, auditor, qa_dir):
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
                         "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 8,
                         "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.6, "pdf.fonttype": 42, "svg.fonttype": "none",
                         "savefig.facecolor": "white", "legend.frameon": False})
    valid = points.loc[~points.missing]
    lower = min(float((groups["mean"] - groups.sd).min()), float(valid.methylation.min()))
    upper = max(float((groups["mean"] + groups.sd).max()), float(valid.methylation.max()))
    if lower < 0 or upper > 25:
        raise ValueError("data or mean ± SD exceeds the specified 0–25% axis; export stopped to avoid clipping")
    fig = plt.figure(figsize=(89 / 25.4, 85 / 25.4))
    ax = fig.add_subplot(111)
    fig.subplots_adjust(left=0.17, right=0.965, bottom=0.40, top=0.86)
    ax.set_gid("a")
    ax.axvspan(22, 26, color="#E8E8E8", alpha=0.75, lw=0, zorder=0)
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#E3E3E3", lw=0.4)
    ax.scatter(valid.x_plot, valid.methylation, s=13, color=BLUE, alpha=0.40,
               edgecolors="none", zorder=3)
    ax.plot(groups.temperature, groups["mean"], color=BLUE, lw=0.8, zorder=2)
    ax.errorbar(groups.temperature, groups["mean"], yerr=groups.sd, fmt="o", ms=3.5,
                color=BLUE, ecolor=BLUE, elinewidth=0.75, capsize=2.5, capthick=0.75,
                markeredgecolor="white", markeredgewidth=0.4, zorder=4)
    ax.set(xlim=(2, 34), ylim=(0, 25), xticks=TEMPERATURES, yticks=np.arange(0, 26, 5),
           xlabel="Temperature (°C)", ylabel="CHH methylation (%)")
    ax.tick_params(length=2.5, width=0.6)
    ax.set_title("Regional CHH methylation", pad=9, fontweight="bold")
    ax.text(0, 1, "(a)", fontsize=8, fontweight="bold", ha="left", va="bottom",
            transform=ax.transAxes + ScaledTranslation(-23 / 72, 7 / 72, fig.dpi_scale_trans))
    if statistics["status"] == "ok":
        p = statistics["p"]
        p_text = f"{p:.3f}" if p >= 0.001 else f"{p:.2e}"
        annotation = f"ρ = {statistics['rho']:.2f}; P = {p_text}"
    else:
        annotation = "ρ = NA; P = NA (constant)"
    ax.text(0.04, 0.94, annotation, transform=ax.transAxes, va="top", fontsize=7)
    n_two = ", ".join(str(int(t)) for t in groups.loc[groups.n_valid.eq(2), "temperature"])
    count_line = f"n = 2 at {n_two} °C; n = 3 otherwise." if n_two else "n = 3 at each temperature."
    notes = [f"{CHROM}: {START}–{END} (1-based)",
             f"{int(groups.n_sites.iloc[0])} CHH sites; equal-weight mean within each sample.",
             "Faint points: samples; solid points: mean ± SD.",
             count_line, "Lines describe trends; gray band marks 24 °C.",
             "Spearman: 8 means; two-sided exact P (unadjusted)."]
    for i, note in enumerate(notes):
        fig.text(0.065, (4 + (len(notes) - 1 - i) * 3.25) / 85, note, fontsize=5.5, va="bottom")
    fig.canvas.draw()
    auditor.require_matplotlib_panel_alignment(fig, json_out=qa_dir / "alignment.json", strict=True)
    return fig, dict(panel="(a)", plotted_replicates=len(valid), means=len(groups),
                     error_bars=len(groups), data_and_sd_extent=[lower, upper], ylim=[0, 25])


def run_audit(command, output):
    result = subprocess.run(command, capture_output=True, text=True)
    output.write_text(result.stdout + result.stderr, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(f"audit failed: {output.name}\n{result.stdout[-2500:]}\n{result.stderr[-1000:]}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--bed", type=Path, default=DEFAULT_BED)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--qa-tools", type=Path, default=DEFAULT_QA_TOOLS)
    parser.add_argument("--validate-only", action="store_true", help="validate data without creating outputs")
    args = parser.parse_args(argv)
    sites = load_sites(args.bed)
    data = load_data(args.input, sites)
    points, groups = summarize(data, sites)
    statistics = exact_spearman(groups["mean"])
    if args.validate_only:
        print(f"PASS: {len(sites)} sites, {len(points)} sample rows, {int((~points.missing).sum())} valid replicates")
        print(groups.to_string(index=False))
        print(statistics)
        return 0
    auditor = load_auditor(args.qa_tools)
    if importlib.util.find_spec("fitz") is None:
        raise RuntimeError("NOT AUDITABLE: PDF collision QA requires PyMuPDF (fitz)")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    qa_dir = args.output_dir / f"{PREFIX}_QA"
    qa_dir.mkdir(exist_ok=True)
    run_audit([sys.executable, str(args.qa_tools / "validate_figure.py"), __file__], qa_dir / "source-preflight.txt")
    fig, panel = make_figure(points, groups, statistics, auditor, qa_dir)
    stem = args.output_dir / PREFIX
    try:
        fig.savefig(stem.with_suffix(".pdf"))
        run_audit([sys.executable, str(args.qa_tools / "audit_pdf_text.py"), str(stem.with_suffix(".pdf")),
                   "--min-pt", "5", "--json"], qa_dir / "pdf-text.json")
        run_audit([sys.executable, str(args.qa_tools / "audit_figure_collisions.py"), str(stem.with_suffix(".pdf")),
                   "--strict", "--json-out", str(qa_dir / "collision-audit.json")], qa_dir / "collision-audit.txt")
        fig.savefig(stem.with_suffix(".svg"))
        fig.savefig(stem.with_suffix(".png"), dpi=600)
        fig.savefig(stem.with_suffix(".tiff"), dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    finally:
        plt.close(fig)
    for name, table in [("sites", sites), ("replicates", points), ("summary", groups),
                        ("statistics", pd.DataFrame([statistics]))]:
        table.to_csv(stem.with_suffix(f".{name}.tsv"), sep="\t", index=False, na_rep="NA", float_format="%.12g")
    counts = ", ".join(f"{r.temperature} °C: n = {r.n_valid}" for r in groups.itertuples())
    caption = ("# Regional CHH methylation and temperature\n\n"
               f"(a) Region {CHROM}:{START}–{END} in 1-based coordinates, including both endpoints, "
               f"including {len(sites)} CHH sites. The interval extends from 16574892 − 150 to 16575292 + 150. "
               "The numbered coordinate file contains CHH ID, chromosome, BED start (0-based), and BED end "
               "in columns 1–4. All sites are single-base intervals; each 1-based site position is BED start + 1, "
               "equal to BED end. The site table retains original BED start/end and adds position_1based. "
               "Each sample value is the equal-weight arithmetic mean of its available site percentages; "
               "no read-depth weighting or percentage rescaling is applied. "
               "Faint points show individual sample means with display-only horizontal offsets. "
               "Solid points and error bars show the mean and sample SD (ddof=1) across replicates. "
               "Connecting lines describe trends only, without fitting a model. "
               "The gray band marks 24 °C.\n\n"
               f"Valid replicate counts: {counts}. Empty measurements are missing, never zero; "
               "samples with no measured sites are omitted from summaries and plotting, but retained in the replicate table. "
               "Partially missing samples use available sites and are flagged with their valid site counts. "
               "Sample rows are treated as independent experimental replicates; biological versus technical status "
               "is not specified by the input table.\n\n"
               "Spearman correlation uses eight temperature means. The two-sided exact P value enumerates "
               "all 8! = 40,320 pairings, including tied-rank multiplicities. No multiple-testing correction "
               "is applied because one regional correlation is tested. These descriptive associations do not establish causality.\n")
    stem.with_suffix(".caption.md").write_text(caption, encoding="utf-8")
    report = dict(status="automated_QA_PASS", panel=panel,
                  region=dict(chromosome=CHROM, start=START, end=END, coordinate_system="1-based inclusive"),
                  equivalent_bed_region=dict(chromosome=CHROM, start=START - 1, end=END),
                  coordinate_convention="single-base BED site: position_1based = column 3 + 1 = column 4; START <= position_1based <= END",
                  input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
                  bed_sha256=hashlib.sha256(args.bed.read_bytes()).hexdigest(),
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  n_sites=len(sites), source_sample_rows=len(points),
                  missing_samples=points.loc[points.missing, "sample"].tolist(),
                  partial_missing_samples=points.loc[points.partial_missing, "sample"].tolist(),
                  excluded_columns="Features outside the selected 1-based inclusive region and non-methylation columns.",
                  correlation_unit="eight temperature means", alignment="not applicable: single panel",
                  visual_review="pending", versions=dict(python=sys.version.split()[0], numpy=np.__version__,
                  pandas=pd.__version__, scipy=scipy.__version__, matplotlib=matplotlib.__version__))
    (qa_dir / "summary.json").write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(pd.DataFrame([statistics]).to_string(index=False))
    print(f"PASS: {len(sites)} sites; {panel['plotted_replicates']} points; 8 means/SD; outputs: {stem}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, RuntimeError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
