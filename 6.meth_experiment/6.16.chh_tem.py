#!/usr/bin/env python3
"""Regional CG, CHG and CHH methylation versus temperature, with exact Spearman.

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
CONTEXTS = ("CG", "CHG", "CHH")
DEFAULT_BEDS = {context: DEFAULT_OUTPUT / f"{context}.CV.before.numbered.tsv" for context in CONTEXTS}
DEFAULT_QA_TOOLS = Path("/datapool/home/2023102768/.agents/skills/nature-figure/scripts")
# Requested 1-based closed interval: [16574892 - 150, 16575292 + 150].
CHROM, START, END = "NC_003071.7", 16574742, 16575442
TEMPERATURES = tuple(range(4, 33, 4))
PREFIX = f"CG_CHG_CHH_region_{START}_{END}_temperature"
BLUE = "#225E7A"
REPLICATE_COLOR = "#85AEC5"
PANEL_WIDTH_MM = 75
FIG_WIDTH_MM = 225
FIG_HEIGHT_MM = 95
DEFAULT_DPI = 600
Y_SCALES = {"CG": (80, 100, 5), "CHG": (35, 70, 5), "CHH": (5, 25, 5)}


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


def load_sites(path: Path, context: str = "CHH") -> pd.DataFrame:
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
            if not feature.startswith(f"{context}_") or not chrom or start < 0 or end <= start:
                raise ValueError(f"coordinate row {line}: invalid {context} ID or BED interval")
            if end != start + 1:
                raise ValueError(f"coordinate row {line}: expected a single-base BED site (end = start + 1)")
            if feature in seen_ids or (chrom, start, end) in seen_positions:
                raise ValueError(f"coordinate row {line}: duplicate site ID or interval")
            seen_ids.add(feature)
            seen_positions.add((chrom, start, end))
            position_1based = start + 1
            if chrom == CHROM and START <= position_1based <= END:
                records.append(dict(context=context, feature=feature, chromosome=chrom, start=start, end=end,
                                    position_1based=position_1based))
    if not records:
        raise ValueError(f"no {context} sites fall within the requested 1-based inclusive region")
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
        points.loc[block.index, "x_plot"] = temperature + np.linspace(-0.42, 0.42, len(block))
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


def analyze_contexts(data, sites):
    """Keep each context's within-sample averaging independent of other sites."""
    point_tables, group_tables, tests = [], [], []
    for context in CONTEXTS:
        selected = sites.loc[sites.context.eq(context)]
        points, groups = summarize(data, selected)
        points.insert(0, "context", context)
        groups.insert(0, "context", context)
        point_tables.append(points)
        group_tables.append(groups)
        tests.append(dict(context=context, n_sites=len(selected),
                          n_valid_samples=int((~points.missing).sum()),
                          **exact_spearman(groups["mean"])))
    statistics = pd.DataFrame(tests)
    # Report the original exact P in panels; provide BH q for the three tests.
    p = statistics.p.to_numpy()
    q = np.full(len(p), np.nan)
    finite = np.flatnonzero(np.isfinite(p))
    order = finite[np.argsort(p[finite], kind="stable")]
    if len(order):
        adjusted = p[order] * len(p) / np.arange(1, len(order) + 1)
        q[order] = np.minimum(1, np.minimum.accumulate(adjusted[::-1])[::-1])
    statistics["q_bh"] = q
    statistics["adjustment"] = "panel P unadjusted; q_bh across three contexts"
    return pd.concat(point_tables, ignore_index=True), pd.concat(group_tables, ignore_index=True), statistics


def load_auditor(directory):
    if not Path(directory / 'audit_panel_alignment.py').is_file():
        _skip_optional_qa(f"Optional QA tool not found: {directory / 'audit_panel_alignment.py'}")
        return None
    spec = importlib.util.spec_from_file_location("chh_region_alignment", directory / "audit_panel_alignment.py")
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load panel alignment auditor")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    try:
        spec.loader.exec_module(module)
    except ImportError as exc:
        _skip_optional_qa(f"Optional alignment dependency unavailable: {exc}")
        return None
    return module


def make_figure(points, groups, statistics, auditor, qa_dir):
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
                         "font.size": 6.6, "axes.labelsize": 7, "axes.titlesize": 8.2,
                         "xtick.labelsize": 6.3, "ytick.labelsize": 6.3,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.75, "pdf.fonttype": 42, "svg.fonttype": "none",
                         "savefig.facecolor": "white", "legend.frameon": False})
    fig = plt.figure(figsize=(FIG_WIDTH_MM / 25.4, FIG_HEIGHT_MM / 25.4))
    # Keep panel width fixed; increase height to reduce the plot width/height ratio.
    plot_width_mm = PANEL_WIDTH_MM * (0.965 - 0.17)
    gutter_mm = PANEL_WIDTH_MM - plot_width_mm
    grid = fig.add_gridspec(1, 3, left=PANEL_WIDTH_MM * 0.17 / FIG_WIDTH_MM,
                           right=1 - PANEL_WIDTH_MM * 0.035 / FIG_WIDTH_MM,
                           bottom=0.40, top=0.86, wspace=gutter_mm / plot_width_mm)
    labels, panels = [], []
    for index, context in enumerate(CONTEXTS):
        valid = points.loc[points.context.eq(context) & ~points.missing]
        sg = groups.loc[groups.context.eq(context)]
        st = statistics.loc[statistics.context.eq(context)].iloc[0]
        lower = min(float((sg["mean"] - sg.sd).min()), float(valid.methylation.min()))
        upper = max(float((sg["mean"] + sg.sd).max()), float(valid.methylation.max()))
        ymin, tick_max, step = Y_SCALES[context]
        # Reserve annotation space as in the reference, relative to the displayed lower bound.
        ymax = max(tick_max, ymin + (upper - ymin) * 1.28)
        if lower < ymin:
            raise ValueError(f"{context}: data or mean ± SD falls below {ymin}%; export stopped")
        ax = fig.add_subplot(grid[0, index])
        ax.set_gid(chr(97 + index))
        ax.axvspan(22, 26, color="#E9E9E9", alpha=0.9, lw=0, zorder=0)
        ax.set_axisbelow(True)
        ax.scatter(valid.x_plot, valid.methylation, s=18, color=REPLICATE_COLOR, alpha=0.78,
                   edgecolors="white", linewidths=0.45, zorder=3)
        ax.plot(sg.temperature, sg["mean"], color=BLUE, lw=1.15, alpha=0.95, zorder=4)
        ax.errorbar(sg.temperature, sg["mean"], yerr=sg.sd, fmt="D", ms=4.0,
                    color=BLUE, ecolor=BLUE, elinewidth=0.9, capsize=2.5, capthick=0.9,
                    markeredgecolor="white", markeredgewidth=0.55, zorder=5)
        ax.set(xlim=(2, 34), ylim=(ymin, ymax), xticks=TEMPERATURES, yticks=np.arange(ymin, min(tick_max, 100) + 1, step),
               xlabel="Temperature (°C)", ylabel=f"{context} methylation (%)")
        ax.tick_params(length=2.7, width=0.65, direction="out")
        ax.set_title(context, pad=9, fontweight="semibold")
        labels.append(ax.text(0, 1, f"({chr(97 + index)})", fontsize=8, fontweight="bold", ha="left", va="bottom",
                              transform=ax.transAxes + ScaledTranslation(-23 / 72, 7 / 72, fig.dpi_scale_trans)))
        if st.status == "ok":
            p_text = f"{st.p:.4f}" if st.p >= 0.0001 else f"{st.p:.2e}"
            rho_text = f"{st.rho:.3f}".replace("-", "\N{MINUS SIGN}")
            annotation = f"ρ = {rho_text}, P = {p_text}"
        else:
            annotation = "ρ = NA; P = NA (constant)"
        ax.text(0.025, 0.965, annotation, transform=ax.transAxes,
                ha="left", va="top", fontsize=5.8, color="#222222", zorder=6)
        panels.append(dict(panel=f"({chr(97 + index)})", context=context, n_sites=int(sg.n_sites.iloc[0]),
                           plotted_replicates=len(valid), means=len(sg), error_bars=len(sg),
                           data_and_sd_extent=[lower, upper], ylim=[ymin, ymax],
                           ytick_max=min(tick_max, 100), annotation_axes_fraction=[0.025, 0.965],
                           uncertainty="sample SD (ddof=1)"))
    count_descriptions = []
    for context in CONTEXTS:
        sg = groups.loc[groups.context.eq(context)]
        n_two = ", ".join(str(int(t)) for t in sg.loc[sg.n_valid.eq(2), "temperature"])
        count_descriptions.append(f"n = 2 at {n_two} °C; n = 3 otherwise." if n_two else "n = 3 at each temperature.")
    count_line = count_descriptions[0] if len(set(count_descriptions)) == 1 else "Replicate counts by context and temperature are reported in the caption."
    notes = [f"{CHROM}: {START}–{END} (1-based, inclusive); equal-weight site mean within each sample.",
             "Light circles: samples; dark diamonds: mean ± sample SD. " + count_line,
             "Lines describe trends only; gray bands mark 24 °C. Note the different y-axis scales.",
             "Spearman: 8 temperature means; two-sided exact P (unadjusted). BH q values are provided in the statistics table."]
    for i, note in enumerate(notes):
        fig.text(0.035, (4 + (len(notes) - 1 - i) * 3.7) / FIG_HEIGHT_MM, note, fontsize=6, va="bottom")
    fig.canvas.draw()
    if auditor is None:
        _skip_optional_qa("audit_panel_alignment.py is unavailable", qa_dir / 'alignment.json')
    else:
        auditor.require_matplotlib_panel_alignment(fig, tolerance_pt=1.5, gutter_tolerance_pt=1.5, strict=True)
        manifest = auditor.matplotlib_layout_manifest(fig)
        for panel, label in zip(manifest["panels"], labels):
            anchor = label.get_transform().transform(label.get_position())
            panel["panel_label"] = label.get_text().strip("()")
            panel["panel_label_anchor_pt"] = list(map(float, fig.dpi_scale_trans.inverted().transform(anchor) * 72))
        report = auditor.audit_layout_manifest(manifest, tolerance_pt=1.5, gutter_tolerance_pt=1.5, require_panel_labels=True)
        auditor.write_json_report(report, qa_dir / "alignment.json")
        if auditor.exit_code(report, strict=True):
            raise RuntimeError(auditor.render_text(report, strict=True))
    return fig, panels


def run_audit(command, output):
    if _skip_missing_qa_command(command, output):
        return
    result = subprocess.run(command, capture_output=True, text=True)
    output.write_text(result.stdout + result.stderr, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(f"audit failed: {output.name}\n{result.stdout[-2500:]}\n{result.stderr[-1000:]}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--bed", "--chh-bed", dest="chh_bed", type=Path, default=DEFAULT_BED,
                        help="numbered CHH BED table (--bed retained as an alias)")
    parser.add_argument("--cg-bed", type=Path, default=DEFAULT_BEDS["CG"])
    parser.add_argument("--chg-bed", type=Path, default=DEFAULT_BEDS["CHG"])
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--qa-tools", type=Path, default=DEFAULT_QA_TOOLS)
    parser.add_argument("--validate-only", action="store_true", help="validate data without creating outputs")
    args = parser.parse_args(argv)
    bed_paths = {"CG": args.cg_bed, "CHG": args.chg_bed, "CHH": args.chh_bed}
    sites = pd.concat([load_sites(bed_paths[c], c) for c in CONTEXTS], ignore_index=True)
    data = load_data(args.input, sites)
    points, groups, statistics = analyze_contexts(data, sites)
    if args.validate_only:
        print(f"PASS: {len(sites)} sites, {len(data)} sample rows, {int((~points.missing).sum())} valid context/sample means")
        print(groups.to_string(index=False))
        print(statistics)
        return 0
    auditor = load_auditor(args.qa_tools)
    # PDF QA dependencies are checked only when the optional tools run.
    args.output_dir.mkdir(parents=True, exist_ok=True)
    qa_dir = args.output_dir / f"{PREFIX}_QA"
    qa_dir.mkdir(exist_ok=True)
    run_audit([sys.executable, str(args.qa_tools / "validate_figure.py"), __file__], qa_dir / "source-preflight.txt")
    fig, panels = make_figure(points, groups, statistics, auditor, qa_dir)
    stem = args.output_dir / PREFIX
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
    for name, table in [("sites", sites), ("replicates", points), ("summary", groups),
                        ("statistics", statistics)]:
        table.to_csv(stem.with_suffix(f".{name}.tsv"), sep="\t", index=False, na_rep="NA", float_format="%.12g")
    counts = "\n\n".join(c + ": " + ", ".join(f"{r.temperature} °C: n = {r.n_valid}" for r in
                           groups.loc[groups.context.eq(c)].itertuples()) for c in CONTEXTS)
    panel_description = "; ".join(f"({chr(97+i)}) {c}, {int(sites.context.eq(c).sum())} sites" for i, c in enumerate(CONTEXTS))
    caption = ("# Regional CG, CHG and CHH methylation and temperature\n\n"
               f"Panels {panel_description}. Region {CHROM}:{START}–{END} uses 1-based coordinates, including both endpoints. "
               "The interval extends from 16574892 − 150 to 16575292 + 150. "
               "The numbered coordinate files contain site ID, chromosome, BED start (0-based), and BED end "
               "in columns 1–4. All sites are single-base intervals; each 1-based site position is BED start + 1, "
               "equal to BED end. The site table retains original BED start/end and adds position_1based. "
               "Each sample value is the equal-weight arithmetic mean of its available site percentages within its methylation context; "
               "no read-depth weighting or percentage rescaling is applied. "
               "Light blue circles show individual sample means with display-only horizontal offsets. "
               "Dark blue diamonds and error bars show the mean and sample SD (ddof=1) across replicates. "
               "Connecting lines describe trends only, without fitting a model. "
               "The gray bands mark 24 °C. Labeled y-axis ticks span 80–100% (CG), 35–70% (CHG), and 5–25% (CHH). "
               "The displayed upper limit leaves at least 28% headroom above the data/SD extent relative to the lower bound; "
               "CG has an unlabeled annotation area above 100%, with no ticks above 100%. "
               f"Each panel is {PANEL_WIDTH_MM} × {FIG_HEIGHT_MM:.3f} mm, with increased height at fixed width. "
               "Correlation annotations are inside the upper-left plot area at axes fraction (0.025, 0.965); "
               "visual slopes should not be compared directly across these different scales.\n\n"
               f"Valid replicate counts: {counts}. Empty measurements are missing, never zero; "
               "samples with no measured sites are omitted from summaries and plotting, but retained in the replicate table. "
               "Partially missing samples use available sites and are flagged with their valid site counts. "
               "Sample rows are treated as independent experimental replicates; biological versus technical status "
               "is not specified by the input table.\n\n"
               "Spearman correlation uses eight temperature means. The two-sided exact P value enumerates "
               "all 8! = 40,320 pairings, including tied-rank multiplicities, separately for each context. "
               "Panel annotations show unadjusted exact P values; BH q values across the three planned tests "
               "are provided in the statistics table. These descriptive associations do not establish causality.\n")
    stem.with_suffix(".caption.md").write_text(caption, encoding="utf-8")
    report = dict(status=("automated_QA_INCOMPLETE" if _OPTIONAL_QA_SKIPPED else "automated_QA_PASS"), panels=panels,
                  figure_contract=dict(question="How does regional mean methylation vary with temperature across CG, CHG and CHH?",
                                       archetype="quantitative grid", panel_order=list(CONTEXTS), backend="python",
                                       size_mm=[FIG_WIDTH_MM, FIG_HEIGHT_MM], panel_size_mm=[PANEL_WIDTH_MM, FIG_HEIGHT_MM],
                                       plot_area_mm=[PANEL_WIDTH_MM * 0.795, FIG_HEIGHT_MM * 0.46],
                                       reuse="structural adaptation"),
                  region=dict(chromosome=CHROM, start=START, end=END, coordinate_system="1-based inclusive"),
                  equivalent_bed_region=dict(chromosome=CHROM, start=START - 1, end=END),
                  coordinate_convention="single-base BED site: position_1based = column 3 + 1 = column 4; START <= position_1based <= END",
                  input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
                  bed_sha256={c: hashlib.sha256(p.read_bytes()).hexdigest() for c, p in bed_paths.items()},
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  n_sites=len(sites), site_counts={c: int(sites.context.eq(c).sum()) for c in CONTEXTS},
                  source_sample_rows=len(data), source_context_sample_rows=len(points),
                  missing_samples=points.loc[points.missing, ["context", "sample"]].to_dict("records"),
                  partial_missing_samples=points.loc[points.partial_missing, ["context", "sample"]].to_dict("records"),
                  excluded_columns="Features outside the selected 1-based inclusive region and non-methylation columns.",
                  correlation_unit="eight temperature means per context", alignment=("SKIPPED" if auditor is None else "PASS: three comparable panels, 1.5 pt tolerance"),
                  visual_review="pending", versions=dict(python=sys.version.split()[0], numpy=np.__version__,
                  pandas=pd.__version__, scipy=scipy.__version__, matplotlib=matplotlib.__version__))
    (qa_dir / "summary.json").write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(statistics[["context", "n_sites", "rho", "p", "q_bh"]].to_string(index=False))
    print(f"PASS: {len(sites)} sites; {sum(p['plotted_replicates'] for p in panels)} points; 24 means/SD; outputs: {stem}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, RuntimeError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
