#!/usr/bin/env python3
"""Four CHH sites versus temperature: replicates, mean ± SD and exact Spearman.

The eight temperature means are the correlation units. Blank methylation is
missing, never zero. Connecting lines describe trends without fitting a model.
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
import re
import subprocess
import sys
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rca-temperature-mpl"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.transforms import ScaledTranslation
import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata

PROJECT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = PROJECT / "list/meth_experiment/Alt.all.snp_meth.tsv"
DEFAULT_OUTPUT = PROJECT / "list/meth_experiment"
DEFAULT_QA_TOOLS = Path("/datapool/home/2023102768/.codex/skills/nature-figure/scripts")
FEATURES = ("CHH_1634", "CHH_1635", "CHH_1694", "CHH_1696")
TEMPERATURES = tuple(range(4, 33, 4))
PREFIX = "RCA_temperature_methylation"
BLUE = "#2F6DAE"
CONTROL = "#E8E8E8"


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


def load_data(path: Path) -> pd.DataFrame:
    """Validate the sample design and selected sites without changing source data."""
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows or not rows[0] or any(not h.strip() for h in rows[0]):
        raise ValueError("missing or blank column header")
    if len(set(rows[0])) != len(rows[0]):
        raise ValueError("duplicate column names")
    if any(len(row) != len(rows[0]) for row in rows[1:]):
        raise ValueError("inconsistent number of TSV fields")
    absent = set(("Temperature", *FEATURES)) - set(rows[0])
    if absent:
        raise ValueError(f"missing required columns: {', '.join(sorted(absent))}")
    data = pd.DataFrame(rows[1:], columns=rows[0])[["Temperature", *FEATURES]].copy()
    data = data.rename(columns={"Temperature": "sample"})
    parsed = data["sample"].str.extract(r"^T(\d+)-(\d+)$")
    if parsed.isna().any().any() or data["sample"].duplicated().any():
        raise ValueError("sample IDs must be unique and match T<temperature>-<replicate>")
    data["temperature"] = parsed[0].astype(int)
    data["replicate"] = parsed[1].astype(int)
    if (data["replicate"] < 1).any() or data.duplicated(["temperature", "replicate"]).any():
        raise ValueError("replicate numbers must be positive and unique per temperature")
    counts = data.groupby("temperature").size().to_dict()
    if counts != {t: 3 for t in TEMPERATURES}:
        raise ValueError(f"expected three sample rows per temperature at 4–32 °C; got {counts}")
    for feature in FEATURES:
        raw = data[feature].str.strip()
        numeric = pd.to_numeric(raw.replace("", np.nan), errors="coerce").astype(float)
        if ((raw != "") & numeric.isna()).any() or np.isinf(numeric).any():
            raise ValueError(f"{feature}: invalid numeric value; only empty fields denote missingness")
        if ((numeric < 0) | (numeric > 100)).any():
            raise ValueError(f"{feature}: methylation must be between 0 and 100 percent")
        data[feature] = numeric
        valid = data.groupby("temperature")[feature].count()
        if (valid == 0).any():
            raise ValueError(f"{feature}: whole temperature group missing: {valid[valid == 0].index.tolist()}")
        if (valid == 1).any():
            raise ValueError(f"{feature}: at least two valid replicates per temperature are required for sample SD")
    return data.sort_values(["temperature", "replicate"], kind="stable").reset_index(drop=True)


def exact_spearman(means, permutations):
    """Exhaustive two-sided pairing test, retaining tied-rank multiplicities."""
    y = np.asarray(means, dtype=float)
    if y.shape != (8,) or not np.isfinite(y).all():
        raise ValueError("Spearman requires eight finite temperature means")
    x = rankdata(TEMPERATURES) - 4.5
    y = rankdata(y) - 4.5
    denominator = np.linalg.norm(x) * np.linalg.norm(y)
    if denominator == 0:
        return np.nan, np.nan, None, "constant_variable"
    observed = float(x @ y)
    null = np.einsum("ij,j->i", x[permutations], y)
    extreme = int(np.count_nonzero(np.abs(null) >= abs(observed) - 1e-12))
    return float(np.clip(observed / denominator, -1, 1)), extreme / len(permutations), extreme, "ok"


def bh_adjust(pvalues):
    """BH over the four planned tests; undefined tests occupy a conservative slot."""
    p = np.asarray(pvalues, dtype=float)
    q = np.full_like(p, np.nan)
    valid = np.flatnonzero(np.isfinite(p))
    order = valid[np.argsort(p[valid], kind="stable")]
    if len(order):
        ranked = p[order] * len(p) / np.arange(1, len(order) + 1)
        q[order] = np.minimum(1, np.minimum.accumulate(ranked[::-1])[::-1])
    return q


def summarize(data):
    permutations = np.array(list(itertools.permutations(range(8))), dtype=np.int16)
    points, groups, tests = [], [], []
    for feature in FEATURES:
        means = []
        for temperature, block in data.groupby("temperature", sort=True):
            values = block[feature]
            mean, sd = float(values.mean()), float(values.std(ddof=1))
            means.append(mean)
            groups.append(dict(feature=feature, temperature=int(temperature),
                               n_valid=int(values.count()), n_missing=int(values.isna().sum()),
                               mean=mean, sd=sd, unit="%"))
            for offset, (_, row) in zip(np.linspace(-0.55, 0.55, len(block)), block.iterrows()):
                points.append(dict(feature=feature, sample=row["sample"],
                                   temperature=int(temperature), replicate=int(row["replicate"]),
                                   methylation=row[feature], missing=bool(pd.isna(row[feature])),
                                   x_plot=float(temperature + offset), unit="%"))
        rho, p, extreme, status = exact_spearman(means, permutations)
        tests.append(dict(feature=feature, n_groups=8, n_valid_samples=int(data[feature].count()),
                          rho=rho, p=p, extreme_permutations=extreme, total_permutations=40320,
                          status=status, test="two-sided exhaustive pairing permutation Spearman"))
    statistics = pd.DataFrame(tests)
    statistics["q_bh"] = bh_adjust(statistics["p"])
    return pd.DataFrame(points), pd.DataFrame(groups), statistics


def load_auditor(directory):
    if not Path(directory / 'audit_panel_alignment.py').is_file():
        _skip_optional_qa(f"Optional QA tool not found: {directory / 'audit_panel_alignment.py'}")
        return None
    spec = importlib.util.spec_from_file_location("temperature_panel_alignment", directory / "audit_panel_alignment.py")
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


def check_alignment(fig, labels, auditor, output):
    if auditor is None:
        return _skip_optional_qa("audit_panel_alignment.py is unavailable", output)
    else:
        auditor.require_matplotlib_panel_alignment(fig, tolerance_pt=1.5, gutter_tolerance_pt=1.5, strict=True)
        # Explicit anchors preserve the requested parenthesized labels in the audit.
        manifest = auditor.matplotlib_layout_manifest(fig)
        for panel, label in zip(manifest["panels"], labels):
            anchor = label.get_transform().transform(label.get_position())
            anchor = fig.dpi_scale_trans.inverted().transform(anchor) * 72
            panel["panel_label"] = label.get_text().strip("()")
            panel["panel_label_anchor_pt"] = list(map(float, anchor))
        report = auditor.audit_layout_manifest(manifest, tolerance_pt=1.5,
                                               gutter_tolerance_pt=1.5, require_panel_labels=True)
        auditor.write_json_report(report, output)
        if auditor.exit_code(report, strict=True):
            raise RuntimeError(auditor.render_text(report, strict=True))


def caption_lines(groups):
    count_descriptions = []
    for feature in FEATURES:
        sg = groups[groups.feature == feature]
        counts = ", ".join(f"{int(r.temperature)} °C: {int(r.n_valid)}" for r in sg.itertuples())
        count_descriptions.append(counts)
    if len(set(count_descriptions)) == 1:
        counts_text = "Valid replicate n at each site (" + count_descriptions[0] + ")."
    else:
        counts_text = "\n".join(f"{f}: replicate n ({c})." for f, c in zip(FEATURES, count_descriptions))
    return [
        "Faint points: individual samples; solid points and error bars: mean ± sample SD. Lines describe trends only.",
        counts_text,
        "Spearman n = 8 temperature means; two-sided exact permutation P (unadjusted). BH q values: statistics table.",
        "Shading marks the 24 °C control. Missing measurements are omitted without imputation.",
    ]


def make_figure(points, groups, statistics, auditor, qa_dir):
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
                         "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 8,
                         "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.6, "pdf.fonttype": 42, "svg.fonttype": "none",
                         "savefig.facecolor": "white", "legend.frameon": False})
    lines = [line for entry in caption_lines(groups) for line in entry.splitlines()]
    extra = max(0, len(lines) - 4) * 4
    height_mm = 132 + extra
    fig = plt.figure(figsize=(183 / 25.4, height_mm / 25.4))
    grid = fig.add_gridspec(2, 2, left=0.105, right=0.975,
                           bottom=(32 + extra) / height_mm, top=1 - 8 / height_mm,
                           hspace=0.55, wspace=0.27)
    labels, panel_records = [], []
    for index, feature in enumerate(FEATURES):
        ax = fig.add_subplot(grid[index // 2, index % 2])
        ax.set_gid(chr(97 + index))
        sp = points[(points.feature == feature) & ~points.missing]
        sg = groups[groups.feature == feature]
        st = statistics[statistics.feature == feature].iloc[0]
        ax.axvspan(22, 26, color=CONTROL, alpha=0.75, lw=0, zorder=0)
        ax.set_axisbelow(True)
        ax.grid(axis="y", color="#E3E3E3", lw=0.4)
        ax.scatter(sp.x_plot, sp.methylation, s=13, color=BLUE, alpha=0.40,
                   edgecolors="none", zorder=3)
        ax.plot(sg.temperature, sg["mean"], color=BLUE, lw=0.8, zorder=2)
        ax.errorbar(sg.temperature, sg["mean"], yerr=sg.sd, fmt="o", ms=3.5,
                    color=BLUE, ecolor=BLUE, elinewidth=0.75, capsize=2.5, capthick=0.75,
                    markeredgecolor="white", markeredgewidth=0.4, zorder=4)
        ymax = 80 if index < 2 else 100
        ax.set(xlim=(2, 34), ylim=(0, ymax), xticks=TEMPERATURES,
               yticks=np.arange(0, ymax + 1, 20), xlabel="Temperature (°C)", ylabel="CHH methylation (%)")
        ax.tick_params(length=2.5, width=0.6)
        ax.set_title(feature, pad=8, fontweight="bold")
        labels.append(ax.text(0, 1, f"({chr(97 + index)})", fontsize=8, fontweight="bold",
                              ha="left", va="bottom", transform=ax.transAxes +
                              ScaledTranslation(-23 / 72, 7 / 72, fig.dpi_scale_trans)))
        if st.status == "ok":
            p_text = f"{st.p:.3f}" if st.p >= 0.001 else f"{st.p:.2e}"
            annotation = f"ρ = {st.rho:.2f}; P = {p_text}"
        else:
            annotation = "ρ = NA; P = NA (constant)"
        ax.text(0.04, 0.94, annotation, transform=ax.transAxes, va="top", fontsize=7)
        bounds = [float((sg["mean"] - sg.sd).min()), float((sg["mean"] + sg.sd).max())]
        if min(bounds[0], sp.methylation.min()) < 0 or max(bounds[1], sp.methylation.max()) > ymax:
            raise ValueError(f"{feature}: data or mean ± SD exceeds the specified 0–{ymax}% axis; export stopped to avoid clipping")
        panel_records.append(dict(panel=f"({chr(97+index)})", feature=feature,
                                  valid_points=len(sp), means=len(sg), error_bars=len(sg),
                                  uncertainty="sample SD (ddof=1)", uncertainty_bounds=bounds,
                                  ylim=[0, ymax], missing=int(points[points.feature == feature].missing.sum())))
    for i, line in enumerate(lines):
        fig.text(0.055, (4 + (len(lines) - 1 - i) * 3.7) / height_mm, line, fontsize=6, va="bottom")
    fig.canvas.draw()
    check_alignment(fig, labels, auditor, qa_dir / "alignment.json")
    return fig, panel_records


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
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--qa-tools", type=Path, default=DEFAULT_QA_TOOLS)
    parser.add_argument("--validate-only", action="store_true", help="validate data without creating outputs")
    args = parser.parse_args(argv)
    data = load_data(args.input)
    if args.validate_only:
        print(f"PASS: {len(data)} sample rows; valid counts: {data[list(FEATURES)].count().to_dict()}")
        return 0
    points, groups, statistics = summarize(data)
    auditor = load_auditor(args.qa_tools)
    # PDF QA dependencies are checked only when the optional tools run.
    args.output_dir.mkdir(parents=True, exist_ok=True)
    qa_dir = args.output_dir / f"{PREFIX}_QA"
    qa_dir.mkdir(exist_ok=True)
    run_audit([sys.executable, str(args.qa_tools / "validate_figure.py"), __file__], qa_dir / "source-preflight.txt")
    fig, panel_records = make_figure(points, groups, statistics, auditor, qa_dir)
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
    for name, table in [("points", points), ("groups", groups), ("statistics", statistics)]:
        table.to_csv(stem.with_suffix(f".{name}.tsv"), sep="\t", index=False, na_rep="NA", float_format="%.12g")
    caption = "# Temperature-associated CHH methylation at four sites in ***RCA***\n\n"
    caption += "Panels (a) CHH_1634, (b) CHH_1635, (c) CHH_1694 and (d) CHH_1696. "
    caption += "\n\n".join(caption_lines(groups))
    caption += ("\n\nPanels (a–b) use a 0–80% y-axis; panels (c–d) use a 0–100% y-axis. "
                "SD uses ddof=1. Sample rows are treated as biological replicates. "
                "Each point retains its measured percentage; horizontal offsets only separate samples visually. "
                "All 8! = 40,320 permutations are enumerated, including tied-rank multiplicities. "
                "BH correction covers the four planned tests; undefined correlations remain NA. "
                "These descriptive associations do not establish causality.\n")
    stem.with_suffix(".caption.md").write_text(caption, encoding="utf-8")
    pd.DataFrame(panel_records).to_csv(qa_dir / "panels.tsv", sep="\t", index=False)
    report = dict(status=("automated_QA_INCOMPLETE" if _OPTIONAL_QA_SKIPPED else "automated_QA_PASS"), panels=panel_records,
                  input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
                  selected_columns=list(FEATURES), source_sample_rows=len(data),
                  missing_measurements=points[points.missing][["feature", "sample", "temperature"]].to_dict("records"),
                  excluded_columns="Unselected sites and expression columns are outside the requested four-site analysis.",
                  source_points_rows=len(points), plotted_points=int((~points.missing).sum()),
                  groups=len(groups), correlation_unit="eight temperature means",
                  versions=dict(python=sys.version.split()[0], numpy=np.__version__, pandas=pd.__version__,
                                scipy=scipy.__version__, matplotlib=matplotlib.__version__))
    (qa_dir / "summary.json").write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(statistics[["feature", "rho", "p", "q_bh"]].to_string(index=False))
    print(f"PASS: {report['plotted_points']} points, {len(groups)} means/SD; PDF/SVG/600 dpi PNG/TIFF exported to {args.output_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, RuntimeError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
