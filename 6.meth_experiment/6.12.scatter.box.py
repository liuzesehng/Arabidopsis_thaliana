#!/usr/bin/env python3
"""CHH methylation versus RCA abundance/composition: 17 descriptive 3+2 grids.

Expression replicates are plotted in eight equally spaced temperature groups,
ordered by mean methylation, with that mean as the numeric tick label.
Methylation and expression use their own available observations; blank
methylation is never replaced by zero. Correlations use eight group means,
not the vertically aligned expression replicates. Lines are descriptive only.
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

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rca-meth-mpl"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.transforms import ScaledTranslation
import numpy as np
import pandas as pd
from scipy.stats import rankdata

PROJECT = Path(__file__).resolve().parents[2]
DEFAULT_FEATURES = PROJECT / "list/RCA/TF_CHH/AT2G39730_CHH_combined_historical.panel-c-features.tsv"
DEFAULT_INPUT = PROJECT / "list/meth_experiment/Alt.all.snp_meth.tsv"
DEFAULT_OUTPUT = PROJECT / "list/meth_experiment"
DEFAULT_QA_TOOLS = Path("/datapool/home/2023102768/.codex/skills/nature-figure/scripts")
TEMPERATURES = tuple(range(4, 33, 4))
METRICS = ("total", "α", "β", "α/%", "β/%")
COLORS = ("#4D4D4D", "#2F6DAE", "#D9822B", "#2F6DAE", "#D9822B")
TITLES = (r"Total $\mathbf{RCA}$", r"$\mathbf{RCA}$ α", r"$\mathbf{RCA}$ β",
          r"$\mathbf{RCA}$ α proportion", r"$\mathbf{RCA}$ β proportion")
PREFIX = "RCA_methylation_expression"
PERMUTATIONS = np.array(list(itertools.permutations(range(8))), dtype=np.int16)


def read_tsv(path: Path) -> pd.DataFrame:
    """Retain blanks as missing, reject duplicate headers and ragged records."""
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows or not rows[0] or any(not name for name in rows[0]):
        raise ValueError(f"{path.name}: missing or blank header")
    if len(set(rows[0])) != len(rows[0]):
        raise ValueError(f"{path.name}: duplicate column names")
    if any(len(row) != len(rows[0]) for row in rows[1:]):
        raise ValueError(f"{path.name}: inconsistent number of fields")
    return pd.DataFrame(rows[1:], columns=rows[0])


def load_inputs(features_path: Path, input_path: Path):
    features = read_tsv(features_path).iloc[:, 0].tolist()
    if len(features) != 17 or len(set(features)) != 17:
        raise ValueError("feature list must contain exactly 17 unique IDs")
    if any(not re.fullmatch(r"CHH_\d+", name) for name in features):
        raise ValueError("feature IDs must have the form CHH_<integer>")
    source = read_tsv(input_path)
    required = ["Temperature", *features, *METRICS]
    absent = [name for name in required if name not in source]
    if absent:
        raise ValueError(f"missing required columns: {', '.join(absent)}")
    data = source[required].copy().rename(columns={"Temperature": "sample"})
    parsed = data["sample"].str.extract(r"^T(\d+)-(\d+)$")
    if parsed.isna().any().any() or data["sample"].duplicated().any():
        raise ValueError("sample IDs must be unique and match T<temperature>-<replicate>")
    data["temperature"] = parsed[0].astype(int)
    data["replicate"] = parsed[1].astype(int)
    if (data["replicate"] < 1).any() or data.duplicated(["temperature", "replicate"]).any():
        raise ValueError("replicate numbers must be positive and unique per temperature")
    counts = data.groupby("temperature").size().to_dict()
    if counts != {t: 3 for t in TEMPERATURES}:
        raise ValueError(f"expected eight temperatures 4–32 °C with 3 rows each; got {counts}")
    for name in [*features, *METRICS]:
        values = data[name].str.strip()
        # Only an empty field is missing. Text such as 'NA' or 'inf' is invalid.
        numeric = pd.to_numeric(values.replace("", np.nan), errors="raise").astype(float)
        if np.isinf(numeric).any() or (name in METRICS and numeric.isna().any()):
            raise ValueError(f"{name}: non-finite or missing expression value")
        if ((numeric < 0) | ((numeric > 100) if name not in METRICS[:3] else False)).any():
            raise ValueError(f"{name}: values outside permitted range")
        data[name] = numeric
    valid_counts = data.groupby("temperature")[features].count()
    if (valid_counts == 0).any().any():
        bad = [(str(c), int(t)) for c in features for t in TEMPERATURES if valid_counts.loc[t, c] == 0]
        raise ValueError(f"whole temperature group has missing methylation: {bad}")
    if not np.allclose(data["total"], data["α"] + data["β"], rtol=1e-6, atol=0.02):
        raise ValueError("total differs from α + β beyond atol=0.02, rtol=1e-6")
    if (abs(data["α/%"] + data["β/%"] - 100) > 0.02 + np.spacing(100.)).any():
        raise ValueError("α/% + β/% differs from 100 beyond 0.02 plus one ULP")
    for name in METRICS[:3]:
        data[f"{name}_ln"] = np.log1p(data[name])
    return features, data.sort_values(["temperature", "replicate"], kind="stable").reset_index(drop=True)


def exact_spearman(x, y):
    """Two-sided exhaustive pairing test; tied ranks retain multiplicities."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.shape != (8,) or y.shape != (8,) or not np.isfinite([x, y]).all():
        raise ValueError("Spearman requires eight finite pairs")
    rx, ry = rankdata(x), rankdata(y)
    rx -= rx.mean()
    ry -= ry.mean()
    denom = np.linalg.norm(rx) * np.linalg.norm(ry)
    if denom == 0:
        return np.nan, np.nan, 0, "constant_variable"
    observed = float(rx @ ry)
    permuted = np.einsum("ij,j->i", rx[PERMUTATIONS], ry)
    extreme = int(np.count_nonzero(np.abs(permuted) >= abs(observed) - 1e-12))
    return float(np.clip(observed / denom, -1, 1)), extreme / len(PERMUTATIONS), extreme, "ok"


def bh_adjust(pvalues):
    p = np.asarray(pvalues, dtype=float)
    q = np.full_like(p, np.nan)
    finite = np.flatnonzero(np.isfinite(p))
    order = finite[np.argsort(p[finite], kind="stable")]
    if len(order):
        adjusted = p[order] * len(order) / np.arange(1, len(order) + 1)
        q[order] = np.minimum(1, np.minimum.accumulate(adjusted[::-1])[::-1])
    return q


def summarize(features, data):
    groups, points, tests = [], [], []
    for feature in features:
        for metric in METRICS:
            value_col = f"{metric}_ln" if metric in METRICS[:3] else metric
            xmeans, ymeans = [], []
            for temperature, block in data.groupby("temperature", sort=True):
                x = block[feature]
                y = block[value_col]
                xm, ym = float(x.mean()), float(y.mean())
                xmeans.append(xm)
                ymeans.append(ym)
                groups.append(dict(feature=feature, metric=metric, temperature=int(temperature),
                                   methylation_mean=xm, methylation_sd=float(x.std(ddof=1)),
                                   n_methylation=int(x.count()), expression_mean=ym,
                                   expression_sd=float(y.std(ddof=1)), n_expression=int(y.count()),
                                   expression_unit="ln(TPM+1)" if metric in METRICS[:3] else "%"))
                for idx, row in block.iterrows():
                    points.append(dict(feature=feature, metric=metric, sample=row["sample"],
                                       temperature=int(temperature), methylation_raw=row[feature],
                                       methylation_missing=bool(pd.isna(row[feature])), x_plot=xm,
                                       expression_raw=row[metric], y_plot=row[value_col]))
            rho, p, extreme, status = exact_spearman(xmeans, ymeans)
            tests.append(dict(feature=feature, metric=metric, n_groups=8, rho=rho, p=p,
                              extreme_permutations=extreme, total_permutations=40320,
                              status=status, test="two-sided exhaustive pairing permutation Spearman"))
    stats = pd.DataFrame(tests)
    stats["q_bh"] = bh_adjust(stats["p"])
    group_frame, point_frame = pd.DataFrame(groups), pd.DataFrame(points)
    point_frame = point_frame.rename(columns={"x_plot": "methylation_mean"})
    for feature in features:
        order = group_frame[(group_frame.feature == feature) & (group_frame.metric == "total")]
        order = order.sort_values(["methylation_mean", "temperature"], kind="stable")
        positions = {int(t): i for i, t in enumerate(order.temperature)}
        for frame in (group_frame, point_frame):
            mask = frame.feature == feature
            frame.loc[mask, "x_plot"] = frame.loc[mask, "temperature"].map(positions)
    group_frame["x_plot"] = group_frame["x_plot"].astype(int)
    point_frame["x_plot"] = point_frame["x_plot"].astype(int)
    return group_frame, point_frame, stats


def configure_style():
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                         "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 8,
                         "xtick.labelsize": 6, "ytick.labelsize": 6,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.6, "pdf.fonttype": 42, "svg.fonttype": "none",
                         "mathtext.fontset": "custom", "mathtext.bf": "DejaVu Sans:bold:italic",
                         "mathtext.rm": "DejaVu Sans", "mathtext.it": "DejaVu Sans:italic",
                         "savefig.facecolor": "white", "legend.frameon": False})


def methylation_ticks(ax, summary):
    """One independent tick per temperature; identical means repeat as labels."""
    order = summary.sort_values("x_plot")
    ax.set_xticks(order.x_plot, labels=[f"{value:.2f}" for value in order.methylation_mean])
    ax.tick_params(axis="x", length=2.5, width=0.5, pad=3)
    for label in ax.get_xticklabels():
        label.set_fontstyle("normal")
        label.set_rotation(60)
        label.set_rotation_mode("anchor")
        label.set_ha("right")
        label.set_va("top")
    return ax.get_xticklabels()


def p_label(p):
    return "NA" if not np.isfinite(p) else (f"{p:.2e}" if p < 0.001 else f"{p:.3f}")


def load_alignment(tools_dir):
    path = tools_dir / "audit_panel_alignment.py"
    spec = importlib.util.spec_from_file_location("rca_panel_alignment", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load alignment auditor: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def alignment_report(fig, labels, auditor, output):
    # The shared helper recognizes bare letters only. Measure parenthesized
    # labels explicitly and run its same geometric audit with labels required.
    first = auditor.require_matplotlib_panel_alignment(fig, tolerance_pt=1.5,
                                                       gutter_tolerance_pt=1.5, strict=True)
    manifest = auditor.matplotlib_layout_manifest(fig)
    for panel, text in zip(manifest["panels"], labels):
        xy = text.get_transform().transform(text.get_position())
        xy = fig.dpi_scale_trans.inverted().transform(xy) * 72
        panel["panel_label"] = text.get_text().strip("()")
        panel["panel_label_anchor_pt"] = list(map(float, xy))
    report = auditor.audit_layout_manifest(manifest, tolerance_pt=1.5,
                                            gutter_tolerance_pt=1.5, require_panel_labels=True)
    auditor.write_json_report(report, output)
    if auditor.exit_code(report, strict=True):
        raise RuntimeError(auditor.render_text(report, strict=True))
    return report


def make_figure(feature, groups, points, stats, auditor, stem):
    width_mm = 183
    figure_width_in = width_mm / 25.4
    figure_height_in = 150 / 25.4
    fig = plt.figure(figsize=(figure_width_in, figure_height_in))
    # Equal two-column spans; lower two panels centered under upper three.
    grid = fig.add_gridspec(2, 6, left=0.10, right=0.975, bottom=0.24, top=0.87,
                           hspace=0.95, wspace=1.10)
    spans = [(0, 0), (0, 2), (0, 4), (1, 1), (1, 3)]
    labels, panel_audit = [], []
    for i, (metric, (row, col)) in enumerate(zip(METRICS, spans)):
        ax = fig.add_subplot(grid[row, col:col + 2])
        ax.set_gid(chr(97 + i))
        sg = groups[(groups.feature == feature) & (groups.metric == metric)]
        sp = points[(points.feature == feature) & (points.metric == metric)]
        st = stats[(stats.feature == feature) & (stats.metric == metric)].iloc[0]
        order = sg.sort_values(["methylation_mean", "temperature"], kind="stable")
        color = COLORS[i]
        # Every one of the 24 expression observations is kept, without x jitter.
        ax.scatter(sp.x_plot, sp.y_plot, s=9, color=color, alpha=0.40,
                   edgecolors="none", zorder=3)
        ax.plot(order.x_plot, order.expression_mean, color=color, lw=0.7,
                alpha=0.80, zorder=2)
        ax.errorbar(sg.x_plot, sg.expression_mean, yerr=sg.expression_sd,
                    fmt="o", ms=3.0, color=color, ecolor=color, elinewidth=0.65,
                    capsize=2, capthick=0.65, markeredgecolor="white", markeredgewidth=0.35,
                    zorder=4)
        all_y = points.loc[points.metric == metric, "y_plot"].to_numpy()
        all_g = groups.loc[groups.metric == metric]
        ymin = min(all_y.min(), (all_g.expression_mean - all_g.expression_sd).min())
        ymax = max(all_y.max(), (all_g.expression_mean + all_g.expression_sd).max())
        span = max(float(ymax - ymin), 0.1)
        ax.set_ylim(ymin - span * 0.10, ymax + span * 0.12)
        ax.set_xlim(-0.5, 7.5)
        ax.set_ylabel("ln(TPM+1)" if i < 3 else "Proportion (%)", labelpad=3)
        ax.tick_params(axis="y", length=2.5, width=0.5, pad=2)
        ax.set_title(TITLES[i], loc="center", pad=24, fontweight="bold")
        label = ax.text(0, 1, f"({chr(97 + i)})", ha="left", va="bottom", fontsize=8,
                        fontweight="bold", transform=ax.transAxes +
                        ScaledTranslation(-16 / 72, 23 / 72, fig.dpi_scale_trans))
        labels.append(label)
        stat_text = (f"ρ = {st.rho:.2f}; P = {p_label(st.p)}"
                     if st.status == "ok" else "ρ = NA; P = NA")
        ax.text(0.5, 1, stat_text, ha="center", va="bottom", fontsize=6,
                transform=ax.transAxes + ScaledTranslation(0, 7 / 72, fig.dpi_scale_trans))
        tick_text = methylation_ticks(ax, sg)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        lowest = min(t.get_window_extent(renderer).y0 for t in tick_text)
        fig_y = (lowest - 10 * fig.dpi / 72) / fig.bbox.height
        fig.text(ax.get_position().x0 + ax.get_position().width / 2, fig_y,
                 "Mean CHH methylation (%)", ha="center", va="top", fontsize=7)
        panel_audit.append(dict(feature=feature, panel=f"({chr(97+i)})", metric=metric,
                                expression_points=len(sp), temperature_means=len(sg),
                                distinct_x_ticks=len(ax.get_xticks()),
                                ymin=ax.get_ylim()[0], ymax=ax.get_ylim()[1],
                                sd_definition="sample SD (ddof=1), expression axis only"))
    fig.suptitle(feature, y=0.98, fontsize=9)
    fig.text(0.5, 0.083, "Faint points: expression replicates; solid points: mean ± SD. Lines describe trends only.",
             ha="center", fontsize=6)
    n_meth = groups[(groups.feature == feature) & (groups.metric == "total")]
    counts = "; ".join(f"{int(r.temperature)}°C: {int(r.n_methylation)}" for r in n_meth.itertuples())
    fig.text(0.5, 0.059, f"Methylation n ({counts}); expression n = 3 per temperature.",
             ha="center", fontsize=6)
    fig.text(0.5, 0.035, "X: equally spaced temperature groups, ordered by mean methylation (ties: increasing temperature).",
             ha="center", fontsize=6)
    fig.text(0.5, 0.011, "ρ: Spearman correlation; n = 8 temperature means; exact two-sided permutation P. BH q: statistics table.",
             ha="center", fontsize=6)
    alignment_report(fig, labels, auditor, stem.with_suffix(".alignment.json"))
    return fig, panel_audit


def run_audit(command, report_path):
    result = subprocess.run(command, text=True, capture_output=True)
    report_path.write_text(result.stdout + result.stderr, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(f"audit failed ({report_path.name}): {result.stdout[-2000:]} {result.stderr[-500:]}")


def write_json(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=Path, default=DEFAULT_FEATURES)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--qa-tools", type=Path, default=DEFAULT_QA_TOOLS)
    parser.add_argument("--validate-only", action="store_true", help="validate inputs, write nothing")
    args = parser.parse_args(argv)
    features, data = load_inputs(args.features, args.input)
    if args.validate_only:
        print(f"PASS: {len(features)} features, {len(data)} samples, eight temperature groups")
        return 0
    groups, points, stats = summarize(features, data)
    auditor = load_alignment(args.qa_tools)
    if importlib.util.find_spec("fitz") is None:
        raise RuntimeError("PDF collision QA requires PyMuPDF (fitz); rendering not started")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    qa_dir = args.output_dir / f"{PREFIX}_QA"
    qa_dir.mkdir(exist_ok=True)
    run_audit([sys.executable, str(args.qa_tools / "validate_figure.py"), __file__],
              qa_dir / "source-preflight.txt")
    configure_style()
    for suffix, frame in [("groups", groups), ("points", points), ("statistics", stats)]:
        frame.to_csv(args.output_dir / f"{PREFIX}.{suffix}.tsv", sep="\t", index=False,
                     na_rep="NA", float_format="%.12g")
    all_panels = []
    for index, feature in enumerate(features, 1):
        name = f"RCA_{feature}_methylation_expression"
        stem = args.output_dir / name
        fig, panels = make_figure(feature, groups, points, stats, auditor, qa_dir / name)
        try:
            fig.savefig(stem.with_suffix(".pdf"))
            run_audit([sys.executable, str(args.qa_tools / "audit_pdf_text.py"),
                       str(stem.with_suffix(".pdf")), "--min-pt", "5", "--json"],
                      qa_dir / f"{name}.pdf-text.json")
            run_audit([sys.executable, str(args.qa_tools / "audit_figure_collisions.py"),
                       str(stem.with_suffix(".pdf")), "--strict", "--json-out",
                       str(qa_dir / f"{name}.collision-audit.json")],
                      qa_dir / f"{name}.collision-audit.txt")
            fig.savefig(stem.with_suffix(".svg"))
            fig.savefig(stem.with_suffix(".png"), dpi=600)
            fig.savefig(stem.with_suffix(".tiff"), dpi=600, pil_kwargs={"compression": "tiff_lzw"})
        finally:
            plt.close(fig)
        all_panels.extend(panels)
        print(f"[{index:02d}/17] {feature}: exported PDF/SVG/PNG/TIFF; rendered QA PASS", flush=True)
    pd.DataFrame(all_panels).to_csv(qa_dir / "panels.tsv", sep="\t", index=False)
    missing = [{"feature": f, "sample": str(row["sample"]), "temperature": int(row["temperature"])}
               for f in features for _, row in data.loc[data[f].isna()].iterrows()]
    report = dict(status="automated_QA_PASS", figure_count=17, panels=85, tests=85,
                  input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
                  features_sha256=hashlib.sha256(args.features.read_bytes()).hexdigest(),
                  missing_methylation=missing, excluded_expression_observations=0,
                  expression_transform="ln1p per replicate before mean and sample SD",
                  missing_policy="separate available observations; no imputation",
                  trend_order="methylation mean ascending, then temperature ascending",
                  tick_policy="eight equally spaced temperature groups; labels are methylation means only; ties remain separate",
                  plot_size_mm=[183, 150],
                  correlation="eight group means, exact 40320 pairing permutations, two-sided abs(rho)",
                  multiplicity="BH across all finite tests; constant-variable tests are NA",
                  alpha_beta_sum_tolerance="atol=0.02, rtol=1e-6; percentage sum=0.02 plus one ULP",
                  versions={"python": sys.version.split()[0], "numpy": np.__version__,
                            "pandas": pd.__version__, "matplotlib": matplotlib.__version__})
    write_json(qa_dir / "summary.json", report)
    (qa_dir / "README.md").write_text(
        "# RCA methylation–expression figures\n\n"
        "Question: how do temperature-group CHH methylation means associate with RCA abundance and composition?\n"
        "Panels a–c show total abundance and its α/β components; d–e show relative composition.\n"
        "These are descriptive associations across eight temperature groups, not evidence of causation.\n\n"
        "All 24 expression observations are displayed per panel in eight equally spaced temperature groups.\n"
        "Missing methylation is excluded only from methylation summaries and recorded in summary.json.\n"
        "Vertical bars are sample SD; methylation SD is in the group table but is not drawn.\n"
        "Lines follow ascending methylation, with temperature breaking exact ties.\n"
        "Tick labels show methylation means only; identical means retain distinct group positions.\n"
        "x_plot in the point and group tables is the zero-based categorical position; methylation_mean retains the measured value.\n"
        "Spearman uses measured methylation means, not categorical positions.\n"
        "Figure titles contain the feature ID only; full panel titles are bold, RCA is bold italic, and x tick labels use upright numerals rotated 60 degrees; all n values are in the caption.\n"
        "Figure dimensions are 183 × 150 mm; raster exports use 600 dpi.\n"
        "Panel-specific y limits are identical across all 17 figures.\n\n"
        "P values are exact exhaustive two-sided pairing permutation probabilities, including observed pairing.\n"
        "Ties use average ranks; BH q values correct the family of 85 tests.\n"
        "α and β proportions are complementary; their tests are not independent evidence.\n"
        "All PDF alignment, text-size and strict collision checks passed; see per-figure reports.\n",
        encoding="utf-8")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, OSError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
