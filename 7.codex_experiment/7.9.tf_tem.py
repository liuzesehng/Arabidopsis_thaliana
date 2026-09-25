#!/usr/bin/env python3
"""SPL8 temperature scatter; outputs default to list/meth_experiment.

Show all 24 samples and mean ± sample SD after ln(TPM+1) transformation.
Spearman uses eight temperature means and all 8! two-sided permutations.
Connecting lines describe means, not a regression fit.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import os
from pathlib import Path
import re
import sys
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "spl8-temperature-mpl"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "spl8-temperature-cache"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import rankdata

PROJECT = Path(__file__).resolve().parents[2]
GENE = "AT1G02065"
TEMPERATURES = tuple(range(4, 33, 4))
PREFIX = "SPL8_temperature_spearman"


def load_expression(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        if not header or header[0] != "Name":
            raise ValueError("first column must be Name")
        if len(set(header)) != len(header) or any(not s for s in header):
            raise ValueError("blank or duplicate column names")
        samples = header[1:]
        if any(re.fullmatch(r"\d+_[1-9]\d*", s) is None for s in samples):
            raise ValueError("sample IDs must be <temperature>_<positive replicate>")
        target = []
        for line, row in enumerate(reader, 2):
            if len(row) != len(header):
                raise ValueError(f"line {line}: inconsistent field count")
            if row[0] == GENE:
                target.append(row[1:])
    if len(target) != 1:
        raise ValueError(f"expected exactly one {GENE} row; found {len(target)}")
    pairs = [tuple(map(int, s.split("_"))) for s in samples]
    if len(set(pairs)) != len(pairs):
        raise ValueError("duplicate normalized sample IDs")
    values = np.asarray(target[0], dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("SPL8 expression must be finite, nonnegative TPM")
    data = pd.DataFrame(dict(gene=GENE, symbol="SPL8", sample=samples,
                             temperature=[p[0] for p in pairs],
                             replicate=[p[1] for p in pairs], tpm=values))
    if data.groupby("temperature").size().to_dict() != {t: 3 for t in TEMPERATURES}:
        raise ValueError("expected eight temperatures 4–32 °C, three samples each")
    data["ln_tpm_plus1"] = np.log1p(values)
    return data.sort_values(["temperature", "replicate"]).reset_index(drop=True)


def exact_spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if x.shape != (8,) or y.shape != (8,) or not np.isfinite([x, y]).all():
        raise ValueError("Spearman requires eight finite pairs")
    rx, ry = rankdata(x), rankdata(y)
    rx, ry = rx - rx.mean(), ry - ry.mean()
    denominator = np.linalg.norm(rx) * np.linalg.norm(ry)
    if denominator == 0:
        return dict(rho=None, p=None, extreme_permutations=0,
                    total_permutations=40320, status="constant_variable")
    permutations = np.array(list(itertools.permutations(range(8))), dtype=np.int16)
    observed = float(rx @ ry)
    permuted = np.einsum("ij,j->i", rx[permutations], ry)
    extreme = int(np.count_nonzero(np.abs(permuted) >= abs(observed) - 1e-12))
    return dict(rho=float(np.clip(observed / denominator, -1, 1)),
                p=extreme / len(permutations), extreme_permutations=extreme,
                total_permutations=len(permutations), status="ok")


def summarize(data):
    groups = data.groupby("temperature", as_index=False).agg(
        n=("ln_tpm_plus1", "size"), mean=("ln_tpm_plus1", "mean"),
        sd=("ln_tpm_plus1", "std"), mean_tpm=("tpm", "mean"))
    result = dict(gene=GENE, symbol="SPL8", n_samples=len(data), n_groups=len(groups),
                  expression_unit="ln(TPM+1)", aggregation="mean after log1p",
                  test="Spearman; exhaustive two-sided pairing permutation",
                  multiple_testing="none (one prespecified test)")
    result.update(exact_spearman(groups.temperature, groups["mean"]))
    return groups, result


def make_figure(data, groups, result):
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                         "font.size": 6.6, "axes.labelsize": 7,
                         "xtick.labelsize": 6.3, "ytick.labelsize": 6.3,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.75, "pdf.fonttype": 42,
                         "svg.fonttype": "none", "legend.frameon": False})
    width_mm, height_mm = 89, 75
    fig, ax = plt.subplots(figsize=(width_mm / 25.4, height_mm / 25.4))
    fig.subplots_adjust(left=0.17, right=0.97, bottom=0.30, top=0.95)
    ax.axvspan(22, 26, color="#E9E9E9", alpha=0.9, linewidth=0, zorder=0)
    # Horizontal offsets separate replicates visually; statistics use true temperatures.
    for temperature in TEMPERATURES:
        rows = data.loc[data.temperature.eq(temperature)].sort_values("replicate")
        ax.scatter(temperature + np.array([-0.42, 0.0, 0.42]), rows.ln_tpm_plus1,
                   s=18, color="#85AEC5", alpha=0.78, edgecolors="white",
                   linewidths=0.45, zorder=3)
    ax.plot(groups.temperature, groups["mean"], lw=1.15, color="#225E7A", alpha=0.95, zorder=4)
    ax.errorbar(groups.temperature, groups["mean"], yerr=groups.sd,
                fmt="D", ms=4, color="#225E7A", elinewidth=0.9,
                capsize=2.5, capthick=0.9, markeredgecolor="white",
                markeredgewidth=0.55, zorder=5)
    high = max(data.ln_tpm_plus1.max(), (groups["mean"] + groups.sd).max())
    ax.set_ylim(0, max(float(high), 0.1) * 1.28)
    ax.set_xlim(2, 34)
    ax.set_xticks(TEMPERATURES)
    ax.set_xlabel("Temperature (°C)", labelpad=5)
    ax.set_ylabel("ln(TPM+1)", labelpad=5)
    ax.tick_params(direction="out", length=2.7, width=0.65)
    annotation = (f"ρ = {result['rho']:.3f}; P = {result['p']:.3f}"
                  if result["status"] == "ok" else "ρ and P undefined\n(constant means)")
    ax.text(0.025, 0.965, annotation, transform=ax.transAxes,
            ha="left", va="top", fontsize=5.8, color="#222222", zorder=6)
    fig.text(0.5, 0.085, "Points: 3 replicates/temperature; diamonds: mean ± SD.\n"
             "Line: descriptive trend; grey band: 24 °C reference.\n"
             "Spearman: 8 means; exact two-sided permutation P.",
             ha="center", va="center", fontsize=5.8, color="#222222", linespacing=1.5)
    fig.canvas.draw()
    return fig


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=PROJECT / "list/coex_experiment/RCA.tran.txt")
    parser.add_argument("--output-dir", type=Path, default=PROJECT / "list/meth_experiment")
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args(argv)
    data = load_expression(args.input)
    groups, result = summarize(data)
    if args.validate_only:
        print(json.dumps(result, ensure_ascii=False, allow_nan=False))
        return 0
    args.output_dir.mkdir(parents=True, exist_ok=True)
    stem = args.output_dir / PREFIX
    data.to_csv(stem.with_suffix(".points.tsv"), sep="\t", index=False)
    groups.to_csv(stem.with_suffix(".summary.tsv"), sep="\t", index=False)
    pd.DataFrame([result]).to_csv(stem.with_suffix(".statistics.tsv"), sep="\t", index=False)
    fig = make_figure(data, groups, result)
    for suffix in (".pdf", ".svg", ".png", ".tiff"):
        fig.savefig(stem.with_suffix(suffix), dpi=600, facecolor="white")
    plt.close(fig)
    notes = dict(question="Does SPL8 expression change monotonically with temperature?",
                 archetype="quantitative single panel", backend="Python/matplotlib",
                 size_mm=[89, 75], alignment="NOT APPLICABLE: one panel",
                 display="replicate x offsets -0.42/0/0.42 °C; diamonds show means; grey band is 24 °C visual reference only",
                 transform="ln(TPM+1) per sample, then mean and sample SD (ddof=1)",
                 observations="24 of 24 requested gene observations retained",
                 replicate_unit="three supplied samples per temperature; biological/technical status not encoded",
                 interpretation="descriptive association; connecting lines are not regression fits",
                 statistics=result)
    stem.with_suffix(".notes.json").write_text(json.dumps(notes, ensure_ascii=False, indent=2,
                                                          allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, allow_nan=False))
    print(f"Outputs: {stem}.*")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
