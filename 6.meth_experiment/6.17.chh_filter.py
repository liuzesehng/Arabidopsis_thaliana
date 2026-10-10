#!/usr/bin/env python3
"""Six CHH candidates: temperature-mean Spearman correlations and thesis heatmap.

Run without arguments for the project defaults. Requires numpy, pandas, scipy,
matplotlib and, for publication audits, PyMuPDF plus the nature-figure QA tools.
All statistics use raw-scale means and exact two-sided pairing permutations.
No sample matching across external files, log transforms or imputation is used.
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

import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata, spearmanr

PROJECT = Path(__file__).resolve().parents[2]
OUTPUT = PROJECT / "list/meth_experiment"
INPUT = OUTPUT / "Alt.all.snp_meth.filtered.tsv"
PREFIX = "CHH_candidates_temperature_RCA"
QA_TOOLS = Path("/datapool/home/2023102768/.codex/skills/nature-figure/scripts")
SITES = ("CHH_1634", "CHH_1635", "CHH_1651", "CHH_1676", "CHH_1698", "CHH_1771")
TEMPERATURES = tuple(range(4, 33, 4))
TRAITS = ("Temperature", "total", "α", "β", "α%", "β%")
SOURCE_TRAITS = ("total", "α", "β", "α/%", "β/%")
METHOD = "two-sided exact Spearman pairing permutation (40320 permutations)"



def stars(p: float) -> str:
    return "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else ""


def load_data(path: Path) -> pd.DataFrame:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if len(rows) < 2 or len(set(rows[0])) != len(rows[0]):
        raise ValueError("empty input or duplicate headers")
    if any(len(r) != len(rows[0]) for r in rows[1:]):
        raise ValueError("inconsistent TSV field counts")
    required = ["Temperature", *SITES, *SOURCE_TRAITS]
    absent = set(required) - set(rows[0])
    if absent:
        raise ValueError(f"missing required columns: {sorted(absent)}")
    data = pd.DataFrame(rows[1:], columns=rows[0])[required].copy()
    data.rename(columns={"Temperature": "Sample", "α/%": "α%", "β/%": "β%"}, inplace=True)
    parsed = data.Sample.str.extract(r"^T(\d+)-(\d+)$")
    if parsed.isna().any().any() or data.Sample.duplicated().any():
        raise ValueError("sample IDs must be unique T<temperature>-<replicate> values")
    data["Temperature"] = parsed[0].astype(int)
    data["Replicate"] = parsed[1].astype(int)
    if (data.Replicate < 1).any() or data.duplicated(["Temperature", "Replicate"]).any():
        raise ValueError("invalid or repeated temperature/replicate identifiers")
    counts = data.groupby("Temperature").size()
    if tuple(counts.index) != TEMPERATURES or not counts.between(2, 3).all():
        raise ValueError(f"expected eight temperatures, each with 2–3 replicates: {counts.to_dict()}")
    for variable in [*SITES, *TRAITS[1:]]:
        values = pd.to_numeric(data[variable], errors="coerce").astype(float)
        if not np.isfinite(values).all():
            raise ValueError(f"{variable}: missing or non-finite measurement; no imputation or reduced n permitted")
        if (values < 0).any() or (variable in (*SITES, "α%", "β%") and (values > 100).any()):
            raise ValueError(f"{variable}: measurement outside the allowed range")
        data[variable] = values
    return data.sort_values(["Temperature", "Replicate"], kind="stable").reset_index(drop=True)


def summarize(data: pd.DataFrame):
    variables = [*SITES, *TRAITS[1:]]
    means = data.groupby("Temperature", sort=True)[variables].mean().reindex(TEMPERATURES)
    if not np.isfinite(means.to_numpy()).all():
        raise ValueError("all eight finite means are required")
    means["Temperature"] = means.index.astype(float)
    records = []
    for temperature, block in data.groupby("Temperature", sort=True):
        for variable in variables:
            records.append(dict(Temperature=int(temperature), Variable=variable,
                                Mean=float(means.loc[temperature, variable]),
                                n_valid=int(block[variable].count()),
                                Unit="%" if variable in (*SITES, "α%", "β%") else "input expression units"))
    return means, pd.DataFrame(records)


def exact_spearman(x, y, permutations):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if x.shape != (8,) or y.shape != (8,) or not np.isfinite([x, y]).all():
        raise ValueError("exact Spearman requires two finite vectors of length eight")
    rx, ry = rankdata(x, method="average") - 4.5, rankdata(y, method="average") - 4.5
    denominator = np.linalg.norm(rx) * np.linalg.norm(ry)
    if denominator == 0:
        raise ValueError("constant temperature-mean variable: Spearman is undefined")
    dot = float(rx @ ry)
    null = np.einsum("ij,j->i", rx[permutations], ry)
    # Average ranks are integer/half-integer, so this tolerance only absorbs roundoff.
    extreme = int(np.count_nonzero(np.abs(null) >= abs(dot) - 1e-12))
    return float(np.clip(dot / denominator, -1, 1)), extreme / len(permutations), extreme


def compute_statistics(means: pd.DataFrame):
    permutations = np.asarray(list(itertools.permutations(range(8))), dtype=np.int16)
    records = []
    for site in SITES:
        for trait in TRAITS:
            rho, p, extreme = exact_spearman(means[site], means[trait], permutations)
            if not np.isclose(rho, spearmanr(means[site], means[trait]).statistic, atol=1e-14):
                raise AssertionError("independent scipy rho check failed")
            records.append(dict(Site=site, Trait=trait, Spearman_rho=rho, P_value=p,
                                n=8, Stars=stars(p), Method=METHOD, Adjustment="none",
                                Extreme_permutations=extreme, Total_permutations=40320))
    return pd.DataFrame(records)


def load_alignment_tools(directory):
    path = directory / "audit_panel_alignment.py"
    spec = importlib.util.spec_from_file_location("chh_candidates_alignment", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load alignment auditor: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def make_figure(statistics, auditor, qa_dir):
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "chh-candidates-mpl"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, Normalize

    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
                         "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7.5,
                         "ytick.labelsize": 7.5, "axes.linewidth": .6,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "svg.fonttype": "none",
                         "savefig.facecolor": "white", "legend.frameon": False})
    fig = plt.figure(figsize=(180 / 25.4, 124 / 25.4), facecolor="white")

    def add_axes(x, y, w, h, gid):
        ax = fig.add_axes([x / 180, y / 124, w / 180, h / 124])
        ax.set_gid(gid)
        return ax

    heatmap = add_axes(29, 27, 128, 87, "heatmap")
    lookup = statistics.set_index(["Site", "Trait"])
    matrix = np.array([[lookup.loc[(s, t), "Spearman_rho"] for t in TRAITS] for s in SITES])
    cmap = LinearSegmentedColormap.from_list("CHH_diverging", ["#35658A", "#FAFAFA", "#AE504B"])
    norm = Normalize(vmin=-1, vmax=1)
    mesh = heatmap.pcolormesh(np.arange(7), np.arange(7), matrix, cmap=cmap, norm=norm,
                             edgecolors="white", linewidth=.6, rasterized=False)
    heatmap.set(xlim=(0, 6), ylim=(6, 0), xticks=np.arange(6) + .5,
                yticks=np.arange(6) + .5, xticklabels=TRAITS, yticklabels=SITES)
    heatmap.tick_params(length=0, pad=5)
    for spine in heatmap.spines.values():
        spine.set_visible(False)
    for i, site in enumerate(SITES):
        for j, trait in enumerate(TRAITS):
            row = lookup.loc[(site, trait)]
            color = "white" if abs(row.Spearman_rho) >= .72 else "#202020"
            heatmap.text(j + .5, i + .5, f"{row.Spearman_rho:.2f}{row.Stars}",
                         ha="center", va="center", fontsize=9, color=color)
    cax = add_axes(161, 45, 3, 52, "colorbar")
    cb = fig.colorbar(mesh, cax=cax, ticks=[-1, -.5, 0, .5, 1])
    cb.solids.set_rasterized(False)
    cb.set_label("Spearman ρ", fontsize=8, labelpad=5)
    cb.ax.tick_params(labelsize=7.5, length=2, width=.6)
    cb.outline.set_linewidth(.6)
    notes = [
        "Cells show Spearman ρ across eight temperature means (n = 8).",
        "Two-sided exact permutation P values; no multiple-testing correction.",
        "* P < 0.05; ** P < 0.01; *** P < 0.001.",
    ]
    for note, bottom_mm in zip(notes, [15, 10.5, 6]):
        fig.text(10 / 180, bottom_mm / 124, note, fontsize=8, ha="left", va="baseline")
    fig.canvas.draw()
    auditor.require_matplotlib_panel_alignment(
        fig, json_out=qa_dir / "alignment.json", tolerance_pt=1.5,
        axes=[heatmap], panel_ids=["heatmap"], strict=True)
    if matrix.shape != (6, 6) or len(heatmap.texts) != 36:
        raise AssertionError("expected exactly 36 annotated heatmap cells")
    return fig, matplotlib.__version__


def audit(command, destination):
    result = subprocess.run(command, capture_output=True, text=True)
    destination.write_text(result.stdout + result.stderr, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(f"QA failed: {destination}\n{result.stdout[-3500:]}\n{result.stderr[-1500:]}")


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_caption(data: pd.DataFrame) -> str:
    counts = "、".join(f"{t} ℃（{n}个）" for t, n in data.groupby("Temperature").size().items())
    return (
        "# 六个候选CHH甲基化位点与温度及RCA表达性状的Spearman相关性\n\n"
        "热图展示六个候选CHH甲基化位点的甲基化水平与温度及RCA表达性状之间的Spearman秩相关关系。"
        "CHH表示非对称胞嘧啶甲基化序列背景，其中H为A、C或T。"
        "行自上而下依次为CHH_1634、CHH_1635、CHH_1651、CHH_1676、CHH_1698和CHH_1771；"
        "列自左向右依次为温度（Temperature）、RCA总表达量（total）、RCA α和β亚型的表达量（α、β），"
        "以及两种亚型各自占RCA总表达量的比例（α%、β%）。表达性状依据RNA-seq数据获得。"
        "图中行列保持预设顺序，未进行层次聚类。\n\n"
        "分析包含4、8、12、16、20、24、28和32 ℃共8个温度处理。"
        f"各温度纳入的生物学重复数为：{counts}，合计{len(data)}个重复。"
        "首先分别计算每个温度下同一CHH位点甲基化水平及各RCA性状的生物学重复算术均值，"
        "再按温度对应这些均值进行相关分析。甲基化水平以百分数表示；"
        "total、α和β保留输入表达量的原始尺度，不作对数转换。"
        "α%和β%直接对输入表中各重复的比例取算术均值，不以α、β和total的组均值重新计算比例。"
        "每次检验均纳入全部8个温度组，每组在相关分析中具有相同权重；"
        "因此统计样本量为n = 8，而不是生物学重复总数，也未将重复作为独立的相关分析点。\n\n"
        "每个单元格对应一个位点与一项性状的相关关系，共进行6 × 6 = 36项检验。"
        "格内数值为Spearman相关系数ρ，显示时保留两位小数；计算检验及显著性标记时使用未舍入的数值。"
        "颜色同样表示ρ，色标统一固定为−1至1：红色表示正相关，蓝色表示负相关，"
        "接近白色表示ρ接近0，颜色越深表示单调相关的绝对强度越大。"
        "其中，正ρ表示较高甲基化水平倾向于对应较高温度或较大的相应RCA性状值，负ρ表示相反趋势。\n\n"
        "P值采用双侧精确配对置换检验获得。对每次检验穷举8! = 40,320种配对，"
        "以置换后相关系数的绝对值不小于观察值绝对值的配对数占全部配对数的比例计算P值；"
        "存在并列观测值时采用平均秩。按候选位点分别开展相关分析，"
        "本图报告未经多重比较校正的原始P值，不进行Benjamini–Hochberg（BH）校正。"
        "星号表示原始P值所在区间：*，0.01 ≤ P < 0.05；**，0.001 ≤ P < 0.01；***，P < 0.001。"
        "P ≥ 0.05时仍显示ρ，但不加星号。完整的相关系数、精确P值及样本量见配套统计结果表。\n\n"
        "本图反映温度组均值层面的相关性，不能据此推断甲基化变化对RCA表达的因果作用，"
        "也不能排除二者共同响应温度所产生的关联。α%与β%作为两种亚型的组成比例，"
        "其相关结果应结合比例之间的互补关系理解，不宜作为相互独立的重复证据。\n"
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT)
    parser.add_argument("--prefix", default=PREFIX)
    parser.add_argument("--validate-only", action="store_true", help="validate and calculate without writing files")
    parser.add_argument("--qa-tools", type=Path, default=QA_TOOLS)
    args = parser.parse_args(argv)
    if not args.prefix or Path(args.prefix).name != args.prefix or args.prefix in (".", ".."):
        parser.error("--prefix must be a filename stem, not a path")
    input_hash = sha256(args.input)
    data = load_data(args.input)
    means, summary = summarize(data)
    statistics = compute_statistics(means)
    if args.validate_only:
        print(f"PASS: {len(data)} biological replicates; 8 temperature means; 36 exact Spearman tests")
        print(data.groupby("Temperature").size().to_string())
        print(statistics.to_string(index=False))
        return 0
    auditor = load_alignment_tools(args.qa_tools)
    if importlib.util.find_spec("pymupdf") is None and importlib.util.find_spec("fitz") is None:
        raise RuntimeError("NOT AUDITABLE: install PyMuPDF for PDF QA")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    qa_dir = args.output_dir / f"{args.prefix}_QA"
    qa_dir.mkdir(exist_ok=True)
    stem = args.output_dir / args.prefix
    audit([sys.executable, str(args.qa_tools / "validate_figure.py"), __file__], qa_dir / "source-preflight.txt")
    fig, mpl_version = make_figure(statistics, auditor, qa_dir)
    try:
        fig.savefig(str(stem) + ".pdf")
        audit([sys.executable, str(args.qa_tools / "audit_pdf_text.py"), str(stem) + ".pdf",
               "--min-pt", "7", "--json"], qa_dir / "pdf-text.json")
        audit([sys.executable, str(args.qa_tools / "audit_figure_collisions.py"), str(stem) + ".pdf",
               "--strict", "--json-out", str(qa_dir / "collision-audit.json")], qa_dir / "collision-audit.txt")
        fig.savefig(str(stem) + ".svg")
        fig.savefig(str(stem) + ".png", dpi=600)
    finally:
        import matplotlib.pyplot as plt
        plt.close(fig)
    for suffix, table in [("statistics", statistics), ("temperature_means", summary)]:
        table.to_csv(str(stem) + f".{suffix}.tsv", sep="\t", index=False, float_format="%.15g")
    caption = build_caption(data)
    Path(str(stem) + ".caption.md").write_text(caption, encoding="utf-8")
    if sha256(args.input) != input_hash:
        raise RuntimeError("input changed during execution")
    report = dict(status="automated_QA_PASS", visual_review="pending", input_sha256=input_hash,
                  script_sha256=sha256(Path(__file__)), parameters={k: str(v) for k, v in vars(args).items()},
                  biological_replicates=len(data), temperature_groups=8, statistical_tests=36,
                  input_rows_used=len(data), excluded_rows=0,
                  columns_used=["Temperature", *SITES, *SOURCE_TRAITS],
                  other_columns="not part of the requested six-site/five-trait analysis",
                  transform="identity before arithmetic mean", test=METHOD, adjustment="none",
                  panels=["heatmap"], heatmap_cells=36, figure_mm=[180, 124], png_dpi=600,
                  versions=dict(python=sys.version.split()[0], numpy=np.__version__, pandas=pd.__version__,
                                scipy=scipy.__version__, matplotlib=mpl_version))
    (qa_dir / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"PASS: {stem}.{{pdf,svg,png}}; 36 tests; input unchanged")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
