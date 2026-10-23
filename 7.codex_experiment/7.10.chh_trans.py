#!/usr/bin/env python3
"""CHH-associated non-coexpressed TFs versus five RCA expression traits.

Run without arguments for project defaults, or use --validate-only to calculate
without writing outputs. Raw replicates are averaged separately by temperature;
no individual pairing, transformation, imputation or P-value adjustment is used.
Requires numpy, pandas, scipy, matplotlib, PyMuPDF and nature-figure QA tools.
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

import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata, spearmanr

PROJECT = Path(__file__).resolve().parents[2]
DATA = PROJECT / "list"
OUTPUT = DATA / "coex_experiment"
TF_LIST = DATA / "meth_experiment/CHH_TF_coexloci/CHH_TF_plot_data.tsv"
EXPRESSION = OUTPUT / "RCA.tran.txt"
RCA = DATA / "RCA/RCA.rnadata.tsv"
PREFIX = "CHH_TF_RCA_correlation"
QA_TOOLS = Path("/datapool/home/2023102768/.codex/skills/nature-figure/scripts")
TEMPERATURES = tuple(range(4, 33, 4))
TF_GENES = {
    "LHY": "AT1G01060", "EPR1": "AT1G18330", "RVE7L": "AT3G10113",
    "RVE5": "AT4G01280", "LCL1": "AT5G02840", "RVE1": "AT5G17300",
    "RVE6": "AT5G52660",
}
TRAITS = ("total", "α", "β", "α%", "β%")
SOURCE_TRAITS = ("total", "α", "β", "α/%", "β/%")
METHOD = "two-sided exact Spearman pairing permutation (40320 permutations)"


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


def stars(p):
    return "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else ""


def read_tsv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if len(rows) < 2 or not all(rows[0]) or len(set(rows[0])) != len(rows[0]):
        raise ValueError(f"{path}: empty table, blank or duplicate headers")
    if any(len(row) != len(rows[0]) for row in rows[1:]):
        raise ValueError(f"{path}: inconsistent TSV field counts")
    return pd.DataFrame(rows[1:], columns=rows[0])


def require_columns(frame, columns, label):
    absent = set(columns) - set(frame.columns)
    if absent:
        raise ValueError(f"{label}: missing required columns {sorted(absent)}")


def numeric(values, label, percentage=False):
    result = pd.to_numeric(values, errors="coerce").astype(float)
    if not np.isfinite(result).all():
        raise ValueError(f"{label}: missing/non-finite measurement; no imputation permitted")
    if (result < 0).any() or (percentage and (result > 100).any()):
        raise ValueError(f"{label}: measurement outside allowed range")
    return result


def check_design(temperatures, label):
    counts = pd.Series(temperatures).value_counts().sort_index().to_dict()
    if counts != {t: 3 for t in TEMPERATURES}:
        raise ValueError(f"{label}: expected eight temperatures with three replicates each; got {counts}")


def load_inputs(tf_list, expression_file, rca_file):
    candidates = read_tsv(tf_list)
    require_columns(candidates, ("symbol", "tf_id", "group"), "TF list")
    excluded = candidates.tf_id.eq("AT3G09600") | candidates.symbol.eq("LCL5")
    selected = candidates.loc[~excluded, ["symbol", "tf_id", "group"]].copy()
    # Multiple motif hits can identify the same gene; preserve first occurrence.
    selected = selected.drop_duplicates().reset_index(drop=True)
    if selected.tf_id.duplicated().any() or selected.symbol.duplicated().any():
        raise ValueError("TF list: conflicting gene/symbol/group mappings")
    if list(selected.symbol) != list(TF_GENES) or dict(zip(selected.symbol, selected.tf_id)) != TF_GENES:
        raise ValueError("TF list: expected LHY, EPR1, RVE7L, RVE5, LCL1, RVE1, RVE6 in source order")
    if not selected.group.eq("non_coexpressed").all():
        raise ValueError("TF list: selected genes must be marked non_coexpressed")

    expression = read_tsv(expression_file)
    if expression.columns[0] != "Name":
        raise ValueError("TF expression: first column must be Name")
    samples = list(expression.columns[1:])
    if any(re.fullmatch(r"\d+_[1-9]\d*", sample) is None for sample in samples):
        raise ValueError("TF expression: sample IDs must be <temperature>_<positive replicate>")
    sample_keys = [tuple(map(int, sample.split("_"))) for sample in samples]
    if len(set(sample_keys)) != len(sample_keys):
        raise ValueError("TF expression: duplicate normalized sample IDs")
    temperatures = [key[0] for key in sample_keys]
    check_design(temperatures, "TF expression")
    tf = pd.DataFrame({"Sample": samples, "Temperature": temperatures})
    for row in selected.itertuples(index=False):
        matches = expression.loc[expression.Name.eq(row.tf_id)]
        if len(matches) != 1:
            raise ValueError(f"{row.symbol} ({row.tf_id}): expected one expression row, found {len(matches)}")
        tf[row.symbol] = numeric(matches.iloc[0, 1:], row.symbol).to_numpy()

    rca = read_tsv(rca_file)
    require_columns(rca, ("Temperature", *SOURCE_TRAITS), "RCA")
    rca = rca[["Temperature", *SOURCE_TRAITS]].copy()
    rca["Temperature"] = numeric(rca.Temperature, "RCA temperature")
    check_design(rca.Temperature, "RCA")
    rca["Temperature"] = rca.Temperature.astype(int)
    rca.rename(columns={"α/%": "α%", "β/%": "β%"}, inplace=True)
    for trait in TRAITS:
        rca[trait] = numeric(rca[trait], f"RCA {trait}", percentage=trait.endswith("%"))
    provenance = dict(candidate_rows=len(candidates), excluded_LCL5_rows=int(excluded.sum()),
                      selected_genes=len(selected), duplicate_motif_rows_collapsed=int((~excluded).sum())-len(selected),
                      expression_rows=len(expression), selected_expression_rows=len(selected),
                      tf_replicates=len(tf), rca_replicates=len(rca))
    return selected, tf, rca, provenance


def summarize(tf, rca):
    tf_means = tf.groupby("Temperature", sort=True)[list(TF_GENES)].mean().reindex(TEMPERATURES)
    rca_means = rca.groupby("Temperature", sort=True)[list(TRAITS)].mean().reindex(TEMPERATURES)
    means = tf_means.join(rca_means, validate="one_to_one")
    records = []
    for source, frame, variables in [("TF", tf, tuple(TF_GENES)), ("RCA", rca, TRAITS)]:
        for variable in variables:
            if means[variable].nunique() < 2:
                raise ValueError(f"{source} {variable}: constant temperature-mean vector; Spearman undefined")
            for temperature in TEMPERATURES:
                values = frame.loc[frame.Temperature.eq(temperature), variable]
                records.append(dict(Source=source, TF_ID=TF_GENES.get(variable, ""), Variable=variable,
                                    Temperature=temperature, Mean=float(means.loc[temperature, variable]),
                                    n_valid=int(values.count()), Unit="%" if variable.endswith("%") else "input expression units"))
    if not np.isfinite(means.to_numpy()).all():
        raise ValueError("all eight means must be finite")
    return means, pd.DataFrame(records)


def exact_spearman(x, y, permutations):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if x.shape != (8,) or y.shape != (8,) or not np.isfinite([x, y]).all():
        raise ValueError("exact Spearman requires two finite vectors of length eight")
    rx, ry = rankdata(x, method="average") - 4.5, rankdata(y, method="average") - 4.5
    denominator = np.linalg.norm(rx) * np.linalg.norm(ry)
    if denominator == 0:
        raise ValueError("constant temperature-mean variable: Spearman undefined")
    observed = float(rx @ ry)
    null = np.einsum("ij,j->i", rx[permutations], ry)
    extreme = int(np.count_nonzero(np.abs(null) >= abs(observed) - 1e-12))
    return float(np.clip(observed / denominator, -1, 1)), extreme / len(permutations), extreme


def compute_statistics(means):
    permutations = np.asarray(list(itertools.permutations(range(8))), dtype=np.int16)
    records = []
    for symbol, gene in TF_GENES.items():
        for trait in TRAITS:
            rho, p, extreme = exact_spearman(means[symbol], means[trait], permutations)
            if not np.isclose(rho, spearmanr(means[symbol], means[trait]).statistic, atol=1e-14, rtol=0):
                raise AssertionError("independent scipy rho check failed")
            records.append(dict(TF_ID=gene, Symbol=symbol, Trait=trait, Spearman_rho=rho,
                                P_value=p, n=8, Stars=stars(p), Method=METHOD, Adjustment="none",
                                Extreme_permutations=extreme, Total_permutations=40320))
    return pd.DataFrame(records)


def load_alignment_tools(directory):
    if not Path(directory / 'audit_panel_alignment.py').is_file():
        _skip_optional_qa(f"Optional QA tool not found: {directory / 'audit_panel_alignment.py'}")
        return None
    spec = importlib.util.spec_from_file_location("chh_tf_rca_alignment", directory / "audit_panel_alignment.py")
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load alignment auditor")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    try:
        spec.loader.exec_module(module)
    except ImportError as exc:
        _skip_optional_qa(f"Optional alignment dependency unavailable: {exc}")
        return None
    return module


def make_figure(statistics, auditor, qa_dir):
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "chh-tf-rca-mpl"))
    os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "chh-tf-rca-cache"))
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
    # Explicit inches preserve the requested 180 x 124 mm page in static QA.
    fig = plt.figure(figsize=(7.086614173228347, 4.881889763779528), facecolor="white")

    def add_axes(x, y, w, h, gid):
        ax = fig.add_axes([x / 180, y / 124, w / 180, h / 124])
        ax.set_gid(gid)
        return ax

    heatmap = add_axes(29, 27, 128, 87, "heatmap")
    lookup = statistics.set_index(["Symbol", "Trait"])
    matrix = np.array([[lookup.loc[(s, t), "Spearman_rho"] for t in TRAITS] for s in TF_GENES])
    cmap = LinearSegmentedColormap.from_list("TF_RCA_diverging", ["#35658A", "#FAFAFA", "#AE504B"])
    mesh = heatmap.pcolormesh(np.arange(6), np.arange(8), matrix, cmap=cmap,
                             norm=Normalize(vmin=-1, vmax=1), edgecolors="white", linewidth=.6, rasterized=False)
    heatmap.set(xlim=(0, 5), ylim=(7, 0), xticks=np.arange(5) + .5, yticks=np.arange(7) + .5,
                xticklabels=TRAITS, yticklabels=list(TF_GENES))
    heatmap.tick_params(length=0, pad=5)
    for spine in heatmap.spines.values():
        spine.set_visible(False)
    for i, symbol in enumerate(TF_GENES):
        for j, trait in enumerate(TRAITS):
            row = lookup.loc[(symbol, trait)]
            heatmap.text(j + .5, i + .5, f"{row.Spearman_rho:.2f}{row.Stars}",
                         ha="center", va="center", fontsize=9,
                         color="white" if abs(row.Spearman_rho) >= .72 else "#202020")
    cax = add_axes(161, 45, 3, 52, "colorbar")
    cb = fig.colorbar(mesh, cax=cax, ticks=[-1, -.5, 0, .5, 1])
    cb.solids.set_rasterized(False)
    cb.set_label("Spearman ρ", fontsize=8, labelpad=5)
    cb.ax.tick_params(labelsize=7.5, length=2, width=.6)
    cb.outline.set_linewidth(.6)
    notes = ["Cells show Spearman ρ across eight temperature means (n = 8).",
             "Two-sided exact permutation P values; no multiple-testing correction.",
             "* P < 0.05; ** P < 0.01; *** P < 0.001."]
    for note, bottom_mm in zip(notes, [15, 10.5, 6]):
        fig.text(10 / 180, bottom_mm / 124, note, fontsize=8, ha="left", va="baseline")
    fig.canvas.draw()
    if auditor is None:
        _skip_optional_qa("audit_panel_alignment.py is unavailable", qa_dir / 'alignment.json')
    else:
        auditor.require_matplotlib_panel_alignment(fig, json_out=qa_dir / "alignment.json",
            tolerance_pt=1.5, axes=[heatmap], panel_ids=["heatmap"], strict=True)
    if matrix.shape != (7, 5) or len(heatmap.texts) != 35:
        raise AssertionError("expected exactly 35 annotated cells")
    for text, row in zip(heatmap.texts, statistics.itertuples(index=False)):
        if text.get_text() != f"{row.Spearman_rho:.2f}{row.Stars}":
            raise AssertionError("heatmap annotation does not match statistics")
    return fig, matplotlib.__version__


def build_caption():
    genes = "、".join(f"{s}（{g}）" for s, g in TF_GENES.items())
    return (
        "# CHH关联候选转录因子编码基因与RCA表达性状的相关性\n\n"
        "热图比较候选转录因子（TF）编码基因的表达量与RCA总表达量、两种亚型表达量及其组成比例的温度组均值相关性。"
        "候选名单共包含8个TF，其中已归为RCA共表达的LCL5（AT3G09600）按预设方案排除；"
        "保留原名单中归为非共表达的7个TF。此分类沿用候选名单，并非由本次相关检验重新定义。"
        f"行自上而下为{genes}；列自左向右为RCA总表达量（total）、α和β亚型表达量（α、β），"
        "以及各亚型占RCA总表达量的百分比（α%、β%）。行列保持预设顺序，不进行聚类。\n\n"
        "TF与RCA数据均覆盖4、8、12、16、20、24、28和32 ℃共8个温度，各温度各有3个生物学重复，"
        "每份表达数据包含24个重复。对两份数据分别按温度计算原始尺度算术均值，再按温度对应均值；"
        "表达量不进行对数转换、标准化或缺失值填补，0值作为有效观测保留。"
        "比例直接对各重复的输入百分比求算术均值，不根据α、β和total的温度组均值重新计算。"
        "RCA输入未提供个体样本编号，因此不推断两份表之间的个体配对。每个温度组权重相同，"
        "每次相关检验的统计样本量为n = 8，而非24个重复。\n\n"
        "共计算7 × 5 = 35项Spearman秩相关。单元格数值为相关系数ρ，显示时保留两位小数；"
        "蓝—白—红色标固定为−1至1，蓝色表示负相关，红色表示正相关，接近白色表示相关系数接近0。"
        "正ρ表示较高的TF表达均值倾向于对应较高的RCA性状均值，负ρ表示相反趋势。"
        "双侧精确配对置换检验穷举8! = 40,320种配对，P值为置换后|ρ|不小于观察|ρ|的配对比例；"
        "并列观测采用平均秩。所有检验报告未经多重比较校正的原始P值，不进行BH校正。"
        "星号依据未舍入的P值：*，0.01 ≤ P < 0.05；**，0.001 ≤ P < 0.01；***，P < 0.001；"
        "P ≥ 0.05时不加星号。完整精度结果及温度均值见配套统计表和均值表。\n\n"
        "本图提供温度组均值层面的探索性关联证据。未经多重比较校正的显著性应谨慎解读；"
        "相关性不证明TF直接调控RCA，也不能排除共同温度响应造成的关联。"
        "α%与β%近似互补，其相反方向的相关结果不构成相互独立的重复证据。\n"
    )


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def audit(command, destination):
    if _skip_missing_qa_command(command, destination):
        return
    result = subprocess.run(command, capture_output=True, text=True)
    destination.write_text(result.stdout + result.stderr, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(f"QA failed: {destination}\n{result.stdout[-3500:]}\n{result.stderr[-1500:]}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tf-list", type=Path, default=TF_LIST)
    parser.add_argument("--expression", type=Path, default=EXPRESSION)
    parser.add_argument("--rca", type=Path, default=RCA)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT)
    parser.add_argument("--prefix", default=PREFIX)
    parser.add_argument("--qa-tools", type=Path, default=QA_TOOLS)
    parser.add_argument("--validate-only", action="store_true", help="validate and calculate without writing files")
    args = parser.parse_args(argv)
    if not args.prefix or Path(args.prefix).name != args.prefix or args.prefix in (".", ".."):
        parser.error("--prefix must be a filename stem")
    inputs = {name: getattr(args, name) for name in ("tf_list", "expression", "rca")}
    hashes = {name: sha256(path) for name, path in inputs.items()}
    _, tf, rca, provenance = load_inputs(args.tf_list, args.expression, args.rca)
    means, summary = summarize(tf, rca)
    statistics = compute_statistics(means)
    if args.validate_only:
        print("PASS: 7 TFs; 24 replicates per input; 8 temperature means; 35 exact Spearman tests")
        print(statistics.to_string(index=False))
        return 0
    auditor = load_alignment_tools(args.qa_tools)
    # PDF QA dependencies are checked only when the optional tools run.
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
    statistics.to_csv(str(stem) + ".statistics.tsv", sep="\t", index=False, float_format="%.15g")
    summary.to_csv(str(stem) + ".temperature_means.tsv", sep="\t", index=False, float_format="%.15g")
    Path(str(stem) + ".caption.md").write_text(build_caption(), encoding="utf-8")
    if any(sha256(path) != hashes[name] for name, path in inputs.items()):
        raise RuntimeError("input changed during execution")
    report = dict(status=("automated_QA_INCOMPLETE" if _OPTIONAL_QA_SKIPPED else "automated_QA_PASS"), visual_review="pending", input_sha256=hashes,
                  script_sha256=sha256(__file__), parameters={k: str(v) for k, v in vars(args).items()},
                  provenance=provenance, temperature_groups=8, statistical_tests=35,
                  transform="identity before arithmetic mean", pairing="temperature means only",
                  test=METHOD, adjustment="none", excluded_gene="LCL5 (AT3G09600), user-specified",
                  archetype="quantitative grid", reuse="structural adaptation of user-provided heatmap",
                  question="How do seven candidate TF expression profiles associate with five RCA traits across temperature means?",
                  panels=["heatmap"], heatmap_cells=35, figure_mm=[180, 124], png_dpi=600,
                  versions=dict(python=sys.version.split()[0], numpy=np.__version__, pandas=pd.__version__,
                                scipy=scipy.__version__, matplotlib=mpl_version))
    (qa_dir / "summary.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"PASS: {stem}.{{pdf,svg,png}}; 35 tests; inputs unchanged")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
