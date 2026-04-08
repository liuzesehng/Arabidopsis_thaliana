from __future__ import annotations

import re
from collections import OrderedDict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

INPUT_FILES = [
    Path("/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana/list/xgboot/TPM_4.5_all/A_per_TPM_import_results.txt"),
    Path("/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana/list/xgboot/TPM_4.5_all/A_TPM_import_results.txt"),
    Path("/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana/list/xgboot/TPM_4.5_all/B_per_TPM_import_results.txt"),
    Path("/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana/list/xgboot/TPM_4.5_all/B_TPM_import_results.txt"),
    Path("/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana/list/xgboot/TPM_4.5_all/total_TPM_import_results.txt"),
]

OUTPUT_ROOT = Path("/datapool/life-gongl/zesheng/Arabidopsis_thaliana/list/RCA")
GROUP_ORDER = ["SNP", "CG", "CHG", "CHH"]
FEATURE_LINE_RE = re.compile(r"^\s{4}([A-Za-z0-9_./%-]+):\s*([0-9]+(?:\.[0-9]+)?)\s*$")

plt.rcParams.update(
    {
        "font.size": 12,
        "axes.titlesize": 14,
        "axes.labelsize": 13,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
    }
)


def extract_sample_name(input_path: Path) -> str:
    suffix = "_import_results.txt"
    if not input_path.name.endswith(suffix):
        raise ValueError(f"Unexpected input filename: {input_path.name}")
    return input_path.name[: -len(suffix)]


def classify_feature(feature_name: str) -> str | None:
    if feature_name.startswith("snp_"):
        return "SNP"
    if feature_name.startswith("CG_"):
        return "CG"
    if feature_name.startswith("CHG_"):
        return "CHG"
    if feature_name.startswith("CHH_"):
        return "CHH"
    return None


def parse_feature_contributions(input_path: Path) -> list[tuple[str, float]]:
    lines = input_path.read_text(encoding="utf-8").splitlines()
    section_header = "  所有特征及其平均|SHAP|值:"
    try:
        start_index = lines.index(section_header) + 1
    except ValueError as exc:
        raise ValueError(f"Could not find SHAP feature section in {input_path}") from exc

    feature_contributions: list[tuple[str, float]] = []
    for line in lines[start_index:]:
        if not line.strip():
            break
        match = FEATURE_LINE_RE.match(line)
        if not match:
            break
        feature_name, shap_value = match.groups()
        feature_contributions.append((feature_name, float(shap_value)))

    if not feature_contributions:
        raise ValueError(f"No SHAP feature entries found in {input_path}")

    return feature_contributions


def summarize_groups(feature_contributions: list[tuple[str, float]]) -> tuple[OrderedDict[str, float], float]:
    group_sums: OrderedDict[str, float] = OrderedDict((group, 0.0) for group in GROUP_ORDER)

    for feature_name, shap_value in feature_contributions:
        group = classify_feature(feature_name)
        if group is None:
            continue
        group_sums[group] += shap_value

    total_shap = sum(value for _, value in feature_contributions)
    return group_sums, total_shap


def calculate_percentages(group_sums: OrderedDict[str, float], total_shap: float) -> OrderedDict[str, float]:
    if total_shap == 0:
        return OrderedDict((group, 0.0) for group in GROUP_ORDER)
    return OrderedDict((group, value / total_shap * 100) for group, value in group_sums.items())


def write_summary(
    output_dir: Path,
    sample_name: str,
    input_path: Path,
    feature_count: int,
    group_sums: OrderedDict[str, float],
    percentages: OrderedDict[str, float],
    total_shap: float,
) -> Path:
    output_path = output_dir / f"{sample_name}_SHAP_group_contribution.txt"
    lines = [
        f"=== {sample_name} SHAP 分组贡献汇总 ===",
        "",
        f"输入文件: {input_path}",
        f"特征数量: {feature_count}",
        f"总 mean(|SHAP|): {total_shap:.6f}",
        "",
        "分组\tSHAP总和\t百分比",
    ]

    for group in GROUP_ORDER:
        lines.append(f"{group}\t{group_sums[group]:.6f}\t{percentages[group]:.2f}%")

    lines.extend(
        [
            "",
            f"百分比合计: {sum(percentages.values()):.2f}%",
        ]
    )
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output_path


def plot_group_percentages(output_dir: Path, sample_name: str, percentages: OrderedDict[str, float]) -> tuple[Path, Path]:
    sorted_items = sorted(percentages.items(), key=lambda item: item[1], reverse=True)
    x_labels = [group for group, _ in sorted_items]
    y_values = [value for _, value in sorted_items]

    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    bars = ax.bar(x_labels, y_values, color=["#4c78a8", "#72b7b2", "#f2cf5b", "#e45756"], width=0.65)

    ax.set_xlabel("Feature Group")
    ax.set_ylabel("Percentage (%)")
    ax.set_ylim(0, 100)
    ax.grid(axis="y", linestyle=":", linewidth=0.8, alpha=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    for bar, value in zip(bars, y_values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            min(value + 1.2, 99.0),
            f"{value:.2f}%",
            ha="center",
            va="bottom",
            fontsize=11,
        )

    fig.tight_layout()

    pdf_path = output_dir / f"{sample_name}_SHAP_group_contribution.pdf"
    png_path = output_dir / f"{sample_name}_SHAP_group_contribution.png"
    fig.savefig(pdf_path, format="pdf", bbox_inches="tight")
    fig.savefig(png_path, format="png", bbox_inches="tight")
    plt.close(fig)
    return pdf_path, png_path


def process_file(input_path: Path) -> None:
    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    sample_name = extract_sample_name(input_path)
    output_dir = OUTPUT_ROOT / sample_name
    output_dir.mkdir(parents=True, exist_ok=True)

    feature_contributions = parse_feature_contributions(input_path)
    group_sums, total_shap = summarize_groups(feature_contributions)
    percentages = calculate_percentages(group_sums, total_shap)

    summary_path = write_summary(
        output_dir=output_dir,
        sample_name=sample_name,
        input_path=input_path,
        feature_count=len(feature_contributions),
        group_sums=group_sums,
        percentages=percentages,
        total_shap=total_shap,
    )
    pdf_path, png_path = plot_group_percentages(output_dir, sample_name, percentages)

    print(f"[{sample_name}] features={len(feature_contributions)} total={total_shap:.6f}")
    print(f"[{sample_name}] summary={summary_path}")
    print(f"[{sample_name}] plots={pdf_path}, {png_path}")


def main() -> None:
    for input_path in INPUT_FILES:
        process_file(input_path)


if __name__ == "__main__":
    main()
