#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, t as student_t


BASE_DIR = Path(
    "/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana"
)
DEFAULT_DATA_FILE = BASE_DIR / "list" / "RCA" / "Alt.all.snp_meth.filtered.tsv"
DEFAULT_IMPORT_DIR = BASE_DIR / "list" / "xgboot" / "TPM_4.5_all"
DEFAULT_OUTPUT_ROOT = BASE_DIR / "list" / "RCA"
PANEL_BASENAME = "{model}_top10_CHH_CHG_feature_correlation_2panel"
SUMMARY_BASENAME = "{model}_feature_correlation_summary.tsv"

MODEL_CONFIGS = [
    ("A_TPM", "A_TPM_import_results.txt", "α"),
    ("A_per_TPM", "A_per_TPM_import_results.txt", "α/%"),
    ("B_TPM", "B_TPM_import_results.txt", "β"),
    ("B_per_TPM", "B_per_TPM_import_results.txt", "β/%"),
    ("total_TPM", "total_TPM_import_results.txt", "total"),
]

FEATURE_LINE_PATTERN = re.compile(r"^\s*([A-Za-z0-9_/%αβ]+):\s*([-+0-9.eE]+)\s*$")
LOG1P_TARGETS = {"total", "α", "β"}

plt.rcParams.update(
    {
        "font.size": 12,
        "axes.labelsize": 14,
        "xtick.labelsize": 11,
        "ytick.labelsize": 11,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "axes.facecolor": "white",
        "figure.facecolor": "white",
    }
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot feature-target correlations for top CHH/CHG features."
    )
    parser.add_argument("--data-file", default=str(DEFAULT_DATA_FILE))
    parser.add_argument("--import-dir", default=str(DEFAULT_IMPORT_DIR))
    parser.add_argument("--output-root", default=str(DEFAULT_OUTPUT_ROOT))
    return parser.parse_args()


def parse_import_results(import_file):
    features = []
    in_feature_block = False

    with open(import_file, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")

            if "所有特征及其平均|SHAP|值:" in line:
                in_feature_block = True
                continue

            if not in_feature_block:
                continue

            match = FEATURE_LINE_PATTERN.match(line)
            if match:
                features.append((match.group(1), float(match.group(2))))
                continue

            if line.strip() == "":
                if features:
                    break
                continue

            if features:
                break

    if not features:
        raise ValueError(f"未能从 {import_file} 解析到特征重要性列表。")

    return features


def select_top_context_features(feature_pairs, top_n=10):
    top_features = feature_pairs[:top_n]
    chh_feature = next((name for name, _ in top_features if name.startswith("CHH_")), None)
    chg_feature = next((name for name, _ in top_features if name.startswith("CHG_")), None)
    return top_features, chh_feature, chg_feature


def compute_correlations(dataframe, feature_name, target_name):
    original_size = len(dataframe)
    subset = dataframe[[feature_name, target_name]].apply(pd.to_numeric, errors="coerce").dropna()

    if target_name in LOG1P_TARGETS:
        subset[target_name] = np.log1p(subset[target_name])

    n_samples = len(subset)

    if n_samples <= 2:
        raise ValueError(
            f"{feature_name} 与 {target_name} 的有效样本数不足，当前仅 {n_samples} 个。"
        )

    spearman_rho, spearman_p = spearmanr(subset[feature_name], subset[target_name])

    return {
        "data": subset,
        "n_samples": n_samples,
        "n_removed": original_size - n_samples,
        "spearman_rho": spearman_rho,
        "spearman_p": spearman_p,
    }


def format_p_value(value):
    if pd.isna(value):
        return "NA"
    if value < 1e-4:
        return f"{value:.2e}"
    return f"{value:.4f}"


def format_target_label(target_name):
    if target_name in LOG1P_TARGETS:
        return f"{target_name}(ln(TPM+1))"
    return target_name


def add_fit_line_with_ci(ax, x_values, y_values):
    if len(x_values) < 3:
        return

    if np.allclose(np.nanstd(x_values), 0):
        return

    coeffs = np.polyfit(x_values, y_values, deg=1)
    x_line = np.linspace(np.nanmin(x_values), np.nanmax(x_values), 200)
    y_line = np.polyval(coeffs, x_line)
    x_mean = np.mean(x_values)
    sxx = np.sum((x_values - x_mean) ** 2)
    if sxx <= 0:
        ax.plot(x_line, y_line, color="#d55e5e", linewidth=2.0)
        return

    residuals = y_values - np.polyval(coeffs, x_values)
    dof = len(x_values) - 2
    if dof <= 0:
        ax.plot(x_line, y_line, color="#d55e5e", linewidth=2.0)
        return

    residual_std = np.sqrt(np.sum(residuals ** 2) / dof)
    t_value = student_t.ppf(0.975, dof)
    se_line = residual_std * np.sqrt((1 / len(x_values)) + ((x_line - x_mean) ** 2 / sxx))
    lower = y_line - (t_value * se_line)
    upper = y_line + (t_value * se_line)

    ax.fill_between(x_line, lower, upper, color="#d55e5e", alpha=0.15, linewidth=0)
    ax.plot(x_line, y_line, color="#d55e5e", linewidth=2.0)


def plot_model_panel(model_name, model_result, output_dir):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)
    feature_types = ["CHH", "CHG"]
    target_name = model_result["target"]
    target_label = format_target_label(target_name)

    for idx, feature_type in enumerate(feature_types):
        ax = axes[idx]
        feature_key = feature_type.lower()
        correlation_key = f"{feature_key}_corr"
        feature_name = model_result.get(feature_key)

        if not feature_name:
            ax.axis("off")
            ax.text(
                0.5,
                0.5,
                f"{model_name}\nNo {feature_type} feature\nin top 10",
                ha="center",
                va="center",
                fontsize=11,
            )
            continue

        corr_result = model_result.get(correlation_key)
        if not corr_result:
            ax.axis("off")
            ax.text(
                0.5,
                0.5,
                f"{model_name}\n{feature_name}\nCorrelation failed",
                ha="center",
                va="center",
                fontsize=11,
            )
            continue

        x_values = corr_result["data"][feature_name].to_numpy(dtype=float)
        y_values = corr_result["data"][target_name].to_numpy(dtype=float)
        ax.scatter(
            x_values,
            y_values,
            s=34,
            alpha=0.78,
            color="#3a789c",
        )
        add_fit_line_with_ci(ax, x_values, y_values)
        ax.set_xlabel(feature_name)
        ax.set_ylabel(target_label)

        stats_text = (
            f"n = {corr_result['n_samples']}\n"
            f"Spearman rho = {corr_result['spearman_rho']:.4f}\n"
            f"Spearman p = {format_p_value(corr_result['spearman_p'])}"
        )
        ax.text(
            0.03,
            0.97,
            stats_text,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=11,
            bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.7, "edgecolor": "#a8a8a8"},
        )
        ax.grid(True, color="#bdbdbd", alpha=0.7, linestyle="-", linewidth=1.0)

    png_path = output_dir / f"{PANEL_BASENAME.format(model=model_name)}.png"
    pdf_path = output_dir / f"{PANEL_BASENAME.format(model=model_name)}.pdf"
    fig.savefig(png_path, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    return png_path, pdf_path


def save_summary(model_name, model_result, output_dir):
    rows = []
    for feature_type in ("CHH", "CHG"):
        feature_key = feature_type.lower()
        corr_key = f"{feature_key}_corr"
        feature_name = model_result.get(feature_key)
        corr_result = model_result.get(corr_key)
        status = "ok" if feature_name and corr_result else "skipped"

        rows.append(
            {
                "model": model_name,
                "target": model_result["target"],
                "target_label": format_target_label(model_result["target"]),
                "feature_type": feature_type,
                "feature_name": feature_name if feature_name else "",
                "status": status,
                "top10_features": ",".join(model_result["top10_features"]),
                "n_samples": corr_result["n_samples"] if corr_result else np.nan,
                "n_removed_missing": corr_result["n_removed"] if corr_result else np.nan,
                "spearman_rho": corr_result["spearman_rho"] if corr_result else np.nan,
                "spearman_p": corr_result["spearman_p"] if corr_result else np.nan,
            }
        )

    summary_df = pd.DataFrame(rows)
    summary_path = output_dir / SUMMARY_BASENAME.format(model=model_name)
    summary_df.to_csv(summary_path, sep="\t", index=False)
    return summary_path


def main():
    args = parse_args()
    data_file = Path(args.data_file)
    import_dir = Path(args.import_dir)
    output_root = Path(args.output_root)

    output_root.mkdir(parents=True, exist_ok=True)

    dataframe = pd.read_csv(data_file, sep="\t", low_memory=False)

    for model_name, import_filename, target_name in MODEL_CONFIGS:
        import_file = import_dir / import_filename
        model_output_dir = output_root / model_name
        model_output_dir.mkdir(parents=True, exist_ok=True)

        feature_pairs = parse_import_results(import_file)
        top10_features, chh_feature, chg_feature = select_top_context_features(feature_pairs, top_n=10)

        if chh_feature is None:
            print(f"Warning: {model_name} 前10特征中未找到 CHH 特征。")
        if chg_feature is None:
            print(f"Warning: {model_name} 前10特征中未找到 CHG 特征。")

        required_columns = [target_name]
        if chh_feature:
            required_columns.append(chh_feature)
        if chg_feature:
            required_columns.append(chg_feature)

        missing_columns = [column for column in required_columns if column not in dataframe.columns]
        if missing_columns:
            raise KeyError(f"{model_name} 缺少数据列: {', '.join(missing_columns)}")

        model_result = {
            "target": target_name,
            "top10_features": [name for name, _ in top10_features],
            "chh": chh_feature,
            "chg": chg_feature,
            "chh_corr": None,
            "chg_corr": None,
        }

        if chh_feature:
            model_result["chh_corr"] = compute_correlations(dataframe, chh_feature, target_name)
        if chg_feature:
            model_result["chg_corr"] = compute_correlations(dataframe, chg_feature, target_name)

        summary_path = save_summary(model_name, model_result, model_output_dir)
        png_path, pdf_path = plot_model_panel(model_name, model_result, model_output_dir)

        print(f"{model_name}: target={target_name}")
        print(f"  top10={', '.join(model_result['top10_features'])}")
        print(f"  CHH={chh_feature if chh_feature else 'None'}")
        print(f"  CHG={chg_feature if chg_feature else 'None'}")
        print(f"  summary={summary_path}")
        print(f"  figure_png={png_path}")
        print(f"  figure_pdf={pdf_path}")


if __name__ == "__main__":
    main()
