import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cache")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xgboost as xgb
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from sklearn import metrics
from sklearn.model_selection import train_test_split

try:
    import shap
except ImportError as exc:
    raise ImportError("The 'shap' package is required to run this import plotting script.") from exc

BASE_FONT_SIZE = 16
LABEL_FONT_SIZE = 18
TITLE_FONT_SIZE = 20
TICK_FONT_SIZE = 15
LEGEND_FONT_SIZE = 15

plt.rcParams.update({
    "font.size": BASE_FONT_SIZE,
    "axes.titlesize": TITLE_FONT_SIZE,
    "axes.labelsize": LABEL_FONT_SIZE,
    "xtick.labelsize": TICK_FONT_SIZE,
    "ytick.labelsize": TICK_FONT_SIZE,
    "legend.fontsize": LEGEND_FONT_SIZE,
})

NAME_TPM = "B_TPM"
TARGET_COLUMN = "β"
MODEL_VERSION = "1.0"
OUTPUT_PREFIX = "TPM_1.5"
FILTER_VALUES = [10, 16, 22]
THRESHOLD_VALUE = 0
TEST_SIZE = 0.2
RANDOM_STATE = 42

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parents[1]
LIST_DIR = PROJECT_DIR / "list"
DATA_PATH = LIST_DIR / "RCA" / "RCA.climate.tsv"


def beeswarm_offsets(values, nbins=60, spread=0.8):
    values = np.asarray(values)
    if len(values) <= 1:
        return np.zeros(len(values))

    bins = np.linspace(values.min(), values.max(), nbins + 1)
    bin_ids = np.digitize(values, bins) - 1
    offsets = np.zeros(len(values), dtype=float)

    for bin_id in np.unique(bin_ids):
        idx = np.where(bin_ids == bin_id)[0]
        if len(idx) <= 1:
            continue

        sorted_idx = idx[np.argsort(values[idx], kind="mergesort")]
        pattern = np.arange(len(sorted_idx), dtype=float)
        pattern = np.where(pattern % 2 == 0, pattern / 2, -(pattern + 1) / 2)
        max_abs = np.max(np.abs(pattern))
        if max_abs > 0:
            pattern = pattern / max_abs * spread
        offsets[sorted_idx] = pattern

    return offsets


def get_filter_column(df):
    for candidate in ("Temperature", "tem"):
        if candidate in df.columns:
            return candidate
    return df.columns[0]


def prepare_group_data(df_all, filter_column, filter_val):
    df = df_all[df_all[filter_column] == filter_val].copy()
    if df.empty:
        return None, None

    drop_columns = [filter_column, "total", "α", "β", "β1", "β2", "α/%", "β/%", "β1/%", "β2/%"]
    x = df.drop(columns=[col for col in drop_columns if col in df.columns])
    x = x.apply(pd.to_numeric, errors="coerce")
    y = df[TARGET_COLUMN]

    _, x_test, _, y_test = train_test_split(
        x,
        y,
        test_size=TEST_SIZE,
        random_state=RANDOM_STATE,
    )
    return x_test, y_test


def align_features(model, x_test):
    model_feature_names = model.get_booster().feature_names
    if not model_feature_names:
        return x_test

    missing = [feature for feature in model_feature_names if feature not in x_test.columns]
    if missing:
        raise ValueError(f"Missing features required by model: {missing}")

    return x_test.loc[:, model_feature_names]


def write_results(
    result_path,
    model_path,
    filter_val,
    mse,
    rmse,
    mae,
    r2,
    top_10_features,
    first_feature_contribution,
    total_contribution,
    first_feature_percentage,
    all_features,
    all_feature_contributions,
    top_80_features,
    top_80_contribution,
):
    with open(result_path, "w", encoding="utf-8") as f:
        f.write(f"=== {NAME_TPM} 导入模型结果 (Group {filter_val}) ===\n\n")
        f.write(f"数据文件: {DATA_PATH}\n")
        f.write(f"模型文件: {model_path}\n")
        f.write(f"测试集划分: test_size={TEST_SIZE}, random_state={RANDOM_STATE}\n\n")

        f.write("测试集预测结果:\n")
        f.write(f"  均方误差 (MSE): {mse:.6f}\n")
        f.write(f"  均方根误差 (RMSE): {rmse:.6f}\n")
        f.write(f"  平均绝对误差 (MAE): {mae:.6f}\n")
        f.write(f"  拟合优度 (R-squared): {r2:.6f}\n")

        f.write("\nSHAP阈值特征分析:\n")
        f.write("  说明: 基于测试集 SHAP\n")
        f.write(f"  所有特征的 SHAP 贡献度总和: {total_contribution:.6f}\n")
        f.write("  所有特征及其平均|SHAP|值:\n")
        for feature, value in zip(all_features, all_feature_contributions):
            f.write(f"    {feature}: {value:.6f}\n")

        f.write("\nSHAP特征重要性分析:\n")
        f.write(f"  前10个最重要特征: {', '.join(top_10_features)}\n")
        f.write(f"  第一个最重要特征: {top_10_features[0]}\n")
        f.write(f"  第一个特征的平均|SHAP|值: {first_feature_contribution:.6f}\n")
        f.write(f"  所有特征的平均|SHAP|值总和: {total_contribution:.6f}\n")
        f.write(f"  第一个特征贡献度百分比: {first_feature_percentage:.2f}%\n")
        f.write("\n  累计贡献度达到80%的特征:\n")
        f.write(f"    特征数量: {len(top_80_features)}\n")
        f.write(f"    累计贡献度: {top_80_contribution:.2f}%\n")
        f.write(f"    特征列表: {', '.join(top_80_features)}\n")


def write_chh_results(
    chh_result_path,
    model_path,
    filter_val,
    chh_features,
    chh_total_contribution,
    chh_percentage,
    total_contribution,
):
    with open(chh_result_path, "w", encoding="utf-8") as f:
        f.write(f"=== {NAME_TPM} CHH SHAP 汇总 (Group {filter_val}) ===\n\n")
        f.write(f"数据文件: {DATA_PATH}\n")
        f.write(f"模型文件: {model_path}\n\n")
        f.write(f"所有特征 mean(|SHAP|) 总和: {total_contribution:.6f}\n")
        f.write(f"以 CHH 开头的特征数量: {len(chh_features)}\n")
        f.write(f"以 CHH 开头的特征 mean(|SHAP|) 总和: {chh_total_contribution:.6f}\n")
        f.write(f"以 CHH 开头的特征占全部 mean(|SHAP|) 的百分比: {chh_percentage:.2f}%\n")
        f.write(f"CHH 特征列表: {', '.join(chh_features) if chh_features else '无'}\n")


def plot_combined_shap(output_dir, x_test, shap_values_numpy, sorted_indices, mean_abs_shap):
    max_display = min(15, x_test.shape[1])
    plot_indices = sorted_indices[:max_display]
    plot_features = x_test.columns[plot_indices].tolist()
    plot_mean_abs_shap = mean_abs_shap[plot_indices]

    threshold_plot_count = int(np.sum(plot_mean_abs_shap > THRESHOLD_VALUE))
    threshold_line_y = threshold_plot_count - 0.5 if 0 < threshold_plot_count < max_display else None

    cmap = plt.get_cmap("cool")
    fig = plt.figure(figsize=(12, 6), dpi=1200)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.9, 0.05], wspace=0.15)
    ax_bar = fig.add_subplot(gs[0, 0])
    ax_bee = fig.add_subplot(gs[0, 1], sharey=ax_bar)
    cax = fig.add_subplot(gs[0, 2])
    y_positions = np.arange(max_display)

    ax_bar.barh(y_positions, plot_mean_abs_shap, color="#2f83d0")
    ax_bar.set_yticks(y_positions)
    ax_bar.set_yticklabels(plot_features, fontsize=TICK_FONT_SIZE)
    ax_bar.invert_yaxis()
    ax_bar.set_xlabel("Mean(|SHAP value|)", fontsize=LABEL_FONT_SIZE)
    ax_bar.grid(axis="y", linestyle=":", linewidth=0.8, alpha=0.4)
    ax_bar.spines["top"].set_visible(False)
    ax_bar.spines["right"].set_visible(False)
    ax_bar.tick_params(axis="x", labelsize=TICK_FONT_SIZE)
    ax_bar.tick_params(axis="y", length=0, labelsize=TICK_FONT_SIZE)

    for row_idx, feature_idx in enumerate(plot_indices):
        feature_name = x_test.columns[feature_idx]
        feature_values = x_test[feature_name].to_numpy()
        shap_row_values = shap_values_numpy[:, feature_idx]
        offsets = beeswarm_offsets(shap_row_values)

        vmin = np.nanpercentile(feature_values, 5)
        vmax = np.nanpercentile(feature_values, 95)
        if np.isclose(vmin, vmax):
            vmin = np.nanmin(feature_values)
            vmax = np.nanmax(feature_values)
        if np.isclose(vmin, vmax):
            vmin -= 1
            vmax += 1
        norm = Normalize(vmin=vmin, vmax=vmax)

        ax_bee.scatter(
            shap_row_values,
            np.full_like(shap_row_values, row_idx, dtype=float) + offsets * 0.32,
            c=feature_values,
            cmap=cmap,
            norm=norm,
            s=16,
            alpha=0.95,
            linewidths=0,
        )

    ax_bee.axvline(0, color="gray", linewidth=1)
    ax_bee.set_yticks(y_positions)
    ax_bee.tick_params(axis="x", labelsize=TICK_FONT_SIZE)
    ax_bee.tick_params(axis="y", left=False, labelleft=False, labelsize=TICK_FONT_SIZE)
    ax_bee.set_xlabel("SHAP value (impact on model output)", fontsize=LABEL_FONT_SIZE)
    ax_bee.grid(axis="y", linestyle=":", linewidth=0.8, alpha=0.4)
    ax_bee.spines["top"].set_visible(False)
    ax_bee.spines["right"].set_visible(False)
    ax_bee.spines["left"].set_visible(False)

    if threshold_line_y is not None:
        ax_bar.axhline(y=threshold_line_y, color="black", linewidth=1.1)
        ax_bee.axhline(y=threshold_line_y, color="black", linewidth=1.1)

    colorbar = fig.colorbar(ScalarMappable(norm=Normalize(vmin=0, vmax=1), cmap=cmap), cax=cax)
    colorbar.set_ticks([0, 1])
    colorbar.set_ticklabels(["Low", "High"])
    colorbar.ax.tick_params(labelsize=TICK_FONT_SIZE)
    colorbar.set_label("Feature value", rotation=270, labelpad=20, fontsize=LABEL_FONT_SIZE)
    colorbar.outline.set_visible(False)

    fig.savefig(
        output_dir / f"SHAP_combined_{NAME_TPM}_{MODEL_VERSION}.pdf",
        format="pdf",
        bbox_inches="tight",
    )
    plt.close(fig)


def plot_dependence(output_dir, top_10_features, shap_values_explanation, x_test):
    for feature in top_10_features:
        shap.dependence_plot(
            feature,
            shap_values_explanation.values,
            x_test,
            interaction_index="auto",
            x_jitter=0,
            dot_size=12,
            alpha=0.9,
            show=False,
        )
        ax = plt.gca()
        feature_min = x_test[feature].min()
        feature_max = x_test[feature].max()
        if feature_min >= 0:
            right_margin = (feature_max - feature_min) * 0.03 if feature_max > feature_min else 0.1
            ax.set_xlim(left=0, right=feature_max + right_margin)
        ax.set_xlabel(ax.get_xlabel(), fontsize=LABEL_FONT_SIZE)
        ax.set_ylabel(ax.get_ylabel(), fontsize=LABEL_FONT_SIZE)
        ax.tick_params(axis="both", labelsize=TICK_FONT_SIZE)
        plt.axhline(y=0, color="black", linestyle="-.", linewidth=1)

        safe_feature_name = feature.replace("/", "_").replace("\\", "_")
        plt.savefig(
            output_dir / f"SHAP_Dependence_{safe_feature_name}_{NAME_TPM}_{MODEL_VERSION}.pdf",
            format="pdf",
            bbox_inches="tight",
            dpi=1200,
        )
        plt.close()


def plot_prediction_scatter(output_dir, filter_val, y_test, y_pred, r2, mae):
    plt.figure(figsize=(8, 6), dpi=1200)
    plt.scatter(y_test, y_pred, color="#f4ba8a", label="Predicted", alpha=0.2)
    max_value = max(max(y_test), max(y_pred))
    plt.plot([0, max_value], [0, max_value], color="black", alpha=0.6, linestyle="--", label="x=y")

    z = np.polyfit(y_test, y_pred, 1)
    p = np.poly1d(z)
    plt.plot(
        y_test,
        p(y_test),
        color="#b4d4e1",
        alpha=0.6,
        label=f"Line of Best Fit\n$R^2$ = {r2:.2f}, MAE = {mae:.2f}",
    )

    plt.title(rf'Expression level of Rca $\beta$ ({filter_val:g}$^\circ$C)', fontsize=TITLE_FONT_SIZE)
    plt.xlabel("Actual Values", fontsize=LABEL_FONT_SIZE)
    plt.ylabel("Predicted Values", fontsize=LABEL_FONT_SIZE)
    plt.xticks(fontsize=TICK_FONT_SIZE)
    plt.yticks(fontsize=TICK_FONT_SIZE)
    plt.legend(loc="upper left", fontsize=LEGEND_FONT_SIZE)
    plt.savefig(
        output_dir / f"{NAME_TPM}_{MODEL_VERSION}.pdf",
        format="pdf",
        bbox_inches="tight",
        dpi=1200,
    )
    plt.close()


def main():
    if not DATA_PATH.exists():
        raise FileNotFoundError(f"Input data not found: {DATA_PATH}")

    df_all = pd.read_csv(DATA_PATH, sep="\t")
    cols_to_log = ["total", "α", "β", "β1", "β2"]
    cols_to_transform = [col for col in cols_to_log if col in df_all.columns]
    df_all[cols_to_transform] = np.log1p(df_all[cols_to_transform])

    filter_column = get_filter_column(df_all)

    for filter_val in FILTER_VALUES:
        output_dir = LIST_DIR / "xgboot" / f"{OUTPUT_PREFIX}_{filter_val}"
        model_path = output_dir / f"my_model_{NAME_TPM}_{MODEL_VERSION}.json"
        result_path = output_dir / f"{NAME_TPM}_import_results.txt"
        chh_result_path = output_dir / f"{NAME_TPM}_CHH_SHAP_summary.txt"

        print(f"Processing group {filter_val}")
        print(f"Using data file: {DATA_PATH}")
        print(f"Using model file: {model_path}")

        if not model_path.exists():
            print(f"Skipping group {filter_val}: model not found.")
            continue

        output_dir.mkdir(parents=True, exist_ok=True)
        x_test, y_test = prepare_group_data(df_all, filter_column, filter_val)
        if x_test is None:
            print(f"Skipping group {filter_val}: no matching rows.")
            continue

        model = xgb.XGBRegressor()
        model.load_model(model_path)
        x_test = align_features(model, x_test)

        y_pred = model.predict(x_test)
        mse = metrics.mean_squared_error(y_test, y_pred)
        rmse = np.sqrt(mse)
        mae = metrics.mean_absolute_error(y_test, y_pred)
        r2 = metrics.r2_score(y_test, y_pred)

        print(f"{NAME_TPM} ({filter_val}) 导入模型预测结果:")
        print("均方误差 (MSE):", mse)
        print("均方根误差 (RMSE):", rmse)
        print("平均绝对误差 (MAE):", mae)
        print("拟合优度 (R-squared):", r2)

        explainer = shap.TreeExplainer(model)
        shap_values_numpy = explainer.shap_values(x_test)
        shap_values_explanation = explainer(x_test)

        mean_abs_shap = np.abs(shap_values_numpy).mean(axis=0)
        top_10_indices = np.argsort(mean_abs_shap)[-10:][::-1]
        top_10_features = x_test.columns[top_10_indices].tolist()

        first_feature_contribution = mean_abs_shap[top_10_indices[0]]
        total_contribution = mean_abs_shap.sum()
        first_feature_percentage = (first_feature_contribution / total_contribution) * 100
        chh_mask = x_test.columns.str.startswith("CHH")
        chh_features = x_test.columns[chh_mask].tolist()
        chh_total_contribution = mean_abs_shap[chh_mask].sum()
        chh_percentage = (chh_total_contribution / total_contribution * 100) if total_contribution else 0.0

        sorted_indices = np.argsort(mean_abs_shap)[::-1]
        all_features = x_test.columns[sorted_indices].tolist()
        all_feature_contributions = mean_abs_shap[sorted_indices]
        sorted_contributions = mean_abs_shap[sorted_indices]
        cumulative_contribution = np.cumsum(sorted_contributions) / total_contribution * 100

        threshold_80_idx = np.where(cumulative_contribution >= 80)[0][0]
        top_80_features = x_test.columns[sorted_indices[: threshold_80_idx + 1]].tolist()
        top_80_contribution = cumulative_contribution[threshold_80_idx]

        print(f"CHH 特征 mean(|SHAP|) 总和: {chh_total_contribution}")
        print(f"CHH 特征贡献度百分比: {chh_percentage:.2f}%")

        write_results(
            result_path,
            model_path,
            filter_val,
            mse,
            rmse,
            mae,
            r2,
            top_10_features,
            first_feature_contribution,
            total_contribution,
            first_feature_percentage,
            all_features,
            all_feature_contributions,
            top_80_features,
            top_80_contribution,
        )
        write_chh_results(
            chh_result_path,
            model_path,
            filter_val,
            chh_features,
            chh_total_contribution,
            chh_percentage,
            total_contribution,
        )
        plot_combined_shap(output_dir, x_test, shap_values_numpy, sorted_indices, mean_abs_shap)
        plot_dependence(output_dir, top_10_features, shap_values_explanation, x_test)
        plot_prediction_scatter(output_dir, filter_val, y_test, y_pred, r2, mae)

        print(f"Finished group {filter_val}. Results written to: {result_path}")
        print(f"CHH summary written to: {chh_result_path}")


if __name__ == "__main__":
    main()
