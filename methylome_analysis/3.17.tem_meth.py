#!/usr/bin/env python3
"""Plot temperature-group methylation/SNP differences for SHAP-selected features."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams
from matplotlib import font_manager
from scipy import stats
from statsmodels.stats.multitest import multipletests


PROJECT_MARKER = "Arabidopsis_thaliana"
INPUT_TSV_NAME = "Alt.all.snp_meth.filtered.tsv"
WEIGHTED_INPUT_TSV_NAME = "Alt.weighted.snp_meth.tsv"
IMPORT_RESULTS_DIRNAME = "TPM_4.5_all"
OUTPUT_SUBDIR = "temperature_meth_plots"
TEMPERATURE_COLUMN = "Temperature"
TEMPERATURE_ORDER = [22, 16, 10]
TEMPERATURE_LABELS = {22: "22°C", 16: "16°C", 10: "10°C"}
TEMPERATURE_COLORS = {22: "#d1495b", 16: "#edae49", 10: "#00798c"}
FEATURE_TYPE_ORDER = ["CHH", "CHG", "CG"]
DATASET_EXPRESSION_COLUMNS = {
    "A_per_TPM": "α/%",
    "A_TPM": "α",
    "B_per_TPM": "β/%",
    "B_TPM": "β",
    "total_TPM": "total",
}
LOG1P_TARGET_COLUMNS = {"total", "α", "β"}
WEIGHTED_ALIGNMENT_COLUMNS = ["Temperature", "total", "α", "β", "β1", "β2", "α/%", "β/%", "β1/%", "β2/%"]
PAIRWISE_TEMPERATURES = [(22, 16), (22, 10), (16, 10)]
FEATURE_PATTERN = re.compile(r"^\s*((?:CHH|CHG|CG|snp)_[^:]+):\s*([-+0-9.eE]+)\s*$")
DATASET_CONFIGS = [
    ("A_per_TPM_import_results.txt", "A_per_TPM"),
    ("A_TPM_import_results.txt", "A_TPM"),
    ("B_per_TPM_import_results.txt", "B_per_TPM"),
    ("B_TPM_import_results.txt", "B_TPM"),
    ("total_TPM_import_results.txt", "total_TPM"),
]
PERMUTATION_ITERATIONS = 4000
GLOBAL_ALPHA = 0.05
METHYLATION_COLORBAR_LABEL = "Methylation level"
SCATTER_TOP_N = 10


@dataclass
class DatasetResult:
    dataset_name: str
    features_by_type: Dict[str, List[str]]
    shap_values: Dict[str, float]


def configure_matplotlib() -> None:
    global METHYLATION_COLORBAR_LABEL
    cjk_candidates = [
        "Noto Sans CJK SC",
        "Noto Sans CJK JP",
        "WenQuanYi Zen Hei",
        "SimHei",
        "Microsoft YaHei",
        "Arial Unicode MS",
    ]
    available_fonts = {font.name for font in font_manager.fontManager.ttflist}
    cjk_font = next((name for name in cjk_candidates if name in available_fonts), None)
    if cjk_font is not None:
        rcParams["font.family"] = cjk_font
        METHYLATION_COLORBAR_LABEL = "甲基化水平"
    else:
        rcParams["font.family"] = "DejaVu Sans"
        METHYLATION_COLORBAR_LABEL = "Methylation level"
    rcParams["axes.unicode_minus"] = False
    rcParams["pdf.fonttype"] = 42
    rcParams["ps.fonttype"] = 42


def resolve_project_root() -> Path:
    script_path = Path(__file__).resolve()
    for parent in script_path.parents:
        if parent.name == PROJECT_MARKER:
            return parent
    raise RuntimeError(f"Unable to locate project root for {script_path}")


def candidate_paths(*paths: Path) -> Iterable[Path]:
    for path in paths:
        if path.exists():
            yield path


def resolve_paths() -> tuple[Path, Path, Path, Path]:
    project_root = resolve_project_root()
    local_tsv = project_root / "list" / "RCA" / INPUT_TSV_NAME
    local_weighted_tsv = project_root / "list" / "RCA" / WEIGHTED_INPUT_TSV_NAME
    local_import_dir = project_root / "list" / "xgboot" / IMPORT_RESULTS_DIRNAME
    external_tsv = Path("/datapool/life-gongl/zesheng/Arabidopsis_thaliana/list/RCA") / INPUT_TSV_NAME
    external_weighted_tsv = Path("/datapool/life-gongl/zesheng/Arabidopsis_thaliana/list/RCA") / WEIGHTED_INPUT_TSV_NAME
    external_import_dir = Path("/datapool/life-gongl/zesheng/Arabidopsis_thaliana/list/xgboot") / IMPORT_RESULTS_DIRNAME

    tsv_path = next(candidate_paths(external_tsv, local_tsv), None)
    weighted_tsv_path = next(candidate_paths(external_weighted_tsv, local_weighted_tsv), None)
    import_dir = next(candidate_paths(external_import_dir, local_import_dir), None)
    if tsv_path is None or weighted_tsv_path is None or import_dir is None:
        raise FileNotFoundError("Failed to resolve the methylation table, weighted methylation table, or import_results directory.")

    output_dir = project_root / "list" / "xgboot" / IMPORT_RESULTS_DIRNAME / OUTPUT_SUBDIR
    output_dir.mkdir(parents=True, exist_ok=True)
    return tsv_path, weighted_tsv_path, import_dir, output_dir


def parse_import_results(result_path: Path) -> DatasetResult:
    features_by_type = {feature_type: [] for feature_type in FEATURE_TYPE_ORDER}
    shap_values: Dict[str, float] = {}
    in_feature_block = False

    with result_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if "所有特征及其平均|SHAP|值" in line:
                in_feature_block = True
                continue
            if not in_feature_block:
                continue

            match = FEATURE_PATTERN.match(line)
            if not match:
                if line.strip():
                    continue
                break

            feature_name, shap_text = match.groups()
            shap_value = float(shap_text)
            if shap_value <= 0:
                continue
            feature_type = feature_name.split("_", 1)[0]
            if feature_type not in features_by_type:
                continue
            features_by_type[feature_type].append(feature_name)
            shap_values[feature_name] = shap_value

    total_features = sum(len(items) for items in features_by_type.values())
    if total_features == 0:
        raise ValueError(f"No SHAP-positive features were found in {result_path}")
    dataset_name = result_path.name.replace("_import_results.txt", "")
    return DatasetResult(dataset_name=dataset_name, features_by_type=features_by_type, shap_values=shap_values)


def select_top_features_for_scatter(dataset_result: DatasetResult, top_n: int = SCATTER_TOP_N) -> Dict[str, List[str]]:
    selected: Dict[str, List[str]] = {}
    for feature_type in FEATURE_TYPE_ORDER:
        ranked = sorted(
            dataset_result.features_by_type[feature_type],
            key=lambda feature: dataset_result.shap_values[feature],
            reverse=True,
        )
        selected[feature_type] = ranked[:top_n]
    return selected


def load_main_table(tsv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    if TEMPERATURE_COLUMN not in df.columns:
        raise KeyError(f"Column {TEMPERATURE_COLUMN!r} not found in {tsv_path}")

    df[TEMPERATURE_COLUMN] = pd.to_numeric(df[TEMPERATURE_COLUMN], errors="coerce")
    observed_temps = sorted(int(temp) for temp in df[TEMPERATURE_COLUMN].dropna().unique())
    if observed_temps != sorted(TEMPERATURE_ORDER):
        raise ValueError(f"Expected temperatures {sorted(TEMPERATURE_ORDER)}, observed {observed_temps}")

    feature_columns = [col for col in df.columns if col != TEMPERATURE_COLUMN]
    df[feature_columns] = df[feature_columns].apply(pd.to_numeric, errors="coerce")
    return df


def load_weighted_table(tsv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    if TEMPERATURE_COLUMN not in df.columns:
        raise KeyError(f"Column {TEMPERATURE_COLUMN!r} not found in {tsv_path}")

    df[TEMPERATURE_COLUMN] = pd.to_numeric(df[TEMPERATURE_COLUMN], errors="coerce")
    observed_temps = sorted(int(temp) for temp in df[TEMPERATURE_COLUMN].dropna().unique())
    if observed_temps != sorted(TEMPERATURE_ORDER):
        raise ValueError(f"Expected temperatures {sorted(TEMPERATURE_ORDER)}, observed {observed_temps}")

    feature_columns = [col for col in df.columns if col != TEMPERATURE_COLUMN]
    df[feature_columns] = df[feature_columns].apply(pd.to_numeric, errors="coerce")
    return df


def benjamini_hochberg(p_values: Sequence[float]) -> list[float]:
    if not p_values:
        return []
    adjusted = multipletests(list(p_values), method="fdr_bh")[1]
    return [float(value) for value in adjusted]


def significance_stars(p_value: float | None) -> str:
    if p_value is None or pd.isna(p_value):
        return "n.s."
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return "n.s."


def draw_significance_bracket(
    ax: plt.Axes,
    x1: float,
    x2: float,
    y: float,
    height: float,
    stars: str,
    text_offset: float,
    text_x_shift: float = 0.0,
) -> None:
    if not stars:
        return
    ax.plot([x1, x1, x2, x2], [y, y + height, y + height, y], color="black", linewidth=0.9, clip_on=False)
    ax.text(
        (x1 + x2) / 2 + text_x_shift,
        y + height + text_offset,
        stars,
        ha="center",
        va="bottom",
        fontsize=11,
        fontweight="bold",
    )


def dunn_pairwise(groups: Dict[int, pd.Series]) -> list[tuple[str, float, float]]:
    arrays = [groups[temp].dropna().to_numpy(dtype=float) for temp in TEMPERATURE_ORDER]
    if any(len(arr) == 0 for arr in arrays):
        return []

    all_values = np.concatenate(arrays)
    if len(np.unique(all_values)) <= 1:
        return []

    ranks = stats.rankdata(all_values)
    lengths = [len(arr) for arr in arrays]
    offsets = np.cumsum([0] + lengths)
    mean_ranks = [np.mean(ranks[offsets[i] : offsets[i + 1]]) for i in range(len(arrays))]

    n_total = len(all_values)
    tie_counts = np.unique(all_values, return_counts=True)[1]
    tie_correction = float(np.sum(tie_counts**3 - tie_counts) / (12 * (n_total - 1))) if n_total > 1 else 0.0
    variance = n_total * (n_total + 1) / 12 - tie_correction
    if variance <= 0:
        return []

    labels = []
    raw_p_values = []
    for temp_a, temp_b in PAIRWISE_TEMPERATURES:
        idx_a = TEMPERATURE_ORDER.index(temp_a)
        idx_b = TEMPERATURE_ORDER.index(temp_b)
        z_score = abs(mean_ranks[idx_a] - mean_ranks[idx_b]) / np.sqrt(
            variance * (1.0 / lengths[idx_a] + 1.0 / lengths[idx_b])
        )
        p_value = float(2 * stats.norm.sf(abs(z_score)))
        labels.append(f"{temp_a}_vs_{temp_b}")
        raw_p_values.append(p_value)

    adjusted = benjamini_hochberg(raw_p_values)
    return list(zip(labels, raw_p_values, adjusted))


def chi_square_permutation_test(values: pd.Series, temperatures: pd.Series, n_permutations: int) -> float:
    clean = pd.DataFrame({"value": values, "temp": temperatures}).dropna()
    if clean.empty:
        return math.nan

    clean["value"] = clean["value"].astype(int)
    observed = pd.crosstab(clean["temp"], clean["value"]).reindex(index=TEMPERATURE_ORDER, columns=[0, 1], fill_value=0)
    observed_stat = float(stats.chi2_contingency(observed.to_numpy(), correction=False)[0])

    temp_array = clean["temp"].to_numpy()
    value_array = clean["value"].to_numpy()
    rng = np.random.default_rng(42)
    exceedances = 0
    for _ in range(n_permutations):
        shuffled = rng.permutation(temp_array)
        perm_table = pd.crosstab(shuffled, value_array).reindex(index=TEMPERATURE_ORDER, columns=[0, 1], fill_value=0)
        perm_stat = float(stats.chi2_contingency(perm_table.to_numpy(), correction=False)[0])
        if perm_stat >= observed_stat - 1e-12:
            exceedances += 1
    return (exceedances + 1) / (n_permutations + 1)


def fisher_pairwise(groups: Dict[int, pd.Series]) -> list[tuple[str, float, float]]:
    labels = []
    raw_p_values = []
    for temp_a, temp_b in PAIRWISE_TEMPERATURES:
        values_a = groups[temp_a].dropna().astype(int)
        values_b = groups[temp_b].dropna().astype(int)
        table = np.array(
            [
                [(values_a == 0).sum(), (values_a == 1).sum()],
                [(values_b == 0).sum(), (values_b == 1).sum()],
            ],
            dtype=int,
        )
        _, p_value = stats.fisher_exact(table)
        labels.append(f"{temp_a}_vs_{temp_b}")
        raw_p_values.append(float(p_value))
    adjusted = benjamini_hochberg(raw_p_values)
    return list(zip(labels, raw_p_values, adjusted))


def snp_pairwise_tests(groups: Dict[int, pd.Series]) -> list[tuple[str, str, float, float]]:
    labels = []
    methods = []
    raw_p_values = []

    for temp_a, temp_b in PAIRWISE_TEMPERATURES:
        values_a = groups[temp_a].dropna().astype(int)
        values_b = groups[temp_b].dropna().astype(int)
        table = np.array(
            [
                [(values_a == 0).sum(), (values_a == 1).sum()],
                [(values_b == 0).sum(), (values_b == 1).sum()],
            ],
            dtype=int,
        )
        _, chi2_p, _dof, expected = stats.chi2_contingency(table, correction=False)
        if np.any(expected < 5):
            method = "Fisher"
            _, p_value = stats.fisher_exact(table)
        else:
            method = "Chi-square"
            p_value = float(chi2_p)
        labels.append(f"{temp_a}_vs_{temp_b}")
        methods.append(method)
        raw_p_values.append(float(p_value))

    adjusted = benjamini_hochberg(raw_p_values)
    return list(zip(labels, methods, raw_p_values, adjusted))


def analyze_continuous_feature(df: pd.DataFrame, feature: str, feature_type: str) -> tuple[dict, list[dict]]:
    groups = {
        temp: pd.to_numeric(
            df.loc[df[TEMPERATURE_COLUMN] == temp, feature],
            errors="coerce",
        ).dropna()
        for temp in TEMPERATURE_ORDER
    }
    if any(group.empty for group in groups.values()):
        overall_p = math.nan
        overall_method = "Kruskal-Wallis"
    else:
        _, overall_p = stats.kruskal(*(groups[temp] for temp in TEMPERATURE_ORDER), nan_policy="omit")
        overall_method = "Kruskal-Wallis"

    overall_row = {
        "feature": feature,
        "type": feature_type,
        "overall_test": overall_method,
        "overall_p": float(overall_p) if pd.notna(overall_p) else math.nan,
        "pairwise_comparison": "",
        "raw_p": math.nan,
        "bh_p": math.nan,
        "significance": significance_stars(float(overall_p)) if pd.notna(overall_p) else "",
    }

    pairwise_rows: list[dict] = []
    if pd.notna(overall_p) and float(overall_p) < GLOBAL_ALPHA:
        for label, raw_p, adj_p in dunn_pairwise(groups):
            pairwise_rows.append(
                {
                    "feature": feature,
                    "type": feature_type,
                    "overall_test": overall_method,
                    "overall_p": float(overall_p),
                    "pairwise_comparison": label,
                    "raw_p": raw_p,
                    "bh_p": adj_p,
                    "significance": significance_stars(adj_p),
                }
            )
    return overall_row, pairwise_rows


def analyze_snp_feature(df: pd.DataFrame, feature: str) -> tuple[dict, list[dict]]:
    values = pd.to_numeric(df[feature], errors="coerce")
    non_na = values.dropna()
    invalid = set(non_na.unique()) - {0, 1}
    if invalid:
        raise ValueError(f"SNP feature {feature} contains non-binary values: {sorted(invalid)}")

    clean = pd.DataFrame({TEMPERATURE_COLUMN: df[TEMPERATURE_COLUMN], feature: values}).dropna()
    overall_row = {
        "feature": feature,
        "type": "snp",
        "overall_test": "Not tested",
        "overall_p": math.nan,
        "pairwise_comparison": "",
        "raw_p": math.nan,
        "bh_p": math.nan,
        "overall_bh_p": math.nan,
        "significance": "",
    }

    groups = {
        temp: clean.loc[clean[TEMPERATURE_COLUMN] == temp, feature].astype(int)
        for temp in TEMPERATURE_ORDER
    }
    pairwise_rows: list[dict] = []
    for label, method, raw_p, adj_p in snp_pairwise_tests(groups):
        pairwise_rows.append(
            {
                "feature": feature,
                "type": "snp",
                "overall_test": method,
                "overall_p": math.nan,
                "pairwise_comparison": label,
                "raw_p": raw_p,
                "bh_p": adj_p,
                "overall_bh_p": math.nan,
                "significance": significance_stars(adj_p),
            }
        )
    return overall_row, pairwise_rows


def summarize_feature_statistics(df: pd.DataFrame, features_by_type: Dict[str, List[str]]) -> pd.DataFrame:
    rows: list[dict] = []
    for feature_type in FEATURE_TYPE_ORDER:
        for feature in features_by_type.get(feature_type, []):
            if feature_type == "snp":
                overall_row, pairwise_rows = analyze_snp_feature(df, feature)
            else:
                overall_row, pairwise_rows = analyze_continuous_feature(df, feature, feature_type)
            rows.append(overall_row)
            rows.extend(pairwise_rows)
    stats_df = pd.DataFrame(rows)
    stats_df["overall_bh_p"] = math.nan

    for feature_type in ("CHH", "CHG", "CG"):
        mask = (stats_df["type"] == feature_type) & (stats_df["pairwise_comparison"] == "") & stats_df["overall_p"].notna()
        if not mask.any():
            continue
        adjusted = benjamini_hochberg(stats_df.loc[mask, "overall_p"].tolist())
        stats_df.loc[mask, "overall_bh_p"] = adjusted
        stats_df.loc[mask, "significance"] = [significance_stars(value) for value in adjusted]
        mapping = stats_df.loc[mask, ["feature", "overall_bh_p"]].set_index("feature")["overall_bh_p"]
        pair_mask = (stats_df["type"] == feature_type) & (stats_df["pairwise_comparison"] != "")
        stats_df.loc[pair_mask, "overall_bh_p"] = stats_df.loc[pair_mask, "feature"].map(mapping)

    column_order = [
        "feature",
        "type",
        "overall_test",
        "overall_p",
        "overall_bh_p",
        "pairwise_comparison",
        "raw_p",
        "bh_p",
        "significance",
    ]
    return stats_df[column_order]


GROUP_ORDER = ["Low", "Medium", "High"]
GROUP_TO_TEMP = {"High": 22, "Medium": 16, "Low": 10}
TEMP_TO_GROUP = {value: key for key, value in GROUP_TO_TEMP.items()}
GROUP_DISPLAY_LABELS = {"High": "22℃", "Medium": "16℃", "Low": "10℃"}
GROUP_COLORS = {"High": "#d1495b", "Medium": "#edae49", "Low": "#00798c"}
BOX_OFFSETS = {"Low": -0.24, "Medium": 0.0, "High": 0.24}
SCATTER_REFERENCE_RATIO = 6267 / 4913


def add_group_labels(df: pd.DataFrame) -> pd.DataFrame:
    labeled = df.copy()
    if "Group" not in labeled.columns:
        labeled["Group"] = labeled[TEMPERATURE_COLUMN].map(TEMP_TO_GROUP)
    return labeled


def align_weighted_inputs(value_df: pd.DataFrame, weighted_df: pd.DataFrame) -> pd.DataFrame:
    merged = value_df.reset_index(names="value_index").merge(
        weighted_df.reset_index(names="weighted_index"),
        on=WEIGHTED_ALIGNMENT_COLUMNS,
        how="left",
        suffixes=("_value", "_weighted"),
    )
    if len(merged) != len(value_df):
        raise ValueError("Weighted methylation table alignment produced duplicate or missing rows.")
    if merged["weighted_index"].isna().any():
        raise ValueError("Weighted methylation table alignment failed for some methylation rows.")
    if (merged.groupby("value_index").size() > 1).any():
        raise ValueError("Weighted methylation table alignment is ambiguous for some methylation rows.")
    aligned_indices = merged.sort_values("value_index")["weighted_index"].astype(int).to_list()
    return weighted_df.iloc[aligned_indices].reset_index(drop=True).copy()


def prepare_dataset_frames(
    value_df: pd.DataFrame,
    weighted_df: pd.DataFrame,
    dataset_name: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr_col = DATASET_EXPRESSION_COLUMNS[dataset_name]
    expr_series = pd.to_numeric(value_df[expr_col], errors="coerce")
    if expr_col in LOG1P_TARGET_COLUMNS:
        eligible_mask = expr_series.notna() & (expr_series > 0)
        transformed = np.log1p(expr_series.loc[eligible_mask])
    else:
        eligible_mask = expr_series.notna()
        transformed = expr_series.loc[eligible_mask]

    if transformed.empty:
        raise ValueError(f"No eligible samples found for dataset {dataset_name} using {expr_col}")

    q1 = float(transformed.quantile(1 / 3))
    q3 = float(transformed.quantile(2 / 3))
    selected_value_df = value_df.loc[eligible_mask].copy()
    selected_weighted_df = weighted_df.loc[eligible_mask].copy()

    def classify(value: float) -> str:
        transformed_value = float(np.log1p(value)) if expr_col in LOG1P_TARGET_COLUMNS else float(value)
        if transformed_value < q1:
            return "Low"
        if transformed_value > q3:
            return "High"
        return "Medium"

    selected_value_df["Group"] = expr_series.loc[eligible_mask].map(classify)
    selected_weighted_df["Group"] = selected_value_df["Group"].to_numpy()
    return selected_value_df, selected_weighted_df


def weighted_average(values: np.ndarray, coverages: np.ndarray) -> np.ndarray:
    valid_mask = (~np.isnan(values)) & (~np.isnan(coverages)) & (coverages > 0)
    safe_values = np.where(valid_mask, values, 0.0)
    safe_coverages = np.where(valid_mask, coverages, 0.0)
    weighted_sum = np.sum(safe_values * safe_coverages, axis=1)
    weight_sum = np.sum(safe_coverages, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        result = weighted_sum / weight_sum
    result[weight_sum == 0] = np.nan
    return result


def build_weighted_type_table(value_df: pd.DataFrame, coverage_df: pd.DataFrame, dataset_result: DatasetResult) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for feature_type in FEATURE_TYPE_ORDER:
        features = dataset_result.features_by_type[feature_type]
        if not features:
            continue
        values = value_df[features].to_numpy(dtype=float)
        coverages = coverage_df[features].to_numpy(dtype=float)
        weighted_values = weighted_average(values, coverages)
        feature_table = pd.DataFrame(
            {
                "sample_index": value_df.index,
                "temperature": value_df[TEMPERATURE_COLUMN],
                "group": value_df["Group"],
                "feature_type": feature_type,
                "feature_value": weighted_values,
            }
        )
        rows.append(feature_table)
    if not rows:
        return pd.DataFrame(columns=["sample_index", "temperature", "group", "feature_type", "feature_value"])
    return pd.concat(rows, ignore_index=True)


def analyze_weighted_continuous(weighted_df: pd.DataFrame, feature_type: str) -> tuple[dict, list[dict]]:
    subset = weighted_df.loc[weighted_df["feature_type"] == feature_type].copy()
    groups = {
        group: pd.to_numeric(subset.loc[subset["group"] == group, "feature_value"], errors="coerce").dropna()
        for group in GROUP_ORDER
    }
    if any(group.empty for group in groups.values()):
        overall_p = math.nan
    else:
        _, overall_p = stats.kruskal(*(groups[group] for group in GROUP_ORDER), nan_policy="omit")

    overall_row = {
        "feature_type": feature_type,
        "overall_test": "Kruskal-Wallis",
        "overall_p": float(overall_p) if pd.notna(overall_p) else math.nan,
        "overall_bh_p": math.nan,
        "pairwise_comparison": "",
        "pairwise_test": "",
        "raw_p": math.nan,
        "bh_p": math.nan,
        "significance": "",
    }
    pairwise_rows: list[dict] = []
    if pd.notna(overall_p) and float(overall_p) < GLOBAL_ALPHA:
        temp_groups = {GROUP_TO_TEMP[group]: groups[group] for group in GROUP_ORDER}
        for label, raw_p, adj_p in dunn_pairwise(temp_groups):
            pairwise_rows.append(
                {
                    "feature_type": feature_type,
                    "overall_test": "Kruskal-Wallis",
                    "overall_p": float(overall_p),
                    "overall_bh_p": math.nan,
                    "pairwise_comparison": (
                        label.replace("22", GROUP_DISPLAY_LABELS["High"])
                        .replace("16", GROUP_DISPLAY_LABELS["Medium"])
                        .replace("10", GROUP_DISPLAY_LABELS["Low"])
                    ),
                    "pairwise_test": "Dunn",
                    "raw_p": raw_p,
                    "bh_p": adj_p,
                    "significance": significance_stars(adj_p),
                }
            )
    return overall_row, pairwise_rows


def summarize_weighted_type_statistics(weighted_df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for feature_type in ("CHH", "CHG", "CG"):
        overall_row, pairwise_rows = analyze_weighted_continuous(weighted_df, feature_type)
        rows.append(overall_row)
        rows.extend(pairwise_rows)

    stats_df = pd.DataFrame(rows)
    cont_mask = (stats_df["feature_type"].isin(["CHH", "CHG", "CG"])) & (stats_df["pairwise_comparison"] == "")
    adjusted = benjamini_hochberg(stats_df.loc[cont_mask, "overall_p"].fillna(1.0).tolist())
    stats_df.loc[cont_mask, "overall_bh_p"] = adjusted
    stats_df.loc[cont_mask, "significance"] = [significance_stars(value) for value in adjusted]
    mapping = stats_df.loc[cont_mask, ["feature_type", "overall_bh_p"]].set_index("feature_type")["overall_bh_p"]
    pair_mask = (stats_df["feature_type"].isin(["CHH", "CHG", "CG"])) & (stats_df["pairwise_comparison"] != "")
    stats_df.loc[pair_mask, "overall_bh_p"] = stats_df.loc[pair_mask, "feature_type"].map(mapping)
    return stats_df


def build_scatter_panel(
    ax: plt.Axes,
    df: pd.DataFrame,
    dataset_result: DatasetResult,
    feature_type: str,
    scatter_features: Dict[str, List[str]],
) -> None:
    features = scatter_features[feature_type]
    display_name = "SNP" if feature_type == "snp" else feature_type
    ax.set_title(f"{display_name} features ({len(features)})", loc="left", fontsize=20, fontweight="bold")
    if not features:
        ax.text(0.5, 0.5, "No features", transform=ax.transAxes, ha="center", va="center", fontsize=22)
        ax.set_xticks([])
        return

    labeled_df = add_group_labels(df)
    norm_max = 100.0
    cmap = "viridis"
    band_height = norm_max
    band_gap = 15.0
    band_base = {
        "Low": 0.0,
        "Medium": band_height + band_gap,
        "High": 2 * (band_height + band_gap),
    }
    median_lines: dict[str, list[float]] = {group: [] for group in GROUP_ORDER}

    for group in GROUP_ORDER:
        subset = labeled_df[labeled_df["Group"] == group]
        x_values: list[float] = []
        y_values: list[float] = []
        c_values: list[float] = []
        for idx, feature in enumerate(features, start=1):
            values = pd.to_numeric(subset[feature], errors="coerce").dropna().to_numpy(dtype=float)
            jitter = np.linspace(-0.12, 0.12, num=len(values)) if len(values) > 1 else np.zeros(len(values))
            x_values.extend((idx + jitter).tolist())
            y_values.extend((band_base[group] + values).tolist())
            c_values.extend(values.tolist())
            median_lines[group].append(float(np.nanmedian(values)) if len(values) else np.nan)
        if x_values:
            ax.scatter(
                np.asarray(x_values, dtype=float),
                np.asarray(y_values, dtype=float),
                c=np.asarray(c_values, dtype=float),
                cmap=cmap,
                vmin=0,
                vmax=norm_max,
                s=18,
                alpha=0.45,
                edgecolors="none",
            )
        valid_points = [
            (idx, band_base[group] + value)
            for idx, value in enumerate(median_lines[group], start=1)
            if not np.isnan(value)
        ]
        if valid_points:
            median_x, median_y = zip(*valid_points)
            ax.plot(median_x, median_y, color=GROUP_COLORS[group], linewidth=2.4, alpha=0.95)

    ax.set_xlim(0.65, len(features) + 0.35)
    ax.set_ylim(-2.0, band_base["High"] + band_height + 2.0)
    tick_values = []
    tick_labels = []
    inner_ticks = [0, 50, 100]
    for group in ["Low", "Medium", "High"]:
        for tick in inner_ticks:
            tick_values.append(band_base[group] + tick)
            tick_labels.append(str(tick))
    ax.axhline(band_base["Medium"], color="#bbbbbb", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.axhline(band_base["High"], color="#bbbbbb", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.set_yticks(tick_values)
    ax.set_yticklabels(tick_labels, fontsize=17)
    ax.set_ylabel("Methylation level", fontsize=24, labelpad=4)
    ax.set_xticks(np.arange(1, len(features) + 1))
    rotation = 0
    x_fontsize = 14
    ax.set_xticklabels(features, rotation=rotation, ha="center", fontsize=x_fontsize)
    ax.set_xlabel("features", fontsize=17)
    ax.tick_params(axis="y", labelsize=17)
    ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.35)
    ax.text(1.008, 0.84, "22℃", transform=ax.transAxes, color=GROUP_COLORS["High"], fontsize=22, va="center")
    ax.text(1.008, 0.50, "16℃", transform=ax.transAxes, color=GROUP_COLORS["Medium"], fontsize=22, va="center")
    ax.text(1.008, 0.16, "10℃", transform=ax.transAxes, color=GROUP_COLORS["Low"], fontsize=22, va="center")


def plot_feature_scatter_figure(df: pd.DataFrame, dataset_result: DatasetResult, output_dir: Path) -> tuple[Path, Path]:
    scatter_features = select_top_features_for_scatter(dataset_result)
    max_features = max(len(features) for features in scatter_features.values()) if scatter_features else 1
    fig_width = max(18, min(24, 12 + max_features * 0.9))
    fig_height = fig_width / SCATTER_REFERENCE_RATIO
    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(fig_width, fig_height), constrained_layout=False)
    fig.subplots_adjust(left=0.07, right=0.985, top=0.985, bottom=0.08, hspace=0.34)

    for axis, feature_type in zip(np.atleast_1d(axes), FEATURE_TYPE_ORDER):
        build_scatter_panel(axis, df, dataset_result, feature_type, scatter_features)

    pdf_path = output_dir / f"{dataset_result.dataset_name}_feature_scatter.pdf"
    png_path = output_dir / f"{dataset_result.dataset_name}_feature_scatter.png"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return pdf_path, png_path


def build_weighted_boxplot_panel(
    ax: plt.Axes,
    weighted_df: pd.DataFrame,
    stats_df: pd.DataFrame,
    feature_types: Sequence[str],
    y_max: float,
    y_ticks: Sequence[float],
) -> None:
    positions = np.arange(len(feature_types), dtype=float)
    type_offsets = BOX_OFFSETS
    bracket_layout = {
        ("High", "Medium"): {"level": 0, "text_shift": -0.03},
        ("Medium", "Low"): {"level": 1, "text_shift": 0.03},
        ("High", "Low"): {"level": 2, "text_shift": 0.0},
    }

    for group in GROUP_ORDER:
        data = []
        for feature_type in feature_types:
            subset = weighted_df.loc[
                (weighted_df["feature_type"] == feature_type) & (weighted_df["group"] == group),
                "feature_value",
            ].dropna()
            data.append(subset.to_numpy(dtype=float))
        bp = ax.boxplot(
            data,
            positions=positions + type_offsets[group],
            widths=0.16,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "black", "linewidth": 1.0},
            whiskerprops={"linewidth": 1.0},
            capprops={"linewidth": 1.0},
            boxprops={"linewidth": 1.0},
        )
        for patch in bp["boxes"]:
            patch.set_facecolor(GROUP_COLORS[group])
            patch.set_alpha(0.8)

    for idx, feature_type in enumerate(feature_types):
        feature_stats = stats_df.loc[stats_df["feature_type"] == feature_type]
        pairwise_rows = feature_stats.loc[feature_stats["pairwise_comparison"] != ""]
        overall_row = feature_stats.loc[feature_stats["pairwise_comparison"] == ""]
        pairwise_map = {row["pairwise_comparison"]: row["significance"] for row in pairwise_rows.to_dict("records")}

        local_ymax = 45.0
        annotation_step = 4.0
        bracket_height = 1.4
        text_offset = 0.8

        for first, second in [("High", "Medium"), ("Medium", "Low"), ("High", "Low")]:
            label = f"{GROUP_DISPLAY_LABELS[first]}_vs_{GROUP_DISPLAY_LABELS[second]}"
            stars = pairwise_map.get(label, "")
            if not stars:
                continue
            layout = bracket_layout[(first, second)]
            y = local_ymax + annotation_step * (0.75 + layout["level"])
            draw_significance_bracket(
                ax,
                idx + type_offsets[first],
                idx + type_offsets[second],
                y,
                bracket_height,
                stars,
                text_offset,
                layout["text_shift"],
            )

        if not pairwise_map and not overall_row.empty and overall_row["significance"].iloc[0]:
            overall_star = overall_row["significance"].iloc[0]
            y = local_ymax + annotation_step * 1.2
            draw_significance_bracket(
                ax,
                idx + type_offsets["High"],
                idx + type_offsets["Low"],
                y,
                bracket_height,
                overall_star,
                text_offset,
            )

    ax.set_xticks(positions)
    ax.set_xticklabels(feature_types, fontsize=18)
    ax.set_ylabel("Methylation level", fontsize=24)
    ax.set_ylim(0, y_max)
    ax.set_yticks(list(y_ticks))
    ax.tick_params(axis="y", labelsize=16)
    ax.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.5)
    ax.legend(
        handles=[
            plt.Line2D(
                [0],
                [0],
                color=GROUP_COLORS[group],
                marker="s",
                linestyle="",
                markersize=11,
                label=GROUP_DISPLAY_LABELS[group],
            )
            for group in GROUP_ORDER
        ],
        title="Temperature",
        fontsize=16,
        title_fontsize=18,
        loc="upper right",
        frameon=False,
    )


def plot_single_weighted_boxplot_figure(
    dataset_result: DatasetResult,
    weighted_df: pd.DataFrame,
    stats_df: pd.DataFrame,
    output_dir: Path,
    feature_type: str,
) -> tuple[Path, Path]:
    fig_main, ax_main = plt.subplots(figsize=(10, 8), constrained_layout=False)
    fig_main.subplots_adjust(left=0.11, right=0.96, top=0.95, bottom=0.16)
    build_weighted_boxplot_panel(
        ax_main,
        weighted_df,
        stats_df,
        [feature_type],
        y_max=60,
        y_ticks=np.arange(0, 61, 10),
    )
    pdf_path = output_dir / f"{dataset_result.dataset_name}_{feature_type}_weighted_boxplot.pdf"
    png_path = output_dir / f"{dataset_result.dataset_name}_{feature_type}_weighted_boxplot.png"
    fig_main.savefig(pdf_path)
    fig_main.savefig(png_path, dpi=300)
    plt.close(fig_main)
    return pdf_path, png_path


def plot_weighted_boxplot_figure(
    dataset_result: DatasetResult,
    weighted_df: pd.DataFrame,
    stats_df: pd.DataFrame,
    output_dir: Path,
) -> tuple[list[tuple[str, Path, Path]], Path]:
    plot_paths: list[tuple[str, Path, Path]] = []
    for feature_type in ["CHH", "CHG", "CG"]:
        pdf_path, png_path = plot_single_weighted_boxplot_figure(
            dataset_result,
            weighted_df,
            stats_df,
            output_dir,
            feature_type,
        )
        plot_paths.append((feature_type, pdf_path, png_path))
    stats_path = output_dir / f"{dataset_result.dataset_name}_weighted_boxplot_stats.tsv"
    stats_df.to_csv(stats_path, sep="\t", index=False)
    return plot_paths, stats_path


def validate_feature_columns(value_df: pd.DataFrame, weighted_df: pd.DataFrame, dataset_results: Sequence[DatasetResult]) -> None:
    available_value_columns = set(value_df.columns)
    available_weighted_columns = set(weighted_df.columns)
    missing = []
    for dataset_result in dataset_results:
        for features in dataset_result.features_by_type.values():
            for feature in features:
                if feature not in available_value_columns or feature not in available_weighted_columns:
                    missing.append(feature)
    if missing:
        preview = ", ".join(sorted(set(missing))[:10])
        raise KeyError(f"Missing features in methylation or weighted methylation table: {preview}")


def run_all() -> None:
    configure_matplotlib()
    tsv_path, weighted_tsv_path, import_dir, output_dir = resolve_paths()
    main_df = load_main_table(tsv_path)
    weighted_meth_df = load_weighted_table(weighted_tsv_path)
    weighted_meth_df = align_weighted_inputs(main_df, weighted_meth_df)

    dataset_results = [parse_import_results(import_dir / filename) for filename, _dataset_name in DATASET_CONFIGS]
    validate_feature_columns(main_df, weighted_meth_df, dataset_results)

    generated_outputs = []
    for dataset_result in dataset_results:
        grouped_value_df, grouped_weighted_source_df = prepare_dataset_frames(main_df, weighted_meth_df, dataset_result.dataset_name)
        scatter_pdf, scatter_png = plot_feature_scatter_figure(grouped_value_df, dataset_result, output_dir)
        weighted_df = build_weighted_type_table(grouped_value_df, grouped_weighted_source_df, dataset_result)
        weighted_stats_df = summarize_weighted_type_statistics(weighted_df)
        box_plots, stats_path = plot_weighted_boxplot_figure(
            dataset_result,
            weighted_df,
            weighted_stats_df,
            output_dir,
        )
        generated_outputs.append((scatter_pdf, scatter_png, box_plots, stats_path))

    print(f"Input table: {tsv_path}")
    print(f"Weighted methylation table: {weighted_tsv_path}")
    print(f"Import-results directory: {import_dir}")
    print(f"Output directory: {output_dir}")
    for scatter_pdf, scatter_png, box_plots, stats_path in generated_outputs:
        box_names = ", ".join(f"{pdf_path.name}, {png_path.name}" for _ft, pdf_path, png_path in box_plots)
        print(f"Generated: {scatter_pdf.name}, {scatter_png.name}, {box_names}, {stats_path.name}")


if __name__ == "__main__":
    run_all()
