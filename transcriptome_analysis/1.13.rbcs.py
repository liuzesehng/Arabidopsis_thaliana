from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr, spearmanr


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[2]
WORK_DIR = PROJECT_ROOT / "list" / "RCA" / "rbcs"
INPUT_MATRIX = PROJECT_ROOT / "list" / "coex2" / "RCA.tran.txt"

RCA_NAME = "RCA"
RCA_ID = "AT2G39730"
RBCS_GENES = {
    "RBCS1A": "AT1G67090",
    "RBCS1B": "AT5G38430",
    "RBCS2B": "AT5G38420",
    "RBCS3B": "AT5G38410",
}


def load_expression_matrix(matrix_path: Path) -> pd.DataFrame:
    matrix = pd.read_csv(matrix_path, sep="\t")
    if "Name" not in matrix.columns:
        raise ValueError(f"Input matrix missing required 'Name' column: {matrix_path}")
    return matrix.set_index("Name")


def extract_gene_vector(matrix: pd.DataFrame, gene_id: str) -> pd.Series:
    if gene_id not in matrix.index:
        raise ValueError(f"Target gene not found in expression matrix: {gene_id}")
    series = pd.to_numeric(matrix.loc[gene_id], errors="coerce")
    if isinstance(series, pd.DataFrame):
        raise ValueError(f"Duplicate entries found for target gene: {gene_id}")
    return series


def calculate_correlations(rca_expr: pd.Series, gene_expr: pd.Series) -> tuple[pd.DataFrame, dict]:
    paired = pd.DataFrame({"RCA": rca_expr, "TARGET": gene_expr}).dropna()
    if paired.shape[0] < 2:
        raise ValueError("Not enough paired non-missing samples to calculate correlation.")

    pearson_r, pearson_p = pearsonr(paired["RCA"], paired["TARGET"])
    spearman_rho, spearman_p = spearmanr(paired["RCA"], paired["TARGET"])
    return paired, {
        "n_samples": paired.shape[0],
        "pearson_r": pearson_r,
        "pearson_p": pearson_p,
        "spearman_rho": spearman_rho,
        "spearman_p": spearman_p,
    }


def create_scatter_plot(
    paired: pd.DataFrame,
    target_gene_name: str,
    target_gene_id: str,
    stats_result: dict,
    output_path: Path,
) -> None:
    plt.figure(figsize=(6.5, 5.5))
    sns.set_theme(style="whitegrid")
    ax = sns.regplot(
        data=paired,
        x="RCA",
        y="TARGET",
        scatter_kws={"s": 28, "alpha": 0.75, "color": "#2f6c8f"},
        line_kws={"color": "#c44e52", "linewidth": 1.5},
    )
    ax.set_xlabel(f"{RCA_NAME} ({RCA_ID}) expression")
    ax.set_ylabel(f"{target_gene_name} ({target_gene_id}) expression")
    ax.set_title(f"{RCA_NAME} vs {target_gene_name}")

    annotation = (
        f"n = {stats_result['n_samples']}\n"
        f"Pearson r = {stats_result['pearson_r']:.4f}\n"
        f"Pearson p = {stats_result['pearson_p']:.3e}\n"
        f"Spearman rho = {stats_result['spearman_rho']:.4f}\n"
        f"Spearman p = {stats_result['spearman_p']:.3e}"
    )
    ax.text(
        0.03,
        0.97,
        annotation,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=10,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85, "edgecolor": "#bbbbbb"},
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def main() -> None:
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    expression_matrix = load_expression_matrix(INPUT_MATRIX)

    rca_expr = extract_gene_vector(expression_matrix, RCA_ID)
    results = []

    for gene_name, gene_id in RBCS_GENES.items():
        gene_expr = extract_gene_vector(expression_matrix, gene_id)
        paired, stats_result = calculate_correlations(rca_expr, gene_expr)

        results.append(
            {
                "gene_name": gene_name,
                "gene_id": gene_id,
                "rca_gene_name": RCA_NAME,
                "rca_gene_id": RCA_ID,
                **stats_result,
            }
        )

        plot_path = WORK_DIR / f"RCA_vs_{gene_name}.scatter.png"
        create_scatter_plot(paired, gene_name, gene_id, stats_result, plot_path)

    result_df = pd.DataFrame(results)
    result_df.to_csv(WORK_DIR / "RCA_RBCS_correlation.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
