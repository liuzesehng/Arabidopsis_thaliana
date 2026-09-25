#!/usr/bin/env python3
"""Plot yellow-module gene expression across treatment temperatures.

Each plotted gene value is the arithmetic mean of three biological samples at
one temperature, transformed as ln(mean TPM + 1). The dark trend connects the
eight module means and is descriptive only; inferential annotation is a
two-sided Spearman correlation across those eight temperature-level means.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Sequence


# Matplotlib/fontconfig cannot write to the default cache in the shared runtime.
_CACHE_ROOT = Path(tempfile.gettempdir()) / "yellow-module-figure-cache"
os.environ.setdefault("MPLCONFIGDIR", str(_CACHE_ROOT / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(_CACHE_ROOT / "xdg"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[2]
DEFAULT_WORK_DIR = PROJECT_ROOT / "list" / "coex_experiment"
DEFAULT_EXPRESSION_FILE = DEFAULT_WORK_DIR / "RCA.tran.txt"
DEFAULT_MODULE_FILE = DEFAULT_WORK_DIR / "yellow_module_genes.txt"
DEFAULT_BASENAME = "yellow_module_temperature_scatter"
DEFAULT_DPI = 600

EXPECTED_TEMPERATURES = np.array([4, 8, 12, 16, 20, 24, 28, 32], dtype=int)
EXPECTED_REPLICATES = 3
JITTER_SEED = 20260905
JITTER_WIDTH_C = 0.65
SAMPLE_COLUMN_PATTERN = re.compile(r"^(?P<temperature>\d+)_(?P<replicate>[^_]+)$")

COLORS = {
    "points": "#7B8794",
    "mean": "#355C7D",
    "control": "#E8E8E8",
    "grid": "#E5E5E5",
    "text": "#222222",
}

FIGURE_WIDTH_MM = 183.0
FIGURE_HEIGHT_MM = 105.0
FIGURE_WIDTH_IN = FIGURE_WIDTH_MM / 25.4
FIGURE_HEIGHT_IN = FIGURE_HEIGHT_MM / 25.4


def configure_style() -> None:
    """Apply a compact publication-oriented Matplotlib style."""
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": 7.0,
            "axes.titlesize": 8.0,
            "axes.labelsize": 7.2,
            "axes.linewidth": 0.75,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.labelsize": 6.8,
            "ytick.labelsize": 6.8,
            "xtick.major.size": 3.0,
            "ytick.major.size": 3.0,
            "xtick.major.width": 0.65,
            "ytick.major.width": 0.65,
            "legend.frameon": False,
            "savefig.facecolor": "white",
        }
    )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot yellow-module gene expression across eight temperatures and "
            "annotate the Spearman correlation of the eight module means."
        )
    )
    parser.add_argument(
        "--expression-file",
        type=Path,
        default=DEFAULT_EXPRESSION_FILE,
        help=f"Expression matrix TSV (default: {DEFAULT_EXPRESSION_FILE})",
    )
    parser.add_argument(
        "--module-file",
        type=Path,
        default=DEFAULT_MODULE_FILE,
        help=f"Yellow-module gene table TSV (default: {DEFAULT_MODULE_FILE})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_WORK_DIR,
        help=f"Output directory (default: {DEFAULT_WORK_DIR})",
    )
    parser.add_argument(
        "--basename",
        default=DEFAULT_BASENAME,
        help=f"Output filename stem (default: {DEFAULT_BASENAME})",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=DEFAULT_DPI,
        help=f"PNG/TIFF resolution in dots per inch (default: {DEFAULT_DPI})",
    )
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> None:
    if args.dpi < 300:
        raise ValueError("--dpi must be at least 300")
    if not args.basename or Path(args.basename).name != args.basename:
        raise ValueError("--basename must be a non-empty filename stem, not a path")


def _validate_gene_ids(table: pd.DataFrame, column: str, source: Path) -> pd.Index:
    if column not in table.columns:
        raise ValueError(f"{source} is missing required column: {column}")
    if table[column].isna().any():
        raise ValueError(f"{source} contains missing gene IDs in column {column}")

    identifiers = table[column].astype(str)
    if identifiers.str.len().eq(0).any():
        raise ValueError(f"{source} contains empty gene IDs in column {column}")
    if identifiers.str.strip().ne(identifiers).any():
        raise ValueError(f"{source} contains gene IDs with leading/trailing whitespace")
    duplicates = identifiers[identifiers.duplicated(keep=False)].unique().tolist()
    if duplicates:
        preview = ", ".join(duplicates[:5])
        raise ValueError(f"{source} contains duplicate gene IDs: {preview}")
    return pd.Index(identifiers, name=column)


def _sample_columns_by_temperature(columns: Sequence[str]) -> dict[int, list[str]]:
    grouped: dict[int, list[str]] = {}
    invalid: list[str] = []
    for column in columns:
        match = SAMPLE_COLUMN_PATTERN.fullmatch(str(column))
        if match is None:
            invalid.append(str(column))
            continue
        temperature = int(match.group("temperature"))
        grouped.setdefault(temperature, []).append(str(column))

    if invalid:
        raise ValueError(
            "expression sample columns must be named '<temperature>_<replicate>'; "
            f"invalid columns: {', '.join(invalid[:5])}"
        )

    expected = set(EXPECTED_TEMPERATURES.tolist())
    observed = set(grouped)
    if observed != expected:
        missing = sorted(expected - observed)
        unexpected = sorted(observed - expected)
        details: list[str] = []
        if missing:
            details.append("missing temperatures=" + ",".join(map(str, missing)))
        if unexpected:
            details.append("unexpected temperatures=" + ",".join(map(str, unexpected)))
        raise ValueError("temperature columns are invalid: " + "; ".join(details))

    bad_counts = {
        temperature: len(grouped[temperature])
        for temperature in EXPECTED_TEMPERATURES
        if len(grouped[temperature]) != EXPECTED_REPLICATES
    }
    if bad_counts:
        observed_text = ", ".join(
            f"{temperature} °C={count}" for temperature, count in bad_counts.items()
        )
        raise ValueError(
            f"each temperature must have exactly {EXPECTED_REPLICATES} replicates; "
            f"observed {observed_text}"
        )
    return {int(t): grouped[int(t)] for t in EXPECTED_TEMPERATURES}


def load_and_validate_inputs(
    expression_file: Path, module_file: Path
) -> tuple[pd.DataFrame, pd.Index, dict[int, list[str]]]:
    """Load both inputs and enforce the gene and experiment contracts."""
    if not expression_file.is_file():
        raise FileNotFoundError(f"expression file not found: {expression_file}")
    if not module_file.is_file():
        raise FileNotFoundError(f"module file not found: {module_file}")

    expression = pd.read_csv(expression_file, sep="\t", encoding="utf-8")
    module = pd.read_csv(module_file, sep="\t", encoding="utf-8")
    if expression.empty:
        raise ValueError(f"expression file contains no data rows: {expression_file}")
    if module.empty:
        raise ValueError(f"module file contains no data rows: {module_file}")

    expression_ids = _validate_gene_ids(expression, "Name", expression_file)
    module_ids = _validate_gene_ids(module, "datSummary", module_file)

    if "colorNEW" not in module.columns:
        raise ValueError(f"{module_file} is missing required column: colorNEW")
    if module["colorNEW"].isna().any():
        raise ValueError(f"{module_file} contains missing colorNEW values")
    module_colors = module["colorNEW"].astype(str)
    if not module_colors.eq("yellow").all():
        observed = ", ".join(sorted(module_colors.unique())[:5])
        raise ValueError(
            "module file must contain only colorNEW='yellow' rows; "
            f"observed: {observed}"
        )
    if len(module_ids) < 2:
        raise ValueError("at least two yellow-module genes are required to compute SD")

    sample_columns = [str(column) for column in expression.columns if column != "Name"]
    if not sample_columns:
        raise ValueError("expression file contains no sample columns")
    columns_by_temperature = _sample_columns_by_temperature(sample_columns)

    numeric = expression.loc[:, sample_columns].apply(pd.to_numeric, errors="raise")
    numeric_array = numeric.to_numpy(dtype=float)
    if not np.isfinite(numeric_array).all():
        raise ValueError("expression matrix contains missing or infinite values")
    if (numeric_array < 0).any():
        raise ValueError("expression matrix contains negative TPM values")

    numeric.index = expression_ids
    missing_genes = module_ids.difference(numeric.index, sort=False)
    if len(missing_genes):
        preview = ", ".join(missing_genes[:5].tolist())
        raise ValueError(
            f"{len(missing_genes)} module genes are absent from the expression matrix: "
            f"{preview}"
        )

    selected = numeric.loc[module_ids].copy()
    return selected, module_ids, columns_by_temperature


def calculate_gene_means(
    expression: pd.DataFrame,
    module_ids: pd.Index,
    columns_by_temperature: dict[int, list[str]],
) -> pd.DataFrame:
    """Return one transformed three-replicate mean per gene and temperature."""
    records: list[pd.DataFrame] = []
    for temperature in EXPECTED_TEMPERATURES:
        columns = columns_by_temperature[int(temperature)]
        mean_tpm = expression.loc[module_ids, columns].mean(axis=1)
        records.append(
            pd.DataFrame(
                {
                    "gene_id": module_ids.to_numpy(),
                    "temperature_c": int(temperature),
                    "replicate_count": len(columns),
                    "mean_tpm": mean_tpm.to_numpy(dtype=float),
                    "ln_mean_tpm_plus1": np.log1p(mean_tpm.to_numpy(dtype=float)),
                }
            )
        )

    gene_means = pd.concat(records, ignore_index=True)
    expected_rows = len(module_ids) * len(EXPECTED_TEMPERATURES)
    if len(gene_means) != expected_rows:
        raise RuntimeError(
            f"gene-temperature table has {len(gene_means)} rows; expected {expected_rows}"
        )
    values = gene_means["ln_mean_tpm_plus1"].to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise RuntimeError("ln(mean TPM + 1) produced non-finite values")
    return gene_means


def summarize_temperatures(gene_means: pd.DataFrame) -> pd.DataFrame:
    """Summarize the plotted gene points at each temperature as mean ± sample SD."""
    summary = (
        gene_means.groupby("temperature_c", sort=True, observed=False)[
            "ln_mean_tpm_plus1"
        ]
        .agg(n_genes="size", mean_ln_mean_tpm_plus1="mean", sd_ln_mean_tpm_plus1="std")
        .reindex(EXPECTED_TEMPERATURES)
        .reset_index()
    )
    if summary.isna().any().any():
        raise RuntimeError("temperature summary contains missing values")
    expected_genes = int(gene_means["gene_id"].nunique())
    if not summary["n_genes"].eq(expected_genes).all():
        raise RuntimeError("one or more temperatures do not contain every module gene")
    if summary["mean_ln_mean_tpm_plus1"].nunique() < 2:
        raise ValueError("Spearman correlation is undefined because module means are constant")
    return summary


def calculate_spearman(summary: pd.DataFrame) -> pd.DataFrame:
    """Calculate the requested two-sided correlation across eight module means."""
    result = spearmanr(
        summary["temperature_c"].to_numpy(dtype=float),
        summary["mean_ln_mean_tpm_plus1"].to_numpy(dtype=float),
        alternative="two-sided",
    )
    rho = float(result.statistic)
    p_value = float(result.pvalue)
    if not np.isfinite(rho) or not np.isfinite(p_value):
        raise RuntimeError("Spearman correlation returned a non-finite result")

    return pd.DataFrame(
        [
            {
                "method": "Spearman rank correlation",
                "alternative": "two-sided",
                "x_variable": "treatment_temperature_c",
                "y_variable": "module_mean_ln_mean_tpm_plus1",
                "n_temperatures": len(summary),
                "rho": rho,
                "p_value": p_value,
                "replicate_summary": "arithmetic mean across 3 samples per gene and temperature",
                "transform": "ln(mean_TPM + 1)",
                "module_summary": "arithmetic mean across yellow-module genes",
                "spread_displayed": "sample SD across yellow-module genes (ddof=1)",
                "multiple_testing_correction": "not applicable (one test)",
                "trend_line_role": "descriptive only",
            }
        ]
    )


def format_p_value(p_value: float) -> str:
    if p_value < 0.0001:
        mantissa, exponent = f"{p_value:.2e}".split("e")
        return f"{mantissa}e{int(exponent):+d}".replace("-", "−")
    if p_value < 0.01:
        return f"{p_value:.4f}"
    return f"{p_value:.3f}"


def deterministic_jitter(gene_ids: Sequence[str], temperature: int) -> np.ndarray:
    """Map gene IDs reproducibly into the permitted horizontal jitter interval."""
    offsets: list[float] = []
    denominator = float((1 << 64) - 1)
    for gene_id in gene_ids:
        key = f"{JITTER_SEED}:{temperature}:{gene_id}".encode("utf-8")
        digest = hashlib.blake2b(key, digest_size=8).digest()
        unit_interval = int.from_bytes(digest, byteorder="big") / denominator
        offsets.append((2.0 * unit_interval - 1.0) * JITTER_WIDTH_C)
    return np.asarray(offsets, dtype=float)


def build_figure(
    gene_means: pd.DataFrame, summary: pd.DataFrame, spearman_table: pd.DataFrame
) -> tuple[Figure, Axes]:
    """Build the single-panel quantitative figure."""
    configure_style()
    fig = plt.figure(figsize=(FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN))
    ax = fig.add_axes([0.09, 0.16, 0.89, 0.66])

    ax.axvspan(22, 26, color=COLORS["control"], alpha=0.85, linewidth=0, zorder=0)
    for temperature in EXPECTED_TEMPERATURES:
        temperature_rows = gene_means.loc[
            gene_means["temperature_c"].eq(int(temperature))
        ]
        values = temperature_rows["ln_mean_tpm_plus1"].to_numpy(dtype=float)
        jitter = deterministic_jitter(
            temperature_rows["gene_id"].astype(str).tolist(), int(temperature)
        )
        ax.scatter(
            np.full(len(values), float(temperature)) + jitter,
            values,
            s=7.0,
            color=COLORS["points"],
            alpha=0.24,
            edgecolors="none",
            zorder=2,
        )

    temperatures = summary["temperature_c"].to_numpy(dtype=float)
    means = summary["mean_ln_mean_tpm_plus1"].to_numpy(dtype=float)
    standard_deviations = summary["sd_ln_mean_tpm_plus1"].to_numpy(dtype=float)
    ax.plot(
        temperatures,
        means,
        color=COLORS["mean"],
        linewidth=1.25,
        alpha=0.95,
        zorder=4,
    )
    ax.errorbar(
        temperatures,
        means,
        yerr=standard_deviations,
        fmt="D",
        markersize=4.0,
        markerfacecolor=COLORS["mean"],
        markeredgecolor="white",
        markeredgewidth=0.55,
        color=COLORS["mean"],
        ecolor=COLORS["mean"],
        elinewidth=1.0,
        capsize=3.0,
        capthick=1.0,
        zorder=5,
    )

    ax.set_xlim(2, 34)
    ax.set_xticks(EXPECTED_TEMPERATURES)
    ax.set_xticklabels([str(int(value)) for value in EXPECTED_TEMPERATURES])
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("ln(TPM + 1)")
    ax.grid(axis="y", color=COLORS["grid"], linewidth=0.5, zorder=0.5)
    ax.set_axisbelow(True)
    ax.tick_params(direction="out")

    plotted_upper = max(
        float(gene_means["ln_mean_tpm_plus1"].max()),
        float(np.max(means + standard_deviations)),
    )
    ax.set_ylim(0, plotted_upper + max(2.0, 0.20 * plotted_upper))
    tick_step = 2.0
    highest_data_tick = np.ceil(plotted_upper / tick_step) * tick_step
    ax.set_yticks(np.arange(0, highest_data_tick + 0.5 * tick_step, tick_step))

    rho = float(spearman_table.loc[0, "rho"])
    p_value = float(spearman_table.loc[0, "p_value"])
    rho_text = f"{rho:.3f}".replace("-", "−")
    ax.text(
        0.018,
        0.965,
        f"ρ = {rho_text}, P = {format_p_value(p_value)}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7.0,
        color=COLORS["text"],
        clip_on=False,
        zorder=10,
    )

    fig.suptitle(
        "Yellow-module gene expression across temperatures",
        x=0.535,
        y=0.955,
        fontsize=8.2,
        fontweight="semibold",
    )
    fig.text(
        0.535,
        0.895,
        "Points = genes (3-replicate mean); dark line and error bars = module mean ± SD "
        "(descriptive); grey band = 24 °C control",
        ha="center",
        va="center",
        fontsize=6.5,
        color=COLORS["text"],
    )
    fig.text(
        0.535,
        0.850,
        "Spearman correlation: n = 8 temperature means; two-sided",
        ha="center",
        va="center",
        fontsize=6.5,
        color=COLORS["text"],
    )
    return fig, ax


def qa_scripts_directory() -> Path:
    return Path(
        os.environ.get(
            "NATURE_FIGURE_QA_DIR",
            Path.home() / ".codex" / "skills" / "nature-figure" / "scripts",
        )
    )


def load_alignment_helper():
    directory = qa_scripts_directory()
    helper_path = directory / "audit_panel_alignment.py"
    if not helper_path.is_file():
        raise RuntimeError(
            "panel-alignment helper not found; set NATURE_FIGURE_QA_DIR to the "
            "directory containing audit_panel_alignment.py"
        )
    sys.path.insert(0, str(directory))
    try:
        from audit_panel_alignment import require_matplotlib_panel_alignment
    except ImportError as exc:
        raise RuntimeError(f"could not import panel-alignment helper: {exc}") from exc
    return require_matplotlib_panel_alignment


def run_collision_audit(pdf_path: Path, json_path: Path) -> None:
    """Require a strict rendered-PDF collision audit before publishing outputs."""
    audit_script = qa_scripts_directory() / "audit_figure_collisions.py"
    if not audit_script.is_file():
        raise RuntimeError(
            "PDF collision-audit helper not found; set NATURE_FIGURE_QA_DIR "
            "to the nature-figure scripts directory"
        )
    completed = subprocess.run(
        [
            sys.executable,
            str(audit_script),
            str(pdf_path),
            "--json-out",
            str(json_path),
            "--strict",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        details = (completed.stdout + "\n" + completed.stderr).strip()
        raise RuntimeError(
            "strict rendered-PDF collision audit failed"
            + (f": {details}" if details else "")
        )
    try:
        report = json.loads(json_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"could not read collision-audit JSON: {exc}") from exc
    if report.get("verdict") != "PASS":
        raise RuntimeError(
            f"collision audit did not pass (verdict: {report.get('verdict', 'unknown')})"
        )


def write_tsv(table: pd.DataFrame, path: Path) -> None:
    table.to_csv(path, sep="\t", index=False, encoding="utf-8", float_format="%.12g")


def export_outputs(
    fig: Figure,
    ax: Axes,
    gene_means: pd.DataFrame,
    summary: pd.DataFrame,
    spearman_table: pd.DataFrame,
    *,
    output_dir: Path,
    basename: str,
    dpi: int,
) -> list[Path]:
    """Stage, audit, and atomically publish the complete output bundle."""
    output_dir.mkdir(parents=True, exist_ok=True)
    suffixes = (
        ".pdf",
        ".svg",
        ".png",
        ".tiff",
        ".gene_means.tsv",
        ".summary.tsv",
        ".spearman.tsv",
        ".alignment.json",
        ".collision-audit.json",
    )
    final_paths = [output_dir / f"{basename}{suffix}" for suffix in suffixes]

    with tempfile.TemporaryDirectory(prefix=f".{basename}.", dir=output_dir) as temp_name:
        temp_dir = Path(temp_name)
        temp_paths = [temp_dir / path.name for path in final_paths]
        (
            pdf_path,
            svg_path,
            png_path,
            tiff_path,
            gene_means_path,
            summary_path,
            spearman_path,
            alignment_path,
            collision_path,
        ) = temp_paths

        require_matplotlib_panel_alignment = load_alignment_helper()
        alignment_report = require_matplotlib_panel_alignment(
            fig,
            axes=[ax],
            panel_ids=["a"],
            json_out=alignment_path,
            tolerance_pt=1.5,
            gutter_tolerance_pt=1.5,
            require_panel_labels=False,
            strict=True,
        )
        if alignment_report.get("verdict") != "NOT APPLICABLE":
            raise RuntimeError(
                "single-panel alignment audit did not return NOT APPLICABLE"
            )

        write_tsv(gene_means, gene_means_path)
        write_tsv(summary, summary_path)
        write_tsv(spearman_table, spearman_path)
        fig.savefig(pdf_path, format="pdf")
        fig.savefig(svg_path, format="svg")
        fig.savefig(png_path, format="png", dpi=dpi)
        fig.savefig(tiff_path, format="tiff", dpi=dpi)
        run_collision_audit(pdf_path, collision_path)

        collision_report = json.loads(collision_path.read_text(encoding="utf-8"))
        collision_report["pdf"] = str(final_paths[0])
        collision_path.write_text(
            json.dumps(collision_report, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )

        for staged, final in zip(temp_paths, final_paths, strict=True):
            os.replace(staged, final)
    return final_paths


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        validate_args(args)
        expression, module_ids, columns_by_temperature = load_and_validate_inputs(
            args.expression_file, args.module_file
        )
        gene_means = calculate_gene_means(
            expression, module_ids, columns_by_temperature
        )
        summary = summarize_temperatures(gene_means)
        spearman_table = calculate_spearman(summary)
        fig, ax = build_figure(gene_means, summary, spearman_table)
        try:
            outputs = export_outputs(
                fig,
                ax,
                gene_means,
                summary,
                spearman_table,
                output_dir=args.output_dir,
                basename=args.basename,
                dpi=args.dpi,
            )
        finally:
            plt.close(fig)
    except (
        FileNotFoundError,
        OSError,
        RuntimeError,
        subprocess.SubprocessError,
        TypeError,
        ValueError,
    ) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    rho = float(spearman_table.loc[0, "rho"])
    p_value = float(spearman_table.loc[0, "p_value"])
    print(
        f"Validated {len(module_ids)} yellow-module genes across "
        f"{len(EXPECTED_TEMPERATURES)} temperatures "
        f"({EXPECTED_REPLICATES} replicates per gene and temperature)."
    )
    print(
        f"Spearman across {len(summary)} temperature means: "
        f"rho={rho:.12g}, p={p_value:.12g} (two-sided)."
    )
    for output in outputs:
        print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
