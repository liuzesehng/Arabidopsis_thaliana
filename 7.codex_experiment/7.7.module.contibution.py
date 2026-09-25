#!/usr/bin/env python3
"""Plot one transcription factor across the RCA temperature series.

The target symbol is resolved through the annotation table.  For each of the
eight temperatures, all three biological replicates are transformed as
ln(TPM + 1), displayed individually, and summarized as mean +/- sample SD.
The line joining temperature means is descriptive only.  The sole inferential
annotation is a two-sided Spearman correlation across the eight means.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Sequence


# Matplotlib/fontconfig cannot write to the default cache in the shared runtime.
_CACHE_ROOT = Path(tempfile.gettempdir()) / "lcl5-temperature-figure-cache"
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
DEFAULT_ANNOTATION_FILE = DEFAULT_WORK_DIR / "total.gene.TF_ids.annotated.txt"
DEFAULT_SYMBOL = "LCL5"
DEFAULT_DPI = 600

EXPECTED_TEMPERATURES = np.array([4, 8, 12, 16, 20, 24, 28, 32], dtype=int)
EXPECTED_REPLICATES = 3
SAMPLE_COLUMN_PATTERN = re.compile(
    r"^(?P<temperature>\d+)_(?P<replicate>[^_]+)$"
)
SAFE_STEM_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
REPLICATE_OFFSETS_C = np.array([-0.42, 0.0, 0.42], dtype=float)

COLORS = {
    "replicate": "#85AEC5",
    "mean": "#225E7A",
    "control": "#E9E9E9",
    "grid": "#E3E5E7",
    "text": "#222222",
}

FIGURE_WIDTH_MM = 89.0
FIGURE_HEIGHT_MM = 75.0
FIGURE_WIDTH_IN = FIGURE_WIDTH_MM / 25.4
FIGURE_HEIGHT_IN = FIGURE_HEIGHT_MM / 25.4


def configure_style() -> None:
    """Apply a compact, publication-oriented Matplotlib style."""
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": 6.6,
            "axes.titlesize": 8.2,
            "axes.labelsize": 7.0,
            "axes.linewidth": 0.75,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.labelsize": 6.3,
            "ytick.labelsize": 6.3,
            "xtick.major.size": 2.7,
            "ytick.major.size": 2.7,
            "xtick.major.width": 0.65,
            "ytick.major.width": 0.65,
            "legend.frameon": False,
            "savefig.facecolor": "white",
        }
    )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot a TF symbol across eight RCA treatment temperatures, showing "
            "all biological replicates, mean +/- SD, and a Spearman correlation "
            "across the eight temperature means."
        )
    )
    parser.add_argument(
        "--expression-file",
        type=Path,
        default=DEFAULT_EXPRESSION_FILE,
        help=f"Expression matrix TSV (default: {DEFAULT_EXPRESSION_FILE})",
    )
    parser.add_argument(
        "--annotation-file",
        type=Path,
        default=DEFAULT_ANNOTATION_FILE,
        help=f"Gene annotation TSV (default: {DEFAULT_ANNOTATION_FILE})",
    )
    parser.add_argument(
        "--symbol",
        default=DEFAULT_SYMBOL,
        help=f"Exact gene symbol to plot (default: {DEFAULT_SYMBOL})",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_WORK_DIR,
        help=f"Output directory (default: {DEFAULT_WORK_DIR})",
    )
    parser.add_argument(
        "--basename",
        default=None,
        help="Output filename stem (default: <symbol>_temperature_scatter)",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=DEFAULT_DPI,
        help=f"PNG resolution in dots per inch (default: {DEFAULT_DPI})",
    )
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> str:
    if not args.symbol or args.symbol.strip() != args.symbol:
        raise ValueError("--symbol must be non-empty and have no surrounding whitespace")
    if not SAFE_STEM_PATTERN.fullmatch(args.symbol):
        raise ValueError(
            "--symbol may contain only letters, digits, periods, underscores, and hyphens"
        )
    if args.dpi < 300:
        raise ValueError("--dpi must be at least 300")

    basename = args.basename or f"{args.symbol}_temperature_scatter"
    if not basename or Path(basename).name != basename:
        raise ValueError("--basename must be a non-empty filename stem, not a path")
    if not SAFE_STEM_PATTERN.fullmatch(basename):
        raise ValueError(
            "--basename may contain only letters, digits, periods, underscores, and hyphens"
        )
    return basename


def resolve_symbol(annotation_file: Path, symbol: str) -> str:
    """Resolve an exact symbol to exactly one unique gene identifier."""
    if not annotation_file.is_file():
        raise FileNotFoundError(f"annotation file not found: {annotation_file}")

    annotation = pd.read_csv(
        annotation_file,
        sep="\t",
        encoding="utf-8",
        dtype=str,
        keep_default_na=False,
    )
    required = {"Gene_ID", "Symbol"}
    missing = sorted(required.difference(annotation.columns))
    if missing:
        raise ValueError(
            f"{annotation_file} is missing required columns: {', '.join(missing)}"
        )
    if annotation.empty:
        raise ValueError(f"annotation file contains no data rows: {annotation_file}")

    matches = annotation.loc[annotation["Symbol"].eq(symbol), "Gene_ID"]
    if matches.empty:
        raise ValueError(f"symbol {symbol!r} was not found in {annotation_file}")
    if matches.eq("").any():
        raise ValueError(f"symbol {symbol!r} has an empty Gene_ID annotation")
    if matches.str.strip().ne(matches).any():
        raise ValueError(f"symbol {symbol!r} has a Gene_ID with surrounding whitespace")

    gene_ids = pd.unique(matches)
    if len(gene_ids) != 1:
        preview = ", ".join(map(str, gene_ids[:5]))
        raise ValueError(
            f"symbol {symbol!r} maps to {len(gene_ids)} unique gene IDs: {preview}"
        )
    return str(gene_ids[0])


def sample_columns_by_temperature(
    columns: Sequence[str],
) -> dict[int, list[tuple[str, str]]]:
    """Validate and order the complete eight-temperature sample design."""
    grouped: dict[int, list[tuple[str, str]]] = {}
    invalid: list[str] = []
    for column in columns:
        column_text = str(column)
        match = SAMPLE_COLUMN_PATTERN.fullmatch(column_text)
        if match is None:
            invalid.append(column_text)
            continue
        temperature = int(match.group("temperature"))
        replicate = match.group("replicate")
        grouped.setdefault(temperature, []).append((column_text, replicate))

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
        if len(grouped[int(temperature)]) != EXPECTED_REPLICATES
    }
    if bad_counts:
        observed_text = ", ".join(
            f"{temperature} C={count}" for temperature, count in bad_counts.items()
        )
        raise ValueError(
            f"each temperature must have exactly {EXPECTED_REPLICATES} replicates; "
            f"observed {observed_text}"
        )

    ordered: dict[int, list[tuple[str, str]]] = {}
    for temperature in EXPECTED_TEMPERATURES:
        pairs = grouped[int(temperature)]
        replicate_ids = [replicate for _, replicate in pairs]
        if len(replicate_ids) != len(set(replicate_ids)):
            raise ValueError(f"duplicate replicate IDs at {int(temperature)} C")
        ordered[int(temperature)] = sorted(pairs, key=lambda pair: pair[1])
    return ordered


def load_replicates(
    expression_file: Path,
    *,
    gene_id: str,
    symbol: str,
) -> pd.DataFrame:
    """Load and validate all 24 target-gene measurements."""
    if not expression_file.is_file():
        raise FileNotFoundError(f"expression file not found: {expression_file}")

    expression = pd.read_csv(expression_file, sep="\t", encoding="utf-8")
    if "Name" not in expression.columns:
        raise ValueError(f"{expression_file} is missing required column: Name")
    if expression.empty:
        raise ValueError(f"expression file contains no data rows: {expression_file}")

    sample_columns = [str(column) for column in expression.columns if column != "Name"]
    if not sample_columns:
        raise ValueError("expression file contains no sample columns")
    grouped = sample_columns_by_temperature(sample_columns)

    identifiers = expression["Name"].astype(str)
    target_rows = expression.loc[identifiers.eq(gene_id)]
    if target_rows.empty:
        raise ValueError(f"gene ID {gene_id!r} was not found in {expression_file}")
    if len(target_rows) != 1:
        raise ValueError(
            f"gene ID {gene_id!r} occurs {len(target_rows)} times in {expression_file}"
        )

    records: list[dict[str, object]] = []
    row = target_rows.iloc[0]
    for temperature in EXPECTED_TEMPERATURES:
        pairs = grouped[int(temperature)]
        columns = [column for column, _ in pairs]
        try:
            numeric = pd.to_numeric(row.loc[columns], errors="raise").to_numpy(
                dtype=float
            )
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"gene ID {gene_id!r} contains non-numeric TPM values at "
                f"{int(temperature)} C"
            ) from exc
        if not np.isfinite(numeric).all():
            raise ValueError(
                f"gene ID {gene_id!r} contains missing or infinite TPM values at "
                f"{int(temperature)} C"
            )
        if (numeric < 0).any():
            raise ValueError(
                f"gene ID {gene_id!r} contains negative TPM values at "
                f"{int(temperature)} C"
            )

        transformed = np.log1p(numeric)
        if not np.isfinite(transformed).all():
            raise RuntimeError("ln(TPM + 1) produced non-finite values")
        for order, ((column, replicate), tpm, log_value) in enumerate(
            zip(pairs, numeric, transformed, strict=True), start=1
        ):
            records.append(
                {
                    "gene_id": gene_id,
                    "symbol": symbol,
                    "temperature_c": int(temperature),
                    "replicate_id": replicate,
                    "sample_column": column,
                    "replicate_order": order,
                    "tpm": float(tpm),
                    "ln_tpm_plus1": float(log_value),
                    "replicate_unit": "biological replicate",
                }
            )

    replicates = pd.DataFrame.from_records(records)
    expected_rows = len(EXPECTED_TEMPERATURES) * EXPECTED_REPLICATES
    if len(replicates) != expected_rows:
        raise RuntimeError(
            f"replicate table has {len(replicates)} rows; expected {expected_rows}"
        )
    return replicates


def summarize_replicates(replicates: pd.DataFrame) -> pd.DataFrame:
    """Compute mean and sample SD on the plotted ln(TPM + 1) scale."""
    summary = (
        replicates.groupby("temperature_c", sort=True, observed=False)["ln_tpm_plus1"]
        .agg(
            n_replicates="size",
            mean_ln_tpm_plus1="mean",
            sd_ln_tpm_plus1="std",
        )
        .reindex(EXPECTED_TEMPERATURES)
        .reset_index()
    )
    if summary.isna().any().any():
        raise RuntimeError("temperature summary contains missing values")
    if not summary["n_replicates"].eq(EXPECTED_REPLICATES).all():
        raise RuntimeError("one or more temperatures do not contain three replicates")
    if summary["mean_ln_tpm_plus1"].nunique() < 2:
        raise ValueError(
            "Spearman correlation is undefined because temperature means are constant"
        )
    return summary


def calculate_spearman(
    summary: pd.DataFrame,
    *,
    gene_id: str,
    symbol: str,
) -> pd.DataFrame:
    """Calculate the requested two-sided test across eight temperature means."""
    result = spearmanr(
        summary["temperature_c"].to_numpy(dtype=float),
        summary["mean_ln_tpm_plus1"].to_numpy(dtype=float),
        alternative="two-sided",
    )
    rho = float(result.statistic)
    p_value = float(result.pvalue)
    if not np.isfinite(rho) or not np.isfinite(p_value):
        raise RuntimeError("Spearman correlation returned a non-finite result")

    return pd.DataFrame(
        [
            {
                "gene_id": gene_id,
                "symbol": symbol,
                "method": "Spearman rank correlation",
                "alternative": "two-sided",
                "x_variable": "treatment_temperature_c",
                "y_variable": "mean_ln_tpm_plus1",
                "n_temperatures": len(summary),
                "rho": rho,
                "p_value": p_value,
                "replicate_summary": "arithmetic mean across 3 biological replicates",
                "transform": "ln(TPM + 1) applied before temperature summary",
                "spread_displayed": "sample SD across 3 biological replicates (ddof=1)",
                "multiple_testing_correction": "not applicable (one test)",
                "trend_line_role": "descriptive only",
                "control_role": "24 C visual reference only; no control comparison test",
            }
        ]
    )


def format_p_value(p_value: float) -> str:
    if p_value < 0.0001:
        mantissa, exponent = f"{p_value:.2e}".split("e")
        return f"{mantissa}e{int(exponent):+d}".replace("-", "\N{MINUS SIGN}")
    if p_value < 0.01:
        return f"{p_value:.4f}"
    return f"{p_value:.4f}"


def build_figure(
    replicates: pd.DataFrame,
    summary: pd.DataFrame,
    spearman_table: pd.DataFrame,
    *,
    gene_id: str,
    symbol: str,
) -> tuple[Figure, Axes]:
    """Build the single-panel quantitative figure."""
    configure_style()
    fig = plt.figure(figsize=(FIGURE_WIDTH_IN, FIGURE_HEIGHT_IN))
    ax = fig.add_axes([0.17, 0.17, 0.80, 0.52])

    ax.axvspan(22, 26, color=COLORS["control"], alpha=0.9, linewidth=0, zorder=0)
    for temperature in EXPECTED_TEMPERATURES:
        rows = replicates.loc[replicates["temperature_c"].eq(int(temperature))]
        rows = rows.sort_values("replicate_order")
        values = rows["ln_tpm_plus1"].to_numpy(dtype=float)
        ax.scatter(
            np.full(EXPECTED_REPLICATES, float(temperature)) + REPLICATE_OFFSETS_C,
            values,
            s=18.0,
            color=COLORS["replicate"],
            alpha=0.78,
            edgecolors="white",
            linewidths=0.45,
            zorder=3,
        )

    temperatures = summary["temperature_c"].to_numpy(dtype=float)
    means = summary["mean_ln_tpm_plus1"].to_numpy(dtype=float)
    standard_deviations = summary["sd_ln_tpm_plus1"].to_numpy(dtype=float)
    ax.plot(
        temperatures,
        means,
        color=COLORS["mean"],
        linewidth=1.15,
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
        elinewidth=0.9,
        capsize=2.5,
        capthick=0.9,
        zorder=5,
    )

    ax.set_xlim(2, 34)
    ax.set_xticks(EXPECTED_TEMPERATURES)
    ax.set_xticklabels([str(int(value)) for value in EXPECTED_TEMPERATURES])
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("ln(TPM + 1)")
    ax.tick_params(direction="out")

    plotted_upper = max(
        float(replicates["ln_tpm_plus1"].max()),
        float(np.max(means + standard_deviations)),
    )
    ax.set_ylim(0, plotted_upper * 1.28)

    rho = float(spearman_table.loc[0, "rho"])
    p_value = float(spearman_table.loc[0, "p_value"])
    rho_text = f"{rho:.3f}".replace("-", "\N{MINUS SIGN}")
    ax.text(
        0.025,
        0.965,
        f"ρ = {rho_text}, P = {format_p_value(p_value)}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=5.8,
        color=COLORS["text"],
        zorder=6,
    )
    fig.suptitle(
        f"{symbol} expression across temperatures",
        x=0.57,
        y=0.965,
        fontsize=8.2,
        fontweight="semibold",
    )
    fig.text(
        0.57,
        0.874,
        f"{gene_id}; points = biological replicates; mean ± SD",
        ha="center",
        va="center",
        fontsize=5.8,
        color=COLORS["text"],
    )
    fig.text(
        0.57,
        0.823,
        "Line = descriptive trend only; grey band = 24 °C control",
        ha="center",
        va="center",
        fontsize=5.8,
        color=COLORS["text"],
    )
    fig.text(
        0.57,
        0.755,
        "Spearman rank correlation (two-sided); n = 8 temperature means",
        ha="center",
        va="center",
        fontsize=5.8,
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


def run_pdf_text_audit(pdf_path: Path, json_path: Path) -> None:
    """Require every rendered PDF text run to be at least 5 pt."""
    audit_script = qa_scripts_directory() / "audit_pdf_text.py"
    if not audit_script.is_file():
        raise RuntimeError(
            "PDF text-audit helper not found; set NATURE_FIGURE_QA_DIR to the "
            "nature-figure scripts directory"
        )
    completed = subprocess.run(
        [
            sys.executable,
            str(audit_script),
            str(pdf_path),
            "--min-pt",
            "5",
            "--json",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        details = (completed.stdout + "\n" + completed.stderr).strip()
        raise RuntimeError(
            "PDF text audit failed" + (f": {details}" if details else "")
        )
    try:
        report = json.loads(completed.stdout)
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"could not parse PDF text-audit JSON: {exc}") from exc
    if not report.get("auditable") or int(report.get("below_minimum_count", 1)) != 0:
        raise RuntimeError("PDF text audit did not satisfy the 5 pt minimum")
    json_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def run_collision_audit(pdf_path: Path, json_path: Path) -> None:
    """Require a strict rendered-PDF collision audit."""
    audit_script = qa_scripts_directory() / "audit_figure_collisions.py"
    if not audit_script.is_file():
        raise RuntimeError(
            "PDF collision-audit helper not found; set NATURE_FIGURE_QA_DIR to the "
            "nature-figure scripts directory"
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
    replicates: pd.DataFrame,
    summary: pd.DataFrame,
    spearman_table: pd.DataFrame,
    *,
    output_dir: Path,
    basename: str,
    dpi: int,
) -> list[Path]:
    """Stage, audit, and publish the requested compact output bundle."""
    output_dir.mkdir(parents=True, exist_ok=True)
    suffixes = (
        ".pdf",
        ".png",
        ".replicates.tsv",
        ".summary.tsv",
        ".spearman.tsv",
        ".alignment.json",
        ".pdf-text-audit.json",
        ".collision-audit.json",
    )
    final_paths = [output_dir / f"{basename}{suffix}" for suffix in suffixes]

    with tempfile.TemporaryDirectory(prefix=f".{basename}.", dir=output_dir) as temp_name:
        temp_dir = Path(temp_name)
        temp_paths = [temp_dir / path.name for path in final_paths]
        (
            pdf_path,
            png_path,
            replicates_path,
            summary_path,
            spearman_path,
            alignment_path,
            pdf_text_path,
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

        write_tsv(replicates, replicates_path)
        write_tsv(summary, summary_path)
        write_tsv(spearman_table, spearman_path)
        fig.savefig(pdf_path, format="pdf", facecolor="white")
        fig.savefig(png_path, format="png", dpi=dpi, facecolor="white")
        run_pdf_text_audit(pdf_path, pdf_text_path)
        run_collision_audit(pdf_path, collision_path)

        pdf_text_report = json.loads(pdf_text_path.read_text(encoding="utf-8"))
        pdf_text_report["pdf"] = str(final_paths[0])
        pdf_text_path.write_text(
            json.dumps(pdf_text_report, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
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
        basename = validate_args(args)
        gene_id = resolve_symbol(args.annotation_file, args.symbol)
        replicates = load_replicates(
            args.expression_file,
            gene_id=gene_id,
            symbol=args.symbol,
        )
        summary = summarize_replicates(replicates)
        spearman_table = calculate_spearman(
            summary,
            gene_id=gene_id,
            symbol=args.symbol,
        )
        fig, ax = build_figure(
            replicates,
            summary,
            spearman_table,
            gene_id=gene_id,
            symbol=args.symbol,
        )
        try:
            outputs = export_outputs(
                fig,
                ax,
                replicates,
                summary,
                spearman_table,
                output_dir=args.output_dir,
                basename=basename,
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
        f"Validated {args.symbol} ({gene_id}) across "
        f"{len(EXPECTED_TEMPERATURES)} temperatures "
        f"({EXPECTED_REPLICATES} biological replicates per temperature)."
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
