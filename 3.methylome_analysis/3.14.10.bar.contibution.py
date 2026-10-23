#!/usr/bin/env python3
"""Create one publication-ready SHAP contribution bar chart per RCA model.

The five model result files are parsed using the same mean(|SHAP|) grouping
method as ``3.14.10.contibution.py``.  SNP, CG, CHG and CHH contributions are
summed within each model and divided by the sum over all parsed features.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import subprocess
import sys
import tempfile
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib-cache"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "font-cache"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


MM_PER_INCH = 25.4
FIGURE_WIDTH_MM = 89.0
FIGURE_HEIGHT_MM = 70.0
SECTION_HEADER = "  所有特征及其平均|SHAP|值:"
GROUP_ORDER = ("SNP", "CG", "CHG", "CHH")
GROUP_COLORS = {
    "SNP": "#A7ADB4",
    "CG": "#79A9B8",
    "CHG": "#D3A35D",
    "CHH": "#4E6F9E",
}
FEATURE_LINE_RE = re.compile(
    r"^\s{4}([A-Za-z0-9_./%+\-]+):\s*"
    r"([+\-]?(?:(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+\-]?\d+)?|"
    r"(?i:nan|inf(?:inity)?)))\s*$"
)
DEFAULT_DPI = 600


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


def _skip_alignment(*args, json_out=None, **kwargs):
    return _skip_optional_qa("audit_panel_alignment.py is unavailable", json_out)


@dataclass(frozen=True)
class ModelSpec:
    response_id: str
    response_label: str
    response_family: str
    filename: str


@dataclass(frozen=True)
class ModelContribution:
    spec: ModelSpec
    input_path: Path
    feature_count: int
    total_shap: float
    group_sums: OrderedDict[str, float]
    percentages: OrderedDict[str, float]


MODEL_SPECS = (
    ModelSpec("total_TPM", "Total expression", "Expression", "total_TPM_import_results.txt"),
    ModelSpec("A_TPM", "α expression", "Expression", "A_TPM_import_results.txt"),
    ModelSpec("B_TPM", "β expression", "Expression", "B_TPM_import_results.txt"),
    ModelSpec("A_per_TPM", "α proportion", "Proportion", "A_per_TPM_import_results.txt"),
    ModelSpec("B_per_TPM", "β proportion", "Proportion", "B_per_TPM_import_results.txt"),
)


def project_root() -> Path:
    """Resolve the Arabidopsis project root independently of the current directory."""
    return Path(__file__).resolve().parents[2]


def default_input_dir() -> Path:
    return project_root() / "list" / "xgboot" / "TPM_4.5_all"


def default_output_dir() -> Path:
    return project_root() / "list" / "RCA"


def configure_matplotlib() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 7.5,
            "axes.labelsize": 8.0,
            "xtick.labelsize": 7.3,
            "ytick.labelsize": 7.0,
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.fontsize": 7.2,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot SNP/CG/CHG/CHH mean(|SHAP|) contributions for five RCA models."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=default_input_dir(),
        help="Directory containing the five *_import_results.txt files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=default_output_dir(),
        help="Root directory containing one output subdirectory per RCA model.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=DEFAULT_DPI,
        help="Raster output resolution; must be at least 300 dpi (default: 600).",
    )
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> None:
    if not args.input_dir.is_dir():
        raise ValueError(f"input directory does not exist: {args.input_dir}")
    if args.output_dir.exists() and not args.output_dir.is_dir():
        raise ValueError(f"output path is not a directory: {args.output_dir}")
    if args.dpi < 300:
        raise ValueError(f"dpi must be at least 300, received {args.dpi}")


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
    if not input_path.is_file():
        raise FileNotFoundError(f"input file not found: {input_path}")

    lines = input_path.read_text(encoding="utf-8").splitlines()
    header_indices = [index for index, line in enumerate(lines) if line == SECTION_HEADER]
    if not header_indices:
        raise ValueError(f"SHAP feature section not found in {input_path}")
    if len(header_indices) != 1:
        raise ValueError(f"multiple SHAP feature sections found in {input_path}")

    contributions: list[tuple[str, float]] = []
    seen_features: set[str] = set()
    start_index = header_indices[0] + 1

    for line_number, line in enumerate(lines[start_index:], start=start_index + 1):
        if not line.strip():
            break
        match = FEATURE_LINE_RE.fullmatch(line)
        if match is None:
            raise ValueError(
                f"malformed SHAP feature entry in {input_path} at line {line_number}: {line!r}"
            )
        feature_name, value_text = match.groups()
        if feature_name in seen_features:
            raise ValueError(f"duplicate feature {feature_name!r} in {input_path}")
        value = float(value_text)
        if not math.isfinite(value):
            raise ValueError(f"non-finite SHAP value for {feature_name!r} in {input_path}")
        if value < 0.0:
            raise ValueError(f"negative mean(|SHAP|) for {feature_name!r} in {input_path}")
        if classify_feature(feature_name) is None:
            raise ValueError(f"unrecognized feature prefix for {feature_name!r} in {input_path}")
        seen_features.add(feature_name)
        contributions.append((feature_name, value))

    if not contributions:
        raise ValueError(f"no SHAP feature entries found in {input_path}")
    return contributions


def summarize_model(spec: ModelSpec, input_path: Path) -> ModelContribution:
    contributions = parse_feature_contributions(input_path)
    grouped_values: dict[str, list[float]] = {group: [] for group in GROUP_ORDER}
    for feature_name, shap_value in contributions:
        group = classify_feature(feature_name)
        if group is None:  # Defensive; the parser already rejects this condition.
            raise ValueError(f"unrecognized feature prefix for {feature_name!r}")
        grouped_values[group].append(shap_value)

    total_shap = math.fsum(value for _, value in contributions)
    if not math.isfinite(total_shap) or total_shap <= 0.0:
        raise ValueError(f"total mean(|SHAP|) must be positive in {input_path}")

    group_sums = OrderedDict(
        (group, math.fsum(grouped_values[group])) for group in GROUP_ORDER
    )
    grouped_total = math.fsum(group_sums.values())
    tolerance = max(1e-12, abs(total_shap) * 1e-12)
    if not math.isclose(grouped_total, total_shap, rel_tol=1e-12, abs_tol=tolerance):
        raise ValueError(
            f"grouped SHAP sum ({grouped_total:.12g}) does not equal the feature total "
            f"({total_shap:.12g}) in {input_path}"
        )

    percentages = OrderedDict(
        (group, group_sums[group] / total_shap * 100.0) for group in GROUP_ORDER
    )
    percentage_total = math.fsum(percentages.values())
    if not math.isclose(percentage_total, 100.0, rel_tol=1e-12, abs_tol=1e-10):
        raise ValueError(
            f"contribution percentages sum to {percentage_total:.12g}, not 100, in {input_path}"
        )

    return ModelContribution(
        spec=spec,
        input_path=input_path,
        feature_count=len(contributions),
        total_shap=total_shap,
        group_sums=group_sums,
        percentages=percentages,
    )


def load_contributions(input_dir: Path) -> list[ModelContribution]:
    return [summarize_model(spec, input_dir / spec.filename) for spec in MODEL_SPECS]


def sorted_groups(model: ModelContribution) -> tuple[str, ...]:
    """Return feature groups in descending contribution order with stable ties."""
    canonical_rank = {group: index for index, group in enumerate(GROUP_ORDER)}
    return tuple(
        sorted(
            GROUP_ORDER,
            key=lambda group: (-model.percentages[group], canonical_rank[group]),
        )
    )


def build_figure(model: ModelContribution):
    configure_matplotlib()
    figure_size = (FIGURE_WIDTH_MM / MM_PER_INCH, FIGURE_HEIGHT_MM / MM_PER_INCH)
    fig, ax = plt.subplots(figsize=figure_size)
    fig.subplots_adjust(left=0.19, right=0.97, bottom=0.20, top=0.96)

    plot_groups = sorted_groups(model)
    x_positions = list(range(len(plot_groups)))
    values = [model.percentages[group] for group in plot_groups]
    bars = ax.bar(
        x_positions,
        values,
        width=0.64,
        color=[GROUP_COLORS[group] for group in plot_groups],
        edgecolor="#30343A",
        linewidth=0.45,
        zorder=2,
    )
    for bar, value in zip(bars, values, strict=True):
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            value + 1.6,
            f"{value:.2f}%",
            ha="center",
            va="bottom",
            color="#202124",
            fontsize=7.0,
            fontweight="medium",
            clip_on=False,
            zorder=3,
        )

    ax.set_xlim(-0.58, len(plot_groups) - 0.42)
    ax.set_ylim(0.0, 80.0)
    ax.set_xticks(x_positions, plot_groups)
    ax.set_yticks([0, 20, 40, 60, 80])
    ax.set_xlabel("Feature group", labelpad=5)
    ax.set_ylabel("Percentage (%)", labelpad=6)
    ax.tick_params(axis="x", direction="out", length=3.0, width=0.8, pad=4)
    ax.tick_params(axis="y", direction="out", length=3.0, width=0.8, pad=3)
    ax.spines["left"].set_bounds(0, 80)
    ax.spines["bottom"].set_linewidth(0.8)

    fig.canvas.draw()
    return fig, ax


def write_source_data(model: ModelContribution, output_path: Path) -> None:
    fieldnames = [
        "response_id",
        "response_label",
        "response_family",
        "feature_group",
        "shap_sum",
        "percentage",
        "feature_count",
        "total_shap",
        "input_file",
    ]
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for group in sorted_groups(model):
            writer.writerow(
                {
                    "response_id": model.spec.response_id,
                    "response_label": model.spec.response_label,
                    "response_family": model.spec.response_family,
                    "feature_group": group,
                    "shap_sum": f"{model.group_sums[group]:.12g}",
                    "percentage": f"{model.percentages[group]:.12g}",
                    "feature_count": model.feature_count,
                    "total_shap": f"{model.total_shap:.12g}",
                    "input_file": str(model.input_path),
                }
            )


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
        _skip_optional_qa("Optional panel-alignment helper not found")
        return _skip_alignment
    sys.path.insert(0, str(directory))
    try:
        from audit_panel_alignment import require_matplotlib_panel_alignment
    except ImportError as exc:
        _skip_optional_qa(str(exc))
        return _skip_alignment
    return require_matplotlib_panel_alignment


def run_collision_audit(pdf_path: Path, json_path: Path) -> None:
    audit_script = qa_scripts_directory() / "audit_figure_collisions.py"
    if not audit_script.is_file():
        _skip_optional_qa(f"Optional QA tool not found: {audit_script}", json_path)
        return
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
        details = "\n".join(
            part.strip() for part in (completed.stdout, completed.stderr) if part.strip()
        )
        raise RuntimeError(
            "strict rendered-PDF collision audit failed"
            + (f":\n{details}" if details else "")
        )
    try:
        report = json.loads(json_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"could not read collision-audit JSON: {exc}") from exc
    if report.get("verdict") != "PASS":
        raise RuntimeError(
            f"collision audit did not pass (verdict: {report.get('verdict', 'unknown')})"
        )


def export_outputs(
    contributions: Sequence[ModelContribution],
    *,
    output_root: Path,
    dpi: int,
) -> list[Path]:
    output_root.mkdir(parents=True, exist_ok=True)
    suffixes = (
        ".svg",
        ".pdf",
        ".png",
        ".tiff",
        ".source_data.tsv",
        ".alignment.json",
        ".collision-audit.json",
    )
    published_paths: list[Path] = []
    staged_to_final: list[tuple[Path, Path]] = []
    require_matplotlib_panel_alignment = load_alignment_helper()

    with tempfile.TemporaryDirectory(
        prefix=".RCA_SHAP_group_contribution.", dir=output_root
    ) as temp_name:
        temp_root = Path(temp_name)

        for model in contributions:
            basename = f"{model.spec.response_id}_SHAP_group_contribution"
            final_dir = output_root / model.spec.response_id
            staged_dir = temp_root / model.spec.response_id
            staged_dir.mkdir(parents=True)

            final_paths = [final_dir / f"{basename}{suffix}" for suffix in suffixes]
            staged_paths = [staged_dir / path.name for path in final_paths]
            svg_path, pdf_path, png_path, tiff_path, data_path, alignment_path, collision_path = (
                staged_paths
            )

            fig, ax = build_figure(model)
            try:
                require_matplotlib_panel_alignment(
                    fig,
                    axes=[ax],
                    panel_ids=["a"],
                    json_out=alignment_path,
                    tolerance_pt=1.5,
                    gutter_tolerance_pt=1.5,
                    require_panel_labels=False,
                    strict=True,
                )

                write_source_data(model, data_path)
                fig.savefig(svg_path, format="svg")
                fig.savefig(pdf_path, format="pdf")
                fig.savefig(png_path, format="png", dpi=dpi)
                fig.savefig(
                    tiff_path,
                    format="tiff",
                    dpi=dpi,
                    pil_kwargs={"compression": "tiff_lzw"},
                )
                run_collision_audit(pdf_path, collision_path)
            finally:
                plt.close(fig)

            collision_report = json.loads(collision_path.read_text(encoding="utf-8"))
            collision_report["pdf"] = str(final_paths[1])
            collision_path.write_text(
                json.dumps(collision_report, indent=2, ensure_ascii=False) + "\n",
                encoding="utf-8",
            )
            staged_to_final.extend(zip(staged_paths, final_paths, strict=True))
            published_paths.extend(final_paths)

        for staged, final in staged_to_final:
            final.parent.mkdir(parents=True, exist_ok=True)
            os.replace(staged, final)
    return published_paths


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        validate_args(args)
        contributions = load_contributions(args.input_dir)
        outputs = export_outputs(
            contributions,
            output_root=args.output_dir,
            dpi=args.dpi,
        )
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    for model in contributions:
        contribution_text = ", ".join(
            f"{group}={model.percentages[group]:.2f}%" for group in GROUP_ORDER
        )
        print(
            f"[{model.spec.response_id}] features={model.feature_count} "
            f"total={model.total_shap:.6f}; {contribution_text}"
        )
    print("Outputs:")
    for output in outputs:
        print(f"  {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
