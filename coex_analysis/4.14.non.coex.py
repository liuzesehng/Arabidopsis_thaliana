#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(
    "/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana"
)
DEFAULT_INPUT = (
    PROJECT_ROOT
    / "list/RCA/AT2G39730_motif_gene_outputs"
    / "total_TPM_negative_AT2G39730_motif_gene_matches.tsv"
)
DEFAULT_OUTPUT = DEFAULT_INPUT.with_name(
    "total_TPM_negative_AT2G39730_motif_gene_matches.in_window_or_in_motif.tsv"
)
DEFAULT_UNIQUE_OUTPUT = DEFAULT_INPUT.with_name(
    "total_TPM_negative_AT2G39730_motif_gene_matches.in_window_or_in_motif.unique_motif_id.tsv"
)
TRUE_VALUES = {"true", "t", "1", "yes", "y"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Select motif-match rows where in_window or in_motif is True."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Input TSV file.")
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Output TSV file for selected rows.",
    )
    parser.add_argument(
        "--unique-output",
        type=Path,
        default=DEFAULT_UNIQUE_OUTPUT,
        help="Output TSV file for unique values from the third column after filtering.",
    )
    return parser.parse_args()


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin(TRUE_VALUES)


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.input, sep="\t")

    required_columns = {"in_window", "in_motif"}
    missing_columns = required_columns.difference(df.columns)
    if missing_columns:
        missing_text = ", ".join(sorted(missing_columns))
        raise ValueError(f"Missing required column(s): {missing_text}")

    selected = df[as_bool(df["in_window"]) | as_bool(df["in_motif"])].copy()
    third_column = selected.columns[2]
    unique_third_column = (
        selected[[third_column]]
        .dropna()
        .drop_duplicates()
        .sort_values(by=third_column)
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.unique_output.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(args.output, sep="\t", index=False)
    unique_third_column.to_csv(args.unique_output, sep="\t", index=False)

    print(f"Input rows: {len(df)}")
    print(f"Selected rows: {len(selected)}")
    print(f"Unique {third_column} rows: {len(unique_third_column)}")
    print(f"Output: {args.output}")
    print(f"Unique output: {args.unique_output}")


if __name__ == "__main__":
    main()
