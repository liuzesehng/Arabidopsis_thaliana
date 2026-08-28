#!/usr/bin/env python3
"""Filter Alt.snp_meth.tsv rows and columns by missingness and support."""

from __future__ import annotations

import csv
import math
from pathlib import Path


SCRIPT_PATH = Path(__file__).absolute()
DEFAULT_INPUT = SCRIPT_PATH.parents[2] / "list" / "RCA" / "Alt.all.snp_meth.tsv"
DEFAULT_OUTPUT = SCRIPT_PATH.parents[2] / "list" / "RCA" / "Alt.all.snp_meth.filtered.tsv"
LAST_FIXED_COLUMN_COUNT = 9
MISSING_VALUES = {"", "NA"}


def is_positive_number(value: str) -> bool:
    if value in MISSING_VALUES:
        return False
    try:
        return float(value) > 0
    except ValueError:
        return False


def row_has_signal(row: list[str]) -> bool:
    return any(value not in MISSING_VALUES for value in row[1:-LAST_FIXED_COLUMN_COUNT])


def main() -> None:
    input_path = DEFAULT_INPUT
    output_path = DEFAULT_OUTPUT

    with input_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        rows = list(reader)

    if len(header) <= LAST_FIXED_COLUMN_COUNT + 1:
        raise ValueError("Input file does not contain enough columns for the requested filtering.")

    filtered_rows = [row for row in rows if row_has_signal(row)]
    removed_row_count = len(rows) - len(filtered_rows)
    kept_row_count = len(filtered_rows)
    threshold = math.ceil(kept_row_count * 0.1) if kept_row_count else 0

    c_columns = [name for name in header if name.startswith("C")]
    s_columns = [name for name in header if name.startswith("s")]
    column_index = {name: idx for idx, name in enumerate(header)}

    kept_c_columns = [
        name
        for name in c_columns
        if sum(is_positive_number(row[column_index[name]]) for row in filtered_rows) >= threshold
    ]
    kept_s_columns = [
        name
        for name in s_columns
        if sum(is_positive_number(row[column_index[name]]) for row in filtered_rows) >= threshold
    ]

    output_columns = [header[0], *kept_c_columns, *kept_s_columns, *header[-LAST_FIXED_COLUMN_COUNT:]]
    output_indices = [column_index[name] for name in output_columns]

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(output_columns)
        for row in filtered_rows:
            writer.writerow([row[idx] for idx in output_indices])

    print(f"Input file: {input_path}")
    print(f"Output file: {output_path}")
    print(f"Original rows: {len(rows)}")
    print(f"Removed empty rows: {removed_row_count}")
    print(f"Kept rows: {kept_row_count}")
    print(f"10% threshold: {threshold}")
    print(f"Original C* columns: {len(c_columns)}")
    print(f"Kept C* columns: {len(kept_c_columns)}")
    print(f"Original s* columns: {len(s_columns)}")
    print(f"Kept s* columns: {len(kept_s_columns)}")


if __name__ == "__main__":
    main()
