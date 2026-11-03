#!/usr/bin/env python3
"""Remove selected tracking columns from the stage 13 CSV files.

This script rewrites each CSV in place. It first validates every CSV header,
then prepares temporary replacements and atomically replaces the originals.
No backup copies are retained.

Run with: python3 preprocess_stage_13.py
"""

import csv
import os
import tempfile
from pathlib import Path


DATA_DIR = Path("/sciclone/home/cabare01/data10/merged_data/stage_13")
COLUMNS_TO_REMOVE = {
    "Heading",
    "Speed",
    "Acceleration",
    "Accel.T",
    "Accel.N",
    "PathLength",
    "PathSpeed",
    "LineLength",
    "LineSpeed",
    "RefLineLength",
    "Perimeter",
    "Width",
    "MaxFeret",
    "MinFeret",
    "Circularity",
    "Elongation",
}


def read_header(path: Path) -> list[str]:
    with path.open("r", newline="", encoding="utf-8-sig") as source:
        try:
            return next(csv.reader(source))
        except StopIteration as error:
            raise ValueError(f"CSV is empty: {path}") from error


def prepare_replacement(path: Path, header: list[str]) -> Path:
    kept_indexes = [i for i, name in enumerate(header) if name not in COLUMNS_TO_REMOVE]
    temp_path = None
    try:
        with path.open("r", newline="", encoding="utf-8-sig") as source:
            reader = csv.reader(source)
            next(reader)  # Header was already validated.
            with tempfile.NamedTemporaryFile(
                mode="w",
                newline="",
                encoding="utf-8",
                dir=path.parent,
                prefix=f".{path.name}.",
                suffix=".tmp",
                delete=False,
            ) as destination:
                temp_path = Path(destination.name)
                writer = csv.writer(destination)
                writer.writerow([header[i] for i in kept_indexes])
                for line_number, row in enumerate(reader, start=2):
                    if len(row) != len(header):
                        raise ValueError(
                            f"{path}:{line_number}: expected {len(header)} fields, got {len(row)}"
                        )
                    writer.writerow([row[i] for i in kept_indexes])
        os.chmod(temp_path, path.stat().st_mode)
        return temp_path
    except Exception:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)
        raise


def main() -> None:
    files = sorted(DATA_DIR.glob("*.csv"))
    if not files:
        raise FileNotFoundError(f"No CSV files found in {DATA_DIR}")

    # Check all headers before creating any replacement files.
    headers = {path: read_header(path) for path in files}
    for path, header in headers.items():
        missing = COLUMNS_TO_REMOVE.difference(header)
        if missing:
            raise ValueError(f"{path.name} is missing columns to remove: {sorted(missing)}")

    replacements = {}
    try:
        # Prepare every rewritten CSV before replacing any original.
        for path in files:
            replacements[path] = prepare_replacement(path, headers[path])
        for path, temp_path in replacements.items():
            os.replace(temp_path, path)
    finally:
        for temp_path in replacements.values():
            temp_path.unlink(missing_ok=True)

    print(f"Updated {len(files)} CSV files in {DATA_DIR}")
    print("Removed columns: " + ", ".join(sorted(COLUMNS_TO_REMOVE)))


if __name__ == "__main__":
    main()
