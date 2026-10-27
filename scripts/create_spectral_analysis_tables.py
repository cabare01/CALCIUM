#!/usr/bin/env python3
"""Create cell- and embryo-level tables from grouped spectral-result CSV files.

The input CSV files are read only.  The lone CSV directly under INPUT_ROOT is
intentionally excluded; only files in the four named treatment directories are
included.

Run on the HPC command line, for example:
    python3 create_spectral_analysis_tables.py

Or choose a different output directory:
    python3 create_spectral_analysis_tables.py --output-dir /path/to/output
"""

import argparse
import csv
import math
import statistics
from collections import Counter
from pathlib import Path


INPUT_ROOT = Path("/sciclone/home/cabare01/data10/embryo_spectral_results")
DEFAULT_OUTPUT_DIR = Path("/sciclone/home/cabare01/CALCIUM/scripts/spectral_analysis_tables")

GROUPS = {
    "stage_13": ("S13", 13, "intact"),
    "stage_13_devitellinated": ("S13_devitellinated", 13, "devitellinated"),
    "stage_14": ("S14", 14, "intact"),
    "stage_14_devitellinated": ("S14_devitellinated", 14, "devitellinated"),
}

FEATURES = [
    "DominantFreqHz",
    "SpectralCentroid",
    "SpectralEntropy",
    "SpectralSpread",
    "TotalBandPower",
    "PeakPowerFraction",
]

INPUT_COLUMNS = [
    "CellID", "RunID", "SourceFile", "DominantFreqHz", "SpectralCentroid",
    "SpectralSpread", "SpectralEntropy", "TotalBandPower", "PeakPowerFraction",
    "Duration_min", "MeanX", "MeanY",
]

CELL_COLUMNS = [
    "Group", "Stage", "MembraneStatus", "EmbryoID", "InputFile", "GlobalID",
    *INPUT_COLUMNS, "EmbryoCellCount", "GroupCellCount",
]

EMBRYO_COLUMNS = [
    "Group", "Stage", "MembraneStatus", "EmbryoID", "InputFile",
    "EmbryoCellCount", "GroupEmbryoCount", "GroupCellCount",
]
for feature in FEATURES:
    EMBRYO_COLUMNS.extend([
        f"{feature}_n_valid", f"{feature}_mean", f"{feature}_sd",
        f"{feature}_median", f"{feature}_iqr",
    ])


def embryo_id_from_filename(path):
    suffix = "_spectral_results.csv"
    if not path.name.endswith(suffix):
        raise ValueError(f"Unexpected result filename: {path.name}")
    return path.name[:-len(suffix)]


def numeric_or_none(value, file_path, column):
    """Convert blank/NA/NaN to None; stop on other non-numeric feature values."""
    value = value.strip()
    if value.lower() in {"", "na", "nan", "null"}:
        return None
    try:
        number = float(value)
    except ValueError as error:
        raise ValueError(f"Non-numeric {column!r} in {file_path}: {value!r}") from error
    return None if math.isnan(number) else number


def percentile(values, proportion):
    """Linearly interpolated percentile for an already sorted numeric list."""
    position = (len(values) - 1) * proportion
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return values[lower]
    return values[lower] + (values[upper] - values[lower]) * (position - lower)


def feature_summary(values):
    if not values:
        return {"n_valid": 0, "mean": None, "sd": None, "median": None, "iqr": None}
    values.sort()
    return {
        "n_valid": len(values),
        "mean": statistics.fmean(values),
        "sd": statistics.stdev(values) if len(values) > 1 else None,
        "median": statistics.median(values),
        "iqr": percentile(values, 0.75) - percentile(values, 0.25),
    }


def grouped_files(input_root):
    files = []
    for folder, metadata in GROUPS.items():
        group_dir = input_root / folder
        if not group_dir.is_dir():
            raise FileNotFoundError(f"Missing required group directory: {group_dir}")
        for path in sorted(group_dir.glob("*_spectral_results.csv")):
            files.append((path, *metadata))
    return files


def read_header(path):
    with path.open(newline="") as handle:
        header = next(csv.reader(handle), None)
    missing_columns = [column for column in INPUT_COLUMNS if column not in (header or [])]
    if missing_columns:
        raise ValueError(f"Missing required columns in {path}: {missing_columns}")


def build_embryo_records(files):
    """First pass: validate input and calculate counts and embryo summaries."""
    records = []
    embryo_ids = set()
    group_cells = Counter()
    group_embryos = Counter()

    for path, group, stage, membrane_status in files:
        read_header(path)
        embryo_id = embryo_id_from_filename(path)
        if embryo_id in embryo_ids:
            raise ValueError(f"Duplicate EmbryoID derived from filenames: {embryo_id}")
        embryo_ids.add(embryo_id)

        cell_ids = set()
        values_by_feature = {feature: [] for feature in FEATURES}
        cell_count = 0
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                cell_count += 1
                cell_id = row["CellID"]
                if cell_id in cell_ids:
                    raise ValueError(f"Duplicate CellID within embryo {embryo_id}: {cell_id}")
                cell_ids.add(cell_id)
                for feature in FEATURES:
                    value = numeric_or_none(row[feature], path, feature)
                    if value is not None:
                        values_by_feature[feature].append(value)

        summaries = {feature: feature_summary(values) for feature, values in values_by_feature.items()}
        records.append({
            "path": path, "Group": group, "Stage": stage,
            "MembraneStatus": membrane_status, "EmbryoID": embryo_id,
            "InputFile": path.name, "EmbryoCellCount": cell_count,
            "summaries": summaries,
        })
        group_cells[group] += cell_count
        group_embryos[group] += 1

    return records, group_cells, group_embryos


def write_cell_master(records, group_cells, output_path):
    """Second pass: write original cell fields plus stable analysis identifiers."""
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CELL_COLUMNS)
        writer.writeheader()
        for record in records:
            with record["path"].open(newline="") as input_handle:
                for row in csv.DictReader(input_handle):
                    output_row = {
                        "Group": record["Group"],
                        "Stage": record["Stage"],
                        "MembraneStatus": record["MembraneStatus"],
                        "EmbryoID": record["EmbryoID"],
                        "InputFile": record["InputFile"],
                        "GlobalID": f"{record['EmbryoID']}::{row['CellID']}",
                        "EmbryoCellCount": record["EmbryoCellCount"],
                        "GroupCellCount": group_cells[record["Group"]],
                    }
                    output_row.update({column: row[column] for column in INPUT_COLUMNS})
                    writer.writerow(output_row)


def write_embryo_table(records, group_cells, group_embryos, output_path):
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=EMBRYO_COLUMNS)
        writer.writeheader()
        for record in records:
            output_row = {
                key: record[key] for key in [
                    "Group", "Stage", "MembraneStatus", "EmbryoID", "InputFile", "EmbryoCellCount"
                ]
            }
            output_row["GroupEmbryoCount"] = group_embryos[record["Group"]]
            output_row["GroupCellCount"] = group_cells[record["Group"]]
            for feature, summary in record["summaries"].items():
                for statistic, value in summary.items():
                    output_row[f"{feature}_{statistic}"] = value
            writer.writerow(output_row)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=INPUT_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()

    files = grouped_files(args.input_root)
    expected_counts = {"S13": 46, "S13_devitellinated": 46, "S14": 46, "S14_devitellinated": 41}
    observed_counts = Counter(group for _, group, _, _ in files)
    if dict(observed_counts) != expected_counts:
        raise ValueError(f"Unexpected embryo file counts: {dict(observed_counts)}")

    records, group_cells, group_embryos = build_embryo_records(files)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    cell_output = args.output_dir / "cell_level_spectral_master.csv"
    embryo_output = args.output_dir / "embryo_level_spectral_summary.csv"
    write_cell_master(records, group_cells, cell_output)
    write_embryo_table(records, group_cells, group_embryos, embryo_output)

    print(f"Wrote {len(records)} embryos and {sum(group_cells.values())} cells.")
    for group in ["S13", "S13_devitellinated", "S14", "S14_devitellinated"]:
        print(f"{group}: {group_embryos[group]} embryos, {group_cells[group]} cells")
    print(f"Cell master: {cell_output}")
    print(f"Embryo summary: {embryo_output}")


if __name__ == "__main__":
    main()
