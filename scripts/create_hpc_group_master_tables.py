#!/usr/bin/env python3
"""Create cell-wise (CW) and embryo-wise (EW) spectral tables for each HPC group.

Only files at <group>/<embryo>/tables/singlecell_results_final.csv are read.
Input files are never modified.  Run non-interactively on the HPC with:

    python3 create_hpc_group_master_tables.py

Reads the v3 spectral outputs (data10/HPC_spectral_v3) by default and writes
data/group_tables. CALCIUM_TABLE_VERSION ("v2", the default) names the column
layout with the extra QC columns, which v2 and v3 runs share. Override paths
with CALCIUM_SPECTRAL_ROOT and CALCIUM_TABLE_DIR, e.g. to rebuild v2 tables:
    CALCIUM_SPECTRAL_ROOT=/sciclone/home/cabare01/data10/HPC_spectral_v2 \
    CALCIUM_TABLE_DIR=/tmp/v2_tables python3 create_hpc_group_master_tables.py
(The v1 run, data10/HPC_spectral, was deleted on 2026-11-05.)
"""

import csv
import math
import os
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import embryo_ref  # noqa: E402


INPUT_ROOT = Path(os.environ.get("CALCIUM_SPECTRAL_ROOT",
                                 "/sciclone/home/cabare01/data10/HPC_spectral_v3"))
OUTPUT_DIR = Path(os.environ.get("CALCIUM_TABLE_DIR",
                                 "/sciclone/home/cabare01/CALCIUM/data/group_tables"))
TABLE_VERSION = os.environ.get("CALCIUM_TABLE_VERSION", "v2")

GROUPS = {
    "stage_13": ("S13", 13, "intact", 46),
    "stage_13_devitellinated": ("S13_devitellinated", 13, "devitellinated", 46),
    "stage_14": ("S14", 14, "intact", 46),
    "stage_14_devitellinated": ("S14_devitellinated", 14, "devitellinated", 41),
}

FEATURES = (
    "DominantFreqHz", "SpectralCentroid", "SpectralEntropy", "SpectralSpread",
    "TotalBandPower", "PeakPowerFraction",
)
INPUT_COLUMNS = (
    "CellID", "RunID", "SourceFile", "DominantFreqHz", "SpectralCentroid",
    "SpectralSpread", "SpectralEntropy", "TotalBandPower", "PeakPowerFraction",
    "Duration_min", "MeanX", "MeanY",
)
if TABLE_VERSION == "v2":
    INPUT_COLUMNS += (
        "UsedDuration_min", "FrameInterval_s", "FramesObserved", "FramesUsed",
        "NSegments", "NBridgedGaps", "NSplitGaps",
    )
CELL_COLUMNS = (
    "Group", "Stage", "MembraneStatus", "EmbryoID", "InputFile", "GlobalID",
    *INPUT_COLUMNS, "EmbryoCellCount", "GroupCellCount",
)
EMBRYO_COLUMNS = [
    "Group", "Stage", "MembraneStatus", "EmbryoID", "InputFile",
    "EmbryoCellCount", "GroupEmbryoCount", "GroupCellCount",
]
for feature in FEATURES:
    EMBRYO_COLUMNS.extend((
        f"{feature}_n_valid", f"{feature}_mean", f"{feature}_sd",
        f"{feature}_median", f"{feature}_iqr",
    ))


def numeric_or_none(value, path, column):
    value = value.strip()
    if value.lower() in {"", "na", "nan", "null"}:
        return None
    try:
        number = float(value)
    except ValueError as error:
        raise ValueError(f"Non-numeric {column} in {path}: {value!r}") from error
    return number if math.isfinite(number) else None


def percentile(sorted_values, proportion):
    position = (len(sorted_values) - 1) * proportion
    lower, upper = math.floor(position), math.ceil(position)
    if lower == upper:
        return sorted_values[lower]
    return sorted_values[lower] + (sorted_values[upper] - sorted_values[lower]) * (position - lower)


def summarize(values):
    if not values:
        return (0, None, None, None, None)
    values.sort()
    return (
        len(values),
        statistics.fmean(values),
        statistics.stdev(values) if len(values) > 1 else None,
        statistics.median(values),
        percentile(values, 0.75) - percentile(values, 0.25),
    )


def group_files(folder):
    group_dir = INPUT_ROOT / folder
    if not group_dir.is_dir():
        raise FileNotFoundError(f"Missing group directory: {group_dir}")
    return sorted(group_dir.glob("*/tables/singlecell_results_final.csv"))


def validate_header(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        header = next(csv.reader(handle), None)
    if header != list(INPUT_COLUMNS):
        raise ValueError(f"Unexpected columns in {path}: {header}")


def read_group(folder, group, stage, membrane_status):
    """Read a group once for its embryo summaries and validation counts."""
    records = []
    embryo_ids = set()
    group_cell_count = 0
    for path in group_files(folder):
        validate_header(path)
        embryo_id = path.parent.parent.name
        if embryo_id in embryo_ids:
            raise ValueError(f"Duplicate embryo directory name in {folder}: {embryo_id}")
        embryo_ids.add(embryo_id)
        cell_ids = set()
        values = {feature: [] for feature in FEATURES}
        cell_count = 0
        with path.open(newline="", encoding="utf-8-sig") as handle:
            for row in csv.DictReader(handle):
                cell_id = row["CellID"]
                if cell_id in cell_ids:
                    raise ValueError(f"Duplicate CellID in {embryo_id}: {cell_id}")
                cell_ids.add(cell_id)
                cell_count += 1
                for feature in FEATURES:
                    value = numeric_or_none(row[feature], path, feature)
                    if value is not None:
                        values[feature].append(value)
        records.append({
            "Group": group, "Stage": stage, "MembraneStatus": membrane_status,
            "EmbryoID": embryo_id, "InputFile": str(path.relative_to(INPUT_ROOT)),
            "path": path, "EmbryoCellCount": cell_count,
            "summaries": {feature: summarize(values[feature]) for feature in FEATURES},
        })
        group_cell_count += cell_count
    return records, group_cell_count


def write_cell_table(records, group_cell_count, output_path):
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CELL_COLUMNS)
        writer.writeheader()
        for record in records:
            with record["path"].open(newline="", encoding="utf-8-sig") as input_handle:
                for row in csv.DictReader(input_handle):
                    output = {column: row[column] for column in INPUT_COLUMNS}
                    output.update({
                        "Group": record["Group"], "Stage": record["Stage"],
                        "MembraneStatus": record["MembraneStatus"],
                        "EmbryoID": record["EmbryoID"], "InputFile": record["InputFile"],
                        "GlobalID": f"{record['EmbryoID']}::{row['CellID']}",
                        "EmbryoCellCount": record["EmbryoCellCount"],
                        "GroupCellCount": group_cell_count,
                    })
                    writer.writerow(output)


def write_embryo_table(records, group_cell_count, output_path):
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=EMBRYO_COLUMNS)
        writer.writeheader()
        for record in records:
            output = {key: record[key] for key in (
                "Group", "Stage", "MembraneStatus", "EmbryoID", "InputFile", "EmbryoCellCount"
            )}
            output["GroupEmbryoCount"] = len(records)
            output["GroupCellCount"] = group_cell_count
            for feature, summary in record["summaries"].items():
                for name, value in zip(("n_valid", "mean", "sd", "median", "iqr"), summary):
                    output[f"{feature}_{name}"] = value
            writer.writerow(output)


def main():
    if OUTPUT_DIR.exists() and any(OUTPUT_DIR.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty table directory: {OUTPUT_DIR}")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    total_embryos = total_cells = 0
    for folder, (group, stage, membrane_status, expected_embryos) in GROUPS.items():
        records, group_cell_count = read_group(folder, group, stage, membrane_status)
        # The reference table (data10/Embryo_Ref_Table.csv) is the definitive grouping.
        wanted = embryo_ref.group_keys(folder)
        found = {embryo_ref.embryo_key(r["EmbryoID"]): r["EmbryoID"] for r in records}
        wrong = sorted(v for k, v in found.items() if k not in wanted)
        if wrong:
            raise ValueError(f"{folder}: embryos not assigned to this group in the reference "
                             f"table: {wrong}")
        missing = sorted(wanted - set(found))
        if missing:
            print(f"{folder}: WARNING reference embryos without spectral output: {missing}")
        write_cell_table(records, group_cell_count, OUTPUT_DIR / f"{folder}_CW.csv")
        write_embryo_table(records, group_cell_count, OUTPUT_DIR / f"{folder}_EW.csv")
        total_embryos += len(records)
        total_cells += group_cell_count
        status = "OK" if len(records) == expected_embryos else f"WARNING: expected {expected_embryos}"
        print(f"{folder}: {len(records)} embryos, {group_cell_count} cells ({status})")
    print(f"Wrote 8 tables to {OUTPUT_DIR}")
    print(f"Total: {total_embryos} embryos, {total_cells} cells")


if __name__ == "__main__":
    main()
