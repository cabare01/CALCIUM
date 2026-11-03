#!/usr/bin/env python3
"""Definitive embryo -> group assignments from data10/Embryo_Ref_Table.csv.

The reference table has one column per group (Stage 13, Stage_13_dev,
Stage_14, Stage_14_dev) listing short IDs such as ``180118SP_E7``. File and
boundary names vary (``180118SPE7_st_13_Q``, ``191020RNSP_E6_st_14``,
``…_st_13_devitelined…``), so everything is matched on a normalised key:
date + initials + embryo number, e.g. ``180118SP_E7``. Stage or membrane
tokens in file names are ignored; the reference table is authoritative.
"""
from __future__ import annotations

import csv
import os
import re
from pathlib import Path

REF_TABLE = Path(os.environ.get("CALCIUM_REF_TABLE",
                                "/sciclone/home/cabare01/data10/Embryo_Ref_Table.csv"))
# reference-table column -> (folder/table name, mapping group label, stage, membrane)
GROUP_COLUMNS = {
    "Stage 13": ("stage_13", "S13", 13, "intact"),
    "Stage_13_dev": ("stage_13_devitellinated", "S13_devitellinated", 13, "devitellinated"),
    "Stage_14": ("stage_14", "S14", 14, "intact"),
    "Stage_14_dev": ("stage_14_devitellinated", "S14_devitellinated", 14, "devitellinated"),
}
FOLDER_TO_COLUMN = {v[0]: k for k, v in GROUP_COLUMNS.items()}
LABEL_TO_COLUMN = {v[1]: k for k, v in GROUP_COLUMNS.items()}
_KEY = re.compile(r"^(\d{6})([A-Za-z]+)_?(E\d+)", re.IGNORECASE)


def embryo_key(name: str) -> str | None:
    """Normalised embryo key from any file/embryo name, or None if unparseable."""
    match = _KEY.match(Path(str(name)).name.strip())
    return f"{match[1]}{match[2].upper()}_{match[3].upper()}" if match else None


def load_reference(path: Path = REF_TABLE) -> dict[str, str]:
    """Return {embryo key: reference column}. Raises on duplicates or bad IDs."""
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.reader(handle))
    header = [h.strip() for h in rows[0]]
    unknown = set(header) - set(GROUP_COLUMNS)
    if unknown:
        raise ValueError(f"Unexpected reference-table columns: {sorted(unknown)}")
    assignment: dict[str, str] = {}
    for row in rows[1:]:
        for column, value in zip(header, row):
            if not value.strip():
                continue
            key = embryo_key(value)
            if key is None:
                raise ValueError(f"Cannot parse embryo ID {value!r} in {path}")
            if key in assignment:
                raise ValueError(f"{value} listed in both {assignment[key]} and {column}")
            assignment[key] = column
    return assignment


def group_keys(group: str, path: Path = REF_TABLE) -> set[str]:
    """Embryo keys for a group given as folder name, mapping label or ref column."""
    column = FOLDER_TO_COLUMN.get(group) or LABEL_TO_COLUMN.get(group) or group
    if column not in GROUP_COLUMNS:
        raise ValueError(f"Unknown group {group!r}")
    return {k for k, c in load_reference(path).items() if c == column}
