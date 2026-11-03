#!/usr/bin/env python3
"""Inventory of every reference embryo's input files (read-only).

For each embryo in data10/Embryo_Ref_Table.csv, lists the merged raw CSV
(and whether it is in the matching group folder), the SuperEmbryo boundary
file, and the previous student's ConvergeSpikeCy5 event file. Also lists
files that match no reference embryo. Writes data/embryo_inventory.csv and
prints a summary, including embryos without boundary files.

    python3 scripts/embryo_inventory.py
"""
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import embryo_ref  # noqa: E402

DATA10 = Path(os.environ.get("CALCIUM_DATA10", "/sciclone/home/cabare01/data10"))
RAW_ROOT = DATA10 / "merged_data"
BOUNDARY_ROOT = DATA10 / "boundary_info" / "Boundaries"
OUTPUT = Path(os.environ.get("CALCIUM_INVENTORY",
                             "/sciclone/home/cabare01/CALCIUM/data/embryo_inventory.csv"))


def index(paths):
    found = defaultdict(list)
    for path in paths:
        found[embryo_ref.embryo_key(path.name)].append(path)
    return found


def main():
    reference = embryo_ref.load_reference()
    raw = index(RAW_ROOT.rglob("*.csv"))
    boundary = index(BOUNDARY_ROOT.glob("*-Boundary_SuperEmbryo.csv"))
    spikes = index(list(BOUNDARY_ROOT.glob("*/ConvergeSpikeCy5/*.csv"))
                   + list(BOUNDARY_ROOT.glob("wrong stage*/*.csv")))
    rows = []
    for key, column in sorted(reference.items(), key=lambda kv: (kv[1], kv[0])):
        folder, label = embryo_ref.GROUP_COLUMNS[column][:2]
        raw_files = raw.get(key, [])
        rows.append({
            "Group": label, "RefColumn": column, "EmbryoKey": key,
            "RawFile": "; ".join(p.name for p in raw_files),
            "RawFolder": "; ".join(p.parent.name for p in raw_files),
            "RawInGroupFolder": bool(raw_files) and all(p.parent.name == folder for p in raw_files),
            "BoundaryFile": "; ".join(p.name for p in boundary.get(key, [])),
            "SpikeCy5File": "; ".join(p.name for p in spikes.get(key, [])),
        })
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    def show(title, items):
        print(f"{title}: {len(items)}")
        for item in items:
            print(f"  {item}")

    counts = defaultdict(int)
    for r in rows:
        counts[r["Group"]] += 1
    print("Reference embryos per group:", dict(counts))
    show("Reference embryos WITHOUT a boundary file",
         [f"{r['Group']}: {r['EmbryoKey']}" for r in rows if not r["BoundaryFile"]])
    show("Reference embryos without a merged raw CSV",
         [f"{r['Group']}: {r['EmbryoKey']}" for r in rows if not r["RawFile"]])
    show("Raw CSVs outside their reference group folder",
         [f"{r['Group']}: {r['RawFile']} (in {r['RawFolder']})" for r in rows
          if r["RawFile"] and not r["RawInGroupFolder"]])
    show("Reference embryos without a ConvergeSpikeCy5 (PDT event) file",
         [f"{r['Group']}: {r['EmbryoKey']}" for r in rows if not r["SpikeCy5File"]])
    show("Raw CSVs matching no reference embryo",
         [p.relative_to(RAW_ROOT) for k, v in raw.items() if k not in reference for p in v])
    show("Boundary files matching no reference embryo",
         [p.name for k, v in boundary.items() if k not in reference for p in v])
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
