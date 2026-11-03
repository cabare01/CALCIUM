#!/usr/bin/env python3
"""Create a visual QC check for one mapped SuperEmbryo boundary CSV.

Plots mapped cell centroids, a subsample of cell outlines, and one selected
cell outline. This is a read-only input check; it writes one PNG.
"""

import argparse
import csv
import re
from pathlib import Path


def num(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("boundary_csv", type=Path)
    parser.add_argument("--output", type=Path, default=Path("boundary_check.png"))
    parser.add_argument("--cell-row", type=int, default=100,
                        help="1-based CSV data-row number to highlight.")
    parser.add_argument("--outline-stride", type=int, default=10,
                        help="Draw every Nth valid cell outline.")
    args = parser.parse_args()

    with args.boundary_csv.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        fields, rows = reader.fieldnames or [], list(reader)
    if not {"Centroid_1", "Centroid_2"} <= set(fields):
        raise ValueError("CSV needs Centroid_1 and Centroid_2 columns.")
    boundary_cols = [f for f in fields if re.fullmatch(r"Boundary_\d+", f)]
    boundary_cols.sort(key=lambda f: int(f.rsplit("_", 1)[1]))
    if not boundary_cols:
        raise ValueError("No Boundary_1, Boundary_2, ... columns found.")
    if args.cell_row < 1 or args.cell_row > len(rows):
        raise ValueError(f"cell-row must be between 1 and {len(rows)}.")
    if args.outline_stride < 1:
        raise ValueError("outline-stride must be >= 1.")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    valid_centroids = []
    for row in rows:
        x, y = num(row["Centroid_1"]), num(row["Centroid_2"])
        if x is not None and y is not None:
            valid_centroids.append((x, y))
    if not valid_centroids:
        raise ValueError("No valid centroid coordinates found.")

    fig, axes = plt.subplots(1, 2, figsize=(13, 6), constrained_layout=True)
    ax = axes[0]
    ax.scatter([p[0] for p in valid_centroids], [p[1] for p in valid_centroids],
               s=8, color="#377eb8")
    ax.set_title(f"Mapped centroids ({len(valid_centroids)} valid)")
    ax.set(xlabel="Template X", ylabel="Template Y")
    ax.set_aspect("equal")
    ax.invert_yaxis()

    ax = axes[1]
    valid_outlines = 0
    selected_has_outline = False
    for i, row in enumerate(rows):
        packed = [num(row[col]) for col in boundary_cols]
        packed = [v for v in packed if v is not None]
        if not packed or len(packed) % 2:
            continue
        valid_outlines += 1
        if i % args.outline_stride:
            continue
        half = len(packed) // 2
        ys, xs = packed[:half], packed[half:]
        ax.plot(xs, ys, color="#777777", linewidth=.4)

    selected = rows[args.cell_row - 1]
    packed = [num(selected[col]) for col in boundary_cols]
    packed = [v for v in packed if v is not None]
    if packed and len(packed) % 2 == 0:
        selected_has_outline = True
        half = len(packed) // 2
        ys, xs = packed[:half], packed[half:]
        ax.plot(xs, ys, color="red", linewidth=1.5, label="Selected cell")
        cx, cy = num(selected["Centroid_1"]), num(selected["Centroid_2"])
        if cx is not None and cy is not None:
            ax.scatter([cx], [cy], s=35, color="red", zorder=3)
        ax.legend()
    outline_title = f"Cell outlines (valid rows: {valid_outlines})"
    if not selected_has_outline:
        outline_title += "; selected row has no valid outline"
    ax.set_title(outline_title)
    ax.set(xlabel="Template X", ylabel="Template Y")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    fig.suptitle(args.boundary_csv.name)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved boundary check: {args.output}")


if __name__ == "__main__":
    main()
