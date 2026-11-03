#!/usr/bin/env python3
"""Map cell-level spectral results to the legacy SuperEmbryo CNPU grid.

This is a Python implementation of the supplied MATLAB Mapping to Template
workflow. It reads a cell-level spectral master table and per-embryo
*-Boundary_SuperEmbryo.csv files, joins cells by CellID/Calcium_ID, assigns
9-pixel CNPUs at origin (120, 1), and creates embryo-balanced group maps.
It does not alter its input files.

The script requires matplotlib for PNG/PDF output.
"""

import argparse
import csv
import math
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import embryo_ref  # noqa: E402

DEFAULT_EXCLUSIONS = Path(__file__).with_name("mapping_exclusions.csv")

FEATURES = (
    "DominantFreqHz", "SpectralCentroid", "SpectralSpread",
    "SpectralEntropy", "TotalBandPower", "PeakPowerFraction",
)
GROUPS = ("S13", "S13_devitellinated", "S14", "S14_devitellinated")
EXPECTED_BY_GROUP = {
    "S13": 46,
    "S13_devitellinated": 46,
    "S14": 46,
    "S14_devitellinated": 41,
}
# Biological Stage 13 groups are overlaid with the Stage 14 Sox2 template;
# biological Stage 14 groups are overlaid with the Stage 15 Sox2 template.
SOX2_TEMPLATE_FOR_GROUP = {
    "S13": "stage14",
    "S13_devitellinated": "stage14",
    "S14": "stage15",
    "S14_devitellinated": "stage15",
}
CNPU_SIZE = 9
ORIGIN_X, ORIGIN_Y = 120, 1
X_RANGE, Y_RANGE = (120, 1122), (1, 1190)


def numeric(value):
    try:
        x = float(value)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def read_table(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames:
            raise ValueError(f"Empty CSV: {path}")
        return reader.fieldnames, list(reader)


def write_table(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def canonical_id(value):
    # Legacy duplicate source convention documented in the MATLAB workflow.
    return re.sub(r"_st13_devitelined_$", "", value, flags=re.IGNORECASE)


def resolve_embryos(master_rows, boundary_root):
    by_id = defaultdict(list)
    for row in master_rows:
        embryo = row["EmbryoID"].strip()
        by_id[canonical_id(embryo)].append(embryo)
    selected = {}
    for biological_id, source_ids in by_id.items():
        candidates = sorted(set(source_ids))
        if len(candidates) == 1:
            selected[biological_id] = candidates[0]
            continue
        dev = [x for x in candidates if re.search(
            r"_st13_devitelined_$", x, re.IGNORECASE)]
        if len(dev) != 1:
            raise ValueError(
                f"Ambiguous mapped source files for {biological_id}: {candidates}")
        selected[biological_id] = dev[0]

    boundary_files = list(Path(boundary_root).rglob("*-Boundary_SuperEmbryo.csv"))
    by_name, by_canonical_name = defaultdict(list), defaultdict(list)
    for path in boundary_files:
        stem = path.name.removesuffix("-Boundary_SuperEmbryo.csv")
        by_name[stem].append(path)
        by_canonical_name[canonical_id(stem)].append(path)
    resolved = {}
    missing = []
    for biological_id, source_id in selected.items():
        # Prefer the exact source-specific boundary when both legacy variants
        # exist; otherwise accept a canonical match or an unambiguous boundary
        # stem at the beginning of the spectral EmbryoID. The latter supports
        # IDs such as 180118SP_E7_st_13_Q paired with
        # 180118SP_E7-Boundary_SuperEmbryo.csv.
        choices = by_name.get(source_id, [])
        if not choices:
            choices = by_canonical_name.get(canonical_id(source_id), [])
        if not choices:
            prefix_matches = [
                (len(stem), path)
                for stem, paths in by_name.items()
                if source_id.startswith(stem + "_")
                for path in paths
            ]
            if prefix_matches:
                longest = max(length for length, _ in prefix_matches)
                choices = [path for length, path in prefix_matches
                           if length == longest]
        if len(choices) != 1:
            if not choices:
                missing.append(source_id)
                continue
            raise ValueError(f"Multiple boundary files for {source_id}: {choices}")
        resolved[source_id] = choices[0]
    return selected, resolved, missing


def load_boundary(path):
    fields, rows = read_table(path)
    id_field = next((x for x in ("Calcium_ID", "CellID") if x in fields), None)
    if id_field is None or not {"Centroid_1", "Centroid_2"} <= set(fields):
        raise ValueError(f"{path} must contain Calcium_ID/CellID and Centroid_1/2")
    found = defaultdict(list)
    for row in rows:
        cell_id = row[id_field].strip()
        x, y = numeric(row["Centroid_1"]), numeric(row["Centroid_2"])
        if x is None or y is None:
            continue
        found[cell_id].append((x, y))
    duplicates = sorted(cell_id for cell_id, coords in found.items()
                        if len(coords) > 1)
    return found, duplicates


def map_cells(master_rows, boundary_root, expected_count):
    selected, boundary_paths, missing_boundary_ids = resolve_embryos(
        master_rows, boundary_root)
    actual_count = len(selected)
    if actual_count > expected_count:
        raise ValueError(
            f"Found {actual_count} biological embryos, above the maximum "
            f"of {expected_count} for this group.")
    if actual_count < expected_count:
        print(
            f"QC warning: found {actual_count} biological embryos; "
            f"the expected maximum is {expected_count}. Proceeding with those found.")
    source_groups = {}
    for row in master_rows:
        source_groups[row["EmbryoID"].strip()] = row["Group"]
    if missing_boundary_ids:
        print(
            f"QC warning: {len(missing_boundary_ids)} of {actual_count} embryos "
            "have no matching boundary CSV; skipping them in spatial mapping.",
            file=sys.stderr, flush=True)
        for source_id in missing_boundary_ids:
            print(f"  Group={source_groups.get(source_id, 'unknown')} "
                  f"EmbryoID={source_id}: no matching boundary CSV",
                  file=sys.stderr, flush=True)
    if not boundary_paths:
        raise FileNotFoundError(
            "No boundary CSVs matched any embryo in this group; cannot map the group.")

    mapped, missing = [], defaultdict(int)
    selected_sources = set(selected.values())
    group_embryos = defaultdict(set)
    for row in master_rows:
        source_id = row["EmbryoID"].strip()
        if source_id in selected_sources:
            group_embryos[row["Group"]].add(canonical_id(source_id))
    for source_id in selected_sources:
        embryo_groups = {r["Group"] for r in master_rows
                         if r["EmbryoID"].strip() == source_id}
        if len(embryo_groups) != 1:
            raise ValueError(f"Embryo {source_id} has inconsistent group labels: {embryo_groups}")
    for biological_id, source_id in selected.items():
        if source_id not in boundary_paths:
            continue
        boundary, duplicate_ids = load_boundary(boundary_paths[source_id])
        if duplicate_ids:
            print(
                f"QC warning: {len(duplicate_ids)} repeated cell IDs in "
                f"{boundary_paths[source_id]}; matching each spectral row to "
                "all corresponding centroids to mirror the MATLAB innerjoin.",
                file=sys.stderr, flush=True)
            print("  Repeated IDs: " + ", ".join(duplicate_ids),
                  file=sys.stderr, flush=True)
        embryo_rows = [r for r in master_rows
                       if r["EmbryoID"].strip() == source_id]
        for row in embryo_rows:
            cell_id = row["CellID"].strip()
            coordinate_matches = boundary.get(cell_id)
            if not coordinate_matches:
                missing[source_id] += 1
                continue
            # MATLAB's innerjoin produces one merged row for every matching
            # boundary row. Preserve that behavior when a Calcium_ID repeats.
            for x, y in coordinate_matches:
                out = dict(row)
                out["BiologicalEmbryoID"] = biological_id
                out["Centroid_1"], out["Centroid_2"] = x, y
                out["CNPU_X"] = math.floor((x - ORIGIN_X) / CNPU_SIZE)
                out["CNPU_Y"] = math.floor((y - ORIGIN_Y) / CNPU_SIZE)
                out["CNPU_ID"] = f"{out['CNPU_X']}_{out['CNPU_Y']}"
                out["CNPU_CenterX"] = ORIGIN_X + out["CNPU_X"] * CNPU_SIZE + CNPU_SIZE / 2
                out["CNPU_CenterY"] = ORIGIN_Y + out["CNPU_Y"] * CNPU_SIZE + CNPU_SIZE / 2
                mapped.append(out)
    if missing:
        print(f"QC: {sum(missing.values())} spectral cells lacked boundary matches "
              f"across {len(missing)} embryos.")
    max_x = math.ceil((X_RANGE[1] - ORIGIN_X) / CNPU_SIZE) - 1
    max_y = math.ceil((Y_RANGE[1] - ORIGIN_Y) / CNPU_SIZE) - 1
    outside = sum(not (0 <= r["CNPU_X"] <= max_x and 0 <= r["CNPU_Y"] <= max_y)
                  for r in mapped)
    if outside:
        print(f"QC: {outside} mapped cells fall outside the complete legacy grid; "
              "they remain in the mapped-cell CSV.")
    return mapped, group_embryos


def complete_grid():
    max_x = math.ceil((X_RANGE[1] - ORIGIN_X) / CNPU_SIZE) - 1
    max_y = math.ceil((Y_RANGE[1] - ORIGIN_Y) / CNPU_SIZE) - 1
    grid = {}
    for ix in range(max_x + 1):
        for iy in range(max_y + 1):
            key = f"{ix}_{iy}"
            grid[key] = {
                "CNPU_X": ix, "CNPU_Y": iy, "CNPU_ID": key,
                "CNPU_StartX": ORIGIN_X + ix * CNPU_SIZE,
                "CNPU_EndX": ORIGIN_X + (ix + 1) * CNPU_SIZE,
                "CNPU_StartY": ORIGIN_Y + iy * CNPU_SIZE,
                "CNPU_EndY": ORIGIN_Y + (iy + 1) * CNPU_SIZE,
                "CNPU_CenterX": ORIGIN_X + ix * CNPU_SIZE + CNPU_SIZE / 2,
                "CNPU_CenterY": ORIGIN_Y + iy * CNPU_SIZE + CNPU_SIZE / 2,
            }
    return grid


def aggregate(mapped, grid, group_embryos, groups):
    embryo_bins = defaultdict(list)
    pooled_bins = defaultdict(list)  # paper method: all cells from all embryos
    for row in mapped:
        group, embryo, cnpu = row["Group"], row["BiologicalEmbryoID"], row["CNPU_ID"]
        for feature in FEATURES:
            value = numeric(row.get(feature))
            if value is not None:
                embryo_bins[group, embryo, cnpu, feature].append(value)
                pooled_bins[group, cnpu, feature].append(value)

    embryo_rows, group_rows = [], []
    per_group = defaultdict(list)
    for (group, embryo, cnpu, feature), values in embryo_bins.items():
        med = statistics.median(values)
        record = {"Group": group, "EmbryoID": embryo, "CNPU_ID": cnpu,
                  "Feature": feature, "EmbryoMedian": med, "NCells": len(values)}
        embryo_rows.append(record)
        per_group[group, cnpu, feature].append((embryo, med, len(values)))

    for group in groups:
        n_group = len(group_embryos[group])
        for feature in FEATURES:
            for cnpu, cell in grid.items():
                contributors = per_group.get((group, cnpu, feature), [])
                vals = [v for _, v, _ in contributors]
                if vals:
                    median = statistics.median(vals)
                    sd = statistics.stdev(vals) if len(vals) > 1 else None
                    n_embryos = len(vals)
                    total_cells = sum(n for _, _, n in contributors)
                else:
                    median = sd = None
                    n_embryos = total_cells = 0
                group_rows.append({
                    "Group": group, "Feature": feature, **cell,
                    "GroupMedian": median, "GroupSD": sd,
                    "NEmbryos": n_embryos, "NGroupEmbryos": n_group,
                    "CoveragePct": 100 * n_embryos / n_group if n_group else None,
                    "NCells": total_cells,
                    # Composite-neural-plate method of the PDT paper (Methods 2.3):
                    # mean over all cells from all embryos whose centroid is in the CNPU.
                    "PooledMean": (statistics.fmean(pooled_bins[group, cnpu, feature])
                                   if pooled_bins.get((group, cnpu, feature)) else None),
                })

    # Legacy 7x7 (49-CNPU) smoother over observed CNPU values, then z-scoring
    # of the observed smoothed values within each group/feature. It runs on
    # (a) the embryo-balanced group medians -> SmoothedMedian/SmoothedZ and
    # (b) the paper's pooled cell means  -> PooledSmoothed/PooledSmoothedZ.
    lookup = {(r["Group"], r["Feature"], r["CNPU_ID"]): r for r in group_rows}
    for group in groups:
        for feature in FEATURES:
            rows_gf = [r for r in group_rows if r["Group"] == group and r["Feature"] == feature]
            for source, smooth_key, z_key in (("GroupMedian", "SmoothedMedian", "SmoothedZ"),
                                              ("PooledMean", "PooledSmoothed", "PooledSmoothedZ")):
                value_by_xy = {(r["CNPU_X"], r["CNPU_Y"]): r[source]
                               for r in rows_gf if r[source] is not None}
                smoothed = {}
                for cnpu, cell in grid.items():
                    x, y = cell["CNPU_X"], cell["CNPU_Y"]
                    # Smooth only CNPUs that contain cells; empty CNPUs stay
                    # empty and are excluded from the z-score mean/SD.
                    if (x, y) not in value_by_xy:
                        smoothed[cnpu] = None
                        continue
                    vals = [value_by_xy[(x + dx, y + dy)]
                            for dx in range(-3, 4) for dy in range(-3, 4)
                            if (x + dx, y + dy) in value_by_xy]
                    smoothed[cnpu] = statistics.fmean(vals) if vals else None
                observed = [v for v in smoothed.values() if v is not None]
                mu = statistics.fmean(observed) if observed else 0
                sd = statistics.stdev(observed) if len(observed) > 1 else 0
                for cnpu in grid:
                    r = lookup[group, feature, cnpu]
                    v = smoothed[cnpu]
                    r[smooth_key] = v
                    r[z_key] = (v - mu) / sd if v is not None and sd else None
    return embryo_rows, group_rows


def load_template(path):
    if not path:
        return None
    fields, rows = read_table(Path(path))
    candidates_x = ("X", "V1", "col", "Centroid_1", "Boundary_1", "Var1")
    candidates_y = ("Y", "V2", "row", "Centroid_2", "Boundary_2", "Var2")
    lower = {f.lower(): f for f in fields}
    xcol = next((lower[x.lower()] for x in candidates_x if x.lower() in lower), None)
    ycol = next((lower[y.lower()] for y in candidates_y if y.lower() in lower), None)
    if (xcol is None or ycol is None) and len(fields) == 2:
        xcol, ycol = fields
    if xcol is None or ycol is None:
        raise ValueError(f"Cannot identify boundary X/Y columns in {path}")
    points = [(numeric(r[xcol]), numeric(r[ycol])) for r in rows]
    points = [(x, y) for x, y in points if x is not None and y is not None]
    if len(points) < 3:
        raise ValueError(f"Fewer than three valid boundary points in {path}")
    if points[0] != points[-1]:
        points.append(points[0])
    return points


def make_plots(mapped, grid, group_rows, output, stage14_template, stage15_template, groups):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import numpy as np

    out = Path(output)
    individual_dir = out / "individual_maps"
    individual_dir.mkdir(parents=True, exist_ok=True)
    lookup = {(r["Group"], r["Feature"], r["CNPU_ID"]): r for r in group_rows}
    # Legacy complete-grid QC: observed and unobserved locations remain visible.
    occupied = {r["CNPU_ID"] for r in mapped}
    fig, ax = plt.subplots(figsize=(10, 8))
    for cid, cell in grid.items():
        ax.scatter(cell["CNPU_CenterX"], cell["CNPU_CenterY"], s=10,
                   c="#4c78a8" if cid in occupied else "#dddddd")
    ax.set(xlim=(X_RANGE[0], X_RANGE[1]), ylim=(Y_RANGE[0], Y_RANGE[1]),
           xlabel="SuperEmbryo X", ylabel="SuperEmbryo Y",
           title="Complete legacy CNPU grid; blue = observed")
    ax.invert_yaxis()
    ax.set_aspect("equal")
    fig.savefig(out / "CompleteCNPUGrid_QC.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    boundaries = {
        group: (stage14_template if SOX2_TEMPLATE_FOR_GROUP[group] == "stage14"
                else stage15_template)
        for group in GROUPS
    }
    all_x = [r["Centroid_1"] for r in mapped]
    all_y = [r["Centroid_2"] for r in mapped]
    x_margin = .02 * (max(all_x) - min(all_x)) if all_x else 0
    y_margin = .02 * (max(all_y) - min(all_y)) if all_y else 0
    individual_xlim = ((min(all_x) - x_margin, max(all_x) + x_margin)
                       if all_x else X_RANGE)
    individual_ylim = ((min(all_y) - y_margin, max(all_y) + y_margin)
                       if all_y else Y_RANGE)
    for feature in FEATURES:
        all_values = [numeric(r.get(feature)) for r in mapped]
        all_values = [v for v in all_values if v is not None]
        limits = tuple(np.percentile(all_values, [1, 99])) if all_values else (0, 1)
        if limits[0] == limits[1]:
            limits = (limits[0] - 1e-6, limits[1] + 1e-6)
        # Group maps use shared limits calculated from group-level metrics,
        # so each color has the same meaning in every panel.
        medians = [r["GroupMedian"] for r in group_rows
                   if r["Feature"] == feature and r["GroupMedian"] is not None]
        sds = [r["GroupSD"] for r in group_rows
               if r["Feature"] == feature and r["GroupSD"] is not None]
        median_limits = ((min(medians), max(medians)) if medians else (0, 1))
        sd_limits = ((min(sds), max(sds)) if sds else (0, 1))
        if median_limits[0] == median_limits[1]:
            median_limits = (median_limits[0] - 1e-6, median_limits[1] + 1e-6)
        if sd_limits[0] == sd_limits[1]:
            sd_limits = (sd_limits[0] - 1e-6, sd_limits[1] + 1e-6)
        pdf_path = out / f"{feature}_IndividualEmbryos.pdf"
        with PdfPages(pdf_path) as pdf:
            for group in groups:
                for embryo in sorted({r["BiologicalEmbryoID"] for r in mapped
                                      if r["Group"] == group}):
                    subset = [r for r in mapped if r["BiologicalEmbryoID"] == embryo]
                    fig, ax = plt.subplots(figsize=(8, 6))
                    vals = [numeric(r.get(feature)) for r in subset]
                    sc = ax.scatter([r["Centroid_1"] for r in subset],
                                    [r["Centroid_2"] for r in subset],
                                    c=vals, s=12, cmap="viridis",
                                    vmin=limits[0], vmax=limits[1])
                    ax.set(xlim=individual_xlim, ylim=individual_ylim,
                           xlabel="Template X", ylabel="Template Y",
                           title=f"{embryo} — {feature}")
                    ax.invert_yaxis()
                    ax.set_aspect("equal")
                    fig.colorbar(sc, ax=ax, label=feature)
                    png = individual_dir / f"{embryo}_Mapped_{feature}.png"
                    fig.savefig(png, dpi=300, bbox_inches="tight")
                    pdf.savefig(fig, bbox_inches="tight")
                    plt.close(fig)

        # Median, SD, coverage, and legacy smoothed-z group maps.
        # The paper's CNP maps (Methods 2.3) are the primary output: pooled cell
        # mean per CNPU, 49-CNPU smoothing, z-score, shown on a fixed -1..1 jet
        # scale that saturates beyond +/-1 (paper Figs 2-4).
        metrics = (
            ("PooledSmoothedZ", f"{feature}_GroupSpatialMaps", (-1.0, 1.0),
             "Smoothed CNPU mean z-score (paper CNP method)"),
            ("SmoothedZ", f"{feature}_EmbryoBalanced_SmoothedZ_GroupSpatialMaps", (-1.0, 1.0),
             "Smoothed embryo-balanced median z-score"),
            ("PooledMean", f"{feature}_PopulationMean_GroupSpatialMaps", None, "CNPU mean (all cells)"),
            ("GroupMedian", f"{feature}_PopulationMedian_GroupSpatialMaps", median_limits, "Population median"),
            ("GroupSD", f"{feature}_PopulationSD_GroupSpatialMaps", sd_limits, "Population SD"),
            ("CoveragePct", f"{feature}_SpatialCoverage_GroupSpatialMaps", (0, 100), "Embryo coverage (%)"),
        )
        for metric, filename, clim, label in metrics:
            all_metric = [lookup[g, feature, cid][metric] for g in groups
                          for cid in grid if lookup[g, feature, cid].get(metric) is not None]
            if clim is None:
                clim = ((min(all_metric), max(all_metric)) if all_metric else (0, 1))
                if clim[0] == clim[1]:
                    clim = (clim[0] - 1e-6, clim[1] + 1e-6)
            ncols = 2 if len(groups) > 1 else 1
            nrows = math.ceil(len(groups) / ncols)
            fig, axes = plt.subplots(nrows, ncols,
                                     figsize=(6 * ncols, 5 * nrows),
                                     squeeze=False, constrained_layout=True)
            last = None
            for ax, group in zip(axes.flat, groups):
                values = [lookup[group, feature, cid].get(metric) for cid in grid]
                xs = [grid[cid]["CNPU_CenterX"] for cid in grid]
                ys = [grid[cid]["CNPU_CenterY"] for cid in grid]
                valid = [(x, y, v) for x, y, v in zip(xs, ys, values) if v is not None]
                if valid:
                    last = ax.scatter([p[0] for p in valid], [p[1] for p in valid],
                                      c=[p[2] for p in valid], s=18,
                                      cmap="jet" if metric.endswith("Z") else "viridis",
                                      vmin=clim[0], vmax=clim[1])
                boundary = boundaries[group]
                if boundary:
                    ax.plot([p[0] for p in boundary], [p[1] for p in boundary],
                            color="black", linestyle=":", linewidth=1.5)
                ax.set(xlim=X_RANGE, ylim=Y_RANGE, title=group,
                       xlabel="SuperEmbryo X", ylabel="SuperEmbryo Y")
                ax.invert_yaxis()
                ax.set_aspect("equal")
            for ax in list(axes.flat)[len(groups):]:
                ax.set_visible(False)
            if last is not None:
                fig.colorbar(last, ax=list(axes.flat)[:len(groups)], label=label, shrink=.8,
                             extend="both" if metric.endswith("Z") else "neither")
            fig.suptitle(f"{feature} — {label}")
            fig.savefig(out / f"{filename}.png", dpi=300, bbox_inches="tight")
            plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path,
                        default=Path("spectral_analysis_tables/cell_level_spectral_master.csv"))
    parser.add_argument("--boundary-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("spatial_spectral_results"))
    parser.add_argument("--group", choices=GROUPS, required=True,
                        help="Analyze one experimental group per run.")
    parser.add_argument("--stage14-sox2", type=Path, required=True,
                        help="Stage 14 Sox2 template used for both Stage 13 groups.")
    parser.add_argument("--stage15-sox2", type=Path, required=True,
                        help="Stage 15 Sox2 template used for both Stage 14 groups.")
    parser.add_argument("--expected-embryos", type=int, default=None,
                        help="Override the maximum count for this group; lower counts proceed.")
    parser.add_argument("--exclusions", type=Path, default=DEFAULT_EXCLUSIONS,
                        help="CSV with EmbryoKey,Reason of embryos left out of mapping "
                             "(default: mapping_exclusions.csv next to this script).")
    args = parser.parse_args()
    fields, master = read_table(args.input)
    required = {"EmbryoID", "Group", "CellID", *FEATURES}
    missing = required - set(fields)
    if missing:
        raise ValueError("Input missing columns: " + ", ".join(sorted(missing)))
    bad_groups = sorted({r["Group"] for r in master} - set(GROUPS))
    if bad_groups:
        raise ValueError(f"Unexpected group labels: {bad_groups}; expected {GROUPS}")
    selected_master = [r for r in master if r["Group"] == args.group]
    excluded = {}
    if args.exclusions and args.exclusions.exists():
        _, exclusion_rows = read_table(args.exclusions)
        excluded = {r["EmbryoKey"].strip(): r.get("Reason", "") for r in exclusion_rows}
    dropped = sorted({r["EmbryoID"].strip() for r in selected_master
                      if embryo_ref.embryo_key(r["EmbryoID"].strip()) in excluded})
    if dropped:
        print(f"Excluded from mapping ({args.exclusions.name}):")
        for embryo in dropped:
            print(f"  {embryo}: {excluded[embryo_ref.embryo_key(embryo)]}")
        selected_master = [r for r in selected_master
                           if r["EmbryoID"].strip() not in dropped]
    expected_count = (args.expected_embryos if args.expected_embryos is not None
                      else EXPECTED_BY_GROUP[args.group])
    mapped, group_embryos = map_cells(selected_master, args.boundary_root, expected_count)
    grid = complete_grid()
    embryo_rows, group_rows = aggregate(mapped, grid, group_embryos, [args.group])

    output = args.output_dir / args.group
    output.mkdir(parents=True, exist_ok=True)
    mapped_columns = list(dict.fromkeys([*fields, "BiologicalEmbryoID",
        "Centroid_1", "Centroid_2", "CNPU_X", "CNPU_Y", "CNPU_ID",
        "CNPU_CenterX", "CNPU_CenterY"]))
    write_table(output / "All_MappedCells_with_CNPU.csv", mapped, mapped_columns)
    write_table(output / "Embryo_Cnpu_Spectral_Summary.csv", embryo_rows,
                ["Group", "EmbryoID", "CNPU_ID", "Feature", "EmbryoMedian", "NCells"])
    group_columns = list(group_rows[0]) if group_rows else []
    write_table(output / "Group_Cnpu_Spectral_Summary.csv", group_rows, group_columns)
    stage14_template = load_template(args.stage14_sox2)
    stage15_template = load_template(args.stage15_sox2)
    make_plots(mapped, grid, group_rows, output,
               stage14_template, stage15_template, [args.group])
    print(f"Mapped {len(mapped)} cells from {len(set(r['BiologicalEmbryoID'] for r in mapped))} embryos.")
    print(f"Outputs: {output}")


if __name__ == "__main__":
    main()
