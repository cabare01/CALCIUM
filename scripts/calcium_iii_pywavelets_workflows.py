#!/usr/bin/env python3
"""Calcium III workflows for Bora/SLURM using PyWavelets.

Select a workflow with CALCIUM_WORKFLOW and set its input variables as described
in the accompanying log. Every run requires a new CALCIUM_OUTPUT_ROOT.
"""
from __future__ import annotations

import csv
import math
import os
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

FEATURES = ["DominantFreqHz", "SpectralCentroid", "SpectralSpread",
            "SpectralEntropy", "TotalBandPower", "PeakPowerFraction"]
EVENT_PREFIXES = ["ActivePeaksPerHr_", "Total_", "CellWiseActiveIncidence_",
                  "OnlySpikes_", "OnlyRosettes_", "OnlyWaves_"]
CNPU_SIZE, GRID_X0, GRID_Y0 = 9, 120, 1
X_RANGE, Y_RANGE = (120, 1122), (1, 1190)


def env_path(name: str, *, file: bool | None = None) -> Path:
    raw = os.environ.get(name)
    if not raw:
        raise ValueError(f"Required environment variable is unset: {name}")
    path = Path(raw).expanduser()
    if file is True and not path.is_file():
        raise FileNotFoundError(path)
    if file is False and not path.is_dir():
        raise NotADirectoryError(path)
    return path


def safe_name(value: object) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(value))


def require_columns(frame: pd.DataFrame, columns: list[str]) -> None:
    missing = [name for name in columns if name not in frame.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")


def output_root() -> Path:
    root = env_path("CALCIUM_OUTPUT_ROOT")
    if root.exists():
        raise FileExistsError(f"Refusing to reuse existing output path: {root}")
    root.mkdir(parents=True, exist_ok=False)
    return root


def main() -> None:
    workflow = os.environ.get("CALCIUM_WORKFLOW", "").lower()
    handlers = {
        "spectral": run_spectral,
        "correlation": run_correlations,
        "corr_summary": run_correlation_summary,
        "individual_maps": run_individual_maps,
        "combine_maps": run_combine_maps,
        "group_spatial": run_group_spatial,
        "boundary_qc": run_boundary_qc,
    }
    if workflow not in handlers:
        raise ValueError(f"Set CALCIUM_WORKFLOW to one of: {', '.join(handlers)}")
    root = output_root()
    print(f"Workflow: {workflow}\nOutput: {root}", flush=True)
    handlers[workflow](root)
    print(f"Workflow finished: {workflow}", flush=True)


def _seconds(time_values: pd.Series, column_name: str) -> np.ndarray:
    if pd.api.types.is_datetime64_any_dtype(time_values):
        return (time_values - time_values.dropna().iloc[0]).dt.total_seconds().to_numpy(float)
    if pd.api.types.is_timedelta64_dtype(time_values):
        return (time_values - time_values.dropna().iloc[0]).dt.total_seconds().to_numpy(float)
    values = pd.to_numeric(time_values, errors="coerce").to_numpy(float)
    if not np.isfinite(values).any():
        parsed = pd.to_datetime(time_values, errors="coerce")
        if parsed.notna().any():
            return (parsed - parsed.dropna().iloc[0]).dt.total_seconds().to_numpy(float)
        durations = pd.to_timedelta(time_values, errors="coerce")
        if durations.notna().any():
            return (durations - durations.dropna().iloc[0]).dt.total_seconds().to_numpy(float)
    unit = os.environ.get("CALCIUM_TIME_UNIT", "s").lower()
    if "ms" in column_name.lower() or unit == "ms":
        values /= 1000.0
    elif "min" in column_name.lower() or unit == "min":
        values *= 60.0
    elif unit not in {"s", "sec", "seconds"}:
        raise ValueError(f"Unsupported CALCIUM_TIME_UNIT: {unit}")
    return values


def _preprocess_trace(rows: pd.DataFrame) -> dict[str, np.ndarray | float]:
    time_candidates = ["Time", "time", "Timestamp", "timestamp", "t", "Time_s", "Time (s)"]
    time_col = next((name for name in time_candidates if name in rows.columns), None)
    if time_col is None:
        raise ValueError("No recognized time column")
    intensity_col = next((name for name in rows.columns if name.lower() == "meanintensity"), None)
    if intensity_col is None:
        intensity_col = next((name for name in rows.columns
                              if "mean" in name.lower() and "intens" in name.lower()), None)
    if intensity_col is None:
        raise ValueError("No MeanIntensity-like column")
    require_columns(rows, ["PositionX", "PositionY"])
    time = _seconds(rows[time_col], time_col)
    intensity = pd.to_numeric(rows[intensity_col], errors="coerce").to_numpy(float)
    x = pd.to_numeric(rows["PositionX"], errors="coerce").to_numpy(float)
    y = pd.to_numeric(rows["PositionY"], errors="coerce").to_numpy(float)
    valid = np.isfinite(time) & np.isfinite(intensity) & np.isfinite(x) & np.isfinite(y)
    time, intensity, x, y = time[valid], intensity[valid], x[valid], y[valid]
    if len(time) < 2:
        raise ValueError("Fewer than two complete observations")
    if np.any(np.diff(time) < 0):
        raise ValueError("Time values are not monotonic; trace was not reordered")
    keep = np.r_[True, np.diff(time) != 0]
    time, intensity, x, y = time[keep], intensity[keep], x[keep], y[keep]
    dt = np.diff(time)
    dt_median = float(np.median(dt))
    if not np.isfinite(dt_median) or dt_median <= 0:
        raise ValueError("Invalid median sampling interval")
    gap_indices = np.flatnonzero(dt > 2 * dt_median)
    starts = np.r_[0, gap_indices + 1]
    ends = np.r_[gap_indices + 1, len(time)]
    t_parts, f_parts, x_parts, y_parts = [], [], [], []
    for start, end in zip(starts, ends):
        ts, fs, xs, ys = time[start:end], intensity[start:end], x[start:end], y[start:end]
        if len(ts) < 2:
            continue
        n = max(2, int(round((ts[-1] - ts[0]) / dt_median)) + 1)
        if n < 10:
            continue
        grid = np.linspace(ts[0], ts[-1], n)
        t_parts.append(grid)
        f_parts.append(np.interp(grid, ts, fs))
        x_parts.append(np.interp(grid, ts, xs))
        y_parts.append(np.interp(grid, ts, ys))
    if not t_parts:
        raise ValueError("No gap-free segment had at least 10 resampled observations")
    t_uniform = np.concatenate(t_parts)
    fluorescence = np.concatenate(f_parts)
    f0 = max(float(np.percentile(fluorescence, 10)), np.finfo(float).eps)
    return {"time": t_uniform, "dff": (fluorescence - f0) / f0,
            "fs": 1.0 / dt_median, "x": np.concatenate(x_parts),
            "y": np.concatenate(y_parts)}


def _wavelet_features(trace: dict, axis_hz: np.ndarray) -> tuple[dict, np.ndarray]:
    try:
        import pywt
    except ImportError as exc:
        raise RuntimeError("PyWavelets is required; install pywt in the SLURM Python environment") from exc
    wavelet_name = os.environ.get("CALCIUM_PYWT_WAVELET", "cmor1.5-1.0")
    wavelet = pywt.ContinuousWavelet(wavelet_name)
    dt = 1.0 / float(trace["fs"])
    scales = pywt.central_frequency(wavelet) / (axis_hz * dt)
    coefficients, _ = pywt.cwt(np.asarray(trace["dff"]), scales, wavelet,
                                sampling_period=dt, method="fft")
    mean_power = np.maximum(np.mean(np.abs(coefficients) ** 2, axis=1), 0.0)
    total = float(np.sum(mean_power))
    result = dict.fromkeys(FEATURES, np.nan)
    if total > 0 and np.isfinite(total):
        probability = mean_power / (total + np.finfo(float).eps)
        peak_index = int(np.argmax(mean_power))
        centroid = float(np.sum(axis_hz * mean_power) / (total + np.finfo(float).eps))
        result.update({
            "DominantFreqHz": float(axis_hz[peak_index]),
            "SpectralCentroid": centroid,
            "SpectralSpread": float(np.sqrt(np.sum((axis_hz-centroid)**2 * mean_power) /
                                            (total + np.finfo(float).eps))),
            "SpectralEntropy": float(-np.sum(probability * np.log2(probability + np.finfo(float).eps))),
            "TotalBandPower": total,
            "PeakPowerFraction": float(mean_power[peak_index] / (total + np.finfo(float).eps)),
        })
    return result, mean_power


def run_spectral(root: Path) -> None:
    source = env_path("CALCIUM_RAW_TRACE_ROOT", file=False)
    pattern = os.environ.get("CALCIUM_RAW_TRACE_PATTERN", "*.csv")
    files = sorted(source.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No raw trace files matching {pattern} in {source}")
    embryo_names = [safe_name(p.stem.removesuffix(" - Not_self_Merged")) for p in files]
    if len(set(embryo_names)) != len(embryo_names):
        raise ValueError("Input files resolve to duplicate embryo output names")
    axis_hz = np.linspace(0.02, 0.20, 80)
    save_plots = os.environ.get("CALCIUM_MAKE_TRACE_PLOTS", "false").lower() == "true"
    import pywt  # fail before processing when dependency is unavailable
    print(f"PyWavelets {pywt.__version__}; wavelet={os.environ.get('CALCIUM_PYWT_WAVELET','cmor1.5-1.0')}")
    for file_index, path in enumerate(files, 1):
        embryo = embryo_names[file_index - 1]
        print(f"Processing {file_index}/{len(files)}: {path.name}", flush=True)
        table = pd.read_csv(path)
        require_columns(table, ["ID"])
        records, power_rows, cell_ids = [], [], []
        figure_dir = root / embryo / "figures"
        if save_plots:
            figure_dir.mkdir(parents=True)
        for cell_id, group in table.groupby("ID", sort=False, dropna=True):
            try:
                trace = _preprocess_trace(group)
                spectral, mean_power = _wavelet_features(trace, axis_hz)
                run_id = pd.Timestamp.now().strftime("%Y-%m-%d_%H%M%S")
                records.append({"CellID": str(cell_id), "RunID": run_id,
                    "SourceFile": path.name, **spectral,
                    "Duration_min": (trace["time"][-1] - trace["time"][0]) / 60.0,
                    "MeanX": float(np.mean(trace["x"])), "MeanY": float(np.mean(trace["y"]))})
                cell_ids.append(str(cell_id)); power_rows.append(mean_power)
                if save_plots:
                    import matplotlib.pyplot as plt
                    fig, ax = plt.subplots(figsize=(9, 3))
                    ax.plot(trace["time"], trace["dff"], lw=1.2)
                    ax.set(xlabel="Time (s)", ylabel="DeltaF/F", title=f"{path.name} | Cell {cell_id}")
                    fig.tight_layout(); fig.savefig(figure_dir / f"{embryo}_Cell_{safe_name(cell_id)}.png", dpi=200); plt.close(fig)
            except Exception as exc:
                warnings.warn(f"Skipping embryo {embryo} cell {cell_id}: {exc}")
        if not records:
            raise ValueError(f"No usable cell traces in {path}")
        embryo_dir = root / embryo
        (embryo_dir / "tables").mkdir(parents=True)
        pd.DataFrame.from_records(records, columns=["CellID", "RunID", "SourceFile", *FEATURES,
            "Duration_min", "MeanX", "MeanY"]).to_csv(embryo_dir / "tables" / "singlecell_results_final.csv", index=False)
        np.savez_compressed(embryo_dir / "SpectralLibrary.npz", CellID=np.asarray(cell_ids),
            Frequency=axis_hz, Power=np.asarray(power_rows))


def _finite_numeric(values) -> np.ndarray:
    return pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(float)


def run_correlations(root: Path) -> None:
    event = pd.read_csv(env_path("CALCIUM_EVENT_FILE", file=True))
    master = pd.read_csv(env_path("CALCIUM_CELL_MASTER", file=True), dtype={"EmbryoID": str, "CellID": str})
    embryos = pd.read_csv(env_path("CALCIUM_EMBRYO_SUMMARY", file=True), dtype={"EmbryoID": str})
    require_columns(event, ["Name", "ID"])
    require_columns(master, ["EmbryoID", "CellID", "Group", *FEATURES])
    require_columns(embryos, ["EmbryoID", "Group"])
    event["Name"] = event["Name"].astype(str); event["CellID"] = event["ID"].astype(str)
    spectral = master.copy(); spectral["CellID"] = spectral["CellID"].astype(str)
    pdt_names = [c.removeprefix("ActivePeaksPerHr_") for c in event.columns if c.startswith("ActivePeaksPerHr_")]
    if not pdt_names:
        raise ValueError("Event file has no ActivePeaksPerHr_* columns")
    records = []
    for _, meta in embryos.iterrows():
        embryo_id = str(meta.EmbryoID)
        event_names = event.Name.dropna().astype(str).unique()
        hits = [name for name in event_names if name.replace("_", "").lower() == embryo_id.replace("_", "").lower()]
        if len(hits) != 1:
            warnings.warn(f"No unique event-data match for {embryo_id}; skipping")
            continue
        merged = event.loc[event.Name.astype(str) == hits[0]].merge(
            spectral.loc[spectral.EmbryoID.astype(str) == embryo_id], on="CellID", how="inner",
            suffixes=("_event", "_spectral"), validate="many_to_many")
        if len(merged) < 3:
            continue
        group_name = str(meta.Group)
        for pdt in pdt_names:
            variables = [prefix + pdt for prefix in EVENT_PREFIXES if prefix + pdt in merged.columns]
            pairs = [(a, b, "Event-Event") for a in variables for b in variables]
            pairs += [(a, b, "Event-Spectral") for a in variables for b in FEATURES]
            pairs += [(a, b, "Spectral-Spectral") for a in FEATURES for b in FEATURES]
            for a, b, kind in pairs:
                x, y = _finite_numeric(merged[a]), _finite_numeric(merged[b])
                valid = np.isfinite(x) & np.isfinite(y); n = int(valid.sum()); r = pval = np.nan
                if n >= 3 and np.unique(x[valid]).size > 1 and np.unique(y[valid]).size > 1:
                    from scipy.stats import pearsonr
                    r, pval = pearsonr(x[valid], y[valid])
                feature1 = a.removesuffix(pdt) if kind in {"Event-Event", "Event-Spectral"} else a
                feature2 = b.removesuffix(pdt) if kind == "Event-Event" else b
                records.append({"Embryo": embryo_id, "Group": group_name, "PDT": pdt,
                    "CorrelationType": kind, "Feature1": feature1, "Feature2": feature2,
                    "N": n, "PearsonR": r, "PValue": pval})
    if not records:
        raise ValueError("No correlation results were generated")
    pd.DataFrame(records).to_csv(root / "All_Embryos_All_PDT_Correlations.csv", index=False)


def run_correlation_summary(root: Path) -> None:
    frame = pd.read_csv(env_path("CALCIUM_CORRELATION_FILE", file=True))
    require_columns(frame, ["Group", "Feature1", "Feature2", "PearsonR"])
    features = EVENT_PREFIXES + FEATURES
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    for group, rows in frame.groupby("Group", sort=False):
        shape = (len(features), len(features)); med=np.full(shape,np.nan); lo=med.copy(); hi=med.copy(); sd=med.copy(); count=np.zeros(shape,int)
        for i, a in enumerate(features):
            for j in range(i, len(features)):
                b=features[j]
                values=pd.to_numeric(rows.loc[(rows.Feature1==a)&(rows.Feature2==b),"PearsonR"],errors="coerce").to_numpy(float)
                values=values[np.isfinite(values)]
                if values.size:
                    med[i,j]=np.median(values); lo[i,j]=np.min(values); hi[i,j]=np.max(values); sd[i,j]=np.std(values,ddof=1) if values.size>1 else np.nan; count[i,j]=values.size
                med[j,i],lo[j,i],hi[j,i],sd[j,i],count[j,i]=med[i,j],lo[i,j],hi[i,j],sd[i,j],count[i,j]
        for label, values in [("Median_R",med),("Min_R",lo),("Max_R",hi),("Std_R",sd),("N_Correlations",count)]:
            pd.DataFrame(values,index=features,columns=features).to_csv(root/f"{safe_name(group)}_{label}.csv",index_label="Feature")
        fig, ax=plt.subplots(figsize=(15,13)); im=ax.imshow(med,vmin=-1,vmax=1,cmap="viridis"); fig.colorbar(im,ax=ax)
        ax.set_xticks(range(len(features)),features,rotation=60,ha="right"); ax.set_yticks(range(len(features)),features)
        ax.set_title(f"{group}: median Pearson r across embryos and PDTs")
        for i in range(len(features)):
            for j in range(len(features)):
                if np.isfinite(med[i,j]): ax.text(j,i,f"{med[i,j]:.2f}\n[{lo[i,j]:.2f}, {hi[i,j]:.2f}]",ha="center",va="center",fontsize=5)
        fig.tight_layout(); fig.savefig(root/f"{safe_name(group)}_Median_Correlation_Matrix.png",dpi=200); plt.close(fig)


def _boundary_xy(frame: pd.DataFrame) -> np.ndarray:
    columns=list(frame.columns)
    xcol=next((c for c in ["X","Centroid_1","Boundary_1","V1","Var1"] if c in columns),None)
    ycol=next((c for c in ["Y","Centroid_2","Boundary_2","V2","Var2"] if c in columns),None)
    if (xcol is None or ycol is None) and len(columns)==2: xcol,ycol=columns
    if xcol is None or ycol is None: raise ValueError("Unable to identify boundary X/Y columns")
    xy=frame[[xcol,ycol]].apply(pd.to_numeric,errors="coerce").to_numpy(float)
    xy=xy[np.isfinite(xy).all(axis=1)]
    if len(xy)<3: raise ValueError("Boundary has fewer than three valid points")
    if not np.array_equal(xy[0],xy[-1]): xy=np.vstack([xy,xy[0]])
    return xy


def run_individual_maps(root: Path) -> None:
    boundary_root=env_path("CALCIUM_BOUNDARY_ROOT",file=False)
    master=pd.read_csv(env_path("CALCIUM_CELL_MASTER",file=True),dtype={"EmbryoID":str,"CellID":str})
    require_columns(master,["EmbryoID","CellID",*FEATURES])
    boundary_files=sorted(boundary_root.rglob("*-Boundary_SuperEmbryo.csv"))
    if not boundary_files: raise FileNotFoundError(f"No boundary CSVs under {boundary_root}")
    merged_list=[]; ids=[]; x_all=[]; y_all=[]; values={f:[] for f in FEATURES}
    for bpath in boundary_files:
        embryo=bpath.name.removesuffix("-Boundary_SuperEmbryo.csv")
        boundary=pd.read_csv(bpath); require_columns(boundary,["Calcium_ID","Centroid_1","Centroid_2"])
        mask=master.EmbryoID.astype(str).eq(embryo)|master.EmbryoID.astype(str).str.startswith(embryo+"_")
        spec=master.loc[mask].copy()
        if spec.empty: warnings.warn(f"No spectral rows for {embryo}"); merged_list.append(None); ids.append(embryo); continue
        spec["CellID"]=spec.CellID.astype(str); boundary["Calcium_ID"]=boundary.Calcium_ID.astype(str)
        merged=boundary.merge(spec,left_on="Calcium_ID",right_on="CellID",how="inner",validate="many_to_many")
        if merged.empty: warnings.warn(f"No matched cells for {embryo}"); merged_list.append(None); ids.append(embryo); continue
        merged_list.append(merged); ids.append(embryo)
        valid=np.isfinite(merged.Centroid_1)&np.isfinite(merged.Centroid_2); x_all.extend(merged.loc[valid,"Centroid_1"]); y_all.extend(merged.loc[valid,"Centroid_2"])
        for feature in FEATURES:
            v=pd.to_numeric(merged[feature],errors="coerce").to_numpy(float); values[feature].extend(v[np.isfinite(v)])
    if not x_all: raise ValueError("No matched rows have finite centroids")
    xmargin=.02*(max(x_all)-min(x_all)); ymargin=.02*(max(y_all)-min(y_all)); xlim=(min(x_all)-xmargin,max(x_all)+xmargin); ylim=(min(y_all)-ymargin,max(y_all)+ymargin)
    limits={f:(np.percentile(v,[1,99]) if len(v) else [0,1]) for f,v in values.items()}
    import matplotlib.pyplot as plt
    for embryo, merged in zip(ids,merged_list):
        if merged is None: continue
        folder=root/safe_name(embryo); folder.mkdir()
        merged.to_csv(folder/f"{embryo}_Spectral_Mapped.csv",index=False)
        for feature in FEATURES:
            vals=pd.to_numeric(merged[feature],errors="coerce").to_numpy(float)
            if not np.isfinite(vals).any(): continue
            lo,hi=limits[feature]
            if lo==hi: lo,hi=lo-max(abs(lo)*.05,1e-6),hi+max(abs(hi)*.05,1e-6)
            fig,ax=plt.subplots(figsize=(8,6)); sc=ax.scatter(merged.Centroid_1,merged.Centroid_2,c=vals,s=25,vmin=lo,vmax=hi)
            ax.set(xlim=xlim,ylim=ylim,aspect="equal",xlabel="Template X",ylabel="Template Y",title=f"{embryo} — {feature}"); fig.colorbar(sc,ax=ax,label=feature)
            fig.tight_layout(); fig.savefig(folder/f"{embryo}_Mapped_{feature}.png",dpi=300); plt.close(fig)


def run_combine_maps(root: Path) -> None:
    map_root=env_path("CALCIUM_MAPPED_ROOT",file=False)
    embryo_summary=pd.read_csv(env_path("CALCIUM_EMBRYO_SUMMARY",file=True),dtype={"EmbryoID":str})
    require_columns(embryo_summary,["EmbryoID","Group"])
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    for feature in FEATURES:
        with PdfPages(root/f"{feature}_AllEmbryos.pdf") as pdf:
            for group, metadata in embryo_summary.groupby("Group",sort=False):
                paths=[]
                available=list(map_root.rglob(f"*_Mapped_{feature}.png"))
                for embryo in metadata.EmbryoID.astype(str):
                    hits=[]
                    for path in available:
                        mapped_id=path.name.removesuffix(f"_Mapped_{feature}.png")
                        if embryo==mapped_id or embryo.startswith(mapped_id+"_"):
                            hits.append(path)
                    if len(hits)==1: paths.append((embryo,hits[0]))
                    elif len(hits)>1: warnings.warn(f"Multiple mapped plots match {embryo} ({feature}); skipping")
                for start in range(0,len(paths),16):
                    page=paths[start:start+16]; fig,axes=plt.subplots(4,4,figsize=(15,11)); fig.suptitle(f"{feature} — {group}")
                    for ax,(embryo,path) in zip(axes.flat,page): ax.imshow(plt.imread(path)); ax.set_title(embryo,fontsize=8); ax.axis("off")
                    for ax in axes.flat[len(page):]: ax.axis("off")
                    fig.tight_layout(); pdf.savefig(fig,dpi=200); plt.close(fig)


def run_group_spatial(root: Path) -> None:
    mapped_root=env_path("CALCIUM_MAPPED_ROOT",file=False)
    assignments=pd.read_csv(env_path("CALCIUM_GROUP_ASSIGNMENTS",file=True),dtype={"Embryo":str,"Group":str})
    require_columns(assignments,["Embryo","Group"])
    if assignments.Embryo.duplicated().any(): raise ValueError("Group assignments contain duplicate embryos")
    mapped_files=sorted(mapped_root.rglob("*_Spectral_Mapped.csv"))
    if not mapped_files: raise FileNotFoundError(f"No mapped CSVs under {mapped_root}")
    tables=[]
    for path in mapped_files:
        embryo=path.name.removesuffix("_Spectral_Mapped.csv")
        assigned=assignments.Embryo.astype(str)
        match=assignments.loc[assigned.str.lower().eq(embryo.lower()) | assigned.str.lower().str.startswith(embryo.lower()+"_") | pd.Series([embryo.lower().startswith(x.lower()+"_") for x in assigned],index=assigned.index)]
        if len(match)!=1: warnings.warn(f"No unique group assignment for {embryo}; skipping"); continue
        frame=pd.read_csv(path); require_columns(frame,["Centroid_1","Centroid_2",*FEATURES])
        frame["Embryo"]=embryo; frame["Group"]=match.Group.iloc[0]; tables.append(frame)
    if not tables: raise ValueError("No mapped rows have unique group assignments")
    cells=pd.concat(tables,ignore_index=True); cells["Centroid_1"]=pd.to_numeric(cells.Centroid_1,errors="coerce"); cells["Centroid_2"]=pd.to_numeric(cells.Centroid_2,errors="coerce")
    cells=cells[np.isfinite(cells.Centroid_1)&np.isfinite(cells.Centroid_2)].copy()
    cells["CNPU_X"]=np.floor((cells.Centroid_1-GRID_X0)/CNPU_SIZE).astype(int); cells["CNPU_Y"]=np.floor((cells.Centroid_2-GRID_Y0)/CNPU_SIZE).astype(int)
    cells["CNPU_ID"]=cells.CNPU_X.astype(str)+"_"+cells.CNPU_Y.astype(str); cells.to_csv(root/"All_MappedCells_with_CNPU.csv",index=False)
    gx=np.arange(0,math.ceil((X_RANGE[1]-GRID_X0)/CNPU_SIZE)); gy=np.arange(0,math.ceil((Y_RANGE[1]-GRID_Y0)/CNPU_SIZE))
    grid=pd.MultiIndex.from_product([gx,gy],names=["CNPU_X","CNPU_Y"]).to_frame(index=False)
    grid["CNPU_StartX"]=GRID_X0+grid.CNPU_X*CNPU_SIZE; grid["CNPU_EndX"]=grid.CNPU_StartX+CNPU_SIZE
    grid["CNPU_StartY"]=GRID_Y0+grid.CNPU_Y*CNPU_SIZE; grid["CNPU_EndY"]=grid.CNPU_StartY+CNPU_SIZE
    grid["CNPU_CenterX"]=grid.CNPU_StartX+CNPU_SIZE/2; grid["CNPU_CenterY"]=grid.CNPU_StartY+CNPU_SIZE/2
    observed=cells.groupby(["CNPU_X","CNPU_Y"]).size().rename("NCells")
    grid=grid.merge(observed,on=["CNPU_X","CNPU_Y"],how="left"); grid.NCells=grid.NCells.fillna(0).astype(int); grid["Observed"]=grid.NCells>0
    b13=_boundary_xy(pd.read_csv(env_path("CALCIUM_STAGE13_SOX2_BOUNDARY",file=True)))
    b14=_boundary_xy(pd.read_csv(env_path("CALCIUM_STAGE14_SOX2_BOUNDARY",file=True)))
    from matplotlib.path import Path as MplPath
    grid["Stage13Template"]=MplPath(b13).contains_points(grid[["CNPU_CenterX","CNPU_CenterY"]].to_numpy())
    grid["Stage14Template"]=MplPath(b14).contains_points(grid[["CNPU_CenterX","CNPU_CenterY"]].to_numpy())
    grid["AnySox2Template"]=grid.Stage13Template|grid.Stage14Template; grid.to_csv(root/"Complete_CNPU_Grid.csv",index=False)
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(15,5))
    for ax,col,title in zip(axes,["Observed","Stage13Template","Stage14Template"],["Observed CNPUs","Stage 13 template","Stage 14 template"]):
        ax.scatter(grid.CNPU_CenterX,grid.CNPU_CenterY,c=grid[col].astype(int),s=8); ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_title(title)
    fig.tight_layout(); fig.savefig(root/"CompleteCNPUGrid_QC.png",dpi=250); plt.close(fig)
    summaries=[]
    for feature in FEATURES:
        valid=pd.to_numeric(cells[feature],errors="coerce"); temp=cells.assign(_value=valid); temp=temp[np.isfinite(temp._value)]
        per=(temp.groupby(["Group","Embryo","CNPU_X","CNPU_Y"],as_index=False)
             .agg(EmbryoMedian=("_value","median"),NCells=("_value","size")))
        group_counts=cells.groupby("Group").Embryo.nunique().to_dict()
        observed_stats=(per.groupby(["Group","CNPU_X","CNPU_Y"],as_index=False)
                       .agg(GroupMedian=("EmbryoMedian","median"),GroupSD=("EmbryoMedian","std"),NEmbryos=("EmbryoMedian","count"),NCells=("NCells","sum")))
        observed_stats["CoveragePct"]=[100*n/group_counts[g] for g,n in zip(observed_stats.Group,observed_stats.NEmbryos)]
        base=grid[["CNPU_X","CNPU_Y","CNPU_CenterX","CNPU_CenterY"]]
        for group in cells.Group.dropna().unique():
            sub=base.merge(observed_stats[observed_stats.Group==group],on=["CNPU_X","CNPU_Y"],how="left")
            sub["Group"]=group; sub["Feature"]=feature; sub["NEmbryos"]=sub.NEmbryos.fillna(0).astype(int); sub["NCells"]=sub.NCells.fillna(0).astype(int); sub["CoveragePct"]=sub.CoveragePct.fillna(0)
            sub["GroupMedian"]=sub.GroupMedian.astype(float); sub["GroupSD"]=sub.GroupSD.astype(float); summaries.append(sub)
    summary=pd.concat(summaries,ignore_index=True)
    summary["SmoothedMedian"]=np.nan; summary["SmoothedZ"]=np.nan
    for (feature,group), idx in summary.groupby(["Feature","Group"]).groups.items():
        block=summary.loc[idx]; coords=list(zip(block.CNPU_X,block.CNPU_Y)); lookup=dict(zip(coords,block.GroupMedian)); smooth=[]
        for x,y in coords:
            vals=[lookup.get((x+dx,y+dy),np.nan) for dx in range(-3,4) for dy in range(-3,4)]
            vals=np.asarray(vals,float); vals=vals[np.isfinite(vals)]; smooth.append(float(np.mean(vals)) if vals.size else np.nan)
        smooth=np.asarray(smooth); summary.loc[idx,"SmoothedMedian"]=smooth; ok=np.isfinite(smooth)
        if ok.any() and np.std(smooth[ok],ddof=1 if ok.sum()>1 else 0)>0:
            summary.loc[np.asarray(idx)[ok],"SmoothedZ"]=(smooth[ok]-np.mean(smooth[ok]))/np.std(smooth[ok],ddof=1 if ok.sum()>1 else 0)
    summary.to_csv(root/"Group_CNPU_Summaries.csv",index=False)
    for feature in FEATURES:
        rows=summary[summary.Feature==feature]
        for metric in ["GroupMedian","GroupSD","CoveragePct"]:
            lo,hi=np.nanmin(rows[metric]),np.nanmax(rows[metric]);
            if not np.isfinite(lo+hi) or lo==hi: lo,hi=0,1
            fig,axes=plt.subplots(2,2,figsize=(10,8))
            for ax,group in zip(axes.flat,rows.Group.drop_duplicates()):
                data=rows[rows.Group==group]; sc=ax.scatter(data.CNPU_CenterX,data.CNPU_CenterY,c=data[metric],s=20,vmin=lo,vmax=hi)
                boundary=b13 if "13" in str(group) else b14; ax.plot(boundary[:,0],boundary[:,1],":k",lw=1); ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_title(str(group)); fig.colorbar(sc,ax=ax)
            fig.suptitle(f"{feature} — {metric}"); fig.tight_layout(); fig.savefig(root/f"{feature}_{metric}_GroupSpatialMaps.png",dpi=250); plt.close(fig)


def run_boundary_qc(root: Path) -> None:
    boundary_root=env_path("CALCIUM_BOUNDARY_ROOT",file=False); records=[]
    for path in sorted(boundary_root.rglob("*-Boundary_SuperEmbryo.csv")):
        frame=pd.read_csv(path); require_columns(frame,["Centroid_1","Centroid_2"])
        boundary_cols=sorted((c for c in frame.columns if c.startswith("Boundary_")),key=lambda c:int(re.search(r"(\d+)$",c).group(1)))
        with_boundary=0
        for _,row in frame[boundary_cols].iterrows():
            vals=pd.to_numeric(row,errors="coerce").to_numpy(float); vals=vals[np.isfinite(vals)]
            if len(vals)>=6 and len(vals)%2==0: with_boundary+=1
        centroids=np.isfinite(pd.to_numeric(frame.Centroid_1,errors="coerce"))&np.isfinite(pd.to_numeric(frame.Centroid_2,errors="coerce"))
        records.append({"BoundaryFile":path.name,"Rows":len(frame),"CellsWithBoundary":with_boundary,"CellsWithCentroid":int(centroids.sum())})
    pd.DataFrame(records).to_csv(root/"Boundary_QC.csv",index=False)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise
