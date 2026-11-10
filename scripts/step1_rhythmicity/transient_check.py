#!/usr/bin/env python3
"""Exploratory check (2026-11-04): are red-noise-"rhythmic" cells oscillating, or do they carry a few slow transients?
Samples up to 150 rhythmic + 150 non-rhythmic cells in 3 median-like embryos per devitellinated group plus 180125SP_E2_st_13,
rebuilds dF/F segments exactly as the v3 pipeline, counts transients with scipy find_peaks (prominence >= max(6 x robust
white-noise SD, 0.1 dF/F)), and reports events/h, half-height width, inter-event interval and CV_ISI.
Usage: python transient_check.py <output_dir>  (writes transients.csv). Quick rule, not a validated event detector.
"""
import numpy as np, pandas as pd
from scipy.signal import find_peaks, peak_widths
sys.path.insert(0,"/sciclone/home/cabare01/CALCIUM/scripts/analysis_common"); import common
RT=common.load_module(common.SCRIPTS/"step1_rhythmicity"/"rhythmicity_test.py","rt")
r=common.load_rhythmicity()
axis=RT.spectral.frequency_axis(); coi=RT.spectral.coi_factor()
pct=r.groupby(["Group","EmbryoID"]).Rhythmic.agg(["mean","sum"])
picks=[]
for g in ["S13_devitellinated","S14_devitellinated"]:
    p=pct.loc[g]; p=p[p["sum"]>=40]; med=pct.loc[g]["mean"].median()
    picks += [(g,e) for e in (p["mean"]-med).abs().sort_values().index[:3]]
picks.append(("S13","180125SP_E2_st_13"))
rng=np.random.default_rng(0); rows=[]
for g,emb in picks:
    gdir=common.SPECTRAL_ROOT/common.GROUP_FOLDER[g]; man=json.loads((gdir/"run_manifest.json").read_text())
    raw={common.safe_name(Path(x).stem.removesuffix(" - Not_self_Merged")):x for x in man["input_files"]}[emb]
    t=pd.read_csv(raw,low_memory=False)
    try: info=RT.preprocess.embryo_frame_times(t)
    except Exception: info=None
    if info is not None: t["_time_s"]=info["time"]
    by={str(k):v for k,v in t.groupby("ID",sort=False)}
    d=r[(r.Group==g)&(r.EmbryoID==emb)]
    for flag in (True,False):
        ids=d[d.Rhythmic==flag].CellID.to_numpy(); ids=rng.choice(ids,min(150,len(ids)),replace=False)
        for cid in ids:
            rowsc=by[cid]; ci=None if info is None else {**info,"time":rowsc["_time_s"].to_numpy()}
            segs=RT.preprocess.preprocess_data(rowsc,float(axis[0]),coi,ci)["segments"]
            times=[];widths=[];dur=0
            for s in segs:
                x=np.asarray(s["dFF"],float); dt=s["dt"]; dur+=len(x)*dt
                noise=1.4826*np.median(np.abs(np.diff(x)))/np.sqrt(2)  # robust white-noise SD
                pk,pr=find_peaks(x,prominence=max(6*noise,0.1))
                if len(pk):
                    w=peak_widths(x,pk,rel_height=0.5)[0]*dt; widths+=list(w); times+=list(s["time"][pk])
            times=np.sort(times); isi=np.diff(times)
            rows.append(dict(Group=g,EmbryoID=emb,Rhythmic=flag,NEvents=len(times),EventsPerHr=len(times)/(dur/3600),
                 MedWidth_s=np.median(widths) if widths else np.nan,
                 MedISI_s=np.median(isi) if len(isi) else np.nan,
                 CV_ISI=isi.std()/isi.mean() if len(isi)>=2 else np.nan))
    print("done",g,emb,flush=True)
df=pd.DataFrame(rows); df.to_csv(sys.argv[1]+"/transients.csv",index=False)
pd.set_option("display.width",250)
print(df.groupby(["Group","EmbryoID","Rhythmic"]).agg(n=("NEvents","size"),medEvents=("NEvents","median"),
   pct0=("NEvents",lambda x:100*(x==0).mean()),pct1_3=("NEvents",lambda x:100*x.between(1,3).mean()),pct_ge6=("NEvents",lambda x:100*(x>=6).mean()),
   medWidth=("MedWidth_s","median"),medISI=("MedISI_s","median"),medCV=("CV_ISI","median")).round(2))
a=df[df.Rhythmic]; print("\nALL rhythmic sampled:", len(a), "median events", a.NEvents.median(), "| % with 1-3 events", round(100*a.NEvents.between(1,3).mean(),1), "| % >=6", round(100*(a.NEvents>=6).mean(),1), "| median width s", round(a.MedWidth_s.median(),1), "| median ISI s", round(a.MedISI_s.median(),1), "| median CV_ISI (>=3 events)", round(a.CV_ISI.median(),2), "| % CV<0.5 among >=6 events", round(100*(a[a.NEvents>=6].CV_ISI<0.5).mean(),1))
