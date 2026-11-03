#!/usr/bin/env python3
"""Step 5a: index PDT event files (*ConvergeSpikeCy5*.csv) against the reference table.

Searches --roots recursively (default: data10/boundary_info plus
data10/pdt_events, the suggested home for event files retrieved from the PDT
analysis), matches each file to a reference embryo by the normalised key
(date + initials + embryo number; folder names are ignored, the reference
table decides the group), and lists the reference embryos still missing.

Outputs in ANALYSIS_ROOT/step5_event_integration/:
  event_file_index.csv    Group, EmbryoKey, EventFile (one row per reference embryo)
  missing_event_files.csv reference embryos with no event file -> request list
  duplicate_event_files.csv  keys with several candidate files (resolve by hand;
                          the index then uses the path given in --prefer, if any)
"""
from __future__ import annotations

import argparse
import hashlib
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis_common"))
import common  # noqa: E402

LABEL = {col: v[1] for col, v in common.embryo_ref.GROUP_COLUMNS.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--roots", nargs="+", type=Path,
                    default=[common.HOME / "data10/boundary_info", common.HOME / "data10/pdt_events"])
    ap.add_argument("--prefer", type=Path, nargs="*", default=[],
                    help="Files to use when a key has several candidates.")
    args = ap.parse_args()
    out = common.out_dir("step5_event_integration")

    found = defaultdict(list)
    for root in args.roots:
        if not root.exists():
            print(f"(root not found, skipped: {root})")
            continue
        for p in root.rglob("*ConvergeSpikeCy5*.csv"):
            found[common.embryo_ref.embryo_key(p.name)].append(p.resolve())
    prefer = {p.resolve() for p in args.prefer}
    ref = common.embryo_ref.load_reference()

    rows, missing, dups = [], [], []
    for key, column in sorted(ref.items()):
        # Byte-identical copies (e.g. the six files duplicated in "wrong stage
        # folder need to cross reference") count as one candidate.
        by_hash = {}
        for c in sorted(set(found.get(key, []))):
            by_hash.setdefault(hashlib.md5(c.read_bytes()).hexdigest(), c)
        cands = sorted(by_hash.values())
        chosen = None
        if len(cands) == 1:
            chosen = cands[0]
        elif len(cands) > 1:
            pick = [c for c in cands if c in prefer]
            chosen = pick[0] if len(pick) == 1 else None
            dups.extend({"EmbryoKey": key, "Candidate": str(c), "Chosen": c == chosen} for c in cands)
        row = {"Group": LABEL[column], "EmbryoKey": key, "EventFile": str(chosen) if chosen else ""}
        rows.append(row)
        if chosen is None:
            missing.append({**row, "NCandidates": len(cands)})
    unref = sorted(k for k in found if k not in ref)
    pd.DataFrame(rows).to_csv(out / "event_file_index.csv", index=False)
    pd.DataFrame(missing).to_csv(out / "missing_event_files.csv", index=False)
    pd.DataFrame(dups).to_csv(out / "duplicate_event_files.csv", index=False)
    have = pd.DataFrame(rows).assign(Has=lambda d: d.EventFile != "")
    print(have.groupby("Group")["Has"].agg(["sum", "size"]).rename(columns={"sum": "WithEvents",
                                                                            "size": "Reference"}))
    if unref:
        print(f"Event files not in the reference table (ignored): {unref}")


if __name__ == "__main__":
    main()
