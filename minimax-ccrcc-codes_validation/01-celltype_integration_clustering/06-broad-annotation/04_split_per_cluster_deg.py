"""Split cluster_top_degs.csv into one CSV per cluster.

User: "每个cluster都保存一个deg的csv" — one DEG CSV per cluster.
"""
from __future__ import annotations

from pathlib import Path
import pandas as pd

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
TAB = ROOT / "epi-cm-core-workflow/tables/01-celltype_integration_clustering/06-broad-annotation"
DEG_DIR = TAB / "per_cluster_degs"
DEG_DIR.mkdir(parents=True, exist_ok=True)

inp = TAB / "cluster_top_degs.csv"
df = pd.read_csv(inp)
clusters = sorted(df["cluster"].astype(str).unique(), key=lambda x: int(x))
print(f"[split] {len(clusters)} clusters", flush=True)
for cl in clusters:
    sub = df[df["cluster"].astype(str) == cl].sort_values("rank").reset_index(drop=True)
    out = DEG_DIR / f"cluster_{cl}_degs.csv"
    sub.to_csv(out, index=False)
    print(f"  C{cl}: {len(sub)} DEGs -> {out.name}", flush=True)
print(f"[split] done", flush=True)
