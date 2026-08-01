"""Per-celltype DEG (rank_genes_groups grouped by cell_type, use_raw=True).

User request: after annotation, run rank_genes_groups once more with the
cell_type column to get per-celltype marker genes (e.g., T markers, B markers).
Save one CSV per cell type.

Outputs:
  tables/06-broad-annotation/per_celltype_degs/celltype_<NAME>_degs.csv
  tables/06-broad-annotation/celltype_top_degs.csv (combined)
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

SEED = 42
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
IN_H5AD = ROOT / "epi-cm-core-workflow/h5ad/01-celltype_integration_clustering/06-broad-annotation/adata_anno.h5ad"
TAB = ROOT / "epi-cm-core-workflow/tables/01-celltype_integration_clustering/06-broad-annotation"
OUT_DIR = TAB / "per_celltype_degs"
OUT_DIR.mkdir(parents=True, exist_ok=True)
N_TOP = 50


def main() -> None:
    t0 = time.time()
    print(f"[ctdeg] reading {IN_H5AD.name}…", flush=True)
    adata = sc.read_h5ad(IN_H5AD)
    print(f"[ctdeg] loaded {adata.shape} in {time.time()-t0:.1f}s", flush=True)
    assert adata.raw is not None, "adata.raw missing"

    print(f"[ctdeg] rank_genes_groups (groupby=cell_type, use_raw=True, t-test)…", flush=True)
    t1 = time.time()
    sc.tl.rank_genes_groups(
        adata,
        groupby="cell_type",
        use_raw=True,
        method="t-test",
        n_genes=N_TOP,
        key_added="rank_genes_celltype",
    )
    print(f"[ctdeg] rank_genes_groups done in {time.time()-t1:.1f}s", flush=True)

    rgg = adata.uns["rank_genes_celltype"]
    groups = list(rgg["names"].dtype.names)
    print(f"[ctdeg] cell types: {groups}", flush=True)

    rows = []
    for g in groups:
        for rank in range(N_TOP):
            rows.append({
                "cell_type": g,
                "rank": rank + 1,
                "gene": rgg["names"][g][rank],
                "logfc": float(rgg["logfoldchanges"][g][rank]),
                "pval": float(rgg["pvals"][g][rank]),
                "padj": float(rgg["pvals_adj"][g][rank]),
                "score": float(rgg["scores"][g][rank]),
            })
    df = pd.DataFrame(rows)
    df.to_csv(TAB / "celltype_top_degs.csv", index=False)
    print(f"[ctdeg] wrote celltype_top_degs.csv ({len(df)} rows)", flush=True)

    # One CSV per cell type
    for g in groups:
        sub = df[df["cell_type"] == g].sort_values("rank").reset_index(drop=True)
        safe = g.replace("/", "_")
        out = OUT_DIR / f"celltype_{safe}_degs.csv"
        sub.to_csv(out, index=False)
        print(f"  {g}: {len(sub)} DEGs -> {out.name}", flush=True)

    # Also print top 10 per cell type for the chat
    print(f"\n[ctdeg] top 10 per cell_type:")
    for g in groups:
        sub = df[df["cell_type"] == g].head(10)
        genes = ", ".join(sub["gene"].tolist())
        print(f"  {g}: {genes}", flush=True)

    print(f"[ctdeg] done. total={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()