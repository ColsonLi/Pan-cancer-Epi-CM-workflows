"""Fix: re-assign leiden_coarse = cell_type and re-export leiden_coarse DEG.

The previous 06_full_degs.py did not replace the existing 'leiden_coarse' column
(branched on the wrong condition), so the leiden_coarse DEG directory still
contains cluster-number CSVs (0..13). This script:
  1. Reassigns adata.obs['leiden_coarse'] = adata.obs['cell_type']
  2. Removes the wrong CSVs from degs_leiden_coarse_pcs25_nn30_res0p4/
  3. Re-runs rank_genes_groups(groupby='leiden_coarse', use_raw=True) with full output
  4. Writes one CSV per cell type label
  5. Updates adata_anno.h5ad
"""
from __future__ import annotations

import time
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

SEED = 42
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
IN_H5AD = ROOT / "epi-cm-core-workflow/h5ad/01-celltype_integration_clustering/06-broad-annotation/adata_anno.h5ad"
TAB = ROOT / "epi-cm-core-workflow/tables/01-celltype_integration_clustering/06-broad-annotation"

PCS = 25
NN = 30
RES = "0p4"
COARSE_COL = "leiden_coarse"
DEG_DIR_COARSE = TAB / f"degs_{COARSE_COL}_pcs{PCS}_nn{NN}_res{RES}"


def main() -> None:
    t0 = time.time()
    print(f"[fix] reading {IN_H5AD.name}…", flush=True)
    adata = sc.read_h5ad(IN_H5AD)
    print(f"[fix] loaded {adata.shape} in {time.time()-t0:.1f}s", flush=True)

    # Reassign leiden_coarse = cell_type
    if "cell_type" not in adata.obs:
        raise RuntimeError("cell_type column missing; cannot re-derive leiden_coarse")
    adata.obs[COARSE_COL] = adata.obs["cell_type"].astype(str).astype("category")
    print(f"[fix] reassigned '{COARSE_COL}' = cell_type ({adata.obs[COARSE_COL].nunique()} unique)")
    print(adata.obs[COARSE_COL].value_counts().to_string(), flush=True)

    # Clear wrong CSVs
    if DEG_DIR_COARSE.exists():
        print(f"[fix] clearing wrong CSVs in {DEG_DIR_COARSE.name}…", flush=True)
        shutil.rmtree(DEG_DIR_COARSE)
    DEG_DIR_COARSE.mkdir(parents=True, exist_ok=True)

    # Re-run full DEG on leiden_coarse
    print(f"[fix] rank_genes_groups (groupby={COARSE_COL}, use_raw=True, t-test, full)…", flush=True)
    t1 = time.time()
    sc.tl.rank_genes_groups(
        adata,
        groupby=COARSE_COL,
        use_raw=True,
        method="t-test",
        n_genes=None,
        key_added=f"rank_genes_{COARSE_COL}",
    )
    print(f"[fix] rank_genes_groups done in {time.time()-t1:.1f}s", flush=True)

    rgg = adata.uns[f"rank_genes_{COARSE_COL}"]
    groups = list(rgg["names"].dtype.names)
    print(f"[fix] {len(groups)} groups: {groups}", flush=True)

    for g in groups:
        n = len(rgg["names"][g])
        df = pd.DataFrame({
            "gene": rgg["names"][g],
            "score": rgg["scores"][g],
            "logfoldchanges": rgg["logfoldchanges"][g],
            "pvals": rgg["pvals"][g],
            "pvals_adj": rgg["pvals_adj"][g],
        })
        safe = str(g).replace("/", "_").replace(" ", "_")
        fname = f"{safe}_degs_{COARSE_COL}_pcs{PCS}_nn{NN}_res{RES}.csv"
        out = DEG_DIR_COARSE / fname
        df.to_csv(out, index=False)
        print(f"  {g}: {n} genes -> {out.name}", flush=True)

    # Save updated h5ad
    print(f"\n[fix] writing updated {IN_H5AD.name}…", flush=True)
    adata.write_h5ad(IN_H5AD, compression="gzip")
    print(f"[fix] done. total={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()