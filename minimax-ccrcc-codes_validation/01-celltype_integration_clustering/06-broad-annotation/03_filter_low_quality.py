"""Filter out LowQuality and Doublets from adata_anno.h5ad.

User instruction: remove LowQuality (C12) and Doublets (C13) cells, keep the rest.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import scanpy as sc

SEED = 42
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
IN_H5AD = ROOT / "epi-cm-core-workflow/h5ad/01-celltype_integration_clustering/06-broad-annotation/adata_anno.h5ad"
OUT_H5AD = ROOT / "epi-cm-core-workflow/h5ad/01-celltype_integration_clustering/06-broad-annotation/adata_anno_filtered.h5ad"
TAB = ROOT / "epi-cm-core-workflow/tables/01-celltype_integration_clustering/06-broad-annotation"

DROP = ["LowQuality", "Doublets"]


def main() -> None:
    t0 = time.time()
    print(f"[filt] reading {IN_H5AD.name}…", flush=True)
    adata = sc.read_h5ad(IN_H5AD)
    print(f"[filt] loaded {adata.shape} in {time.time()-t0:.1f}s", flush=True)

    n_before = adata.n_obs
    keep_mask = ~adata.obs["cell_type"].astype(str).isin(DROP)
    adata = adata[keep_mask].copy()
    n_after = adata.n_obs
    print(f"[filt] dropped {DROP}: {n_before - n_after} cells; kept {n_after}", flush=True)

    # cell_type counts after filter
    print(f"[filt] cell_type counts after filter:")
    print(adata.obs["cell_type"].value_counts().to_string(), flush=True)

    # Per-celltype per-sample summary
    import pandas as pd
    grp = adata.obs.groupby(["cell_type", "sample"], observed=True).size().unstack(fill_value=0)
    grp.to_csv(TAB / "celltype_sample_counts_filtered.csv")
    grp2 = adata.obs.groupby(["cell_type", "status"], observed=True).size().unstack(fill_value=0)
    grp2.to_csv(TAB / "celltype_status_counts_filtered.csv")

    print(f"[filt] writing {OUT_H5AD.name}…", flush=True)
    adata.write_h5ad(OUT_H5AD, compression="gzip")
    print(f"[filt] done. total={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()