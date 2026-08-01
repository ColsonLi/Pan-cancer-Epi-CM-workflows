"""Load merged atlas from adata_merge_raw_barcode_* files and write adata_merge.h5ad.

Per user instruction: SKIP sample-merge step. The triplet
(adata_merge_raw_barcode_counts_csr.h5, obs.tsv.gz, var.tsv.gz) is already a
cross-dataset merged atlas. This script only loads it and writes the canonical
Module 01 h5ad/02-merge-metadata/adata_merge.h5ad so downstream tasks have a
stable path.
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

SEED = 42
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
COUNTS = ROOT / "adata_merge_raw_barcode_counts_csr.h5"
OBS = ROOT / "adata_merge_raw_barcode_obs.tsv.gz"
VAR = ROOT / "adata_merge_raw_barcode_var.tsv.gz"
OUT_H5AD = ROOT / "epi-cm-core-workflow/h5ad/01-celltype_integration_clustering/02-merge-metadata/adata_merge.h5ad"
OUT_H5AD.parent.mkdir(parents=True, exist_ok=True)


def main() -> None:
    t0 = time.time()
    print(f"[load_merge] reading {COUNTS.name} (CSR sparse, backed)…", flush=True)
    # The file is the CSR triplet already on disk; load via h5py-backed sparse
    from scipy.sparse import csr_matrix
    import h5py

    with h5py.File(COUNTS, "r") as f:
        g = f["X"]
        keys = list(g.keys())
        print(f"[load_merge] h5 /X keys: {keys}", flush=True)
        data = g["data"][:]
        indices = g["indices"][:]
        indptr = g["indptr"][:]
    # Use the obs/var-derived shape (CSR indptr has n_rows+1 entries)
    n_rows = indptr.shape[0] - 1
    n_cols = pd.read_csv(VAR, sep="\t", index_col=0).shape[0]
    X = csr_matrix((data, indices, indptr), shape=(n_rows, n_cols))
    print(f"[load_merge] CSR loaded: {X.shape}, nnz={X.nnz}, dtype={X.dtype}", flush=True)

    obs = pd.read_csv(OBS, sep="\t", index_col=0)
    var = pd.read_csv(VAR, sep="\t", index_col=0)
    print(f"[load_merge] obs={obs.shape}, var={var.shape}", flush=True)

    # Align obs/var order to matrix axes
    assert X.shape[0] == obs.shape[0], "obs row count mismatch with X"
    assert X.shape[1] == var.shape[0], "var row count mismatch with X"
    obs.index = obs.index.astype(str)
    var.index = var.index.astype(str)

    adata = ad.AnnData(
        X=X,
        obs=obs.copy(),
        var=var.copy(),
    )
    adata.var_names_make_unique()
    print(f"[load_merge] AnnData: {adata}", flush=True)
    print(f"[load_merge] status counts:\n{adata.obs['status'].value_counts()}", flush=True)
    print(f"[load_merge] series counts:\n{adata.obs['series'].value_counts()}", flush=True)

    # Write h5ad (compressed)
    t1 = time.time()
    print(f"[load_merge] writing {OUT_H5AD} …", flush=True)
    adata.write_h5ad(OUT_H5AD, compression="gzip")
    t2 = time.time()
    print(f"[load_merge] done. load={t1-t0:.1f}s, write={t2-t1:.1f}s, total={t2-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()