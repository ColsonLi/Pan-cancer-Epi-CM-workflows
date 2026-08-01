"""Broad clustering — user-specified (pcs=25, nn=30, res=0.4). SKIP grid.

Per SKILL.md: "If the user explicitly specifies n_neighbors, n_pcs, and
optionally resolution, run only the user-specified combination or resolution-
search path and record it as user-specified." User explicitly chose to skip the
default 64-candidate grid and run a single combination.

Steps:
  1. Read h5ad/04-integration-harmony/adata_harmony.h5ad (clean, fresh process).
  2. rsc.get.anndata_to_GPU(adata)
  3. rsc.pp.neighbors(n_neighbors=30, n_pcs=25, use_rep='X_pca_inte')
  4. rsc.tl.leiden(resolution=0.4, key_added='leiden_coarse')
  5. rsc.tl.umap
  6. Save h5ad/05-clustering-parameter-search/selected/adata_inte.h5ad

Outputs:
  - h5ad/01-celltype_integration_clustering/05-clustering-parameter-search/selected/adata_inte.h5ad
  - tables/.../selected_clustering.csv
  - tables/.../cluster_counts.csv
  - figures/.../umap_leiden_res0p4.pdf
  - figures/.../umap_sample.pdf (optional diagnostic)
"""
from __future__ import annotations

import json
import sys
import time
import gc
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

SEED = 42
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
IN_H5AD = ROOT / "epi-cm-core-workflow/h5ad/01-celltype_integration_clustering/04-integration-harmony/adata_harmony.h5ad"
OUT_H5AD = ROOT / "epi-cm-core-workflow/h5ad/01-celltype_integration_clustering/05-clustering-parameter-search/selected/adata_inte.h5ad"
TAB = ROOT / "epi-cm-core-workflow/tables/01-celltype_integration_clustering/05-clustering-parameter-search"
TAB.mkdir(parents=True, exist_ok=True)
FIG = ROOT / "epi-cm-core-workflow/figures/01-celltype_integration_clustering/05-clustering-parameter-search"
FIG.mkdir(parents=True, exist_ok=True)
OUT_H5AD.parent.mkdir(parents=True, exist_ok=True)

# User-specified parameters
N_PCS = 25
N_NEIGHBORS = 30
LEIDEN_RES = 0.4
LEIDEN_KEY = "leiden_coarse"


def main() -> None:
    t0 = time.time()
    print(f"[clust] reading {IN_H5AD.name}…", flush=True)
    adata = sc.read_h5ad(IN_H5AD)
    print(f"[clust] loaded {adata.shape} in {time.time()-t0:.1f}s", flush=True)

    # Sanity check that X_pca_inte is present
    assert "X_pca_inte" in adata.obsm, "X_pca_inte missing from adata_harmony.h5ad"
    print(f"[clust] X_pca_inte shape: {adata.obsm['X_pca_inte'].shape}", flush=True)

    # GPU
    import rapids_singlecell as rsc
    print("[clust] moving to GPU…", flush=True)
    rsc.get.anndata_to_GPU(adata)

    # Neighbors
    print(f"[clust] neighbors (n_pcs={N_PCS}, n_neighbors={N_NEIGHBORS})…", flush=True)
    t1 = time.time()
    rsc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS, use_rep="X_pca_inte", random_state=SEED)
    print(f"[clust] neighbors done in {time.time()-t1:.1f}s", flush=True)

    # Leiden
    print(f"[clust] leiden (resolution={LEIDEN_RES})…", flush=True)
    t1 = time.time()
    rsc.tl.leiden(adata, resolution=LEIDEN_RES, key_added=LEIDEN_KEY, random_state=SEED)
    n_clusters = int(adata.obs[LEIDEN_KEY].nunique())
    print(f"[clust] leiden done in {time.time()-t1:.1f}s; n_clusters={n_clusters}", flush=True)

    # UMAP
    print("[clust] UMAP…", flush=True)
    t1 = time.time()
    rsc.tl.umap(adata, random_state=SEED)
    print(f"[clust] UMAP done in {time.time()-t1:.1f}s", flush=True)

    # Move back to CPU before write
    rsc.get.anndata_to_CPU(adata)
    gc.collect()

    # Save the selected object
    print(f"[clust] writing {OUT_H5AD.name}…", flush=True)
    adata.write_h5ad(OUT_H5AD, compression="gzip")

    # Cluster counts
    cc = adata.obs[LEIDEN_KEY].value_counts().sort_index()
    pd.DataFrame({"cluster": cc.index.astype(str), "n_cells": cc.values}).to_csv(
        TAB / "cluster_counts.csv", index=False
    )
    print(f"[clust] cluster counts:\n{cc}", flush=True)

    # selected_clustering.csv
    pd.DataFrame([{
        "selected_h5ad": str(OUT_H5AD),
        "source_h5ad": str(IN_H5AD),
        "n_pcs": N_PCS,
        "n_neighbors": N_NEIGHBORS,
        "resolution": LEIDEN_RES,
        "leiden_key": LEIDEN_KEY,
        "n_clusters": n_clusters,
        "use_rep": "X_pca_inte",
        "umap_key": "X_umap",
        "random_seed": SEED,
        "backend": "rapids_singlecell (GPU)",
        "manual_selection_note": "user-specified; SKIP grid (64 candidates not run)",
    }]).to_csv(TAB / "selected_clustering.csv", index=False)
    print(f"[clust] selected_clustering.csv written", flush=True)

    # Plot UMAP
    sc.settings.figdir = str(FIG)
    sc.set_figure_params(figsize=(5, 5), dpi=150)
    sc.pl.umap(adata, color=LEIDEN_KEY, save="_leiden_res0p4.pdf", show=False)
    sc.pl.umap(adata, color="sample", save="_sample_diagnostic.pdf", show=False)
    sc.pl.umap(adata, color="status", save="_status_diagnostic.pdf", show=False)
    sc.pl.umap(adata, color="series", save="_series_diagnostic.pdf", show=False)
    print(f"[clust] UMAP figures saved under {FIG}", flush=True)

    print(f"[clust] done. total={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()