"""Redraw umap_leiden_coarse.pdf with current leiden_coarse = cell_type.

The previous umap_leiden_coarse.pdf was drawn at 11:45 when leiden_coarse
contained cluster numbers 0-13. After 07_fix_leiden_coarse_deg.py, the
leiden_coarse column was reassigned to the broad cell type label.
"""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import scanpy as sc

SEED = 42
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
IN_H5AD = ROOT / "epi-cm-core-workflow/h5ad/01-celltype_integration_clustering/06-broad-annotation/adata_anno.h5ad"
FIG = ROOT / "epi-cm-core-workflow/figures/01-celltype_integration_clustering/06-broad-annotation"


def main() -> None:
    t0 = time.time()
    print(f"[redraw] reading {IN_H5AD.name}…", flush=True)
    adata = sc.read_h5ad(IN_H5AD)
    print(f"[redraw] loaded {adata.shape} in {time.time()-t0:.1f}s", flush=True)
    print(f"[redraw] leiden_coarse value_counts:\n{adata.obs['leiden_coarse'].value_counts()}", flush=True)

    sc.settings.figdir = str(FIG)
    sc.set_figure_params(figsize=(5, 5), dpi=150)
    sc.pl.umap(adata, color="leiden_coarse", save="_leiden_coarse.pdf", show=False)
    print(f"[redraw] umap_leiden_coarse.pdf redrawn", flush=True)
    print(f"[redraw] done. total={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()