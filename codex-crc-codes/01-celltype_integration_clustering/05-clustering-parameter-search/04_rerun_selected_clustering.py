#!/usr/bin/env python3
"""Rerun and save the manually selected broad clustering."""

from __future__ import annotations

import gc
import platform
import random
import traceback
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc

matplotlib.use("Agg")

SEED = 42
random.seed(SEED)
np.random.seed(SEED)

WORKFLOW_ROOT = Path("/mnt/disk18t/lr_xcy/riku/crc_val/epi-cm-core-workflow")
STEP = "01-celltype_integration_clustering/05-clustering-parameter-search"
INPUT_H5AD = WORKFLOW_ROOT / "h5ad/01-celltype_integration_clustering/04-integration-harmony/adata_harmony.h5ad"
OUTPUT_H5AD = WORKFLOW_ROOT / "h5ad" / STEP / "selected/adata_inte.h5ad"
TABLE_DIR = WORKFLOW_ROOT / "tables" / STEP
SELECTED_TABLE_DIR = TABLE_DIR / "selected"
FIGURE_DIR = WORKFLOW_ROOT / "figures" / STEP / "selected"
CODE_FILE = Path(__file__)

N_PCS = 30
N_NEIGHBORS = 30
RESOLUTION = 0.3
BASIS = "X_pca_inte"
NEIGHBORS_KEY = "neighbors"
UMAP_KEY = "X_umap"
LEIDEN_KEY = "leiden_res0p3"
SELECTED_LABEL = "pcs-30_nn-30_res-0p3"


def assert_no_overwrite() -> None:
    outputs = [
        OUTPUT_H5AD,
        TABLE_DIR / "selected_clustering.csv",
        SELECTED_TABLE_DIR / "selected_cluster_counts.csv",
        SELECTED_TABLE_DIR / "selected_clustering_parameters.csv",
        SELECTED_TABLE_DIR / "package_versions.txt",
        SELECTED_TABLE_DIR / "readme.txt",
        FIGURE_DIR / "umap_selected_leiden_res0p3.pdf",
    ]
    existing = [str(path) for path in outputs if path.exists()]
    if existing:
        raise FileExistsError(
            "Selected clustering output already exists; refusing to overwrite:\n"
            + "\n".join(existing)
        )


def main() -> None:
    assert_no_overwrite()
    OUTPUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
    SELECTED_TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    search = None

    import cupy as cp
    import rapids_singlecell as rsc

    sc.settings.autoshow = False
    sc.settings.figdir = str(FIGURE_DIR)
    sc.settings.set_figure_params(figsize=(4, 4), dpi=180)

    status = {
        "step": "selected_clustering_start",
        "input_h5ad": str(INPUT_H5AD),
        "selected_h5ad": str(OUTPUT_H5AD),
        "selected_label": SELECTED_LABEL,
        "n_pcs": N_PCS,
        "n_neighbors": N_NEIGHBORS,
        "resolution": RESOLUTION,
        "basis": BASIS,
        "neighbors_key": NEIGHBORS_KEY,
        "umap_key": UMAP_KEY,
        "raw_leiden_cluster_column": LEIDEN_KEY,
        "backend": "rapids_singlecell",
        "random_seed": SEED,
        "code_file": str(CODE_FILE),
    }
    pd.DataFrame([status]).to_csv(SELECTED_TABLE_DIR / "selected_clustering_status.csv", index=False)

    try:
        adata = ad.read_h5ad(INPUT_H5AD)
        if BASIS not in adata.obsm:
            raise KeyError(f"Missing obsm[{BASIS!r}] in {INPUT_H5AD}")

        search = ad.AnnData(
            X=np.zeros((adata.n_obs, 1), dtype=np.float32),
            obs=adata.obs.copy(),
        )
        search.obsm[BASIS] = np.asarray(adata.obsm[BASIS]).astype(np.float32)

        rsc.get.anndata_to_GPU(search)
        rsc.pp.neighbors(
            search,
            n_neighbors=N_NEIGHBORS,
            n_pcs=N_PCS,
            use_rep=BASIS,
            random_state=SEED,
            key_added=NEIGHBORS_KEY,
        )
        rsc.tl.umap(
            search,
            random_state=SEED,
            neighbors_key=NEIGHBORS_KEY,
            key_added=UMAP_KEY,
        )
        rsc.tl.leiden(
            search,
            resolution=RESOLUTION,
            key_added=LEIDEN_KEY,
            random_state=SEED,
            neighbors_key=NEIGHBORS_KEY,
        )
        rsc.get.anndata_to_CPU(search)

        adata.obsm[UMAP_KEY] = np.asarray(search.obsm[UMAP_KEY])
        adata.obs[LEIDEN_KEY] = search.obs[LEIDEN_KEY].copy()
        adata.uns[NEIGHBORS_KEY] = search.uns[NEIGHBORS_KEY].copy()
        for key in search.obsp.keys():
            adata.obsp[key] = search.obsp[key].copy()

        counts = (
            adata.obs[LEIDEN_KEY]
            .astype(str)
            .value_counts()
            .sort_index()
            .rename_axis("cluster")
            .reset_index(name="n_cells")
        )
        n_clusters = int(counts.shape[0])
        counts.to_csv(SELECTED_TABLE_DIR / "selected_cluster_counts.csv", index=False)

        params = {
            **status,
            "step": "selected_clustering_complete",
            "n_obs": int(adata.n_obs),
            "n_vars": int(adata.n_vars),
            "n_clusters": n_clusters,
            "candidate_grid_label": SELECTED_LABEL,
            "candidate_h5ad_saved": False,
            "manual_selection_note": "User selected pcs=30, n_neighbors=30, resolution=0.3 after UMAP grid review.",
            "selected_object_note": "Graph, UMAP, and Leiden were rerun from the clean Harmony h5ad before saving this selected object.",
        }
        pd.DataFrame([params]).to_csv(SELECTED_TABLE_DIR / "selected_clustering_parameters.csv", index=False)

        selected = pd.DataFrame(
            [
                {
                    "selected_h5ad": str(OUTPUT_H5AD),
                    "input_h5ad": str(INPUT_H5AD),
                    "selected_label": SELECTED_LABEL,
                    "final_raw_leiden_cluster_column": LEIDEN_KEY,
                    "n_pcs": N_PCS,
                    "n_neighbors": N_NEIGHBORS,
                    "resolution": RESOLUTION,
                    "basis": BASIS,
                    "neighbors_key": NEIGHBORS_KEY,
                    "umap_key": UMAP_KEY,
                    "n_clusters": n_clusters,
                    "manual_selection_note": "User selected pcs=30, n_neighbors=30, resolution=0.3.",
                    "code_file": str(CODE_FILE),
                    "random_seed": SEED,
                }
            ]
        )
        selected.to_csv(TABLE_DIR / "selected_clustering.csv", index=False)

        sc.pl.umap(adata, color=LEIDEN_KEY, save="_selected_leiden_res0p3.pdf", show=False)
        adata.write_h5ad(OUTPUT_H5AD)

        with (SELECTED_TABLE_DIR / "package_versions.txt").open("w") as fh:
            fh.write(f"python: {platform.python_version()}\n")
            fh.write(f"anndata: {ad.__version__}\n")
            fh.write(f"scanpy: {sc.__version__}\n")
            fh.write(f"rapids_singlecell: {rsc.__version__}\n")
            fh.write(f"cupy: {cp.__version__}\n")
            fh.write(f"code_file: {CODE_FILE}\n")

        with (SELECTED_TABLE_DIR / "readme.txt").open("w") as fh:
            fh.write("Selected broad clustering completed.\n")
            fh.write(f"Input h5ad: {INPUT_H5AD}\n")
            fh.write(f"Selected h5ad: {OUTPUT_H5AD}\n")
            fh.write(f"Selected parameters: n_pcs={N_PCS}, n_neighbors={N_NEIGHBORS}, resolution={RESOLUTION}\n")
            fh.write(f"Raw Leiden cluster column: obs['{LEIDEN_KEY}']\n")
            fh.write("Note: obs['leiden_coarse'] is not assigned here; it belongs to 06-broad-annotation after marker review.\n")
            fh.write("Figure: figures/01-celltype_integration_clustering/05-clustering-parameter-search/selected/umap_selected_leiden_res0p3.pdf\n")
            fh.write("Tables: selected_clustering.csv, selected/selected_clustering_parameters.csv, selected/selected_cluster_counts.csv\n")

        status.update(
            {
                "step": "selected_clustering_complete",
                "n_obs": int(adata.n_obs),
                "n_vars": int(adata.n_vars),
                "n_clusters": n_clusters,
            }
        )
        pd.DataFrame([status]).to_csv(SELECTED_TABLE_DIR / "selected_clustering_status.csv", index=False)
        print(pd.DataFrame([status]).to_string(index=False))

    except Exception as exc:
        status.update(
            {
                "step": "selected_clustering_failed",
                "error_summary": "".join(traceback.format_exception_only(type(exc), exc)).strip(),
            }
        )
        pd.DataFrame([status]).to_csv(SELECTED_TABLE_DIR / "selected_clustering_status.csv", index=False)
        raise
    finally:
        try:
            if search is not None:
                rsc.get.anndata_to_CPU(search)
        except Exception:
            pass
        try:
            cp.get_default_memory_pool().free_all_blocks()
        except Exception:
            pass
        gc.collect()


if __name__ == "__main__":
    main()
