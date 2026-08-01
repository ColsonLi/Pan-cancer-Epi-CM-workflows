#!/usr/bin/env python3
"""Rerun the user-selected BRCA clustering from the clean Harmony handoff."""

from __future__ import annotations

import gc
import importlib.metadata
import json
import os
import random
import sys
import time
from pathlib import Path

# Required on this host for cudf/cugraph host transfers with the current
# NVIDIA driver; set before importing RAPIDS/Numba-dependent packages.
os.environ.setdefault("NUMBA_CUDA_USE_NVIDIA_BINDING", "1")

import cupy as cp
import cupyx.scipy.sparse as cpsparse
import numpy as np
import pandas as pd
import rapids_singlecell as rsc
import scanpy as sc


SEED = 42
N_PCS = 20
N_NEIGHBORS = 30
RESOLUTION = 0.8
CLUSTER_KEY = "leiden_res0p8"
USE_REP = "X_pca_inte"
UMAP_MIN_DIST = 0.5
UMAP_SPREAD = 1.0

random.seed(SEED)
np.random.seed(SEED)
cp.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/brca_cm-epi_codex")
WORKFLOW = ROOT / "epi-cm-core-workflow"
BLOCK = "01-celltype_integration_clustering"
TASK = "05-clustering-parameter-search"
INPUT_H5AD = WORKFLOW / "h5ad" / BLOCK / "04-integration-harmony" / "adata_harmony.h5ad"
OUT_H5AD = WORKFLOW / "h5ad" / BLOCK / TASK / "selected" / "adata_inte.h5ad"
TABLE_ROOT = WORKFLOW / "tables" / BLOCK / TASK
SELECTED_TABLE_DIR = TABLE_ROOT / "selected"
SELECTED_RECORD = TABLE_ROOT / "selected_clustering.csv"
FIGURE_DIR = WORKFLOW / "figures" / BLOCK / TASK / "selected"
MANIFEST = TABLE_ROOT / "clustering_grid_manifest.csv"
GRID_COUNTS = (
    TABLE_ROOT
    / "pcs-20_nn-30_res-0p1-1p0"
    / "cluster_counts.csv"
)
CODE_PATH = WORKFLOW / "codes" / BLOCK / TASK / "07_run_selected_clustering.py"


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def cpuify(value):
    if isinstance(value, cp.ndarray):
        return cp.asnumpy(value)
    if cpsparse.isspmatrix(value):
        return value.get()
    if isinstance(value, dict):
        return {key: cpuify(item) for key, item in value.items()}
    if isinstance(value, list):
        return [cpuify(item) for item in value]
    if isinstance(value, tuple):
        return tuple(cpuify(item) for item in value)
    return value


def convert_aux_to_cpu(adata) -> None:
    for slot in [adata.obsm, adata.varm, adata.obsp, adata.varp, adata.layers, adata.uns]:
        for key in list(slot.keys()):
            slot[key] = cpuify(slot[key])


def release_gpu() -> None:
    cp.get_default_memory_pool().free_all_blocks()
    cp.get_default_pinned_memory_pool().free_all_blocks()
    gc.collect()


def main() -> None:
    started = time.time()
    OUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
    SELECTED_TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(MANIFEST)
    manifest_row = manifest[
        manifest["n_pcs"].eq(N_PCS)
        & manifest["n_neighbors"].eq(N_NEIGHBORS)
        & np.isclose(manifest["resolution"].astype(float), RESOLUTION)
    ]
    if len(manifest_row) != 1:
        raise ValueError("Selected candidate is not unique in the grid manifest.")
    manifest_row = manifest_row.iloc[0]
    if str(manifest_row["status"]) != "completed" or not bool(manifest_row["completed"]):
        raise ValueError("User-selected candidate is not a completed grid candidate.")
    grid_counts = pd.read_csv(GRID_COUNTS)
    grid_row = grid_counts[np.isclose(grid_counts["resolution"], RESOLUTION)]
    if len(grid_row) != 1 or str(grid_row.iloc[0]["status"]) != "completed":
        raise ValueError("Selected resolution is not completed in its graph table.")
    grid_cluster_count = int(grid_row.iloc[0]["cluster_count"])

    source = sc.read_h5ad(INPUT_H5AD)
    if source.n_obs != 99866:
        raise ValueError(f"Unexpected source cell count: {source.n_obs}")
    if USE_REP not in source.obsm or source.obsm[USE_REP].shape[1] < N_PCS:
        raise ValueError(f"Missing or insufficient {USE_REP} dimensions.")
    if not source.obs_names.is_unique:
        raise ValueError("Harmony source cell IDs are not unique.")
    if {"neighbors", "umap", "leiden"} & set(source.uns):
        raise ValueError("Harmony source is not a clean graph-free handoff.")
    if "X_umap" in source.obsm:
        raise ValueError("Harmony source unexpectedly contains X_umap.")
    source_names = source.obs_names.copy()

    # The skill requires an explicit selected-run copy before mutation.
    adata = source.copy()
    del source
    print(
        f"[Selected] clean source copied: cells={adata.n_obs}, genes={adata.n_vars}, "
        f"raw_genes={adata.raw.n_vars if adata.raw is not None else 0}",
        flush=True,
    )

    step_times: dict[str, float] = {}
    t0 = time.time()
    rsc.get.anndata_to_GPU(adata)
    step_times["anndata_to_GPU"] = time.time() - t0
    t0 = time.time()
    rsc.pp.neighbors(
        adata,
        n_neighbors=N_NEIGHBORS,
        n_pcs=N_PCS,
        use_rep=USE_REP,
        random_state=SEED,
        algorithm="brute",
        metric="euclidean",
    )
    step_times["neighbors"] = time.time() - t0
    t0 = time.time()
    rsc.tl.umap(
        adata,
        min_dist=UMAP_MIN_DIST,
        spread=UMAP_SPREAD,
        random_state=SEED,
    )
    step_times["umap"] = time.time() - t0
    t0 = time.time()
    rsc.tl.leiden(
        adata,
        resolution=RESOLUTION,
        random_state=SEED,
        key_added=CLUSTER_KEY,
        n_iterations=100,
    )
    step_times["leiden"] = time.time() - t0
    selected_cluster_count = int(adata.obs[CLUSTER_KEY].nunique())
    print(
        f"[Selected] graph complete: PCs={N_PCS}, NN={N_NEIGHBORS}, "
        f"res={RESOLUTION:.1f}, clusters={selected_cluster_count}",
        flush=True,
    )

    t0 = time.time()
    rsc.get.anndata_to_CPU(adata)
    convert_aux_to_cpu(adata)
    release_gpu()
    step_times["anndata_to_CPU"] = time.time() - t0

    if not adata.obs_names.equals(source_names):
        raise ValueError("Selected rerun changed cell IDs or their order.")
    if adata.n_obs != 99866:
        raise ValueError("Selected rerun changed the cell count.")
    if CLUSTER_KEY not in adata.obs or adata.obs[CLUSTER_KEY].isna().any():
        raise ValueError("Selected raw Leiden labels are missing or incomplete.")
    if "leiden_coarse" in adata.obs or "cell_type" in adata.obs:
        raise ValueError("Pure selected clustering must not create annotation columns.")
    if "X_umap" not in adata.obsm or adata.obsm["X_umap"].shape != (adata.n_obs, 2):
        raise ValueError("Selected UMAP is missing or has the wrong shape.")
    if not np.isfinite(np.asarray(adata.obsm["X_umap"])).all():
        raise ValueError("Selected UMAP contains non-finite values.")
    if adata.raw is None or adata.raw.n_vars != 27716:
        raise ValueError("Normalized/log raw expression was not preserved.")

    adata.uns["selected_clustering"] = {
        "manual_selection": True,
        "manual_selection_note": (
            "User explicitly selected PCs=20, n_neighbors=30, resolution=0.8 "
            "after reviewing the completed grid."
        ),
        "source_h5ad": str(INPUT_H5AD),
        "integration_source": "Harmony",
        "use_rep": USE_REP,
        "n_pcs": N_PCS,
        "n_neighbors": N_NEIGHBORS,
        "resolution": RESOLUTION,
        "algorithm": "leiden",
        "raw_cluster_column": CLUSTER_KEY,
        "seed": SEED,
        "umap_min_dist": UMAP_MIN_DIST,
        "umap_spread": UMAP_SPREAD,
        "backend": "rapids_singlecell_gpu",
        "grid_cluster_count": grid_cluster_count,
        "selected_rerun_cluster_count": selected_cluster_count,
    }

    sc.settings.figdir = FIGURE_DIR
    sc.settings.file_format_figs = "pdf"
    sc.set_figure_params(figsize=(3, 3), dpi=150, fontsize=8)
    sc.pl.umap(
        adata,
        color=CLUSTER_KEY,
        show=False,
        save="_selected_leiden_res0p8.pdf",
    )
    sc.pl.umap(
        adata,
        color=[CLUSTER_KEY, "series", "status"],
        ncols=3,
        wspace=0.5,
        show=False,
        save="_selected_cluster_series_status.pdf",
    )
    selected_figure = FIGURE_DIR / "umap_selected_leiden_res0p8.pdf"
    diagnostic_figure = FIGURE_DIR / "umap_selected_cluster_series_status.pdf"
    for figure in [selected_figure, diagnostic_figure]:
        if not figure.is_file() or figure.stat().st_size == 0:
            raise FileNotFoundError(f"Expected selected figure was not created: {figure}")

    cluster_sizes = (
        adata.obs[CLUSTER_KEY]
        .value_counts(sort=False)
        .rename_axis(CLUSTER_KEY)
        .rename("n_cells")
        .reset_index()
    )
    cluster_sizes.to_csv(SELECTED_TABLE_DIR / "selected_cluster_sizes.csv", index=False)
    sample_fractions = pd.crosstab(
        adata.obs[CLUSTER_KEY].astype(str),
        adata.obs["sample"].astype(str),
        normalize="index",
    )
    sample_summary = pd.DataFrame(
        {
            CLUSTER_KEY: sample_fractions.index,
            "max_sample_fraction": sample_fractions.max(axis=1).to_numpy(),
            "dominant_sample": sample_fractions.idxmax(axis=1).to_numpy(),
        }
    )
    sample_summary.to_csv(
        SELECTED_TABLE_DIR / "selected_cluster_sample_mixing.csv", index=False
    )
    pd.DataFrame(
        [{"step": key, "elapsed_seconds": value} for key, value in step_times.items()]
    ).to_csv(SELECTED_TABLE_DIR / "selected_step_times.csv", index=False)

    selected_relative = "h5ad/01-celltype_integration_clustering/05-clustering-parameter-search/selected/adata_inte.h5ad"
    selected_record = pd.DataFrame(
        [
            {
                "selected_h5ad": selected_relative,
                "selected_h5ad_absolute": str(OUT_H5AD),
                "source_h5ad": str(INPUT_H5AD),
                "integration_source": "Harmony",
                "use_rep": USE_REP,
                "algorithm": "leiden",
                "final_raw_cluster_column": CLUSTER_KEY,
                "n_pcs": N_PCS,
                "n_neighbors": N_NEIGHBORS,
                "resolution": RESOLUTION,
                "grid_candidate_status": str(manifest_row["status"]),
                "grid_candidate_completed": bool(manifest_row["completed"]),
                "grid_cluster_count": grid_cluster_count,
                "selected_rerun_cluster_count": selected_cluster_count,
                "selected_figure": str(selected_figure),
                "diagnostic_figure": str(diagnostic_figure),
                "manual_selection": True,
                "manual_selection_note": (
                    "User explicitly selected PCs=20, n_neighbors=30, resolution=0.8."
                ),
                "selected_run_started_from_clean_harmony": True,
                "selected_run_created_by_copy_before_mutation": True,
                "qc_or_harmony_rerun_inside_selected_clustering": False,
                "all_cells_preserved": True,
                "seed": SEED,
                "backend": "rapids_singlecell_gpu",
            }
        ]
    )
    selected_record.to_csv(SELECTED_RECORD, index=False)

    versions = {
        "python": sys.version.split()[0],
        "scanpy": package_version("scanpy"),
        "anndata": package_version("anndata"),
        "rapids-singlecell": package_version("rapids-singlecell"),
        "cuml-cu12": package_version("cuml-cu12"),
        "cugraph-cu12": package_version("cugraph-cu12"),
        "cupy-cuda12x": package_version("cupy-cuda12x"),
        "cuda-python": package_version("cuda-python"),
        "numba-cuda": package_version("numba-cuda"),
        "NUMBA_CUDA_USE_NVIDIA_BINDING": os.environ.get(
            "NUMBA_CUDA_USE_NVIDIA_BINDING", ""
        ),
        "code": str(CODE_PATH),
        "seed": str(SEED),
    }
    (SELECTED_TABLE_DIR / "package_versions.txt").write_text(
        "\n".join(f"{key}={value}" for key, value in versions.items()) + "\n",
        encoding="utf-8",
    )
    (SELECTED_TABLE_DIR / "readme.txt").write_text(
        f"""BRCA selected broad clustering

User-selected parameters: PCs={N_PCS}, neighbors={N_NEIGHBORS}, Leiden resolution={RESOLUTION:.1f}
Source: {INPUT_H5AD}
Output: {OUT_H5AD}
Raw technical cluster column: {CLUSTER_KEY}

The graph, UMAP, and Leiden result were rerun from a clean read of the saved
Harmony handoff after an explicit selected-run copy was created. All 99,866
cells were retained. QC, normalization, PCA, and Harmony were not rerun. The
selected object intentionally does not contain leiden_coarse or cell_type;
those biological columns belong to the subsequent broad-annotation task.
Candidate-grid figures, tables, manifests, and code were preserved.
""",
        encoding="utf-8",
    )

    print(f"[Selected] writing {OUT_H5AD}", flush=True)
    adata.write_h5ad(OUT_H5AD, compression="gzip")
    report = {
        "n_cells": int(adata.n_obs),
        "n_hvg": int(adata.n_vars),
        "n_raw_genes": int(adata.raw.n_vars),
        "raw_cluster_column": CLUSTER_KEY,
        "grid_cluster_count": grid_cluster_count,
        "selected_rerun_cluster_count": selected_cluster_count,
        "selected_h5ad": str(OUT_H5AD),
        "selected_figure": str(selected_figure),
        "elapsed_seconds": time.time() - started,
    }
    (SELECTED_TABLE_DIR / "selected_completion.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
