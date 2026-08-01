#!/usr/bin/env python3
"""Rerun one user-selected BRCA lineage subtype clustering from clean Harmony."""

from __future__ import annotations

import gc
import importlib.metadata
import json
import os
import random
import subprocess
import sys
import time
from pathlib import Path

os.environ.setdefault("NUMBA_CUDA_USE_NVIDIA_BINDING", "1")

import anndata as ad
import cupy as cp
import numpy as np
import pandas as pd
import rapids_singlecell as rsc
import scanpy as sc
from scipy import sparse


SEED = 42
random.seed(SEED)
np.random.seed(SEED)
cp.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/brca_cm-epi_codex")
WORKFLOW = ROOT / "epi-cm-core-workflow"
BLOCK = "02-cell_subtype_integration_clustering"
INTEGRATED_BASIS = "X_pca_inte"
UMAP_MIN_DIST = 0.5
UMAP_SPREAD = 1.0
SELECTION_NOTE = "User confirmed all primary shortlist candidates on 2026-07-14."

LINEAGE_CONFIG = {
    "epithelial": {
        "label": "Epithelial Cells",
        "abbrev": "epi",
        "n_pcs": 25,
        "n_neighbors": 35,
        "resolution": 0.5,
        "expected_clusters": 10,
    },
    "t_cells": {
        "label": "T Cells",
        "abbrev": "t",
        "n_pcs": 20,
        "n_neighbors": 30,
        "resolution": 0.7,
        "expected_clusters": 10,
    },
    "myeloid": {
        "label": "Myeloid Cells",
        "abbrev": "mye",
        "n_pcs": 25,
        "n_neighbors": 25,
        "resolution": 0.5,
        "expected_clusters": 9,
    },
    "b_cells": {
        "label": "B Cells",
        "abbrev": "b",
        "n_pcs": 25,
        "n_neighbors": 35,
        "resolution": 0.6,
        "expected_clusters": 7,
    },
    "plasma": {
        "label": "Plasma Cells",
        "abbrev": "plasma",
        "n_pcs": 25,
        "n_neighbors": 35,
        "resolution": 0.4,
        "expected_clusters": 7,
    },
    "endothelial": {
        "label": "Endothelial Cells",
        "abbrev": "endo",
        "n_pcs": 20,
        "n_neighbors": 35,
        "resolution": 0.5,
        "expected_clusters": 10,
    },
    "stromal": {
        "label": "Stromal Cells",
        "abbrev": "stromal",
        "n_pcs": 20,
        "n_neighbors": 35,
        "resolution": 0.7,
        "expected_clusters": 8,
    },
    "perivascular": {
        "label": "Perivascular Cells",
        "abbrev": "pvl",
        "n_pcs": 20,
        "n_neighbors": 25,
        "resolution": 0.6,
        "expected_clusters": 7,
    },
}

LINEAGE_SLUG = "epithelial"
CODE_PATH = Path(__file__).resolve()


def configure_lineage(lineage_slug: str, code_path: Path | None = None) -> None:
    global LINEAGE_SLUG, CODE_PATH
    if lineage_slug not in LINEAGE_CONFIG:
        raise ValueError(f"Unsupported lineage slug: {lineage_slug}")
    LINEAGE_SLUG = lineage_slug
    if code_path is not None:
        CODE_PATH = code_path


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def res_token(resolution: float) -> str:
    return f"{resolution:.1f}".replace(".", "p")


def pdf_is_readable(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    result = subprocess.run(["pdfinfo", str(path)], capture_output=True, text=True)
    return result.returncode == 0 and "Pages:" in result.stdout


def gpu_memory_snapshot() -> dict[str, int | str]:
    command = [
        "nvidia-smi",
        "--query-gpu=index,memory.total,memory.used,memory.free",
        "--format=csv,noheader,nounits",
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0 or not result.stdout.strip():
        return {
            "gpu_index": "unavailable",
            "gpu_total_mb_preflight": -1,
            "gpu_used_mb_preflight": -1,
            "gpu_free_mb_preflight": -1,
        }
    fields = [field.strip() for field in result.stdout.splitlines()[0].split(",")]
    return {
        "gpu_index": fields[0],
        "gpu_total_mb_preflight": int(fields[1]),
        "gpu_used_mb_preflight": int(fields[2]),
        "gpu_free_mb_preflight": int(fields[3]),
    }


def release_gpu() -> None:
    cp.get_default_memory_pool().free_all_blocks()
    cp.get_default_pinned_memory_pool().free_all_blocks()
    gc.collect()


def selected_paths(config: dict[str, object]) -> dict[str, Path]:
    lineage = LINEAGE_SLUG
    table_dir = WORKFLOW / "tables" / BLOCK / "04-subtype-selected-clustering" / lineage
    figure_dir = WORKFLOW / "figures" / BLOCK / "04-subtype-selected-clustering" / lineage
    h5ad_dir = WORKFLOW / "h5ad" / BLOCK / "04-subtype-selected-clustering" / lineage
    return {
        "input_h5ad": (
            WORKFLOW
            / "h5ad"
            / BLOCK
            / "02-subtype-harmony"
            / lineage
            / f"adata_{lineage}_harmony.h5ad"
        ),
        "selected_ids": (
            WORKFLOW
            / "tables"
            / BLOCK
            / "01-lineage-selection"
            / lineage
            / "selected_cell_ids.csv"
        ),
        "grid_manifest": (
            WORKFLOW
            / "tables"
            / BLOCK
            / "03-subtype-clustering-grid"
            / lineage
            / "clustering_grid_manifest.csv"
        ),
        "table_dir": table_dir,
        "figure_dir": figure_dir,
        "h5ad_dir": h5ad_dir,
        "output_h5ad": h5ad_dir / f"adata_{config['abbrev']}_selected_clustered.h5ad",
        "figure": figure_dir / "umap_selected_raw_cluster.pdf",
        "parameters": table_dir / "selected_clustering.csv",
        "cluster_sizes": table_dir / "selected_cluster_sizes.csv",
        "sample_counts": table_dir / "selected_cluster_by_sample_counts.csv",
        "completion": table_dir / "selected_clustering_completion.json",
        "readme": table_dir / "readme.txt",
        "versions": table_dir / "package_versions.txt",
    }


def verify_grid_selection(paths: dict[str, Path], config: dict[str, object]) -> dict[str, object]:
    manifest = pd.read_csv(paths["grid_manifest"])
    mask = (
        manifest["n_pcs"].eq(int(config["n_pcs"]))
        & manifest["n_neighbors"].eq(int(config["n_neighbors"]))
        & np.isclose(manifest["resolution"].astype(float), float(config["resolution"]))
    )
    if int(mask.sum()) != 1:
        raise ValueError(f"Selected candidate is not unique in {paths['grid_manifest']}")
    row = manifest.loc[mask].iloc[0].to_dict()
    if str(row["status"]) != "completed" or str(row["completed"]).lower() != "true":
        raise ValueError(f"Selected candidate is not completed: {row}")
    if str(row["algorithm"]) != "leiden":
        raise ValueError(f"Selected candidate is not Leiden: {row['algorithm']}")
    if int(float(row["cluster_count"])) != int(config["expected_clusters"]):
        raise ValueError(
            "Selected grid cluster count changed: "
            f"manifest={row['cluster_count']}, expected={config['expected_clusters']}"
        )
    return row


def existing_run_is_valid(paths: dict[str, Path], config: dict[str, object]) -> bool:
    required = [
        paths["output_h5ad"],
        paths["figure"],
        paths["parameters"],
        paths["cluster_sizes"],
        paths["sample_counts"],
        paths["completion"],
        paths["readme"],
        paths["versions"],
    ]
    existing = [path.exists() for path in required]
    if not any(existing):
        return False
    if not all(existing):
        missing = [str(path) for path, present in zip(required, existing) if not present]
        raise FileExistsError(
            "Partial selected-clustering outputs already exist; refusing overwrite. "
            f"Missing: {missing}"
        )
    params = pd.read_csv(paths["parameters"]).iloc[0]
    expected = (
        int(params["n_pcs"]) == int(config["n_pcs"])
        and int(params["n_neighbors"]) == int(config["n_neighbors"])
        and np.isclose(float(params["resolution"]), float(config["resolution"]))
        and str(params["manual_selection_confirmed"]).lower() == "true"
    )
    raw_key = f"leiden_res{res_token(float(config['resolution']))}"
    saved = ad.read_h5ad(paths["output_h5ad"], backed="r")
    try:
        expected &= raw_key in saved.obs.columns
        expected &= int(saved.obs[raw_key].nunique()) == int(config["expected_clusters"])
        expected &= "X_umap" in saved.obsm
    finally:
        saved.file.close()
    expected &= pdf_is_readable(paths["figure"])
    if not expected:
        raise FileExistsError("Existing selected-clustering outputs do not match this selection.")
    return True


def sample_mixing_summary(obs: pd.DataFrame, raw_key: str) -> tuple[float, float]:
    fractions = pd.crosstab(
        obs[raw_key].astype(str), obs["sample"].astype(str), normalize="index"
    )
    dominance = fractions.max(axis=1)
    return float(dominance.mean()), float(dominance.max())


def write_selected_umap(source: ad.AnnData, raw_key: str, figure_dir: Path) -> Path:
    figure_dir.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = figure_dir
    sc.settings.file_format_figs = "pdf"
    sc.set_figure_params(figsize=(3, 3), dpi=150, fontsize=8)
    plot_adata = ad.AnnData(
        X=sparse.csr_matrix((source.n_obs, 0), dtype=np.float32),
        obs=source.obs[[raw_key]].copy(),
    )
    plot_adata.obsm["X_umap"] = np.asarray(source.obsm["X_umap"], dtype=np.float32)
    sc.pl.umap(
        plot_adata,
        color=raw_key,
        show=False,
        save="_selected_raw_cluster.pdf",
    )
    figure = figure_dir / "umap_selected_raw_cluster.pdf"
    if not pdf_is_readable(figure):
        raise FileNotFoundError(f"Scanpy did not create a readable selected UMAP: {figure}")
    color_key = f"{raw_key}_colors"
    if color_key in plot_adata.uns:
        source.uns[color_key] = list(plot_adata.uns[color_key])
    return figure


def main() -> None:
    started = time.time()
    config = LINEAGE_CONFIG[LINEAGE_SLUG]
    paths = selected_paths(config)
    grid_row = verify_grid_selection(paths, config)
    if existing_run_is_valid(paths, config):
        print(
            json.dumps(
                {
                    "lineage": LINEAGE_SLUG,
                    "status": "valid_existing_selected_outputs_reused",
                    "output_h5ad": str(paths["output_h5ad"]),
                    "selected_figure": str(paths["figure"]),
                },
                indent=2,
            )
        )
        return

    source = sc.read_h5ad(paths["input_h5ad"])
    selected_ids = pd.read_csv(paths["selected_ids"], usecols=["cell_id"], dtype=str)
    if source.n_obs != len(selected_ids) or not source.obs_names.is_unique:
        raise ValueError(
            f"Unexpected source cell IDs: observed={source.n_obs}, expected={len(selected_ids)}"
        )
    if INTEGRATED_BASIS not in source.obsm:
        raise ValueError(f"Missing {INTEGRATED_BASIS} in {paths['input_h5ad']}")
    if source.obsm[INTEGRATED_BASIS].shape[1] < int(config["n_pcs"]):
        raise ValueError(f"Harmony basis has too few PCs: {source.obsm[INTEGRATED_BASIS].shape}")
    if not source.obs["leiden_coarse"].astype(str).eq(str(config["label"])).all():
        raise ValueError(f"Off-lineage cells found in {paths['input_h5ad']}")
    required_obs = {"sample", "series", "status"}
    if not required_obs.issubset(source.obs.columns):
        raise ValueError(f"Source lacks obs columns: {sorted(required_obs - set(source.obs))}")

    gpu = gpu_memory_snapshot()
    graph = ad.AnnData(
        X=sparse.csr_matrix((source.n_obs, 0), dtype=np.float32),
        obs=source.obs[["sample", "series", "status"]].copy(),
    )
    graph.obsm[INTEGRATED_BASIS] = np.asarray(
        source.obsm[INTEGRATED_BASIS], dtype=np.float32
    )
    rsc.get.anndata_to_GPU(graph)
    print(
        f"[{LINEAGE_SLUG} selected] cells={source.n_obs} "
        f"pcs={config['n_pcs']} nn={config['n_neighbors']} "
        f"res={config['resolution']}",
        flush=True,
    )

    neighbors_started = time.time()
    rsc.pp.neighbors(
        graph,
        n_neighbors=int(config["n_neighbors"]),
        n_pcs=int(config["n_pcs"]),
        use_rep=INTEGRATED_BASIS,
        random_state=SEED,
        algorithm="brute",
        metric="euclidean",
    )
    neighbors_seconds = time.time() - neighbors_started
    umap_started = time.time()
    rsc.tl.umap(
        graph,
        min_dist=UMAP_MIN_DIST,
        spread=UMAP_SPREAD,
        random_state=SEED,
    )
    umap_seconds = time.time() - umap_started
    raw_key = f"leiden_res{res_token(float(config['resolution']))}"
    clustering_started = time.time()
    rsc.tl.leiden(
        graph,
        resolution=float(config["resolution"]),
        random_state=SEED,
        key_added=raw_key,
        n_iterations=100,
    )
    clustering_seconds = time.time() - clustering_started
    observed_clusters = int(graph.obs[raw_key].nunique())
    if observed_clusters != int(config["expected_clusters"]):
        raise RuntimeError(
            "Selected rerun did not reproduce the reviewed grid cluster count: "
            f"observed={observed_clusters}, expected={config['expected_clusters']}"
        )
    rsc.get.anndata_to_CPU(graph)

    required_obsp = {"connectivities", "distances"}
    if not required_obsp.issubset(graph.obsp.keys()):
        raise ValueError(f"Selected graph lacks obsp keys: {sorted(required_obsp - set(graph.obsp))}")
    if "neighbors" not in graph.uns or "X_umap" not in graph.obsm:
        raise ValueError("Selected graph lacks neighbors metadata or X_umap")

    categories = sorted(graph.obs[raw_key].astype(str).unique(), key=lambda value: int(value))
    source.obs[raw_key] = pd.Categorical(
        graph.obs[raw_key].astype(str), categories=categories, ordered=True
    )
    source.obsm["X_umap"] = np.asarray(graph.obsm["X_umap"], dtype=np.float32)
    source.obsp["connectivities"] = graph.obsp["connectivities"].copy()
    source.obsp["distances"] = graph.obsp["distances"].copy()
    source.uns["neighbors"] = graph.uns["neighbors"]
    if "umap" in graph.uns:
        source.uns["umap"] = graph.uns["umap"]
    if raw_key in graph.uns:
        source.uns[raw_key] = graph.uns[raw_key]
    source.uns["selected_subtype_clustering"] = {
        "lineage": LINEAGE_SLUG,
        "target_leiden_coarse": str(config["label"]),
        "n_pcs": int(config["n_pcs"]),
        "n_neighbors": int(config["n_neighbors"]),
        "resolution": float(config["resolution"]),
        "algorithm": "leiden",
        "raw_cluster_key": raw_key,
        "seed": SEED,
        "manual_selection_confirmed": True,
        "selection_note": SELECTION_NOTE,
    }
    release_gpu()

    paths["table_dir"].mkdir(parents=True, exist_ok=True)
    paths["h5ad_dir"].mkdir(parents=True, exist_ok=True)
    figure = write_selected_umap(source, raw_key, paths["figure_dir"])
    mean_dominance, max_dominance = sample_mixing_summary(source.obs, raw_key)

    cluster_sizes = (
        source.obs[raw_key]
        .value_counts(sort=False)
        .rename_axis(raw_key)
        .reset_index(name="n_cells")
    )
    cluster_sizes["fraction_of_lineage"] = cluster_sizes["n_cells"] / source.n_obs
    cluster_sizes.to_csv(paths["cluster_sizes"], index=False)
    sample_counts = (
        pd.crosstab(source.obs[raw_key].astype(str), source.obs["sample"].astype(str))
        .stack()
        .rename("n_cells")
        .reset_index()
    )
    sample_counts.to_csv(paths["sample_counts"], index=False)

    parameters = {
        "lineage": LINEAGE_SLUG,
        "target_leiden_coarse": config["label"],
        "cell_abbrev": config["abbrev"],
        "source_h5ad": str(paths["input_h5ad"]),
        "selected_h5ad": str(paths["output_h5ad"]),
        "grid_manifest": str(paths["grid_manifest"]),
        "grid_candidate_label": grid_row["candidate_label"],
        "n_cells": int(source.n_obs),
        "n_pcs": int(config["n_pcs"]),
        "n_neighbors": int(config["n_neighbors"]),
        "resolution": float(config["resolution"]),
        "algorithm": "leiden",
        "raw_cluster_key": raw_key,
        "cluster_count": observed_clusters,
        "grid_expected_cluster_count": int(config["expected_clusters"]),
        "mean_cluster_max_sample_fraction": mean_dominance,
        "max_cluster_max_sample_fraction": max_dominance,
        "use_rep": INTEGRATED_BASIS,
        "umap_min_dist": UMAP_MIN_DIST,
        "umap_spread": UMAP_SPREAD,
        "neighbors_algorithm": "brute",
        "neighbors_metric": "euclidean",
        "leiden_n_iterations": 100,
        "backend": "rapids_singlecell_gpu",
        "seed": SEED,
        "manual_selection_confirmed": True,
        "manual_selection_note": SELECTION_NOTE,
        "selected_figure": str(figure),
        "cluster_sizes_table": str(paths["cluster_sizes"]),
        "sample_counts_table": str(paths["sample_counts"]),
        "code_file": str(CODE_PATH),
        "concurrency_limit": 2,
        "gpu_assignment": os.environ.get("CUDA_VISIBLE_DEVICES", "0"),
        "gpu_safety_margin_mb": 4096,
        "neighbors_seconds": neighbors_seconds,
        "umap_seconds": umap_seconds,
        "clustering_seconds": clustering_seconds,
        **gpu,
    }
    pd.DataFrame([parameters]).to_csv(paths["parameters"], index=False)

    source.write_h5ad(paths["output_h5ad"], compression="gzip")
    saved = ad.read_h5ad(paths["output_h5ad"], backed="r")
    try:
        if raw_key not in saved.obs or "X_umap" not in saved.obsm:
            raise ValueError(f"Saved selected h5ad failed validation: {paths['output_h5ad']}")
        if int(saved.obs[raw_key].nunique()) != observed_clusters:
            raise ValueError(f"Saved selected h5ad cluster count mismatch: {paths['output_h5ad']}")
    finally:
        saved.file.close()

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
        "CUDA_VISIBLE_DEVICES": os.environ.get("CUDA_VISIBLE_DEVICES", ""),
        "code": str(CODE_PATH),
        "seed": str(SEED),
    }
    paths["versions"].write_text(
        "\n".join(f"{key}={value}" for key, value in versions.items()) + "\n",
        encoding="utf-8",
    )
    paths["readme"].write_text(
        f"""BRCA {config['label']} selected subtype clustering

Source: {paths['input_h5ad']}
Grid manifest: {paths['grid_manifest']}
Code: {CODE_PATH}
Manual choice: n_pcs={config['n_pcs']}; n_neighbors={config['n_neighbors']}; resolution={config['resolution']}
Raw cluster key: {raw_key}; clusters={observed_clusters}
Selected h5ad: {paths['output_h5ad']}
Selected UMAP: {figure}
All {source.n_obs} strict-consistent lineage cells were retained. The selected
neighbors, UMAP, and Leiden result were rerun from a clean read of the saved
lineage Harmony h5ad. No cell_subtype labels or DEG results are assigned here.
""",
        encoding="utf-8",
    )
    completion = {
        "lineage": LINEAGE_SLUG,
        "status": "completed",
        "n_cells": int(source.n_obs),
        "selected_h5ad": str(paths["output_h5ad"]),
        "selected_figure": str(figure),
        "raw_cluster_key": raw_key,
        "cluster_count": observed_clusters,
        "manual_selection_confirmed": True,
        "deg_annotation_started": False,
        "elapsed_seconds": time.time() - started,
    }
    paths["completion"].write_text(
        json.dumps(completion, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(completion, indent=2), flush=True)


if __name__ == "__main__":
    main()
