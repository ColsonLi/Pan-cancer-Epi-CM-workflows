#!/usr/bin/env python3
"""Run the user-specified BRCA Harmony clustering grid on all retained cells."""

from __future__ import annotations

import gc
import importlib.metadata
import json
import os
import random
import sys
import time
from pathlib import Path

# numba-cuda's legacy ctypes driver binding segfaults in cuCtxGetDevice on the
# host's newer NVIDIA driver when cudf transfers Leiden labels to host memory.
# The supported NVIDIA Python binding was verified on both a minimal cugraph
# and the full 99,866-cell graph and must be selected before importing RAPIDS.
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

ANALYSIS_ROOT = Path("/mnt/disk18t/lr_xcy/riku/brca_cm-epi_codex")
WORKFLOW_ROOT = ANALYSIS_ROOT / "epi-cm-core-workflow"
BLOCK = "01-celltype_integration_clustering"
INPUT_H5AD = (
    WORKFLOW_ROOT
    / "h5ad"
    / BLOCK
    / "04-integration-harmony"
    / "adata_harmony.h5ad"
)
TABLE_ROOT = WORKFLOW_ROOT / "tables" / BLOCK / "05-clustering-parameter-search"
FIGURE_ROOT = WORKFLOW_ROOT / "figures" / BLOCK / "05-clustering-parameter-search"
CODE_PATH = (
    WORKFLOW_ROOT
    / "codes"
    / BLOCK
    / "05-clustering-parameter-search"
    / "04_run_full_grid.py"
)
MANIFEST = TABLE_ROOT / "clustering_grid_manifest.csv"
REVIEW = TABLE_ROOT / "candidate_review_manifest.csv"
COMPLETION = TABLE_ROOT / "clustering_grid_completion_check.csv"

PCS_VALUES = list(range(10, 51, 5))
NN_VALUES = list(range(10, 51, 5))
RESOLUTIONS = [round(x / 10, 1) for x in range(1, 11)]
MAX_CLUSTER_COUNT = 20
INTEGRATED_BASIS = "X_pca_inte"
ALGORITHM = "leiden"
UMAP_MIN_DIST = 0.5
UMAP_SPREAD = 1.0


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def res_token(resolution: float) -> str:
    return f"{resolution:.1f}".replace(".", "p")


def graph_label(n_pcs: int, n_neighbors: int) -> str:
    return f"pcs-{n_pcs}_nn-{n_neighbors}_res-0p1-1p0"


def candidate_label(n_pcs: int, n_neighbors: int, resolution: float) -> str:
    return f"pcs-{n_pcs}_nn-{n_neighbors}_res-{res_token(resolution)}"


def cluster_key(resolution: float) -> str:
    return f"leiden_res{res_token(resolution)}"


def release_gpu() -> None:
    cp.get_default_memory_pool().free_all_blocks()
    cp.get_default_pinned_memory_pool().free_all_blocks()
    gc.collect()


def initial_manifest() -> pd.DataFrame:
    rows = []
    for n_pcs in PCS_VALUES:
        for n_neighbors in NN_VALUES:
            label = graph_label(n_pcs, n_neighbors)
            figure_dir = FIGURE_ROOT / label
            table_dir = TABLE_ROOT / label
            for resolution in RESOLUTIONS:
                rows.append(
                    {
                        "n_pcs": n_pcs,
                        "n_neighbors": n_neighbors,
                        "resolution": resolution,
                        "candidate_label": candidate_label(n_pcs, n_neighbors, resolution),
                        "graph_label": label,
                        "expected_cluster_count_table": str(table_dir / "cluster_counts.csv"),
                        "expected_parameter_table": str(table_dir / "clustering_parameters.csv"),
                        "expected_figure_dir": str(figure_dir),
                        "expected_figure": str(figure_dir / "umap_leiden_grid.pdf"),
                        "status": "planned",
                        "completed": False,
                        "reason_if_skipped": "",
                        "cluster_count": np.nan,
                    }
                )
    return pd.DataFrame(rows)


def load_or_initialize_manifest() -> pd.DataFrame:
    if MANIFEST.exists():
        manifest = pd.read_csv(MANIFEST)
        expected = len(PCS_VALUES) * len(NN_VALUES) * len(RESOLUTIONS)
        if len(manifest) != expected:
            raise ValueError(
                f"Existing manifest has {len(manifest)} rows; expected {expected}."
            )
        return manifest
    manifest = initial_manifest()
    TABLE_ROOT.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(MANIFEST, index=False)
    return manifest


def restore_completed_graphs(manifest: pd.DataFrame) -> tuple[pd.DataFrame, set[str]]:
    existing: set[str] = set()
    for n_pcs in PCS_VALUES:
        for n_neighbors in NN_VALUES:
            label = graph_label(n_pcs, n_neighbors)
            table_dir = TABLE_ROOT / label
            figure = FIGURE_ROOT / label / "umap_leiden_grid.pdf"
            counts_csv = table_dir / "cluster_counts.csv"
            params_csv = table_dir / "clustering_parameters.csv"
            if not (figure.exists() and counts_csv.exists() and params_csv.exists()):
                continue
            counts = pd.read_csv(counts_csv)
            required = {"resolution", "status", "cluster_count"}
            if len(counts) != len(RESOLUTIONS) or not required.issubset(counts.columns):
                continue
            for row in counts.itertuples(index=False):
                mask = (
                    manifest["n_pcs"].eq(n_pcs)
                    & manifest["n_neighbors"].eq(n_neighbors)
                    & np.isclose(manifest["resolution"].astype(float), float(row.resolution))
                )
                manifest.loc[mask, "status"] = str(row.status)
                manifest.loc[mask, "completed"] = str(row.status) == "completed"
                manifest.loc[mask, "cluster_count"] = row.cluster_count
                if str(row.status) == "skipped_user_approved":
                    manifest.loc[mask, "reason_if_skipped"] = (
                        "higher resolution after first cluster_count > 20 for the same graph"
                    )
            existing.add(label)
    return manifest, existing


def sample_mixing_summary(obs: pd.DataFrame, key: str) -> tuple[float, float]:
    fractions = pd.crosstab(obs[key].astype(str), obs["sample"].astype(str), normalize="index")
    dominance = fractions.max(axis=1)
    return float(dominance.mean()), float(dominance.max())


def write_combined_umap(plot_adata: ad.AnnData, keys: list[str], figure_dir: Path) -> Path:
    figure_dir.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = figure_dir
    sc.settings.file_format_figs = "pdf"
    sc.set_figure_params(figsize=(3, 3), dpi=150, fontsize=8)
    save_name = "_leiden_grid.pdf"
    if len(keys) == 1:
        sc.pl.umap(plot_adata, color=keys, show=False, save=save_name)
    else:
        sc.pl.umap(
            plot_adata,
            color=keys,
            ncols=min(5, len(keys)),
            wspace=0.35,
            show=False,
            save=save_name,
        )
    out = figure_dir / "umap_leiden_grid.pdf"
    if not out.exists():
        raise FileNotFoundError(f"Scanpy did not create expected figure: {out}")
    return out


def run_graph(template: ad.AnnData, n_pcs: int, n_neighbors: int) -> tuple[pd.DataFrame, dict[str, object]]:
    label = graph_label(n_pcs, n_neighbors)
    table_dir = TABLE_ROOT / label
    figure_dir = FIGURE_ROOT / label
    table_dir.mkdir(parents=True, exist_ok=True)

    graph_start = time.time()
    print(
        f"[Grid] start {label} cells={template.n_obs}",
        flush=True,
    )
    adata_run = template.copy()
    step_start = time.time()
    rsc.pp.neighbors(
        adata_run,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        use_rep=INTEGRATED_BASIS,
        random_state=SEED,
        algorithm="brute",
        metric="euclidean",
    )
    neighbors_seconds = time.time() - step_start
    step_start = time.time()
    rsc.tl.umap(
        adata_run,
        min_dist=UMAP_MIN_DIST,
        spread=UMAP_SPREAD,
        random_state=SEED,
    )
    umap_seconds = time.time() - step_start
    print(
        f"[Grid] {label} graph+UMAP completed "
        f"neighbors={neighbors_seconds:.1f}s umap={umap_seconds:.1f}s",
        flush=True,
    )
    rows: list[dict[str, object]] = []
    executed_keys: list[str] = []
    cutoff_resolution: float | None = None
    for resolution in RESOLUTIONS:
        key = cluster_key(resolution)
        if cutoff_resolution is not None:
            rows.append(
                {
                    "n_pcs": n_pcs,
                    "n_neighbors": n_neighbors,
                    "resolution": resolution,
                    "raw_cluster_key": key,
                    "algorithm": ALGORITHM,
                    "cluster_count": np.nan,
                    "status": "skipped_user_approved",
                    "reason": (
                        f"higher resolution after res={cutoff_resolution:.1f} "
                        f"produced >{MAX_CLUSTER_COUNT} clusters"
                    ),
                    "mean_cluster_max_sample_fraction": np.nan,
                    "max_cluster_max_sample_fraction": np.nan,
                }
            )
            continue
        rsc.tl.leiden(
            adata_run,
            resolution=resolution,
            random_state=SEED,
            key_added=key,
            n_iterations=100,
        )
        n_clusters = int(adata_run.obs[key].nunique())
        print(
            f"[Grid] {label} res={resolution:.1f} clusters={n_clusters}",
            flush=True,
        )
        mean_dom, max_dom = sample_mixing_summary(adata_run.obs, key)
        rows.append(
            {
                "n_pcs": n_pcs,
                "n_neighbors": n_neighbors,
                "resolution": resolution,
                "raw_cluster_key": key,
                "algorithm": ALGORITHM,
                "cluster_count": n_clusters,
                "status": "completed",
                "reason": "",
                "mean_cluster_max_sample_fraction": mean_dom,
                "max_cluster_max_sample_fraction": max_dom,
            }
        )
        executed_keys.append(key)
        if n_clusters > MAX_CLUSTER_COUNT:
            cutoff_resolution = resolution
            print(
                f"[Grid] {label} cutoff reached at res={resolution:.1f}; "
                "higher resolutions will be skipped",
                flush=True,
            )

    counts = pd.DataFrame(rows)
    counts.to_csv(table_dir / "cluster_counts.csv", index=False)

    umap = adata_run.obsm["X_umap"]
    if isinstance(umap, cp.ndarray):
        umap = cp.asnumpy(umap)
    plot_obs = adata_run.obs[executed_keys + ["sample", "series", "status"]].copy()
    plot_adata = ad.AnnData(
        X=sparse.csr_matrix((adata_run.n_obs, 0), dtype=np.float32),
        obs=plot_obs,
    )
    plot_adata.obsm["X_umap"] = np.asarray(umap, dtype=np.float32)
    figure = write_combined_umap(plot_adata, executed_keys, figure_dir)

    params = counts.copy()
    params["source_h5ad"] = str(INPUT_H5AD)
    params["integration_source"] = "Harmony"
    params["use_rep"] = INTEGRATED_BASIS
    params["source_template_transferred_to_gpu_before_loop"] = True
    params["candidate_object_created_by_copy"] = True
    params["graph_template_expression_omitted"] = True
    params["graph_template_all_cells_preserved"] = True
    params["umap_min_dist"] = UMAP_MIN_DIST
    params["umap_spread"] = UMAP_SPREAD
    params["figure"] = str(figure)
    params["cluster_count_table"] = str(table_dir / "cluster_counts.csv")
    params["code_file"] = str(CODE_PATH)
    params["seed"] = SEED
    params["backend"] = "rapids_singlecell_gpu"
    params["candidate_h5ad_saved"] = False
    params["candidate_object_deleted_after_outputs"] = True
    params["gpu_memory_release_after_candidate"] = True
    params["user_cluster_cutoff_rule"] = (
        "execute and draw the first resolution with cluster_count > 20; "
        "skip all higher resolutions for the same PC/NN graph"
    )
    params.to_csv(table_dir / "clustering_parameters.csv", index=False)
    (table_dir / "readme.txt").write_text(
        f"""BRCA clustering grid candidate graph

Source: {INPUT_H5AD}
Code: {CODE_PATH}
n_pcs={n_pcs}; n_neighbors={n_neighbors}; resolutions=0.1..1.0 step 0.1
All retained cells are used. The first resolution with more than 20 clusters is
executed and drawn; later resolutions for this graph are skipped by explicit
user instruction. All executed resolutions share one neighbor graph and UMAP.
No candidate h5ad is saved.
""",
        encoding="utf-8",
    )
    review = {
        "n_pcs": n_pcs,
        "n_neighbors": n_neighbors,
        "graph_label": label,
        "n_executed_resolutions": int(counts["status"].eq("completed").sum()),
        "n_skipped_resolutions": int(counts["status"].eq("skipped_user_approved").sum()),
        "first_resolution_over_20_clusters": cutoff_resolution,
        "max_executed_resolution": float(
            counts.loc[counts["status"].eq("completed"), "resolution"].max()
        ),
        "max_executed_cluster_count": int(
            counts.loc[counts["status"].eq("completed"), "cluster_count"].max()
        ),
        "cluster_counts": ";".join(
            f"{row.resolution:.1f}:{int(row.cluster_count)}"
            for row in counts.loc[counts["status"].eq("completed")].itertuples(index=False)
        ),
        "figure": str(figure),
        "cluster_count_table": str(table_dir / "cluster_counts.csv"),
        "elapsed_seconds": time.time() - graph_start,
    }
    del plot_adata, adata_run
    release_gpu()
    print(
        f"[Grid] completed {label} executed={review['n_executed_resolutions']} "
        f"skipped={review['n_skipped_resolutions']} "
        f"elapsed={review['elapsed_seconds']:.1f}s",
        flush=True,
    )
    return counts, review


def update_manifest_for_graph(
    manifest: pd.DataFrame, n_pcs: int, n_neighbors: int, counts: pd.DataFrame
) -> pd.DataFrame:
    for row in counts.itertuples(index=False):
        mask = (
            manifest["n_pcs"].eq(n_pcs)
            & manifest["n_neighbors"].eq(n_neighbors)
            & np.isclose(manifest["resolution"].astype(float), float(row.resolution))
        )
        manifest.loc[mask, "status"] = row.status
        manifest.loc[mask, "completed"] = row.status == "completed"
        manifest.loc[mask, "cluster_count"] = row.cluster_count
        manifest.loc[mask, "reason_if_skipped"] = row.reason
    manifest.to_csv(MANIFEST, index=False)
    return manifest


def finalize(manifest: pd.DataFrame, review_rows: list[dict[str, object]]) -> None:
    review = pd.DataFrame(review_rows).sort_values(["n_pcs", "n_neighbors"])
    review.to_csv(REVIEW, index=False)
    terminal = manifest["status"].isin(["completed", "skipped_user_approved"])
    expected_graphs = len(PCS_VALUES) * len(NN_VALUES)
    completed_graphs = int(review["graph_label"].nunique())
    all_outputs = True
    for row in review.itertuples(index=False):
        all_outputs &= Path(row.figure).exists()
        all_outputs &= Path(row.cluster_count_table).exists()
    check = pd.DataFrame(
        [
            {
                "expected_graph_candidates": expected_graphs,
                "completed_graph_candidates": completed_graphs,
                "expected_resolution_candidates": len(manifest),
                "completed_resolution_candidates": int(manifest["status"].eq("completed").sum()),
                "skipped_resolution_candidates": int(
                    manifest["status"].eq("skipped_user_approved").sum()
                ),
                "failed_candidates": int(manifest["status"].eq("failed").sum()),
                "missing_candidates": int((~terminal).sum()),
                "all_expected_lightweight_outputs_exist": bool(all_outputs),
                "manual_selection_ready": bool(
                    completed_graphs == expected_graphs
                    and terminal.all()
                    and all_outputs
                ),
                "dynamic_cutoff": "first cluster_count>20 is drawn; later resolutions skipped",
            }
        ]
    )
    check.to_csv(COMPLETION, index=False)


def main() -> None:
    start = time.time()
    manifest = load_or_initialize_manifest()
    manifest, existing_graphs = restore_completed_graphs(manifest)
    manifest.to_csv(MANIFEST, index=False)

    source = sc.read_h5ad(INPUT_H5AD)
    if INTEGRATED_BASIS not in source.obsm:
        raise ValueError(f"Missing {INTEGRATED_BASIS} in {INPUT_H5AD}")
    if source.obsm[INTEGRATED_BASIS].shape[1] < max(PCS_VALUES):
        raise ValueError(
            f"Harmony basis has only {source.obsm[INTEGRATED_BASIS].shape[1]} dimensions; "
            f"grid requires {max(PCS_VALUES)}."
        )
    if source.n_obs == 0 or not source.obs_names.is_unique:
        raise ValueError("Harmony source has no cells or non-unique cell IDs.")

    # The clustering methods use only X_pca_inte. Build a graph-only template
    # that preserves every cell and all required obs/embedding values while
    # omitting unused expression arrays from candidate copies. This is not cell
    # downsampling and does not alter the Harmony source on disk.
    template = ad.AnnData(
        X=sparse.csr_matrix((source.n_obs, 0), dtype=np.float32),
        obs=source.obs[["sample", "series", "status"]].copy(),
    )
    template.obsm[INTEGRATED_BASIS] = np.asarray(
        source.obsm[INTEGRATED_BASIS], dtype=np.float32
    )
    del source
    rsc.get.anndata_to_GPU(template)
    print(
        f"[Grid] source ready: cells={template.n_obs}, "
        f"graphs={len(PCS_VALUES) * len(NN_VALUES)}, "
        f"resolutions_per_graph={len(RESOLUTIONS)}, "
        f"restored_complete_graphs={len(existing_graphs)}",
        flush=True,
    )

    review_by_label: dict[str, dict[str, object]] = {}
    if REVIEW.exists():
        for row in pd.read_csv(REVIEW).to_dict(orient="records"):
            review_by_label[str(row["graph_label"])] = row

    for n_pcs in PCS_VALUES:
        for n_neighbors in NN_VALUES:
            label = graph_label(n_pcs, n_neighbors)
            if label in existing_graphs:
                print(f"[Grid] existing valid outputs skipped: {label}", flush=True)
                counts = pd.read_csv(TABLE_ROOT / label / "cluster_counts.csv")
                if label not in review_by_label:
                    completed = counts[counts["status"].eq("completed")]
                    cutoff = completed.loc[
                        completed["cluster_count"].astype(float) > MAX_CLUSTER_COUNT,
                        "resolution",
                    ]
                    review_by_label[label] = {
                        "n_pcs": n_pcs,
                        "n_neighbors": n_neighbors,
                        "graph_label": label,
                        "n_executed_resolutions": len(completed),
                        "n_skipped_resolutions": int(
                            counts["status"].eq("skipped_user_approved").sum()
                        ),
                        "first_resolution_over_20_clusters": (
                            float(cutoff.iloc[0]) if len(cutoff) else np.nan
                        ),
                        "max_executed_resolution": float(completed["resolution"].max()),
                        "max_executed_cluster_count": int(completed["cluster_count"].max()),
                        "cluster_counts": ";".join(
                            f"{row.resolution:.1f}:{int(row.cluster_count)}"
                            for row in completed.itertuples(index=False)
                        ),
                        "figure": str(FIGURE_ROOT / label / "umap_leiden_grid.pdf"),
                        "cluster_count_table": str(TABLE_ROOT / label / "cluster_counts.csv"),
                        "elapsed_seconds": np.nan,
                    }
                continue
            mask = manifest["graph_label"].eq(label)
            manifest.loc[mask, "status"] = "running"
            manifest.to_csv(MANIFEST, index=False)
            try:
                counts, review = run_graph(template, n_pcs, n_neighbors)
            except Exception as exc:
                manifest.loc[mask, "status"] = "failed"
                manifest.loc[mask, "reason_if_skipped"] = repr(exc)
                manifest.to_csv(MANIFEST, index=False)
                raise
            manifest = update_manifest_for_graph(
                manifest, n_pcs, n_neighbors, counts
            )
            review_by_label[label] = review
            finalize(manifest, list(review_by_label.values()))

    release_gpu()
    finalize(manifest, list(review_by_label.values()))
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
    (TABLE_ROOT / "package_versions.txt").write_text(
        "\n".join(f"{k}={v}" for k, v in versions.items()) + "\n",
        encoding="utf-8",
    )
    (TABLE_ROOT / "readme.txt").write_text(
        f"""BRCA Harmony clustering parameter grid

Source: {INPUT_H5AD}
Code: {CODE_PATH}
Grid: PCs 10..50 step 5; neighbors 10..50 step 5; resolution 0.1..1.0 step 0.1.
For each PC/NN graph, the first resolution producing more than 20 clusters is
executed and drawn, and higher resolutions are skipped by explicit user request.
Already complete aggregate figures/tables are validated and skipped on resume.
Every candidate uses all retained cells, seed {SEED}, and {INTEGRATED_BASIS}.
No candidate h5ad is saved and no parameter is automatically selected.
""",
        encoding="utf-8",
    )
    completion = pd.read_csv(COMPLETION).iloc[0].to_dict()
    completion["elapsed_seconds_this_invocation"] = time.time() - start
    print(json.dumps(completion, indent=2, default=str))


if __name__ == "__main__":
    main()
