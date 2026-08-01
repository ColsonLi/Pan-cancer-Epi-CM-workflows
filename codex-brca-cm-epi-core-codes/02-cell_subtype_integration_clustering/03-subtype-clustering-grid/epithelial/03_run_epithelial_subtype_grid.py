#!/usr/bin/env python3
"""Run the complete BRCA epithelial subtype clustering grid on Harmony PCs."""

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
INPUT_H5AD = (
    WORKFLOW
    / "h5ad"
    / BLOCK
    / "02-subtype-harmony"
    / "epithelial"
    / "adata_epithelial_harmony.h5ad"
)
TABLE_ROOT = WORKFLOW / "tables" / BLOCK / "03-subtype-clustering-grid" / "epithelial"
FIGURE_ROOT = WORKFLOW / "figures" / BLOCK / "03-subtype-clustering-grid" / "epithelial"
CODE_PATH = (
    WORKFLOW
    / "codes"
    / BLOCK
    / "03-subtype-clustering-grid"
    / "epithelial"
    / "03_run_epithelial_subtype_grid.py"
)
MANIFEST = TABLE_ROOT / "clustering_grid_manifest.csv"
REVIEW = TABLE_ROOT / "candidate_review_manifest.csv"
COMPLETION = TABLE_ROOT / "clustering_grid_completion_check.csv"
LINEAGE_SLUG = "epithelial"
TARGET_LABEL = "Epithelial Cells"
SELECTED_IDS = (
    WORKFLOW
    / "tables"
    / BLOCK
    / "01-lineage-selection"
    / LINEAGE_SLUG
    / "selected_cell_ids.csv"
)
LINEAGE_LABELS = {
    "epithelial": "Epithelial Cells",
    "t_cells": "T Cells",
    "myeloid": "Myeloid Cells",
    "b_cells": "B Cells",
    "plasma": "Plasma Cells",
    "endothelial": "Endothelial Cells",
    "stromal": "Stromal Cells",
    "perivascular": "Perivascular Cells",
}

PCS_VALUES = list(range(10, 51, 5))
NN_VALUES = list(range(10, 51, 5))
RESOLUTIONS = [round(value / 10, 1) for value in range(1, 16)]
INTEGRATED_BASIS = "X_pca_inte"
UMAP_MIN_DIST = 0.5
UMAP_SPREAD = 1.0
ABNORMAL_CLUSTER_COUNT = 100


def configure_lineage(lineage_slug: str, code_path: Path | None = None) -> None:
    """Configure this validated grid implementation for another broad lineage."""
    global LINEAGE_SLUG, TARGET_LABEL, INPUT_H5AD, SELECTED_IDS
    global TABLE_ROOT, FIGURE_ROOT, CODE_PATH, MANIFEST, REVIEW, COMPLETION
    if lineage_slug not in LINEAGE_LABELS:
        raise ValueError(f"Unsupported lineage slug: {lineage_slug}")
    LINEAGE_SLUG = lineage_slug
    TARGET_LABEL = LINEAGE_LABELS[lineage_slug]
    INPUT_H5AD = (
        WORKFLOW
        / "h5ad"
        / BLOCK
        / "02-subtype-harmony"
        / lineage_slug
        / f"adata_{lineage_slug}_harmony.h5ad"
    )
    SELECTED_IDS = (
        WORKFLOW
        / "tables"
        / BLOCK
        / "01-lineage-selection"
        / lineage_slug
        / "selected_cell_ids.csv"
    )
    TABLE_ROOT = WORKFLOW / "tables" / BLOCK / "03-subtype-clustering-grid" / lineage_slug
    FIGURE_ROOT = WORKFLOW / "figures" / BLOCK / "03-subtype-clustering-grid" / lineage_slug
    MANIFEST = TABLE_ROOT / "clustering_grid_manifest.csv"
    REVIEW = TABLE_ROOT / "candidate_review_manifest.csv"
    COMPLETION = TABLE_ROOT / "clustering_grid_completion_check.csv"
    if code_path is not None:
        CODE_PATH = code_path


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def res_token(resolution: float) -> str:
    return f"{resolution:.1f}".replace(".", "p")


def graph_label(n_pcs: int, n_neighbors: int) -> str:
    return f"pcs-{n_pcs}_nn-{n_neighbors}_res-0p1-1p5"


def candidate_label(n_pcs: int, n_neighbors: int, resolution: float) -> str:
    return f"pcs-{n_pcs}_nn-{n_neighbors}_res-{res_token(resolution)}"


def cluster_key(algorithm: str, resolution: float) -> str:
    return f"{algorithm}_res{res_token(resolution)}"


def release_gpu() -> None:
    cp.get_default_memory_pool().free_all_blocks()
    cp.get_default_pinned_memory_pool().free_all_blocks()
    gc.collect()


def initial_manifest() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for n_pcs in PCS_VALUES:
        for n_neighbors in NN_VALUES:
            graph = graph_label(n_pcs, n_neighbors)
            table_dir = TABLE_ROOT / graph
            figure_dir = FIGURE_ROOT / graph
            for resolution in RESOLUTIONS:
                rows.append(
                    {
                        "n_pcs": n_pcs,
                        "n_neighbors": n_neighbors,
                        "resolution": resolution,
                        "candidate_label": candidate_label(n_pcs, n_neighbors, resolution),
                        "graph_label": graph,
                        "expected_cluster_count_table": str(table_dir / "cluster_counts.csv"),
                        "expected_parameter_table": str(table_dir / "clustering_parameters.csv"),
                        "expected_figure_dir": str(figure_dir),
                        "expected_figure": str(figure_dir / "umap_subtype_grid.pdf"),
                        "status": "planned",
                        "completed": False,
                        "reason_if_failed": "",
                        "algorithm": "leiden_planned",
                        "raw_cluster_key": "",
                        "cluster_count": np.nan,
                    }
                )
    return pd.DataFrame(rows)


def load_or_initialize_manifest() -> pd.DataFrame:
    expected = len(PCS_VALUES) * len(NN_VALUES) * len(RESOLUTIONS)
    if MANIFEST.exists():
        manifest = pd.read_csv(MANIFEST)
        if len(manifest) != expected:
            raise ValueError(f"Existing manifest has {len(manifest)} rows; expected {expected}.")
        for column in ["reason_if_failed", "raw_cluster_key", "algorithm"]:
            manifest[column] = manifest[column].fillna("").astype(str)
        return manifest
    TABLE_ROOT.mkdir(parents=True, exist_ok=True)
    manifest = initial_manifest()
    manifest.to_csv(MANIFEST, index=False)
    return manifest


def pdf_is_readable(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    result = subprocess.run(["pdfinfo", str(path)], capture_output=True, text=True)
    return result.returncode == 0 and "Pages:" in result.stdout


def restore_completed_graphs(manifest: pd.DataFrame) -> tuple[pd.DataFrame, set[str]]:
    existing: set[str] = set()
    for n_pcs in PCS_VALUES:
        for n_neighbors in NN_VALUES:
            graph = graph_label(n_pcs, n_neighbors)
            table_dir = TABLE_ROOT / graph
            figure = FIGURE_ROOT / graph / "umap_subtype_grid.pdf"
            counts_csv = table_dir / "cluster_counts.csv"
            params_csv = table_dir / "clustering_parameters.csv"
            if not (pdf_is_readable(figure) and counts_csv.exists() and params_csv.exists()):
                continue
            counts = pd.read_csv(counts_csv)
            required = {"resolution", "status", "cluster_count", "algorithm", "raw_cluster_key"}
            observed_resolutions = sorted(counts["resolution"].astype(float).round(1).tolist())
            if (
                len(counts) != len(RESOLUTIONS)
                or not required.issubset(counts.columns)
                or observed_resolutions != RESOLUTIONS
                or not counts["status"].eq("completed").all()
            ):
                continue
            for row in counts.itertuples(index=False):
                mask = (
                    manifest["n_pcs"].eq(n_pcs)
                    & manifest["n_neighbors"].eq(n_neighbors)
                    & np.isclose(manifest["resolution"].astype(float), float(row.resolution))
                )
                manifest.loc[mask, "status"] = "completed"
                manifest.loc[mask, "completed"] = True
                manifest.loc[mask, "cluster_count"] = int(row.cluster_count)
                manifest.loc[mask, "algorithm"] = str(row.algorithm)
                manifest.loc[mask, "raw_cluster_key"] = str(row.raw_cluster_key)
                manifest.loc[mask, "reason_if_failed"] = ""
            existing.add(graph)
    return manifest, existing


def sample_mixing_summary(obs: pd.DataFrame, key: str) -> tuple[float, float]:
    fractions = pd.crosstab(
        obs[key].astype(str), obs["sample"].astype(str), normalize="index"
    )
    dominance = fractions.max(axis=1)
    return float(dominance.mean()), float(dominance.max())


def write_combined_umap(
    plot_adata: ad.AnnData,
    diagnostic_keys: list[str],
    cluster_keys: list[str],
    figure_dir: Path,
) -> Path:
    figure_dir.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = figure_dir
    sc.settings.file_format_figs = "pdf"
    sc.set_figure_params(figsize=(3, 3), dpi=150, fontsize=8)
    keys = diagnostic_keys + cluster_keys
    sc.pl.umap(
        plot_adata,
        color=keys,
        ncols=5,
        wspace=0.4,
        show=False,
        save="_subtype_grid.pdf",
    )
    output = figure_dir / "umap_subtype_grid.pdf"
    if not pdf_is_readable(output):
        raise FileNotFoundError(f"Scanpy did not create a readable figure: {output}")
    return output


def run_graph(
    template: ad.AnnData, n_pcs: int, n_neighbors: int
) -> tuple[pd.DataFrame, dict[str, object]]:
    graph = graph_label(n_pcs, n_neighbors)
    table_dir = TABLE_ROOT / graph
    figure_dir = FIGURE_ROOT / graph
    table_dir.mkdir(parents=True, exist_ok=True)
    graph_started = time.time()
    print(f"[{LINEAGE_SLUG} grid] start {graph} cells={template.n_obs}", flush=True)
    adata_run = template.copy()

    t0 = time.time()
    rsc.pp.neighbors(
        adata_run,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        use_rep=INTEGRATED_BASIS,
        random_state=SEED,
        algorithm="brute",
        metric="euclidean",
    )
    neighbors_seconds = time.time() - t0
    t0 = time.time()
    rsc.tl.umap(
        adata_run,
        min_dist=UMAP_MIN_DIST,
        spread=UMAP_SPREAD,
        random_state=SEED,
    )
    umap_seconds = time.time() - t0
    print(
        f"[{LINEAGE_SLUG} grid] {graph} neighbors={neighbors_seconds:.1f}s "
        f"umap={umap_seconds:.1f}s",
        flush=True,
    )

    rows: list[dict[str, object]] = []
    cluster_keys: list[str] = []
    for resolution in RESOLUTIONS:
        algorithm = "leiden"
        key = cluster_key(algorithm, resolution)
        rsc.tl.leiden(
            adata_run,
            resolution=resolution,
            random_state=SEED,
            key_added=key,
            n_iterations=100,
        )
        n_clusters = int(adata_run.obs[key].nunique())
        fallback_reason = ""
        if n_clusters > ABNORMAL_CLUSTER_COUNT:
            fallback_reason = (
                f"Leiden produced {n_clusters} clusters (> {ABNORMAL_CLUSTER_COUNT}); "
                "used required abnormal-Leiden Louvain fallback."
            )
            del adata_run.obs[key]
            algorithm = "louvain"
            key = cluster_key(algorithm, resolution)
            rsc.tl.louvain(
                adata_run,
                resolution=resolution,
                key_added=key,
                n_iterations=100,
            )
            n_clusters = int(adata_run.obs[key].nunique())
        mean_dominance, max_dominance = sample_mixing_summary(adata_run.obs, key)
        rows.append(
            {
                "n_pcs": n_pcs,
                "n_neighbors": n_neighbors,
                "resolution": resolution,
                "raw_cluster_key": key,
                "algorithm": algorithm,
                "cluster_count": n_clusters,
                "status": "completed",
                "fallback_reason": fallback_reason,
                "mean_cluster_max_sample_fraction": mean_dominance,
                "max_cluster_max_sample_fraction": max_dominance,
            }
        )
        cluster_keys.append(key)
        print(
            f"[{LINEAGE_SLUG} grid] {graph} res={resolution:.1f} "
            f"algorithm={algorithm} clusters={n_clusters}",
            flush=True,
        )

    counts = pd.DataFrame(rows)
    counts.to_csv(table_dir / "cluster_counts.csv", index=False)
    umap = adata_run.obsm["X_umap"]
    if isinstance(umap, cp.ndarray):
        umap = cp.asnumpy(umap)
    diagnostic_keys = ["sample"]
    for column in ["series", "status"]:
        if adata_run.obs[column].astype(str).nunique() > 1:
            diagnostic_keys.append(column)
    plot_columns = diagnostic_keys + cluster_keys
    plot_adata = ad.AnnData(
        X=sparse.csr_matrix((adata_run.n_obs, 0), dtype=np.float32),
        obs=adata_run.obs[plot_columns].copy(),
    )
    plot_adata.obsm["X_umap"] = np.asarray(umap, dtype=np.float32)
    figure = write_combined_umap(plot_adata, diagnostic_keys, cluster_keys, figure_dir)

    params = counts.copy()
    params["source_h5ad"] = str(INPUT_H5AD)
    params["lineage"] = LINEAGE_SLUG
    params["target_leiden_coarse"] = TARGET_LABEL
    params["integration_source"] = "Harmony"
    params["use_rep"] = INTEGRATED_BASIS
    params["source_template_transferred_to_gpu_before_loop"] = True
    params["candidate_object_created_by_copy"] = True
    params["graph_template_expression_omitted"] = True
    params["graph_template_all_cells_preserved"] = True
    params["n_cells"] = int(adata_run.n_obs)
    params["umap_min_dist"] = UMAP_MIN_DIST
    params["umap_spread"] = UMAP_SPREAD
    params["diagnostic_umap_keys"] = ";".join(diagnostic_keys)
    params["figure"] = str(figure)
    params["cluster_count_table"] = str(table_dir / "cluster_counts.csv")
    params["code_file"] = str(CODE_PATH)
    params["seed"] = SEED
    params["backend"] = "rapids_singlecell_gpu"
    params["candidate_h5ad_saved"] = False
    params["candidate_object_deleted_after_outputs"] = True
    params["gpu_memory_release_after_candidate"] = True
    params["all_default_resolutions_executed"] = True
    params.to_csv(table_dir / "clustering_parameters.csv", index=False)
    (table_dir / "readme.txt").write_text(
        f"""BRCA {TARGET_LABEL} subtype clustering candidate graph

Source: {INPUT_H5AD}
Code: {CODE_PATH}
n_pcs={n_pcs}; n_neighbors={n_neighbors}; resolutions=0.1..1.5 step 0.1
All {adata_run.n_obs} strict-consistent {TARGET_LABEL} cells are used. All 15
resolutions are executed on one graph and one UMAP. The combined review PDF
contains the sample diagnostic and all raw clustering labels. No candidate h5ad
is saved.
""",
        encoding="utf-8",
    )
    review = {
        "n_pcs": n_pcs,
        "n_neighbors": n_neighbors,
        "graph_label": graph,
        "n_completed_resolutions": int(len(counts)),
        "min_cluster_count": int(counts["cluster_count"].min()),
        "max_cluster_count": int(counts["cluster_count"].max()),
        "n_louvain_fallbacks": int(counts["algorithm"].eq("louvain").sum()),
        "cluster_counts": ";".join(
            f"{row.resolution:.1f}:{row.algorithm}:{int(row.cluster_count)}"
            for row in counts.itertuples(index=False)
        ),
        "figure": str(figure),
        "cluster_count_table": str(table_dir / "cluster_counts.csv"),
        "elapsed_seconds": time.time() - graph_started,
    }
    del plot_adata, adata_run
    release_gpu()
    print(
        f"[{LINEAGE_SLUG} grid] completed {graph} resolutions={len(counts)} "
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
        manifest.loc[mask, "status"] = "completed"
        manifest.loc[mask, "completed"] = True
        manifest.loc[mask, "cluster_count"] = int(row.cluster_count)
        manifest.loc[mask, "algorithm"] = str(row.algorithm)
        manifest.loc[mask, "raw_cluster_key"] = str(row.raw_cluster_key)
        manifest.loc[mask, "reason_if_failed"] = ""
    manifest.to_csv(MANIFEST, index=False)
    return manifest


def finalize(manifest: pd.DataFrame, review_rows: list[dict[str, object]]) -> None:
    review = pd.DataFrame(review_rows).sort_values(["n_pcs", "n_neighbors"])
    review.to_csv(REVIEW, index=False)
    expected_graphs = len(PCS_VALUES) * len(NN_VALUES)
    expected_candidates = expected_graphs * len(RESOLUTIONS)
    completed_graphs = int(review["graph_label"].nunique()) if len(review) else 0
    all_outputs = len(review) == expected_graphs
    for row in review.itertuples(index=False):
        all_outputs &= pdf_is_readable(Path(row.figure))
        all_outputs &= Path(row.cluster_count_table).exists()
    completed_candidates = int(manifest["status"].eq("completed").sum())
    failed_candidates = int(manifest["status"].eq("failed").sum())
    missing_candidates = expected_candidates - completed_candidates
    check = pd.DataFrame(
        [
            {
                "expected_graph_candidates": expected_graphs,
                "completed_graph_candidates": completed_graphs,
                "expected_resolution_candidates": expected_candidates,
                "completed_resolution_candidates": completed_candidates,
                "failed_candidates": failed_candidates,
                "missing_candidates": missing_candidates,
                "all_expected_lightweight_outputs_exist": bool(all_outputs),
                "manual_selection_ready": bool(
                    completed_graphs == expected_graphs
                    and completed_candidates == expected_candidates
                    and failed_candidates == 0
                    and all_outputs
                ),
                "candidate_h5ad_files_saved": 0,
                "grid": "PCs 10..50 step 5; NN 10..50 step 5; res 0.1..1.5 step 0.1",
            }
        ]
    )
    check.to_csv(COMPLETION, index=False)


def main() -> None:
    started = time.time()
    manifest = load_or_initialize_manifest()
    manifest, existing_graphs = restore_completed_graphs(manifest)
    manifest.to_csv(MANIFEST, index=False)

    source = sc.read_h5ad(INPUT_H5AD)
    if INTEGRATED_BASIS not in source.obsm:
        raise ValueError(f"Missing {INTEGRATED_BASIS} in {INPUT_H5AD}")
    if source.obsm[INTEGRATED_BASIS].shape != (source.n_obs, 50):
        raise ValueError(f"Unexpected Harmony basis shape: {source.obsm[INTEGRATED_BASIS].shape}")
    selected_ids = pd.read_csv(SELECTED_IDS, usecols=["cell_id"], dtype={"cell_id": str})
    expected_cells = int(len(selected_ids))
    if source.n_obs != expected_cells or not source.obs_names.is_unique:
        raise ValueError(
            f"Unexpected {LINEAGE_SLUG} Harmony cells/IDs: "
            f"observed={source.n_obs}, expected={expected_cells}"
        )
    if not source.obs["leiden_coarse"].astype(str).eq(TARGET_LABEL).all():
        raise ValueError(f"Off-lineage cells found in {LINEAGE_SLUG} Harmony source.")
    required_obs = {"sample", "series", "status"}
    if not required_obs.issubset(source.obs.columns):
        raise ValueError(f"Harmony source lacks obs columns: {sorted(required_obs - set(source.obs))}")

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
        f"[{LINEAGE_SLUG} grid] source ready cells={template.n_obs} "
        f"graphs={len(PCS_VALUES) * len(NN_VALUES)} "
        f"resolutions_per_graph={len(RESOLUTIONS)} "
        f"restored_graphs={len(existing_graphs)}",
        flush=True,
    )

    review_by_graph: dict[str, dict[str, object]] = {}
    if REVIEW.exists():
        for row in pd.read_csv(REVIEW).to_dict(orient="records"):
            review_by_graph[str(row["graph_label"])] = row

    for n_pcs in PCS_VALUES:
        for n_neighbors in NN_VALUES:
            graph = graph_label(n_pcs, n_neighbors)
            if graph in existing_graphs:
                print(f"[{LINEAGE_SLUG} grid] valid existing outputs skipped: {graph}", flush=True)
                if graph not in review_by_graph:
                    counts = pd.read_csv(TABLE_ROOT / graph / "cluster_counts.csv")
                    review_by_graph[graph] = {
                        "n_pcs": n_pcs,
                        "n_neighbors": n_neighbors,
                        "graph_label": graph,
                        "n_completed_resolutions": int(len(counts)),
                        "min_cluster_count": int(counts["cluster_count"].min()),
                        "max_cluster_count": int(counts["cluster_count"].max()),
                        "n_louvain_fallbacks": int(counts["algorithm"].eq("louvain").sum()),
                        "cluster_counts": ";".join(
                            f"{row.resolution:.1f}:{row.algorithm}:{int(row.cluster_count)}"
                            for row in counts.itertuples(index=False)
                        ),
                        "figure": str(FIGURE_ROOT / graph / "umap_subtype_grid.pdf"),
                        "cluster_count_table": str(TABLE_ROOT / graph / "cluster_counts.csv"),
                        "elapsed_seconds": np.nan,
                    }
                continue
            mask = manifest["graph_label"].eq(graph)
            manifest.loc[mask, "status"] = "running"
            manifest.to_csv(MANIFEST, index=False)
            try:
                counts, review = run_graph(template, n_pcs, n_neighbors)
            except Exception as exc:
                manifest.loc[mask, "status"] = "failed"
                manifest.loc[mask, "completed"] = False
                manifest.loc[mask, "reason_if_failed"] = repr(exc)
                manifest.to_csv(MANIFEST, index=False)
                raise
            manifest = update_manifest_for_graph(manifest, n_pcs, n_neighbors, counts)
            review_by_graph[graph] = review
            finalize(manifest, list(review_by_graph.values()))

    release_gpu()
    finalize(manifest, list(review_by_graph.values()))
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
        "\n".join(f"{key}={value}" for key, value in versions.items()) + "\n",
        encoding="utf-8",
    )
    (TABLE_ROOT / "readme.txt").write_text(
        f"""BRCA {TARGET_LABEL} subtype Harmony clustering grid

Source: {INPUT_H5AD}
Code: {CODE_PATH}
Grid: PCs 10..50 step 5; neighbors 10..50 step 5; resolution 0.1..1.5 step 0.1.
All 81 graph candidates and all 1215 resolution candidates use all {expected_cells}
strict-consistent {TARGET_LABEL} cells, seed {SEED}, and {INTEGRATED_BASIS}. Each
PC/NN graph stores one combined PDF with a sample diagnostic and all 15 raw
clustering resolutions. No candidate h5ad is saved and no candidate is selected
automatically. Valid completed graphs are detected and skipped on resume.
""",
        encoding="utf-8",
    )
    completion = pd.read_csv(COMPLETION).iloc[0].to_dict()
    completion["elapsed_seconds_this_invocation"] = time.time() - started
    print(json.dumps(completion, indent=2, default=str))


if __name__ == "__main__":
    main()
