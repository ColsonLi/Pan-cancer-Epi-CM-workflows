#!/usr/bin/env python3
"""Run one lineage subtype clustering grid with RAPIDS, saving figures/tables only."""

from __future__ import annotations

import argparse
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
MODULE = "02-cell_subtype_integration_clustering"
STEP = "03-subtype-clustering-grid"
CODE_FILE = Path(__file__)

PCS_VALUES = list(range(10, 51, 5))
NN_VALUES = list(range(10, 51, 5))
RES_VALUES = [round(x, 1) for x in np.arange(0.1, 1.51, 0.1)]
BASIS = "X_pca_inte"
UMAP_KEY = "X_umap"
NEIGHBORS_KEY = "neighbors"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lineage", required=True)
    parser.add_argument("--allow-existing", action="store_true")
    return parser.parse_args()


def graph_label(n_pcs: int, n_neighbors: int) -> str:
    return f"pcs-{n_pcs}_nn-{n_neighbors}_res-0p1-1p5"


def cluster_key(resolution: float) -> str:
    return f"leiden_res{str(resolution).replace('.', 'p')}"


def make_manifest(table_dir: Path, fig_dir: Path) -> pd.DataFrame:
    rows = []
    for n_pcs in PCS_VALUES:
        for n_neighbors in NN_VALUES:
            glabel = graph_label(n_pcs, n_neighbors)
            for resolution in RES_VALUES:
                rows.append(
                    {
                        "n_pcs": n_pcs,
                        "n_neighbors": n_neighbors,
                        "resolution": resolution,
                        "graph_label": glabel,
                        "candidate_label": f"pcs-{n_pcs}_nn-{n_neighbors}_res-{str(resolution).replace('.', 'p')}",
                        "expected_cluster_count_table": str(table_dir / glabel / "cluster_counts_by_resolution.csv"),
                        "expected_parameter_table": str(table_dir / glabel / "clustering_parameters_by_resolution.csv"),
                        "expected_figure_pdf": str(fig_dir / glabel / f"umap_{glabel}.pdf"),
                        "status": "planned",
                        "completed": False,
                        "reason_if_skipped": "",
                    }
                )
    return pd.DataFrame(rows)


def existing_graph_complete(table_dir: Path, fig_dir: Path, n_pcs: int, n_neighbors: int) -> bool:
    glabel = graph_label(n_pcs, n_neighbors)
    count_path = table_dir / glabel / "cluster_counts_by_resolution.csv"
    param_path = table_dir / glabel / "clustering_parameters_by_resolution.csv"
    fig_path = fig_dir / glabel / f"umap_{glabel}.pdf"
    if not (count_path.exists() and param_path.exists() and fig_path.exists()):
        return False
    params = pd.read_csv(param_path)
    return set(params["resolution"].round(1)) == set(RES_VALUES)


def main() -> None:
    args = parse_args()
    input_h5ad = WORKFLOW_ROOT / "h5ad" / MODULE / "02-subtype-harmony" / args.lineage / f"adata_{args.lineage}_harmony.h5ad"
    table_dir = WORKFLOW_ROOT / "tables" / MODULE / STEP / args.lineage
    fig_dir = WORKFLOW_ROOT / "figures" / MODULE / STEP / args.lineage
    table_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    sc.settings.autoshow = False
    sc.settings.set_figure_params(figsize=(3, 3), dpi=150)

    import cupy as cp
    import rapids_singlecell as rsc

    manifest_path = table_dir / "subtype_clustering_grid_manifest.csv"
    if manifest_path.exists():
        manifest = pd.read_csv(manifest_path)
    else:
        manifest = make_manifest(table_dir, fig_dir)
        manifest.to_csv(manifest_path, index=False)

    preflight = {
        "step": "subtype_rsc_clustering_preflight",
        "lineage": args.lineage,
        "input_h5ad": str(input_h5ad),
        "planned_backend": "rapids_singlecell",
        "basis": BASIS,
        "status": "not_run",
        "code_file": str(CODE_FILE),
        "random_seed": SEED,
    }
    try:
        tiny = ad.AnnData(X=np.ones((100, 5), dtype=np.float32))
        tiny.obsm[BASIS] = np.random.normal(size=(100, 10)).astype(np.float32)
        rsc.get.anndata_to_GPU(tiny)
        rsc.pp.neighbors(tiny, n_neighbors=10, n_pcs=10, use_rep=BASIS, random_state=SEED)
        rsc.tl.umap(tiny, random_state=SEED)
        rsc.tl.leiden(tiny, resolution=0.5, random_state=SEED)
        rsc.get.anndata_to_CPU(tiny)
        preflight["status"] = "ok"
    except Exception as exc:
        preflight["status"] = "failed"
        preflight["error_summary"] = "".join(traceback.format_exception_only(type(exc), exc)).strip()
        pd.DataFrame([preflight]).to_csv(table_dir / "gpu_backend_capability_summary.csv", index=False)
        raise
    pd.DataFrame([preflight]).to_csv(table_dir / "gpu_backend_capability_summary.csv", index=False)

    source = ad.read_h5ad(input_h5ad, backed="r")
    if BASIS not in source.obsm:
        raise KeyError(f"Missing obsm[{BASIS!r}] in {input_h5ad}")
    obs = source.obs.copy()
    basis = np.asarray(source.obsm[BASIS]).astype(np.float32)
    source.file.close()

    search = ad.AnnData(X=np.zeros((obs.shape[0], 1), dtype=np.float32), obs=obs)
    search.obsm[BASIS] = basis
    del source, obs, basis
    gc.collect()

    review_rows: list[dict] = []
    review_path = table_dir / "candidate_review_manifest.csv"
    if review_path.exists():
        review_rows = pd.read_csv(review_path).to_dict("records")

    rsc.get.anndata_to_GPU(search)
    try:
        for n_pcs in PCS_VALUES:
            for n_neighbors in NN_VALUES:
                glabel = graph_label(n_pcs, n_neighbors)
                graph_table_dir = table_dir / glabel
                graph_fig_dir = fig_dir / glabel
                graph_table_dir.mkdir(parents=True, exist_ok=True)
                graph_fig_dir.mkdir(parents=True, exist_ok=True)

                if args.allow_existing and existing_graph_complete(table_dir, fig_dir, n_pcs, n_neighbors):
                    manifest.loc[manifest["graph_label"] == glabel, ["status", "completed"]] = ["completed_existing", True]
                    manifest.to_csv(manifest_path, index=False)
                    continue

                adata_run = search.copy()
                count_rows: list[dict] = []
                param_rows: list[dict] = []
                try:
                    rsc.pp.neighbors(
                        adata_run,
                        n_neighbors=n_neighbors,
                        n_pcs=n_pcs,
                        use_rep=BASIS,
                        random_state=SEED,
                        key_added=NEIGHBORS_KEY,
                    )
                    rsc.tl.umap(
                        adata_run,
                        random_state=SEED,
                        neighbors_key=NEIGHBORS_KEY,
                        key_added=UMAP_KEY,
                    )
                    keys = []
                    for resolution in RES_VALUES:
                        key = cluster_key(resolution)
                        rsc.tl.leiden(
                            adata_run,
                            resolution=resolution,
                            key_added=key,
                            random_state=SEED,
                            neighbors_key=NEIGHBORS_KEY,
                        )
                        keys.append(key)

                    rsc.get.anndata_to_CPU(adata_run)
                    for resolution, key in zip(RES_VALUES, keys):
                        counts = adata_run.obs[key].astype(str).value_counts().rename_axis("cluster").reset_index(name="n_cells")
                        n_clusters = int(counts.shape[0])
                        for _, row in counts.iterrows():
                            count_rows.append(
                                {
                                    "lineage": args.lineage,
                                    "graph_label": glabel,
                                    "n_pcs": n_pcs,
                                    "n_neighbors": n_neighbors,
                                    "resolution": resolution,
                                    "cluster_key": key,
                                    "cluster": row["cluster"],
                                    "n_cells": int(row["n_cells"]),
                                }
                            )
                        param_rows.append(
                            {
                                "lineage": args.lineage,
                                "source_h5ad": str(input_h5ad),
                                "graph_label": glabel,
                                "candidate_label": f"pcs-{n_pcs}_nn-{n_neighbors}_res-{str(resolution).replace('.', 'p')}",
                                "n_pcs": n_pcs,
                                "n_neighbors": n_neighbors,
                                "resolution": resolution,
                                "basis": BASIS,
                                "neighbors_key": NEIGHBORS_KEY,
                                "umap_key": UMAP_KEY,
                                "cluster_key": key,
                                "n_clusters": n_clusters,
                                "backend_package": "rapids_singlecell",
                                "candidate_h5ad_saved": False,
                                "candidate_object_created_by_copy": True,
                                "candidate_object_deleted_after_outputs": True,
                                "code_file": str(CODE_FILE),
                                "random_seed": SEED,
                            }
                        )

                    count_path = graph_table_dir / "cluster_counts_by_resolution.csv"
                    param_path = graph_table_dir / "clustering_parameters_by_resolution.csv"
                    pd.DataFrame(count_rows).to_csv(count_path, index=False)
                    pd.DataFrame(param_rows).to_csv(param_path, index=False)

                    sc.settings.figdir = str(graph_fig_dir)
                    sc.pl.umap(
                        adata_run,
                        color=keys,
                        ncols=5,
                        wspace=0.35,
                        save=f"_{glabel}.pdf",
                        show=False,
                    )
                    fig_path = graph_fig_dir / f"umap_{glabel}.pdf"

                    for row in param_rows:
                        row = dict(row)
                        row["figure_pdf"] = str(fig_path)
                        row["cluster_count_table"] = str(count_path)
                        row["parameter_table"] = str(param_path)
                        review_rows.append(row)
                    pd.DataFrame(review_rows).to_csv(review_path, index=False)
                    manifest.loc[manifest["graph_label"] == glabel, ["status", "completed"]] = ["completed", True]
                    manifest.to_csv(manifest_path, index=False)
                except Exception as exc:
                    manifest.loc[manifest["graph_label"] == glabel, ["status", "reason_if_skipped"]] = [
                        "failed",
                        "".join(traceback.format_exception_only(type(exc), exc)).strip(),
                    ]
                    manifest.to_csv(manifest_path, index=False)
                    raise
                finally:
                    try:
                        rsc.get.anndata_to_CPU(adata_run)
                    except Exception:
                        pass
                    del adata_run
                    cp.get_default_memory_pool().free_all_blocks()
                    gc.collect()
    finally:
        try:
            rsc.get.anndata_to_CPU(search)
        except Exception:
            pass
        cp.get_default_memory_pool().free_all_blocks()

    manifest = pd.read_csv(manifest_path)
    expected = len(manifest)
    completed = int(manifest["completed"].astype(bool).sum())
    failed = int((manifest["status"] == "failed").sum())
    missing = int((manifest["status"] == "planned").sum())
    graph_expected = len(PCS_VALUES) * len(NN_VALUES)
    completed_graphs = int(manifest.loc[manifest["completed"].astype(bool), "graph_label"].nunique())
    all_outputs_exist = True
    for glabel in sorted(manifest.loc[manifest["completed"].astype(bool), "graph_label"].unique()):
        if not (
            (table_dir / glabel / "cluster_counts_by_resolution.csv").exists()
            and (table_dir / glabel / "clustering_parameters_by_resolution.csv").exists()
            and (fig_dir / glabel / f"umap_{glabel}.pdf").exists()
        ):
            all_outputs_exist = False
            break
    completion = pd.DataFrame(
        [
            {
                "lineage": args.lineage,
                "expected_resolution_candidates": expected,
                "completed_resolution_candidates": completed,
                "expected_graph_candidates": graph_expected,
                "completed_graph_candidates": completed_graphs,
                "failed_candidates": failed,
                "missing_candidates": missing,
                "all_expected_lightweight_outputs_exist": bool(all_outputs_exist),
                "manual_selection_ready": bool(completed == expected and failed == 0 and missing == 0 and all_outputs_exist),
                "candidate_h5ad_saved": False,
            }
        ]
    )
    completion.to_csv(table_dir / "clustering_grid_completion_check.csv", index=False)

    with (table_dir / "package_versions.txt").open("w") as fh:
        fh.write(f"python: {platform.python_version()}\n")
        fh.write(f"anndata: {ad.__version__}\n")
        fh.write(f"scanpy: {sc.__version__}\n")
        fh.write(f"numpy: {np.__version__}\n")
        fh.write(f"pandas: {pd.__version__}\n")
        try:
            import rapids_singlecell as rsc_version
            fh.write(f"rapids_singlecell: {rsc_version.__version__}\n")
        except Exception:
            fh.write("rapids_singlecell: unavailable_version\n")
        try:
            import cupy
            fh.write(f"cupy: {cupy.__version__}\n")
        except Exception:
            fh.write("cupy: unavailable_version\n")
        fh.write(f"code_file: {CODE_FILE}\n")

    with (table_dir / "readme.txt").open("w") as fh:
        fh.write(f"{args.lineage} subtype clustering grid completed with RAPIDS.\n")
        fh.write(f"Input Harmony h5ad: {input_h5ad}\n")
        fh.write("Grid: pcs 10-50 step 5; n_neighbors 10-50 step 5; Leiden resolution 0.1-1.5 step 0.1.\n")
        fh.write("Candidate h5ad files were not saved. Outputs are lightweight figures and tables only.\n")

    print(completion.to_string(index=False))


if __name__ == "__main__":
    main()
