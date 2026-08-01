#!/usr/bin/env python3
"""Disjoint RAPIDS clustering-grid worker.

Each worker owns a non-overlapping n_pcs range to avoid CSV write races.
Candidate h5ad files are never saved.
"""

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
INPUT_H5AD = (
    WORKFLOW_ROOT
    / "h5ad/01-celltype_integration_clustering/04-integration-harmony/adata_harmony.h5ad"
)
TABLE_DIR = WORKFLOW_ROOT / "tables/01-celltype_integration_clustering/05-clustering-parameter-search"
FIG_DIR = WORKFLOW_ROOT / "figures/01-celltype_integration_clustering/05-clustering-parameter-search"
CODE_FILE = Path(__file__)

ALL_PCS_VALUES = list(range(10, 51, 5))
NN_VALUES = list(range(10, 51, 5))
RES_VALUES = [round(x, 1) for x in np.arange(0.1, 1.21, 0.1)]
MAX_CLUSTERS_BEFORE_STOP_HIGHER_RES = 20
BASIS = "X_pca_inte"
UMAP_KEY = "X_umap"


def label_for(n_pcs: int, n_neighbors: int, resolution: float) -> str:
    return f"pcs-{n_pcs}_nn-{n_neighbors}_res-{str(resolution).replace('.', 'p')}"


def graph_label(n_pcs: int, n_neighbors: int) -> str:
    return f"pcs-{n_pcs}_nn-{n_neighbors}_res-range"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pcs-min", type=int, required=True)
    parser.add_argument("--pcs-max", type=int, required=True)
    parser.add_argument("--worker-id", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    pcs_values = [p for p in ALL_PCS_VALUES if args.pcs_min <= p <= args.pcs_max]
    if not pcs_values:
        raise ValueError(f"No pcs values in range {args.pcs_min}-{args.pcs_max}")

    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    sc.settings.autoshow = False
    sc.settings.set_figure_params(figsize=(3, 3), dpi=150)

    import cupy as cp
    import rapids_singlecell as rsc

    worker_rows = []
    review_rows = []
    worker_manifest = TABLE_DIR / f"worker_{args.worker_id}_grid_manifest.csv"
    worker_review = TABLE_DIR / f"worker_{args.worker_id}_candidate_review_manifest.csv"
    worker_log = TABLE_DIR / f"worker_{args.worker_id}_status.csv"

    preflight = {
        "worker_id": args.worker_id,
        "step": "worker_start",
        "pcs_values": ";".join(map(str, pcs_values)),
        "input_h5ad": str(INPUT_H5AD),
        "backend": "rapids_singlecell",
        "status": "running",
        "code_file": str(CODE_FILE),
        "random_seed": SEED,
    }
    pd.DataFrame([preflight]).to_csv(worker_log, index=False)

    source = ad.read_h5ad(INPUT_H5AD, backed="r")
    obs = source.obs.copy()
    basis = np.asarray(source.obsm[BASIS]).astype(np.float32)
    source.file.close()
    search = ad.AnnData(
        X=np.zeros((obs.shape[0], 1), dtype=np.float32),
        obs=obs,
    )
    search.obsm[BASIS] = basis
    del source, obs, basis
    gc.collect()

    rsc.get.anndata_to_GPU(search)
    for n_pcs in pcs_values:
        for n_neighbors in NN_VALUES:
            glabel = graph_label(n_pcs, n_neighbors)
            graph_table_dir = TABLE_DIR / glabel
            graph_fig_dir = FIG_DIR / glabel
            graph_table_dir.mkdir(parents=True, exist_ok=True)
            graph_fig_dir.mkdir(parents=True, exist_ok=True)

            adata_run = search.copy()
            try:
                rsc.pp.neighbors(
                    adata_run,
                    n_neighbors=n_neighbors,
                    n_pcs=n_pcs,
                    use_rep=BASIS,
                    random_state=SEED,
                    key_added="neighbors",
                )
                rsc.tl.umap(
                    adata_run,
                    random_state=SEED,
                    neighbors_key="neighbors",
                    key_added=UMAP_KEY,
                )
                stop_higher_res = False
                for resolution in RES_VALUES:
                    cand = label_for(n_pcs, n_neighbors, resolution)
                    param_path = graph_table_dir / f"{cand}_parameters.csv"
                    count_path = graph_table_dir / f"{cand}_cluster_counts.csv"
                    fig_path = graph_fig_dir / f"umap_{cand}.pdf"

                    if stop_higher_res:
                        worker_rows.append(
                            {
                                "worker_id": args.worker_id,
                                "candidate_label": cand,
                                "n_pcs": n_pcs,
                                "n_neighbors": n_neighbors,
                                "resolution": resolution,
                                "status": "skipped_user_approved",
                                "completed": False,
                                "reason_if_skipped": f"previous resolution exceeded {MAX_CLUSTERS_BEFORE_STOP_HIGHER_RES} clusters for same graph",
                            }
                        )
                        continue

                    if param_path.exists() and count_path.exists() and fig_path.exists():
                        old = pd.read_csv(param_path).iloc[0].to_dict()
                        n_clusters = int(old["n_clusters"])
                        status = "completed_existing"
                    else:
                        cluster_key = f"leiden_res{str(resolution).replace('.', 'p')}"
                        try:
                            rsc.tl.leiden(
                                adata_run,
                                resolution=resolution,
                                key_added=cluster_key,
                                random_state=SEED,
                                neighbors_key="neighbors",
                            )
                            rsc.get.anndata_to_CPU(adata_run)
                            counts = (
                                adata_run.obs[cluster_key]
                                .astype(str)
                                .value_counts()
                                .rename_axis("cluster")
                                .reset_index(name="n_cells")
                            )
                            n_clusters = int(counts.shape[0])
                            counts.to_csv(count_path, index=False)
                            params = pd.DataFrame(
                                [
                                    {
                                        "source_h5ad": str(INPUT_H5AD),
                                        "candidate_label": cand,
                                        "n_pcs": n_pcs,
                                        "n_neighbors": n_neighbors,
                                        "resolution": resolution,
                                        "basis": BASIS,
                                        "neighbors_key": "neighbors",
                                        "umap_key": UMAP_KEY,
                                        "cluster_key": cluster_key,
                                        "n_clusters": n_clusters,
                                        "backend_package": "rapids_singlecell",
                                        "candidate_h5ad_saved": False,
                                        "candidate_object_created_by_copy": True,
                                        "candidate_object_deleted_after_outputs": False,
                                        "stop_rule": f"stop higher resolutions after n_clusters > {MAX_CLUSTERS_BEFORE_STOP_HIGHER_RES}",
                                        "code_file": str(CODE_FILE),
                                        "random_seed": SEED,
                                        "worker_id": args.worker_id,
                                    }
                                ]
                            )
                            params.to_csv(param_path, index=False)
                            sc.settings.figdir = str(graph_fig_dir)
                            sc.pl.umap(adata_run, color=cluster_key, save=f"_{cand}.pdf", show=False)
                            status = "completed"
                            rsc.get.anndata_to_GPU(adata_run)
                        except Exception as exc:
                            worker_rows.append(
                                {
                                    "worker_id": args.worker_id,
                                    "candidate_label": cand,
                                    "n_pcs": n_pcs,
                                    "n_neighbors": n_neighbors,
                                    "resolution": resolution,
                                    "status": "failed",
                                    "completed": False,
                                    "reason_if_skipped": "",
                                    "error_summary": "".join(traceback.format_exception_only(type(exc), exc)).strip(),
                                }
                            )
                            pd.DataFrame(worker_rows).to_csv(worker_manifest, index=False)
                            raise

                    worker_rows.append(
                        {
                            "worker_id": args.worker_id,
                            "candidate_label": cand,
                            "n_pcs": n_pcs,
                            "n_neighbors": n_neighbors,
                            "resolution": resolution,
                            "n_clusters": n_clusters,
                            "status": status,
                            "completed": True,
                            "reason_if_skipped": "",
                        }
                    )
                    review_rows.append(
                        {
                            "worker_id": args.worker_id,
                            "candidate_label": cand,
                            "n_pcs": n_pcs,
                            "n_neighbors": n_neighbors,
                            "resolution": resolution,
                            "n_clusters": n_clusters,
                            "figure_pdf": str(fig_path),
                            "cluster_count_table": str(count_path),
                            "parameter_table": str(param_path),
                            "candidate_h5ad_saved": False,
                        }
                    )
                    pd.DataFrame(worker_rows).to_csv(worker_manifest, index=False)
                    pd.DataFrame(review_rows).to_csv(worker_review, index=False)
                    if n_clusters > MAX_CLUSTERS_BEFORE_STOP_HIGHER_RES:
                        stop_higher_res = True
            finally:
                try:
                    rsc.get.anndata_to_CPU(adata_run)
                except Exception:
                    pass
                del adata_run
                cp.get_default_memory_pool().free_all_blocks()
                gc.collect()

    try:
        rsc.get.anndata_to_CPU(search)
    except Exception:
        pass
    cp.get_default_memory_pool().free_all_blocks()

    pd.DataFrame(worker_rows).to_csv(worker_manifest, index=False)
    pd.DataFrame(review_rows).to_csv(worker_review, index=False)
    status = {
        "worker_id": args.worker_id,
        "step": "worker_complete",
        "pcs_values": ";".join(map(str, pcs_values)),
        "completed_candidates": int(sum(r.get("completed", False) for r in worker_rows)),
        "skipped_candidates": int(sum(str(r.get("status", "")).startswith("skipped") for r in worker_rows)),
        "candidate_h5ad_saved": False,
        "backend": "rapids_singlecell",
        "python": platform.python_version(),
        "rapids_singlecell": rsc.__version__,
        "cupy": cp.__version__,
        "code_file": str(CODE_FILE),
        "random_seed": SEED,
    }
    pd.DataFrame([status]).to_csv(worker_log, index=False)
    print(pd.DataFrame([status]).to_string(index=False))


if __name__ == "__main__":
    main()
