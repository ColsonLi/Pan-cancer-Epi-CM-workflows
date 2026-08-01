#!/usr/bin/env python3
"""RAPIDS clustering grid search, UMAP/table outputs only."""

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
INPUT_H5AD = (
    WORKFLOW_ROOT
    / "h5ad/01-celltype_integration_clustering/04-integration-harmony/adata_harmony.h5ad"
)
TABLE_DIR = WORKFLOW_ROOT / "tables/01-celltype_integration_clustering/05-clustering-parameter-search"
FIG_DIR = WORKFLOW_ROOT / "figures/01-celltype_integration_clustering/05-clustering-parameter-search"
CODE_FILE = Path(__file__)

PCS_VALUES = list(range(10, 51, 5))
NN_VALUES = list(range(10, 51, 5))
RES_VALUES = [round(x, 1) for x in np.arange(0.1, 1.21, 0.1)]
MAX_CLUSTERS_BEFORE_STOP_HIGHER_RES = 20
BASIS = "X_pca_inte"
LEIDEN_PREFIX = "leiden"
UMAP_KEY = "X_umap"


def label_for(n_pcs: int, n_neighbors: int, resolution: float) -> str:
    return f"pcs-{n_pcs}_nn-{n_neighbors}_res-{str(resolution).replace('.', 'p')}"


def graph_label(n_pcs: int, n_neighbors: int) -> str:
    return f"pcs-{n_pcs}_nn-{n_neighbors}_res-range"


def write_manifest(rows: list[dict]) -> None:
    pd.DataFrame(rows).to_csv(TABLE_DIR / "clustering_grid_manifest.csv", index=False)


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    sc.settings.autoshow = False
    sc.settings.figdir = str(FIG_DIR)
    sc.settings.set_figure_params(figsize=(3, 3), dpi=150)

    import cupy as cp
    import rapids_singlecell as rsc

    expected_rows = []
    for n_pcs in PCS_VALUES:
        for n_neighbors in NN_VALUES:
            for resolution in RES_VALUES:
                cand = label_for(n_pcs, n_neighbors, resolution)
                expected_rows.append(
                    {
                        "n_pcs": n_pcs,
                        "n_neighbors": n_neighbors,
                        "resolution": resolution,
                        "candidate_label": cand,
                        "expected_cluster_count_table": str(TABLE_DIR / graph_label(n_pcs, n_neighbors) / f"{cand}_cluster_counts.csv"),
                        "expected_parameter_table": str(TABLE_DIR / graph_label(n_pcs, n_neighbors) / f"{cand}_parameters.csv"),
                        "expected_figure_dir": str(FIG_DIR / graph_label(n_pcs, n_neighbors)),
                        "status": "planned",
                        "completed": False,
                        "reason_if_skipped": "",
                    }
                )
    write_manifest(expected_rows)

    preflight = {
        "step": "rsc_clustering_preflight",
        "input_h5ad": str(INPUT_H5AD),
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
        pd.DataFrame([preflight]).to_csv(TABLE_DIR / "gpu_backend_capability_summary.csv", index=False)
        raise
    pd.DataFrame([preflight]).to_csv(TABLE_DIR / "gpu_backend_capability_summary.csv", index=False)

    # Keep only metadata and the integrated PCA basis for the parameter search.
    source = ad.read_h5ad(INPUT_H5AD)
    search = ad.AnnData(
        X=np.zeros((source.n_obs, 1), dtype=np.float32),
        obs=source.obs.copy(),
    )
    search.obsm[BASIS] = source.obsm[BASIS].astype(np.float32)
    del source
    gc.collect()

    manifest = pd.read_csv(TABLE_DIR / "clustering_grid_manifest.csv")
    completed_candidate_rows = []
    review_rows = []

    rsc.get.anndata_to_GPU(search)
    for n_pcs in PCS_VALUES:
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
                    if stop_higher_res:
                        manifest.loc[manifest["candidate_label"] == cand, ["status", "reason_if_skipped"]] = [
                            "skipped_user_approved",
                            f"previous resolution exceeded {MAX_CLUSTERS_BEFORE_STOP_HIGHER_RES} clusters for same n_pcs/n_neighbors",
                        ]
                        continue

                    cluster_key = f"{LEIDEN_PREFIX}_res{str(resolution).replace('.', 'p')}"
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
                    counts.to_csv(graph_table_dir / f"{cand}_cluster_counts.csv", index=False)
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
                            }
                        ]
                    )
                    params.to_csv(graph_table_dir / f"{cand}_parameters.csv", index=False)

                    sc.settings.figdir = str(graph_fig_dir)
                    sc.pl.umap(
                        adata_run,
                        color=cluster_key,
                        save=f"_{cand}.pdf",
                        show=False,
                    )

                    manifest.loc[manifest["candidate_label"] == cand, ["status", "completed"]] = [
                        "completed",
                        True,
                    ]
                    completed_candidate_rows.append(params.iloc[0].to_dict())
                    review_rows.append(
                        {
                            "candidate_label": cand,
                            "n_pcs": n_pcs,
                            "n_neighbors": n_neighbors,
                            "resolution": resolution,
                            "n_clusters": n_clusters,
                            "figure_pdf": str(graph_fig_dir / f"umap_{cand}.pdf"),
                            "cluster_count_table": str(graph_table_dir / f"{cand}_cluster_counts.csv"),
                            "parameter_table": str(graph_table_dir / f"{cand}_parameters.csv"),
                            "candidate_h5ad_saved": False,
                        }
                    )
                    pd.DataFrame(review_rows).to_csv(TABLE_DIR / "candidate_review_manifest.csv", index=False)
                    write_manifest(manifest.to_dict("records"))
                    if n_clusters > MAX_CLUSTERS_BEFORE_STOP_HIGHER_RES:
                        stop_higher_res = True
                    rsc.get.anndata_to_GPU(adata_run)
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

    manifest = pd.read_csv(TABLE_DIR / "clustering_grid_manifest.csv")
    expected = len(manifest)
    completed = int((manifest["status"] == "completed").sum())
    skipped = int(manifest["status"].astype(str).str.startswith("skipped").sum())
    failed = int((manifest["status"] == "failed").sum())
    missing = int((manifest["status"] == "planned").sum())
    all_outputs_exist = True
    for _, row in manifest[manifest["status"] == "completed"].iterrows():
        if not Path(row["expected_cluster_count_table"]).exists() or not Path(row["expected_parameter_table"]).exists():
            all_outputs_exist = False
            break
    completion = pd.DataFrame(
        [
            {
                "expected_candidates_full_grid": expected,
                "completed_candidates": completed,
                "failed_candidates": failed,
                "skipped_candidates_by_cluster_stop_rule": skipped,
                "missing_candidates": missing,
                "all_expected_lightweight_outputs_exist": all_outputs_exist,
                "manual_selection_ready": failed == 0 and missing == 0 and all_outputs_exist,
                "candidate_h5ad_saved": False,
                "user_grid": "pcs 10-50 step5; nn 10-50 step5; res 0.1-1.2; stop higher res when clusters >20",
            }
        ]
    )
    completion.to_csv(TABLE_DIR / "clustering_grid_completion_check.csv", index=False)

    with (TABLE_DIR / "readme.txt").open("w", encoding="utf-8") as fh:
        fh.write("05-clustering-parameter-search completed.\n")
        fh.write(f"Input: {INPUT_H5AD}\n")
        fh.write("Candidate h5ad files were not saved, by user request.\n")
        fh.write("Final selected h5ad will be saved only after user chooses final parameters.\n")

    versions = [
        f"python={platform.python_version()}",
        f"platform={platform.platform()}",
        f"anndata={ad.__version__}",
        f"scanpy={sc.__version__}",
        f"numpy={np.__version__}",
        f"pandas={pd.__version__}",
        f"rapids_singlecell={rsc.__version__}",
        f"cupy={cp.__version__}",
        "environment=/mnt/disk18t/lr_xcy/riku/crc_val/uv_envs/rapids/.venv",
        f"code_file={CODE_FILE}",
        f"random_seed={SEED}",
    ]
    (TABLE_DIR / "package_versions.txt").write_text("\n".join(versions) + "\n", encoding="utf-8")
    print(completion.to_string(index=False))


if __name__ == "__main__":
    main()
