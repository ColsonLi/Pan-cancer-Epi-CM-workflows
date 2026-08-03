#!/usr/bin/env python3
"""Finalize broad-clustering grid-search manifests without writing h5ad files."""

from __future__ import annotations

import platform
from pathlib import Path

import pandas as pd


WORKFLOW_ROOT = Path("/mnt/disk18t/lr_xcy/riku/crc_val/epi-cm-core-workflow")
STEP = "01-celltype_integration_clustering/05-clustering-parameter-search"
TABLE_DIR = WORKFLOW_ROOT / "tables" / STEP
FIGURE_DIR = WORKFLOW_ROOT / "figures" / STEP
H5AD_DIR = WORKFLOW_ROOT / "h5ad/01-celltype_integration_clustering/05-clustering-parameter-search"
CODE_DIR = WORKFLOW_ROOT / "codes" / STEP
INPUT_H5AD = (
    WORKFLOW_ROOT
    / "h5ad/01-celltype_integration_clustering/04-integration-harmony/adata_harmony.h5ad"
)


def read_many(pattern: str) -> pd.DataFrame:
    frames = [pd.read_csv(path) for path in sorted(TABLE_DIR.glob(pattern))]
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def main() -> None:
    review = read_many("worker_*_candidate_review_manifest.csv")
    grid = read_many("worker_*_grid_manifest.csv")
    status = read_many("worker_*_status.csv")

    if review.empty or grid.empty or status.empty:
        raise RuntimeError("Missing worker manifest/status files; cannot finalize grid outputs.")

    sort_cols = ["n_pcs", "n_neighbors", "resolution"]
    review = review.sort_values(sort_cols).reset_index(drop=True)
    grid = grid.sort_values(sort_cols).reset_index(drop=True)
    status = status.sort_values("worker_id").reset_index(drop=True)

    review.to_csv(TABLE_DIR / "candidate_review_manifest.csv", index=False)
    grid.to_csv(TABLE_DIR / "clustering_grid_manifest.csv", index=False)
    status.to_csv(TABLE_DIR / "worker_status_summary.csv", index=False)

    pdf_files = sorted(FIGURE_DIR.glob("pcs-*_nn-*_res-range/umap_*.pdf"))
    h5ad_files = sorted(H5AD_DIR.rglob("*.h5ad")) if H5AD_DIR.exists() else []

    missing_figures = [row.figure_pdf for row in review.itertuples() if not Path(row.figure_pdf).exists()]
    missing_cluster_tables = [
        row.cluster_count_table for row in review.itertuples() if not Path(row.cluster_count_table).exists()
    ]
    missing_parameter_tables = [
        row.parameter_table for row in review.itertuples() if not Path(row.parameter_table).exists()
    ]

    expected_total_grid = 9 * 9 * 12
    completed = int(grid["completed"].fillna(False).astype(bool).sum())
    skipped = int((~grid["completed"].fillna(False).astype(bool)).sum())
    completion = pd.DataFrame(
        [
            {
                "step": "clustering_parameter_search",
                "input_h5ad": str(INPUT_H5AD),
                "n_pcs_values": "10;15;20;25;30;35;40;45;50",
                "n_neighbors_values": "10;15;20;25;30;35;40;45;50",
                "resolution_values": "0.1-1.2 step 0.1",
                "stop_rule": "for each n_pcs/n_neighbors graph, stop increasing resolution after n_clusters > 20",
                "expected_grid_rows": expected_total_grid,
                "grid_rows_recorded": int(len(grid)),
                "completed_candidates": completed,
                "skipped_candidates": skipped,
                "review_manifest_rows": int(len(review)),
                "pdf_umap_files": int(len(pdf_files)),
                "candidate_h5ad_files": int(len(h5ad_files)),
                "missing_figures": int(len(missing_figures)),
                "missing_cluster_count_tables": int(len(missing_cluster_tables)),
                "missing_parameter_tables": int(len(missing_parameter_tables)),
                "all_workers_complete": bool((status["step"] == "worker_complete").all()),
                "candidate_h5ad_saved": bool(review["candidate_h5ad_saved"].fillna(False).astype(bool).any()),
                "backend": "rapids_singlecell",
                "random_seed": 42,
                "code_file": str(CODE_DIR / "02_rsc_clustering_grid_worker.py"),
                "finalize_code_file": str(CODE_DIR / "03_finalize_clustering_grid_outputs.py"),
            }
        ]
    )
    completion.to_csv(TABLE_DIR / "clustering_grid_completion_check.csv", index=False)

    pd.DataFrame({"missing_figure_pdf": missing_figures}).to_csv(
        TABLE_DIR / "missing_candidate_figures.csv", index=False
    )
    pd.DataFrame({"candidate_h5ad_file": [str(p) for p in h5ad_files]}).to_csv(
        TABLE_DIR / "candidate_h5ad_file_check.csv", index=False
    )

    package_versions = {
        "python": platform.python_version(),
        "pandas": pd.__version__,
        "finalize_code_file": str(CODE_DIR / "03_finalize_clustering_grid_outputs.py"),
        "worker_code_file": str(CODE_DIR / "02_rsc_clustering_grid_worker.py"),
    }
    with (TABLE_DIR / "package_versions.txt").open("w") as fh:
        for key, value in package_versions.items():
            fh.write(f"{key}: {value}\n")

    with (TABLE_DIR / "readme.txt").open("w") as fh:
        fh.write("05-clustering-parameter-search completed.\n")
        fh.write(f"Input h5ad: {INPUT_H5AD}\n")
        fh.write("Basis: obsm['X_pca_inte'] from Harmony output.\n")
        fh.write("Grid: n_pcs and n_neighbors 10-50 step 5; resolution 0.1-1.2 step 0.1.\n")
        fh.write("Stop rule: stop higher resolutions for a graph after n_clusters > 20.\n")
        fh.write("Candidate outputs: UMAP PDFs and per-candidate CSV tables only.\n")
        fh.write("No candidate h5ad files were written; selected h5ad is deferred until manual parameter choice.\n")
        fh.write("Main tables: candidate_review_manifest.csv, clustering_grid_manifest.csv, clustering_grid_completion_check.csv.\n")

    print(completion.to_string(index=False))


if __name__ == "__main__":
    main()
