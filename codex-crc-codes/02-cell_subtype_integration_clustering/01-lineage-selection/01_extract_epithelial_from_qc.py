#!/usr/bin/env python3
"""Extract score/rank-consistent epithelial cells from adata_qc."""

from __future__ import annotations

import platform
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

WORKFLOW_ROOT = Path("/mnt/disk18t/lr_xcy/riku/crc_val/epi-cm-core-workflow")
MODULE = "02-cell_subtype_integration_clustering"
STEP = "01-lineage-selection"
LINEAGE = "epithelial"
TARGET_LABEL = "Epithelial Cells"

QC_H5AD = WORKFLOW_ROOT / "h5ad/01-celltype_integration_clustering/03-qc/adata_qc.h5ad"
ANNOTATION_H5AD = (
    WORKFLOW_ROOT
    / "h5ad/01-celltype_integration_clustering/07-score-rank-qc/adata_anno_score_genes_rank_consistent.h5ad"
)
OUTPUT_H5AD = WORKFLOW_ROOT / "h5ad" / MODULE / STEP / LINEAGE / "adata_epithelial_qc.h5ad"
TABLE_DIR = WORKFLOW_ROOT / "tables" / MODULE / STEP / LINEAGE
CODE_FILE = Path(__file__)

TRANSFER_COLUMNS = [
    "leiden_res0p3",
    "leiden_coarse",
    "cell_type",
    "best_rank_type_global",
    "best_rank_score_rank_pct",
    "best_rank_score_col",
    "score_rank_consistent",
]


def assert_no_overwrite() -> None:
    outputs = [
        OUTPUT_H5AD,
        TABLE_DIR / "lineage_selection_parameters.csv",
        TABLE_DIR / "lineage_selection_summary.csv",
        TABLE_DIR / "selected_cell_ids.csv",
        TABLE_DIR / "readme.txt",
        TABLE_DIR / "package_versions.txt",
    ]
    existing = [str(path) for path in outputs if path.exists()]
    if existing:
        raise FileExistsError("Lineage selection output already exists; refusing to overwrite:\n" + "\n".join(existing))


def main() -> None:
    assert_no_overwrite()
    OUTPUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    anno = ad.read_h5ad(ANNOTATION_H5AD, backed="r")
    if "leiden_coarse" not in anno.obs.columns:
        raise KeyError("Missing obs['leiden_coarse'] in annotation h5ad")
    if "score_rank_consistent" not in anno.obs.columns:
        raise KeyError("Missing obs['score_rank_consistent'] in annotation h5ad")

    mask = (anno.obs["leiden_coarse"].astype(str) == TARGET_LABEL) & anno.obs["score_rank_consistent"].astype(bool)
    selected_ids = anno.obs_names[mask].copy()
    transfer_cols = [col for col in TRANSFER_COLUMNS if col in anno.obs.columns]
    transfer = anno.obs.loc[selected_ids, transfer_cols].copy()
    anno.file.close()

    qc = ad.read_h5ad(QC_H5AD)
    missing_ids = selected_ids.difference(qc.obs_names)
    if len(missing_ids) > 0:
        raise RuntimeError(f"{len(missing_ids)} selected IDs are absent from adata_qc.")
    epi = qc[selected_ids].copy()
    del qc

    for col in transfer_cols:
        epi.obs[col] = transfer.loc[epi.obs_names, col].values
    epi.obs["lineage_selection_source"] = str(ANNOTATION_H5AD)
    epi.obs["lineage_selection_label"] = TARGET_LABEL
    epi.obs["lineage_selection_rank_filter"] = "score_rank_consistent"

    epi.write_h5ad(OUTPUT_H5AD)

    pd.DataFrame({"cell_id": selected_ids}).to_csv(TABLE_DIR / "selected_cell_ids.csv", index=False)
    summary = pd.DataFrame(
        [
            {
                "lineage": LINEAGE,
                "target_leiden_coarse": TARGET_LABEL,
                "annotation_source_h5ad": str(ANNOTATION_H5AD),
                "expression_source_h5ad": str(QC_H5AD),
                "output_h5ad": str(OUTPUT_H5AD),
                "n_selected_cells": int(epi.n_obs),
                "n_genes": int(epi.n_vars),
                "rank_filter_mode": "required_consistent_object",
                "all_selected_score_rank_consistent": bool(epi.obs["score_rank_consistent"].astype(bool).all()),
                "obs_names_unique": bool(epi.obs_names.is_unique),
            }
        ]
    )
    summary.to_csv(TABLE_DIR / "lineage_selection_summary.csv", index=False)

    params = summary.copy()
    params["transfer_columns"] = ";".join(transfer_cols)
    params["code_file"] = str(CODE_FILE)
    params.to_csv(TABLE_DIR / "lineage_selection_parameters.csv", index=False)

    with (TABLE_DIR / "package_versions.txt").open("w") as fh:
        fh.write(f"python: {platform.python_version()}\n")
        fh.write(f"anndata: {ad.__version__}\n")
        fh.write(f"numpy: {np.__version__}\n")
        fh.write(f"pandas: {pd.__version__}\n")
        fh.write(f"code_file: {CODE_FILE}\n")

    with (TABLE_DIR / "readme.txt").open("w") as fh:
        fh.write("Epithelial lineage selection completed.\n")
        fh.write(f"Annotation/ID source: {ANNOTATION_H5AD}\n")
        fh.write(f"Expression/reclustering source: {QC_H5AD}\n")
        fh.write(f"Target label: {TARGET_LABEL}\n")
        fh.write("Rank filter: selected IDs must come from adata_anno_score_genes_rank_consistent.h5ad.\n")
        fh.write(f"Output h5ad: {OUTPUT_H5AD}\n")

    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
