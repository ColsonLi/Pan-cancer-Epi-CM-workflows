#!/usr/bin/env python3
"""Select strict-consistent BRCA epithelial IDs for joint subtype analysis."""

from __future__ import annotations

import importlib.metadata
import json
import random
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


SEED = 42
random.seed(SEED)
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/brca_cm-epi_codex")
WORKFLOW = ROOT / "epi-cm-core-workflow"
BLOCK01 = "01-celltype_integration_clustering"
BLOCK02 = "02-cell_subtype_integration_clustering"
ANNOTATION_H5AD = (
    WORKFLOW
    / "h5ad"
    / BLOCK01
    / "07-score-rank-qc"
    / "adata_anno_score_genes_rank_consistent.h5ad"
)
QC_H5AD = WORKFLOW / "h5ad" / BLOCK01 / "03-qc" / "adata_qc.h5ad"
TABLE_DIR = WORKFLOW / "tables" / BLOCK02 / "01-lineage-selection" / "epithelial"
CODE_PATH = (
    WORKFLOW
    / "codes"
    / BLOCK02
    / "01-lineage-selection"
    / "epithelial"
    / "01_select_epithelial_ids.py"
)

TARGET_LABEL = "Epithelial Cells"
COARSE_KEY = "leiden_coarse"
BEST_KEY = "best_rank_type_global"
CONSISTENT_KEY = "score_rank_consistent"


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    annotation = ad.read_h5ad(ANNOTATION_H5AD, backed="r")
    required_annotation = {
        COARSE_KEY,
        "cell_type",
        BEST_KEY,
        CONSISTENT_KEY,
        "sample",
        "series",
        "status",
        "original_barcode",
        "leiden_res0p8",
    }
    missing = required_annotation - set(annotation.obs.columns)
    if missing:
        raise ValueError(f"Strict annotation source lacks columns: {sorted(missing)}")
    if not annotation.obs_names.is_unique:
        raise ValueError("Strict annotation source has non-unique cell IDs.")
    if not annotation.obs[CONSISTENT_KEY].astype(bool).all():
        raise ValueError("Strict annotation source contains score/rank-inconsistent cells.")
    if not annotation.obs[COARSE_KEY].astype(str).equals(
        annotation.obs[BEST_KEY].astype(str)
    ):
        raise ValueError("Strict annotation source violates broad best-rank equality.")

    all_samples = sorted(annotation.obs["sample"].astype(str).unique())
    target_mask = annotation.obs[COARSE_KEY].astype(str).eq(TARGET_LABEL)
    selected = annotation.obs.loc[target_mask].copy()
    if selected.empty:
        raise ValueError(f"No strict-consistent cells found for {TARGET_LABEL}.")
    selected_ids = pd.Index(selected.index.astype(str))
    selected_samples = sorted(selected["sample"].astype(str).unique())
    absent_samples = sorted(set(all_samples) - set(selected_samples))

    qc = ad.read_h5ad(QC_H5AD, backed="r")
    if not qc.obs_names.is_unique:
        raise ValueError("QC expression base has non-unique cell IDs.")
    missing_ids = selected_ids.difference(qc.obs_names.astype(str))
    if len(missing_ids):
        raise ValueError(f"Selected epithelial IDs missing from adata_qc: {len(missing_ids)}")

    export_columns = [
        "sample",
        "series",
        "status",
        "original_barcode",
        "leiden_res0p8",
        COARSE_KEY,
        "cell_type",
        BEST_KEY,
        CONSISTENT_KEY,
    ]
    score_columns = sorted(
        column
        for column in selected.columns
        if column.endswith("_score") or column.endswith("_score_rank_pct")
    )
    selected_table = selected[export_columns + score_columns].copy()
    selected_table.insert(0, "cell_id", selected_ids)
    selected_table.to_csv(TABLE_DIR / "selected_cell_ids.csv", index=False)

    sample_counts = (
        selected.assign(sample=selected["sample"].astype(str))
        .groupby("sample", observed=True)
        .size()
        .rename("n_selected_epithelial_cells")
        .reindex(all_samples, fill_value=0)
        .rename_axis("sample")
        .reset_index()
    )
    sample_counts["has_selected_epithelial_cells"] = (
        sample_counts["n_selected_epithelial_cells"] > 0
    )
    sample_counts.to_csv(TABLE_DIR / "selected_epithelial_cells_by_sample.csv", index=False)

    pd.DataFrame(
        [
            {
                "annotation_source_h5ad": str(ANNOTATION_H5AD),
                "expression_base_h5ad": str(QC_H5AD),
                "target_leiden_coarse": TARGET_LABEL,
                "rank_filter_mode": "strict_consistent",
                "selection_rule": (
                    "IDs from strict-consistent object where "
                    "leiden_coarse == Epithelial Cells"
                ),
                "n_annotation_source_cells": int(annotation.n_obs),
                "n_selected_cells": int(len(selected_ids)),
                "n_samples_in_annotation_source": int(len(all_samples)),
                "n_samples_with_selected_cells": int(len(selected_samples)),
                "n_samples_without_selected_cells": int(len(absent_samples)),
                "samples_without_selected_cells": ";".join(absent_samples),
                "n_ids_missing_from_qc": int(len(missing_ids)),
                "all_selected_ids_unique": bool(selected_ids.is_unique),
                "all_selected_ids_match_qc": len(missing_ids) == 0,
                "seed": SEED,
                "code_file": str(CODE_PATH),
            }
        ]
    ).to_csv(TABLE_DIR / "lineage_selection_parameters.csv", index=False)

    report = {
        "pass": True,
        "target_leiden_coarse": TARGET_LABEL,
        "rank_filter_mode": "strict_consistent",
        "n_selected_cells": int(len(selected_ids)),
        "n_samples_in_strict_atlas": int(len(all_samples)),
        "n_samples_with_selected_cells": int(len(selected_samples)),
        "samples_without_selected_cells": absent_samples,
        "n_ids_missing_from_qc": int(len(missing_ids)),
        "selected_ids_csv": str(TABLE_DIR / "selected_cell_ids.csv"),
    }
    (TABLE_DIR / "lineage_selection_completion.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    (TABLE_DIR / "readme.txt").write_text(
        f"""BRCA strict epithelial lineage selection

Annotation/ID source: {ANNOTATION_H5AD}
Expression/reclustering base: {QC_H5AD}
Selection: {COARSE_KEY} == {TARGET_LABEL!r} in the strict-consistent object.
Selected cells: {len(selected_ids)}
Samples represented: {len(selected_samples)} of {len(all_samples)}
Samples with zero strict-consistent epithelial cells: {', '.join(absent_samples)}

All selected IDs match adata_qc exactly. No cells were downsampled. The absent
samples were not manually restored because the user accepted strict filtering.
""",
        encoding="utf-8",
    )
    (TABLE_DIR / "package_versions.txt").write_text(
        "\n".join(
            [
                f"python={sys.version.split()[0]}",
                f"anndata={package_version('anndata')}",
                f"pandas={package_version('pandas')}",
                f"code={CODE_PATH}",
                f"seed={SEED}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    annotation.file.close()
    qc.file.close()
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
