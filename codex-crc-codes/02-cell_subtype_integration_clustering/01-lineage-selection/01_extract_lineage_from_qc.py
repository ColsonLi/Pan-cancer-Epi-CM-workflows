#!/usr/bin/env python3
"""Extract one score/rank-consistent broad lineage from adata_qc."""

from __future__ import annotations

import argparse
import platform
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

WORKFLOW_ROOT = Path("/mnt/disk18t/lr_xcy/riku/crc_val/epi-cm-core-workflow")
MODULE = "02-cell_subtype_integration_clustering"
STEP = "01-lineage-selection"
QC_H5AD = WORKFLOW_ROOT / "h5ad/01-celltype_integration_clustering/03-qc/adata_qc.h5ad"
ANNOTATION_H5AD = (
    WORKFLOW_ROOT
    / "h5ad/01-celltype_integration_clustering/07-score-rank-qc/adata_anno_score_genes_rank_consistent.h5ad"
)
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lineage", required=True, help="Output lineage slug, e.g. myeloid.")
    parser.add_argument("--target-label", required=True, help="Exact obs['leiden_coarse'] label.")
    parser.add_argument("--allow-existing", action="store_true")
    return parser.parse_args()


def output_paths(lineage: str) -> tuple[Path, Path]:
    output_h5ad = WORKFLOW_ROOT / "h5ad" / MODULE / STEP / lineage / f"adata_{lineage}_qc.h5ad"
    table_dir = WORKFLOW_ROOT / "tables" / MODULE / STEP / lineage
    return output_h5ad, table_dir


def assert_no_overwrite(output_h5ad: Path, table_dir: Path, allow_existing: bool) -> None:
    outputs = [
        output_h5ad,
        table_dir / "lineage_selection_parameters.csv",
        table_dir / "lineage_selection_summary.csv",
        table_dir / "selected_cell_ids.csv",
        table_dir / "readme.txt",
        table_dir / "package_versions.txt",
    ]
    existing = [str(path) for path in outputs if path.exists()]
    if existing and not allow_existing:
        raise FileExistsError("Lineage selection output already exists; refusing to overwrite:\n" + "\n".join(existing))


def validate_existing(output_h5ad: Path, table_dir: Path, lineage: str, target_label: str) -> None:
    if not output_h5ad.exists():
        raise FileNotFoundError(output_h5ad)
    adata = ad.read_h5ad(output_h5ad, backed="r")
    ok = "leiden_coarse" in adata.obs.columns and set(adata.obs["leiden_coarse"].astype(str).unique()) == {target_label}
    n_obs, n_vars = adata.shape
    adata.file.close()
    if not ok:
        raise RuntimeError(f"Existing lineage h5ad failed label validation: {output_h5ad}")
    summary = pd.DataFrame(
        [
            {
                "lineage": lineage,
                "target_leiden_coarse": target_label,
                "output_h5ad": str(output_h5ad),
                "n_selected_cells": int(n_obs),
                "n_genes": int(n_vars),
                "status": "completed_existing",
            }
        ]
    )
    table_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(table_dir / "lineage_selection_existing_validation.csv", index=False)
    print(summary.to_string(index=False))


def main() -> None:
    args = parse_args()
    output_h5ad, table_dir = output_paths(args.lineage)
    if args.allow_existing and output_h5ad.exists():
        validate_existing(output_h5ad, table_dir, args.lineage, args.target_label)
        return
    assert_no_overwrite(output_h5ad, table_dir, args.allow_existing)

    output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    anno = ad.read_h5ad(ANNOTATION_H5AD, backed="r")
    if "leiden_coarse" not in anno.obs.columns:
        raise KeyError("Missing obs['leiden_coarse'] in annotation h5ad")
    if "score_rank_consistent" not in anno.obs.columns:
        raise KeyError("Missing obs['score_rank_consistent'] in annotation h5ad")

    mask = (anno.obs["leiden_coarse"].astype(str) == args.target_label) & anno.obs["score_rank_consistent"].astype(bool)
    selected_ids = anno.obs_names[mask].copy()
    transfer_cols = [col for col in TRANSFER_COLUMNS if col in anno.obs.columns]
    transfer = anno.obs.loc[selected_ids, transfer_cols].copy()
    anno.file.close()

    qc = ad.read_h5ad(QC_H5AD)
    missing_ids = selected_ids.difference(qc.obs_names)
    if len(missing_ids) > 0:
        raise RuntimeError(f"{len(missing_ids)} selected IDs are absent from adata_qc.")
    lineage_adata = qc[selected_ids].copy()
    del qc

    for col in transfer_cols:
        lineage_adata.obs[col] = transfer.loc[lineage_adata.obs_names, col].values
    lineage_adata.obs["lineage_selection_source"] = str(ANNOTATION_H5AD)
    lineage_adata.obs["lineage_selection_label"] = args.target_label
    lineage_adata.obs["lineage_selection_rank_filter"] = "score_rank_consistent"
    lineage_adata.write_h5ad(output_h5ad)

    pd.DataFrame({"cell_id": selected_ids}).to_csv(table_dir / "selected_cell_ids.csv", index=False)
    summary = pd.DataFrame(
        [
            {
                "lineage": args.lineage,
                "target_leiden_coarse": args.target_label,
                "annotation_source_h5ad": str(ANNOTATION_H5AD),
                "expression_source_h5ad": str(QC_H5AD),
                "output_h5ad": str(output_h5ad),
                "n_selected_cells": int(lineage_adata.n_obs),
                "n_genes": int(lineage_adata.n_vars),
                "rank_filter_mode": "required_consistent_object",
                "all_selected_score_rank_consistent": bool(lineage_adata.obs["score_rank_consistent"].astype(bool).all()),
                "obs_names_unique": bool(lineage_adata.obs_names.is_unique),
                "status": "completed",
            }
        ]
    )
    summary.to_csv(table_dir / "lineage_selection_summary.csv", index=False)

    params = summary.copy()
    params["transfer_columns"] = ";".join(transfer_cols)
    params["code_file"] = str(CODE_FILE)
    params.to_csv(table_dir / "lineage_selection_parameters.csv", index=False)

    with (table_dir / "package_versions.txt").open("w") as fh:
        fh.write(f"python: {platform.python_version()}\n")
        fh.write(f"anndata: {ad.__version__}\n")
        fh.write(f"numpy: {np.__version__}\n")
        fh.write(f"pandas: {pd.__version__}\n")
        fh.write(f"code_file: {CODE_FILE}\n")

    with (table_dir / "readme.txt").open("w") as fh:
        fh.write(f"{args.lineage} lineage selection completed.\n")
        fh.write(f"Annotation/ID source: {ANNOTATION_H5AD}\n")
        fh.write(f"Expression/reclustering source: {QC_H5AD}\n")
        fh.write(f"Target label: {args.target_label}\n")
        fh.write("Rank filter: selected IDs must come from adata_anno_score_genes_rank_consistent.h5ad.\n")
        fh.write(f"Output h5ad: {output_h5ad}\n")

    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
