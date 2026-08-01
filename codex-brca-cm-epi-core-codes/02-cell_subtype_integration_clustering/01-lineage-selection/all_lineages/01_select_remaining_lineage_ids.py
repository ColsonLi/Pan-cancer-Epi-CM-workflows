#!/usr/bin/env python3
"""Select strict-consistent IDs for all remaining BRCA broad lineages."""

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
TABLE_ROOT = WORKFLOW / "tables" / BLOCK02 / "01-lineage-selection"
CODE_PATH = (
    WORKFLOW
    / "codes"
    / BLOCK02
    / "01-lineage-selection"
    / "all_lineages"
    / "01_select_remaining_lineage_ids.py"
)
LINEAGES = {
    "t_cells": "T Cells",
    "myeloid": "Myeloid Cells",
    "b_cells": "B Cells",
    "plasma": "Plasma Cells",
    "endothelial": "Endothelial Cells",
    "stromal": "Stromal Cells",
    "perivascular": "Perivascular Cells",
}
COARSE_KEY = "leiden_coarse"
BEST_KEY = "best_rank_type_global"
CONSISTENT_KEY = "score_rank_consistent"


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def main() -> None:
    annotation = ad.read_h5ad(ANNOTATION_H5AD, backed="r")
    qc = ad.read_h5ad(QC_H5AD, backed="r")
    required = {
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
    missing = required - set(annotation.obs.columns)
    if missing:
        raise ValueError(f"Strict annotation source lacks columns: {sorted(missing)}")
    if not annotation.obs_names.is_unique or not qc.obs_names.is_unique:
        raise ValueError("Annotation or QC expression source has non-unique cell IDs.")
    if not annotation.obs[CONSISTENT_KEY].astype(bool).all():
        raise ValueError("Strict annotation source contains inconsistent cells.")
    if not annotation.obs[COARSE_KEY].astype(str).equals(
        annotation.obs[BEST_KEY].astype(str)
    ):
        raise ValueError("Strict annotation source violates best-rank equality.")

    all_samples = sorted(annotation.obs["sample"].astype(str).unique())
    qc_names = pd.Index(qc.obs_names.astype(str))
    score_columns = sorted(
        column
        for column in annotation.obs.columns
        if column.endswith("_score") or column.endswith("_score_rank_pct")
    )
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
        *score_columns,
    ]
    summary_rows: list[dict[str, object]] = []
    for slug, label in LINEAGES.items():
        table_dir = TABLE_ROOT / slug
        table_dir.mkdir(parents=True, exist_ok=True)
        mask = annotation.obs[COARSE_KEY].astype(str).eq(label)
        selected = annotation.obs.loc[mask, export_columns].copy()
        if selected.empty:
            raise ValueError(f"No strict-consistent cells found for {label}.")
        selected_ids = pd.Index(selected.index.astype(str))
        missing_ids = selected_ids.difference(qc_names)
        if len(missing_ids):
            raise ValueError(f"{label}: {len(missing_ids)} selected IDs missing from adata_qc.")
        selected_samples = sorted(selected["sample"].astype(str).unique())
        absent_samples = sorted(set(all_samples) - set(selected_samples))

        selected.insert(0, "cell_id", selected_ids)
        selected.to_csv(table_dir / "selected_cell_ids.csv", index=False)
        sample_counts = (
            selected.assign(sample=selected["sample"].astype(str))
            .groupby("sample", observed=True)
            .size()
            .rename("n_selected_lineage_cells")
            .reindex(all_samples, fill_value=0)
            .rename_axis("sample")
            .reset_index()
        )
        sample_counts["has_selected_lineage_cells"] = (
            sample_counts["n_selected_lineage_cells"] > 0
        )
        sample_counts.to_csv(table_dir / "selected_lineage_cells_by_sample.csv", index=False)
        parameter_row = {
            "lineage_slug": slug,
            "target_leiden_coarse": label,
            "annotation_source_h5ad": str(ANNOTATION_H5AD),
            "expression_base_h5ad": str(QC_H5AD),
            "rank_filter_mode": "strict_consistent",
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
        pd.DataFrame([parameter_row]).to_csv(
            table_dir / "lineage_selection_parameters.csv", index=False
        )
        report = {
            "pass": True,
            **parameter_row,
            "selected_ids_csv": str(table_dir / "selected_cell_ids.csv"),
        }
        (table_dir / "lineage_selection_completion.json").write_text(
            json.dumps(report, indent=2) + "\n", encoding="utf-8"
        )
        (table_dir / "readme.txt").write_text(
            f"""BRCA strict lineage selection: {label}

Annotation/ID source: {ANNOTATION_H5AD}
Expression/reclustering base: {QC_H5AD}
Selected cells: {len(selected_ids)}
Samples represented: {len(selected_samples)} of {len(all_samples)}
Samples with zero strict-consistent cells: {', '.join(absent_samples) or 'none'}

All selected IDs match adata_qc exactly. No cells were downsampled or restored.
""",
            encoding="utf-8",
        )
        (table_dir / "package_versions.txt").write_text(
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
        summary_rows.append(parameter_row)
        print(
            f"[Lineage selection] {slug}: cells={len(selected_ids)} "
            f"samples={len(selected_samples)}/{len(all_samples)}",
            flush=True,
        )

    summary_dir = TABLE_ROOT / "all_lineages"
    summary_dir.mkdir(parents=True, exist_ok=True)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(summary_dir / "remaining_lineages_selection_summary.csv", index=False)
    completion = {
        "pass": True,
        "n_remaining_lineages": len(summary_rows),
        "all_ids_match_qc": bool(summary["all_selected_ids_match_qc"].all()),
        "lineage_cell_counts": {
            row["lineage_slug"]: int(row["n_selected_cells"]) for row in summary_rows
        },
    }
    (summary_dir / "remaining_lineages_selection_completion.json").write_text(
        json.dumps(completion, indent=2) + "\n", encoding="utf-8"
    )
    annotation.file.close()
    qc.file.close()
    print(json.dumps(completion, indent=2))


if __name__ == "__main__":
    main()
