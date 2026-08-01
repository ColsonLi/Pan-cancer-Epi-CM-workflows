#!/usr/bin/env python3
"""Final validation for BRCA selected clustering, broad annotation, and score/rank QC."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


ROOT = Path("/mnt/disk18t/lr_xcy/riku/brca_cm-epi_codex")
WORKFLOW = ROOT / "epi-cm-core-workflow"
BLOCK = "01-celltype_integration_clustering"
SELECTED = WORKFLOW / "h5ad" / BLOCK / "05-clustering-parameter-search" / "selected" / "adata_inte.h5ad"
ANNOTATED = WORKFLOW / "h5ad" / BLOCK / "06-broad-annotation" / "adata_anno.h5ad"
SCORED = WORKFLOW / "h5ad" / BLOCK / "07-score-rank-qc" / "adata_anno_score_genes_rank.h5ad"
CONSISTENT = WORKFLOW / "h5ad" / BLOCK / "07-score-rank-qc" / "adata_anno_score_genes_rank_consistent.h5ad"
TABLE06 = WORKFLOW / "tables" / BLOCK / "06-broad-annotation"
TABLE07 = WORKFLOW / "tables" / BLOCK / "07-score-rank-qc"
FIG06 = WORKFLOW / "figures" / BLOCK / "06-broad-annotation"
FIG07 = WORKFLOW / "figures" / BLOCK / "07-score-rank-qc"
ROUND1 = TABLE06 / "degs_leiden_res0p8_pcs20_nn30_res0p8"
ROUND2 = TABLE06 / "degs_leiden_coarse_pcs20_nn30_res0p8_myo_merged"
MAPPING = TABLE06 / "broad_annotation_mapping.csv"
OUTPUT = TABLE07 / "module01_final_validation.json"
FLAGS = TABLE07 / "score_rank_retention_flags.csv"

RAW_KEY = "leiden_res0p8"
COARSE_KEY = "leiden_coarse"
BEST_KEY = "best_rank_type_global"
CONSISTENT_KEY = "score_rank_consistent"
EXPECTED_LABELS = [
    "Epithelial Cells",
    "T Cells",
    "Myeloid Cells",
    "B Cells",
    "Plasma Cells",
    "Endothelial Cells",
    "Stromal Cells",
    "Perivascular Cells",
]
EXPECTED_COUNTS = {
    "Epithelial Cells": 30658,
    "T Cells": 32884,
    "Myeloid Cells": 9239,
    "B Cells": 3422,
    "Plasma Cells": 3681,
    "Endothelial Cells": 7643,
    "Stromal Cells": 6935,
    "Perivascular Cells": 5404,
}


def validate_full_deg_dir(path: Path, expected_files: int) -> None:
    files = sorted(path.glob("*.csv"))
    if len(files) != expected_files:
        raise AssertionError(f"Expected {expected_files} DEG CSVs in {path}, found {len(files)}")
    for file in files:
        frame = pd.read_csv(file)
        if len(frame) != 27716:
            raise AssertionError(f"DEG table is not full length: {file}")
        canonical = ["gene", "score", "logfoldchanges", "pvals", "pvals_adj"]
        if frame.columns[:5].tolist() != canonical:
            raise AssertionError(f"DEG canonical column order is wrong: {file}")
        if frame["gene"].isna().any() or frame["gene"].duplicated().any():
            raise AssertionError(f"DEG genes are missing or duplicated: {file}")


def validate_pdfs(directory: Path, expected: int) -> list[str]:
    pdfs = sorted(directory.glob("*.pdf"))
    if len(pdfs) != expected:
        raise AssertionError(f"Expected {expected} PDFs in {directory}, found {len(pdfs)}")
    for pdf in pdfs:
        info = subprocess.run(
            ["pdfinfo", str(pdf)], check=True, capture_output=True, text=True
        ).stdout
        if "Pages:" not in info:
            raise AssertionError(f"Unreadable PDF: {pdf}")
    return [pdf.name for pdf in pdfs]


def main() -> None:
    validate_full_deg_dir(ROUND1, 15)
    validate_full_deg_dir(ROUND2, 8)

    selected = ad.read_h5ad(SELECTED, backed="r")
    annotated = ad.read_h5ad(ANNOTATED, backed="r")
    scored = ad.read_h5ad(SCORED, backed="r")
    filtered = ad.read_h5ad(CONSISTENT, backed="r")
    if selected.shape != annotated.shape or annotated.shape != scored.shape:
        raise AssertionError("Selected, annotated, and unfiltered scored shapes differ.")
    if scored.shape != (99866, 3000) or not (0 < filtered.n_obs <= scored.n_obs):
        raise AssertionError(f"Unexpected scored/filtered shapes: {scored.shape}, {filtered.shape}")
    if not selected.obs_names.equals(annotated.obs_names):
        raise AssertionError("Annotation changed selected cell IDs/order.")
    if not annotated.obs_names.equals(scored.obs_names):
        raise AssertionError("Scoring changed unfiltered cell IDs/order.")
    for obj, label in [(selected, "selected"), (annotated, "annotated"), (scored, "scored"), (filtered, "filtered")]:
        if not obj.obs_names.is_unique:
            raise AssertionError(f"{label} object cell IDs are not unique.")
        if obj.raw is None or obj.raw.n_vars != 27716:
            raise AssertionError(f"{label} object lacks the preserved full raw gene space.")
        for key in ["sample", "series", "status", "original_barcode", RAW_KEY]:
            if key not in obj.obs:
                raise AssertionError(f"{label} object lacks required obs column {key}")

    labels = annotated.obs[COARSE_KEY].cat.categories.astype(str).tolist()
    if labels != EXPECTED_LABELS:
        raise AssertionError(f"Broad category order mismatch: {labels}")
    counts = annotated.obs[COARSE_KEY].value_counts(sort=False).to_dict()
    if {str(k): int(v) for k, v in counts.items()} != EXPECTED_COUNTS:
        raise AssertionError(f"Broad category counts mismatch: {counts}")
    if not annotated.obs[COARSE_KEY].astype(str).equals(annotated.obs["cell_type"].astype(str)):
        raise AssertionError("Annotated cell_type is not initialized from leiden_coarse.")
    mapping = pd.read_csv(MAPPING, dtype={"raw_cluster": str})
    if len(mapping) != 15 or mapping["raw_cluster"].duplicated().any():
        raise AssertionError("Broad mapping does not cover 15 raw clusters exactly once.")
    if mapping.set_index("raw_cluster").loc["11", COARSE_KEY] != "Epithelial Cells":
        raise AssertionError("BRCA myoepithelial raw cluster 11 was not merged into Epithelial Cells.")
    if mapping.set_index("raw_cluster").loc["8", COARSE_KEY] != "Epithelial Cells":
        raise AssertionError("BRCA luminal-progenitor cluster 8 mapping unexpectedly changed.")

    score_mapping = pd.read_csv(TABLE07 / "broad_score_column_mapping.csv")
    gene_sets = pd.read_csv(TABLE07 / "broad_score_gene_sets.csv")
    if len(score_mapping) != 8 or set(score_mapping[COARSE_KEY]) != set(EXPECTED_LABELS):
        raise AssertionError("Score mapping does not match the eight observed broad labels.")
    gene_counts = gene_sets.groupby(COARSE_KEY, observed=True).size().to_dict()
    if gene_counts != {label: 100 for label in EXPECTED_LABELS}:
        raise AssertionError(f"Score gene-set sizes are not exactly 100: {gene_counts}")
    if not gene_sets["source_full_deg_csv"].map(lambda value: Path(value).parent == ROUND2).all():
        raise AssertionError("Score gene sets are not sourced only from saved round-2 full DEGs.")
    if not score_mapping["use_raw"].astype(bool).all():
        raise AssertionError("Score mapping does not record use_raw=True for every label.")
    score_cols = score_mapping["score_column"].tolist()
    rank_cols = score_mapping["rank_column"].tolist()
    for key in [COARSE_KEY, "cell_type", BEST_KEY, CONSISTENT_KEY, *score_cols, *rank_cols]:
        if key not in scored.obs:
            raise AssertionError(f"Unfiltered scored object lacks {key}")
    if scored.obs[[*score_cols, *rank_cols]].isna().any().any():
        raise AssertionError("Score/rank columns contain missing values.")
    expected_consistency = (
        scored.obs[BEST_KEY].astype(str) == scored.obs[COARSE_KEY].astype(str)
    )
    if not expected_consistency.equals(scored.obs[CONSISTENT_KEY].astype(bool)):
        raise AssertionError("Stored score/rank consistency calls do not match the rule.")
    expected_filtered_ids = scored.obs_names[expected_consistency]
    if not filtered.obs_names.equals(expected_filtered_ids):
        raise AssertionError("Filtered handoff IDs/order do not equal consistent cells.")
    if not filtered.obs[CONSISTENT_KEY].astype(bool).all():
        raise AssertionError("Filtered handoff contains inconsistent cells.")

    broad_summary = pd.read_csv(TABLE07 / "score_rank_consistency_by_leiden_coarse.csv")
    raw_summary = pd.read_csv(TABLE07 / "score_rank_consistency_by_raw_cluster.csv")
    if int(broad_summary["n_input"].sum()) != 99866 or int(broad_summary["n_kept"].sum()) != filtered.n_obs:
        raise AssertionError("Broad retention summary totals are inconsistent.")
    if int(raw_summary["n_input"].sum()) != 99866 or int(raw_summary["n_kept"].sum()) != filtered.n_obs:
        raise AssertionError("Raw-cluster retention summary totals are inconsistent.")
    flag_rows: list[dict[str, object]] = []
    for row in raw_summary.itertuples(index=False):
        if float(row.retained_fraction) < 0.5:
            raw_cluster = str(getattr(row, RAW_KEY))
            if raw_cluster == "8":
                note = (
                    "Luminal-progenitor/basal epithelial program has low broad score/rank agreement; "
                    "retain the unfiltered object for audit and request downstream input choice."
                )
            elif raw_cluster == "11":
                note = (
                    "Contractile myoepithelial raw cluster was merged into broad Epithelial Cells "
                    "by user request but has low merged-epithelial score/rank agreement; use the "
                    "unfiltered scored object to preserve it for epithelial subtype analysis."
                )
            else:
                note = (
                    "Small selected T-cell cluster has low score/rank agreement; "
                    "retained in the unfiltered object."
                )
            flag_rows.append(
                {
                    "level": "raw_cluster",
                    "value": str(getattr(row, RAW_KEY)),
                    "n_input": int(row.n_input),
                    "n_kept": int(row.n_kept),
                    "retained_fraction": float(row.retained_fraction),
                    "severity": "high",
                    "note": note,
                }
            )
    pd.DataFrame(flag_rows).to_csv(FLAGS, index=False)

    broad_pdfs = validate_pdfs(FIG06, 4)
    score_pdfs = validate_pdfs(FIG07, 4)
    raster_files = [
        path
        for path in WORKFLOW.rglob("*")
        if path.is_file() and path.suffix.lower() in {".png", ".jpg", ".jpeg"}
    ]
    if raster_files:
        raise AssertionError(f"Raster figures are forbidden: {raster_files}")

    report = {
        "pass": True,
        "dataset": "BRCA Breast_Wu2021 GSE176078",
        "selected_parameters": {"n_pcs": 20, "n_neighbors": 30, "resolution": 0.8},
        "selected_cells": int(selected.n_obs),
        "selected_raw_clusters": int(selected.obs[RAW_KEY].nunique()),
        "broad_labels": EXPECTED_COUNTS,
        "round1_full_deg_csvs": 15,
        "round2_full_deg_csvs": 8,
        "score_gene_sets": 8,
        "score_genes_per_label": 100,
        "score_genes_use_raw": True,
        "unfiltered_scored_cells": int(scored.n_obs),
        "consistent_handoff_cells": int(filtered.n_obs),
        "consistent_handoff_fraction": int(filtered.n_obs) / int(scored.n_obs),
        "retention_qc_flags": flag_rows,
        "unfiltered_object_preserved": True,
        "broad_annotation_pdfs": broad_pdfs,
        "score_rank_pdfs": score_pdfs,
        "raster_files": 0,
        "downstream_input_choice_required_before_subtype_work": bool(flag_rows),
    }
    OUTPUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
