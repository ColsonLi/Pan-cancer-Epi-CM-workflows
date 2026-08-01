#!/usr/bin/env python3
"""Validate the manually selected BRCA broad-clustering handoff."""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


ROOT = Path("/mnt/disk18t/lr_xcy/riku/brca_cm-epi_codex")
WORKFLOW = ROOT / "epi-cm-core-workflow"
BLOCK = "01-celltype_integration_clustering"
TASK = "05-clustering-parameter-search"
HARMONY = WORKFLOW / "h5ad" / BLOCK / "04-integration-harmony" / "adata_harmony.h5ad"
SELECTED = WORKFLOW / "h5ad" / BLOCK / TASK / "selected" / "adata_inte.h5ad"
TABLE_ROOT = WORKFLOW / "tables" / BLOCK / TASK
SELECTED_DIR = TABLE_ROOT / "selected"
SELECTED_RECORD = TABLE_ROOT / "selected_clustering.csv"
GRID_MANIFEST = TABLE_ROOT / "clustering_grid_manifest.csv"
GRID_COUNTS = TABLE_ROOT / "pcs-20_nn-30_res-0p1-1p0" / "cluster_counts.csv"
PDF_HASH_AUDIT = TABLE_ROOT / "independent_pdf_panel_audit.csv"
FIGURE_ROOT = WORKFLOW / "figures" / BLOCK / TASK
SELECTED_FIGURE = FIGURE_ROOT / "selected" / "umap_selected_leiden_res0p8.pdf"
DIAGNOSTIC_FIGURE = FIGURE_ROOT / "selected" / "umap_selected_cluster_series_status.pdf"
OUTPUT = SELECTED_DIR / "selected_validation.json"
CLUSTER_KEY = "leiden_res0p8"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    record = pd.read_csv(SELECTED_RECORD)
    if len(record) != 1:
        raise AssertionError("selected_clustering.csv must contain exactly one row.")
    row = record.iloc[0]
    expected = {
        "integration_source": "Harmony",
        "use_rep": "X_pca_inte",
        "algorithm": "leiden",
        "final_raw_cluster_column": CLUSTER_KEY,
    }
    for key, value in expected.items():
        if str(row[key]) != value:
            raise AssertionError(f"Selected record {key} != {value}")
    if int(row["n_pcs"]) != 20 or int(row["n_neighbors"]) != 30:
        raise AssertionError("Selected PC/NN values are not 20/30.")
    if not np.isclose(float(row["resolution"]), 0.8):
        raise AssertionError("Selected resolution is not 0.8.")
    if str(row["grid_candidate_status"]) != "completed":
        raise AssertionError("Selected candidate is not recorded as completed.")
    if not bool(row["manual_selection"]):
        raise AssertionError("Selected record is not marked as a manual user selection.")

    manifest = pd.read_csv(GRID_MANIFEST)
    manifest_row = manifest[
        manifest["n_pcs"].eq(20)
        & manifest["n_neighbors"].eq(30)
        & np.isclose(manifest["resolution"], 0.8)
    ]
    if len(manifest_row) != 1 or str(manifest_row.iloc[0]["status"]) != "completed":
        raise AssertionError("Selected candidate is not completed in the live manifest.")
    grid_row = pd.read_csv(GRID_COUNTS)
    grid_row = grid_row[np.isclose(grid_row["resolution"], 0.8)]
    if len(grid_row) != 1 or int(grid_row.iloc[0]["cluster_count"]) != 15:
        raise AssertionError("Selected grid candidate does not have the expected 15 clusters.")

    harmony = ad.read_h5ad(HARMONY, backed="r")
    selected = ad.read_h5ad(SELECTED, backed="r")
    if selected.shape != harmony.shape or selected.shape != (99866, 3000):
        raise AssertionError(f"Selected shape mismatch: {selected.shape}")
    if not selected.obs_names.equals(harmony.obs_names):
        raise AssertionError("Selected cell IDs/order differ from the Harmony source.")
    if not selected.obs_names.is_unique:
        raise AssertionError("Selected cell IDs are not unique.")
    if selected.raw is None or selected.raw.shape != (99866, 27716):
        raise AssertionError("Selected raw normalized/log expression is missing or changed.")
    for key in ["sample", "series", "status", "original_barcode", CLUSTER_KEY]:
        if key not in selected.obs:
            raise AssertionError(f"Required selected obs column is missing: {key}")
    if selected.obs[CLUSTER_KEY].isna().any():
        raise AssertionError("Selected raw cluster column contains missing values.")
    n_clusters = int(selected.obs[CLUSTER_KEY].nunique())
    if n_clusters != 15:
        raise AssertionError(f"Selected rerun has {n_clusters} rather than 15 clusters.")
    if "leiden_coarse" in selected.obs or "cell_type" in selected.obs:
        raise AssertionError("Selected pure clustering prematurely contains annotation columns.")
    if "X_umap" not in selected.obsm or selected.obsm["X_umap"].shape != (99866, 2):
        raise AssertionError("Selected UMAP is missing or malformed.")
    if not np.isfinite(np.asarray(selected.obsm["X_umap"])).all():
        raise AssertionError("Selected UMAP has non-finite values.")
    if not {"distances", "connectivities"}.issubset(selected.obsp.keys()):
        raise AssertionError("Selected neighbor graph matrices are missing.")
    if "neighbors" not in selected.uns or CLUSTER_KEY not in selected.uns:
        raise AssertionError("Selected graph/Leiden provenance is missing from .uns.")

    size_table = pd.read_csv(SELECTED_DIR / "selected_cluster_sizes.csv")
    observed_sizes = (
        selected.obs[CLUSTER_KEY]
        .value_counts(sort=False)
        .rename_axis(CLUSTER_KEY)
        .rename("n_cells")
        .reset_index()
    )
    size_table[CLUSTER_KEY] = size_table[CLUSTER_KEY].astype(str)
    observed_sizes[CLUSTER_KEY] = observed_sizes[CLUSTER_KEY].astype(str)
    size_table = size_table.sort_values(CLUSTER_KEY).reset_index(drop=True)
    observed_sizes = observed_sizes.sort_values(CLUSTER_KEY).reset_index(drop=True)
    if not size_table.equals(observed_sizes):
        raise AssertionError("Saved selected cluster sizes do not match the h5ad.")

    for figure in [SELECTED_FIGURE, DIAGNOSTIC_FIGURE]:
        info = subprocess.run(
            ["pdfinfo", str(figure)], check=True, capture_output=True, text=True
        ).stdout
        if "Pages:           1" not in info:
            raise AssertionError(f"Selected PDF is invalid: {figure}")
    selected_text = subprocess.run(
        ["pdftotext", str(SELECTED_FIGURE), "-"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    if CLUSTER_KEY not in selected_text:
        raise AssertionError("Selected raw cluster PDF lacks the cluster title.")

    # Prove that adding the selected outputs did not alter any candidate PDF.
    old_hashes = pd.read_csv(PDF_HASH_AUDIT).set_index("graph_label")["sha256"].to_dict()
    candidate_pdfs = sorted(
        path
        for path in FIGURE_ROOT.rglob("umap_leiden_grid.pdf")
        if "selected" not in path.relative_to(FIGURE_ROOT).parts
    )
    if len(candidate_pdfs) != 81:
        raise AssertionError(f"Expected 81 preserved candidate PDFs, found {len(candidate_pdfs)}")
    for figure in candidate_pdfs:
        label = figure.parent.name
        if old_hashes.get(label) != sha256(figure):
            raise AssertionError(f"Candidate PDF changed after selected rerun: {label}")

    raster_files = [
        path
        for path in WORKFLOW.rglob("*")
        if path.is_file() and path.suffix.lower() in {".png", ".jpg", ".jpeg"}
    ]
    if raster_files:
        raise AssertionError(f"Raster files are forbidden: {raster_files}")

    report = {
        "pass": True,
        "manual_selection": True,
        "n_pcs": 20,
        "n_neighbors": 30,
        "resolution": 0.8,
        "raw_cluster_column": CLUSTER_KEY,
        "n_cells": int(selected.n_obs),
        "n_clusters": n_clusters,
        "raw_genes": int(selected.raw.n_vars),
        "all_cell_ids_and_order_match_harmony": True,
        "annotation_columns_absent_before_broad_annotation": True,
        "selected_pdf_valid": True,
        "candidate_pdfs_preserved_with_identical_hashes": 81,
        "raster_files": 0,
        "selected_h5ad_bytes": SELECTED.stat().st_size,
    }
    OUTPUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
