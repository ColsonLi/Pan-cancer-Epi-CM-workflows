#!/usr/bin/env python3
"""Independently validate the complete BRCA clustering grid and its artifacts."""

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
GRID = "05-clustering-parameter-search"
TABLE_ROOT = WORKFLOW / "tables" / BLOCK / GRID
FIGURE_ROOT = WORKFLOW / "figures" / BLOCK / GRID
HARMONY = (
    WORKFLOW
    / "h5ad"
    / BLOCK
    / "04-integration-harmony"
    / "adata_harmony.h5ad"
)
SELECTED_H5AD = (
    WORKFLOW
    / "h5ad"
    / BLOCK
    / GRID
    / "selected"
    / "adata_inte.h5ad"
)
MANIFEST = TABLE_ROOT / "clustering_grid_manifest.csv"
REVIEW = TABLE_ROOT / "candidate_review_manifest.csv"
COMPLETION = TABLE_ROOT / "clustering_grid_completion_check.csv"
SELECTED_RECORD = TABLE_ROOT / "selected_clustering.csv"
GRAPH_AUDIT = TABLE_ROOT / "independent_graph_rule_audit.csv"
PDF_AUDIT = TABLE_ROOT / "independent_pdf_panel_audit.csv"
SUMMARY = TABLE_ROOT / "independent_grid_audit.json"

PCS_VALUES = list(range(10, 51, 5))
NN_VALUES = list(range(10, 51, 5))
RESOLUTIONS = [round(x / 10, 1) for x in range(1, 11)]


def res_token(value: float) -> str:
    return f"{value:.1f}".replace(".", "p")


def graph_label(n_pcs: int, n_neighbors: int) -> str:
    return f"pcs-{n_pcs}_nn-{n_neighbors}_res-0p1-1p0"


def assert_same_float_or_nan(left: pd.Series, right: pd.Series, label: str) -> None:
    a = pd.to_numeric(left, errors="coerce").to_numpy(dtype=float)
    b = pd.to_numeric(right, errors="coerce").to_numpy(dtype=float)
    if not np.allclose(a, b, equal_nan=True):
        raise AssertionError(f"Numeric mismatch: {label}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    expected_combos = pd.MultiIndex.from_product(
        [PCS_VALUES, NN_VALUES, RESOLUTIONS],
        names=["n_pcs", "n_neighbors", "resolution"],
    )
    manifest = pd.read_csv(MANIFEST).sort_values(
        ["n_pcs", "n_neighbors", "resolution"]
    ).reset_index(drop=True)
    observed_combos = pd.MultiIndex.from_frame(
        manifest[["n_pcs", "n_neighbors", "resolution"]]
    )
    if len(manifest) != 810 or not observed_combos.equals(expected_combos):
        raise AssertionError("Manifest does not contain the exact 9 x 9 x 10 grid.")
    if manifest.duplicated(["n_pcs", "n_neighbors", "resolution"]).any():
        raise AssertionError("Manifest contains duplicate parameter combinations.")
    terminal = {"completed", "skipped_user_approved"}
    if not set(manifest["status"]).issubset(terminal):
        raise AssertionError("Manifest contains non-terminal or failed statuses.")

    harmony = ad.read_h5ad(HARMONY, backed="r")
    source_n_cells = int(harmony.n_obs)
    if source_n_cells != 99866:
        raise AssertionError(f"Unexpected Harmony cell count: {source_n_cells}")

    graph_rows: list[dict[str, object]] = []
    pdf_rows: list[dict[str, object]] = []
    for n_pcs in PCS_VALUES:
        for n_neighbors in NN_VALUES:
            label = graph_label(n_pcs, n_neighbors)
            table_dir = TABLE_ROOT / label
            figure = FIGURE_ROOT / label / "umap_leiden_grid.pdf"
            counts_path = table_dir / "cluster_counts.csv"
            params_path = table_dir / "clustering_parameters.csv"
            readme_path = table_dir / "readme.txt"
            for path in [figure, counts_path, params_path, readme_path]:
                if not path.is_file() or path.stat().st_size == 0:
                    raise AssertionError(f"Missing or empty graph artifact: {path}")

            counts = pd.read_csv(counts_path).sort_values("resolution").reset_index(drop=True)
            params = pd.read_csv(params_path).sort_values("resolution").reset_index(drop=True)
            if len(counts) != 10 or len(params) != 10:
                raise AssertionError(f"Expected 10 rows in graph tables: {label}")
            if not np.allclose(counts["resolution"], RESOLUTIONS):
                raise AssertionError(f"Resolution sequence mismatch: {label}")
            if not counts["status"].isin(terminal).all():
                raise AssertionError(f"Non-terminal graph status: {label}")
            if not counts["status"].equals(params["status"]):
                raise AssertionError(f"Counts/parameters status mismatch: {label}")
            assert_same_float_or_nan(
                counts["cluster_count"], params["cluster_count"], f"params {label}"
            )

            statuses = counts["status"].tolist()
            completed_idx = [i for i, status in enumerate(statuses) if status == "completed"]
            skipped_idx = [
                i for i, status in enumerate(statuses) if status == "skipped_user_approved"
            ]
            if not completed_idx or completed_idx != list(range(len(completed_idx))):
                raise AssertionError(f"Completed resolutions are not a nonempty prefix: {label}")
            if skipped_idx != list(range(len(completed_idx), 10)):
                raise AssertionError(f"Skipped resolutions are not the suffix: {label}")

            completed = counts.iloc[completed_idx].copy()
            completed_counts = completed["cluster_count"].to_numpy(dtype=float)
            if not np.isfinite(completed_counts).all():
                raise AssertionError(f"Completed cluster count is missing: {label}")
            if not np.equal(completed_counts, completed_counts.astype(int)).all():
                raise AssertionError(f"Completed cluster count is not integral: {label}")
            if skipped_idx and counts.iloc[skipped_idx]["cluster_count"].notna().any():
                raise AssertionError(f"Skipped cluster count should be empty: {label}")

            over = np.flatnonzero(completed_counts > 20)
            if len(over):
                if int(over[0]) != completed_idx[-1]:
                    raise AssertionError(
                        f"First >20 result is not the last executed resolution: {label}"
                    )
                if (completed_counts[: over[0]] > 20).any():
                    raise AssertionError(f"Earlier >20 result was not cut off: {label}")
                cutoff = float(completed.iloc[int(over[0])]["resolution"])
            else:
                if len(completed_idx) != 10:
                    raise AssertionError(f"Graph skipped without an executed >20 result: {label}")
                cutoff = np.nan
            if skipped_idx and counts.iloc[skipped_idx]["reason"].fillna("").eq("").any():
                raise AssertionError(f"Skipped rows lack reasons: {label}")

            manifest_graph = manifest[
                manifest["n_pcs"].eq(n_pcs)
                & manifest["n_neighbors"].eq(n_neighbors)
            ].sort_values("resolution").reset_index(drop=True)
            if not manifest_graph["status"].equals(counts["status"]):
                raise AssertionError(f"Manifest/graph status mismatch: {label}")
            assert_same_float_or_nan(
                manifest_graph["cluster_count"],
                counts["cluster_count"],
                f"manifest {label}",
            )

            if not params["graph_template_all_cells_preserved"].astype(bool).all():
                raise AssertionError(f"All-cell preservation flag failed: {label}")
            if params["candidate_h5ad_saved"].astype(bool).any():
                raise AssertionError(f"Candidate h5ad was marked as saved: {label}")
            if not params["backend"].eq("rapids_singlecell_gpu").all():
                raise AssertionError(f"Unexpected grid backend: {label}")

            pdf_info = subprocess.run(
                ["pdfinfo", str(figure)],
                check=True,
                capture_output=True,
                text=True,
            ).stdout
            if "Pages:           1" not in pdf_info:
                raise AssertionError(f"PDF is unreadable or not one page: {figure}")
            pdf_text = subprocess.run(
                ["pdftotext", str(figure), "-"],
                check=True,
                capture_output=True,
                text=True,
            ).stdout
            expected_titles = {
                f"leiden_res{res_token(value)}" for value in completed["resolution"]
            }
            observed_titles = {
                token
                for token in pdf_text.split()
                if token.startswith("leiden_res")
            }
            if observed_titles != expected_titles:
                raise AssertionError(
                    f"PDF panel titles do not equal executed resolutions: {label}; "
                    f"expected={sorted(expected_titles)}, observed={sorted(observed_titles)}"
                )

            graph_rows.append(
                {
                    "n_pcs": n_pcs,
                    "n_neighbors": n_neighbors,
                    "graph_label": label,
                    "n_executed": len(completed_idx),
                    "n_skipped": len(skipped_idx),
                    "first_resolution_over_20": cutoff,
                    "last_executed_resolution": float(completed["resolution"].iloc[-1]),
                    "last_executed_cluster_count": int(completed_counts[-1]),
                    "prefix_suffix_rule_valid": True,
                    "dynamic_cutoff_rule_valid": True,
                    "all_cells_preserved": True,
                    "gpu_backend": True,
                }
            )
            pdf_rows.append(
                {
                    "graph_label": label,
                    "pdf": str(figure),
                    "bytes": figure.stat().st_size,
                    "sha256": sha256(figure),
                    "n_expected_panels": len(expected_titles),
                    "n_observed_panel_titles": len(observed_titles),
                    "panel_set_exact": True,
                    "pdfinfo_valid": True,
                }
            )

    graph_audit = pd.DataFrame(graph_rows)
    pdf_audit = pd.DataFrame(pdf_rows)
    if len(graph_audit) != 81 or len(pdf_audit) != 81:
        raise AssertionError("Independent graph/PDF audit did not cover 81 graphs.")
    if int(graph_audit["n_executed"].sum()) != 701:
        raise AssertionError("Independent executed-resolution total is not 701.")
    if int(graph_audit["n_skipped"].sum()) != 109:
        raise AssertionError("Independent skipped-resolution total is not 109.")

    review = pd.read_csv(REVIEW)
    completion = pd.read_csv(COMPLETION)
    if len(review) != 81 or review["graph_label"].nunique() != 81:
        raise AssertionError("Candidate review manifest does not cover 81 unique graphs.")
    if len(completion) != 1:
        raise AssertionError("Completion check should contain one summary row.")
    check = completion.iloc[0]
    required_completion = {
        "expected_graph_candidates": 81,
        "completed_graph_candidates": 81,
        "expected_resolution_candidates": 810,
        "completed_resolution_candidates": 701,
        "skipped_resolution_candidates": 109,
        "failed_candidates": 0,
        "missing_candidates": 0,
    }
    for key, expected in required_completion.items():
        if int(check[key]) != expected:
            raise AssertionError(f"Completion field {key} != {expected}")
    if not bool(check["all_expected_lightweight_outputs_exist"]):
        raise AssertionError("Completion check reports missing lightweight outputs.")
    if not bool(check["manual_selection_ready"]):
        raise AssertionError("Completion check is not ready for manual review.")

    # The manual selected rerun is a valid post-grid output under a selected/
    # subdirectory.  Candidate counts must exclude those selected figures.
    grid_pdfs = sorted(
        path
        for path in FIGURE_ROOT.rglob("*.pdf")
        if "selected" not in path.relative_to(FIGURE_ROOT).parts
    )
    raster_files = sorted(
        path
        for path in WORKFLOW.rglob("*")
        if path.is_file() and path.suffix.lower() in {".png", ".jpg", ".jpeg"}
    )
    grid_h5ads = sorted(
        path
        for path in WORKFLOW.rglob("*.h5ad")
        if GRID in path.parts and path != SELECTED_H5AD
    )
    selected_record_exists = SELECTED_RECORD.is_file()
    selected_h5ad_exists = SELECTED_H5AD.is_file()
    if selected_record_exists != selected_h5ad_exists:
        raise AssertionError(
            "Manual selected record and selected h5ad must either both exist or both be absent."
        )
    manual_selection_present = selected_record_exists and selected_h5ad_exists
    if len(grid_pdfs) != 81:
        raise AssertionError(f"Expected exactly 81 grid PDFs, found {len(grid_pdfs)}")
    if raster_files:
        raise AssertionError(f"Raster figures are forbidden: {raster_files}")
    if grid_h5ads:
        raise AssertionError(f"Candidate grid h5ad files are forbidden: {grid_h5ads}")

    graph_audit.to_csv(GRAPH_AUDIT, index=False)
    pdf_audit.to_csv(PDF_AUDIT, index=False)
    summary = {
        "pass": True,
        "source_harmony_cells": source_n_cells,
        "expected_pc_values": PCS_VALUES,
        "expected_nn_values": NN_VALUES,
        "expected_resolutions": RESOLUTIONS,
        "manifest_rows": len(manifest),
        "unique_graphs": len(graph_audit),
        "executed_resolutions": int(graph_audit["n_executed"].sum()),
        "skipped_resolutions": int(graph_audit["n_skipped"].sum()),
        "failed_or_missing_resolutions": 0,
        "valid_pdfs": len(pdf_audit),
        "pdf_panel_sets_exact": bool(pdf_audit["panel_set_exact"].all()),
        "dynamic_cutoff_rules_valid": bool(
            graph_audit["dynamic_cutoff_rule_valid"].all()
        ),
        "all_cells_preserved": bool(graph_audit["all_cells_preserved"].all()),
        "gpu_backend_all_graphs": bool(graph_audit["gpu_backend"].all()),
        "raster_files": len(raster_files),
        "candidate_grid_h5ads": len(grid_h5ads),
        "manual_selected_record_present": manual_selection_present,
        "manual_selected_h5ad_present": manual_selection_present,
        "automatic_final_candidate_selected": False,
    }
    SUMMARY.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
