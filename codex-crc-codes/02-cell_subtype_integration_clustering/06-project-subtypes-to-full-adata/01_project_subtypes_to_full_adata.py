#!/usr/bin/env python3
"""Project final lineage subtype annotations back to the consistent full atlas."""

from __future__ import annotations

import gc
import platform
import random
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
MODULE = "02-cell_subtype_integration_clustering"
STEP = "06-project-subtypes-to-full-adata"
CODE_FILE = Path(__file__)
TARGET_H5AD = (
    WORKFLOW_ROOT
    / "h5ad/01-celltype_integration_clustering/07-score-rank-qc/adata_anno_score_genes_rank_consistent.h5ad"
)
OUTPUT_H5AD = WORKFLOW_ROOT / "h5ad" / MODULE / STEP / "adata_anno_cellsubtype.h5ad"
TABLE_DIR = WORKFLOW_ROOT / "tables" / MODULE / STEP
FIGURE_DIR = WORKFLOW_ROOT / "figures" / MODULE / STEP

LINEAGES = {
    "b_cells": {"abbr": "b", "broad": "B Cells", "file": "adata_b.h5ad"},
    "cycling_immune": {"abbr": "cycling", "broad": "Cycling Immune Cells", "file": "adata_cycling.h5ad"},
    "endothelial": {"abbr": "endo", "broad": "Endothelial Cells", "file": "adata_endo.h5ad"},
    "epithelial": {"abbr": "epi", "broad": "Epithelial Cells", "file": "adata_epi.h5ad"},
    "mast": {"abbr": "mast", "broad": "Mast Cells", "file": "adata_mast.h5ad"},
    "myeloid": {"abbr": "mye", "broad": "Myeloid Cells", "file": "adata_mye.h5ad"},
    "plasma_cells": {"abbr": "plasma", "broad": "Plasma Cells", "file": "adata_plasma.h5ad"},
    "schwann": {"abbr": "schwann", "broad": "Schwann Cells", "file": "adata_schwann.h5ad"},
    "stromal": {"abbr": "stromal", "broad": "Stromal Cells", "file": "adata_stromal.h5ad"},
    "t_cells": {"abbr": "t", "broad": "T Cells", "file": "adata_t.h5ad"},
}

TRANSFER_COLUMNS = [
    "cell_subtype",
    "functional_state",
    "annotation_source",
    "annotation_confidence",
    "marker_evidence",
    "gene_selection_rationale",
    "annotation_note",
]

STATUS_MAP = {
    "primary_crc_no_metastasis_no_lymph_node": "tumor",
    "normal_colon_no_metastasis_no_lymph_node": "normal",
}


def main() -> None:
    if OUTPUT_H5AD.exists():
        raise FileExistsError(f"Projection output already exists: {OUTPUT_H5AD}")
    OUTPUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    full = ad.read_h5ad(TARGET_H5AD)
    full.obs["cell_type_before_subtype_projection"] = full.obs["cell_type"].astype(str).values if "cell_type" in full.obs else ""
    full.obs["cell_subtype_before_subtype_projection"] = (
        full.obs["cell_subtype"].astype(str).values if "cell_subtype" in full.obs else ""
    )
    full.obs["functional_state_before_subtype_projection"] = (
        full.obs["functional_state"].astype(str).values if "functional_state" in full.obs else ""
    )
    if "status" in full.obs:
        full.obs["status_original"] = full.obs["status"].astype(str).values
        full.obs["status"] = full.obs["status_original"].map(STATUS_MAP).fillna(full.obs["status_original"]).astype("category")

    projected = pd.DataFrame(index=full.obs_names)
    projected["cell_type"] = full.obs["leiden_coarse"].astype(str).values
    for col in TRANSFER_COLUMNS:
        projected[col] = pd.NA
    projected["annotation_source_h5ad"] = pd.NA
    projected["projection_lineage"] = pd.NA

    all_reports = []
    for lineage, cfg in LINEAGES.items():
        lineage_h5ad = WORKFLOW_ROOT / "h5ad" / MODULE / "05-subtype-deg-annotation" / lineage / cfg["file"]
        table_subdir = TABLE_DIR / lineage
        table_subdir.mkdir(parents=True, exist_ok=True)
        lineage_adata = ad.read_h5ad(lineage_h5ad, backed="r")
        obs = lineage_adata.obs.copy()
        lineage_adata.file.close()
        duplicated = int(obs.index.duplicated().sum())
        matched = obs.index.intersection(full.obs_names)
        unmatched = obs.index.difference(full.obs_names)
        target_mask = full.obs["leiden_coarse"].astype(str) == cfg["broad"]
        expected_ids = full.obs_names[target_mask]
        missing_expected = expected_ids.difference(obs.index)
        extra_off_broad = matched.difference(expected_ids)

        for col in TRANSFER_COLUMNS:
            if col in obs.columns:
                projected.loc[matched, col] = obs.loc[matched, col].astype(str).values
        projected.loc[matched, "cell_type"] = cfg["broad"]
        projected.loc[matched, "annotation_source_h5ad"] = str(lineage_h5ad)
        projected.loc[matched, "projection_lineage"] = lineage

        report = pd.DataFrame(
            [
                {
                    "lineage": lineage,
                    "broad_label": cfg["broad"],
                    "source_h5ad": str(lineage_h5ad),
                    "target_h5ad": str(TARGET_H5AD),
                    "source_cells": int(obs.shape[0]),
                    "target_expected_cells": int(target_mask.sum()),
                    "matched_cells": int(len(matched)),
                    "unmatched_source_cells": int(len(unmatched)),
                    "missing_expected_target_cells": int(len(missing_expected)),
                    "duplicated_source_obs_names": duplicated,
                    "matched_off_broad_cells": int(len(extra_off_broad)),
                }
            ]
        )
        report.to_csv(table_subdir / "projection_match_report.csv", index=False)
        pd.DataFrame({"cell_id": unmatched}).to_csv(table_subdir / "unmatched_source_cell_ids.csv", index=False)
        pd.DataFrame({"cell_id": missing_expected}).to_csv(table_subdir / "missing_expected_target_cell_ids.csv", index=False)
        pd.DataFrame({"cell_id": extra_off_broad}).to_csv(table_subdir / "matched_off_broad_cell_ids.csv", index=False)
        all_reports.append(report)

        line_fig_dir = FIGURE_DIR / lineage
        line_fig_dir.mkdir(parents=True, exist_ok=True)
        lineage_plot = ad.read_h5ad(lineage_h5ad)
        sc.settings.autoshow = False
        sc.settings.figdir = str(line_fig_dir)
        sc.settings.set_figure_params(figsize=(3, 3), dpi=150)
        sc.pl.umap(lineage_plot, color="cell_subtype", save="_cell_subtype.pdf", show=False)
        del lineage_plot, obs
        gc.collect()

    for col in ["cell_type", *TRANSFER_COLUMNS, "annotation_source_h5ad", "projection_lineage"]:
        full.obs[col] = projected[col].astype("category")
    full.obs["cell_type"] = full.obs["cell_type"].astype("category")
    full.write_h5ad(OUTPUT_H5AD)

    all_report_df = pd.concat(all_reports, ignore_index=True)
    all_report_df.to_csv(TABLE_DIR / "projection_match_report_all_lineages.csv", index=False)
    subtype_counts = full.obs["cell_subtype"].astype(str).value_counts(dropna=False).rename_axis("cell_subtype").reset_index(name="n_cells")
    subtype_counts.to_csv(TABLE_DIR / "projected_cell_subtype_counts.csv", index=False)
    broad_counts = full.obs[["leiden_coarse", "cell_type"]].astype(str).value_counts(dropna=False).reset_index(name="n_cells")
    broad_counts.to_csv(TABLE_DIR / "projected_broad_label_comparison_counts.csv", index=False)

    sc.settings.autoshow = False
    sc.settings.figdir = str(FIGURE_DIR)
    sc.settings.set_figure_params(figsize=(3, 3), dpi=150)
    sc.pl.umap(
        full,
        color=["leiden_coarse", "cell_type"],
        ncols=2,
        wspace=0.4,
        save="_leiden_coarse_vs_projected_cell_type.pdf",
        show=False,
    )
    sc.pl.umap(full, color="cell_subtype", save="_projected_cell_subtype.pdf", show=False)

    with (TABLE_DIR / "package_versions.txt").open("w") as fh:
        fh.write(f"python: {platform.python_version()}\n")
        fh.write(f"anndata: {ad.__version__}\n")
        fh.write(f"scanpy: {sc.__version__}\n")
        fh.write(f"numpy: {np.__version__}\n")
        fh.write(f"pandas: {pd.__version__}\n")
        fh.write(f"code_file: {CODE_FILE}\n")

    with (TABLE_DIR / "readme.txt").open("w") as fh:
        fh.write("Projected final lineage subtype annotations back to the score/rank-consistent full atlas.\n")
        fh.write(f"Target full atlas: {TARGET_H5AD}\n")
        fh.write(f"Projected output h5ad: {OUTPUT_H5AD}\n")
        fh.write("Projection matched by exact obs_names.\n")
        fh.write("Original cell_type/cell_subtype/functional_state columns were preserved in backup columns before projection.\n")
        fh.write("status was standardized to tumor/normal; original values preserved in status_original.\n")
        fh.write("SVG output for full-atlas UMAPs was not emitted to avoid very large vector files for 643k cells; PDF QC figures were saved.\n")

    print(all_report_df.to_string(index=False))
    print(f"output_h5ad={OUTPUT_H5AD}")


if __name__ == "__main__":
    main()
