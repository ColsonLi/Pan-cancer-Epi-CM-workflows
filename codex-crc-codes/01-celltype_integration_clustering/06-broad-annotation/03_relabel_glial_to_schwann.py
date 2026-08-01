#!/usr/bin/env python3
"""Relabel broad annotation from Glial Cells to Schwann Cells."""

from __future__ import annotations

import platform
from pathlib import Path

import anndata as ad
import matplotlib
import pandas as pd
import scanpy as sc

matplotlib.use("Agg")

WORKFLOW_ROOT = Path("/mnt/disk18t/lr_xcy/riku/crc_val/epi-cm-core-workflow")
STEP = "01-celltype_integration_clustering/06-broad-annotation"
H5AD = WORKFLOW_ROOT / "h5ad" / STEP / "adata_anno.h5ad"
TABLE_DIR = WORKFLOW_ROOT / "tables" / STEP
FIGURE_DIR = WORKFLOW_ROOT / "figures" / STEP
ROUND2_DEG_DIR = TABLE_DIR / "degs_leiden_coarse_pcs30_nn30_res0p3"
CODE_FILE = Path(__file__)

OLD = "Glial Cells"
NEW = "Schwann Cells"
CATEGORY_ORDER = [
    "Epithelial Cells",
    "T Cells",
    "Cycling Immune Cells",
    "Myeloid Cells",
    "B Cells",
    "Plasma Cells",
    "Endothelial Cells",
    "Stromal Cells",
    "Schwann Cells",
    "Mast Cells",
]

MARKER_PANEL = {
    "Epithelial Cells": ["EPCAM", "KRT8", "KRT18"],
    "T Cells": ["CD3D", "CD3E", "CD2"],
    "Cycling Immune Cells": ["MKI67", "TOP2A", "STMN1"],
    "Myeloid Cells": ["LYZ", "LST1", "TYROBP"],
    "B Cells": ["MS4A1", "CD79A", "CD74"],
    "Plasma Cells": ["JCHAIN", "MZB1", "XBP1"],
    "Endothelial Cells": ["PECAM1", "VWF", "CLDN5"],
    "Stromal Cells": ["DCN", "COL1A1", "COL1A2"],
    "Schwann Cells": ["S100B", "PLP1", "PMP22"],
    "Mast Cells": ["TPSAB1", "CPA3", "KIT"],
}


def replace_label_series(series: pd.Series) -> pd.Categorical:
    values = series.astype(str).replace({OLD: NEW})
    return pd.Categorical(values, categories=CATEGORY_ORDER, ordered=True)


def update_csv_labels(path: Path) -> None:
    if not path.exists():
        return
    df = pd.read_csv(path)
    changed = False
    for col in ("leiden_coarse", "cell_type", "marker_class"):
        if col in df.columns:
            df[col] = df[col].astype(str).replace({OLD: NEW})
            changed = True
    if "evidence" in df.columns:
        df["evidence"] = df["evidence"].astype(str).str.replace(
            "glial/Schwann-like markers", "Schwann-cell markers", regex=False
        )
        changed = True
    if "annotation_note" in df.columns:
        df["annotation_note"] = df["annotation_note"].fillna("")
        df.loc[df.get("raw_cluster", pd.Series(dtype=str)).astype(str) == "15", "annotation_note"] = (
            "Schwann-cell label supported by S100B, PLP1, PMP22, MPZ, and SOX2."
        )
        changed = True
    if changed:
        df.to_csv(path, index=False)


def main() -> None:
    adata = ad.read_h5ad(H5AD)
    for col in ("leiden_coarse", "cell_type"):
        if col not in adata.obs.columns:
            raise KeyError(f"Missing obs[{col!r}]")
        adata.obs[col] = replace_label_series(adata.obs[col])

    update_csv_labels(TABLE_DIR / "02_raw_cluster_to_leiden_coarse_mapping.csv")
    update_csv_labels(TABLE_DIR / "02_broad_annotation_summary.csv")
    update_csv_labels(TABLE_DIR / "02_raw_cluster_broad_marker_summary.csv")
    update_csv_labels(TABLE_DIR / "02_round2_leiden_coarse_deg_summary.csv")

    old_deg = ROUND2_DEG_DIR / "Glial_Cells_degs_leiden_coarse_pcs30_nn30_res0p3.csv"
    new_deg = ROUND2_DEG_DIR / "Schwann_Cells_degs_leiden_coarse_pcs30_nn30_res0p3.csv"
    if old_deg.exists() and not new_deg.exists():
        df = pd.read_csv(old_deg)
        df.to_csv(new_deg, index=False)

    summary = (
        adata.obs["leiden_coarse"]
        .value_counts()
        .rename_axis("leiden_coarse")
        .reset_index(name="n_cells")
        .sort_values("leiden_coarse")
    )
    summary.to_csv(TABLE_DIR / "02_broad_annotation_summary.csv", index=False)

    marker_symbols = {
        label: [gene for gene in genes if gene in set(adata.raw.var["gene_symbol"].dropna().astype(str))]
        for label, genes in MARKER_PANEL.items()
    }
    sc.settings.autoshow = False
    sc.settings.figdir = str(FIGURE_DIR)
    sc.settings.set_figure_params(figsize=(3, 3), dpi=150)
    sc.pl.umap(
        adata,
        color=["leiden_res0p3", "leiden_coarse"],
        ncols=2,
        wspace=0.4,
        save="_raw_leiden_and_leiden_coarse.pdf",
        show=False,
    )
    sc.pl.dotplot(
        adata,
        var_names=marker_symbols,
        groupby="leiden_coarse",
        categories_order=CATEGORY_ORDER,
        standard_scale="var",
        use_raw=True,
        gene_symbols="gene_symbol",
        save="_broad_markers_leiden_coarse.pdf",
        show=False,
    )
    adata.write_h5ad(H5AD)

    pd.DataFrame(
        [
            {
                "step": "relabel_glial_to_schwann",
                "h5ad": str(H5AD),
                "old_label": OLD,
                "new_label": NEW,
                "n_schwann_cells": int((adata.obs["leiden_coarse"].astype(str) == NEW).sum()),
                "reason": "User approved renaming Glial Cells to Schwann Cells based on S100B, PLP1, PMP22, MPZ, and SOX2 marker evidence.",
                "code_file": str(CODE_FILE),
                "python": platform.python_version(),
                "scanpy": sc.__version__,
                "anndata": ad.__version__,
            }
        ]
    ).to_csv(TABLE_DIR / "03_relabel_glial_to_schwann.csv", index=False)


if __name__ == "__main__":
    main()
