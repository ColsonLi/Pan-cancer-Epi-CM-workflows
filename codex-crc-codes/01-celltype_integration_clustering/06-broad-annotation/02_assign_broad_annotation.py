#!/usr/bin/env python3
"""Assign broad labels, export Round 2 DEGs, and save annotated atlas."""

from __future__ import annotations

import platform
import random
import traceback
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
EXTERNAL_ROOT = Path("/mnt/disk18t/lr_xcy/riku/external_atlas_for_cm-epi_validation/CRC_S-BIAD2208")
STEP = "01-celltype_integration_clustering/06-broad-annotation"
INPUT_H5AD = (
    WORKFLOW_ROOT
    / "h5ad/01-celltype_integration_clustering/05-clustering-parameter-search/selected/adata_inte.h5ad"
)
OUTPUT_H5AD = WORKFLOW_ROOT / "h5ad" / STEP / "adata_anno.h5ad"
GTF_VAR_TABLE = (
    EXTERNAL_ROOT
    / "input_datasets_extracted/load_datasets/harmonize_datasets/artifacts/adata_var_gtf.csv"
)
TABLE_DIR = WORKFLOW_ROOT / "tables" / STEP
FIGURE_DIR = WORKFLOW_ROOT / "figures" / STEP
ROUND2_DEG_DIR = TABLE_DIR / "degs_leiden_coarse_pcs30_nn30_res0p3"
CODE_FILE = Path(__file__)

RAW_CLUSTER_COL = "leiden_res0p3"
ANNOTATION_COL = "leiden_coarse"
CELL_TYPE_COL = "cell_type"
METHOD = "t-test"
USE_RAW = True

CLUSTER_TO_LABEL = {
    "0": "Stromal Cells",
    "1": "Epithelial Cells",
    "2": "T Cells",
    "3": "Epithelial Cells",
    "4": "B Cells",
    "5": "T Cells",
    "6": "Plasma Cells",
    "7": "Endothelial Cells",
    "8": "Epithelial Cells",
    "9": "Myeloid Cells",
    "10": "Epithelial Cells",
    "11": "Mast Cells",
    "12": "Stromal Cells",
    "13": "Epithelial Cells",
    "14": "Epithelial Cells",
    "15": "Schwann Cells",
    "16": "Cycling Immune Cells",
}

ANNOTATION_EVIDENCE = {
    "0": "CALD1, DCN, COL1A2, COL6A2, LUM, COL3A1 fibroblast/stromal markers.",
    "1": "FABP1, SLC26A3, KRT8, EPCAM, KRT20, CLDN4 epithelial/colon epithelial markers.",
    "2": "NKG7, GZMA, CD3D, CD3E, CD8A cytotoxic T-cell markers.",
    "3": "EPCAM, KRT8, KRT18, KRT19, CLDN4 epithelial markers.",
    "4": "HLA-DRA, CD74, CD79A, MS4A1, BANK1 B-cell markers.",
    "5": "CD3D, IL7R, CD2, CD3E, CD3G T-cell markers.",
    "6": "JCHAIN, MZB1, XBP1, TNFRSF17 plasma-cell markers.",
    "7": "VWF, PLVAP, SPARCL1, RAMP2, CLDN5 endothelial markers.",
    "8": "PHGR1, KRT8, KRT18, EPCAM, TSPAN8 epithelial markers.",
    "9": "TYROBP, FCER1G, LYZ, LST1, AIF1 myeloid markers.",
    "10": "TFF3, AGR2, KRT18, KRT8, CLCA1, EPCAM epithelial/goblet markers.",
    "11": "TPSAB1, CPA3, MS4A2, KIT mast-cell markers.",
    "12": "ACTA2, TAGLN, MYL9, CALD1, COL1A2 smooth muscle/pericyte stromal markers.",
    "13": "CA7, BEST4, FABP1, KRT8, EPCAM epithelial markers.",
    "14": "POU2F3, TRPM5, SH2D6 with KRT8/KRT18 tuft-like epithelial evidence.",
    "15": "S100B, PLP1, PMP22, MPZ, SOX2 Schwann-cell markers.",
    "16": "MKI67, TOP2A, STMN1 cycling genes with immune markers PTPRC, CD52, CORO1A.",
}

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


def assert_no_overwrite() -> None:
    status_path = TABLE_DIR / "02_broad_annotation_status.csv"
    failed_resume = False
    if status_path.exists():
        try:
            status = pd.read_csv(status_path)
            failed_resume = str(status.iloc[0].get("step", "")) == "broad_annotation_failed"
        except Exception:
            failed_resume = False
    outputs = [
        OUTPUT_H5AD,
        TABLE_DIR / "02_raw_cluster_to_leiden_coarse_mapping.csv",
        TABLE_DIR / "02_broad_annotation_parameters.csv",
        TABLE_DIR / "02_broad_annotation_summary.csv",
        TABLE_DIR / "package_versions_broad_annotation.txt",
        TABLE_DIR / "readme_broad_annotation_complete.txt",
        FIGURE_DIR / "umap_raw_leiden_and_leiden_coarse.pdf",
        FIGURE_DIR / "dotplot__broad_markers_leiden_coarse.pdf",
    ]
    resumable_after_failure = {
        TABLE_DIR / "02_raw_cluster_to_leiden_coarse_mapping.csv",
        TABLE_DIR / "02_broad_annotation_status.csv",
        FIGURE_DIR / "umap_raw_leiden_and_leiden_coarse.pdf",
        TABLE_DIR / "02_missing_broad_marker_genes.csv",
    }
    existing = [
        str(path)
        for path in outputs
        if path.exists() and not (failed_resume and path in resumable_after_failure)
    ]
    if ROUND2_DEG_DIR.exists() and any(ROUND2_DEG_DIR.glob("*.csv")):
        existing.append(str(ROUND2_DEG_DIR))
    if existing:
        raise FileExistsError(
            "Broad annotation output already exists; refusing to overwrite:\n"
            + "\n".join(existing)
        )


def load_gene_maps() -> tuple[dict[str, str], dict[str, str]]:
    gtf = pd.read_csv(GTF_VAR_TABLE, usecols=["ensembl", "GeneSymbol"])
    gtf = gtf.dropna(subset=["ensembl", "GeneSymbol"])
    ensembl_to_symbol = dict(
        zip(gtf.drop_duplicates("ensembl")["ensembl"].astype(str), gtf.drop_duplicates("ensembl")["GeneSymbol"].astype(str))
    )
    symbol_to_ensembl = dict(
        zip(gtf.drop_duplicates("GeneSymbol")["GeneSymbol"].astype(str), gtf.drop_duplicates("GeneSymbol")["ensembl"].astype(str))
    )
    return ensembl_to_symbol, symbol_to_ensembl


def add_gene_symbols(adata: ad.AnnData, ensembl_to_symbol: dict[str, str]) -> None:
    adata.var["gene_symbol"] = pd.Series(adata.var_names.astype(str), index=adata.var_names).map(ensembl_to_symbol)
    if adata.raw is not None:
        adata.raw.var["gene_symbol"] = pd.Series(adata.raw.var_names.astype(str), index=adata.raw.var_names).map(
            ensembl_to_symbol
        )


def export_rank_genes_groups(adata: ad.AnnData, ensembl_to_symbol: dict[str, str]) -> pd.DataFrame:
    result = adata.uns["rank_genes_groups"]
    groups = list(result["names"].dtype.names)
    summary_rows = []
    for group in groups:
        ensembl = pd.Series(result["names"][group].astype(str), name="ensembl")
        symbols = ensembl.map(ensembl_to_symbol)
        genes = symbols.fillna(ensembl)
        out = pd.DataFrame(
            {
                "gene": genes,
                "score": result["scores"][group],
                "logfoldchanges": result["logfoldchanges"][group],
                "pvals": result["pvals"][group],
                "pvals_adj": result["pvals_adj"][group],
                "ensembl": ensembl,
                "gene_symbol": symbols,
            }
        )
        safe_group = str(group).replace("/", "_").replace(" ", "_")
        out.to_csv(ROUND2_DEG_DIR / f"{safe_group}_degs_leiden_coarse_pcs30_nn30_res0p3.csv", index=False)
        summary_rows.append(
            {
                "leiden_coarse": group,
                "n_genes": int(out.shape[0]),
                "top_gene": str(out.loc[0, "gene"]),
                "top_gene_ensembl": str(out.loc[0, "ensembl"]),
                "top_score": float(out.loc[0, "score"]),
                "top_logfoldchanges": float(out.loc[0, "logfoldchanges"]),
                "top_pvals_adj": float(out.loc[0, "pvals_adj"]),
            }
        )
    return pd.DataFrame(summary_rows)


def marker_panel_symbols(adata: ad.AnnData) -> dict[str, list[str]]:
    raw_symbols = set(adata.raw.var["gene_symbol"].dropna().astype(str))
    converted = {}
    missing_rows = []
    for label, genes in MARKER_PANEL.items():
        present_symbols = []
        for gene in genes:
            if gene in raw_symbols:
                present_symbols.append(gene)
            else:
                missing_rows.append({"leiden_coarse": label, "marker_gene": gene})
        converted[label] = present_symbols
    pd.DataFrame(missing_rows).to_csv(TABLE_DIR / "02_missing_broad_marker_genes.csv", index=False)
    return converted


def main() -> None:
    assert_no_overwrite()
    OUTPUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    ROUND2_DEG_DIR.mkdir(parents=True, exist_ok=True)

    sc.settings.autoshow = False
    sc.settings.figdir = str(FIGURE_DIR)
    sc.settings.set_figure_params(figsize=(3, 3), dpi=150)

    status = {
        "step": "broad_annotation_start",
        "input_h5ad": str(INPUT_H5AD),
        "output_h5ad": str(OUTPUT_H5AD),
        "raw_cluster_col": RAW_CLUSTER_COL,
        "annotation_col": ANNOTATION_COL,
        "cell_type_col": CELL_TYPE_COL,
        "round2_deg_method": METHOD,
        "use_raw": USE_RAW,
        "gtf_var_table": str(GTF_VAR_TABLE),
        "code_file": str(CODE_FILE),
        "random_seed": SEED,
    }
    pd.DataFrame([status]).to_csv(TABLE_DIR / "02_broad_annotation_status.csv", index=False)

    try:
        ensembl_to_symbol, symbol_to_ensembl = load_gene_maps()
        adata = ad.read_h5ad(INPUT_H5AD)
        if adata.raw is None:
            raise RuntimeError("adata.raw is absent; broad annotation requires raw-normalized expression.")
        if RAW_CLUSTER_COL not in adata.obs.columns:
            raise KeyError(f"Missing obs[{RAW_CLUSTER_COL!r}] in selected h5ad.")
        add_gene_symbols(adata, ensembl_to_symbol)

        raw_clusters = adata.obs[RAW_CLUSTER_COL].astype(str)
        missing = sorted(set(raw_clusters.unique()) - set(CLUSTER_TO_LABEL))
        if missing:
            raise ValueError(f"Missing broad-label mapping for raw clusters: {missing}")

        adata.obs[ANNOTATION_COL] = raw_clusters.map(CLUSTER_TO_LABEL)
        adata.obs[ANNOTATION_COL] = pd.Categorical(adata.obs[ANNOTATION_COL], categories=CATEGORY_ORDER, ordered=True)
        adata.obs[CELL_TYPE_COL] = adata.obs[ANNOTATION_COL].copy()

        cluster_sizes = raw_clusters.value_counts().rename_axis("raw_cluster").reset_index(name="n_cells")
        mapping = cluster_sizes.sort_values("raw_cluster", key=lambda s: s.astype(int)).copy()
        mapping[ANNOTATION_COL] = mapping["raw_cluster"].map(CLUSTER_TO_LABEL)
        mapping[CELL_TYPE_COL] = mapping[ANNOTATION_COL]
        mapping["evidence"] = mapping["raw_cluster"].map(ANNOTATION_EVIDENCE)
        mapping["selected_review_depth"] = 120
        mapping["confidence"] = "high"
        mapping.loc[mapping["raw_cluster"].isin(["14", "16"]), "confidence"] = "medium"
        mapping.loc[mapping["raw_cluster"] == "14", "annotation_note"] = "Rare tuft-like epithelial cluster assigned to broad epithelial label."
        mapping.loc[mapping["raw_cluster"] == "16", "annotation_note"] = "Cycling cluster with mixed immune markers kept as cycling immune broad label."
        mapping.to_csv(TABLE_DIR / "02_raw_cluster_to_leiden_coarse_mapping.csv", index=False)

        sc.pl.umap(
            adata,
            color=[RAW_CLUSTER_COL, ANNOTATION_COL],
            ncols=2,
            wspace=0.4,
            save="_raw_leiden_and_leiden_coarse.pdf",
            show=False,
        )

        marker_symbols = marker_panel_symbols(adata)
        sc.pl.dotplot(
            adata,
            var_names=marker_symbols,
            groupby=ANNOTATION_COL,
            categories_order=CATEGORY_ORDER,
            standard_scale="var",
            use_raw=True,
            gene_symbols="gene_symbol",
            save="_broad_markers_leiden_coarse.pdf",
            show=False,
        )

        sc.tl.rank_genes_groups(
            adata,
            groupby=ANNOTATION_COL,
            method=METHOD,
            use_raw=USE_RAW,
        )
        round2_summary = export_rank_genes_groups(adata, ensembl_to_symbol)
        round2_summary.to_csv(TABLE_DIR / "02_round2_leiden_coarse_deg_summary.csv", index=False)

        summary = (
            adata.obs[ANNOTATION_COL]
            .value_counts()
            .rename_axis(ANNOTATION_COL)
            .reset_index(name="n_cells")
            .sort_values(ANNOTATION_COL)
        )
        summary.to_csv(TABLE_DIR / "02_broad_annotation_summary.csv", index=False)

        params = {
            **status,
            "step": "broad_annotation_complete",
            "n_obs": int(adata.n_obs),
            "n_vars": int(adata.n_vars),
            "raw_n_vars": int(adata.raw.n_vars),
            "n_raw_clusters": int(adata.obs[RAW_CLUSTER_COL].nunique()),
            "n_leiden_coarse": int(adata.obs[ANNOTATION_COL].nunique()),
            "round2_deg_dir": str(ROUND2_DEG_DIR),
            "expression_source_for_degs_and_dotplot": "adata.raw",
        }
        pd.DataFrame([params]).to_csv(TABLE_DIR / "02_broad_annotation_parameters.csv", index=False)
        pd.DataFrame([params]).to_csv(TABLE_DIR / "02_broad_annotation_status.csv", index=False)

        with (TABLE_DIR / "package_versions_broad_annotation.txt").open("w") as fh:
            fh.write(f"python: {platform.python_version()}\n")
            fh.write(f"anndata: {ad.__version__}\n")
            fh.write(f"scanpy: {sc.__version__}\n")
            fh.write(f"pandas: {pd.__version__}\n")
            fh.write(f"numpy: {np.__version__}\n")
            fh.write(f"code_file: {CODE_FILE}\n")

        with (TABLE_DIR / "readme_broad_annotation_complete.txt").open("w") as fh:
            fh.write("06-broad-annotation completed.\n")
            fh.write(f"Input selected h5ad: {INPUT_H5AD}\n")
            fh.write(f"Output annotated h5ad: {OUTPUT_H5AD}\n")
            fh.write(f"Raw cluster column: obs['{RAW_CLUSTER_COL}']\n")
            fh.write(f"Broad annotation column: obs['{ANNOTATION_COL}']\n")
            fh.write(f"Cell type column initialized from broad annotation: obs['{CELL_TYPE_COL}']\n")
            fh.write("Round 1 DEGs were exported before annotation; Round 2 DEGs were exported after leiden_coarse assignment.\n")
            fh.write("DEG and marker plotting used adata.raw with gene symbols mapped from adata_var_gtf.csv.\n")
            fh.write("Cluster 14 was assigned to Epithelial Cells with tuft-like epithelial evidence.\n")
            fh.write("Cluster 15 was assigned to Schwann Cells based on S100B, PLP1, PMP22, MPZ, and SOX2.\n")
            fh.write("Cluster 16 was assigned to Cycling Immune Cells due cycling genes plus immune markers.\n")

        adata.write_h5ad(OUTPUT_H5AD)
        print(pd.DataFrame([params]).to_string(index=False))

    except Exception as exc:
        status.update(
            {
                "step": "broad_annotation_failed",
                "error_summary": "".join(traceback.format_exception_only(type(exc), exc)).strip(),
            }
        )
        pd.DataFrame([status]).to_csv(TABLE_DIR / "02_broad_annotation_status.csv", index=False)
        raise


if __name__ == "__main__":
    main()
