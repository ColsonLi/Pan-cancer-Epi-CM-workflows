#!/usr/bin/env python3
"""Build the BRCA merged count AnnData from the published GSE176078 matrix."""

from __future__ import annotations

import importlib.metadata
import json
import random
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy
from scipy import sparse
from scipy.io import mmread


SEED = 42
random.seed(SEED)
np.random.seed(SEED)

ANALYSIS_ROOT = Path("/mnt/disk18t/lr_xcy/riku/brca_cm-epi_codex")
DATA_ROOT = Path(
    "/mnt/disk18t/lr_xcy/riku/external_atlas_for_cm-epi_validation/"
    "Breast_Wu2021_Zenodo4739739/single_cell/merged"
)
WORKFLOW_ROOT = ANALYSIS_ROOT / "epi-cm-core-workflow"
BLOCK = "01-celltype_integration_clustering"
TABLE_DIR = WORKFLOW_ROOT / "tables" / BLOCK / "02-merge-metadata"
H5AD_DIR = WORKFLOW_ROOT / "h5ad" / BLOCK / "02-merge-metadata"
CODE_PATH = (
    WORKFLOW_ROOT
    / "codes"
    / BLOCK
    / "02-merge-metadata"
    / "01_prepare_merge.py"
)
OUT_H5AD = H5AD_DIR / "adata_merge.h5ad"

MATRIX = DATA_ROOT / "count_matrix_sparse.mtx"
GENES = DATA_ROOT / "count_matrix_genes.tsv"
BARCODES = DATA_ROOT / "count_matrix_barcodes.tsv"
METADATA = DATA_ROOT / "metadata.csv"


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def write_readme() -> None:
    text = f"""BRCA GSE176078 merge and metadata audit

Input count matrix: {MATRIX}
Input genes: {GENES}
Input barcodes: {BARCODES}
Input published metadata: {METADATA}
Code: {CODE_PATH}
Output: {OUT_H5AD}

The supplied combined processed count matrix is used directly because it already
contains all 26 biological samples from one GEO series. No cells are sampled or
downsampled. The matrix remains raw integer counts in .X. Published cell labels
are preserved as metadata but are not used to generate technical clusters.
All biological samples are recorded as Primary tumor; there is no independent
adjacent-tissue or healthy-normal sample in this dataset.
"""
    (TABLE_DIR / "readme.txt").write_text(text, encoding="utf-8")


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    H5AD_DIR.mkdir(parents=True, exist_ok=True)

    genes = pd.read_csv(GENES, sep="\t", header=None, dtype=str)[0]
    barcodes = pd.read_csv(BARCODES, sep="\t", header=None, dtype=str)[0]
    obs = pd.read_csv(METADATA, index_col=0)
    obs.index = obs.index.astype(str)

    if not genes.is_unique:
        raise ValueError("Published gene symbols are not unique.")
    if not barcodes.is_unique:
        raise ValueError("Published cell IDs are not unique.")
    if not obs.index.is_unique:
        raise ValueError("Published metadata cell IDs are not unique.")
    if not np.array_equal(barcodes.to_numpy(), obs.index.to_numpy()):
        raise ValueError("Barcode order does not exactly match metadata row order.")

    matrix_gene_by_cell = mmread(MATRIX)
    if matrix_gene_by_cell.shape != (len(genes), len(barcodes)):
        raise ValueError(
            f"Matrix shape {matrix_gene_by_cell.shape} does not match "
            f"genes/barcodes {(len(genes), len(barcodes))}."
        )
    X = sparse.csr_matrix(matrix_gene_by_cell.T, dtype=np.int32)
    del matrix_gene_by_cell

    total_counts = np.asarray(X.sum(axis=1)).ravel().astype(np.int64)
    n_features = X.getnnz(axis=1).astype(np.int64)
    expected_counts = obs["nCount_RNA"].to_numpy(dtype=np.int64)
    expected_features = obs["nFeature_RNA"].to_numpy(dtype=np.int64)
    if not np.array_equal(total_counts, expected_counts):
        bad = int(np.count_nonzero(total_counts != expected_counts))
        raise ValueError(f"Count totals disagree with published metadata for {bad} cells.")
    if not np.array_equal(n_features, expected_features):
        bad = int(np.count_nonzero(n_features != expected_features))
        raise ValueError(f"Detected-gene totals disagree with published metadata for {bad} cells.")

    obs["sample"] = obs["orig.ident"].astype(str)
    obs["series"] = "GSE176078"
    obs["status"] = "Primary tumor"
    obs["status_source"] = "study_design_and_user_confirmed_primary_tumor_only"
    obs["original_barcode"] = [
        cell_id[len(sample) + 1 :]
        if cell_id.startswith(f"{sample}_")
        else cell_id
        for cell_id, sample in zip(obs.index, obs["sample"], strict=True)
    ]
    if obs["original_barcode"].eq(obs.index).any():
        raise ValueError("At least one merged cell ID could not be split into sample and barcode.")
    obs["source_matrix_path"] = str(MATRIX)

    var = pd.DataFrame(index=pd.Index(genes.astype(str), name="gene"))
    var["gene_name"] = var.index.astype(str)
    var["gene_identifier_type"] = "gene_name"

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.uns["dataset"] = "GSE176078_Wu_etal_2021_BRCA_scRNASeq"
    adata.uns["expression_contract"] = {
        "X": "raw integer UMI counts",
        "raw_counts_preserved_in": ".X of adata_merge.h5ad",
        "matrix_orientation_source": "gene_by_cell MatrixMarket transposed to cell_by_gene",
    }
    adata.uns["random_seed"] = SEED

    sample_meta = (
        obs.groupby("sample", observed=True)
        .agg(n_cells=("sample", "size"), subtype=("subtype", "first"))
        .reset_index()
    )
    sample_meta.insert(1, "series", "GSE176078")
    sample_meta["status"] = "Primary tumor"
    sample_meta["tissue_site"] = "Breast"
    sample_meta["disease_context"] = "Primary breast cancer"
    sample_meta["treatment_context"] = "not_available_in_local_processed_metadata"
    sample_meta["include"] = True
    sample_meta["source_url_or_file"] = str(METADATA)
    sample_meta["notes"] = "No independent AT or healthy-normal biological sample."
    sample_meta.to_csv(TABLE_DIR / "sample_metadata.csv", index=False)

    mapping = pd.DataFrame(
        {
            "merged_cell_id": obs.index,
            "sample": obs["sample"].to_numpy(),
            "original_barcode": obs["original_barcode"].to_numpy(),
            "source_matrix_path": str(MATRIX),
        }
    )
    mapping.to_csv(
        TABLE_DIR / "merged_cell_id_mapping.csv.gz",
        index=False,
        compression="gzip",
    )

    gene_audit = pd.DataFrame(
        [
            {
                "input_h5ad": str(OUT_H5AD),
                "object_layer": "var_names",
                "var_names_type": "gene_name",
                "n_features": len(genes),
                "n_unique": int(genes.nunique()),
                "n_duplicates": int(genes.duplicated().sum()),
                "ensembl_like_n": int(
                    genes.str.match(r"^ENS[A-Z]*G[0-9]+(?:\.[0-9]+)?$").sum()
                ),
                "ensembl_like_fraction": float(
                    genes.str.match(r"^ENS[A-Z]*G[0-9]+(?:\.[0-9]+)?$").mean()
                ),
                "candidate_gene_id_columns": "",
                "candidate_gene_name_columns": "gene_name",
                "selected_gene_id_source": "",
                "selected_gene_name_source": "count_matrix_genes.tsv",
                "action_taken": "preserved_unique_gene_symbols_as_var_names",
            }
        ]
    )
    gene_audit.to_csv(TABLE_DIR / "gene_identifier_audit.csv", index=False)
    pd.DataFrame(
        {"object_layer": "var_names", "example": genes.head(20).to_numpy()}
    ).to_csv(TABLE_DIR / "gene_identifier_examples.csv", index=False)

    validation = {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_nonzero": int(adata.X.nnz),
        "n_samples": int(obs["sample"].nunique()),
        "obs_names_unique": bool(adata.obs_names.is_unique),
        "var_names_unique": bool(adata.var_names.is_unique),
        "barcode_metadata_order_exact": True,
        "published_count_totals_exact": True,
        "published_feature_totals_exact": True,
        "raw_counts_dtype": str(adata.X.dtype),
        "random_seed": SEED,
    }
    (TABLE_DIR / "input_validation.json").write_text(
        json.dumps(validation, indent=2), encoding="utf-8"
    )

    versions = {
        "python": sys.version.split()[0],
        "anndata": package_version("anndata"),
        "numpy": package_version("numpy"),
        "pandas": package_version("pandas"),
        "scipy": scipy.__version__,
        "code": str(CODE_PATH),
        "seed": str(SEED),
    }
    (TABLE_DIR / "package_versions.txt").write_text(
        "\n".join(f"{k}={v}" for k, v in versions.items()) + "\n",
        encoding="utf-8",
    )
    write_readme()
    adata.write_h5ad(OUT_H5AD, compression="gzip")
    print(json.dumps(validation, indent=2))


if __name__ == "__main__":
    main()
