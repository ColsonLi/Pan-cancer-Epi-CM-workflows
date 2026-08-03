#!/usr/bin/env python3
"""Round 1 broad-annotation DEGs for the selected raw Leiden clusters."""

from __future__ import annotations

import platform
import random
import traceback
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

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
GTF_VAR_TABLE = (
    EXTERNAL_ROOT
    / "input_datasets_extracted/load_datasets/harmonize_datasets/artifacts/adata_var_gtf.csv"
)
TABLE_DIR = WORKFLOW_ROOT / "tables" / STEP
DEG_DIR = TABLE_DIR / "degs_leiden_res0p3_pcs30_nn30_res0p3"
CODE_FILE = Path(__file__)

GROUPBY = "leiden_res0p3"
METHOD = "t-test"
USE_RAW = True


def assert_no_overwrite() -> None:
    outputs = [
        TABLE_DIR / "01_round1_raw_leiden_deg_parameters.csv",
        TABLE_DIR / "01_round1_raw_leiden_deg_summary.csv",
        TABLE_DIR / "package_versions_round1_deg.txt",
        TABLE_DIR / "readme.txt",
    ]
    existing = [str(path) for path in outputs if path.exists()]
    if DEG_DIR.exists() and any(DEG_DIR.glob("*.csv")):
        existing.append(str(DEG_DIR))
    if existing:
        raise FileExistsError(
            "Round 1 broad-annotation DEG output already exists; refusing to overwrite:\n"
            + "\n".join(existing)
        )


def load_gene_symbol_map() -> dict[str, str]:
    gtf = pd.read_csv(GTF_VAR_TABLE, usecols=["ensembl", "GeneSymbol"])
    gtf = gtf.dropna(subset=["ensembl", "GeneSymbol"]).drop_duplicates("ensembl")
    return dict(zip(gtf["ensembl"].astype(str), gtf["GeneSymbol"].astype(str)))


def export_rank_genes_groups(adata: ad.AnnData, symbol_map: dict[str, str]) -> pd.DataFrame:
    result = adata.uns["rank_genes_groups"]
    groups = list(result["names"].dtype.names)
    summary_rows = []
    for group in groups:
        ensembl = pd.Series(result["names"][group].astype(str), name="ensembl")
        symbols = ensembl.map(symbol_map)
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
        safe_group = str(group).replace("/", "_")
        out.to_csv(DEG_DIR / f"{safe_group}_degs_leiden_res0p3_pcs30_nn30_res0p3.csv", index=False)
        summary_rows.append(
            {
                "group": group,
                "n_genes": int(out.shape[0]),
                "top_gene": str(out.loc[0, "gene"]),
                "top_gene_ensembl": str(out.loc[0, "ensembl"]),
                "top_score": float(out.loc[0, "score"]),
                "top_logfoldchanges": float(out.loc[0, "logfoldchanges"]),
                "top_pvals_adj": float(out.loc[0, "pvals_adj"]),
            }
        )
    return pd.DataFrame(summary_rows)


def main() -> None:
    assert_no_overwrite()
    DEG_DIR.mkdir(parents=True, exist_ok=True)

    status = {
        "step": "round1_raw_leiden_deg_start",
        "input_h5ad": str(INPUT_H5AD),
        "groupby": GROUPBY,
        "method": METHOD,
        "use_raw": USE_RAW,
        "gtf_var_table": str(GTF_VAR_TABLE),
        "code_file": str(CODE_FILE),
        "random_seed": SEED,
    }
    pd.DataFrame([status]).to_csv(TABLE_DIR / "01_round1_raw_leiden_deg_status.csv", index=False)

    try:
        symbol_map = load_gene_symbol_map()
        adata = ad.read_h5ad(INPUT_H5AD)
        if adata.raw is None:
            raise RuntimeError("adata.raw is absent; broad annotation DEG requires use_raw=True.")
        if GROUPBY not in adata.obs.columns:
            raise KeyError(f"Missing obs[{GROUPBY!r}] in selected clustering h5ad.")

        raw_symbols = pd.Series(adata.raw.var_names.astype(str), index=adata.raw.var_names).map(symbol_map)
        adata.raw.var["gene_symbol"] = raw_symbols
        var_symbols = pd.Series(adata.var_names.astype(str), index=adata.var_names).map(symbol_map)
        adata.var["gene_symbol"] = var_symbols

        sc.tl.rank_genes_groups(
            adata,
            groupby=GROUPBY,
            method=METHOD,
            use_raw=USE_RAW,
        )
        summary = export_rank_genes_groups(adata, symbol_map)
        summary.to_csv(TABLE_DIR / "01_round1_raw_leiden_deg_summary.csv", index=False)

        params = {
            **status,
            "step": "round1_raw_leiden_deg_complete",
            "n_obs": int(adata.n_obs),
            "n_vars": int(adata.n_vars),
            "raw_n_vars": int(adata.raw.n_vars),
            "n_groups": int(adata.obs[GROUPBY].nunique()),
            "deg_dir": str(DEG_DIR),
            "canonical_deg_columns": "gene;score;logfoldchanges;pvals;pvals_adj",
            "additional_columns": "ensembl;gene_symbol",
            "expression_source": "adata.raw",
        }
        pd.DataFrame([params]).to_csv(TABLE_DIR / "01_round1_raw_leiden_deg_parameters.csv", index=False)
        pd.DataFrame([params]).to_csv(TABLE_DIR / "01_round1_raw_leiden_deg_status.csv", index=False)

        with (TABLE_DIR / "package_versions_round1_deg.txt").open("w") as fh:
            fh.write(f"python: {platform.python_version()}\n")
            fh.write(f"anndata: {ad.__version__}\n")
            fh.write(f"scanpy: {sc.__version__}\n")
            fh.write(f"pandas: {pd.__version__}\n")
            fh.write(f"numpy: {np.__version__}\n")
            fh.write(f"code_file: {CODE_FILE}\n")

        with (TABLE_DIR / "readme.txt").open("w") as fh:
            fh.write("06-broad-annotation started with Round 1 raw Leiden DEG export.\n")
            fh.write(f"Input h5ad: {INPUT_H5AD}\n")
            fh.write(f"Raw cluster column: obs['{GROUPBY}']\n")
            fh.write("DEG method: scanpy.tl.rank_genes_groups(method='t-test', use_raw=True)\n")
            fh.write(f"DEG output directory: {DEG_DIR}\n")
            fh.write(f"Gene symbol mapping: {GTF_VAR_TABLE}\n")
            fh.write("No leiden_coarse labels are assigned in this script.\n")

        print(pd.DataFrame([params]).to_string(index=False))

    except Exception as exc:
        status.update(
            {
                "step": "round1_raw_leiden_deg_failed",
                "error_summary": "".join(traceback.format_exception_only(type(exc), exc)).strip(),
            }
        )
        pd.DataFrame([status]).to_csv(TABLE_DIR / "01_round1_raw_leiden_deg_status.csv", index=False)
        raise


if __name__ == "__main__":
    main()
