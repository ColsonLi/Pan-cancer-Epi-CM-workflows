#!/usr/bin/env python3
"""MT/Ribo metrics, gene filtering, and final cell QC."""

from __future__ import annotations

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
INPUT_H5AD = (
    WORKFLOW_ROOT
    / "h5ad/01-celltype_integration_clustering/03-qc/adata_doublet_filtered.h5ad"
)
OUTPUT_H5AD = WORKFLOW_ROOT / "h5ad/01-celltype_integration_clustering/03-qc/adata_qc.h5ad"
TABLE_DIR = WORKFLOW_ROOT / "tables/01-celltype_integration_clustering/03-qc"
FIG_DIR = WORKFLOW_ROOT / "figures/01-celltype_integration_clustering/03-qc"
CODE_FILE = Path(__file__)

MIN_CELLS_GENE = 3
MIN_COUNTS_GENE = 1
MIN_GENES_FINAL = 200
MAX_GENES_FINAL = 5000
MAX_PCT_MT = 15.0
MT_PREFIX = "MT-"
RIBO_PREFIX = "RPS"


def _write_params(name: str, row: dict) -> None:
    pd.DataFrame([row]).to_csv(TABLE_DIR / name, index=False)


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    if OUTPUT_H5AD.exists():
        raise FileExistsError(f"{OUTPUT_H5AD} exists; cleanup must be explicit before rerun.")

    sc.settings.figdir = str(FIG_DIR)
    sc.settings.set_figure_params(figsize=(3, 3), dpi=150)
    sc.settings.autoshow = False
    sc.settings.verbosity = 2

    adata = ad.read_h5ad(INPUT_H5AD)
    n_obs0, n_vars0 = adata.n_obs, adata.n_vars

    adata.var["MT"] = adata.var_names.astype(str).str.startswith(MT_PREFIX)
    adata.var["RIBO"] = adata.var_names.astype(str).str.startswith(RIBO_PREFIX)
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["MT", "RIBO"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    adata.obs["pct_counts_MT"] = adata.obs["pct_counts_MT"].astype(float)
    adata.obs["pct_counts_RIBO"] = adata.obs["pct_counts_RIBO"].astype(float)
    _write_params(
        "04_mt_ribo_qc_metric_parameters.csv",
        {
            "step": "mt_ribo_qc_metrics_before_gene_filter",
            "input_h5ad": str(INPUT_H5AD),
            "output_h5ad_or_object": "in_memory_adata",
            "n_obs_before": n_obs0,
            "n_obs_after": adata.n_obs,
            "n_vars_before": n_vars0,
            "n_vars_after": adata.n_vars,
            "mt_gene_prefix": MT_PREFIX,
            "ribo_gene_prefix": RIBO_PREFIX,
            "backend_package": "scanpy",
            "code_file": str(CODE_FILE),
            "random_seed": SEED,
            "notes": "QC metrics computed before gene filtering and recalculated after gene filtering.",
        },
    )

    try:
        sc.pl.violin(
            adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_MT", "pct_counts_RIBO", "doublet_score"],
            groupby=None,
            jitter=0.4,
            multi_panel=True,
            save="_qc_metrics_prefilter.pdf",
            show=False,
        )
    except Exception as exc:
        (TABLE_DIR / "qc_violin_prefilter_error.txt").write_text(repr(exc) + "\n", encoding="utf-8")

    n_obs_before_gene, n_vars_before_gene = adata.n_obs, adata.n_vars
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS_GENE)
    after_min_cells_vars = adata.n_vars
    sc.pp.filter_genes(adata, min_counts=MIN_COUNTS_GENE)
    _write_params(
        "05_gene_filter_parameters.csv",
        {
            "step": "gene_filter",
            "input_h5ad": str(INPUT_H5AD),
            "output_h5ad_or_object": "in_memory_adata_after_gene_filter",
            "n_obs_before": n_obs_before_gene,
            "n_obs_after": adata.n_obs,
            "n_vars_before": n_vars_before_gene,
            "n_vars_after": adata.n_vars,
            "min_cells_gene": MIN_CELLS_GENE,
            "n_vars_after_min_cells": after_min_cells_vars,
            "min_counts_gene": MIN_COUNTS_GENE,
            "backend_package": "scanpy",
            "code_file": str(CODE_FILE),
            "random_seed": SEED,
        },
    )

    adata.var["MT"] = adata.var_names.astype(str).str.startswith(MT_PREFIX)
    adata.var["RIBO"] = adata.var_names.astype(str).str.startswith(RIBO_PREFIX)
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["MT", "RIBO"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    n_obs_before_final, n_vars_before_final = adata.n_obs, adata.n_vars
    keep = (
        (adata.obs["n_genes_by_counts"] >= MIN_GENES_FINAL)
        & (adata.obs["n_genes_by_counts"] <= MAX_GENES_FINAL)
        & (adata.obs["pct_counts_MT"] <= MAX_PCT_MT)
    )
    adata = adata[keep.to_numpy()].copy()
    adata.write_h5ad(OUTPUT_H5AD, compression="lzf")
    _write_params(
        "06_final_cell_filter_parameters.csv",
        {
            "step": "final_cell_filter",
            "input_h5ad": str(INPUT_H5AD),
            "output_h5ad_or_object": str(OUTPUT_H5AD),
            "n_obs_before": n_obs_before_final,
            "n_obs_after": adata.n_obs,
            "n_vars_before": n_vars_before_final,
            "n_vars_after": adata.n_vars,
            "min_genes_final": MIN_GENES_FINAL,
            "max_genes_final": MAX_GENES_FINAL,
            "max_pct_mt": MAX_PCT_MT,
            "backend_package": "scanpy",
            "code_file": str(CODE_FILE),
            "random_seed": SEED,
        },
    )

    qc_report = pd.DataFrame(
        [
            {
                "metric": "input_after_doublet_filter",
                "n_obs": n_obs0,
                "n_vars": n_vars0,
            },
            {
                "metric": "after_gene_filter_before_final_cell_filter",
                "n_obs": n_obs_before_final,
                "n_vars": n_vars_before_final,
            },
            {
                "metric": "final_qc",
                "n_obs": adata.n_obs,
                "n_vars": adata.n_vars,
            },
            {
                "metric": "n_samples_final",
                "n_obs": adata.obs["sample"].nunique(),
                "n_vars": np.nan,
            },
            {
                "metric": "n_series_final",
                "n_obs": adata.obs["series"].nunique(),
                "n_vars": np.nan,
            },
        ]
    )
    qc_report.to_csv(TABLE_DIR / "qc_report.csv", index=False)
    adata.obs.groupby(["series", "sample"], observed=True).size().reset_index(name="n_cells_after_qc").to_csv(
        TABLE_DIR / "qc_cells_by_sample.csv", index=False
    )

    try:
        sc.pl.violin(
            adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_MT", "pct_counts_RIBO", "doublet_score"],
            groupby=None,
            jitter=0.4,
            multi_panel=True,
            save="_qc_metrics_final.pdf",
            show=False,
        )
    except Exception as exc:
        (TABLE_DIR / "qc_violin_final_error.txt").write_text(repr(exc) + "\n", encoding="utf-8")

    with (TABLE_DIR / "readme.txt").open("a", encoding="utf-8") as fh:
        fh.write("\nMT/Ribo metrics, gene filtering, and final cell QC completed.\n")
        fh.write(f"Input: {INPUT_H5AD}\n")
        fh.write(f"Output: {OUTPUT_H5AD}\n")
        fh.write(
            f"Final thresholds: min_genes={MIN_GENES_FINAL}, max_genes={MAX_GENES_FINAL}, max_pct_mt={MAX_PCT_MT}.\n"
        )

    versions = [
        f"python={platform.python_version()}",
        f"platform={platform.platform()}",
        f"anndata={ad.__version__}",
        f"scanpy={sc.__version__}",
        f"numpy={np.__version__}",
        f"pandas={pd.__version__}",
        "environment=/mnt/disk18t/lr_xcy/riku/crc_val/uv_envs/main/.venv",
        f"code_file={CODE_FILE}",
        f"random_seed={SEED}",
    ]
    (TABLE_DIR / "package_versions_final_qc.txt").write_text(
        "\n".join(versions) + "\n", encoding="utf-8"
    )
    print(qc_report.to_string(index=False))


if __name__ == "__main__":
    main()
