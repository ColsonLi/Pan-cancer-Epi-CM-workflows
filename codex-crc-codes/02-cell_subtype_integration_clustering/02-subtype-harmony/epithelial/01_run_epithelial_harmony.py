#!/usr/bin/env python3
"""Run epithelial subtype normalization, PCA, and Harmony integration."""

from __future__ import annotations

import gc
import platform
import random
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

SEED = 42
random.seed(SEED)
np.random.seed(SEED)

WORKFLOW_ROOT = Path("/mnt/disk18t/lr_xcy/riku/crc_val/epi-cm-core-workflow")
MODULE = "02-cell_subtype_integration_clustering"
STEP = "02-subtype-harmony"
LINEAGE = "epithelial"
INPUT_H5AD = WORKFLOW_ROOT / "h5ad" / MODULE / "01-lineage-selection" / LINEAGE / "adata_epithelial_qc.h5ad"
OUTPUT_H5AD = WORKFLOW_ROOT / "h5ad" / MODULE / STEP / LINEAGE / "adata_epithelial_harmony.h5ad"
TABLE_DIR = WORKFLOW_ROOT / "tables" / MODULE / STEP / LINEAGE
CODE_FILE = Path(__file__)

TARGET_SUM = 10000.0
N_TOP_GENES = 3000
HVG_FLAVOR = "seurat"
REGRESS_KEYS = ["total_counts", "pct_counts_MT"]
SCALE_MAX_VALUE = 10
N_COMPS = 50
BATCH_KEY = "sample"
PCA_BASIS = "X_pca"
HARMONY_BASIS = "X_pca_inte"


def assert_no_overwrite() -> None:
    outputs = [
        OUTPUT_H5AD,
        TABLE_DIR / "01_normalize_total_parameters.csv",
        TABLE_DIR / "02_log1p_parameters.csv",
        TABLE_DIR / "03_highly_variable_genes_parameters.csv",
        TABLE_DIR / "04_raw_assignment_and_hvg_subset_record.csv",
        TABLE_DIR / "05_regress_out_parameters.csv",
        TABLE_DIR / "06_scale_parameters.csv",
        TABLE_DIR / "07_pca_parameters.csv",
        TABLE_DIR / "08_harmony_integrate_parameters.csv",
        TABLE_DIR / "subtype_harmony_summary.csv",
        TABLE_DIR / "package_versions.txt",
        TABLE_DIR / "readme.txt",
    ]
    # Allow retry after a failed run that wrote early parameter tables but did not
    # produce the final Harmony h5ad. The stable final object remains protected.
    retry_partial = not OUTPUT_H5AD.exists()
    partial_retry_allowed = {
        TABLE_DIR / "01_normalize_total_parameters.csv",
        TABLE_DIR / "02_log1p_parameters.csv",
        TABLE_DIR / "03_highly_variable_genes_parameters.csv",
        TABLE_DIR / "04_raw_assignment_and_hvg_subset_record.csv",
        TABLE_DIR / "05_regress_out_parameters.csv",
        TABLE_DIR / "06_scale_parameters.csv",
        TABLE_DIR / "07_pca_parameters.csv",
    }
    existing = [
        str(path)
        for path in outputs
        if path.exists() and not (retry_partial and path in partial_retry_allowed)
    ]
    if existing:
        raise FileExistsError("Subtype Harmony output already exists; refusing to overwrite:\n" + "\n".join(existing))


def one_row(path: Path, row: dict) -> None:
    pd.DataFrame([row]).to_csv(path, index=False)


def main() -> None:
    assert_no_overwrite()
    OUTPUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(INPUT_H5AD)
    n_obs0, n_vars0 = adata.n_obs, adata.n_vars
    if BATCH_KEY not in adata.obs.columns:
        raise KeyError(f"Missing obs[{BATCH_KEY!r}] for Harmony batch correction.")
    for key in REGRESS_KEYS:
        if key not in adata.obs.columns:
            raise KeyError(f"Missing obs[{key!r}] for regress_out.")

    sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
    one_row(
        TABLE_DIR / "01_normalize_total_parameters.csv",
        {
            "step": "normalize_total",
            "input_h5ad": str(INPUT_H5AD),
            "output_h5ad_or_object": "in_memory_adata",
            "n_obs_before": n_obs0,
            "n_obs_after": adata.n_obs,
            "n_vars_before": n_vars0,
            "n_vars_after": adata.n_vars,
            "target_sum": TARGET_SUM,
            "backend_package": "scanpy_cpu",
            "code_file": str(CODE_FILE),
            "random_seed": SEED,
        },
    )

    sc.pp.log1p(adata)
    one_row(
        TABLE_DIR / "02_log1p_parameters.csv",
        {
            "step": "log1p",
            "input_h5ad": str(INPUT_H5AD),
            "output_h5ad_or_object": "in_memory_adata",
            "n_obs_before": n_obs0,
            "n_obs_after": adata.n_obs,
            "n_vars_before": n_vars0,
            "n_vars_after": adata.n_vars,
            "backend_package": "scanpy_cpu",
            "code_file": str(CODE_FILE),
            "random_seed": SEED,
        },
    )

    sc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_GENES, flavor=HVG_FLAVOR)
    actual_hvg = int(adata.var["highly_variable"].sum())
    one_row(
        TABLE_DIR / "03_highly_variable_genes_parameters.csv",
        {
            "step": "highly_variable_genes",
            "input_h5ad": str(INPUT_H5AD),
            "output_h5ad_or_object": "in_memory_adata_with_hvg_flags",
            "n_obs_before": n_obs0,
            "n_obs_after": adata.n_obs,
            "n_vars_before": n_vars0,
            "n_vars_after": adata.n_vars,
            "n_top_genes": N_TOP_GENES,
            "actual_hvg_count": actual_hvg,
            "flavor": HVG_FLAVOR,
            "backend_package": "scanpy_cpu",
            "code_file": str(CODE_FILE),
            "random_seed": SEED,
        },
    )

    adata.raw = adata
    adata = adata[:, adata.var["highly_variable"].to_numpy()].copy()
    one_row(
        TABLE_DIR / "04_raw_assignment_and_hvg_subset_record.csv",
        {
            "step": "raw_assignment_and_hvg_subset",
            "input_h5ad": str(INPUT_H5AD),
            "output_h5ad_or_object": "in_memory_hvg_subset",
            "n_obs_before": n_obs0,
            "n_obs_after": adata.n_obs,
            "n_vars_before": n_vars0,
            "n_vars_after": adata.n_vars,
            "raw_assigned": True,
            "raw_n_vars": int(adata.raw.n_vars),
            "hvg_subset": True,
            "code_file": str(CODE_FILE),
            "random_seed": SEED,
        },
    )

    sc.pp.regress_out(adata, keys=REGRESS_KEYS)
    one_row(
        TABLE_DIR / "05_regress_out_parameters.csv",
        {
            "step": "regress_out",
            "input_h5ad": str(INPUT_H5AD),
            "output_h5ad_or_object": "in_memory_hvg_subset",
            "n_obs_before": n_obs0,
            "n_obs_after": adata.n_obs,
            "n_vars_before": N_TOP_GENES,
            "n_vars_after": adata.n_vars,
            "regress_keys": ";".join(REGRESS_KEYS),
            "backend_package": "scanpy_cpu",
            "code_file": str(CODE_FILE),
            "random_seed": SEED,
        },
    )

    sc.pp.scale(adata, max_value=SCALE_MAX_VALUE)
    one_row(
        TABLE_DIR / "06_scale_parameters.csv",
        {
            "step": "scale",
            "input_h5ad": str(INPUT_H5AD),
            "output_h5ad_or_object": "in_memory_hvg_subset",
            "n_obs_before": n_obs0,
            "n_obs_after": adata.n_obs,
            "n_vars_before": N_TOP_GENES,
            "n_vars_after": adata.n_vars,
            "max_value": SCALE_MAX_VALUE,
            "backend_package": "scanpy_cpu",
            "code_file": str(CODE_FILE),
            "random_seed": SEED,
        },
    )

    sc.tl.pca(adata, n_comps=N_COMPS, random_state=SEED)
    one_row(
        TABLE_DIR / "07_pca_parameters.csv",
        {
            "step": "pca",
            "input_h5ad": str(INPUT_H5AD),
            "output_h5ad_or_object": f"obsm_{PCA_BASIS}",
            "n_obs_before": n_obs0,
            "n_obs_after": adata.n_obs,
            "n_vars_before": N_TOP_GENES,
            "n_vars_after": adata.n_vars,
            "n_comps": N_COMPS,
            "basis": PCA_BASIS,
            "backend_package": "scanpy_cpu",
            "code_file": str(CODE_FILE),
            "random_seed": SEED,
        },
    )

    import harmonypy as hm

    harmony_out = hm.run_harmony(
        adata.obsm[PCA_BASIS],
        adata.obs,
        BATCH_KEY,
        max_iter_harmony=10,
        random_state=SEED,
    )
    corrected = np.asarray(harmony_out.Z_corr)
    if corrected.shape == (N_COMPS, adata.n_obs):
        corrected = corrected.T
    if corrected.shape != (adata.n_obs, N_COMPS):
        raise RuntimeError(
            f"Harmony corrected matrix has shape {corrected.shape}; expected {(adata.n_obs, N_COMPS)}."
        )
    adata.obsm[HARMONY_BASIS] = corrected.astype(np.float32)
    one_row(
        TABLE_DIR / "08_harmony_integrate_parameters.csv",
        {
            "step": "harmony_integrate",
            "input_h5ad": str(INPUT_H5AD),
            "output_h5ad_or_object": str(OUTPUT_H5AD),
            "n_obs_before": n_obs0,
            "n_obs_after": adata.n_obs,
            "n_vars_before": N_TOP_GENES,
            "n_vars_after": adata.n_vars,
            "key": BATCH_KEY,
            "basis": PCA_BASIS,
            "adjusted_basis": HARMONY_BASIS,
            "backend_package": "harmonypy_cpu",
            "code_file": str(CODE_FILE),
            "random_seed": SEED,
        },
    )

    adata.write_h5ad(OUTPUT_H5AD)
    summary = pd.DataFrame(
        [
            {
                "step": "epithelial_subtype_harmony",
                "lineage": LINEAGE,
                "input_h5ad": str(INPUT_H5AD),
                "output_h5ad": str(OUTPUT_H5AD),
                "n_obs": int(adata.n_obs),
                "n_vars": int(adata.n_vars),
                "raw_n_vars": int(adata.raw.n_vars),
                "hvg_count": actual_hvg,
                "batch_key": BATCH_KEY,
                "pca_basis": PCA_BASIS,
                "harmony_basis": HARMONY_BASIS,
                "code_file": str(CODE_FILE),
                "random_seed": SEED,
            }
        ]
    )
    summary.to_csv(TABLE_DIR / "subtype_harmony_summary.csv", index=False)

    with (TABLE_DIR / "package_versions.txt").open("w") as fh:
        fh.write(f"python: {platform.python_version()}\n")
        fh.write(f"anndata: {ad.__version__}\n")
        fh.write(f"scanpy: {sc.__version__}\n")
        fh.write(f"numpy: {np.__version__}\n")
        fh.write(f"pandas: {pd.__version__}\n")
        try:
            import harmonypy

            fh.write(f"harmonypy: {harmonypy.__version__}\n")
        except Exception:
            fh.write("harmonypy: unavailable_version\n")
        fh.write(f"code_file: {CODE_FILE}\n")

    with (TABLE_DIR / "readme.txt").open("w") as fh:
        fh.write("Epithelial subtype Harmony integration completed.\n")
        fh.write(f"Input lineage h5ad: {INPUT_H5AD}\n")
        fh.write(f"Output Harmony h5ad: {OUTPUT_H5AD}\n")
        fh.write("Parameters matched Module 01 defaults: normalize_total/log1p/HVG/raw/HVG subset/regress/scale/PCA/Harmony.\n")
        fh.write(f"Harmony batch key: {BATCH_KEY}\n")
        fh.write(f"Harmony basis: obsm['{HARMONY_BASIS}']\n")

    print(summary.to_string(index=False))
    del adata
    gc.collect()


if __name__ == "__main__":
    main()
