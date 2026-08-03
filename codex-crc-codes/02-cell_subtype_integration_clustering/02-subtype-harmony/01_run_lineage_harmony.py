#!/usr/bin/env python3
"""Run one lineage-specific subtype Harmony template."""

from __future__ import annotations

import argparse
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
MAX_ITER_HARMONY = 10


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lineage", required=True)
    parser.add_argument("--allow-existing", action="store_true")
    return parser.parse_args()


def paths(lineage: str) -> tuple[Path, Path, Path]:
    input_h5ad = WORKFLOW_ROOT / "h5ad" / MODULE / "01-lineage-selection" / lineage / f"adata_{lineage}_qc.h5ad"
    output_h5ad = WORKFLOW_ROOT / "h5ad" / MODULE / STEP / lineage / f"adata_{lineage}_harmony.h5ad"
    table_dir = WORKFLOW_ROOT / "tables" / MODULE / STEP / lineage
    return input_h5ad, output_h5ad, table_dir


def one_row(path: Path, row: dict) -> None:
    pd.DataFrame([row]).to_csv(path, index=False)


def validate_existing(output_h5ad: Path, table_dir: Path, lineage: str) -> None:
    adata = ad.read_h5ad(output_h5ad, backed="r")
    ok = adata.raw is not None and PCA_BASIS in adata.obsm and HARMONY_BASIS in adata.obsm
    n_obs, n_vars = adata.shape
    adata.file.close()
    if not ok:
        raise RuntimeError(f"Existing Harmony h5ad failed validation: {output_h5ad}")
    table_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [
            {
                "lineage": lineage,
                "output_h5ad": str(output_h5ad),
                "n_obs": int(n_obs),
                "n_vars": int(n_vars),
                "raw_present": True,
                "pca_basis": PCA_BASIS,
                "harmony_basis": HARMONY_BASIS,
                "status": "completed_existing",
            }
        ]
    ).to_csv(table_dir / "subtype_harmony_existing_validation.csv", index=False)


def assert_no_overwrite(output_h5ad: Path, table_dir: Path, allow_existing: bool) -> None:
    outputs = [
        output_h5ad,
        table_dir / "01_normalize_total_parameters.csv",
        table_dir / "02_log1p_parameters.csv",
        table_dir / "03_highly_variable_genes_parameters.csv",
        table_dir / "04_raw_assignment_and_hvg_subset_record.csv",
        table_dir / "05_regress_out_parameters.csv",
        table_dir / "06_scale_parameters.csv",
        table_dir / "07_pca_parameters.csv",
        table_dir / "08_harmony_integrate_parameters.csv",
        table_dir / "subtype_harmony_summary.csv",
        table_dir / "package_versions.txt",
        table_dir / "readme.txt",
    ]
    existing = [str(path) for path in outputs if path.exists()]
    if existing and not allow_existing:
        raise FileExistsError("Subtype Harmony output already exists; refusing to overwrite:\n" + "\n".join(existing))


def main() -> None:
    args = parse_args()
    input_h5ad, output_h5ad, table_dir = paths(args.lineage)
    if args.allow_existing and output_h5ad.exists():
        validate_existing(output_h5ad, table_dir, args.lineage)
        return
    assert_no_overwrite(output_h5ad, table_dir, args.allow_existing)

    output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(input_h5ad)
    n_obs0, n_vars0 = adata.n_obs, adata.n_vars
    if BATCH_KEY not in adata.obs.columns:
        raise KeyError(f"Missing obs[{BATCH_KEY!r}] for Harmony batch correction.")
    for key in REGRESS_KEYS:
        if key not in adata.obs.columns:
            raise KeyError(f"Missing obs[{key!r}] for regress_out.")

    sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
    one_row(table_dir / "01_normalize_total_parameters.csv", {
        "step": "normalize_total", "input_h5ad": str(input_h5ad), "output_h5ad_or_object": "in_memory_adata",
        "n_obs_before": n_obs0, "n_obs_after": adata.n_obs, "n_vars_before": n_vars0, "n_vars_after": adata.n_vars,
        "target_sum": TARGET_SUM, "backend_package": "scanpy_cpu", "code_file": str(CODE_FILE), "random_seed": SEED,
    })

    sc.pp.log1p(adata)
    one_row(table_dir / "02_log1p_parameters.csv", {
        "step": "log1p", "input_h5ad": str(input_h5ad), "output_h5ad_or_object": "in_memory_adata",
        "n_obs_before": n_obs0, "n_obs_after": adata.n_obs, "n_vars_before": n_vars0, "n_vars_after": adata.n_vars,
        "backend_package": "scanpy_cpu", "code_file": str(CODE_FILE), "random_seed": SEED,
    })

    sc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_GENES, flavor=HVG_FLAVOR)
    actual_hvg = int(adata.var["highly_variable"].sum())
    one_row(table_dir / "03_highly_variable_genes_parameters.csv", {
        "step": "highly_variable_genes", "input_h5ad": str(input_h5ad),
        "output_h5ad_or_object": "in_memory_adata_with_hvg_flags", "n_obs_before": n_obs0, "n_obs_after": adata.n_obs,
        "n_vars_before": n_vars0, "n_vars_after": adata.n_vars, "n_top_genes": N_TOP_GENES,
        "actual_hvg_count": actual_hvg, "flavor": HVG_FLAVOR, "backend_package": "scanpy_cpu",
        "code_file": str(CODE_FILE), "random_seed": SEED,
    })

    adata.raw = adata
    adata = adata[:, adata.var["highly_variable"].to_numpy()].copy()
    one_row(table_dir / "04_raw_assignment_and_hvg_subset_record.csv", {
        "step": "raw_assignment_and_hvg_subset", "input_h5ad": str(input_h5ad), "output_h5ad_or_object": "in_memory_hvg_subset",
        "n_obs_before": n_obs0, "n_obs_after": adata.n_obs, "n_vars_before": n_vars0, "n_vars_after": adata.n_vars,
        "raw_assigned": True, "raw_n_vars": int(adata.raw.n_vars), "hvg_subset": True,
        "code_file": str(CODE_FILE), "random_seed": SEED,
    })

    sc.pp.regress_out(adata, keys=REGRESS_KEYS)
    one_row(table_dir / "05_regress_out_parameters.csv", {
        "step": "regress_out", "input_h5ad": str(input_h5ad), "output_h5ad_or_object": "in_memory_hvg_subset",
        "n_obs_before": n_obs0, "n_obs_after": adata.n_obs, "n_vars_before": N_TOP_GENES, "n_vars_after": adata.n_vars,
        "regress_keys": ";".join(REGRESS_KEYS), "backend_package": "scanpy_cpu", "code_file": str(CODE_FILE),
        "random_seed": SEED,
    })

    sc.pp.scale(adata, max_value=SCALE_MAX_VALUE)
    one_row(table_dir / "06_scale_parameters.csv", {
        "step": "scale", "input_h5ad": str(input_h5ad), "output_h5ad_or_object": "in_memory_hvg_subset",
        "n_obs_before": n_obs0, "n_obs_after": adata.n_obs, "n_vars_before": N_TOP_GENES, "n_vars_after": adata.n_vars,
        "max_value": SCALE_MAX_VALUE, "backend_package": "scanpy_cpu", "code_file": str(CODE_FILE), "random_seed": SEED,
    })

    sc.tl.pca(adata, n_comps=N_COMPS, random_state=SEED)
    one_row(table_dir / "07_pca_parameters.csv", {
        "step": "pca", "input_h5ad": str(input_h5ad), "output_h5ad_or_object": f"obsm_{PCA_BASIS}",
        "n_obs_before": n_obs0, "n_obs_after": adata.n_obs, "n_vars_before": N_TOP_GENES, "n_vars_after": adata.n_vars,
        "n_comps": N_COMPS, "basis": PCA_BASIS, "backend_package": "scanpy_cpu",
        "code_file": str(CODE_FILE), "random_seed": SEED,
    })

    import harmonypy as hm

    harmony_out = hm.run_harmony(
        adata.obsm[PCA_BASIS],
        adata.obs,
        BATCH_KEY,
        max_iter_harmony=MAX_ITER_HARMONY,
        random_state=SEED,
    )
    corrected = np.asarray(harmony_out.Z_corr)
    if corrected.shape == (N_COMPS, adata.n_obs):
        corrected = corrected.T
    if corrected.shape != (adata.n_obs, N_COMPS):
        raise RuntimeError(f"Harmony corrected matrix has shape {corrected.shape}; expected {(adata.n_obs, N_COMPS)}.")
    adata.obsm[HARMONY_BASIS] = corrected.astype(np.float32)
    one_row(table_dir / "08_harmony_integrate_parameters.csv", {
        "step": "harmony_integrate", "input_h5ad": str(input_h5ad), "output_h5ad_or_object": str(output_h5ad),
        "n_obs_before": n_obs0, "n_obs_after": adata.n_obs, "n_vars_before": N_TOP_GENES, "n_vars_after": adata.n_vars,
        "key": BATCH_KEY, "basis": PCA_BASIS, "adjusted_basis": HARMONY_BASIS, "max_iter_harmony": MAX_ITER_HARMONY,
        "backend_package": "harmonypy_cpu", "code_file": str(CODE_FILE), "random_seed": SEED,
    })

    adata.write_h5ad(output_h5ad)
    summary = pd.DataFrame([{
        "step": "subtype_harmony", "lineage": args.lineage, "input_h5ad": str(input_h5ad), "output_h5ad": str(output_h5ad),
        "n_obs": int(adata.n_obs), "n_vars": int(adata.n_vars), "raw_n_vars": int(adata.raw.n_vars),
        "hvg_count": actual_hvg, "batch_key": BATCH_KEY, "pca_basis": PCA_BASIS, "harmony_basis": HARMONY_BASIS,
        "code_file": str(CODE_FILE), "random_seed": SEED,
    }])
    summary.to_csv(table_dir / "subtype_harmony_summary.csv", index=False)

    with (table_dir / "package_versions.txt").open("w") as fh:
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

    with (table_dir / "readme.txt").open("w") as fh:
        fh.write(f"{args.lineage} subtype Harmony integration completed.\n")
        fh.write(f"Input lineage h5ad: {input_h5ad}\n")
        fh.write(f"Output Harmony h5ad: {output_h5ad}\n")
        fh.write("Parameters matched Module 01 defaults: normalize_total/log1p/HVG/raw/HVG subset/regress/scale/PCA/Harmony.\n")
        fh.write(f"Harmony batch key: {BATCH_KEY}\n")
        fh.write(f"Harmony basis: obsm['{HARMONY_BASIS}']\n")

    print(summary.to_string(index=False))
    del adata
    gc.collect()


if __name__ == "__main__":
    main()
