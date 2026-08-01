#!/usr/bin/env python3
"""Normalize, select HVGs, regress, PCA, and Harmony-integrate full BRCA data."""

from __future__ import annotations

import gc
import importlib.metadata
import json
import random
import sys
import time
import warnings
from pathlib import Path

import cupy as cp
import cupyx.scipy.sparse as cpsparse
import numpy as np
import pandas as pd
import rapids_singlecell as rsc
import scanpy as sc


SEED = 42
random.seed(SEED)
np.random.seed(SEED)
cp.random.seed(SEED)

ANALYSIS_ROOT = Path("/mnt/disk18t/lr_xcy/riku/brca_cm-epi_codex")
WORKFLOW_ROOT = ANALYSIS_ROOT / "epi-cm-core-workflow"
BLOCK = "01-celltype_integration_clustering"
INPUT_H5AD = WORKFLOW_ROOT / "h5ad" / BLOCK / "03-qc" / "adata_qc.h5ad"
OUT_H5AD = (
    WORKFLOW_ROOT / "h5ad" / BLOCK / "04-integration-harmony" / "adata_harmony.h5ad"
)
TABLE_DIR = WORKFLOW_ROOT / "tables" / BLOCK / "04-integration-harmony"
CODE_PATH = (
    WORKFLOW_ROOT
    / "codes"
    / BLOCK
    / "04-integration-harmony"
    / "03_integrate_harmony.py"
)

TARGET_SUM = 1e4
N_TOP_GENES = 3000
HVG_FLAVOR = "seurat"
HVG_BATCH_KEY = "sample"
REGRESSION_KEYS = ["total_counts", "pct_counts_MT"]
REGRESSION_BATCHSIZE = 10000
SCALE_MAX_VALUE = 10.0
PCA_N_COMPS = 50
BATCH_KEY = "sample"
PCA_BASIS = "X_pca"
INTEGRATED_BASIS = "X_pca_inte"
HARMONY_INITIAL_MAX_ITER = 10
HARMONY_RETRY_MAX_ITER = 20


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def release_gpu() -> None:
    cp.get_default_memory_pool().free_all_blocks()
    cp.get_default_pinned_memory_pool().free_all_blocks()
    gc.collect()


def cpuify(value):
    if isinstance(value, cp.ndarray):
        return cp.asnumpy(value)
    if cpsparse.isspmatrix(value):
        return value.get()
    if isinstance(value, dict):
        return {k: cpuify(v) for k, v in value.items()}
    if isinstance(value, list):
        return [cpuify(v) for v in value]
    if isinstance(value, tuple):
        return tuple(cpuify(v) for v in value)
    return value


def convert_adata_aux_to_cpu(adata) -> None:
    for key in list(adata.obsm.keys()):
        adata.obsm[key] = cpuify(adata.obsm[key])
    for key in list(adata.varm.keys()):
        adata.varm[key] = cpuify(adata.varm[key])
    for key in list(adata.obsp.keys()):
        adata.obsp[key] = cpuify(adata.obsp[key])
    for key in list(adata.varp.keys()):
        adata.varp[key] = cpuify(adata.varp[key])
    for key in list(adata.layers.keys()):
        adata.layers[key] = cpuify(adata.layers[key])
    for key in list(adata.uns.keys()):
        adata.uns[key] = cpuify(adata.uns[key])


def main() -> None:
    start = time.time()
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    OUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(INPUT_H5AD)
    required_obs = {"sample", "series", "status", "original_barcode", "total_counts", "pct_counts_MT"}
    missing = required_obs - set(adata.obs.columns)
    if missing:
        raise ValueError(f"Missing required obs columns: {sorted(missing)}")
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError("QC input is empty.")
    if not adata.obs_names.is_unique:
        raise ValueError("QC input obs_names are not unique.")

    source_count_dtype = str(adata.X.dtype)
    # CuPy 14 sparse matrices do not accept integer dtypes.  The authoritative
    # integer counts remain in INPUT_H5AD; normalization starts from an exact
    # float32 representation in this derived integration object.
    if adata.X.dtype not in (np.float32, np.float64):
        adata.X = adata.X.astype(np.float32)

    backend_rows: list[dict[str, object]] = []
    print(
        f"[Harmony] input shape={adata.shape} source_dtype={source_count_dtype} "
        f"working_dtype={adata.X.dtype}",
        flush=True,
    )
    t0 = time.time()
    rsc.get.anndata_to_GPU(adata)
    backend_rows.append(
        {
            "step": "anndata_to_GPU",
            "planned_backend": "rapids_gpu",
            "attempted_backend": "rapids_gpu",
            "status": "completed",
            "error_summary": "",
            "fallback_backend": "",
            "clean_input_reloaded": True,
            "final_backend_for_rerun": "rapids_gpu",
            "elapsed_seconds": time.time() - t0,
        }
    )

    step_times: dict[str, float] = {}
    t0 = time.time()
    rsc.pp.normalize_total(adata, target_sum=int(TARGET_SUM))
    rsc.pp.log1p(adata)
    step_times["normalize_log1p"] = time.time() - t0
    print(
        f"[Harmony] normalize_log1p completed in "
        f"{step_times['normalize_log1p']:.1f}s",
        flush=True,
    )

    t0 = time.time()
    rsc.pp.highly_variable_genes(
        adata,
        n_top_genes=N_TOP_GENES,
        flavor=HVG_FLAVOR,
        batch_key=HVG_BATCH_KEY,
    )
    step_times["highly_variable_genes"] = time.time() - t0
    n_hvg = int(adata.var["highly_variable"].sum())
    print(
        f"[Harmony] highly_variable_genes completed: n_hvg={n_hvg} "
        f"elapsed={step_times['highly_variable_genes']:.1f}s",
        flush=True,
    )
    if n_hvg < 2000:
        raise ValueError(f"Only {n_hvg} HVGs selected; at least 2000 are required.")

    rsc.get.anndata_to_CPU(adata)
    convert_adata_aux_to_cpu(adata)
    adata.raw = adata.copy()
    hvg_mask = adata.var["highly_variable"].to_numpy(dtype=bool)
    adata = adata[:, hvg_mask].copy()
    print(
        f"[Harmony] raw stored with {adata.raw.n_vars} genes; "
        f"working HVG matrix={adata.shape}",
        flush=True,
    )

    rsc.get.anndata_to_GPU(adata)
    t0 = time.time()
    rsc.pp.regress_out(
        adata,
        keys=REGRESSION_KEYS,
        batchsize=REGRESSION_BATCHSIZE,
        verbose=True,
    )
    step_times["regress_out"] = time.time() - t0
    print(
        f"[Harmony] regress_out completed in {step_times['regress_out']:.1f}s",
        flush=True,
    )
    t0 = time.time()
    rsc.pp.scale(adata, max_value=SCALE_MAX_VALUE)
    step_times["scale"] = time.time() - t0
    print(f"[Harmony] scale completed in {step_times['scale']:.1f}s", flush=True)
    t0 = time.time()
    rsc.tl.pca(
        adata,
        n_comps=PCA_N_COMPS,
        random_state=SEED,
        dtype="float32",
    )
    step_times["pca"] = time.time() - t0
    print(
        f"[Harmony] PCA completed: shape={adata.obsm[PCA_BASIS].shape} "
        f"elapsed={step_times['pca']:.1f}s",
        flush=True,
    )

    harmony_warnings: list[str] = []
    harmony_max_iter_used = HARMONY_INITIAL_MAX_ITER
    t0 = time.time()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        rsc.pp.harmony_integrate(
            adata,
            key=BATCH_KEY,
            basis=PCA_BASIS,
            adjusted_basis=INTEGRATED_BASIS,
            max_iter_harmony=HARMONY_INITIAL_MAX_ITER,
            random_state=SEED,
            verbose=True,
        )
    harmony_warnings.extend(str(w.message) for w in caught)
    did_not_converge = any("did not converge" in text.lower() for text in harmony_warnings)
    if did_not_converge:
        harmony_max_iter_used = HARMONY_RETRY_MAX_ITER
        with warnings.catch_warnings(record=True) as caught_retry:
            warnings.simplefilter("always")
            rsc.pp.harmony_integrate(
                adata,
                key=BATCH_KEY,
                basis=PCA_BASIS,
                adjusted_basis=INTEGRATED_BASIS,
                max_iter_harmony=HARMONY_RETRY_MAX_ITER,
                random_state=SEED,
                verbose=True,
            )
        retry_messages = [str(w.message) for w in caught_retry]
        harmony_warnings.extend(f"retry: {text}" for text in retry_messages)
        if any("did not converge" in text.lower() for text in retry_messages):
            raise RuntimeError(
                "Harmony did not converge after the documented increased max_iter_harmony retry."
            )
    step_times["harmony"] = time.time() - t0
    print(
        f"[Harmony] integration completed with max_iter={harmony_max_iter_used} "
        f"in {step_times['harmony']:.1f}s",
        flush=True,
    )

    rsc.get.anndata_to_CPU(adata)
    convert_adata_aux_to_cpu(adata)
    release_gpu()

    if PCA_BASIS not in adata.obsm or INTEGRATED_BASIS not in adata.obsm:
        raise ValueError("PCA or Harmony basis missing after integration.")
    if adata.obsm[PCA_BASIS].shape != (adata.n_obs, PCA_N_COMPS):
        raise ValueError(f"Unexpected PCA shape: {adata.obsm[PCA_BASIS].shape}")
    if adata.obsm[INTEGRATED_BASIS].shape != (adata.n_obs, PCA_N_COMPS):
        raise ValueError(f"Unexpected Harmony shape: {adata.obsm[INTEGRATED_BASIS].shape}")
    if not np.isfinite(adata.obsm[INTEGRATED_BASIS]).all():
        raise ValueError("Harmony embedding contains non-finite values.")
    forbidden = {"neighbors", "umap", "leiden"}
    if forbidden & set(adata.uns):
        raise ValueError("Harmony handoff unexpectedly contains graph/clustering outputs.")
    if "X_umap" in adata.obsm:
        raise ValueError("Harmony handoff unexpectedly contains X_umap.")

    adata.uns["integration_parameters"] = {
        "target_sum": TARGET_SUM,
        "log1p": True,
        "n_top_genes": N_TOP_GENES,
        "hvg_flavor": HVG_FLAVOR,
        "hvg_batch_key": HVG_BATCH_KEY,
        "raw_assignment_point": "after normalize_total/log1p/HVG marking before HVG subsetting",
        "regression_keys": REGRESSION_KEYS,
        "regression_batchsize": REGRESSION_BATCHSIZE,
        "scale_max_value": SCALE_MAX_VALUE,
        "pca_n_comps": PCA_N_COMPS,
        "batch_key": BATCH_KEY,
        "pca_basis": PCA_BASIS,
        "integrated_basis": INTEGRATED_BASIS,
        "harmony_max_iter_used": harmony_max_iter_used,
        "source_count_dtype": source_count_dtype,
        "working_dtype_after_count_conversion": "float32",
        "seed": SEED,
    }
    adata.uns["expression_contract"] = {
        "X": "HVG-subset regressed and scaled expression",
        "raw": "full-gene normalized/log1p expression",
        "raw_counts": str(INPUT_H5AD) + "::.X",
    }

    pd.DataFrame(
        [
            {
                "input_h5ad": str(INPUT_H5AD),
                "output_h5ad": str(OUT_H5AD),
                "n_cells": adata.n_obs,
                "n_genes_input": adata.raw.n_vars,
                "n_hvg": n_hvg,
                "target_sum": TARGET_SUM,
                "hvg_flavor": HVG_FLAVOR,
                "hvg_batch_key": HVG_BATCH_KEY,
                "regression_keys": ";".join(REGRESSION_KEYS),
                "regression_batchsize": REGRESSION_BATCHSIZE,
                "scale_max_value": SCALE_MAX_VALUE,
                "pca_n_comps": PCA_N_COMPS,
                "batch_key": BATCH_KEY,
                "pca_basis": PCA_BASIS,
                "integrated_basis": INTEGRATED_BASIS,
                "harmony_initial_max_iter": HARMONY_INITIAL_MAX_ITER,
                "harmony_final_max_iter": harmony_max_iter_used,
                "seed": SEED,
                "backend": "rapids_singlecell_gpu",
            }
        ]
    ).to_csv(TABLE_DIR / "harmony_parameters.csv", index=False)
    pd.DataFrame(
        [{"step": key, "elapsed_seconds": value} for key, value in step_times.items()]
    ).to_csv(TABLE_DIR / "integration_step_times.csv", index=False)
    (TABLE_DIR / "harmony_warnings.txt").write_text(
        "\n".join(harmony_warnings) + ("\n" if harmony_warnings else ""),
        encoding="utf-8",
    )
    backend_rows.extend(
        {
            "step": step,
            "planned_backend": "rapids_gpu",
            "attempted_backend": "rapids_gpu",
            "status": "completed",
            "error_summary": "",
            "fallback_backend": "",
            "clean_input_reloaded": True,
            "final_backend_for_rerun": "rapids_gpu",
            "elapsed_seconds": elapsed,
        }
        for step, elapsed in step_times.items()
    )
    pd.DataFrame(backend_rows).to_csv(
        TABLE_DIR / "gpu_backend_capability_summary.csv", index=False
    )

    ens_re = r"^ENS[A-Z]*G[0-9]+(?:\.[0-9]+)?$"
    gene_rows = []
    for object_layer, names in [
        ("var_names", pd.Index(adata.var_names.astype(str))),
        ("raw.var_names", pd.Index(adata.raw.var_names.astype(str))),
    ]:
        ensembl = names.to_series(index=range(len(names))).str.match(ens_re)
        fraction = float(ensembl.mean()) if len(names) else 0.0
        gene_rows.append(
            {
                "input_h5ad": str(OUT_H5AD),
                "object_layer": object_layer,
                "var_names_type": "gene_id" if fraction >= 0.9 else "gene_name" if fraction <= 0.1 else "mixed_or_unknown",
                "n_features": len(names),
                "n_unique": int(names.nunique()),
                "n_duplicates": int(names.duplicated().sum()),
                "ensembl_like_n": int(ensembl.sum()),
                "ensembl_like_fraction": fraction,
                "candidate_gene_id_columns": "",
                "candidate_gene_name_columns": "gene_name",
                "selected_gene_id_source": "",
                "selected_gene_name_source": "var_names",
                "action_taken": "preserved_unique_gene_symbols",
            }
        )
    pd.DataFrame(gene_rows).to_csv(TABLE_DIR / "gene_identifier_audit.csv", index=False)

    versions = {
        "python": sys.version.split()[0],
        "scanpy": package_version("scanpy"),
        "anndata": package_version("anndata"),
        "rapids-singlecell": package_version("rapids-singlecell"),
        "cuml-cu12": package_version("cuml-cu12"),
        "cugraph-cu12": package_version("cugraph-cu12"),
        "cupy-cuda12x": package_version("cupy-cuda12x"),
        "code": str(CODE_PATH),
        "seed": str(SEED),
    }
    (TABLE_DIR / "package_versions.txt").write_text(
        "\n".join(f"{k}={v}" for k, v in versions.items()) + "\n",
        encoding="utf-8",
    )
    (TABLE_DIR / "readme.txt").write_text(
        f"""BRCA post-QC Harmony integration

Input: {INPUT_H5AD}
Code: {CODE_PATH}
Output: {OUT_H5AD}

Full cells are normalized to {TARGET_SUM:g}, log1p transformed, marked for
{N_TOP_GENES} batch-aware {HVG_FLAVOR} HVGs, and saved to adata.raw before HVG
subsetting. total_counts and pct_counts_MT are regressed, values are scaled with
max_value={SCALE_MAX_VALUE:g}, PCA uses {PCA_N_COMPS} components, and Harmony
uses sample as the batch key. No neighbor graph, UMAP, or clustering is stored
in this handoff object. No cells are downsampled.
""",
        encoding="utf-8",
    )

    adata.write_h5ad(OUT_H5AD, compression="gzip")
    report = {
        "n_cells": int(adata.n_obs),
        "n_hvg": int(adata.n_vars),
        "n_raw_genes": int(adata.raw.n_vars),
        "pca_shape": list(adata.obsm[PCA_BASIS].shape),
        "harmony_shape": list(adata.obsm[INTEGRATED_BASIS].shape),
        "harmony_max_iter_used": harmony_max_iter_used,
        "elapsed_seconds": time.time() - start,
    }
    (TABLE_DIR / "harmony_completion.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
