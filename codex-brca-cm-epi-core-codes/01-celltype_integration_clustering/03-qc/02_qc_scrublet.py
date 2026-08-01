#!/usr/bin/env python3
"""Apply sample-wise Scrublet and preserve the published BRCA QC contract."""

from __future__ import annotations

import gc
import importlib.metadata
import json
import random
import sys
import time
from pathlib import Path

import anndata as ad
import cupy as cp
import numpy as np
import pandas as pd
import rapids_singlecell as rsc
import scanpy as sc
from scipy import sparse


SEED = 42
random.seed(SEED)
np.random.seed(SEED)
cp.random.seed(SEED)

ANALYSIS_ROOT = Path("/mnt/disk18t/lr_xcy/riku/brca_cm-epi_codex")
WORKFLOW_ROOT = ANALYSIS_ROOT / "epi-cm-core-workflow"
BLOCK = "01-celltype_integration_clustering"
INPUT_H5AD = (
    WORKFLOW_ROOT / "h5ad" / BLOCK / "02-merge-metadata" / "adata_merge.h5ad"
)
OUT_H5AD = WORKFLOW_ROOT / "h5ad" / BLOCK / "03-qc" / "adata_qc.h5ad"
TABLE_DIR = WORKFLOW_ROOT / "tables" / BLOCK / "03-qc"
FIGURE_DIR = WORKFLOW_ROOT / "figures" / BLOCK / "03-qc"
CODE_PATH = WORKFLOW_ROOT / "codes" / BLOCK / "03-qc" / "02_qc_scrublet.py"

MIN_GENES_INITIAL = 200
MIN_CELLS_GENE = 3
MIN_COUNTS_GENE = 3
MIN_GENES_FINAL = 200
MAX_GENES_FINAL = 11000
MAX_PCT_MT = 20.0
MT_PREFIX = "MT-"
RIBO_PREFIX = "RPS"


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def release_gpu() -> None:
    cp.get_default_memory_pool().free_all_blocks()
    cp.get_default_pinned_memory_pool().free_all_blocks()
    gc.collect()


def validate_scrublet_result(sample_adata: ad.AnnData, input_names: pd.Index) -> None:
    required = {"doublet_score", "predicted_doublet"}
    if not required.issubset(sample_adata.obs.columns):
        raise ValueError("missing required Scrublet output columns")
    if not sample_adata.obs_names.equals(input_names):
        raise ValueError("Scrublet output cell IDs do not match input cell IDs")
    if sample_adata.obs["doublet_score"].isna().any():
        raise ValueError("doublet_score contains missing values")
    if sample_adata.obs["predicted_doublet"].isna().any():
        raise ValueError("predicted_doublet contains missing values")


def run_gpu_scrublet(sample_adata: ad.AnnData) -> ad.AnnData:
    input_names = sample_adata.obs_names.copy()
    # CuPy 14 sparse matrices intentionally reject integer dtypes.  Scrublet
    # receives an exact float32 copy; the canonical QC object's raw integer
    # counts are never replaced.
    if sample_adata.X.dtype not in (np.float32, np.float64):
        sample_adata.X = sample_adata.X.astype(np.float32)
    rsc.get.anndata_to_GPU(sample_adata)
    rsc.pp.scrublet(sample_adata, random_state=SEED, verbose=False)
    rsc.get.anndata_to_CPU(sample_adata)
    validate_scrublet_result(sample_adata, input_names)
    return sample_adata


def run_cpu_scrublet(sample_adata: ad.AnnData) -> ad.AnnData:
    input_names = sample_adata.obs_names.copy()
    sc.pp.scrublet(sample_adata, random_state=SEED, verbose=False)
    validate_scrublet_result(sample_adata, input_names)
    return sample_adata


def write_qc_figures(adata: ad.AnnData) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = FIGURE_DIR
    sc.settings.file_format_figs = "pdf"
    sc.set_figure_params(figsize=(9, 3), dpi=150, fontsize=8)
    metrics = [
        "n_genes_by_counts",
        "total_counts",
        "pct_counts_MT",
        "pct_counts_RIBO",
        "doublet_score",
    ]
    for metric in metrics:
        sc.pl.violin(
            adata,
            keys=[metric],
            groupby="sample",
            jitter=0,
            rotation=90,
            show=False,
            save=f"_{metric}_by_sample.pdf",
        )


def main() -> None:
    start = time.time()
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    OUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(INPUT_H5AD)
    if not sparse.isspmatrix_csr(adata.X):
        adata.X = sparse.csr_matrix(adata.X)
    if not adata.obs_names.is_unique:
        raise ValueError("Input obs_names are not unique.")
    required_obs = {"sample", "series", "status", "original_barcode"}
    missing = required_obs - set(adata.obs.columns)
    if missing:
        raise ValueError(f"Missing required obs columns: {sorted(missing)}")

    preflight_rows = []
    try:
        props = cp.cuda.runtime.getDeviceProperties(0)
        test = adata[:2, :3].copy()
        if test.X.dtype not in (np.float32, np.float64):
            test.X = test.X.astype(np.float32)
        rsc.get.anndata_to_GPU(test)
        gpu_type = type(test.X).__name__
        rsc.get.anndata_to_CPU(test)
        preflight_rows.append(
            {
                "gpu_visible": True,
                "device": props["name"].decode(),
                "rapids_singlecell_import": True,
                "minimal_anndata_transfer": True,
                "gpu_matrix_type": gpu_type,
                "cuda_runtime_version": cp.cuda.runtime.runtimeGetVersion(),
                "cupy": package_version("cupy-cuda12x"),
                "rapids_singlecell": package_version("rapids-singlecell"),
                "cuml": package_version("cuml-cu12"),
                "cugraph": package_version("cugraph-cu12"),
            }
        )
    except Exception as exc:
        preflight_rows.append(
            {
                "gpu_visible": False,
                "device": "",
                "rapids_singlecell_import": True,
                "minimal_anndata_transfer": False,
                "gpu_matrix_type": "",
                "cuda_runtime_version": "",
                "cupy": package_version("cupy-cuda12x"),
                "rapids_singlecell": package_version("rapids-singlecell"),
                "cuml": package_version("cuml-cu12"),
                "cugraph": package_version("cugraph-cu12"),
                "error": repr(exc),
            }
        )
        raise
    pd.DataFrame(preflight_rows).to_csv(TABLE_DIR / "gpu_preflight.csv", index=False)

    before_counts = adata.obs.groupby("sample", observed=True).size().rename("n_input")
    initial_mask = adata.obs["nFeature_RNA"].astype(float) > MIN_GENES_INITIAL
    adata = adata[initial_mask].copy()
    adata.obs["doublet_score"] = np.nan
    adata.obs["predicted_doublet"] = pd.Series(
        False, index=adata.obs_names, dtype=bool
    )

    scrublet_rows: list[dict[str, object]] = []
    failed_rows: list[dict[str, object]] = []
    failed_samples: list[str] = []
    sample_order = adata.obs["sample"].astype(str).drop_duplicates().tolist()
    for order, sample_id in enumerate(sample_order, start=1):
        sample_mask = adata.obs["sample"].astype(str).eq(sample_id)
        sample_input = adata[sample_mask].copy()
        input_n = sample_input.n_obs
        t0 = time.time()
        backend = "rapids_singlecell_gpu"
        attempt = "gpu_initial"
        error_summary = ""
        print(
            f"[Scrublet {order:02d}/{len(sample_order):02d}] "
            f"sample={sample_id} n_cells={input_n}",
            flush=True,
        )
        try:
            result = run_gpu_scrublet(sample_input)
        except (cp.cuda.memory.OutOfMemoryError, MemoryError) as exc:
            release_gpu()
            backend = "scanpy_cpu"
            attempt = "gpu_oom_then_cpu_same_method"
            error_summary = repr(exc)
            fresh = adata[sample_mask].copy()
            try:
                result = run_cpu_scrublet(fresh)
            except Exception as cpu_exc:
                failed_samples.append(sample_id)
                failed_rows.append(
                    {
                        "sample": sample_id,
                        "n_cells_excluded": input_n,
                        "failure_reason": f"GPU OOM: {exc!r}; CPU Scrublet: {cpu_exc!r}",
                        "action": "exclude_entire_sample",
                    }
                )
                print(
                    f"[Scrublet {order:02d}/{len(sample_order):02d}] "
                    f"FAILED sample={sample_id}: {failed_rows[-1]['failure_reason']}",
                    flush=True,
                )
                continue
        except Exception as exc:
            release_gpu()
            attempt = "gpu_initial_then_clean_retry"
            error_summary = repr(exc)
            fresh = adata[sample_mask].copy()
            try:
                result = run_gpu_scrublet(fresh)
            except Exception as retry_exc:
                failed_samples.append(sample_id)
                failed_rows.append(
                    {
                        "sample": sample_id,
                        "n_cells_excluded": input_n,
                        "failure_reason": f"GPU initial: {exc!r}; clean GPU retry: {retry_exc!r}",
                        "action": "exclude_entire_sample",
                    }
                )
                print(
                    f"[Scrublet {order:02d}/{len(sample_order):02d}] "
                    f"FAILED sample={sample_id}: {failed_rows[-1]['failure_reason']}",
                    flush=True,
                )
                continue

        validate_scrublet_result(result, sample_input.obs_names)
        adata.obs.loc[result.obs_names, "doublet_score"] = result.obs[
            "doublet_score"
        ].to_numpy(dtype=float)
        adata.obs.loc[result.obs_names, "predicted_doublet"] = result.obs[
            "predicted_doublet"
        ].to_numpy(dtype=bool)
        scrublet_rows.append(
            {
                "sample_order": order,
                "sample": sample_id,
                "n_input": input_n,
                "n_output": result.n_obs,
                "n_predicted_doublets": int(result.obs["predicted_doublet"].sum()),
                "backend": backend,
                "attempt_route": attempt,
                "fallback_error_summary": error_summary,
                "seed": SEED,
                "elapsed_seconds": time.time() - t0,
                "valid_complete_calls": True,
            }
        )
        print(
            f"[Scrublet {order:02d}/{len(sample_order):02d}] completed "
            f"backend={backend} doublets={scrublet_rows[-1]['n_predicted_doublets']} "
            f"elapsed={scrublet_rows[-1]['elapsed_seconds']:.1f}s",
            flush=True,
        )
        del sample_input, result
        release_gpu()

    if failed_samples:
        adata = adata[~adata.obs["sample"].astype(str).isin(failed_samples)].copy()
    if adata.n_obs == 0:
        raise RuntimeError("No samples remain after excluding failed Scrublet samples.")
    if adata.obs[["doublet_score", "predicted_doublet"]].isna().any().any():
        raise RuntimeError("Retained samples do not all have complete Scrublet outputs.")

    pd.DataFrame(scrublet_rows).to_csv(
        TABLE_DIR / "02_scrublet_parameters.csv", index=False
    )
    pd.DataFrame(
        failed_rows,
        columns=["sample", "n_cells_excluded", "failure_reason", "action"],
    ).to_csv(TABLE_DIR / "02_scrublet_failed_sample_exclusions.csv", index=False)

    n_before_doublet = adata.n_obs
    adata = adata[~adata.obs["predicted_doublet"].astype(bool)].copy()
    n_after_doublet = adata.n_obs

    genes_before = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS_GENE)
    after_min_cells = adata.n_vars
    sc.pp.filter_genes(adata, min_counts=MIN_COUNTS_GENE)
    after_min_counts = adata.n_vars

    adata.var["MT"] = adata.var_names.str.startswith(MT_PREFIX)
    adata.var["RIBO"] = adata.var_names.str.startswith(RIBO_PREFIX)
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["MT", "RIBO"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    final_mask = (
        (adata.obs["n_genes_by_counts"] > MIN_GENES_FINAL)
        & (adata.obs["n_genes_by_counts"] <= MAX_GENES_FINAL)
        & (adata.obs["pct_counts_MT"] < MAX_PCT_MT)
    )
    n_before_final = adata.n_obs
    adata = adata[final_mask].copy()
    n_after_final = adata.n_obs
    adata.uns["expression_contract"] = {
        "X": "raw integer UMI counts after published QC plus sample-wise Scrublet",
        "raw_counts_preserved_in": ".X of adata_qc.h5ad",
    }
    adata.uns["qc_parameters"] = {
        "min_genes_initial_strict_gt": MIN_GENES_INITIAL,
        "min_cells_gene": MIN_CELLS_GENE,
        "min_counts_gene": MIN_COUNTS_GENE,
        "min_genes_final_strict_gt": MIN_GENES_FINAL,
        "max_genes_final": MAX_GENES_FINAL,
        "max_pct_mt_strict_lt": MAX_PCT_MT,
        "mt_prefix": MT_PREFIX,
        "ribo_prefix": RIBO_PREFIX,
        "seed": SEED,
    }

    after_counts = adata.obs.groupby("sample", observed=True).size().rename("n_final")
    doublet_counts = (
        pd.DataFrame(scrublet_rows)
        .set_index("sample")["n_predicted_doublets"]
        .rename("n_predicted_doublets")
    )
    sample_summary = pd.concat(
        [before_counts, doublet_counts, after_counts], axis=1
    ).fillna(0)
    sample_summary["n_removed_total"] = (
        sample_summary["n_input"] - sample_summary["n_final"]
    )
    sample_summary.to_csv(TABLE_DIR / "qc_cell_counts_by_sample.csv")

    pd.DataFrame(
        [
            {"parameter": "min_genes_initial_strict_gt", "value": MIN_GENES_INITIAL},
            {"parameter": "min_cells_gene", "value": MIN_CELLS_GENE},
            {"parameter": "min_counts_gene", "value": MIN_COUNTS_GENE},
            {"parameter": "min_genes_final_strict_gt", "value": MIN_GENES_FINAL},
            {"parameter": "max_genes_final", "value": MAX_GENES_FINAL},
            {"parameter": "max_pct_mt_strict_lt", "value": MAX_PCT_MT},
            {"parameter": "seed", "value": SEED},
        ]
    ).to_csv(TABLE_DIR / "qc_thresholds.csv", index=False)

    pd.DataFrame(
        [
            {
                "genes_before": genes_before,
                "genes_after_min_cells": after_min_cells,
                "genes_after_min_counts": after_min_counts,
                "cells_before_doublet_filter": n_before_doublet,
                "cells_after_doublet_filter": n_after_doublet,
                "cells_before_final_filter": n_before_final,
                "cells_after_final_filter": n_after_final,
                "failed_samples": ";".join(failed_samples),
                "n_failed_samples": len(failed_samples),
            }
        ]
    ).to_csv(TABLE_DIR / "qc_summary.csv", index=False)

    backend_rows = [
        {
            "step": "minimal_anndata_gpu_transfer",
            "planned_backend": "rapids_gpu",
            "attempted_backend": "rapids_gpu",
            "status": "completed",
            "error_summary": "",
            "fallback_backend": "",
            "clean_input_reloaded": True,
            "final_backend_for_rerun": "rapids_gpu",
        },
        {
            "step": "sample_wise_scrublet",
            "planned_backend": "rapids_gpu",
            "attempted_backend": "rapids_gpu",
            "status": "completed" if not failed_samples else "completed_with_failed_sample_exclusions",
            "error_summary": "",
            "fallback_backend": "per-sample CPU only after recorded GPU OOM",
            "clean_input_reloaded": True,
            "final_backend_for_rerun": "see 02_scrublet_parameters.csv",
        },
    ]
    pd.DataFrame(backend_rows).to_csv(
        TABLE_DIR / "gpu_backend_capability_summary.csv", index=False
    )

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
        "package_mirror_note": "base packages from TUNA; RAPIDS 25.04 repair from official PyPI after TUNA stalled",
    }
    (TABLE_DIR / "package_versions.txt").write_text(
        "\n".join(f"{k}={v}" for k, v in versions.items()) + "\n",
        encoding="utf-8",
    )
    (TABLE_DIR / "readme.txt").write_text(
        f"""BRCA QC and sample-wise Scrublet

Input: {INPUT_H5AD}
Code: {CODE_PATH}
Output: {OUT_H5AD}

The local processed matrix already passes its published cell-level boundaries
(nCount_RNA > 250, nFeature_RNA > 200, percent.mito < 20). The workflow keeps
those boundaries, adds sample-wise Scrublet, removes predicted doublets, filters
genes using min_cells={MIN_CELLS_GENE} and min_counts={MIN_COUNTS_GENE}, and
recalculates MT-/RPS QC metrics. No cell downsampling is performed.
""",
        encoding="utf-8",
    )

    write_qc_figures(adata)
    adata.write_h5ad(OUT_H5AD, compression="gzip")
    report = {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_samples": int(adata.obs["sample"].nunique()),
        "n_predicted_doublets_removed": int(n_before_doublet - n_after_doublet),
        "n_final_qc_removed": int(n_before_final - n_after_final),
        "failed_samples": failed_samples,
        "elapsed_seconds": time.time() - start,
    }
    (TABLE_DIR / "qc_completion.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
