#!/usr/bin/env python3
"""Recompute the BRCA epithelial subtype Harmony template from adata_qc."""

from __future__ import annotations

import gc
import importlib.metadata
import json
import os
import random
import sys
import time
import warnings
from pathlib import Path

os.environ.setdefault("NUMBA_CUDA_USE_NVIDIA_BINDING", "1")

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

ROOT = Path("/mnt/disk18t/lr_xcy/riku/brca_cm-epi_codex")
WORKFLOW = ROOT / "epi-cm-core-workflow"
BLOCK01 = "01-celltype_integration_clustering"
BLOCK02 = "02-cell_subtype_integration_clustering"
QC_H5AD = WORKFLOW / "h5ad" / BLOCK01 / "03-qc" / "adata_qc.h5ad"
SELECTED_IDS = (
    WORKFLOW
    / "tables"
    / BLOCK02
    / "01-lineage-selection"
    / "epithelial"
    / "selected_cell_ids.csv"
)
OUT_H5AD = (
    WORKFLOW
    / "h5ad"
    / BLOCK02
    / "02-subtype-harmony"
    / "epithelial"
    / "adata_epithelial_harmony.h5ad"
)
TABLE_DIR = WORKFLOW / "tables" / BLOCK02 / "02-subtype-harmony" / "epithelial"
CODE_PATH = (
    WORKFLOW
    / "codes"
    / BLOCK02
    / "02-subtype-harmony"
    / "epithelial"
    / "02_run_epithelial_harmony.py"
)
LINEAGE_SLUG = "epithelial"
TARGET_LABEL = "Epithelial Cells"
LINEAGE_LABELS = {
    "epithelial": "Epithelial Cells",
    "t_cells": "T Cells",
    "myeloid": "Myeloid Cells",
    "b_cells": "B Cells",
    "plasma": "Plasma Cells",
    "endothelial": "Endothelial Cells",
    "stromal": "Stromal Cells",
    "perivascular": "Perivascular Cells",
}

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
HARMONY_RETRY_MAX_ITER = 50


def configure_lineage(lineage_slug: str, code_path: Path | None = None) -> None:
    """Configure this validated Harmony implementation for another broad lineage."""
    global LINEAGE_SLUG, TARGET_LABEL, SELECTED_IDS, OUT_H5AD, TABLE_DIR, CODE_PATH
    if lineage_slug not in LINEAGE_LABELS:
        raise ValueError(f"Unsupported lineage slug: {lineage_slug}")
    LINEAGE_SLUG = lineage_slug
    TARGET_LABEL = LINEAGE_LABELS[lineage_slug]
    SELECTED_IDS = (
        WORKFLOW
        / "tables"
        / BLOCK02
        / "01-lineage-selection"
        / lineage_slug
        / "selected_cell_ids.csv"
    )
    OUT_H5AD = (
        WORKFLOW
        / "h5ad"
        / BLOCK02
        / "02-subtype-harmony"
        / lineage_slug
        / f"adata_{lineage_slug}_harmony.h5ad"
    )
    TABLE_DIR = WORKFLOW / "tables" / BLOCK02 / "02-subtype-harmony" / lineage_slug
    if code_path is not None:
        CODE_PATH = code_path


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
        return {key: cpuify(item) for key, item in value.items()}
    if isinstance(value, list):
        return [cpuify(item) for item in value]
    if isinstance(value, tuple):
        return tuple(cpuify(item) for item in value)
    return value


def convert_adata_aux_to_cpu(adata) -> None:
    for mapping_name in ["obsm", "varm", "obsp", "varp", "layers"]:
        mapping = getattr(adata, mapping_name)
        for key in list(mapping.keys()):
            mapping[key] = cpuify(mapping[key])
    for key in list(adata.uns.keys()):
        adata.uns[key] = cpuify(adata.uns[key])


def main() -> None:
    started = time.time()
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    OUT_H5AD.parent.mkdir(parents=True, exist_ok=True)

    selection = pd.read_csv(SELECTED_IDS, dtype={"cell_id": str})
    if selection.empty or selection["cell_id"].duplicated().any():
        raise ValueError(f"Selected {LINEAGE_SLUG} ID table is empty or duplicated.")
    selected_ids = pd.Index(selection["cell_id"].astype(str))

    adata_qc = sc.read_h5ad(QC_H5AD)
    required_qc = {
        "sample",
        "series",
        "status",
        "original_barcode",
        "total_counts",
        "pct_counts_MT",
    }
    missing = required_qc - set(adata_qc.obs.columns)
    if missing:
        raise ValueError(f"QC expression base lacks columns: {sorted(missing)}")
    if not adata_qc.obs_names.is_unique:
        raise ValueError("QC expression base has non-unique cell IDs.")
    missing_ids = selected_ids.difference(adata_qc.obs_names.astype(str))
    if len(missing_ids):
        raise ValueError(f"Selected {LINEAGE_SLUG} IDs missing from adata_qc: {len(missing_ids)}")
    adata = adata_qc[selected_ids.tolist()].copy()
    del adata_qc
    if not adata.obs_names.equals(selected_ids):
        raise ValueError(f"QC subset order does not exactly match selected {LINEAGE_SLUG} IDs.")

    core_metadata = ["sample", "series", "status", "original_barcode"]
    selection_by_id = selection.set_index("cell_id").loc[selected_ids]
    transfer_audit: list[dict[str, object]] = []
    for column in core_metadata:
        matches = adata.obs[column].astype(str).equals(
            selection_by_id[column].astype(str)
        )
        transfer_audit.append(
            {"column": column, "role": "identity_validation", "all_values_match": matches}
        )
        if not matches:
            raise ValueError(f"Metadata mismatch between selection source and adata_qc: {column}")
    annotation_columns = [column for column in selection.columns if column not in {"cell_id", *core_metadata}]
    for column in annotation_columns:
        values = selection_by_id[column]
        if column in {"leiden_res0p8", "leiden_coarse", "cell_type", "best_rank_type_global"}:
            adata.obs[column] = pd.Categorical(values.astype(str).to_numpy())
        elif column == "score_rank_consistent":
            adata.obs[column] = values.astype(bool).to_numpy()
        else:
            adata.obs[column] = pd.to_numeric(values, errors="raise").to_numpy()
        transfer_audit.append(
            {"column": column, "role": "annotation_transfer", "all_values_match": True}
        )
    pd.DataFrame(transfer_audit).to_csv(
        TABLE_DIR / "annotation_transfer_audit.csv", index=False
    )
    if not adata.obs["leiden_coarse"].astype(str).eq(TARGET_LABEL).all():
        raise ValueError(f"Off-lineage broad labels entered the {LINEAGE_SLUG} subset.")
    if not adata.obs["score_rank_consistent"].astype(bool).all():
        raise ValueError(f"Score/rank-inconsistent cells entered the strict {LINEAGE_SLUG} subset.")

    source_count_dtype = str(adata.X.dtype)
    if adata.X.dtype not in (np.float32, np.float64):
        adata.X = adata.X.astype(np.float32)
    n_samples = int(adata.obs["sample"].astype(str).nunique())
    print(
        f"[{LINEAGE_SLUG} Harmony] input={adata.shape} samples={n_samples} "
        f"source_dtype={source_count_dtype} working_dtype={adata.X.dtype}",
        flush=True,
    )

    backend_rows: list[dict[str, object]] = []
    step_times: dict[str, float] = {}
    t0 = time.time()
    gpu_count = int(cp.cuda.runtime.getDeviceCount())
    if gpu_count < 1:
        raise RuntimeError("No CUDA GPU is visible to CuPy.")
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

    t0 = time.time()
    rsc.pp.normalize_total(adata, target_sum=int(TARGET_SUM))
    rsc.pp.log1p(adata)
    step_times["normalize_log1p"] = time.time() - t0

    t0 = time.time()
    hvg_batch_counts = adata.obs[HVG_BATCH_KEY].astype(str).value_counts().sort_index()
    hvg_ineligible_samples = hvg_batch_counts[hvg_batch_counts < 2].index.tolist()
    hvg_eligible_mask = ~adata.obs[HVG_BATCH_KEY].astype(str).isin(hvg_ineligible_samples)
    hvg_eligibility = hvg_batch_counts.rename("n_lineage_cells").rename_axis("sample").reset_index()
    hvg_eligibility["used_for_hvg_estimation"] = hvg_eligibility["n_lineage_cells"] >= 2
    hvg_eligibility["reason_if_excluded"] = np.where(
        hvg_eligibility["used_for_hvg_estimation"],
        "",
        "singleton batch cannot estimate within-sample variance",
    )
    hvg_eligibility.to_csv(TABLE_DIR / "hvg_batch_eligibility.csv", index=False)
    hvg_source = adata if not hvg_ineligible_samples else adata[hvg_eligible_mask].copy()
    original_var_columns = set(adata.var.columns)
    rsc.pp.highly_variable_genes(
        hvg_source,
        n_top_genes=N_TOP_GENES,
        flavor=HVG_FLAVOR,
        batch_key=HVG_BATCH_KEY,
    )
    if hvg_source is not adata:
        expected_hvg_columns = [
            "highly_variable",
            "means",
            "dispersions",
            "dispersions_norm",
            "highly_variable_nbatches",
            "highly_variable_intersection",
        ]
        transfer_columns = [
            column
            for column in hvg_source.var.columns
            if column in expected_hvg_columns or column not in original_var_columns
        ]
        for column in transfer_columns:
            adata.var[column] = hvg_source.var[column].to_numpy()
        del hvg_source
        gc.collect()
    step_times["highly_variable_genes"] = time.time() - t0
    n_hvg = int(adata.var["highly_variable"].sum())
    if n_hvg < 2000:
        raise ValueError(f"Only {n_hvg} {LINEAGE_SLUG} HVGs were selected.")

    rsc.get.anndata_to_CPU(adata)
    convert_adata_aux_to_cpu(adata)
    adata.raw = adata.copy()
    hvg_mask = adata.var["highly_variable"].to_numpy(dtype=bool)
    adata = adata[:, hvg_mask].copy()
    print(
        f"[{LINEAGE_SLUG} Harmony] raw_genes={adata.raw.n_vars} hvg_matrix={adata.shape}",
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
    t0 = time.time()
    rsc.pp.scale(adata, max_value=SCALE_MAX_VALUE)
    step_times["scale"] = time.time() - t0
    t0 = time.time()
    rsc.tl.pca(
        adata,
        n_comps=PCA_N_COMPS,
        random_state=SEED,
        dtype="float32",
    )
    step_times["pca"] = time.time() - t0

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
    harmony_warnings.extend(str(item.message) for item in caught)
    if any("did not converge" in message.lower() for message in harmony_warnings):
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
        retry_messages = [str(item.message) for item in caught_retry]
        harmony_warnings.extend(f"retry: {message}" for message in retry_messages)
        if any("did not converge" in message.lower() for message in retry_messages):
            raise RuntimeError(f"{LINEAGE_SLUG} Harmony did not converge after max_iter retry.")
    step_times["harmony"] = time.time() - t0

    rsc.get.anndata_to_CPU(adata)
    convert_adata_aux_to_cpu(adata)
    release_gpu()
    if adata.obsm[PCA_BASIS].shape != (adata.n_obs, PCA_N_COMPS):
        raise ValueError(f"Unexpected {LINEAGE_SLUG} PCA shape: {adata.obsm[PCA_BASIS].shape}")
    if adata.obsm[INTEGRATED_BASIS].shape != (adata.n_obs, PCA_N_COMPS):
        raise ValueError(
            f"Unexpected {LINEAGE_SLUG} Harmony shape: {adata.obsm[INTEGRATED_BASIS].shape}"
        )
    if not np.isfinite(adata.obsm[INTEGRATED_BASIS]).all():
        raise ValueError(f"{LINEAGE_SLUG} Harmony embedding contains non-finite values.")
    if adata.raw is None or adata.raw.n_vars != 27716:
        raise ValueError(f"Full normalized/log {LINEAGE_SLUG} raw gene space was not preserved.")
    if {"neighbors", "umap", "leiden"} & set(adata.uns) or "X_umap" in adata.obsm:
        raise ValueError("Harmony template unexpectedly contains graph/clustering outputs.")

    adata.uns["integration_parameters"] = {
        "lineage": LINEAGE_SLUG,
        "target_leiden_coarse": TARGET_LABEL,
        "selection_ids_csv": str(SELECTED_IDS),
        "expression_base_h5ad": str(QC_H5AD),
        "rank_filter_mode": "strict_consistent",
        "target_sum": TARGET_SUM,
        "log1p": True,
        "n_top_genes": N_TOP_GENES,
        "hvg_flavor": HVG_FLAVOR,
        "hvg_batch_key": HVG_BATCH_KEY,
        "hvg_estimation_min_cells_per_batch": 2,
        "hvg_estimation_n_cells": int(hvg_eligible_mask.sum()),
        "hvg_estimation_n_samples": int((hvg_batch_counts >= 2).sum()),
        "hvg_estimation_excluded_singleton_samples": hvg_ineligible_samples,
        "hvg_estimation_cell_retention_rule": (
            "singleton samples excluded only from HVG variance estimation; "
            "all selected cells retained for scaling/PCA/Harmony"
        ),
        "raw_assignment_point": "after normalize_total/log1p/HVG marking before HVG subsetting",
        "regression_keys": REGRESSION_KEYS,
        "regression_batchsize": REGRESSION_BATCHSIZE,
        "scale_max_value": SCALE_MAX_VALUE,
        "pca_n_comps": PCA_N_COMPS,
        "batch_key": BATCH_KEY,
        "pca_basis": PCA_BASIS,
        "integrated_basis": INTEGRATED_BASIS,
        "harmony_max_iter_used": harmony_max_iter_used,
        "seed": SEED,
    }
    adata.uns["expression_contract"] = {
        "X": f"{LINEAGE_SLUG}-HVG regressed and scaled expression",
        "raw": f"full-gene {LINEAGE_SLUG} normalized/log1p expression",
        "raw_counts": str(QC_H5AD) + "::.X for selected IDs",
    }

    pd.DataFrame(
        [
            {
                "input_qc_h5ad": str(QC_H5AD),
                "lineage_slug": LINEAGE_SLUG,
                "target_leiden_coarse": TARGET_LABEL,
                "selected_ids_csv": str(SELECTED_IDS),
                "output_h5ad": str(OUT_H5AD),
                "n_cells": int(adata.n_obs),
                "n_samples": n_samples,
                "n_genes_input": int(adata.raw.n_vars),
                "n_hvg": n_hvg,
                "target_sum": TARGET_SUM,
                "hvg_flavor": HVG_FLAVOR,
                "hvg_batch_key": HVG_BATCH_KEY,
                "hvg_estimation_min_cells_per_batch": 2,
                "hvg_estimation_n_cells": int(hvg_eligible_mask.sum()),
                "hvg_estimation_n_samples": int((hvg_batch_counts >= 2).sum()),
                "hvg_estimation_excluded_singleton_samples": ";".join(hvg_ineligible_samples),
                "regression_keys": ";".join(REGRESSION_KEYS),
                "regression_batchsize": REGRESSION_BATCHSIZE,
                "scale_max_value": SCALE_MAX_VALUE,
                "pca_n_comps": PCA_N_COMPS,
                "batch_key": BATCH_KEY,
                "pca_basis": PCA_BASIS,
                "integrated_basis": INTEGRATED_BASIS,
                "harmony_initial_max_iter": HARMONY_INITIAL_MAX_ITER,
                "harmony_final_max_iter": harmony_max_iter_used,
                "rank_filter_mode": "strict_consistent",
                "seed": SEED,
                "backend": "rapids_singlecell_gpu",
                "code_file": str(CODE_PATH),
            }
        ]
    ).to_csv(TABLE_DIR / "harmony_parameters.csv", index=False)
    pd.DataFrame(
        [{"step": step, "elapsed_seconds": elapsed} for step, elapsed in step_times.items()]
    ).to_csv(TABLE_DIR / "integration_step_times.csv", index=False)
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
    (TABLE_DIR / "harmony_warnings.txt").write_text(
        "\n".join(harmony_warnings) + ("\n" if harmony_warnings else ""),
        encoding="utf-8",
    )

    ens_re = r"^ENS[A-Z]*G[0-9]+(?:\.[0-9]+)?$"
    gene_rows: list[dict[str, object]] = []
    for layer, names in [
        ("var_names", pd.Index(adata.var_names.astype(str))),
        ("raw.var_names", pd.Index(adata.raw.var_names.astype(str))),
    ]:
        ensembl = names.to_series(index=range(len(names))).str.match(ens_re)
        fraction = float(ensembl.mean()) if len(names) else 0.0
        gene_rows.append(
            {
                "input_h5ad": str(OUT_H5AD),
                "object_layer": layer,
                "var_names_type": (
                    "gene_id" if fraction >= 0.9 else "gene_name" if fraction <= 0.1 else "mixed_or_unknown"
                ),
                "n_features": len(names),
                "n_unique": int(names.nunique()),
                "n_duplicates": int(names.duplicated().sum()),
                "ensembl_like_n": int(ensembl.sum()),
                "ensembl_like_fraction": fraction,
                "action_taken": "none",
            }
        )
    pd.DataFrame(gene_rows).to_csv(TABLE_DIR / "gene_identifier_audit.csv", index=False)

    print(f"[{LINEAGE_SLUG} Harmony] writing {OUT_H5AD}", flush=True)
    adata.write_h5ad(OUT_H5AD, compression="gzip")
    report = {
        "pass": True,
        "lineage_slug": LINEAGE_SLUG,
        "target_leiden_coarse": TARGET_LABEL,
        "n_cells": int(adata.n_obs),
        "n_samples": n_samples,
        "n_hvg": n_hvg,
        "hvg_estimation_n_cells": int(hvg_eligible_mask.sum()),
        "hvg_estimation_excluded_singleton_samples": hvg_ineligible_samples,
        "n_raw_genes": int(adata.raw.n_vars),
        "pca_shape": list(adata.obsm[PCA_BASIS].shape),
        "harmony_shape": list(adata.obsm[INTEGRATED_BASIS].shape),
        "harmony_max_iter_used": harmony_max_iter_used,
        "rank_filter_mode": "strict_consistent",
        "elapsed_seconds": time.time() - started,
        "output_h5ad": str(OUT_H5AD),
    }
    (TABLE_DIR / "harmony_completion.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    (TABLE_DIR / "readme.txt").write_text(
        f"""BRCA {TARGET_LABEL} subtype Harmony integration

QC expression base: {QC_H5AD}
Strict selected IDs: {SELECTED_IDS}
Output: {OUT_H5AD}

All {adata.n_obs} strict-consistent {TARGET_LABEL} cells from {n_samples} represented
samples are jointly normalized to 10000, log1p transformed, marked for 3000
sample-aware seurat HVGs, stored in adata.raw before HVG subsetting, regressed
for total_counts and pct_counts_MT, scaled to max_value=10, reduced to 50 PCs,
and Harmony-corrected by sample. No cells are downsampled and no graph, UMAP, or
clustering result is stored in this template.
Samples with fewer than two lineage cells are excluded only from batch-aware
HVG variance estimation because within-sample variance is undefined; every
selected cell remains in the normalized/scaled/PCA/Harmony object. Excluded
HVG-estimation samples: {', '.join(hvg_ineligible_samples) or 'none'}.
""",
        encoding="utf-8",
    )
    (TABLE_DIR / "package_versions.txt").write_text(
        "\n".join(
            [
                f"python={sys.version.split()[0]}",
                f"scanpy={package_version('scanpy')}",
                f"anndata={package_version('anndata')}",
                f"rapids-singlecell={package_version('rapids-singlecell')}",
                f"cupy-cuda12x={package_version('cupy-cuda12x')}",
                f"cuml-cu12={package_version('cuml-cu12')}",
                f"NUMBA_CUDA_USE_NVIDIA_BINDING={os.environ.get('NUMBA_CUDA_USE_NVIDIA_BINDING', '')}",
                f"code={CODE_PATH}",
                f"seed={SEED}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
