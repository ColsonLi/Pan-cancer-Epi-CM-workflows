#!/usr/bin/env python3
"""Normalize/log/HVG/regress/scale/PCA/Harmony integration."""

from __future__ import annotations

import gc
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
INPUT_H5AD = WORKFLOW_ROOT / "h5ad/01-celltype_integration_clustering/03-qc/adata_qc.h5ad"
OUTPUT_H5AD = (
    WORKFLOW_ROOT
    / "h5ad/01-celltype_integration_clustering/04-integration-harmony/adata_harmony.h5ad"
)
TABLE_DIR = WORKFLOW_ROOT / "tables/01-celltype_integration_clustering/04-integration-harmony"
CODE_FILE = Path(__file__)

TARGET_SUM = 1e4
N_TOP_GENES = 3000
REGRESS_KEYS = ["total_counts", "pct_counts_MT"]
SCALE_MAX_VALUE = 10
N_PCS = 50
HARMONY_KEY = "sample"
PCA_BASIS = "X_pca"
HARMONY_BASIS = "X_pca_inte"
FORCE_CPU_SPARSE_PREPROCESS_AFTER_RECORDED_GPU_OOM = True
RECORDED_GPU_OOM_STEP = "log1p"
RECORDED_GPU_OOM_ERROR = (
    "RAPIDS log1p OOM on full sparse matrix after normalize_total; failed to allocate "
    "4809910648 bytes on RTX 2080 Ti."
)
FORCE_CPU_DENSE_STEPS_AFTER_RECORDED_GPU_OOM = True
RECORDED_DENSE_GPU_OOM_STEP = "regress_out"
RECORDED_DENSE_GPU_OOM_ERROR = (
    "RAPIDS regress_out OOM on 896869 x 3000 HVG matrix; failed to allocate "
    "954532736 bytes during cuML SVD on RTX 2080 Ti."
)


def _write_param(name: str, row: dict) -> None:
    pd.DataFrame([row]).to_csv(TABLE_DIR / name, index=False)


def _is_oom(exc: BaseException) -> bool:
    text = f"{type(exc).__name__}: {exc}".lower()
    return "out of memory" in text or "cuda_error_out_of_memory" in text or "memoryerror" in text


def _preflight() -> dict:
    row = {
        "step": "integration_gpu_preflight",
        "planned_backend": "rapids_singlecell",
        "code_file": str(CODE_FILE),
        "random_seed": SEED,
    }
    try:
        import cupy as cp
        import rapids_singlecell as rsc
        import scipy.sparse as sp

        tiny = ad.AnnData(X=sp.csr_matrix(np.array([[1, 0, 2], [0, 3, 0]], dtype=np.float32)))
        rsc.get.anndata_to_GPU(tiny)
        gpu_type = type(tiny.X).__name__
        rsc.pp.normalize_total(tiny, target_sum=TARGET_SUM)
        rsc.pp.log1p(tiny)
        rsc.get.anndata_to_CPU(tiny)
        row.update(
            {
                "status": "ok",
                "rapids_singlecell_version": rsc.__version__,
                "cupy_version": cp.__version__,
                "minimal_transfer_result": gpu_type,
                "fallback_backend": "",
                "error_summary": "",
            }
        )
    except Exception as exc:
        row.update(
            {
                "status": "failed",
                "fallback_backend": "none",
                "error_summary": "".join(traceback.format_exception_only(type(exc), exc)).strip(),
            }
        )
    return row


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
    if OUTPUT_H5AD.exists():
        raise FileExistsError(f"{OUTPUT_H5AD} exists; cleanup must be explicit before rerun.")

    preflight = _preflight()
    backend_rows = [preflight]
    if preflight["status"] != "ok":
        pd.DataFrame(backend_rows).to_csv(TABLE_DIR / "gpu_backend_capability_summary.csv", index=False)
        raise RuntimeError("RAPIDS preflight failed; stopping before integration.")

    import cupy as cp
    import rapids_singlecell as rsc

    adata = ad.read_h5ad(INPUT_H5AD)
    adata.X = adata.X.astype(np.float32)
    n_obs0, n_vars0 = adata.n_obs, adata.n_vars

    final_backend_plan = {}
    gpu_loaded = False
    if FORCE_CPU_SPARSE_PREPROCESS_AFTER_RECORDED_GPU_OOM:
        backend_rows.append(
            {
                "step": RECORDED_GPU_OOM_STEP,
                "planned_backend": "rapids_singlecell",
                "attempted_backend": "rapids_singlecell",
                "status": "gpu_oom_recorded_from_profiling_pass",
                "fallback_backend": "scanpy_cpu",
                "error_summary": RECORDED_GPU_OOM_ERROR,
                "clean_input_reloaded": True,
                "final_backend_for_rerun": "scanpy_cpu_for_sparse_preprocess",
            }
        )
    else:
        try:
            rsc.get.anndata_to_GPU(adata)
            gpu_loaded = True
            backend_rows.append(
                {
                    "step": "anndata_to_GPU",
                    "planned_backend": "rapids_singlecell",
                    "attempted_backend": "rapids_singlecell",
                    "status": "completed",
                    "fallback_backend": "",
                    "error_summary": "",
                    "clean_input_reloaded": True,
                    "final_backend_for_rerun": "rapids_singlecell",
                }
            )
        except Exception as exc:
            backend_rows.append(
                {
                    "step": "anndata_to_GPU",
                    "planned_backend": "rapids_singlecell",
                    "attempted_backend": "rapids_singlecell",
                    "status": "gpu_oom" if _is_oom(exc) else "failed_non_oom",
                    "fallback_backend": "scanpy_cpu" if _is_oom(exc) else "none",
                    "error_summary": "".join(traceback.format_exception_only(type(exc), exc)).strip(),
                    "clean_input_reloaded": True,
                    "final_backend_for_rerun": "scanpy_cpu" if _is_oom(exc) else "blocked",
                }
            )
            if not _is_oom(exc):
                pd.DataFrame(backend_rows).to_csv(TABLE_DIR / "gpu_backend_capability_summary.csv", index=False)
                raise

    if not gpu_loaded:
        # CPU fallback for the full sparse normalize/log/HVG path after recorded GPU OOM.
        sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
        _write_param(
            "01_normalize_total_parameters.csv",
            {
                "step": "normalize_total",
                "input_h5ad": str(INPUT_H5AD),
                "output_h5ad_or_object": "in_memory_adata",
                "n_obs_before": n_obs0,
                "n_obs_after": adata.n_obs,
                "n_vars_before": n_vars0,
                "n_vars_after": adata.n_vars,
                "target_sum": TARGET_SUM,
                "backend_package": "scanpy_cpu_after_gpu_oom",
                "code_file": str(CODE_FILE),
                "random_seed": SEED,
            },
        )
        sc.pp.log1p(adata)
        _write_param(
            "02_log1p_parameters.csv",
            {
                "step": "log1p",
                "input_h5ad": str(INPUT_H5AD),
                "output_h5ad_or_object": "in_memory_adata",
                "n_obs_before": adata.n_obs,
                "n_obs_after": adata.n_obs,
                "n_vars_before": n_vars0,
                "n_vars_after": adata.n_vars,
                "base": "natural_log",
                "backend_package": "scanpy_cpu_after_gpu_oom",
                "code_file": str(CODE_FILE),
                "random_seed": SEED,
            },
        )
        sc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_GENES, flavor="seurat")
        backend_name = "scanpy_cpu_after_gpu_oom"
    else:
        rsc.pp.normalize_total(adata, target_sum=TARGET_SUM)
        _write_param(
            "01_normalize_total_parameters.csv",
            {
                "step": "normalize_total",
                "input_h5ad": str(INPUT_H5AD),
                "output_h5ad_or_object": "in_memory_gpu_adata",
                "n_obs_before": n_obs0,
                "n_obs_after": adata.n_obs,
                "n_vars_before": n_vars0,
                "n_vars_after": adata.n_vars,
                "target_sum": TARGET_SUM,
                "backend_package": "rapids_singlecell",
                "code_file": str(CODE_FILE),
                "random_seed": SEED,
            },
        )
        rsc.pp.log1p(adata)
        _write_param(
            "02_log1p_parameters.csv",
            {
                "step": "log1p",
                "input_h5ad": str(INPUT_H5AD),
                "output_h5ad_or_object": "in_memory_gpu_adata",
                "n_obs_before": adata.n_obs,
                "n_obs_after": adata.n_obs,
                "n_vars_before": n_vars0,
                "n_vars_after": adata.n_vars,
                "base": "natural_log",
                "backend_package": "rapids_singlecell",
                "code_file": str(CODE_FILE),
                "random_seed": SEED,
            },
        )
        rsc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_GENES, flavor="seurat")
        backend_name = "rapids_singlecell"

    n_hvg = int(np.asarray(adata.var["highly_variable"]).sum())
    _write_param(
        "03_highly_variable_genes_parameters.csv",
        {
            "step": "highly_variable_genes",
            "input_h5ad": str(INPUT_H5AD),
            "output_h5ad_or_object": "in_memory_adata_with_hvg_flags",
            "n_obs_before": adata.n_obs,
            "n_obs_after": adata.n_obs,
            "n_vars_before": n_vars0,
            "n_vars_after": adata.n_vars,
            "n_top_genes": N_TOP_GENES,
            "actual_hvg_count": n_hvg,
            "flavor": "seurat",
            "backend_package": backend_name,
            "code_file": str(CODE_FILE),
            "random_seed": SEED,
        },
    )

    if gpu_loaded:
        rsc.get.anndata_to_CPU(adata)
        cp.get_default_memory_pool().free_all_blocks()
    adata.raw = adata.copy()
    adata = adata[:, adata.var["highly_variable"].to_numpy()].copy()
    _write_param(
        "04_raw_assignment_and_hvg_subset_record.csv",
        {
            "step": "raw_assignment_and_hvg_subset",
            "input_h5ad": str(INPUT_H5AD),
            "output_h5ad_or_object": "in_memory_hvg_subset_with_raw",
            "n_obs_before": n_obs0,
            "n_obs_after": adata.n_obs,
            "n_vars_before": n_vars0,
            "n_vars_after": adata.n_vars,
            "raw_assigned_after_normalize_log_hvg_marking": True,
            "raw_n_vars": n_vars0,
            "subset_to_hvg": True,
            "backend_package": "anndata_cpu",
            "code_file": str(CODE_FILE),
            "random_seed": SEED,
        },
    )

    # Downstream dense operations use CPU here after recorded GPU OOM, from a clean rerun.
    if FORCE_CPU_DENSE_STEPS_AFTER_RECORDED_GPU_OOM:
        backend_rows.append(
            {
                "step": RECORDED_DENSE_GPU_OOM_STEP,
                "planned_backend": "rapids_singlecell",
                "attempted_backend": "rapids_singlecell",
                "status": "gpu_oom_recorded_from_profiling_pass",
                "fallback_backend": "scanpy_cpu",
                "error_summary": RECORDED_DENSE_GPU_OOM_ERROR,
                "clean_input_reloaded": True,
                "final_backend_for_rerun": "scanpy_cpu_for_dense_steps",
            }
        )
        sc.pp.regress_out(adata, REGRESS_KEYS)
        regress_backend = "scanpy_cpu_after_gpu_oom"
        sc.pp.scale(adata, max_value=SCALE_MAX_VALUE)
        scale_backend = "scanpy_cpu_after_gpu_oom"
        sc.tl.pca(adata, n_comps=N_PCS, random_state=SEED, svd_solver="arpack")
        pca_backend = "scanpy_cpu_after_gpu_oom"
        import harmonypy as hm

        ho = hm.run_harmony(
            adata.obsm[PCA_BASIS],
            adata.obs,
            vars_use=[HARMONY_KEY],
            random_state=SEED,
            verbose=False,
        )
        z_corr = np.asarray(ho.Z_corr)
        if z_corr.shape == (N_PCS, adata.n_obs):
            z_corr = z_corr.T
        if z_corr.shape != (adata.n_obs, N_PCS):
            raise ValueError(
                f"Unexpected Harmony output shape {z_corr.shape}; expected {(adata.n_obs, N_PCS)}"
            )
        adata.obsm[HARMONY_BASIS] = z_corr.astype(np.float32)
        harmony_backend = "harmonypy_cpu_after_gpu_oom"
    else:
        try:
            rsc.get.anndata_to_GPU(adata)
            rsc.pp.regress_out(adata, REGRESS_KEYS, batchsize=100000)
            regress_backend = "rapids_singlecell"
            rsc.pp.scale(adata, max_value=SCALE_MAX_VALUE)
            scale_backend = "rapids_singlecell"
            rsc.pp.pca(adata, n_comps=N_PCS, random_state=SEED)
            pca_backend = "rapids_singlecell"
            rsc.pp.harmony_integrate(
                adata,
                key=HARMONY_KEY,
                basis=PCA_BASIS,
                adjusted_basis=HARMONY_BASIS,
                random_state=SEED,
                verbose=False,
            )
            harmony_backend = "rapids_singlecell"
            rsc.get.anndata_to_CPU(adata)
            cp.get_default_memory_pool().free_all_blocks()
        except Exception as exc:
            err = "".join(traceback.format_exception_only(type(exc), exc)).strip()
            try:
                rsc.get.anndata_to_CPU(adata)
            except Exception:
                pass
            cp.get_default_memory_pool().free_all_blocks()
            if not _is_oom(exc):
                pd.DataFrame(backend_rows + [{"step": "downstream_gpu", "status": "failed_non_oom", "error_summary": err}]).to_csv(
                    TABLE_DIR / "gpu_backend_capability_summary.csv", index=False
                )
                raise
            backend_rows.append(
                {
                    "step": "downstream_gpu_dense_steps",
                    "planned_backend": "rapids_singlecell",
                    "attempted_backend": "rapids_singlecell",
                    "status": "gpu_oom",
                    "fallback_backend": "scanpy_cpu",
                    "error_summary": err,
                    "clean_input_reloaded": False,
                    "final_backend_for_rerun": "scanpy_cpu_for_dense_steps",
                }
            )
            sc.pp.regress_out(adata, REGRESS_KEYS)
            regress_backend = "scanpy_cpu_after_gpu_oom"
            sc.pp.scale(adata, max_value=SCALE_MAX_VALUE)
            scale_backend = "scanpy_cpu_after_gpu_oom"
            sc.tl.pca(adata, n_comps=N_PCS, random_state=SEED, svd_solver="arpack")
            pca_backend = "scanpy_cpu_after_gpu_oom"
            import harmonypy as hm

            ho = hm.run_harmony(
                adata.obsm[PCA_BASIS],
                adata.obs,
                vars_use=[HARMONY_KEY],
                random_state=SEED,
                verbose=False,
            )
            z_corr = np.asarray(ho.Z_corr)
            if z_corr.shape == (N_PCS, adata.n_obs):
                z_corr = z_corr.T
            if z_corr.shape != (adata.n_obs, N_PCS):
                raise ValueError(
                    f"Unexpected Harmony output shape {z_corr.shape}; expected {(adata.n_obs, N_PCS)}"
                )
            adata.obsm[HARMONY_BASIS] = z_corr.astype(np.float32)
            harmony_backend = "harmonypy_cpu_after_gpu_oom"

    _write_param(
        "05_regress_out_parameters.csv",
        {
            "step": "regress_out",
            "input_h5ad": str(INPUT_H5AD),
            "output_h5ad_or_object": "in_memory_hvg_subset",
            "n_obs_before": adata.n_obs,
            "n_obs_after": adata.n_obs,
            "n_vars_before": adata.n_vars,
            "n_vars_after": adata.n_vars,
            "regress_keys": ";".join(REGRESS_KEYS),
            "backend_package": regress_backend,
            "code_file": str(CODE_FILE),
            "random_seed": SEED,
        },
    )
    _write_param(
        "06_scale_parameters.csv",
        {
            "step": "scale",
            "input_h5ad": str(INPUT_H5AD),
            "output_h5ad_or_object": "in_memory_hvg_subset",
            "n_obs_before": adata.n_obs,
            "n_obs_after": adata.n_obs,
            "n_vars_before": adata.n_vars,
            "n_vars_after": adata.n_vars,
            "max_value": SCALE_MAX_VALUE,
            "backend_package": scale_backend,
            "code_file": str(CODE_FILE),
            "random_seed": SEED,
        },
    )
    _write_param(
        "07_pca_parameters.csv",
        {
            "step": "pca",
            "input_h5ad": str(INPUT_H5AD),
            "output_h5ad_or_object": "obsm_X_pca",
            "n_obs_before": adata.n_obs,
            "n_obs_after": adata.n_obs,
            "n_vars_before": adata.n_vars,
            "n_vars_after": adata.n_vars,
            "n_comps": N_PCS,
            "basis": PCA_BASIS,
            "backend_package": pca_backend,
            "code_file": str(CODE_FILE),
            "random_seed": SEED,
        },
    )
    _write_param(
        "08_harmony_integrate_parameters.csv",
        {
            "step": "harmony_integrate",
            "input_h5ad": str(INPUT_H5AD),
            "output_h5ad_or_object": str(OUTPUT_H5AD),
            "n_obs_before": adata.n_obs,
            "n_obs_after": adata.n_obs,
            "n_vars_before": adata.n_vars,
            "n_vars_after": adata.n_vars,
            "key": HARMONY_KEY,
            "basis": PCA_BASIS,
            "adjusted_basis": HARMONY_BASIS,
            "backend_package": harmony_backend,
            "code_file": str(CODE_FILE),
            "random_seed": SEED,
        },
    )

    adata.write_h5ad(OUTPUT_H5AD, compression="lzf")
    summary = pd.DataFrame(
        [
            {
                "step": "integration_harmony",
                "input_h5ad": str(INPUT_H5AD),
                "output_h5ad": str(OUTPUT_H5AD),
                "n_obs": adata.n_obs,
                "n_vars_hvg": adata.n_vars,
                "raw_present": adata.raw is not None,
                "raw_n_vars": n_vars0,
                "x_pca_present": PCA_BASIS in adata.obsm,
                "x_pca_inte_present": HARMONY_BASIS in adata.obsm,
                "n_pcs": N_PCS,
                "harmony_key": HARMONY_KEY,
                "code_file": str(CODE_FILE),
                "random_seed": SEED,
            }
        ]
    )
    summary.to_csv(TABLE_DIR / "integration_harmony_summary.csv", index=False)
    pd.DataFrame(backend_rows).to_csv(TABLE_DIR / "gpu_backend_capability_summary.csv", index=False)

    with (TABLE_DIR / "readme.txt").open("w", encoding="utf-8") as fh:
        fh.write("04-integration-harmony completed.\n")
        fh.write(f"Input: {INPUT_H5AD}\n")
        fh.write(f"Output: {OUTPUT_H5AD}\n")
        fh.write("No UMAP, neighbors, or Leiden clustering was run in this object.\n")

    versions = [
        f"python={platform.python_version()}",
        f"platform={platform.platform()}",
        f"anndata={ad.__version__}",
        f"scanpy={sc.__version__}",
        f"numpy={np.__version__}",
        f"pandas={pd.__version__}",
        f"rapids_singlecell={rsc.__version__}",
        f"cupy={cp.__version__}",
        "environment=/mnt/disk18t/lr_xcy/riku/crc_val/uv_envs/rapids/.venv",
        f"code_file={CODE_FILE}",
        f"random_seed={SEED}",
    ]
    (TABLE_DIR / "package_versions.txt").write_text("\n".join(versions) + "\n", encoding="utf-8")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
