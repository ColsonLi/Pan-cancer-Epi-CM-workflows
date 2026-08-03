#!/usr/bin/env python3
"""Run Scrublet per sample after the initial cell filter."""

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
INPUT_H5AD = (
    WORKFLOW_ROOT
    / "h5ad/01-celltype_integration_clustering/03-qc/adata_initial_cell_filtered.h5ad"
)
H5AD_DIR = WORKFLOW_ROOT / "h5ad/01-celltype_integration_clustering/03-qc"
TABLE_DIR = WORKFLOW_ROOT / "tables/01-celltype_integration_clustering/03-qc"
OUTPUT_SCORED = H5AD_DIR / "adata_scrublet_scored.h5ad"
OUTPUT_FILTERED = H5AD_DIR / "adata_doublet_filtered.h5ad"
CODE_FILE = Path(__file__)


def _gpu_preflight() -> tuple[bool, dict[str, str]]:
    info: dict[str, str] = {
        "step": "scrublet_preflight",
        "planned_backend": "rapids_singlecell",
        "random_seed": str(SEED),
        "code_file": str(CODE_FILE),
    }
    try:
        import cupy as cp
        import rapids_singlecell as rsc
        import scipy.sparse as sp

        tiny = ad.AnnData(X=sp.csr_matrix(np.array([[1, 0, 2], [0, 1, 0]], dtype=np.float32)))
        rsc.get.anndata_to_GPU(tiny)
        gpu_type = type(tiny.X).__name__
        rsc.get.anndata_to_CPU(tiny)
        info.update(
            {
                "status": "ok",
                "rapids_singlecell_version": rsc.__version__,
                "cupy_version": cp.__version__,
                "minimal_transfer_result": gpu_type,
            }
        )
        return True, info
    except Exception as exc:  # pragma: no cover - recorded for runtime provenance
        info.update({"status": "failed", "error_summary": repr(exc)})
        return False, info


def _is_oom(exc: BaseException) -> bool:
    text = f"{type(exc).__name__}: {exc}".lower()
    return "out of memory" in text or "cuda_error_out_of_memory" in text or "memoryerror" in text


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    H5AD_DIR.mkdir(parents=True, exist_ok=True)
    if OUTPUT_SCORED.exists() or OUTPUT_FILTERED.exists():
        raise FileExistsError("Scrublet output exists; cleanup must be explicit before rerun.")

    gpu_ok, preflight = _gpu_preflight()
    pd.DataFrame([preflight]).to_csv(TABLE_DIR / "gpu_backend_capability_summary.csv", index=False)
    if not gpu_ok:
        raise RuntimeError(
            "rapids_singlecell preflight failed with visible GPU context; stopping before CPU fallback."
        )

    import cupy as cp
    import rapids_singlecell as rsc

    adata = ad.read_h5ad(INPUT_H5AD)
    adata.obs["doublet_score"] = np.nan
    adata.obs["predicted_doublet"] = False
    adata.obs["scrublet_backend"] = "not_run"
    adata.obs["scrublet_status"] = "not_run"

    sample_order = adata.obs["sample"].astype(str).value_counts().index.tolist()
    rows = []
    failed_samples: list[str] = []

    for order_i, sample_id in enumerate(sample_order, start=1):
        mask = adata.obs["sample"].astype(str).to_numpy() == sample_id
        n_in = int(mask.sum())
        sample_adata = adata[mask].copy()
        backend = "rsc"
        status = "completed"
        error_summary = ""
        n_doublet = 0
        try:
            sample_adata.X = sample_adata.X.astype(np.float32)
            rsc.get.anndata_to_GPU(sample_adata)
            rsc.pp.scrublet(sample_adata, random_state=SEED, verbose=False)
            rsc.get.anndata_to_CPU(sample_adata)
        except Exception as exc:
            error_summary = "".join(traceback.format_exception_only(type(exc), exc)).strip()
            try:
                rsc.get.anndata_to_CPU(sample_adata)
            except Exception:
                pass
            cp.get_default_memory_pool().free_all_blocks()
            if _is_oom(exc):
                backend = "scanpy_cpu_after_gpu_oom"
                try:
                    sample_adata = adata[mask].copy()
                    sample_adata.X = sample_adata.X.astype(np.float32)
                    sc.pp.scrublet(sample_adata, random_state=SEED, verbose=False)
                    status = "completed_cpu_fallback_after_gpu_oom"
                except Exception as cpu_exc:
                    status = "failed"
                    error_summary += " | CPU fallback: " + "".join(
                        traceback.format_exception_only(type(cpu_exc), cpu_exc)
                    ).strip()
            else:
                status = "failed_non_oom_rsc"

        if status.startswith("completed") and {
            "doublet_score",
            "predicted_doublet",
        }.issubset(sample_adata.obs.columns):
            scores = sample_adata.obs["doublet_score"].astype(float)
            calls = sample_adata.obs["predicted_doublet"].astype(bool)
            adata.obs.loc[sample_adata.obs_names, "doublet_score"] = scores.to_numpy()
            adata.obs.loc[sample_adata.obs_names, "predicted_doublet"] = calls.to_numpy()
            adata.obs.loc[sample_adata.obs_names, "scrublet_backend"] = backend
            adata.obs.loc[sample_adata.obs_names, "scrublet_status"] = status
            n_doublet = int(calls.sum())
        else:
            failed_samples.append(sample_id)
            adata.obs.loc[sample_adata.obs_names, "scrublet_status"] = status
            adata.obs.loc[sample_adata.obs_names, "scrublet_backend"] = backend

        rows.append(
            {
                "step": "scrublet_per_sample",
                "sample": sample_id,
                "sample_order": order_i,
                "input_h5ad": str(INPUT_H5AD),
                "n_obs_before": n_in,
                "n_obs_after": n_in if status.startswith("completed") else 0,
                "n_vars_before": adata.n_vars,
                "n_vars_after": adata.n_vars,
                "batch_key": "sample",
                "backend": backend,
                "status": status,
                "n_predicted_doublet": n_doublet,
                "sim_doublet_ratio": 2.0,
                "expected_doublet_rate": 0.05,
                "stdev_doublet_rate": 0.02,
                "n_prin_comps": 30,
                "threshold": "auto",
                "code_file": str(CODE_FILE),
                "random_seed": SEED,
                "error_summary": error_summary,
            }
        )
        pd.DataFrame(rows).to_csv(TABLE_DIR / "02_scrublet_parameters.csv", index=False)
        del sample_adata
        cp.get_default_memory_pool().free_all_blocks()
        gc.collect()

    adata.write_h5ad(OUTPUT_SCORED, compression="lzf")
    if failed_samples:
        keep = ~adata.obs["sample"].astype(str).isin(failed_samples)
    else:
        keep = np.ones(adata.n_obs, dtype=bool)
    keep = keep & ~adata.obs["predicted_doublet"].astype(bool).to_numpy()
    filtered = adata[keep].copy()
    filtered.write_h5ad(OUTPUT_FILTERED, compression="lzf")

    doublet_params = pd.DataFrame(
        [
            {
                "step": "doublet_filter",
                "input_h5ad": str(OUTPUT_SCORED),
                "output_h5ad_or_object": str(OUTPUT_FILTERED),
                "n_obs_before": int(adata.n_obs),
                "n_obs_after": int(filtered.n_obs),
                "n_vars_before": int(adata.n_vars),
                "n_vars_after": int(filtered.n_vars),
                "filter_rule": "predicted_doublet == False and scrublet_failed_samples_removed",
                "failed_samples_removed": ";".join(failed_samples),
                "n_failed_samples_removed": len(failed_samples),
                "backend_package": "rapids_singlecell_scrublet_per_sample",
                "code_file": str(CODE_FILE),
                "random_seed": SEED,
            }
        ]
    )
    doublet_params.to_csv(TABLE_DIR / "03_doublet_filter_parameters.csv", index=False)

    with (TABLE_DIR / "readme.txt").open("a", encoding="utf-8") as fh:
        fh.write("\nScrublet per-sample doublet scoring completed.\n")
        fh.write(f"Input: {INPUT_H5AD}\n")
        fh.write(f"Scored output: {OUTPUT_SCORED}\n")
        fh.write(f"Doublet-filtered output: {OUTPUT_FILTERED}\n")
        fh.write(f"Failed samples removed: {len(failed_samples)}\n")

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
    (TABLE_DIR / "package_versions_scrublet.txt").write_text(
        "\n".join(versions) + "\n", encoding="utf-8"
    )
    print(doublet_params.to_string(index=False))


if __name__ == "__main__":
    main()
