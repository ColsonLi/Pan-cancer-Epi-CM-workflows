"""Integration step: normalize/log/HVG/regress/scale/PCA/Harmony.

Per SKILL.md default order:
  1. Read adata_qc.h5ad
  2. rsc.get.anndata_to_GPU
  3. normalize_total (target_sum=1e4)
  4. log1p
  5. highly_variable_genes (n_top_genes=3000, user-specified)
  6. adata.raw = adata (preserve all genes)
  7. Subset to HVGs
  8. regress_out (total_counts, pct_counts_mt)
  9. scale (max_value=10)
 10. PCA (n_comps=50)
 11. harmony_integrate (key='sample', basis='X_pca', adjusted_basis='X_pca_inte')
 12. Save adata_harmony.h5ad

Outputs:
  - h5ad/04-integration-harmony/adata_harmony.h5ad
  - tables/04-integration-harmony/{01..08}_*_parameters.csv
  - tables/04-integration-harmony/integration_harmony_summary.csv
  - tables/04-integration-harmony/package_versions.txt
"""
from __future__ import annotations

import json
import sys
import time
import gc
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

SEED = 42
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
IN_H5AD = ROOT / "epi-cm-core-workflow/h5ad/01-celltype_integration_clustering/03-qc/adata_qc.h5ad"
OUT_H5AD = ROOT / "epi-cm-core-workflow/h5ad/01-celltype_integration_clustering/04-integration-harmony/adata_harmony.h5ad"
TAB = ROOT / "epi-cm-core-workflow/tables/01-celltype_integration_clustering/04-integration-harmony"
TAB.mkdir(parents=True, exist_ok=True)

# User-specified / standard params
N_TOP_GENES = 3000
TARGET_SUM = 1e4
SCALE_MAX = 10
N_COMPS = 50
REGRESS_KEYS = ["total_counts", "pct_counts_mt"]
HARMONY_KEY = "sample"
PCA_BASIS = "X_pca"
INTE_BASIS = "X_pca_inte"


def write_param(step: str, before: tuple, after: tuple, params: dict, code_file: str, notes: str = ""):
    df = pd.DataFrame([{
        "step": step,
        "input_h5ad": str(IN_H5AD) if "in-memory" not in params.get("input", "") else "in-memory",
        "output_h5ad_or_object": str(OUT_H5AD) if step in ("08_harmony_integrate",) else "in-memory adata",
        "n_obs_before": before[0],
        "n_obs_after": after[0],
        "n_vars_before": before[1],
        "n_vars_after": after[1],
        "parameters": json.dumps(params, default=str),
        "backend": params.get("backend", "rsc (GPU)"),
        "code_file": code_file,
        "random_seed": SEED,
        "notes": notes,
    }])
    df.to_csv(TAB / f"{step}_parameters.csv", index=False)
    print(f"[integ] wrote {step}_parameters.csv", flush=True)


def main() -> None:
    t0 = time.time()
    print(f"[integ] reading {IN_H5AD.name}…", flush=True)
    adata = sc.read_h5ad(IN_H5AD)
    print(f"[integ] loaded {adata.shape} in {time.time()-t0:.1f}s", flush=True)
    n_start, v_start = adata.shape

    # Move to GPU
    import rapids_singlecell as rsc
    print("[integ] moving to GPU…", flush=True)
    rsc.get.anndata_to_GPU(adata)
    print(f"[integ] on GPU. shape={adata.shape}", flush=True)

    # 1. normalize_total
    rsc.pp.normalize_total(adata, target_sum=TARGET_SUM)
    write_param("01_normalize_total", (n_start, v_start), adata.shape,
                {"target_sum": TARGET_SUM, "backend": "rsc.pp.normalize_total (GPU)"},
                "epi-cm-core-workflow/codes/01-celltype_integration_clustering/04-integration-harmony/01_integrate.py",
                "")

    # 2. log1p
    rsc.pp.log1p(adata)
    write_param("02_log1p", adata.shape, adata.shape,
                {"backend": "rsc.pp.log1p (GPU)"},
                "epi-cm-core-workflow/codes/01-celltype_integration_clustering/04-integration-harmony/01_integrate.py",
                "")

    # 3. HVG
    rsc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_GENES, flavor="seurat", batch_key=HARMONY_KEY)
    n_hvg = int(adata.var["highly_variable"].sum())
    print(f"[integ] HVG selected: {n_hvg} genes", flush=True)
    write_param("03_highly_variable_genes", adata.shape, adata.shape,
                {"n_top_genes": N_TOP_GENES, "flavor": "seurat", "batch_key": HARMONY_KEY,
                 "actual_n_hvg": n_hvg, "backend": "rsc.pp.highly_variable_genes (GPU)"},
                "epi-cm-core-workflow/codes/01-celltype_integration_clustering/04-integration-harmony/01_integrate.py",
                "user-specified n_top_genes=3000")

    # 4. store raw (before HVG subset)
    adata.raw = adata
    raw_n_vars = adata.raw.n_vars
    write_param("04_raw_assignment_and_hvg_subset_record", adata.shape, adata.shape,
                {"raw_n_vars": raw_n_vars, "note": "adata.raw = adata set before HVG subset"},
                "epi-cm-core-workflow/codes/01-celltype_integration_clustering/04-integration-harmony/01_integrate.py",
                "")

    # 5. subset to HVGs
    adata = adata[:, adata.var["highly_variable"]].copy()
    print(f"[integ] after HVG subset: {adata.shape}", flush=True)
    write_param("04b_hvg_subset", (n_start, v_start), adata.shape,
                {"subset_rule": "highly_variable == True"},
                "epi-cm-core-workflow/codes/01-celltype_integration_clustering/04-integration-harmony/01_integrate.py",
                "")

    # 6. regress_out
    print("[integ] regress_out (total_counts, pct_counts_mt)…", flush=True)
    t1 = time.time()
    rsc.pp.regress_out(adata, keys=REGRESS_KEYS)
    print(f"[integ] regress_out done in {time.time()-t1:.1f}s", flush=True)
    write_param("05_regress_out", adata.shape, adata.shape,
                {"keys": REGRESS_KEYS, "backend": "rsc.pp.regress_out (GPU)"},
                "epi-cm-core-workflow/codes/01-celltype_integration_clustering/04-integration-harmony/01_integrate.py",
                "")

    # 7. scale
    print("[integ] scale…", flush=True)
    t1 = time.time()
    rsc.pp.scale(adata, max_value=SCALE_MAX)
    print(f"[integ] scale done in {time.time()-t1:.1f}s", flush=True)
    write_param("06_scale", adata.shape, adata.shape,
                {"max_value": SCALE_MAX, "backend": "rsc.pp.scale (GPU)"},
                "epi-cm-core-workflow/codes/01-celltype_integration_clustering/04-integration-harmony/01_integrate.py",
                "")

    # 8. PCA (CPU scanpy to avoid RMM 6.67GB contiguous alloc OOM)
    print("[integ] moving to CPU for PCA (avoid RMM frag OOM)…", flush=True)
    rsc.get.anndata_to_CPU(adata)
    gc.collect()
    t1 = time.time()
    sc.tl.pca(adata, n_comps=N_COMPS, random_state=SEED, svd_solver="arpack")
    print(f"[integ] PCA (CPU scanpy) done in {time.time()-t1:.1f}s; X_pca shape={adata.obsm['X_pca'].shape}", flush=True)
    write_param("07_pca", adata.shape, adata.shape,
                {"n_comps": N_COMPS, "random_state": SEED, "svd_solver": "arpack",
                 "backend": "sc.tl.pca (CPU, fallback from rsc OOM)"},
                "epi-cm-core-workflow/codes/01-celltype_integration_clustering/04-integration-harmony/01_integrate.py",
                "rsc.tl.pca OOM (RMM frag); switched to sc.tl.pca arpack on CPU")

    # 9. Harmony (back on GPU)
    print("[integ] moving to GPU for Harmony…", flush=True)
    rsc.get.anndata_to_GPU(adata)
    print("[integ] Harmony…", flush=True)
    t1 = time.time()
    rsc.pp.harmony_integrate(adata, key=HARMONY_KEY, basis=PCA_BASIS, adjusted_basis=INTE_BASIS, random_state=SEED)
    print(f"[integ] Harmony done in {time.time()-t1:.1f}s; X_pca_inte shape={adata.obsm[INTE_BASIS].shape}", flush=True)
    write_param("08_harmony_integrate", adata.shape, adata.shape,
                {"key": HARMONY_KEY, "basis": PCA_BASIS, "adjusted_basis": INTE_BASIS,
                 "random_state": SEED, "backend": "rsc.pp.harmony_integrate (GPU)"},
                "epi-cm-core-workflow/codes/01-celltype_integration_clustering/04-integration-harmony/01_integrate.py",
                "")

    # 10. write output
    print(f"[integ] writing {OUT_H5AD.name}…", flush=True)
    adata.write_h5ad(OUT_H5AD, compression="gzip")
    print(f"[integ] done. final shape={adata.shape}, total={time.time()-t0:.1f}s", flush=True)

    # Summary
    summary = {
        "n_obs_start": n_start,
        "n_vars_start": v_start,
        "n_obs_final": adata.n_obs,
        "n_vars_final_hvg": adata.n_vars,
        "n_raw_vars_preserved": raw_n_vars,
        "n_hvg_actual": n_hvg,
        "n_pcs": N_COMPS,
        "harmony_key": HARMONY_KEY,
        "harmony_basis": PCA_BASIS,
        "harmony_adjusted_basis": INTE_BASIS,
        "regress_keys": REGRESS_KEYS,
        "scale_max": SCALE_MAX,
        "normalize_target_sum": TARGET_SUM,
        "hvg_n_top_genes_requested": N_TOP_GENES,
        "total_runtime_sec": round(time.time() - t0, 1),
        "random_seed": SEED,
    }
    pd.DataFrame([summary]).to_csv(TAB / "integration_harmony_summary.csv", index=False)
    print(f"[integ] integration_harmony_summary.csv written", flush=True)
    print(json.dumps(summary, indent=2), flush=True)

    # package_versions
    pkg = TAB / "package_versions.txt"
    with open(pkg, "w") as f:
        f.write("Integration step — Harmony\n")
        f.write(f"python: {sys.version.split()[0]}\n")
        f.write(f"scanpy: {sc.__version__}\n")
        f.write(f"rapids_singlecell: {rsc.__version__}\n")
        f.write(f"numpy: {np.__version__}\n")
        f.write(f"pandas: {pd.__version__}\n")
        f.write(f"random_seed: {SEED}\n")
        f.write(f"hvg_n_top_genes: {N_TOP_GENES}\n")
        f.write(f"n_pcs: {N_COMPS}\n")
        f.write(f"regress_keys: {REGRESS_KEYS}\n")
        f.write(f"harmony_key: {HARMONY_KEY}\n")
    print(f"[integ] package_versions.txt written", flush=True)


if __name__ == "__main__":
    main()