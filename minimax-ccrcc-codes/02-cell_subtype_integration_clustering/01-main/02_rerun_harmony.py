"""Re-run B / Mast / Epi subclustering with max_iter_harmony=20.

These three lineages showed 'Harmony did not converge' in the first pass.
This script re-runs only those three with max_iter_harmony=20 to ensure
convergence, overwriting the corresponding h5ad and figures.
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

SEED = 42
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
CONSISTENT = ROOT / "epi-cm-core-workflow/h5ad/01-celltype_integration_clustering/07-score-rank-qc/adata_anno_score_genes_rank_consistent.h5ad"
QC = ROOT / "epi-cm-core-workflow/h5ad/01-celltype_integration_clustering/03-qc/adata_qc.h5ad"
H5AD = ROOT / "epi-cm-core-workflow/h5ad/02-cell_subtype_integration_clustering/01-main"
TAB = ROOT / "epi-cm-core-workflow/tables/02-cell_subtype_integration_clustering/01-main"
FIG = ROOT / "epi-cm-core-workflow/figures/02-cell_subtype_integration_clustering/01-main"

LINEAGE_PARAMS = {
    "B":   dict(pcs=10, nn=25, res=0.2),
    "Mast": dict(pcs=15, nn=40, res=0.8),
    "Epi": dict(pcs=15, nn=30, res=0.6),
}

N_TOP_HVG = 3000
TARGET_SUM = 1e4
SCALE_MAX = 10
REGRESS_KEYS = ["total_counts", "pct_counts_mt"]
MAX_ITER_HARMONY = 20


def process_lineage(lineage: str, ids: pd.Index, params: dict,
                    adata_qc: sc.AnnData) -> None:
    print(f"\n========== {lineage} (pcs={params['pcs']}, nn={params['nn']}, res={params['res']}, max_iter={MAX_ITER_HARMONY}) ==========", flush=True)
    print(f"[{lineage}] n_input_ids: {len(ids)}", flush=True)

    sub = adata_qc[ids].copy()
    print(f"[{lineage}] sub shape: {sub.shape}", flush=True)

    import rapids_singlecell as rsc
    rsc.get.anndata_to_GPU(sub)
    rsc.pp.normalize_total(sub, target_sum=TARGET_SUM)
    rsc.pp.log1p(sub)
    rsc.pp.highly_variable_genes(sub, n_top_genes=N_TOP_HVG)
    sub.raw = sub
    sub = sub[:, sub.var["highly_variable"]].copy()
    rsc.pp.regress_out(sub, keys=REGRESS_KEYS)
    rsc.pp.scale(sub, max_value=SCALE_MAX)
    rsc.tl.pca(sub, n_comps=50, random_state=SEED)
    rsc.pp.harmony_integrate(
        sub, key="sample", basis="X_pca", adjusted_basis="X_pca_inte",
        random_state=SEED, max_iter_harmony=MAX_ITER_HARMONY,
    )
    rsc.pp.neighbors(sub, n_neighbors=params["nn"], n_pcs=params["pcs"], use_rep="X_pca_inte", random_state=SEED)
    leiden_key = "leiden_sub"
    rsc.tl.leiden(sub, resolution=params["res"], key_added=leiden_key, random_state=SEED)
    n_clusters = int(sub.obs[leiden_key].nunique())
    print(f"[{lineage}] leiden n_clusters: {n_clusters}", flush=True)
    rsc.tl.umap(sub, random_state=SEED)
    rsc.get.anndata_to_CPU(sub)
    gc.collect()

    out_h5ad = H5AD / f"adata_{lineage}_subclustered.h5ad"
    sub.write_h5ad(out_h5ad, compression="gzip")
    print(f"[{lineage}] wrote {out_h5ad.name}", flush=True)

    sc.tl.rank_genes_groups(
        sub, groupby=leiden_key, use_raw=True, method="t-test",
        n_genes=None, key_added="rank_genes_sub",
    )
    rgg = sub.uns["rank_genes_sub"]
    groups = list(rgg["names"].dtype.names)
    deg_dir = TAB / f"degs_{lineage}"
    deg_dir.mkdir(parents=True, exist_ok=True)
    for g in groups:
        df = pd.DataFrame({
            "gene": rgg["names"][g],
            "score": rgg["scores"][g],
            "logfoldchanges": rgg["logfoldchanges"][g],
            "pvals": rgg["pvals"][g],
            "pvals_adj": rgg["pvals_adj"][g],
        })
        df.to_csv(deg_dir / f"{g}_degs_{lineage}.csv", index=False)

    cc = sub.obs[leiden_key].value_counts().sort_index()
    pd.DataFrame({"cluster": cc.index.astype(str), "n_cells": cc.values}).to_csv(
        TAB / f"cluster_counts_{lineage}.csv", index=False
    )

    sc.settings.figdir = str(FIG)
    sc.set_figure_params(figsize=(4, 4), dpi=150)
    sc.pl.umap(sub, color=leiden_key, save=f"_{lineage}_leiden_res{str(params['res']).replace('.', 'p')}.pdf", show=False)
    sc.pl.umap(sub, color="status", save=f"_{lineage}_status.pdf", show=False)
    sc.pl.umap(sub, color="sample", save=f"_{lineage}_sample_diagnostic.pdf", show=False)
    print(f"[{lineage}] UMAP figures saved", flush=True)

    with open(TAB / f"run_summary_{lineage}.json", "w") as f:
        json.dump({
            "lineage": lineage,
            "params": params,
            "max_iter_harmony": MAX_ITER_HARMONY,
            "n_input_ids": int(len(ids)),
            "n_sub_cells": int(sub.n_obs),
            "n_sub_vars_hvg": int(sub.n_vars),
            "n_clusters": n_clusters,
            "hvg_n_top": N_TOP_HVG,
            "regress_keys": REGRESS_KEYS,
            "scale_max": SCALE_MAX,
            "normalize_target_sum": TARGET_SUM,
            "harmony_key": "sample",
            "random_seed": SEED,
            "rerun_reason": "harmony did not converge in initial run; re-run with max_iter_harmony=20",
            "input_consistent_h5ad": str(CONSISTENT),
            "input_qc_h5ad": str(QC),
        }, f, indent=2)
    print(f"[{lineage}] provenance saved", flush=True)


def main() -> None:
    t0 = time.time()
    print(f"[main] reading {CONSISTENT.name} (for cell IDs)…", flush=True)
    adata_c = sc.read_h5ad(CONSISTENT)
    print(f"[main] consistent shape: {adata_c.shape}, time={time.time()-t0:.1f}s", flush=True)

    ids_per_lineage = {}
    for lineage in LINEAGE_PARAMS:
        mask = adata_c.obs["leiden_coarse"].astype(str) == lineage
        ids = adata_c.obs_names[mask]
        ids_per_lineage[lineage] = ids
        print(f"[main] {lineage}: {len(ids)} cells", flush=True)

    print(f"\n[main] reading {QC.name} (for expression base)…", flush=True)
    adata_qc = sc.read_h5ad(QC)
    print(f"[main] qc shape: {adata_qc.shape}, time={time.time()-t0:.1f}s", flush=True)

    for lineage, params in LINEAGE_PARAMS.items():
        process_lineage(lineage, ids_per_lineage[lineage], params, adata_qc)
        gc.collect()

    print(f"\n[main] all 3 rerun done. total={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()