"""Redo Epi: re-run subclustering + annotation + DEG round 2.

The Epi h5ad was corrupted when the previous script was interrupted.
Re-run from clean adata_qc with max_iter_harmony=20.
"""
from __future__ import annotations

import json
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

LINEAGE = "Epi"
PARAMS = dict(pcs=15, nn=30, res=0.6)
N_TOP_HVG = 3000
TARGET_SUM = 1e4
SCALE_MAX = 10
REGRESS_KEYS = ["total_counts", "pct_counts_mt"]
MAX_ITER_HARMONY = 20


def main() -> None:
    t0 = time.time()
    print(f"[epi] reading consistent (for IDs)…", flush=True)
    adata_c = sc.read_h5ad(CONSISTENT)
    mask = adata_c.obs["leiden_coarse"].astype(str) == LINEAGE
    ids = adata_c.obs_names[mask]
    print(f"[epi] n_ids: {len(ids)}, time={time.time()-t0:.1f}s", flush=True)

    print(f"[epi] reading adata_qc (expression base)…", flush=True)
    adata_qc = sc.read_h5ad(QC)
    print(f"[epi] qc shape: {adata_qc.shape}, time={time.time()-t0:.1f}s", flush=True)

    sub = adata_qc[ids].copy()
    print(f"[epi] sub shape: {sub.shape}", flush=True)

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
    rsc.pp.neighbors(sub, n_neighbors=PARAMS["nn"], n_pcs=PARAMS["pcs"], use_rep="X_pca_inte", random_state=SEED)
    leiden_key = "leiden_sub"
    rsc.tl.leiden(sub, resolution=PARAMS["res"], key_added=leiden_key, random_state=SEED)
    n_clusters = int(sub.obs[leiden_key].nunique())
    print(f"[epi] leiden n_clusters: {n_clusters}", flush=True)
    rsc.tl.umap(sub, random_state=SEED)
    rsc.get.anndata_to_CPU(sub)
    gc.collect()

    # DEG round 1 (cluster)
    sc.tl.rank_genes_groups(
        sub, groupby=leiden_key, use_raw=True, method="t-test",
        n_genes=None, key_added="rank_genes_sub",
    )
    rgg = sub.uns["rank_genes_sub"]
    groups = list(rgg["names"].dtype.names)
    deg_dir = TAB / f"degs_{LINEAGE}"
    deg_dir.mkdir(parents=True, exist_ok=True)
    sub_map = {}
    for g in groups:
        df = pd.DataFrame({
            "gene": rgg["names"][g],
            "score": rgg["scores"][g],
            "logfoldchanges": rgg["logfoldchanges"][g],
            "pvals": rgg["pvals"][g],
            "pvals_adj": rgg["pvals_adj"][g],
        })
        df.to_csv(deg_dir / f"{g}_degs_{LINEAGE}.csv", index=False)
        sub_map[g] = f"{LINEAGE}_{rgg['names'][g][0]}"
        print(f"  C{g} -> {sub_map[g]} (first DEG = {rgg['names'][g][0]})", flush=True)

    # cluster counts
    cc = sub.obs[leiden_key].value_counts().sort_index()
    pd.DataFrame({"cluster": cc.index.astype(str), "n_cells": cc.values}).to_csv(
        TAB / f"cluster_counts_{LINEAGE}.csv", index=False
    )

    # cell_subtype
    sub.obs["cell_subtype"] = sub.obs[leiden_key].astype(str).map(sub_map).astype("category")
    print(f"[epi] cell_subtype counts:\n{sub.obs['cell_subtype'].value_counts().to_string()}", flush=True)

    # Save cluster_to_subtype
    pd.DataFrame([
        {"cluster": cl, "first_deg_gene": rgg['names'][cl][0], "cell_subtype": sub_map[cl]}
        for cl in sub_map
    ]).to_csv(TAB / f"cluster_to_subtype_{LINEAGE}.csv", index=False)

    # UMAP figures
    sc.settings.figdir = str(FIG)
    sc.set_figure_params(figsize=(4, 4), dpi=150)
    sc.pl.umap(sub, color=leiden_key, save=f"_{LINEAGE}_leiden_res0p6.pdf", show=False)
    sc.pl.umap(sub, color="status", save=f"_{LINEAGE}_status.pdf", show=False)
    print(f"[epi] UMAP figures saved", flush=True)

    # DEG round 2 (cell_subtype)
    sc.tl.rank_genes_groups(
        sub, groupby="cell_subtype", use_raw=True, method="t-test",
        n_genes=None, key_added="rank_genes_subtype",
    )
    rgg2 = sub.uns["rank_genes_subtype"]
    sub2_groups = list(rgg2["names"].dtype.names)
    sub_deg_dir = TAB / f"degs_subtype_{LINEAGE}"
    sub_deg_dir.mkdir(parents=True, exist_ok=True)
    for g in sub2_groups:
        df = pd.DataFrame({
            "gene": rgg2["names"][g],
            "score": rgg2["scores"][g],
            "logfoldchanges": rgg2["logfoldchanges"][g],
            "pvals": rgg2["pvals"][g],
            "pvals_adj": rgg2["pvals_adj"][g],
        })
        safe = g.replace("/", "_").replace(" ", "_")
        df.to_csv(sub_deg_dir / f"{safe}_degs_subtype_{LINEAGE}.csv", index=False)
    print(f"[epi] wrote {len(sub2_groups)} subtype DEG CSVs to {sub_deg_dir.name}", flush=True)

    # Save h5ad
    out_h5ad = H5AD / f"adata_{LINEAGE}_subclustered.h5ad"
    sub.write_h5ad(out_h5ad, compression="gzip")
    print(f"[epi] wrote {out_h5ad.name}", flush=True)

    with open(TAB / f"run_summary_{LINEAGE}.json", "w") as f:
        json.dump({
            "lineage": LINEAGE, "params": PARAMS, "max_iter_harmony": MAX_ITER_HARMONY,
            "n_input_ids": int(len(ids)), "n_sub_cells": int(sub.n_obs),
            "n_sub_vars_hvg": int(sub.n_vars), "n_clusters": n_clusters,
            "n_subtypes": len(sub2_groups), "hvg_n_top": N_TOP_HVG,
            "regress_keys": REGRESS_KEYS, "scale_max": SCALE_MAX,
            "normalize_target_sum": TARGET_SUM, "harmony_key": "sample",
            "random_seed": SEED, "rerun_reason": "Epi h5ad was corrupted during previous annotation run; full redo",
            "input_consistent_h5ad": str(CONSISTENT), "input_qc_h5ad": str(QC),
        }, f, indent=2)
    print(f"[epi] done. total={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()