"""Reorganize DEG CSVs per updated SKILL.md.

Per skill:
  Round 1: DEG on raw cluster column 'leiden_res0p4' (cluster numbers 0-13).
    Path: tables/degs_leiden_res0p4_pcs25_nn30_res0p4/
    File: {group}_degs_leiden_res0p4_pcs25_nn30_res0p4.csv
    Columns: gene, score, logfoldchanges, pvals, pvals_adj (+ extras preserved)
    Full DEG (no topXX), use_raw=True.

  Round 2: DEG on annotated column 'leiden_coarse' (broad cell type label).
    Path: tables/degs_leiden_coarse_pcs25_nn30_res0p4/
    File: {group}_degs_leiden_coarse_pcs25_nn30_res0p4.csv
    Same canonical columns.

Also: rename obs columns so 'leiden_res0p4' = raw cluster number, 'leiden_coarse'
= broad cell type, 'cell_type' = copy of leiden_coarse.
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
IN_H5AD = ROOT / "epi-cm-core-workflow/h5ad/01-celltype_integration_clustering/06-broad-annotation/adata_anno.h5ad"
TAB = ROOT / "epi-cm-core-workflow/tables/01-celltype_integration_clustering/06-broad-annotation"

PCS = 25
NN = 30
RES = "0p4"
RAW_COL = f"leiden_res{RES}"  # the raw cluster number column
COARSE_COL = "leiden_coarse"  # the broad annotation column (cell type)
CELLTYPE_COL = "cell_type"

DEG_DIR_RAW = TAB / f"degs_{RAW_COL}_pcs{PCS}_nn{NN}_res{RES}"
DEG_DIR_COARSE = TAB / f"degs_{COARSE_COL}_pcs{PCS}_nn{NN}_res{RES}"
DEG_DIR_RAW.mkdir(parents=True, exist_ok=True)
DEG_DIR_COARSE.mkdir(parents=True, exist_ok=True)


def export_full_degs(adata: sc.AnnData, groupby: str, out_dir: Path, suffix: str) -> None:
    """Run rank_genes_groups (full) and export one CSV per group."""
    key = f"rank_genes_{groupby}"
    print(f"\n[deg] === {groupby} ===", flush=True)
    print(f"[deg] rank_genes_groups (groupby={groupby}, use_raw=True, t-test, n_genes=None=full)…", flush=True)
    t1 = time.time()
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        use_raw=True,
        method="t-test",
        n_genes=None,
        key_added=key,
    )
    print(f"[deg] rank_genes_groups done in {time.time()-t1:.1f}s", flush=True)

    rgg = adata.uns[key]
    groups = list(rgg["names"].dtype.names)
    print(f"[deg] {len(groups)} groups in {groupby}", flush=True)

    for g in groups:
        n = len(rgg["names"][g])
        # canonical columns: gene, score, logfoldchanges, pvals, pvals_adj
        # rgg also has pts (if method supports it). We only have the 5 standard for t-test.
        df = pd.DataFrame({
            "gene": rgg["names"][g],
            "score": rgg["scores"][g],
            "logfoldchanges": rgg["logfoldchanges"][g],
            "pvals": rgg["pvals"][g],
            "pvals_adj": rgg["pvals_adj"][g],
        })
        # Sanitize group name for filename (only for the COARSE column; raw is "0".."13")
        safe = str(g).replace("/", "_").replace(" ", "_")
        fname = f"{safe}_degs_{suffix}_pcs{PCS}_nn{NN}_res{RES}.csv"
        out = out_dir / fname
        df.to_csv(out, index=False)
        print(f"  {g}: {n} genes -> {out.name}", flush=True)


def main() -> None:
    t0 = time.time()
    print(f"[deg] reading {IN_H5AD.name}…", flush=True)
    adata = sc.read_h5ad(IN_H5AD)
    print(f"[deg] loaded {adata.shape} in {time.time()-t0:.1f}s", flush=True)

    # Step 1: ensure column layout matches skill
    # The current adata has:
    #   obs['leiden_coarse'] = cluster number (e.g., "0", "1", ...)
    #   obs['cell_type'] = annotation (e.g., "T", "B", ...)
    # Rename to: leiden_res0p4 = raw cluster, leiden_coarse = broad cell type
    if "leiden_coarse" in adata.obs and RAW_COL not in adata.obs:
        adata.obs[RAW_COL] = adata.obs["leiden_coarse"].astype(str).astype("category")
        print(f"[deg] copied 'leiden_coarse' -> '{RAW_COL}' as raw cluster column", flush=True)
    if CELLTYPE_COL in adata.obs and COARSE_COL not in adata.obs:
        adata.obs[COARSE_COL] = adata.obs[CELLTYPE_COL].astype(str).astype("category")
        print(f"[deg] copied 'cell_type' -> '{COARSE_COL}' as broad annotation", flush=True)
    elif COARSE_COL not in adata.obs:
        # Should not happen, but safety
        adata.obs[COARSE_COL] = adata.obs[CELLTYPE_COL].astype(str).astype("category")

    print(f"[deg] adata.obs[leiden_res0p4] n_cats: {adata.obs[RAW_COL].nunique()}", flush=True)
    print(f"[deg] adata.obs[leiden_coarse] n_cats: {adata.obs[COARSE_COL].nunique()}", flush=True)
    print(f"[deg] adata.obs[leiden_coarse] value counts:\n{adata.obs[COARSE_COL].value_counts()}", flush=True)

    assert adata.raw is not None, "adata.raw missing — required for use_raw=True DEG"

    # Step 2: full DEG per group for both columns
    export_full_degs(adata, RAW_COL, DEG_DIR_RAW, suffix=RAW_COL)
    gc.collect()
    export_full_degs(adata, COARSE_COL, DEG_DIR_COARSE, suffix=COARSE_COL)

    # Save updated h5ad
    print(f"\n[deg] writing updated {IN_H5AD.name}…", flush=True)
    adata.write_h5ad(IN_H5AD, compression="gzip")
    print(f"[deg] done. total={time.time()-t0:.1f}s", flush=True)

    # Provenance
    with open(TAB / "degs_export_readme.txt", "w") as f:
        f.write("DEG export (per updated SKILL.md)\n\n")
        f.write("Round 1 (raw cluster DEGs):\n")
        f.write(f"  column: {RAW_COL}\n")
        f.write(f"  dir: {DEG_DIR_RAW}\n")
        f.write(f"  method: t-test, use_raw=True, n_genes=None (full)\n\n")
        f.write("Round 2 (broad annotation DEGs):\n")
        f.write(f"  column: {COARSE_COL}\n")
        f.write(f"  dir: {DEG_DIR_COARSE}\n")
        f.write(f"  method: t-test, use_raw=True, n_genes=None (full)\n\n")
        f.write("Per-group CSV columns: gene, score, logfoldchanges, pvals, pvals_adj\n")
        f.write(f"params: pcs={PCS}, nn={NN}, res=0.4\n")
        f.write(f"random_seed: {SEED}\n")


if __name__ == "__main__":
    main()