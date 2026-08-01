"""Subtype annotation (simple): cell_subtype = <leiden_coarse>_<first_DEG_gene>.

User instruction: annotation = leiden_coarse prefix + underscore + first DEG gene.
Then run a second DEG round grouped by cell_subtype.

Steps per lineage:
  1. Load adata_<lineage>_subclustered.h5ad (has leiden_sub + leiden_coarse columns)
  2. For each leiden_sub cluster, read its DEG CSV, get the first gene (rank 1)
  3. Build cell_subtype = f"{leiden_coarse}_{first_gene}"
  4. Save annotated h5ad
  5. Run DEG round 2: rank_genes_groups(groupby='cell_subtype', use_raw=True, n_genes=None)
  6. Save per-subtype DEG CSVs in degs_subtype_<lineage>/
"""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

SEED = 42
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
H5AD = ROOT / "epi-cm-core-workflow/h5ad/02-cell_subtype_integration_clustering/01-main"
TAB = ROOT / "epi-cm-core-workflow/tables/02-cell_subtype_integration_clustering/01-main"

LINEAGES = ["T", "B", "NK", "Mye", "Mast", "Epi", "Endo", "S"]


def process_lineage(lineage: str) -> None:
    in_h5ad = H5AD / f"adata_{lineage}_subclustered.h5ad"
    deg_dir = TAB / f"degs_{lineage}"
    print(f"\n========== {lineage} ==========", flush=True)
    if not in_h5ad.exists():
        print(f"[{lineage}] missing {in_h5ad.name}, skip", flush=True)
        return

    adata = sc.read_h5ad(in_h5ad)
    leiden_key = "leiden_sub"
    clusters = sorted(adata.obs[leiden_key].astype(str).unique().tolist(),
                      key=lambda x: int(x))
    print(f"[{lineage}] n_clusters: {len(clusters)}", flush=True)

    # 1. Build subtype map
    sub_map = {}
    for cl in clusters:
        fpath = deg_dir / f"{cl}_degs_{lineage}.csv"
        if not fpath.exists():
            print(f"  C{cl}: missing DEG, skip", flush=True)
            continue
        deg_df = pd.read_csv(fpath)
        first_gene = str(deg_df.iloc[0]["gene"])
        sub_type = f"{lineage}_{first_gene}"
        sub_map[cl] = sub_type
        print(f"  C{cl} -> {sub_type} (first DEG = {first_gene})", flush=True)

    # 2. Apply
    adata.obs["cell_subtype"] = adata.obs[leiden_key].astype(str).map(sub_map).astype("category")
    print(f"[{lineage}] cell_subtype counts:")
    print(adata.obs["cell_subtype"].value_counts().to_string(), flush=True)

    # 3. Save annotated h5ad
    adata.write_h5ad(in_h5ad, compression="gzip")
    print(f"[{lineage}] wrote {in_h5ad.name}", flush=True)

    # 4. Save cluster -> subtype map
    pd.DataFrame([
        {"cluster": cl, "first_deg_gene": pd.read_csv(deg_dir / f"{cl}_degs_{lineage}.csv").iloc[0]["gene"],
         "cell_subtype": sub_map[cl]}
        for cl in sub_map
    ]).to_csv(TAB / f"cluster_to_subtype_{lineage}.csv", index=False)
    print(f"[{lineage}] wrote cluster_to_subtype_{lineage}.csv", flush=True)

    # 5. DEG round 2: groupby=cell_subtype
    print(f"[{lineage}] DEG round 2 (groupby=cell_subtype)…", flush=True)
    assert adata.raw is not None
    sc.tl.rank_genes_groups(
        adata, groupby="cell_subtype", use_raw=True, method="t-test",
        n_genes=None, key_added="rank_genes_subtype",
    )
    rgg = adata.uns["rank_genes_subtype"]
    groups = list(rgg["names"].dtype.names)
    sub_deg_dir = TAB / f"degs_subtype_{lineage}"
    sub_deg_dir.mkdir(parents=True, exist_ok=True)
    for g in groups:
        df = pd.DataFrame({
            "gene": rgg["names"][g],
            "score": rgg["scores"][g],
            "logfoldchanges": rgg["logfoldchanges"][g],
            "pvals": rgg["pvals"][g],
            "pvals_adj": rgg["pvals_adj"][g],
        })
        # sanitize filename
        safe = g.replace("/", "_").replace(" ", "_")
        df.to_csv(sub_deg_dir / f"{safe}_degs_subtype_{lineage}.csv", index=False)
    print(f"[{lineage}] wrote {len(groups)} subtype DEG CSVs to {sub_deg_dir.name}", flush=True)

    # 6. Save again with subtype DEG stored
    adata.write_h5ad(in_h5ad, compression="gzip")


def main() -> None:
    t0 = time.time()
    for lineage in LINEAGES:
        process_lineage(lineage)
    print(f"\n[main] all 8 lineage annotation + DEG round 2 done. total={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()