"""Fix: if multiple clusters in a lineage share the same top1 DEG gene,
the 2nd+ occurrence uses top2 instead of top1.

Per user: avoid duplicate cell_subtype names within a lineage.
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
    map_csv = TAB / f"cluster_to_subtype_{lineage}.csv"
    sub_deg_dir = TAB / f"degs_subtype_{lineage}"
    print(f"\n========== {lineage} ==========", flush=True)
    if not map_csv.exists():
        print(f"[{lineage}] no map csv, skip", flush=True)
        return

    m = pd.read_csv(map_csv)
    used_genes = set()
    new_map = {}
    changes = []
    for _, row in m.iterrows():
        cl = str(row["cluster"])
        deg_csv = deg_dir / f"{cl}_degs_{lineage}.csv"
        if not deg_csv.exists():
            new_map[cl] = row["cell_subtype"]
            continue
        deg = pd.read_csv(deg_csv)
        # find first gene not yet used
        chosen_gene = None
        for g in deg["gene"].tolist():
            if g not in used_genes:
                chosen_gene = str(g)
                break
        if chosen_gene is None:
            chosen_gene = str(deg.iloc[0]["gene"])  # fallback
        old_gene = str(row["first_deg_gene"])
        if chosen_gene != old_gene:
            changes.append((cl, old_gene, chosen_gene))
        new_map[cl] = f"{lineage}_{chosen_gene}"
        used_genes.add(chosen_gene)

    if changes:
        print(f"[{lineage}] changes:")
        for cl, old, new in changes:
            print(f"  C{cl}: {old} -> {new}", flush=True)
    else:
        print(f"[{lineage}] no changes", flush=True)
        return

    # Update map csv
    new_rows = []
    for _, row in m.iterrows():
        cl = str(row["cluster"])
        new_rows.append({
            "cluster": cl,
            "first_deg_gene": new_map[cl].split("_", 1)[1],
            "cell_subtype": new_map[cl],
        })
    pd.DataFrame(new_rows).to_csv(map_csv, index=False)
    print(f"[{lineage}] updated {map_csv.name}", flush=True)

    # Update h5ad
    adata = sc.read_h5ad(in_h5ad)
    adata.obs["cell_subtype"] = adata.obs["leiden_sub"].astype(str).map(new_map).astype("category")
    print(f"[{lineage}] cell_subtype counts:\n{adata.obs['cell_subtype'].value_counts().to_string()}", flush=True)

    # DEG round 2 (re-run since subtypes changed)
    print(f"[{lineage}] DEG round 2 (groupby=cell_subtype, n_genes=None)…", flush=True)
    assert adata.raw is not None
    sc.tl.rank_genes_groups(
        adata, groupby="cell_subtype", use_raw=True, method="t-test",
        n_genes=None, key_added="rank_genes_subtype",
    )
    rgg = adata.uns["rank_genes_subtype"]
    groups = list(rgg["names"].dtype.names)
    sub_deg_dir.mkdir(parents=True, exist_ok=True)
    # remove old CSVs
    for f in sub_deg_dir.glob(f"*_degs_subtype_{lineage}.csv"):
        f.unlink()
    for g in groups:
        df = pd.DataFrame({
            "gene": rgg["names"][g],
            "score": rgg["scores"][g],
            "logfoldchanges": rgg["logfoldchanges"][g],
            "pvals": rgg["pvals"][g],
            "pvals_adj": rgg["pvals_adj"][g],
        })
        safe = g.replace("/", "_").replace(" ", "_")
        df.to_csv(sub_deg_dir / f"{safe}_degs_subtype_{lineage}.csv", index=False)
    print(f"[{lineage}] wrote {len(groups)} subtype DEG CSVs to {sub_deg_dir.name}", flush=True)

    adata.write_h5ad(in_h5ad, compression="gzip")
    print(f"[{lineage}] wrote {in_h5ad.name}", flush=True)


def main() -> None:
    t0 = time.time()
    for lin in LINEAGES:
        process_lineage(lin)
    print(f"\n[main] done. total={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()