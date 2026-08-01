"""Re-run figure plotting for 07-score-rank-qc (UMAP + dotplot)."""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

SEED = 42
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
IN_H5AD = ROOT / "epi-cm-core-workflow/h5ad/01-celltype_integration_clustering/07-score-rank-qc/adata_anno_score_genes_rank.h5ad"
FIG = ROOT / "epi-cm-core-workflow/figures/01-celltype_integration_clustering/07-score-rank-qc"
FIG.mkdir(parents=True, exist_ok=True)
TAB = ROOT / "epi-cm-core-workflow/tables/01-celltype_integration_clustering/07-score-rank-qc"
DEG_DIR = ROOT / "epi-cm-core-workflow/tables/01-celltype_integration_clustering/06-broad-annotation/degs_leiden_coarse_pcs25_nn30_res0p4"


def sanitize(name: str) -> str:
    return name.replace(" ", "_").replace("/", "_").replace("-", "_")


def main() -> None:
    t0 = time.time()
    print(f"[fig] reading {IN_H5AD.name}…", flush=True)
    adata = sc.read_h5ad(IN_H5AD)
    print(f"[fig] loaded {adata.shape} in {time.time()-t0:.1f}s", flush=True)

    sc.settings.figdir = str(FIG)
    sc.set_figure_params(figsize=(5, 5), dpi=150)

    print("[fig] UMAP leiden_coarse…", flush=True)
    sc.pl.umap(adata, color="leiden_coarse", save="_major_lineage_umap.pdf", show=False)
    print("[fig] UMAP best_rank_type_global…", flush=True)
    sc.pl.umap(adata, color="best_rank_type_global", save="_best_rank_type_global_umap.pdf", show=False)
    print("[fig] UMAP rank_consistent…", flush=True)
    sc.pl.umap(adata, color="rank_consistent", save="_rank_consistent_umap.pdf", show=False)

    # canonical dotplot: top 3 per leiden_coarse label
    labels = sorted(adata.obs["leiden_coarse"].astype(str).unique().tolist())
    marker_dict = {}
    for lab in labels:
        safe = sanitize(lab)
        fpath = DEG_DIR / f"{safe}_degs_leiden_coarse_pcs25_nn30_res0p4.csv"
        if not fpath.exists():
            print(f"[fig] missing DEG for {lab}: {fpath}", flush=True)
            continue
        df = pd.read_csv(fpath)
        genes = []
        for g in df["gene"].tolist():
            if g in adata.raw.var_names:
                genes.append(g)
            if len(genes) >= 3:
                break
        marker_dict[lab] = genes
    print(f"[fig] marker dict:", flush=True)
    for k, v in marker_dict.items():
        print(f"  {k}: {v}", flush=True)

    if marker_dict:
        # Build flat var_names for older scanpy
        flat = []
        for lab in labels:
            if lab in marker_dict:
                flat.extend(marker_dict[lab])
        print(f"[fig] flat markers: {flat}", flush=True)
        if flat:
            sc.pl.dotplot(
                adata,
                var_names=flat,
                groupby="leiden_coarse",
                standard_scale="var",
                use_raw=True,
                save="_all_leiden_coarse.pdf",
                show=False,
            )
            print(f"[fig] dotplot saved", flush=True)

    # Save summary
    summary = {
        "n_cells_input": int(adata.n_obs),
        "n_obs_score": int(adata.n_obs),
        "n_obs_consistent": int(adata.obs["rank_consistent"].sum()),
        "n_labels": len(labels),
        "labels": labels,
        "marker_dict": marker_dict,
    }
    with open(TAB / "score_rank_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[fig] done. total={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()