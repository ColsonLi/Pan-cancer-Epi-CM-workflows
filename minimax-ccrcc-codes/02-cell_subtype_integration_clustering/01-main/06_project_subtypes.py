"""06-ccrcc-subtypes-to-full-adata.

For each cell in the consistent filtered atlas, look up its lineage's
cell_subtype and project it back. Build cell_type from cell_subtype prefix.

Outputs:
  h5ad/02-cell_subtype_integration_clustering/06-ccrcc-subtypes-to-full-adata/
    adata_anno_cellsubtype.h5ad
  tables/02-cell_subtype_integration_clustering/06-ccrcc-subtypes-to-full-adata/
    celltype_cellsubtype_counts.csv
    unmatched_cells.csv (if any)
  figures/02-cell_subtype_integration_clustering/06-ccrcc-subtypes-to-full-adata/
    umap_cell_type.pdf
    umap_cell_subtype.pdf
    umap_leiden_coarse_vs_cell_type.pdf
"""
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
ATLAS = ROOT / "epi-cm-core-workflow/h5ad/01-celltype_integration_clustering/07-score-rank-qc/adata_anno_score_genes_rank_consistent.h5ad"
H5AD_IN = ROOT / "epi-cm-core-workflow/h5ad/02-cell_subtype_integration_clustering/01-main"
H5AD_OUT = ROOT / "epi-cm-core-workflow/h5ad/02-cell_subtype_integration_clustering/06-ccrcc-subtypes-to-full-adata"
TAB = ROOT / "epi-cm-core-workflow/tables/02-cell_subtype_integration_clustering/06-ccrcc-subtypes-to-full-adata"
FIG = ROOT / "epi-cm-core-workflow/figures/02-cell_subtype_integration_clustering/06-ccrcc-subtypes-to-full-adata"

LINEAGES = ["T", "B", "NK", "Mye", "Mast", "Epi", "Endo", "S"]
H5AD_OUT.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)


def main() -> None:
    t0 = time.time()
    print(f"[proj] reading atlas {ATLAS.name}…", flush=True)
    atlas = sc.read_h5ad(ATLAS)
    print(f"[proj] atlas shape: {atlas.shape}, time={time.time()-t0:.1f}s", flush=True)
    print(f"[proj] leiden_coarse cats: {atlas.obs['leiden_coarse'].nunique()}", flush=True)

    # Build per-cell cell_subtype lookup from each lineage h5ad
    cell_subtype_map = {}
    cell_type_overrides = {}  # cell -> override cell_type
    for lineage in LINEAGES:
        fpath = H5AD_IN / f"adata_{lineage}_subclustered.h5ad"
        if not fpath.exists():
            print(f"[proj] WARN: {fpath.name} missing, skip", flush=True)
            continue
        sub = sc.read_h5ad(fpath, backed="r")
        if "cell_subtype" not in sub.obs:
            print(f"[proj] WARN: cell_subtype missing in {fpath.name}, skip", flush=True)
            continue
        # cell_subtype index
        cs = sub.obs["cell_subtype"].astype(str)
        overlap = set(sub.obs_names) & set(atlas.obs_names)
        print(f"[proj] {lineage}: {len(sub.obs_names)} sub cells, {len(overlap)} in atlas", flush=True)
        for cid, csub in zip(sub.obs_names, cs):
            if cid in atlas.obs_names:
                cell_subtype_map[cid] = csub
                # cell_type override = lineage prefix (broad category)
                # In our naming, cell_subtype is like "T_CD3D", "NK_NKG7"
                # The lineage prefix already equals the cell_type label
                cell_type_overrides[cid] = lineage

    print(f"[proj] mapped {len(cell_subtype_map)} / {atlas.n_obs} cells", flush=True)

    # Project back to atlas
    atlas.obs["cell_subtype"] = pd.Categorical(
        [cell_subtype_map.get(c, "Unknown") for c in atlas.obs_names],
        categories=sorted(set(cell_subtype_map.values()) | {"Unknown"}),
    )
    # cell_type = lineage prefix (if mapped) or leiden_coarse (if unknown)
    atlas.obs["cell_type"] = pd.Categorical(
        [cell_type_overrides.get(c, str(atlas.obs.loc[c, "leiden_coarse"])) for c in atlas.obs_names]
    )
    print(f"[proj] cell_subtype counts (top 10):")
    print(atlas.obs["cell_subtype"].value_counts().head(10).to_string(), flush=True)
    print(f"[proj] Unknown cells: {(atlas.obs['cell_subtype'] == 'Unknown').sum()}", flush=True)
    print(f"[proj] cell_type counts:")
    print(atlas.obs["cell_type"].value_counts().to_string(), flush=True)

    # Save
    out_h5ad = H5AD_OUT / "adata_anno_cellsubtype.h5ad"
    print(f"[proj] writing {out_h5ad.name}…", flush=True)
    atlas.write_h5ad(out_h5ad, compression="gzip")
    print(f"[proj] wrote {out_h5ad.name}", flush=True)

    # Tables
    counts = pd.crosstab(atlas.obs["cell_type"], atlas.obs["cell_subtype"], dropna=False)
    counts.to_csv(TAB / "celltype_cellsubtype_counts.csv")
    print(f"[proj] wrote celltype_cellsubtype_counts.csv", flush=True)

    unmatched = atlas.obs[atlas.obs["cell_subtype"] == "Unknown"]
    if len(unmatched) > 0:
        unmatched.to_csv(TAB / "unmatched_cells.csv")
        print(f"[proj] {len(unmatched)} unmatched cells -> unmatched_cells.csv", flush=True)

    pd.crosstab(atlas.obs["cell_type"], atlas.obs["status"]).to_csv(TAB / "celltype_status_counts.csv")
    pd.crosstab(atlas.obs["cell_type"], atlas.obs["series"]).to_csv(TAB / "celltype_series_counts.csv")
    print(f"[proj] wrote celltype_status/series counts", flush=True)

    # Figures
    sc.settings.figdir = str(FIG)
    sc.set_figure_params(figsize=(6, 6), dpi=150)
    sc.pl.umap(atlas, color="cell_type", save="_cell_type.pdf", show=False)
    sc.pl.umap(atlas, color="cell_subtype", save="_cell_subtype.pdf", show=False)
    sc.pl.umap(atlas, color=["leiden_coarse", "cell_type"], ncols=2, wspace=0.3, save="_leiden_coarse_vs_projected.pdf", show=False)
    print(f"[proj] UMAP figures saved", flush=True)

    # Summary
    summary = {
        "n_atlas_cells": int(atlas.n_obs),
        "n_mapped_cells": int(len(cell_subtype_map)),
        "n_unmatched_cells": int((atlas.obs["cell_subtype"] == "Unknown").sum()),
        "n_lineages_used": len(LINEAGES),
        "lineages": LINEAGES,
        "input_atlas": str(ATLAS),
        "input_subclustered_dir": str(H5AD_IN),
    }
    with open(TAB / "projection_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[proj] wrote projection_summary.json", flush=True)
    print(f"[proj] done. total={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()