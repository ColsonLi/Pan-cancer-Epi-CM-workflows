"""Additional diagnostic figures for 07-score-rank-qc:
  - UMAP on consistent-filtered object (leiden_coarse / status / series)
  - UMAP on unfiltered (leiden_res0p4 / sample)
  - Cell-type composition stacked bar (leiden_coarse x status, x series)
"""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

SEED = 42
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
IN_UNFILT = ROOT / "epi-cm-core-workflow/h5ad/01-celltype_integration_clustering/07-score-rank-qc/adata_anno_score_genes_rank.h5ad"
IN_FILT = ROOT / "epi-cm-core-workflow/h5ad/01-celltype_integration_clustering/07-score-rank-qc/adata_anno_score_genes_rank_consistent.h5ad"
FIG = ROOT / "epi-cm-core-workflow/figures/01-celltype_integration_clustering/07-score-rank-qc"
TAB = ROOT / "epi-cm-core-workflow/tables/01-celltype_integration_clustering/07-score-rank-qc"


def main() -> None:
    t0 = time.time()
    sc.settings.figdir = str(FIG)
    sc.set_figure_params(figsize=(5, 5), dpi=150)

    # ---- 1. UMAP on consistent-filtered (additional views) ----
    print("[fig] loading filtered consistent…", flush=True)
    adata_f = sc.read_h5ad(IN_FILT)
    print(f"[fig] filtered shape={adata_f.shape}", flush=True)

    sc.pl.umap(adata_f, color="leiden_coarse", save="_filtered_leiden_coarse_umap.pdf", show=False)
    sc.pl.umap(adata_f, color="status", save="_filtered_status_umap.pdf", show=False)
    sc.pl.umap(adata_f, color="series", save="_filtered_series_umap.pdf", show=False)
    sc.pl.umap(adata_f, color="sample", save="_filtered_sample_diagnostic.pdf", show=False)
    print("[fig] filtered UMAPs saved", flush=True)

    # ---- 2. UMAP on unfiltered with raw leiden_res0p4 ----
    print("[fig] loading unfiltered…", flush=True)
    adata = sc.read_h5ad(IN_UNFILT)
    print(f"[fig] unfiltered shape={adata.shape}", flush=True)

    if "leiden_res0p4" in adata.obs:
        sc.pl.umap(adata, color="leiden_res0p4", save="_raw_leiden_res0p4_umap.pdf", show=False)
        print("[fig] raw leiden_res0p4 UMAP saved", flush=True)

    # ---- 3. Composition stacked bar ----
    # leiden_coarse x status
    ct_order = ["T", "B", "NK", "Mye", "Mast", "Epi", "Endo", "S"]
    avail = [c for c in ct_order if c in adata.obs["leiden_coarse"].unique()]
    other = [c for c in adata.obs["leiden_coarse"].unique() if c not in avail]
    final_order = avail + sorted(other)

    cross = pd.crosstab(adata.obs["leiden_coarse"], adata.obs["status"], normalize="index") * 100
    cross = cross.reindex(final_order)
    cross.to_csv(TAB / "composition_leiden_coarse_x_status_pct.csv")
    print("[fig] composition leiden_coarse x status:", flush=True)
    print(cross.round(2).to_string(), flush=True)

    fig, ax = plt.subplots(figsize=(7, 4))
    cross.plot(kind="bar", stacked=True, ax=ax, colormap="tab10")
    ax.set_ylabel("% of cells")
    ax.set_xlabel("leiden_coarse")
    ax.set_title("Cell type composition by status")
    ax.legend(title="status", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    out = FIG / "composition_leiden_coarse_x_status.pdf"
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"[fig] wrote {out.name}", flush=True)

    # leiden_coarse x series
    cross2 = pd.crosstab(adata.obs["leiden_coarse"], adata.obs["series"], normalize="index") * 100
    cross2 = cross2.reindex(final_order)
    cross2.to_csv(TAB / "composition_leiden_coarse_x_series_pct.csv")
    fig, ax = plt.subplots(figsize=(9, 4))
    cross2.plot(kind="bar", stacked=True, ax=ax, colormap="tab20")
    ax.set_ylabel("% of cells")
    ax.set_xlabel("leiden_coarse")
    ax.set_title("Cell type composition by series")
    ax.legend(title="series", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    out2 = FIG / "composition_leiden_coarse_x_series.pdf"
    plt.savefig(out2, bbox_inches="tight")
    plt.close()
    print(f"[fig] wrote {out2.name}", flush=True)

    print(f"[fig] done. total={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()