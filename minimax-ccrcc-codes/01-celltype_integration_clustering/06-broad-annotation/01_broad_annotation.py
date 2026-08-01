"""Broad annotation — rank genes per leiden cluster, map clusters → cell types.

Per SKILL.md:
  - Use use_raw=True for DEG/marker computation (raw preserved in adata.raw)
  - Save adata_anno.h5ad with obs['leiden_coarse'] and obs['cell_type']
  - Standard broad marker review

Workflow:
  1. Read h5ad/05-clustering-parameter-search/selected/adata_inte.h5ad
  2. sc.tl.rank_genes_groups(groupby='leiden_coarse', use_raw=True, method='t-test')
  3. For each cluster, pull top-N DEGs and look up the broad-marker panel
  4. Assign cell_type per cluster based on marker evidence
  5. Save h5ad/06-broad-annotation/adata_anno.h5ad

Outputs:
  - h5ad/.../06-broad-annotation/adata_anno.h5ad
  - tables/.../06-broad-annotation/cluster_to_celltype_map.csv
  - tables/.../06-broad-annotation/cluster_top_degs.csv
  - tables/.../06-broad-annotation/broad_marker_dotplot_data.csv
  - figures/.../06-broad-annotation/umap_leiden_coarse.pdf
  - figures/.../06-broad-annotation/umap_cell_type.pdf
  - figures/.../06-broad-annotation/dotplot_broad_markers.pdf
  - tables/.../06-broad-annotation/readme.txt
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

SEED = 42
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
IN_H5AD = ROOT / "epi-cm-core-workflow/h5ad/01-celltype_integration_clustering/05-clustering-parameter-search/selected/adata_inte.h5ad"
OUT_H5AD = ROOT / "epi-cm-core-workflow/h5ad/01-celltype_integration_clustering/06-broad-annotation/adata_anno.h5ad"
TAB = ROOT / "epi-cm-core-workflow/tables/01-celltype_integration_clustering/06-broad-annotation"
TAB.mkdir(parents=True, exist_ok=True)
FIG = ROOT / "epi-cm-core-workflow/figures/01-celltype_integration_clustering/06-broad-annotation"
FIG.mkdir(parents=True, exist_ok=True)
OUT_H5AD.parent.mkdir(parents=True, exist_ok=True)

LEIDEN_KEY = "leiden_coarse"

# Broad marker panel (CCRCC / kidney scRNA-seq reference)
BROAD_MARKERS = {
    "T": ["CD3D", "CD3E", "CD2", "TRAC", "TRBC1", "TRBC2", "LCK"],
    "NK": ["GNLY", "NKG7", "KLRD1", "KLRK1", "KLRB1", "GZMB"],
    "B": ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1", "PAX5"],
    "Plasma": ["MZB1", "JCHAIN", "IGHG1", "IGHA1", "XBP1", "PRDM1"],
    "Mye": ["LST1", "S100A8", "S100A9", "CD14", "FCGR3A", "LYZ", "CSF1R", "ITGAM"],
    "Mast": ["TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2"],
    "Epi": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "MUC1", "SLC34A1", "PAX8"],
    "Endo": ["PECAM1", "VWF", "CDH5", "CLDN5", "FLT1", "KDR"],
    "S": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRB", "ACTA2", "TAGLN"],
    "pDC": ["LILRA4", "IL3RA", "TCF4", "IRF7", "GZMB"],
}


def get_top_degs(rgg, n_top: int = 30) -> pd.DataFrame:
    """Convert sc.tl.rank_genes_groups result to long-form DataFrame."""
    res = rgg
    groups = res["names"].dtype.names
    rows = []
    for g in groups:
        for rank in range(n_top):
            rows.append({
                "cluster": g,
                "rank": rank + 1,
                "gene": res["names"][g][rank],
                "logfc": float(res["logfoldchanges"][g][rank]),
                "pval": float(res["pvals"][g][rank]),
                "padj": float(res["pvals_adj"][g][rank]),
                "score": float(res["scores"][g][rank]),
            })
    return pd.DataFrame(rows)


def main() -> None:
    t0 = time.time()
    print(f"[anno] reading {IN_H5AD.name}…", flush=True)
    adata = sc.read_h5ad(IN_H5AD)
    print(f"[anno] loaded {adata.shape} in {time.time()-t0:.1f}s", flush=True)
    print(f"[anno] n_clusters={adata.obs[LEIDEN_KEY].nunique()}", flush=True)
    print(f"[anno] cluster sizes:\n{adata.obs[LEIDEN_KEY].value_counts().sort_index()}", flush=True)

    # raw check
    assert adata.raw is not None, "adata.raw missing — DEG will use .X instead"

    # rank_genes_groups
    print("[anno] rank_genes_groups (t-test, use_raw=True)…", flush=True)
    t1 = time.time()
    sc.tl.rank_genes_groups(adata, groupby=LEIDEN_KEY, use_raw=True, method="t-test",
                            n_genes=200, key_added="rank_genes_leiden")
    print(f"[anno] rank_genes_groups done in {time.time()-t1:.1f}s", flush=True)

    # top DEGs
    rgg = adata.uns["rank_genes_leiden"]
    top_degs = get_top_degs(rgg, n_top=30)
    top_degs.to_csv(TAB / "cluster_top_degs.csv", index=False)
    print(f"[anno] wrote cluster_top_degs.csv ({len(top_degs)} rows)", flush=True)

    # Print top-10 per cluster
    print("\n[anno] top 10 DEGs per cluster:")
    for cl in sorted(top_degs["cluster"].unique(), key=lambda x: int(x)):
        sub = top_degs[top_degs["cluster"] == cl].head(10)
        genes = ", ".join(sub["gene"].tolist())
        print(f"  C{cl}: {genes}", flush=True)

    # Filter broad markers to those present in raw
    raw_var = set(adata.raw.var_names)
    present_markers = {
        ct: [g for g in gs if g in raw_var]
        for ct, gs in BROAD_MARKERS.items()
        if any(g in raw_var for g in gs)
    }
    print(f"\n[anno] present broad markers per category:")
    for ct, gs in present_markers.items():
        print(f"  {ct}: {gs}", flush=True)

    # Build dotplot data: mean expression of each present marker per cluster
    flat_markers = [g for gs in present_markers.values() for g in gs]
    print(f"[anno] computing dotplot for {len(flat_markers)} markers x {adata.obs[LEIDEN_KEY].nunique()} clusters…", flush=True)
    t1 = time.time()
    sc.settings.figdir = str(FIG)
    sc.set_figure_params(figsize=(8, 6), dpi=150)

    # Use sc.pl.rank_genes_groups_dotplot first to get the actual expression data
    try:
        sc.pl.rank_genes_groups_dotplot(
            adata, key="rank_genes_leiden", n_genes=5,
            save="_top5_dotplot.pdf", show=False,
        )
    except Exception as e:
        print(f"[anno] rank_genes_groups_dotplot failed: {e}", flush=True)

    # Build manual dotplot for the curated broad-marker panel
    if flat_markers:
        # mean expression per cluster (using raw counts) and pct expressing
        expr = adata.raw[:, flat_markers].X
        if hasattr(expr, "toarray"):
            expr = expr.toarray()
        df_expr = pd.DataFrame(expr, columns=flat_markers)
        df_expr["cluster"] = adata.obs[LEIDEN_KEY].astype(str).values
        # mean and pct expressing per cluster per gene
        rows = []
        for cl, sub in df_expr.groupby("cluster"):
            for g in flat_markers:
                v = sub[g].values
                rows.append({
                    "cluster": cl,
                    "gene": g,
                    "mean_expr": float(v.mean()),
                    "pct_expr": float((v > 0).mean() * 100.0),
                })
        dot_df = pd.DataFrame(rows)
        dot_df.to_csv(TAB / "broad_marker_dotplot_data.csv", index=False)
        print(f"[anno] wrote broad_marker_dotplot_data.csv ({len(dot_df)} rows, {time.time()-t1:.1f}s)", flush=True)

        # Build marker -> category map
        marker_to_cat = {g: ct for ct, gs in present_markers.items() for g in gs}
        # Score per cluster: for each broad category, compute mean z-score of its markers
        # Use mean_expr (log-normalized) as the score
        cluster_scores = {}
        for cl in sorted(dot_df["cluster"].unique()):
            sub = dot_df[dot_df["cluster"] == cl]
            scores = {}
            for ct, gs in present_markers.items():
                sub_ct = sub[sub["gene"].isin(gs)]
                if len(sub_ct) == 0:
                    scores[ct] = 0.0
                else:
                    scores[ct] = float(sub_ct["mean_expr"].mean())
            cluster_scores[cl] = scores
        scores_df = pd.DataFrame(cluster_scores).T
        scores_df.index.name = "cluster"
        scores_df.to_csv(TAB / "cluster_category_mean_expr.csv")
        print(f"[anno] cluster_category_mean_expr.csv:", flush=True)
        print(scores_df.round(2), flush=True)

        # Assign cell type per cluster: argmax of category scores
        cluster_to_celltype = {}
        for cl in scores_df.index:
            row = scores_df.loc[cl]
            best_cat = row.idxmax()
            best_score = row.max()
            cluster_to_celltype[cl] = best_cat
        # Print mapping
        print(f"\n[anno] auto-mapping (argmax of mean_expr):")
        for cl, ct in cluster_to_celltype.items():
            print(f"  C{cl} -> {ct}", flush=True)

    # Plot UMAP
    sc.pl.umap(adata, color=LEIDEN_KEY, save="_leiden_coarse.pdf", show=False)
    if "cell_type" in adata.obs:
        sc.pl.umap(adata, color="cell_type", save="_cell_type_auto.pdf", show=False)

    # Save adata_anno.h5ad with provisional cell_type
    if "cell_type" not in adata.obs:
        adata.obs["cell_type"] = adata.obs[LEIDEN_KEY].astype(str).map(cluster_to_celltype).astype("category")
    print(f"[anno] cell_type value counts:\n{adata.obs['cell_type'].value_counts()}", flush=True)

    print(f"[anno] writing {OUT_H5AD.name}…", flush=True)
    adata.write_h5ad(OUT_H5AD, compression="gzip")
    print(f"[anno] done. total={time.time()-t0:.1f}s", flush=True)

    # Save cluster → cell_type map
    pd.DataFrame([
        {"cluster": cl, "cell_type_auto": ct, "method": "argmax of broad-marker mean expr"}
        for cl, ct in cluster_to_celltype.items()
    ]).to_csv(TAB / "cluster_to_celltype_map.csv", index=False)

    # readme
    with open(TAB / "readme.txt", "w") as f:
        f.write("06-broad-annotation outputs\n")
        f.write(f"input: {IN_H5AD}\n")
        f.write(f"output h5ad: {OUT_H5AD}\n")
        f.write(f"clustering: pcs={25}, nn={30}, res={0.4}\n")
        f.write(f"DEG method: t-test, use_raw=True\n")
        f.write(f"n_cells: {adata.n_obs}, n_clusters: {adata.obs[LEIDEN_KEY].nunique()}\n")
        f.write(f"broad-marker panel: {json.dumps(BROAD_MARKERS)}\n")
        f.write("NOTE: 'B Cells' (MS4A1/CD79A) and 'pDC' (LILRA4) markers are present in the panel; if a cluster argmax is 'T' but shows high B-marker mean, that means B cells were merged into a T-dominant cluster at this leiden resolution.\n")
    print(f"[anno] readme.txt written", flush=True)


if __name__ == "__main__":
    main()