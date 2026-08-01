"""07-score-rank-qc: build score_genes + best_rank_type_global per skill.

Per updated SKILL.md:
  1. Read saved full DEG CSVs from tables/degs_leiden_coarse_pcs25_nn30_res0p4/
  2. For each observed leiden_coarse label, take first 100 usable genes
     (restricted to genes present in adata).
  3. sc.tl.score_genes per label -> score column.
  4. rank(ascending=False, pct=True) -> score_rank_pct column.
  5. best_rank_type_global = label with smallest score_rank_pct.
  6. Save unfiltered scored object:
     h5ad/07-score-rank-qc/adata_anno_score_genes_rank.h5ad
  7. Filter: keep only cells where best_rank_type_global == leiden_coarse.
     Save consistent object:
     h5ad/07-score-rank-qc/adata_anno_score_genes_rank_consistent.h5ad
  8. Per-cell table, consistent/inconsistent summary tables.

Note: input is adata_anno_filtered.h5ad (LowQuality/Doublets already removed).
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

SEED = 42
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
IN_H5AD = ROOT / "epi-cm-core-workflow/h5ad/01-celltype_integration_clustering/06-broad-annotation/adata_anno_filtered.h5ad"
OUT_DIR = ROOT / "epi-cm-core-workflow/h5ad/01-celltype_integration_clustering/07-score-rank-qc"
OUT_DIR.mkdir(parents=True, exist_ok=True)
TAB = ROOT / "epi-cm-core-workflow/tables/01-celltype_integration_clustering/07-score-rank-qc"
TAB.mkdir(parents=True, exist_ok=True)
FIG = ROOT / "epi-cm-core-workflow/figures/01-celltype_integration_clustering/07-score-rank-qc"
FIG.mkdir(parents=True, exist_ok=True)

OUT_UNFILT = OUT_DIR / "adata_anno_score_genes_rank.h5ad"
OUT_FILT = OUT_DIR / "adata_anno_score_genes_rank_consistent.h5ad"

DEG_DIR = ROOT / "epi-cm-core-workflow/tables/01-celltype_integration_clustering/06-broad-annotation/degs_leiden_coarse_pcs25_nn30_res0p4"

N_TOP_GENES = 100


def sanitize(name: str) -> str:
    return name.replace(" ", "_").replace("/", "_").replace("-", "_")


def main() -> None:
    t0 = time.time()
    print(f"[qcrk] reading {IN_H5AD.name}…", flush=True)
    adata = sc.read_h5ad(IN_H5AD)
    print(f"[qcrk] loaded {adata.shape} in {time.time()-t0:.1f}s", flush=True)
    assert adata.raw is not None, "adata.raw missing"
    raw_var = set(adata.raw.var_names)
    print(f"[qcrk] raw var: {len(raw_var)}", flush=True)

    # observed leiden_coarse labels
    labels = sorted(adata.obs["leiden_coarse"].astype(str).unique().tolist())
    print(f"[qcrk] {len(labels)} leiden_coarse labels: {labels}", flush=True)

    # 1. read per-label DEG CSVs and pick first 100 usable genes
    score_gene_sets: dict[str, list[str]] = {}
    for lab in labels:
        safe = sanitize(lab)
        fpath = DEG_DIR / f"{safe}_degs_leiden_coarse_pcs25_nn30_res0p4.csv"
        assert fpath.exists(), f"missing DEG CSV: {fpath}"
        df = pd.read_csv(fpath)
        # Take first 100 rows (already sorted by rank; gene column is 'gene')
        usable = []
        for g in df["gene"].tolist():
            if g in raw_var:
                usable.append(g)
            if len(usable) >= N_TOP_GENES:
                break
        score_gene_sets[lab] = usable
        print(f"[qcrk] {lab}: {len(usable)} usable genes", flush=True)

    # Save gene-set mapping
    with open(TAB / "score_gene_sets.json", "w") as f:
        json.dump(score_gene_sets, f, indent=2)
    print(f"[qcrk] wrote score_gene_sets.json", flush=True)

    # Save score column name mapping
    score_col_map = {lab: f"{sanitize(lab)}_score" for lab in labels}
    rank_col_map = {lab: f"{sanitize(lab)}_score_rank_pct" for lab in labels}
    with open(TAB / "score_column_mapping.json", "w") as f:
        json.dump({"score_columns": score_col_map,
                   "rank_pct_columns": rank_col_map}, f, indent=2)
    print(f"[qcrk] wrote score_column_mapping.json", flush=True)

    # 2. score_genes per label (use_raw=True by default in scanpy)
    for lab in labels:
        sc.tl.score_genes(adata, gene_list=score_gene_sets[lab],
                          score_name=score_col_map[lab], random_state=SEED)
        print(f"[qcrk] score_genes done for {lab} -> {score_col_map[lab]}", flush=True)

    # 3. rank percentile per score column (rank(ascending=False, pct=True))
    # 1.0 = top; 0.0 = bottom
    for lab in labels:
        s = adata.obs[score_col_map[lab]].astype(float)
        rpct = s.rank(ascending=False, pct=True)
        adata.obs[rank_col_map[lab]] = rpct
        print(f"[qcrk] rank_pct done for {lab} -> {rank_col_map[lab]}", flush=True)

    # 4. best_rank_type_global = label with smallest rank_pct
    rank_df = pd.DataFrame({lab: adata.obs[rank_col_map[lab]] for lab in labels},
                           index=adata.obs_names)
    best_label = rank_df.idxmin(axis=1)
    adata.obs["best_rank_type_global"] = pd.Categorical(best_label, categories=labels)
    print(f"[qcrk] best_rank_type_global value counts:")
    print(adata.obs["best_rank_type_global"].value_counts().to_string(), flush=True)

    # 5. consistency: best_rank_type_global == leiden_coarse
    is_consistent = (adata.obs["best_rank_type_global"].astype(str) ==
                     adata.obs["leiden_coarse"].astype(str))
    adata.obs["rank_consistent"] = is_consistent
    print(f"[qcrk] consistent cells: {is_consistent.sum()} / {len(is_consistent)} "
          f"({100*is_consistent.mean():.2f}%)", flush=True)

    # Per-cell table
    per_cell = pd.DataFrame({
        "leiden_coarse": adata.obs["leiden_coarse"].astype(str).values,
        "best_rank_type_global": adata.obs["best_rank_type_global"].astype(str).values,
        "rank_consistent": is_consistent.values,
    })
    for lab in labels:
        per_cell[score_col_map[lab]] = adata.obs[score_col_map[lab]].values
        per_cell[rank_col_map[lab]] = adata.obs[rank_col_map[lab]].values
    per_cell.to_csv(TAB / "per_cell_score_rank.csv", index=True)
    print(f"[qcrk] wrote per_cell_score_rank.csv ({len(per_cell)} rows)", flush=True)

    # Per-celltype summary
    summary = per_cell.groupby(["leiden_coarse", "rank_consistent"]).size().unstack(fill_value=0)
    summary.columns = [f"{c}" for c in summary.columns]
    summary["total"] = summary.sum(axis=1)
    summary["pct_consistent"] = (summary[True] / summary["total"] * 100).round(2) if True in summary.columns else 0
    summary.to_csv(TAB / "consistency_summary_by_leiden_coarse.csv")
    print(f"[qcrk] consistency summary by leiden_coarse:")
    print(summary.to_string(), flush=True)

    # 6. Save unfiltered scored object
    print(f"[qcrk] writing {OUT_UNFILT.name}…", flush=True)
    adata.write_h5ad(OUT_UNFILT, compression="gzip")

    # 7. Save filtered consistent object
    adata_filt = adata[is_consistent].copy()
    print(f"[qcrk] filtered consistent: {adata_filt.shape}", flush=True)
    adata_filt.write_h5ad(OUT_FILT, compression="gzip")
    print(f"[qcrk] wrote {OUT_FILT.name}", flush=True)

    # 8. Major-lineage UMAP and dotplot
    sc.settings.figdir = str(FIG)
    sc.set_figure_params(figsize=(5, 5), dpi=150)
    sc.pl.umap(adata, color="leiden_coarse", save="_major_lineage_umap.pdf", show=False)
    sc.pl.umap(adata, color="best_rank_type_global", save="_best_rank_type_global_umap.pdf", show=False)
    print(f"[qcrk] UMAP figures saved", flush=True)

    # Canonical marker dotplot: top 3 per leiden_coarse label from full DEG tables
    top3_markers = []
    for lab in labels:
        safe = sanitize(lab)
        fpath = DEG_DIR / f"{safe}_degs_leiden_coarse_pcs25_nn30_res0p4.csv"
        df = pd.read_csv(fpath)
        # take top 3 by rank (already sorted)
        for g in df["gene"].head(3).tolist():
            if g not in top3_markers:
                top3_markers.append((lab, g))
    # Build marker dict: {label: [3 markers]} preserving label order
    marker_dict = {}
    for lab in labels:
        genes = [g for (l, g) in top3_markers if l == lab]
        marker_dict[lab] = genes
    print(f"[qcrk] marker dict (3 per label):")
    for lab, gs in marker_dict.items():
        print(f"  {lab}: {gs}", flush=True)
    if marker_dict:
        sc.pl.dotplot(
            adata,
            var_names=marker_dict,
            groupby="leiden_coarse",
            standard_scale="var",
            use_raw=True,
            save="_all_leiden_coarse.pdf",
            show=False,
        )
        print(f"[qcrk] canonical dotplot saved", flush=True)

    # Save summary json
    summary_json = {
        "n_cells_input": int(adata.n_obs),
        "n_obs_score": int(adata.n_obs),
        "n_obs_consistent": int(is_consistent.sum()),
        "n_labels": len(labels),
        "labels": labels,
        "n_top_genes_per_label": N_TOP_GENES,
        "score_columns": score_col_map,
        "rank_pct_columns": rank_col_map,
        "per_label_consistent_pct": {lab: float(per_cell[(per_cell["leiden_coarse"]==lab) & per_cell["rank_consistent"]].shape[0] /
                                                  max(1, per_cell[per_cell["leiden_coarse"]==lab].shape[0]) * 100)
                                      for lab in labels},
    }
    with open(TAB / "score_rank_summary.json", "w") as f:
        json.dump(summary_json, f, indent=2)
    print(f"[qcrk] score_rank_summary.json written", flush=True)
    print(f"[qcrk] done. total={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()