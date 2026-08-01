#!/usr/bin/env python3
"""Run BRCA broad-lineage DEG-derived score/rank consistency QC."""

from __future__ import annotations

import importlib.metadata
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


ROOT = Path("/mnt/disk18t/lr_xcy/riku/brca_cm-epi_codex")
WORKFLOW = ROOT / "epi-cm-core-workflow"
BLOCK = "01-celltype_integration_clustering"
INPUT_H5AD = WORKFLOW / "h5ad" / BLOCK / "06-broad-annotation" / "adata_anno.h5ad"
OUT_UNFILTERED = (
    WORKFLOW
    / "h5ad"
    / BLOCK
    / "07-score-rank-qc"
    / "adata_anno_score_genes_rank.h5ad"
)
OUT_FILTERED = (
    WORKFLOW
    / "h5ad"
    / BLOCK
    / "07-score-rank-qc"
    / "adata_anno_score_genes_rank_consistent.h5ad"
)
TABLE_DIR = WORKFLOW / "tables" / BLOCK / "07-score-rank-qc"
FIGURE_DIR = WORKFLOW / "figures" / BLOCK / "07-score-rank-qc"
DEG_DIR = (
    WORKFLOW
    / "tables"
    / BLOCK
    / "06-broad-annotation"
    / "degs_leiden_coarse_pcs20_nn30_res0p8_myo_merged"
)
CODE_PATH = (
    WORKFLOW
    / "codes"
    / BLOCK
    / "07-score-rank-qc"
    / "11_score_rank_qc.py"
)

COARSE_KEY = "leiden_coarse"
RAW_KEY = "leiden_res0p8"
BEST_KEY = "best_rank_type_global"
CONSISTENT_KEY = "score_rank_consistent"
N_SCORE_GENES = 100
CTRL_SIZE = 50
N_BINS = 25
SEED = 42


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_").lower()


def find_deg_file(label: str) -> Path:
    files = sorted(DEG_DIR.glob(f"{re.sub(r'[^A-Za-z0-9]+', '_', label).strip('_')}_degs_*.csv"))
    if len(files) != 1:
        raise ValueError(f"Expected one saved full DEG table for {label}, found {files}")
    return files[0]


def choose_score_genes(deg: pd.DataFrame, genes_present: set[str]) -> list[str]:
    required = {"gene", "logfoldchanges", "pvals_adj"}
    if not required.issubset(deg.columns):
        raise ValueError(f"DEG table lacks required columns: {sorted(required - set(deg.columns))}")
    usable: list[str] = []
    seen: set[str] = set()
    for row in deg.itertuples(index=False):
        gene = str(row.gene)
        if gene in seen or gene not in genes_present:
            continue
        if not np.isfinite(float(row.logfoldchanges)) or float(row.logfoldchanges) <= 0:
            continue
        if not np.isfinite(float(row.pvals_adj)) or float(row.pvals_adj) >= 0.05:
            continue
        if gene.startswith(("MT-", "RPL", "RPS")):
            continue
        usable.append(gene)
        seen.add(gene)
        if len(usable) == N_SCORE_GENES:
            return usable
    raise ValueError(
        f"Only {len(usable)} usable positive non-MT/non-ribosomal genes found; "
        f"{N_SCORE_GENES} are required."
    )


def main() -> None:
    started = time.time()
    OUT_UNFILTERED.parent.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(INPUT_H5AD)
    if adata.raw is None:
        raise ValueError("score_genes(use_raw=True) requires adata.raw.")
    for key in [COARSE_KEY, "cell_type", RAW_KEY, "sample", "series", "status"]:
        if key not in adata.obs:
            raise ValueError(f"Required annotation/metadata column is missing: {key}")
    if not adata.obs[COARSE_KEY].astype(str).equals(adata.obs["cell_type"].astype(str)):
        raise ValueError("Initial cell_type must equal leiden_coarse before subtype analysis.")
    labels = adata.obs[COARSE_KEY].cat.categories.astype(str).tolist()
    if set(labels) != set(adata.obs[COARSE_KEY].astype(str).unique()):
        raise ValueError("leiden_coarse contains unused or missing categories.")
    deg_files = sorted(DEG_DIR.glob("*.csv"))
    if len(deg_files) != len(labels):
        raise ValueError("The saved round-2 DEG directory does not match observed broad labels.")

    genes_present = set(adata.raw.var_names.astype(str))
    score_to_label: dict[str, str] = {}
    score_to_rank: dict[str, str] = {}
    gene_set_rows: list[dict[str, object]] = []
    print(
        f"[Score/rank] cells={adata.n_obs}, broad_labels={len(labels)}, "
        f"score_genes_per_label={N_SCORE_GENES}, use_raw=True",
        flush=True,
    )
    for label in labels:
        source = find_deg_file(label)
        deg = pd.read_csv(source)
        if len(deg) != adata.raw.n_vars:
            raise ValueError(f"Saved DEG table is not full length: {source}")
        genes = choose_score_genes(deg, genes_present)
        score_col = f"{safe_name(label)}_score"
        rank_col = f"{safe_name(label)}_score_rank_pct"
        if score_col in adata.obs or rank_col in adata.obs:
            raise ValueError(f"Score/rank column already exists: {score_col} / {rank_col}")
        sc.tl.score_genes(
            adata,
            gene_list=genes,
            ctrl_size=CTRL_SIZE,
            gene_pool=adata.raw.var_names,
            n_bins=N_BINS,
            score_name=score_col,
            random_state=SEED,
            use_raw=True,
        )
        score_to_label[score_col] = label
        score_to_rank[score_col] = rank_col
        gene_set_rows.extend(
            {
                "leiden_coarse": label,
                "score_column": score_col,
                "rank_column": rank_col,
                "gene_rank": rank,
                "gene": gene,
                "source_full_deg_csv": str(source),
                "selection_rule": (
                    "first 100 positive padj<0.05 genes present in adata.raw, "
                    "excluding MT-/RPL/RPS"
                ),
            }
            for rank, gene in enumerate(genes, start=1)
        )
        print(f"[Score/rank] scored {label} with {len(genes)} DEG-derived genes", flush=True)

    score_cols = list(score_to_label)
    score_frame = adata.obs[score_cols].astype(float)
    if not np.isfinite(score_frame.to_numpy()).all():
        raise ValueError("Broad score columns contain non-finite values.")
    rank_frame = score_frame.rank(
        axis=0, method="average", ascending=False, pct=True
    )
    for score_col, rank_col in score_to_rank.items():
        adata.obs[rank_col] = rank_frame[score_col].to_numpy(dtype=float)
    rank_cols = [score_to_rank[column] for column in score_cols]
    best_rank_col = adata.obs[rank_cols].idxmin(axis=1)
    rank_to_label = {
        rank_col: score_to_label[score_col]
        for score_col, rank_col in score_to_rank.items()
    }
    best_by_rank = best_rank_col.map(rank_to_label)
    adata.obs[BEST_KEY] = pd.Categorical(
        best_by_rank, categories=labels, ordered=True
    )
    adata.obs[CONSISTENT_KEY] = (
        adata.obs[BEST_KEY].astype(str) == adata.obs[COARSE_KEY].astype(str)
    )
    consistent = adata.obs[CONSISTENT_KEY].astype(bool)
    n_kept = int(consistent.sum())
    n_removed = int((~consistent).sum())
    if n_kept == 0 or n_removed == adata.n_obs:
        raise ValueError("Score/rank consistency filtering retained no cells.")

    pd.DataFrame(gene_set_rows).to_csv(
        TABLE_DIR / "broad_score_gene_sets.csv", index=False
    )
    pd.DataFrame(
        [
            {
                "leiden_coarse": score_to_label[score_col],
                "score_column": score_col,
                "rank_column": score_to_rank[score_col],
                "n_score_genes": N_SCORE_GENES,
                "use_raw": True,
            }
            for score_col in score_cols
        ]
    ).to_csv(TABLE_DIR / "broad_score_column_mapping.csv", index=False)

    per_cell_cols = [
        COARSE_KEY,
        "cell_type",
        RAW_KEY,
        BEST_KEY,
        CONSISTENT_KEY,
        *score_cols,
        *rank_cols,
    ]
    per_cell = adata.obs[per_cell_cols].copy()
    per_cell.index.name = "cell_id"
    per_cell.to_csv(TABLE_DIR / "broad_score_rank_per_cell.csv.gz", compression="gzip")
    per_cell.loc[consistent].to_csv(
        TABLE_DIR / "broad_score_rank_consistent_cells.csv.gz", compression="gzip"
    )
    per_cell.loc[~consistent].to_csv(
        TABLE_DIR / "broad_score_rank_inconsistent_cells.csv.gz", compression="gzip"
    )

    summary = (
        adata.obs.groupby(COARSE_KEY, observed=True)[CONSISTENT_KEY]
        .agg(n_input="size", n_kept="sum")
        .reindex(labels)
    )
    summary["n_removed"] = summary["n_input"] - summary["n_kept"]
    summary["retained_fraction"] = summary["n_kept"] / summary["n_input"]
    summary.to_csv(TABLE_DIR / "score_rank_consistency_by_leiden_coarse.csv")
    raw_summary = (
        adata.obs.groupby(RAW_KEY, observed=True)[CONSISTENT_KEY]
        .agg(n_input="size", n_kept="sum")
    )
    raw_summary["n_removed"] = raw_summary["n_input"] - raw_summary["n_kept"]
    raw_summary["retained_fraction"] = raw_summary["n_kept"] / raw_summary["n_input"]
    raw_summary.to_csv(TABLE_DIR / "score_rank_consistency_by_raw_cluster.csv")
    pd.crosstab(
        adata.obs[COARSE_KEY].astype(str),
        adata.obs[BEST_KEY].astype(str),
    ).reindex(index=labels, columns=labels, fill_value=0).to_csv(
        TABLE_DIR / "leiden_coarse_vs_best_rank_counts.csv"
    )

    adata.uns["score_rank_qc_parameters"] = {
        "source_h5ad": str(INPUT_H5AD),
        "source_full_deg_dir": str(DEG_DIR),
        "observed_broad_labels": labels,
        "n_score_genes_per_label": N_SCORE_GENES,
        "score_gene_selection": (
            "first 100 positive padj<0.05 genes present in adata.raw, "
            "excluding MT-/RPL/RPS"
        ),
        "score_genes_use_raw": True,
        "score_genes_ctrl_size": CTRL_SIZE,
        "score_genes_n_bins": N_BINS,
        "rank_method": "per-score-column descending percentile; smaller is better",
        "best_label_rule": "label with smallest rank percentile",
        "consistency_rule": "best_rank_type_global exactly equals leiden_coarse",
        "seed": SEED,
    }

    sc.settings.figdir = FIGURE_DIR
    sc.settings.file_format_figs = "pdf"
    sc.set_figure_params(figsize=(3, 3), dpi=150, fontsize=8)
    sc.pl.umap(
        adata,
        color=[COARSE_KEY, BEST_KEY],
        ncols=2,
        wspace=0.5,
        show=False,
        save="_leiden_coarse_vs_best_rank.pdf",
    )
    sc.pl.umap(
        adata,
        color=CONSISTENT_KEY,
        show=False,
        save="_score_rank_consistency.pdf",
    )

    print(f"[Score/rank] writing unfiltered scored object: {OUT_UNFILTERED}", flush=True)
    adata.write_h5ad(OUT_UNFILTERED, compression="gzip")
    filtered = adata[consistent].copy()
    if not filtered.obs[CONSISTENT_KEY].astype(bool).all():
        raise ValueError("Filtered object still contains inconsistent cells.")
    sc.pl.umap(
        filtered,
        color=[COARSE_KEY, BEST_KEY],
        ncols=2,
        wspace=0.5,
        show=False,
        save="_consistent_cells_leiden_coarse_vs_best_rank.pdf",
    )
    if n_removed:
        inconsistent_adata = adata[~consistent].copy()
        sc.pl.umap(
            inconsistent_adata,
            color=[COARSE_KEY, BEST_KEY],
            ncols=2,
            wspace=0.5,
            show=False,
            save="_inconsistent_cells_leiden_coarse_vs_best_rank.pdf",
        )
        del inconsistent_adata
    print(f"[Score/rank] writing consistent handoff: {OUT_FILTERED}", flush=True)
    filtered.write_h5ad(OUT_FILTERED, compression="gzip")

    versions = {
        "python": sys.version.split()[0],
        "scanpy": package_version("scanpy"),
        "anndata": package_version("anndata"),
        "code": str(CODE_PATH),
        "seed": str(SEED),
    }
    (TABLE_DIR / "package_versions.txt").write_text(
        "\n".join(f"{key}={value}" for key, value in versions.items()) + "\n",
        encoding="utf-8",
    )
    (TABLE_DIR / "readme.txt").write_text(
        f"""BRCA broad score/rank consistency QC

Input: {INPUT_H5AD}
Unfiltered scored output: {OUT_UNFILTERED}
Default downstream handoff: {OUT_FILTERED}
Source full broad-label DEG tables: {DEG_DIR}

For each of the {len(labels)} observed BRCA broad labels, the score gene set is the
first {N_SCORE_GENES} positive adjusted-P<0.05 genes present in adata.raw from
that label's saved full-length DEG table, excluding MT-/RPL/RPS genes. No
preset lineage list or marker list is used to build scores. score_genes is run
with use_raw=True. Per-score descending rank percentiles are stored; smaller
rank percentile is better. A cell is retained only when best_rank_type_global
equals its current leiden_coarse label. Both unfiltered and filtered objects
are preserved, and the filtered object is the default downstream handoff.
""",
        encoding="utf-8",
    )
    report = {
        "n_input_cells": int(adata.n_obs),
        "n_consistent_cells": n_kept,
        "n_inconsistent_cells": n_removed,
        "retained_fraction": n_kept / adata.n_obs,
        "n_observed_broad_labels": len(labels),
        "n_score_genes_per_label": N_SCORE_GENES,
        "use_raw": True,
        "unfiltered_h5ad": str(OUT_UNFILTERED),
        "consistent_h5ad": str(OUT_FILTERED),
        "elapsed_seconds": time.time() - started,
    }
    (TABLE_DIR / "score_rank_completion.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
