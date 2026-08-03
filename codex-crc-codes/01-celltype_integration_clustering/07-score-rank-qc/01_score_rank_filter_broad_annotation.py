#!/usr/bin/env python3
"""Broad annotation score/rank QC and consistent-cell handoff."""

from __future__ import annotations

import platform
import random
import re
import traceback
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc

matplotlib.use("Agg")

SEED = 42
random.seed(SEED)
np.random.seed(SEED)

WORKFLOW_ROOT = Path("/mnt/disk18t/lr_xcy/riku/crc_val/epi-cm-core-workflow")
STEP = "01-celltype_integration_clustering/07-score-rank-qc"
INPUT_H5AD = WORKFLOW_ROOT / "h5ad/01-celltype_integration_clustering/06-broad-annotation/adata_anno.h5ad"
OUTPUT_H5AD = WORKFLOW_ROOT / "h5ad" / STEP / "adata_anno_score_genes_rank.h5ad"
OUTPUT_CONSISTENT_H5AD = (
    WORKFLOW_ROOT / "h5ad" / STEP / "adata_anno_score_genes_rank_consistent.h5ad"
)
TABLE_DIR = WORKFLOW_ROOT / "tables" / STEP
FIGURE_DIR = WORKFLOW_ROOT / "figures" / STEP
DEG_DIR = (
    WORKFLOW_ROOT
    / "tables/01-celltype_integration_clustering/06-broad-annotation/degs_leiden_coarse_pcs30_nn30_res0p3"
)
CODE_FILE = Path(__file__)

GROUPBY = "leiden_coarse"
CELL_TYPE_COL = "cell_type"
RAW_CLUSTER_COL = "leiden_res0p3"
N_TOP_DEG_GENES = 100
CTRL_SIZE = 50
N_BINS = 25


def sanitize_label(label: str) -> str:
    text = re.sub(r"[^A-Za-z0-9]+", "_", label).strip("_")
    return text or "label"


def assert_no_overwrite() -> None:
    outputs = [
        OUTPUT_H5AD,
        OUTPUT_CONSISTENT_H5AD,
        TABLE_DIR / "01_score_gene_set_manifest.csv",
        TABLE_DIR / "01_score_column_label_mapping.csv",
        TABLE_DIR / "01_score_rank_qc_per_cell.csv",
        TABLE_DIR / "01_score_rank_qc_consistent_cells.csv",
        TABLE_DIR / "01_score_rank_qc_inconsistent_cells.csv",
        TABLE_DIR / "01_score_rank_qc_summary_by_leiden_coarse.csv",
        TABLE_DIR / "01_score_rank_qc_parameters.csv",
        TABLE_DIR / "package_versions.txt",
        TABLE_DIR / "readme.txt",
        FIGURE_DIR / "umap_leiden_coarse_best_rank_consistency.pdf",
        FIGURE_DIR / "umap_leiden_coarse_consistent_cells.pdf",
    ]
    existing = [str(path) for path in outputs if path.exists()]
    if existing:
        raise FileExistsError(
            "Score/rank QC output already exists; refusing to overwrite:\n"
            + "\n".join(existing)
        )


def deg_path_for_label(label: str) -> Path:
    return DEG_DIR / f"{sanitize_label(label)}_degs_leiden_coarse_pcs30_nn30_res0p3.csv"


def build_gene_sets(adata: ad.AnnData, labels: list[str]) -> tuple[pd.DataFrame, dict[str, list[str]]]:
    raw_var = set(adata.raw.var_names.astype(str))
    rows = []
    gene_sets = {}
    for label in labels:
        path = deg_path_for_label(label)
        if not path.exists():
            raise FileNotFoundError(f"Missing DEG table for observed leiden_coarse label {label!r}: {path}")
        deg = pd.read_csv(path)
        if "ensembl" not in deg.columns:
            raise ValueError(f"DEG table lacks required 'ensembl' column: {path}")
        usable = []
        seen = set()
        for row in deg.itertuples(index=False):
            ensembl = str(getattr(row, "ensembl"))
            if ensembl in raw_var and ensembl not in seen:
                usable.append(ensembl)
                seen.add(ensembl)
            if len(usable) >= N_TOP_DEG_GENES:
                break
        if len(usable) < N_TOP_DEG_GENES:
            raise RuntimeError(
                f"Only {len(usable)} usable DEG genes found for {label}; expected {N_TOP_DEG_GENES}."
            )
        gene_sets[label] = usable
        gene_symbol_map = (
            adata.raw.var["gene_symbol"].astype(str).to_dict()
            if "gene_symbol" in adata.raw.var.columns
            else {}
        )
        for rank, ensembl in enumerate(usable, start=1):
            rows.append(
                {
                    "leiden_coarse": label,
                    "rank_in_deg_table": rank,
                    "ensembl": ensembl,
                    "gene_symbol": gene_symbol_map.get(ensembl, ""),
                    "deg_table": str(path),
                }
            )
    return pd.DataFrame(rows), gene_sets


def main() -> None:
    assert_no_overwrite()
    OUTPUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    sc.settings.autoshow = False
    sc.settings.figdir = str(FIGURE_DIR)
    sc.settings.set_figure_params(figsize=(3, 3), dpi=150)

    status = {
        "step": "score_rank_qc_start",
        "input_h5ad": str(INPUT_H5AD),
        "output_h5ad": str(OUTPUT_H5AD),
        "output_consistent_h5ad": str(OUTPUT_CONSISTENT_H5AD),
        "groupby": GROUPBY,
        "deg_dir": str(DEG_DIR),
        "n_top_deg_genes_per_label": N_TOP_DEG_GENES,
        "score_genes_use_raw": True,
        "ctrl_size": CTRL_SIZE,
        "n_bins": N_BINS,
        "random_seed": SEED,
        "code_file": str(CODE_FILE),
    }
    pd.DataFrame([status]).to_csv(TABLE_DIR / "01_score_rank_qc_status.csv", index=False)

    try:
        adata = ad.read_h5ad(INPUT_H5AD)
        if adata.raw is None:
            raise RuntimeError("adata.raw is absent; score_genes QC requires use_raw=True.")
        for col in (GROUPBY, CELL_TYPE_COL, RAW_CLUSTER_COL):
            if col not in adata.obs.columns:
                raise KeyError(f"Missing obs[{col!r}] in annotated h5ad.")

        labels = list(adata.obs[GROUPBY].cat.categories) if hasattr(adata.obs[GROUPBY].dtype, "categories") else []
        labels = [str(label) for label in labels if (adata.obs[GROUPBY].astype(str) == str(label)).any()]
        if not labels:
            labels = sorted(adata.obs[GROUPBY].astype(str).unique())

        gene_manifest, gene_sets = build_gene_sets(adata, labels)
        gene_manifest.to_csv(TABLE_DIR / "01_score_gene_set_manifest.csv", index=False)

        mapping_rows = []
        score_cols = []
        rank_cols = []
        for label in labels:
            safe = sanitize_label(label)
            score_col = f"{safe}_score"
            rank_col = f"{safe}_score_rank_pct"
            sc.tl.score_genes(
                adata,
                gene_list=gene_sets[label],
                score_name=score_col,
                ctrl_size=CTRL_SIZE,
                n_bins=N_BINS,
                random_state=SEED,
                use_raw=True,
            )
            adata.obs[rank_col] = adata.obs[score_col].rank(ascending=False, pct=True)
            score_cols.append(score_col)
            rank_cols.append(rank_col)
            mapping_rows.append(
                {
                    "leiden_coarse": label,
                    "score_col": score_col,
                    "rank_pct_col": rank_col,
                    "n_genes": len(gene_sets[label]),
                }
            )
        mapping = pd.DataFrame(mapping_rows)
        mapping.to_csv(TABLE_DIR / "01_score_column_label_mapping.csv", index=False)

        rank_matrix = adata.obs[rank_cols]
        best_rank_col = rank_matrix.idxmin(axis=1)
        rank_to_label = dict(zip(mapping["rank_pct_col"], mapping["leiden_coarse"]))
        score_to_label = dict(zip(mapping["score_col"], mapping["leiden_coarse"]))
        best_labels = best_rank_col.map(rank_to_label)
        adata.obs["best_rank_type_global"] = pd.Categorical(best_labels, categories=labels)
        adata.obs["best_rank_score_rank_pct"] = rank_matrix.min(axis=1).astype(float)
        best_score_cols = best_labels.map({v: k for k, v in score_to_label.items()})
        adata.obs["best_rank_score_col"] = best_score_cols.astype(str)
        adata.obs["score_rank_consistent"] = (
            adata.obs["best_rank_type_global"].astype(str) == adata.obs[GROUPBY].astype(str)
        )

        per_cell_cols = [
            GROUPBY,
            CELL_TYPE_COL,
            RAW_CLUSTER_COL,
            "best_rank_type_global",
            "best_rank_score_rank_pct",
            "best_rank_score_col",
            "score_rank_consistent",
        ]
        per_cell = adata.obs[per_cell_cols].copy()
        per_cell.insert(0, "cell_id", adata.obs_names)
        per_cell.to_csv(TABLE_DIR / "01_score_rank_qc_per_cell.csv", index=False)
        per_cell[per_cell["score_rank_consistent"]].to_csv(
            TABLE_DIR / "01_score_rank_qc_consistent_cells.csv", index=False
        )
        per_cell[~per_cell["score_rank_consistent"]].to_csv(
            TABLE_DIR / "01_score_rank_qc_inconsistent_cells.csv", index=False
        )

        summary = (
            per_cell.groupby(GROUPBY, observed=True)["score_rank_consistent"]
            .agg(n_cells="size", n_consistent="sum")
            .reset_index()
        )
        summary["n_inconsistent"] = summary["n_cells"] - summary["n_consistent"]
        summary["pct_consistent"] = summary["n_consistent"] / summary["n_cells"] * 100
        summary.to_csv(TABLE_DIR / "01_score_rank_qc_summary_by_leiden_coarse.csv", index=False)

        sc.pl.umap(
            adata,
            color=[GROUPBY, "best_rank_type_global", "score_rank_consistent"],
            ncols=3,
            wspace=0.45,
            save="_leiden_coarse_best_rank_consistency.pdf",
            show=False,
        )

        consistent_mask = adata.obs["score_rank_consistent"].astype(bool).to_numpy()
        consistent = adata[consistent_mask].copy()
        sc.pl.umap(
            consistent,
            color=GROUPBY,
            save="_leiden_coarse_consistent_cells.pdf",
            show=False,
        )

        adata.write_h5ad(OUTPUT_H5AD)
        consistent.write_h5ad(OUTPUT_CONSISTENT_H5AD)

        params = {
            **status,
            "step": "score_rank_qc_complete",
            "n_obs_input": int(adata.n_obs),
            "n_obs_consistent": int(consistent.n_obs),
            "n_obs_inconsistent": int(adata.n_obs - consistent.n_obs),
            "pct_consistent": float(consistent.n_obs / adata.n_obs * 100),
            "n_labels_scored": len(labels),
            "score_columns": ";".join(score_cols),
            "rank_pct_columns": ";".join(rank_cols),
            "filter_rule": "keep cells where best_rank_type_global equals leiden_coarse",
        }
        pd.DataFrame([params]).to_csv(TABLE_DIR / "01_score_rank_qc_parameters.csv", index=False)
        pd.DataFrame([params]).to_csv(TABLE_DIR / "01_score_rank_qc_status.csv", index=False)

        with (TABLE_DIR / "package_versions.txt").open("w") as fh:
            fh.write(f"python: {platform.python_version()}\n")
            fh.write(f"anndata: {ad.__version__}\n")
            fh.write(f"scanpy: {sc.__version__}\n")
            fh.write(f"pandas: {pd.__version__}\n")
            fh.write(f"numpy: {np.__version__}\n")
            fh.write(f"code_file: {CODE_FILE}\n")

        with (TABLE_DIR / "readme.txt").open("w") as fh:
            fh.write("07-score-rank-qc completed.\n")
            fh.write(f"Input annotated h5ad: {INPUT_H5AD}\n")
            fh.write(f"Unfiltered scored h5ad: {OUTPUT_H5AD}\n")
            fh.write(f"Consistent-cell handoff h5ad: {OUTPUT_CONSISTENT_H5AD}\n")
            fh.write("Score gene sets were selected from current leiden_coarse Round2 DEG tables only.\n")
            fh.write(f"Top usable DEG genes per label: {N_TOP_DEG_GENES}\n")
            fh.write("Filter rule: best_rank_type_global must match leiden_coarse.\n")
            fh.write("Expression source: adata.raw via sc.tl.score_genes(use_raw=True).\n")

        print(pd.DataFrame([params]).to_string(index=False))

    except Exception as exc:
        status.update(
            {
                "step": "score_rank_qc_failed",
                "error_summary": "".join(traceback.format_exception_only(type(exc), exc)).strip(),
            }
        )
        pd.DataFrame([status]).to_csv(TABLE_DIR / "01_score_rank_qc_status.csv", index=False)
        raise


if __name__ == "__main__":
    main()
