#!/usr/bin/env python3
"""Export full-length round-1 DEGs for the selected raw Leiden clusters."""

from __future__ import annotations

import importlib.metadata
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


ROOT = Path("/mnt/disk18t/lr_xcy/riku/brca_cm-epi_codex")
WORKFLOW = ROOT / "epi-cm-core-workflow"
BLOCK = "01-celltype_integration_clustering"
INPUT_H5AD = (
    WORKFLOW
    / "h5ad"
    / BLOCK
    / "05-clustering-parameter-search"
    / "selected"
    / "adata_inte.h5ad"
)
TABLE_DIR = WORKFLOW / "tables" / BLOCK / "06-broad-annotation"
DEG_DIR = TABLE_DIR / "degs_leiden_res0p8_pcs20_nn30_res0p8"
CODE_PATH = (
    WORKFLOW
    / "codes"
    / BLOCK
    / "06-broad-annotation"
    / "09_round1_raw_cluster_degs.py"
)
CLUSTER_KEY = "leiden_res0p8"
METHOD = "t-test"
SEED = 42


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def main() -> None:
    started = time.time()
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    DEG_DIR.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(INPUT_H5AD)
    if adata.raw is None:
        raise ValueError("Round-1 broad DEGs require adata.raw normalized/log expression.")
    if CLUSTER_KEY not in adata.obs or adata.obs[CLUSTER_KEY].isna().any():
        raise ValueError(f"Missing or incomplete raw cluster column: {CLUSTER_KEY}")
    if "leiden_coarse" in adata.obs:
        raise ValueError("Round-1 raw-cluster DEGs must precede broad annotation.")
    groups = adata.obs[CLUSTER_KEY].cat.categories.astype(str).tolist()
    if len(groups) != 15:
        raise ValueError(f"Expected 15 selected raw clusters, found {len(groups)}")

    print(
        f"[Broad DEG round 1] cells={adata.n_obs}, raw_genes={adata.raw.n_vars}, "
        f"groups={len(groups)}, method={METHOD}, use_raw=True",
        flush=True,
    )
    sc.tl.rank_genes_groups(
        adata,
        groupby=CLUSTER_KEY,
        method=METHOD,
        use_raw=True,
        n_genes=adata.raw.n_vars,
        key_added="rank_genes_groups_raw_cluster",
    )
    result = adata.uns["rank_genes_groups_raw_cluster"]
    required = ["names", "scores", "logfoldchanges", "pvals", "pvals_adj"]
    missing = [key for key in required if key not in result]
    if missing:
        raise ValueError(f"rank_genes_groups result is missing fields: {missing}")

    audit_rows: list[dict[str, object]] = []
    for group in groups:
        frame = pd.DataFrame(
            {
                "gene": np.asarray(result["names"][group]).astype(str),
                "score": np.asarray(result["scores"][group], dtype=float),
                "logfoldchanges": np.asarray(
                    result["logfoldchanges"][group], dtype=float
                ),
                "pvals": np.asarray(result["pvals"][group], dtype=float),
                "pvals_adj": np.asarray(result["pvals_adj"][group], dtype=float),
            }
        )
        if len(frame) != adata.raw.n_vars:
            raise ValueError(f"DEG table for cluster {group} is not full length.")
        if frame["gene"].isna().any() or frame["gene"].duplicated().any():
            raise ValueError(f"DEG genes for cluster {group} are missing or duplicated.")
        out = DEG_DIR / f"{group}_degs_leiden_res0p8_pcs20_nn30_res0p8.csv"
        frame.to_csv(out, index=False)
        audit_rows.append(
            {
                "raw_cluster": group,
                "n_cells": int((adata.obs[CLUSTER_KEY].astype(str) == group).sum()),
                "n_deg_rows": len(frame),
                "full_length": True,
                "use_raw": True,
                "method": METHOD,
                "deg_csv": str(out),
            }
        )
        print(f"[Broad DEG round 1] wrote cluster={group} rows={len(frame)}", flush=True)

    pd.DataFrame(audit_rows).to_csv(
        TABLE_DIR / "round1_raw_cluster_deg_audit.csv", index=False
    )
    author_counts = pd.crosstab(
        adata.obs[CLUSTER_KEY].astype(str),
        adata.obs["celltype_major"].astype(str),
    )
    author_counts.to_csv(TABLE_DIR / "raw_cluster_by_published_major_counts.csv")
    author_fractions = author_counts.div(author_counts.sum(axis=1), axis=0)
    author_fractions.to_csv(TABLE_DIR / "raw_cluster_by_published_major_fractions.csv")
    pd.crosstab(
        adata.obs[CLUSTER_KEY].astype(str),
        adata.obs["celltype_minor"].astype(str),
    ).to_csv(TABLE_DIR / "raw_cluster_by_published_minor_counts.csv")

    versions = {
        "python": sys.version.split()[0],
        "scanpy": package_version("scanpy"),
        "anndata": package_version("anndata"),
        "code": str(CODE_PATH),
        "input_h5ad": str(INPUT_H5AD),
        "groupby": CLUSTER_KEY,
        "method": METHOD,
        "use_raw": "True",
        "review_depth_initial": "50",
        "seed": str(SEED),
    }
    (TABLE_DIR / "round1_package_versions.txt").write_text(
        "\n".join(f"{key}={value}" for key, value in versions.items()) + "\n",
        encoding="utf-8",
    )
    (TABLE_DIR / "round1_readme.txt").write_text(
        f"""BRCA broad annotation round-1 raw-cluster DEGs

Input: {INPUT_H5AD}
Groupby: {CLUSTER_KEY}
Method: {METHOD}
Expression: adata.raw normalized/log expression via use_raw=True
Output: {DEG_DIR}

One full-length DEG CSV with {adata.raw.n_vars} rows is written per raw cluster.
No top-only DEG CSV or top-N DEG manifest is produced. Marker review starts from
the first 50 rows of each already-saved full table and expands consistently only
if broad-lineage evidence is insufficient. Published author labels are exported
only as independent validation crosstabs; they do not replace DEG/marker review.
""",
        encoding="utf-8",
    )
    report = {
        "n_cells": int(adata.n_obs),
        "n_raw_genes": int(adata.raw.n_vars),
        "n_raw_clusters": len(groups),
        "n_full_length_deg_csvs": len(audit_rows),
        "rows_per_deg_csv": int(adata.raw.n_vars),
        "method": METHOD,
        "use_raw": True,
        "elapsed_seconds": time.time() - started,
    }
    (TABLE_DIR / "round1_completion.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
