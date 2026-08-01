"""QC stage A — initial inspection & threshold selection.

Loads adata_merge.h5ad, runs an initial min_genes cell filter, then prints
per-sample count/genes summary so the agent can pick dataset-appropriate
thresholds for the final QC pass.

Per SKILL.md:
  - Filter cells using an agent-selected initial cell-level minimum gene threshold.
  - Compute QC metrics (n_genes, total_counts, pct_counts_mt, pct_counts_ribo).
  - Agent picks thresholds from observed distributions.

Outputs:
  - tables/03-qc/01_initial_cell_filter_parameters.csv
  - tables/03-qc/qc_report_stageA.csv (per-sample count/genes summary)
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

SEED = 42
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
IN_H5AD = ROOT / "epi-cm-core-workflow/h5ad/01-celltype_integration_clustering/02-merge-metadata/adata_merge.h5ad"
TAB = ROOT / "epi-cm-core-workflow/tables/01-celltype_integration_clustering/03-qc"
TAB.mkdir(parents=True, exist_ok=True)

# Agent-selected initial min_genes (post-QC 591k cells × 21.5k genes was the prior
# run's output; using 200 here as a permissive floor to drop empty barcodes only)
INITIAL_MIN_GENES = 200


def main() -> None:
    t0 = time.time()
    print(f"[qc-A] reading {IN_H5AD.name}…", flush=True)
    adata = sc.read_h5ad(IN_H5AD)
    print(f"[qc-A] loaded {adata.shape} in {time.time()-t0:.1f}s", flush=True)

    n_before, n_var_before = adata.shape
    print(f"[qc-A] before initial cell filter: n_obs={n_before}, n_vars={n_var_before}", flush=True)

    # Initial cell filter
    sc.pp.filter_cells(adata, min_genes=INITIAL_MIN_GENES)
    n_after = adata.n_obs
    print(f"[qc-A] after initial cell filter (min_genes={INITIAL_MIN_GENES}): {n_after} (dropped {n_before-n_after})", flush=True)

    # Write initial-cell-filter parameter table
    p1 = TAB / "01_initial_cell_filter_parameters.csv"
    pd.DataFrame([{
        "step": "01_initial_cell_filter",
        "input_h5ad": str(IN_H5AD),
        "output_h5ad_or_object": "in-memory adata",
        "n_obs_before": n_before,
        "n_obs_after": n_after,
        "n_vars_before": n_var_before,
        "n_vars_after": adata.n_vars,
        "parameters": json.dumps({"min_genes": INITIAL_MIN_GENES}),
        "backend": "scanpy",
        "code_file": "epi-cm-core-workflow/codes/01-celltype_integration_clustering/03-qc/02_inspect_raw.py",
        "random_seed": SEED,
        "notes": "permissive floor to drop empty barcodes only; final thresholds selected in stage C",
    }]).to_csv(p1, index=False)
    print(f"[qc-A] wrote {p1.name}", flush=True)

    # Mark MT/RIBO var FIRST so calculate_qc_metrics can compute their % columns
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    n_mt = int(adata.var["mt"].sum())
    n_ribo = int(adata.var["ribo"].sum())
    print(f"[qc-A] MT genes (MT-*): {n_mt}; ribo genes (RPS*/RPL*): {n_ribo}", flush=True)

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    print(f"[qc-A] obs cols after qc_metrics: {list(adata.obs.columns)}", flush=True)

    # Per-sample summary (median + IQR)
    grp = adata.obs.groupby("sample", observed=True)
    rows = []
    for sample_id, sub in grp:
        rows.append({
            "sample": sample_id,
            "series": sub["series"].iloc[0],
            "status": sub["status"].iloc[0],
            "n_cells": sub.shape[0],
            "median_total_counts": float(sub["total_counts"].median()),
            "p25_total_counts": float(sub["total_counts"].quantile(0.25)),
            "p75_total_counts": float(sub["total_counts"].quantile(0.75)),
            "median_n_genes_by_counts": float(sub["n_genes_by_counts"].median()),
            "p25_n_genes": float(sub["n_genes_by_counts"].quantile(0.25)),
            "p75_n_genes": float(sub["n_genes_by_counts"].quantile(0.75)),
            "median_pct_counts_mt": float(sub["pct_counts_mt"].median()),
            "p95_pct_counts_mt": float(sub["pct_counts_mt"].quantile(0.95)),
            "median_pct_counts_ribo": float(sub["pct_counts_ribo"].median()),
            "p95_pct_counts_ribo": float(sub["pct_counts_ribo"].quantile(0.95)),
        })
    rep = pd.DataFrame(rows).sort_values("n_cells", ascending=False)
    rep_path = TAB / "qc_report_stageA.csv"
    rep.to_csv(rep_path, index=False)
    print(f"[qc-A] per-sample summary written to {rep_path.name}", flush=True)
    print(rep.head(15).to_string(), flush=True)

    # Save threshold-selection hints
    summary = {
        "n_obs_after_initial_filter": n_after,
        "n_vars": adata.n_vars,
        "n_mt_genes": n_mt,
        "n_ribo_genes": n_ribo,
        "median_total_counts_overall": float(adata.obs["total_counts"].median()),
        "p05_total_counts_overall": float(adata.obs["total_counts"].quantile(0.05)),
        "p95_total_counts_overall": float(adata.obs["total_counts"].quantile(0.95)),
        "median_n_genes_overall": float(adata.obs["n_genes_by_counts"].median()),
        "p05_n_genes_overall": float(adata.obs["n_genes_by_counts"].quantile(0.05)),
        "p95_n_genes_overall": float(adata.obs["n_genes_by_counts"].quantile(0.95)),
        "median_pct_MT_overall": float(adata.obs["pct_counts_mt"].median()),
        "p95_pct_MT_overall": float(adata.obs["pct_counts_mt"].quantile(0.95)),
        "median_pct_ribo_overall": float(adata.obs["pct_counts_ribo"].median()),
        "p95_pct_ribo_overall": float(adata.obs["pct_counts_ribo"].quantile(0.95)),
    }
    with open(TAB / "qc_stageA_thresholds.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=2), flush=True)

    # Cache the post-stage-A adata for stage B/C
    cache = ROOT / "epi-cm-core-workflow/h5ad/01-celltype_integration_clustering/03-qc/adata_qc_stageA.h5ad"
    adata.write_h5ad(cache, compression="gzip")
    print(f"[qc-A] cached to {cache.name}", flush=True)
    print(f"[qc-A] done. total={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()