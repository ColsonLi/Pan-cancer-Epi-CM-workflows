"""QC stage C — Scrublet per sample + gene/cell filters → adata_qc.h5ad.

Per SKILL.md default protocol:
  2. Run Scrublet per sample using batch_key = 'sample'.
  3. Remove cells with obs['predicted_doublet'] == True.
  6. Filter genes using min_cells threshold.
  7. Filter genes using min_counts threshold.
  8. Recalculate QC metrics after gene filtering before final cell filtering.
  9. Keep cells using final thresholds.
  10. Save h5ad/03-qc/adata_qc.h5ad.

User-selected thresholds (recorded here):
  - Scrublet: full rsc.pp.scrublet per sample, 84 samples
  - gene filter: min_cells=3
  - final cell filter: min_genes=200, max_pct_mt=15%, max_n_genes=6000

Outputs:
  - h5ad/01-celltype_integration_clustering/03-qc/adata_qc.h5ad
  - tables/01-celltype_integration_clustering/03-qc/02_scrublet_parameters.csv
  - tables/01-celltype_integration_clustering/03-qc/03_doublet_filter_parameters.csv
  - tables/01-celltype_integration_clustering/03-qc/05_gene_filter_parameters.csv
  - tables/01-celltype_integration_clustering/03-qc/06_final_cell_filter_parameters.csv
  - tables/01-celltype_integration_clustering/03-qc/qc_report.csv
  - tables/01-celltype_integration_clustering/03-qc/package_versions.txt
"""
from __future__ import annotations

import json
import os
import sys
import time
import gc
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

SEED = 42
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
STAGEA = ROOT / "epi-cm-core-workflow/h5ad/01-celltype_integration_clustering/03-qc/adata_qc_stageA.h5ad"
OUT_H5AD = ROOT / "epi-cm-core-workflow/h5ad/01-celltype_integration_clustering/03-qc/adata_qc.h5ad"
TAB = ROOT / "epi-cm-core-workflow/tables/01-celltype_integration_clustering/03-qc"
TAB.mkdir(parents=True, exist_ok=True)

# User-selected thresholds
MIN_CELLS_GENE = 3
MIN_GENES_FINAL = 200
MAX_PCT_MT_FINAL = 15.0
MAX_GENES_FINAL = 6000

# Scrublet params (defaults from rsc.pp.scrublet; we set random_state=SEED)
SCRUBLET_KW = dict(
    sim_doublet_ratio=2.0,
    expected_doublet_rate=0.05,
    n_prin_comps=30,
    random_state=SEED,
    verbose=False,
)


def run_scrublet_per_sample(adata: ad.AnnData) -> tuple[ad.AnnData, list[dict]]:
    """Run rsc.pp.scrublet on each sample subset, write scores back to adata.obs.

    Returns a list of per-sample records (for the parameters table).
    Samples that fail Scrublet on both rsc and CPU fallbacks are dropped from adata.
    """
    import rapids_singlecell as rsc

    adata.obs["doublet_score"] = np.nan
    adata.obs["predicted_doublet"] = False

    sample_records: list[dict] = []
    failed_samples: list[str] = []

    samples = adata.obs["sample"].astype(str).unique()
    print(f"[scrublet] {len(samples)} samples to process", flush=True)

    for i, sample_id in enumerate(samples, 1):
        mask = adata.obs["sample"].astype(str) == sample_id
        n_in = int(mask.sum())
        rec = {
            "sample": sample_id,
            "n_cells_in": n_in,
            "n_cells_out_after_scrublet": n_in,
            "n_predicted_doublets": 0,
            "backend": "rsc.pp.scrublet (GPU)",
            "status": "pending",
            "fallback_attempted": False,
            "failure_reason": "",
        }
        t0 = time.time()
        sub = adata[mask].copy()
        try:
            rsc.get.anndata_to_GPU(sub)
            rsc.pp.scrublet(sub, **SCRUBLET_KW)
            rsc.get.anndata_to_CPU(sub)
            # Sanity check
            if "doublet_score" not in sub.obs or "predicted_doublet" not in sub.obs:
                raise RuntimeError("scrublet did not write doublet_score/predicted_doublet")
            scores = sub.obs["doublet_score"].to_numpy()
            if not np.isfinite(scores).any():
                raise RuntimeError("scrublet produced all-NaN doublet_score")
            adata.obs.loc[sub.obs_names, "doublet_score"] = scores
            adata.obs.loc[sub.obs_names, "predicted_doublet"] = sub.obs["predicted_doublet"].astype(bool).to_numpy()
            n_pred = int(sub.obs["predicted_doublet"].sum())
            rec["n_predicted_doublets"] = n_pred
            rec["status"] = "ok"
        except Exception as e:
            print(f"[scrublet][{i}/{len(samples)}] {sample_id}: rsc failed: {e}; trying CPU fallback", flush=True)
            rec["fallback_attempted"] = True
            try:
                # CPU fallback via scanpy.scrublet (wrap with try/except)
                sc.pp.scrublet(sub, **SCRUBLET_KW)
                scores = sub.obs["doublet_score"].to_numpy()
                if not np.isfinite(scores).any():
                    raise RuntimeError("CPU scrublet produced all-NaN doublet_score")
                adata.obs.loc[sub.obs_names, "doublet_score"] = scores
                adata.obs.loc[sub.obs_names, "predicted_doublet"] = sub.obs["predicted_doublet"].astype(bool).to_numpy()
                n_pred = int(sub.obs["predicted_doublet"].sum())
                rec["n_predicted_doublets"] = n_pred
                rec["backend"] = "sc.pp.scrublet (CPU fallback)"
                rec["status"] = "ok"
            except Exception as e2:
                print(f"[scrublet][{i}/{len(samples)}] {sample_id}: BOTH rsc and CPU failed: {e2}", flush=True)
                rec["status"] = "failed"
                rec["failure_reason"] = f"rsc: {e}; cpu: {e2}"
                failed_samples.append(sample_id)
        finally:
            try:
                del sub
            except Exception:
                pass
            gc.collect()

        rec["elapsed_sec"] = round(time.time() - t0, 1)
        sample_records.append(rec)
        if i % 5 == 0 or i == len(samples):
            print(f"[scrublet] processed {i}/{len(samples)} samples", flush=True)

    # Drop failed samples entirely
    if failed_samples:
        keep = ~adata.obs["sample"].astype(str).isin(failed_samples)
        print(f"[scrublet] dropping {len(failed_samples)} failed samples: {failed_samples}", flush=True)
        adata = adata[keep].copy()
        # mark dropped in records
        for r in sample_records:
            if r["sample"] in failed_samples:
                r["note"] = "sample dropped from cleaned adata"
    return adata, sample_records, failed_samples


def main() -> None:
    t0 = time.time()
    print(f"[qc-C] reading {STAGEA.name}…", flush=True)
    adata = sc.read_h5ad(STAGEA)
    print(f"[qc-C] loaded {adata.shape} in {time.time()-t0:.1f}s", flush=True)
    n_start, v_start = adata.shape

    # ---- Step 2: Scrublet per sample ----
    print(f"[qc-C] === step 2-3: Scrublet per sample ===", flush=True)
    adata, scrub_records, failed_samples = run_scrublet_per_sample(adata)
    n_after_scrublet = adata.n_obs
    n_dropped_doublets = int(adata.obs["predicted_doublet"].fillna(False).astype(bool).sum())
    print(f"[qc-C] after scrublet: {n_after_scrublet} cells; {n_dropped_doublets} predicted_doublets marked", flush=True)

    # Apply doublet filter
    keep = ~adata.obs["predicted_doublet"].fillna(False).astype(bool)
    adata = adata[keep].copy()
    n_after_doublet = adata.n_obs
    print(f"[qc-C] after doublet filter: {n_after_doublet} cells (dropped {n_after_scrublet - n_after_doublet} doublets, {len(failed_samples)} failed samples dropped)", flush=True)

    # Write scrublet parameters table
    pd.DataFrame(scrub_records).to_csv(TAB / "02_scrublet_parameters.csv", index=False)
    pd.DataFrame([{
        "step": "03_doublet_filter",
        "input_h5ad": str(STAGEA),
        "output_h5ad_or_object": "in-memory adata",
        "n_obs_before": n_after_scrublet,
        "n_obs_after": n_after_doublet,
        "n_vars_before": v_start,
        "n_vars_after": adata.n_vars,
        "parameters": json.dumps({
            "filter_rule": "~predicted_doublet.astype(bool)",
            "n_failed_samples_dropped": len(failed_samples),
            "failed_samples": failed_samples,
        }),
        "backend": "scanpy",
        "code_file": "epi-cm-core-workflow/codes/01-celltype_integration_clustering/03-qc/03_scrublet_doublet.py",
        "random_seed": SEED,
        "notes": "predicted_doublet assigned per sample by rsc.pp.scrublet; samples failing both rsc and CPU are dropped",
    }]).to_csv(TAB / "03_doublet_filter_parameters.csv", index=False)

    # ---- Step 4-5: QC metrics (already computed in stageA, keep) ----
    # (qc metrics are inherited from stageA; adata.var['mt'], adata.var['ribo'] preserved)

    # ---- Step 6: gene filter (min_cells) ----
    n_pre_gene = adata.shape
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS_GENE)
    n_post_gene = adata.n_vars
    print(f"[qc-C] gene filter min_cells={MIN_CELLS_GENE}: {n_pre_gene[1]} -> {n_post_gene} genes", flush=True)
    pd.DataFrame([{
        "step": "05_gene_filter",
        "input_h5ad": "in-memory after doublet filter",
        "output_h5ad_or_object": "in-memory adata",
        "n_obs_before": n_pre_gene[0],
        "n_obs_after": n_pre_gene[0],
        "n_vars_before": n_pre_gene[1],
        "n_vars_after": n_post_gene,
        "parameters": json.dumps({"min_cells": MIN_CELLS_GENE}),
        "backend": "scanpy",
        "code_file": "epi-cm-core-workflow/codes/01-celltype_integration_clustering/03-qc/03_scrublet_doublet.py",
        "random_seed": SEED,
        "notes": "",
    }]).to_csv(TAB / "05_gene_filter_parameters.csv", index=False)

    # ---- Step 8: recalculate QC metrics after gene filtering ----
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    # ---- Step 9: final cell filter ----
    n_pre = adata.n_obs
    cell_mask = (
        (adata.obs["n_genes_by_counts"] >= MIN_GENES_FINAL) &
        (adata.obs["pct_counts_mt"] <= MAX_PCT_MT_FINAL) &
        (adata.obs["n_genes_by_counts"] <= MAX_GENES_FINAL)
    )
    adata = adata[cell_mask].copy()
    n_post = adata.n_obs
    print(f"[qc-C] final cell filter (min_genes>={MIN_GENES_FINAL}, max_pct_mt<={MAX_PCT_MT_FINAL}, max_n_genes<={MAX_GENES_FINAL}): {n_pre} -> {n_post}", flush=True)
    pd.DataFrame([{
        "step": "06_final_cell_filter",
        "input_h5ad": "in-memory after gene filter",
        "output_h5ad_or_object": str(OUT_H5AD),
        "n_obs_before": n_pre,
        "n_obs_after": n_post,
        "n_vars_before": adata.n_vars,
        "n_vars_after": adata.n_vars,
        "parameters": json.dumps({
            "min_genes_final": MIN_GENES_FINAL,
            "max_pct_mt_final": MAX_PCT_MT_FINAL,
            "max_genes_final": MAX_GENES_FINAL,
        }),
        "backend": "scanpy",
        "code_file": "epi-cm-core-workflow/codes/01-celltype_integration_clustering/03-qc/03_scrublet_doublet.py",
        "random_seed": SEED,
        "notes": "",
    }]).to_csv(TAB / "06_final_cell_filter_parameters.csv", index=False)

    # ---- Step 10: write adata_qc.h5ad ----
    print(f"[qc-C] writing {OUT_H5AD}…", flush=True)
    adata.write_h5ad(OUT_H5AD, compression="gzip")
    print(f"[qc-C] done. final shape={adata.shape}, total={time.time()-t0:.1f}s", flush=True)

    # QC summary
    summary = {
        "n_obs_start": n_start,
        "n_vars_start": v_start,
        "n_obs_after_initial_filter": n_start,  # stageA = after initial cell filter
        "n_obs_after_scrublet_drop_failures": n_after_scrublet,
        "n_obs_after_doublet_filter": n_after_doublet,
        "n_obs_final": n_post,
        "n_vars_final": adata.n_vars,
        "n_failed_samples_dropped": len(failed_samples),
        "failed_samples": failed_samples,
        "thresholds": {
            "initial_min_genes": 200,
            "min_cells_gene": MIN_CELLS_GENE,
            "min_genes_final": MIN_GENES_FINAL,
            "max_pct_mt_final": MAX_PCT_MT_FINAL,
            "max_genes_final": MAX_GENES_FINAL,
        },
    }
    with open(TAB / "qc_report.json", "w") as f:
        json.dump(summary, f, indent=2)
    # Also csv form
    pd.DataFrame([{
        "stage": "initial_cell_filter",
        "n_obs": n_start,
        "n_vars": v_start,
        "notes": "stageA: min_genes=200",
    }, {
        "stage": "scrublet_drop_failed_samples",
        "n_obs": n_after_scrublet,
        "n_vars": v_start,
        "notes": f"failed_samples={len(failed_samples)}",
    }, {
        "stage": "doublet_filter",
        "n_obs": n_after_doublet,
        "n_vars": v_start,
        "notes": "~predicted_doublet",
    }, {
        "stage": "gene_filter",
        "n_obs": n_after_doublet,
        "n_vars": n_post_gene,
        "notes": f"min_cells={MIN_CELLS_GENE}",
    }, {
        "stage": "final_cell_filter",
        "n_obs": n_post,
        "n_vars": adata.n_vars,
        "notes": f"min_genes={MIN_GENES_FINAL}, max_pct_mt={MAX_PCT_MT_FINAL}, max_genes={MAX_GENES_FINAL}",
    }]).to_csv(TAB / "qc_report.csv", index=False)
    print(f"[qc-C] qc_report.csv written", flush=True)

    # Per-sample post-QC summary
    grp = adata.obs.groupby("sample", observed=True)
    rows = []
    for sid, sub in grp:
        rows.append({
            "sample": sid,
            "series": sub["series"].iloc[0],
            "status": sub["status"].iloc[0],
            "n_cells": sub.shape[0],
            "median_total_counts": float(sub["total_counts"].median()),
            "median_n_genes": float(sub["n_genes_by_counts"].median()),
            "median_pct_mt": float(sub["pct_counts_mt"].median()),
        })
    pd.DataFrame(rows).sort_values("n_cells", ascending=False).to_csv(TAB / "qc_report_per_sample_post.csv", index=False)

    # package_versions
    import rapids_singlecell as rsc
    pkg = TAB / "package_versions.txt"
    with open(pkg, "w") as f:
        f.write("QC stage C — Scrublet + final filters\n")
        f.write(f"python: {sys.version.split()[0]}\n")
        f.write(f"scanpy: {sc.__version__}\n")
        f.write(f"rapids_singlecell: {rsc.__version__}\n")
        f.write(f"numpy: {np.__version__}\n")
        f.write(f"pandas: {pd.__version__}\n")
        f.write(f"random_seed: {SEED}\n")
        f.write(f"scrublet_kwargs: {json.dumps(SCRUBLET_KW)}\n")
        f.write(f"thresholds: {json.dumps(summary['thresholds'])}\n")
    print(f"[qc-C] package_versions.txt written", flush=True)


if __name__ == "__main__":
    main()