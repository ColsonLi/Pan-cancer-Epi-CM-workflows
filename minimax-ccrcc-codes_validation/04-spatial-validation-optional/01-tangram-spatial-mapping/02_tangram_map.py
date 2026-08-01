"""Block 04 step 2: per-sample canonical Tangram mapping.

Per SKILL.md (Module 05 Method Integrity Rules):
  Required Tangram preflight gate before any spatial mapping.
  Canonical Tangram call (FIXED):
      tg.pp_adatas(asc, asp, genes=common)
      tg.map_cells_to_space(mode='cells', device='cuda:0', num_epochs=350,
                            learning_rate=0.05, random_state=42, verbose=False)
      tg.project_cell_annotations(mapper, asp, annotation='cell_subtype')

Inputs (produced by step 1):
  - gene_intersection_list.csv     (Tangram gene set, 2488 genes)
  - subtype_mean_reference.csv     (66 subtypes x 2488-gene-set mean)
  - sample_scope.csv

Outputs:
  - h5ad/04-spatial-validation-optional/01-tangram-spatial-mapping/<sample>_tangram_mapped.h5ad
  - tables/.../tangram_run.log
"""
from __future__ import annotations

import gc
import importlib.util
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import torch

SEED = 42
np.random.seed(SEED)
torch.manual_seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
SP_DIR = ROOT / "spatial_h5ad"
BLOCK = ROOT / "epi-cm-core-workflow"
H5AD_OUT = BLOCK / "h5ad/04-spatial-validation-optional/01-tangram-spatial-mapping"
TAB = BLOCK / "tables/04-spatial-validation-optional/01-tangram-spatial-mapping"
LOG = TAB / "tangram_run.log"
H5AD_OUT.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

if importlib.util.find_spec("tangram") is None:
    raise RuntimeError(
        "Tangram is required for canonical spatial validation. "
        "Install tangram-sc/PyTorch and rerun; do not substitute score_genes, "
        "marker-based scores, NNLS, regression, or nearest-neighbor matching."
    )
import tangram as tg
print(f"[tangram] version {tg.__version__} loaded", flush=True)


def normalize_spatial_sample(adata, target_sum=1e4):
    adata.var_names_make_unique()
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    else:
        adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    return adata


def standardize_spatial_obs_names(adata, sample):
    if "barcode" not in adata.obs.columns:
        adata.obs["barcode"] = adata.obs_names.astype(str)
    adata.obs["sample"] = sample
    adata.obs_names = [f"{b}_{sample}" for b in adata.obs["barcode"].astype(str)]
    if not adata.obs_names.is_unique:
        raise ValueError(f"Non-unique obs_names: {sample}")
    return adata


def main():
    t0 = time.time()
    log_lines = []

    def log(msg):
        print(msg, flush=True)
        log_lines.append(msg)

    common_genes = pd.read_csv(TAB / "gene_intersection_list.csv")["gene"].astype(str).tolist()
    log(f"[tangram] common_genes (markers ∩ ref ∩ all spatial): {len(common_genes)}")

    log("[ref] reloading subtype_mean_reference from step 1 output…")
    subtype_mean = pd.read_csv(TAB / "subtype_mean_reference.csv", index_col=0)
    subtype_mean.index = subtype_mean.index.astype(str)
    log(f"[ref] subtype_mean: {subtype_mean.shape}")
    common = sorted(set(common_genes) & set(subtype_mean.columns))
    log(f"[tangram] canonical Tangram gene set: {len(common)}")

    sp_files = sorted(SP_DIR.glob("*.h5ad"))
    log(f"[tangram] {len(sp_files)} spatial samples to map")

    for fpath in sp_files:
        sample = fpath.stem
        log(f"\n[tangram] === {sample} ===")
        adata_sp = sc.read_h5ad(fpath)
        adata_sp = standardize_spatial_obs_names(adata_sp, sample)
        adata_sp = normalize_spatial_sample(adata_sp)
        sp_genes = sorted(set(common) & set(adata_sp.var_names))
        log(f"[tangram] {sample}: per-sample intersect = {len(sp_genes)}")
        asc = ad.AnnData(
            X=subtype_mean[sp_genes].values.astype(np.float32),
            obs=pd.DataFrame(index=subtype_mean.index),
            var=pd.DataFrame(index=sp_genes),
        )
        asc.obs["cell_subtype"] = asc.obs.index.astype(str)
        asp = adata_sp[:, sp_genes].copy()
        try:
            tg.pp_adatas(asc, asp, genes=sp_genes)
            mapper = tg.map_cells_to_space(
                asc, asp, mode="cells", device="cuda:0",
                num_epochs=350, learning_rate=0.05,
                random_state=SEED, verbose=False,
            )
            tg.project_cell_annotations(mapper, asp, annotation="cell_subtype")
        except Exception as e:
            log(f"[tangram] {sample} failed: {e}")
            raise
        if "tangram_ct_pred" not in asp.obsm:
            raise RuntimeError(f"{sample}: tangram_ct_pred not in asp.obsm")
        ct = asp.obsm["tangram_ct_pred"]
        log(f"[tangram] {sample}: tangram_ct_pred shape = {ct.shape}")
        asp.write_h5ad(H5AD_OUT / f"{sample}_tangram_mapped.h5ad", compression="gzip")
        log(f"[tangram] {sample} saved h5ad")
        del mapper, asc, asp, adata_sp
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

    log(f"\n[tangram] all done. total={time.time()-t0:.1f}s")
    LOG.write_text("\n".join(log_lines) + "\n")


if __name__ == "__main__":
    main()