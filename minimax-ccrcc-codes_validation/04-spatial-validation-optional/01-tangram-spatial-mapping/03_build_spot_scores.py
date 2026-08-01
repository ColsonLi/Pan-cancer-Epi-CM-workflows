"""Block 04 step 3: per-sample spot-score CSV."""
from __future__ import annotations

import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

SEED = 42
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
H5AD_DIR = ROOT / "epi-cm-core-workflow/h5ad/04-spatial-validation-optional/01-tangram-spatial-mapping"
TAB = ROOT / "epi-cm-core-workflow/tables/04-spatial-validation-optional/01-tangram-spatial-mapping"
H_CANON = ROOT / "epi-cm-core-workflow/tables/03-epi-cm-discovery/H_df.csv"


def save_spot_scores(sample, obs_df, cm_df, epi_frac_df):
    score = obs_df.copy()
    score.insert(0, "sample", sample)
    score.insert(1, "spot_id", score.index.astype(str))
    if "barcode" not in score.columns:
        score["barcode"] = score["spot_id"].astype(str)
    if "array_row" not in score.columns and "spatial_y" in score.columns:
        score["array_row"] = score["spatial_y"]
    if "array_col" not in score.columns and "spatial_x" in score.columns:
        score["array_col"] = score["spatial_x"]
    if "spatial_x" not in score.columns and "array_col" in score.columns:
        score["spatial_x"] = score["array_col"]
    if "spatial_y" not in score.columns and "array_row" in score.columns:
        score["spatial_y"] = score["array_row"]
    score = score.join(cm_df.add_prefix("CMact__"), how="left")
    score = score.join(epi_frac_df.add_prefix("EPIfrac__"), how="left")
    required = ["sample", "spot_id", "array_row", "array_col"]
    missing = [c for c in required if c not in score.columns]
    if missing:
        raise KeyError(f"{sample}: spot-score CSV missing required columns: {missing}")
    out = TAB / f"{sample}_tangram_pseudobulk_epi_cm_spot_scores.csv"
    score.to_csv(out, index=False)
    return out


def main():
    t0 = time.time()
    H = pd.read_csv(H_CANON, index_col=0)
    canonical_cms = H.index.astype(str).tolist()
    print(f"[fix] H: {H.shape}; canonical_cms={len(canonical_cms)}", flush=True)
    files = sorted(H5AD_DIR.glob("*_tangram_mapped.h5ad"))
    print(f"[fix] {len(files)} mapped h5ad files", flush=True)
    n_rows_log = []
    for fpath in files:
        sample = fpath.stem.replace("_tangram_mapped", "")
        print(f"\n[fix] === {sample} ===", flush=True)
        asp = ad.read_h5ad(fpath)
        ct = asp.obsm["tangram_ct_pred"]
        print(f"[fix] ct shape: {ct.shape}", flush=True)
        row_sums = ct.sum(axis=1).replace(0, np.nan)
        ab_norm = ct.div(row_sums, axis=0).fillna(0.0)
        epi_subs = [c for c in ab_norm.columns if str(c).startswith("Epi_")]
        if not epi_subs:
            raise ValueError(f"{sample}: no Epi_* columns in tangram_ct_pred")
        epi_abundance = ab_norm.loc[:, epi_subs].astype(float)
        epi_total = epi_abundance.sum(axis=1)
        epi_frac = epi_abundance.div(epi_total.replace(0, np.nan), axis=0).fillna(0.0)
        positive = epi_total.to_numpy() > 0
        if positive.any():
            row_sums_check = epi_frac.loc[positive].sum(axis=1).to_numpy()
            if not np.allclose(row_sums_check, 1.0, rtol=1e-6, atol=1e-8):
                raise ValueError(f"{sample}: EPIfrac row-normalization failed")
        cm_act = pd.DataFrame(index=ab_norm.index, columns=canonical_cms, dtype=float)
        for cm in canonical_cms:
            h_row = H.loc[cm].reindex(ab_norm.columns).fillna(0.0)
            cm_act[cm] = (ab_norm * h_row).sum(axis=1)
        obs_df = asp.obs[["array_row", "array_col"]].copy()
        if "spatial" in asp.obsm:
            obs_df["spatial_x"] = asp.obsm["spatial"][:, 0]
            obs_df["spatial_y"] = asp.obsm["spatial"][:, 1]
        else:
            obs_df["spatial_x"] = obs_df["array_col"]
            obs_df["spatial_y"] = obs_df["array_row"]
        out = save_spot_scores(sample, obs_df, cm_act, epi_frac)
        n_rows_log.append((sample, len(obs_df), out.stat().st_size))
        print(f"[fix] {sample}: spot-score written -> {out.name} ({len(obs_df)} spots)", flush=True)
    print(f"\n[fix] done. total={time.time()-t0:.1f}s", flush=True)
    for s, n, sz in n_rows_log:
        print(f"  {s}: {n} spots, {sz/1024:.1f} KB")


if __name__ == "__main__":
    main()