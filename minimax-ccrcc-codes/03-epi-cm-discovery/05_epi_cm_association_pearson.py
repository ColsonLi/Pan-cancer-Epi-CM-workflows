"""Block 03 step 5: Epi-CM Association - Pearson branch (per SKILL.md).

Separate branch from Spearman, parallel output inventory.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests

SEED = 42
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
TAB = ROOT / "epi-cm-core-workflow/tables/03-epi-cm-discovery"
EPI_FREQ = TAB / "epi_subtype_frequency.csv"
W_CANON = TAB / "W_df.csv"
SAMPLE_INCL = TAB / "sample_inclusion_exclusion.csv"
SAMPLE_STATUS = TAB / "sample_status.csv"
OUT = TAB / "association-pearson"
OUT.mkdir(parents=True, exist_ok=True)


def epi_cm_correlations_pearson(epi_freq: pd.DataFrame, cm_activity: pd.DataFrame):
    common = epi_freq.index.intersection(cm_activity.index)
    E = epi_freq.loc[common].astype(float)
    C = cm_activity.loc[common].astype(float)
    corr = pd.DataFrame(index=E.columns, columns=C.columns, dtype=float)
    pval = pd.DataFrame(index=E.columns, columns=C.columns, dtype=float)
    for epi in E.columns:
        for cm in C.columns:
            x, y = E[epi], C[cm]
            ok = x.notna() & y.notna()
            if ok.sum() < 3 or x[ok].nunique() < 2 or y[ok].nunique() < 2:
                corr.loc[epi, cm] = np.nan
                pval.loc[epi, cm] = np.nan
                continue
            r, p = pearsonr(x[ok], y[ok])
            corr.loc[epi, cm] = r
            pval.loc[epi, cm] = p
    flat = pval.to_numpy().ravel()
    valid = np.isfinite(flat)
    qflat = np.full(flat.shape, np.nan, dtype=float)
    qflat[valid] = multipletests(flat[valid], method="fdr_bh")[1]
    qval = pd.DataFrame(qflat.reshape(pval.shape), index=pval.index, columns=pval.columns)
    return corr, pval, qval


def main() -> None:
    t0 = time.time()
    epi_freq = pd.read_csv(EPI_FREQ, index_col=0)
    W = pd.read_csv(W_CANON, index_col=0)
    status = pd.read_csv(SAMPLE_STATUS, index_col=0)["status"]
    incl = pd.read_csv(SAMPLE_INCL, index_col=0)
    keep_samples = incl.index[incl["keep_for_epi_cm"]].astype(str).tolist()
    keep_samples = [s for s in keep_samples if s in epi_freq.index and s in W.index and s in status.index]

    contexts = {
        "overall": keep_samples,
        "tumor": [s for s in keep_samples if status.loc[s] == "tumor"],
        "normal-like": [s for s in keep_samples if status.loc[s] == "normal-like"],
    }
    for ctx, samples in contexts.items():
        epi_c = epi_freq.loc[samples]
        W_c = W.loc[samples]
        rmat, pmat, qmat = epi_cm_correlations_pearson(epi_c, W_c)
        rmat.to_csv(OUT / f"epi_cm_association_{ctx}_pearson_r_matrix.csv")
        qmat.to_csv(OUT / f"epi_cm_association_{ctx}_pearson_q_matrix.csv")
        rows = []
        for e in rmat.index:
            for c in rmat.columns:
                rows.append({
                    "epi_subtype": e, "CM": c,
                    "r": float(rmat.loc[e, c]),
                    "pval": float(pmat.loc[e, c]),
                    "qval": float(qmat.loc[e, c]),
                    "n_samples": len(samples),
                })
        pd.DataFrame(rows).to_csv(OUT / f"epi_cm_association_{ctx}_pearson_long.csv", index=False)

    with open(OUT / "association_pearson_summary.json", "w") as f:
        json.dump({"n_samples": len(keep_samples),
                   "n_tumor": len(contexts["tumor"]),
                   "n_normal_like": len(contexts["normal-like"]),
                   "n_epi_subtypes": epi_freq.shape[1],
                   "n_cm": W.shape[1]}, f, indent=2)
    print(f"[assoc-pe] done. total={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()