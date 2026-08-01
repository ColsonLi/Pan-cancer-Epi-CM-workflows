"""Block 03 step 4: Epi-CM Association - Spearman branch (per SKILL.md).

Strictly canonical:
  - Run all pairwise (epi_subtype x CM) combinations
  - Use multipletests(fdr_bh) for q values
  - For each status context: tumor, normal-like, and overall
  - Save to a Spearman-specific task subdirectory

Outputs (under tables/03-epi-cm-discovery/association-spearman/):
  epi_cm_association_overall_rho_matrix.csv
  epi_cm_association_overall_q_matrix.csv
  epi_cm_association_tumor_rho_matrix.csv
  epi_cm_association_tumor_q_matrix.csv
  epi_cm_association_normal-like_rho_matrix.csv
  epi_cm_association_normal-like_q_matrix.csv
  epi_cm_association_overall_long.csv
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

SEED = 42
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
TAB = ROOT / "epi-cm-core-workflow/tables/03-epi-cm-discovery"
EPI_FREQ = TAB / "epi_subtype_frequency.csv"
W_CANON = TAB / "W_df.csv"
SAMPLE_INCL = TAB / "sample_inclusion_exclusion.csv"
SAMPLE_STATUS = TAB / "sample_status.csv"
OUT = TAB / "association-spearman"
OUT.mkdir(parents=True, exist_ok=True)


def epi_cm_correlations(epi_freq: pd.DataFrame, cm_activity: pd.DataFrame, method: str = "spearman"):
    common = epi_freq.index.intersection(cm_activity.index)
    E = epi_freq.loc[common].astype(float)
    C = cm_activity.loc[common].astype(float)
    corr = pd.DataFrame(index=E.columns, columns=C.columns, dtype=float)
    pval = pd.DataFrame(index=E.columns, columns=C.columns, dtype=float)
    test = spearmanr if method == "spearman" else None
    for epi in E.columns:
        for cm in C.columns:
            x, y = E[epi], C[cm]
            ok = x.notna() & y.notna()
            if ok.sum() < 3 or x[ok].nunique() < 2 or y[ok].nunique() < 2:
                corr.loc[epi, cm] = np.nan
                pval.loc[epi, cm] = np.nan
                continue
            r, p = test(x[ok], y[ok])
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
    # Use keep_for_epi_cm samples for overall / per-status
    keep_samples = incl.index[incl["keep_for_epi_cm"]].astype(str).tolist()
    keep_samples = [s for s in keep_samples if s in epi_freq.index and s in W.index and s in status.index]
    print(f"[assoc-sp] n_samples (keep_for_epi_cm): {len(keep_samples)}", flush=True)

    contexts = {
        "overall": keep_samples,
        "tumor": [s for s in keep_samples if status.loc[s] == "tumor"],
        "normal-like": [s for s in keep_samples if status.loc[s] == "normal-like"],
    }
    for ctx, samples in contexts.items():
        epi_c = epi_freq.loc[samples]
        W_c = W.loc[samples]
        # CM activity is the W columns (no status column for Epi-CM assoc)
        # activity_df_sample_by_CM has 'status' and 'non_epi_cells' as extra cols
        cm_act = W_c
        rho, pval, qval = epi_cm_correlations(epi_c, cm_act, method="spearman")
        rho.to_csv(OUT / f"epi_cm_association_{ctx}_rho_matrix.csv")
        qval.to_csv(OUT / f"epi_cm_association_{ctx}_q_matrix.csv")
        # long format
        rows = []
        for e in rho.index:
            for c in rho.columns:
                rows.append({
                    "epi_subtype": e, "CM": c,
                    "rho": float(rho.loc[e, c]),
                    "pval": float(pval.loc[e, c]),
                    "qval": float(qval.loc[e, c]),
                    "n_samples": len(samples),
                })
        ld = pd.DataFrame(rows)
        ld.to_csv(OUT / f"epi_cm_association_{ctx}_long.csv", index=False)
        sig = ld[ld["qval"] < 0.05].sort_values("rho", key=abs, ascending=False)
        sig.to_csv(OUT / f"epi_cm_association_{ctx}_significant_q0.05.csv", index=False)
        print(f"[assoc-sp] {ctx}: {len(sig)}/{len(ld)} sig (q<0.05)", flush=True)

    # Summary
    summary = {
        "n_samples": len(keep_samples),
        "n_tumor": len(contexts["tumor"]),
        "n_normal_like": len(contexts["normal-like"]),
        "n_epi_subtypes": epi_freq.shape[1],
        "n_cm": W.shape[1],
    }
    with open(OUT / "association_spearman_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[assoc-sp] done. total={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()