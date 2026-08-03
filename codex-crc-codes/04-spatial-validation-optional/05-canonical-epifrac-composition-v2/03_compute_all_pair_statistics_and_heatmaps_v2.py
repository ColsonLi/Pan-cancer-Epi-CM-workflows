#!/usr/bin/env python3
"""Compute complete CSV-driven Spearman/Fisher/Stouffer spatial statistics."""

from __future__ import annotations

import json
import math
import platform
import random
import re
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
from scipy import stats

SEED = 42
random.seed(SEED)
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/crc_codex/epi-cm-core-workflow")
TASK = "05-canonical-epifrac-composition-v2"
BLOCK = "04-spatial-validation-optional-sig_genes"
BASE_TABLE = ROOT / f"tables/{BLOCK}/{TASK}"
MAPPING_TABLES = BASE_TABLE / "02-tangram-mapping"
SCORE_DIR = MAPPING_TABLES / "spot_scores"
TABLE_DIR = BASE_TABLE / "03-all-pair-statistics-and-plots"
FIG_DIR = ROOT / f"figures/{BLOCK}/{TASK}/03-all-pair-statistics-and-plots/stouffer_heatmaps"
SAMPLE_SCOPE = BASE_TABLE / "01-reference-and-manifest/spatial_sample_scope_manifest.csv"

plt.rcParams.update(
    {
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "font.family": "DejaVu Sans",
        "font.size": 8,
        "axes.titlesize": 9,
    }
)


def bh(pvalues: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvalues, errors="coerce").to_numpy(dtype=float)
    out = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    values = p[valid]
    if len(values):
        order = np.argsort(values)
        ranked = values[order]
        adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
        adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
        restored = np.empty_like(adjusted)
        restored[order] = np.minimum(adjusted, 1.0)
        out[valid] = restored
    return pd.Series(out, index=pvalues.index)


def signed_z_from_p(p: float, direction: float) -> float:
    if not np.isfinite(p) or not np.isfinite(direction) or direction == 0:
        return np.nan
    p = max(min(float(p), 1.0), np.finfo(float).tiny)
    return float(np.sign(direction) * stats.norm.isf(p / 2.0))


def percentile_rank(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.rank(method="average", pct=True)


def sample_statistics(path: Path, expected_cms: list[str], expected_epis: list[str]) -> tuple[list[dict], list[dict]]:
    score = pd.read_csv(path)
    sample_values = score["sample"].astype(str).unique()
    if len(sample_values) != 1:
        raise ValueError(f"{path}: expected one sample")
    sample = sample_values[0]
    status = str(score["status"].iloc[0])
    cm_cols = [f"CMact__{x}" for x in expected_cms]
    epi_cols = [f"EPIfrac__{x}" for x in expected_epis]
    features = cm_cols + epi_cols
    x = score[features].apply(pd.to_numeric, errors="coerce")
    corr, pmat = stats.spearmanr(x.to_numpy(dtype=float), axis=0, nan_policy="omit")
    pct = percentile_rank(x)
    spearman_rows = []
    fisher_rows = []
    for cm_i, cm in enumerate(expected_cms):
        cm_col = f"CMact__{cm}"
        cm_high = pct[cm_col].ge(0.5)
        for epi_i, epi in enumerate(expected_epis):
            epi_col = f"EPIfrac__{epi}"
            epi_high = pct[epi_col].ge(0.5)
            j = len(expected_cms) + epi_i
            rho = float(corr[cm_i, j])
            p_s = float(pmat[cm_i, j])
            spearman_rows.append(
                {
                    "sample": sample,
                    "status": status,
                    "CM": cm,
                    "epi_subtype": epi,
                    "n_spatial_cells": len(score),
                    "spearman_rho": rho,
                    "spearman_p": p_s,
                }
            )
            both = int((cm_high & epi_high).sum())
            cm_only = int((cm_high & ~epi_high).sum())
            epi_only = int((~cm_high & epi_high).sum())
            both_low = int((~cm_high & ~epi_high).sum())
            odds, p_f = stats.fisher_exact([[both, cm_only], [epi_only, both_low]])
            corrected_or = ((both + 0.5) * (both_low + 0.5)) / ((cm_only + 0.5) * (epi_only + 0.5))
            log2_or = float(np.log2(corrected_or))
            fisher_rows.append(
                {
                    "sample": sample,
                    "status": status,
                    "CM": cm,
                    "epi_subtype": epi,
                    "percentile_cutoff": 0.5,
                    "both_high": both,
                    "cm_high_only": cm_only,
                    "epi_high_only": epi_only,
                    "both_low": both_low,
                    "fisher_odds_ratio": float(odds),
                    "haldane_log2_odds_ratio": log2_or,
                    "fisher_p": float(p_f),
                    "signed_z": signed_z_from_p(float(p_f), log2_or),
                }
            )
    return spearman_rows, fisher_rows


def combine_stouffer(fisher: pd.DataFrame, scope: str) -> pd.DataFrame:
    rows = []
    for (cm, epi), group in fisher.groupby(["CM", "epi_subtype"], sort=False):
        z = pd.to_numeric(group["signed_z"], errors="coerce").dropna().to_numpy(dtype=float)
        if len(z):
            combined = float(z.sum() / math.sqrt(len(z)))
            p = float(2.0 * stats.norm.sf(abs(combined)))
        else:
            combined, p = np.nan, np.nan
        rows.append(
            {
                "scope": scope,
                "CM": cm,
                "epi_subtype": epi,
                "n_samples": len(z),
                "combined_signed_z": combined,
                "p_value": p,
            }
        )
    out = pd.DataFrame(rows)
    out["q_value_bh"] = bh(out["p_value"])
    return out


def cm_sort_key(value: str):
    match = re.search(r"CM(\d+)", str(value))
    return int(match.group(1)) if match else 999


def q_label(q: float) -> str:
    if not np.isfinite(q):
        return "ns"
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return "ns"


def plot_heatmap(table: pd.DataFrame, stem: Path, title: str) -> None:
    z = table.pivot(index="epi_subtype", columns="CM", values="combined_signed_z")
    q = table.pivot(index="epi_subtype", columns="CM", values="q_value_bh")
    rows = sorted(z.index)
    cols = sorted(z.columns, key=cm_sort_key)
    z, q = z.reindex(index=rows, columns=cols), q.reindex(index=rows, columns=cols)
    annot = q.map(q_label)
    finite = z.to_numpy()[np.isfinite(z.to_numpy())]
    vmax = float(np.max(np.abs(finite))) if len(finite) else 1.0
    vmax = vmax if vmax > 0 else 1.0
    fig, ax = plt.subplots(figsize=(max(8.0, 0.52 * len(cols) + 2.6), max(5.0, 0.42 * len(rows) + 1.9)))
    sns.heatmap(
        z,
        cmap="coolwarm",
        center=0,
        vmin=-vmax,
        vmax=vmax,
        linewidths=0.25,
        linecolor="white",
        annot=annot,
        fmt="",
        annot_kws={"fontsize": 6.5, "color": "black"},
        square=True,
        cbar_kws={"label": "Combined signed Stouffer Z"},
        ax=ax,
    )
    ax.set_title(title, fontweight="bold")
    ax.set_xlabel("CM")
    ax.set_ylabel("Epithelial subtype")
    ax.tick_params(axis="x", rotation=45)
    ax.tick_params(axis="y", rotation=0)
    fig.tight_layout()
    fig.savefig(stem.with_suffix(".pdf"), dpi=300, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".svg"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    started = time.time()
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    scope = pd.read_csv(SAMPLE_SCOPE)
    cm_names = pd.read_csv(MAPPING_TABLES / "all_sample_cm_epi_pair_manifest.csv")["CM"].drop_duplicates().tolist()
    epi_names = pd.read_csv(MAPPING_TABLES / "all_sample_cm_epi_pair_manifest.csv")["epi_subtype"].drop_duplicates().tolist()
    score_paths = sorted(SCORE_DIR.glob("*_tangram_pseudobulk_epi_cm_spot_scores.csv"))
    if len(score_paths) != 37:
        raise ValueError(f"Expected 37 score tables, found {len(score_paths)}")
    all_spearman, all_fisher = [], []
    for i, path in enumerate(score_paths, 1):
        spearman, fisher = sample_statistics(path, cm_names, epi_names)
        all_spearman.extend(spearman)
        all_fisher.extend(fisher)
        print(f"DONE statistics {i}/37 {path.name}", flush=True)
    spearman = pd.DataFrame(all_spearman)
    fisher = pd.DataFrame(all_fisher)
    spearman["spearman_q_within_sample"] = spearman.groupby("sample")["spearman_p"].transform(bh)
    fisher["fisher_q_within_sample"] = fisher.groupby("sample")["fisher_p"].transform(bh)
    spearman.to_csv(TABLE_DIR / "per_sample_spearman.csv", index=False)
    fisher.to_csv(TABLE_DIR / "percentile_quadrant_fisher_per_sample.csv", index=False)

    all_samples = combine_stouffer(fisher, "all-samples")
    tumor_samples = set(scope.loc[scope["include_tumor_only"].astype(bool), "sample"].astype(str))
    tumor_only = combine_stouffer(fisher[fisher["sample"].isin(tumor_samples)], "tumor-only")
    all_path = TABLE_DIR / "statistics/all-samples/percentile_quadrant_fisher_sample_stouffer_all_samples.csv"
    tumor_path = TABLE_DIR / "statistics/tumor-only/percentile_quadrant_fisher_sample_stouffer_tumor_only.csv"
    all_path.parent.mkdir(parents=True, exist_ok=True)
    tumor_path.parent.mkdir(parents=True, exist_ok=True)
    all_samples.to_csv(all_path, index=False)
    tumor_only.to_csv(tumor_path, index=False)
    all_samples.to_csv(TABLE_DIR / "percentile_quadrant_fisher_sample_stouffer_all_samples.csv", index=False)
    tumor_only.to_csv(TABLE_DIR / "percentile_quadrant_fisher_sample_stouffer_tumor_only.csv", index=False)
    plot_heatmap(all_samples, FIG_DIR / "all_samples_stouffer_signedZ_qstars", "All spatial samples")
    plot_heatmap(tumor_only, FIG_DIR / "tumor_only_stouffer_signedZ_qstars", "Tumor spatial samples")

    summary = {
        "n_samples": int(fisher["sample"].nunique()),
        "n_tumor_samples": len(tumor_samples),
        "n_cm": len(cm_names),
        "n_epi_subtypes": len(epi_names),
        "n_pair_sample_rows": len(fisher),
        "n_stouffer_pairs_all_samples": len(all_samples),
        "n_stouffer_pairs_tumor_only": len(tumor_only),
        "stouffer_weighting": "unweighted sample-level signed Z; spatial samples are replicates",
        "elapsed_seconds": round(time.time() - started, 3),
        "python": platform.python_version(),
        "scipy": scipy.__version__,
        "seed": SEED,
    }
    (TABLE_DIR / "statistics_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
