#!/usr/bin/env python3
"""Render every canonical CM-Epi pair for one spatial sample from its score CSV."""

from __future__ import annotations

import argparse
import json
import random
import time
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

SEED = 42
random.seed(SEED)
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/crc_codex/epi-cm-core-workflow")
TASK = "05-canonical-epifrac-composition-v2"
BLOCK = "04-spatial-validation-optional-sig_genes"
SCORE_DIR = ROOT / f"tables/{BLOCK}/{TASK}/02-tangram-mapping/spot_scores"
PAIR_MANIFEST = ROOT / f"tables/{BLOCK}/{TASK}/02-tangram-mapping/all_sample_cm_epi_pair_manifest.csv"
TABLE_DIR = ROOT / f"tables/{BLOCK}/{TASK}/03-all-pair-statistics-and-plots/plot_manifests_by_sample"
FIG_ROOT = ROOT / f"figures/{BLOCK}/{TASK}/03-all-pair-statistics-and-plots"

plt.rcParams.update(
    {
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "font.family": "DejaVu Sans",
    }
)

CM_CMAP = LinearSegmentedColormap.from_list(
    "cm_percentile_light_to_blue", ["#f2f2f2", "#b7d7ea", "#4c78a8", "#2b6cb0"]
)
EPI_CMAP = LinearSegmentedColormap.from_list(
    "epi_percentile_light_to_orange", ["#f2f2f2", "#fde6b3", "#fbbf24", "#f59e0b"]
)
CATEGORY_ORDER = ["both_low", "epi_high_only", "cm_high_only", "both_high"]
CATEGORY_LABELS = {
    "both_low": "CM<0.5, Epi<0.5",
    "epi_high_only": "Epi>=0.5 only",
    "cm_high_only": "CM>=0.5 only",
    "both_high": "Both>=0.5",
}
CATEGORY_COLORS = {
    "both_low": "#d9d9d9",
    "epi_high_only": "#f59e0b",
    "cm_high_only": "#2b6cb0",
    "both_high": "#b91c1c",
}


def args():
    p = argparse.ArgumentParser()
    p.add_argument("--sample", required=True)
    return p.parse_args()


def robust_limits(values):
    arr = pd.to_numeric(values, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna().to_numpy()
    if not len(arr):
        return 0.0, 1.0
    lo, hi = np.nanquantile(arr, [0.01, 0.99])
    if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
        lo, hi = float(np.nanmin(arr)), float(np.nanmax(arr))
    return (float(lo), float(hi if hi != lo else lo + 1e-9))


def set_axes(ax):
    ax.invert_yaxis()
    ax.set_box_aspect(1)
    ax.set_xlabel("array_col")
    ax.set_ylabel("array_row")


def save_pdf_svg(fig, stem: Path):
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", dpi=300)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", dpi=300)
    plt.close(fig)


def raw_plot(score, sample, cm, epi, out_dir, size):
    cm_raw = score[f"CMact__{cm}"]
    epi_raw = score[f"EPIfrac__{epi}"]
    x, y = score["array_col"].to_numpy(), score["array_row"].to_numpy()
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6), constrained_layout=True)
    a = axes[0].scatter(x, y, c=cm_raw, s=size, cmap="magma", norm=Normalize(*robust_limits(cm_raw)), linewidths=0, rasterized=True)
    set_axes(axes[0]); axes[0].set_title(f"CMact: {cm}")
    fig.colorbar(a, ax=axes[0], fraction=0.046, pad=0.02).set_label("CMact")
    b = axes[1].scatter(x, y, c=epi_raw, s=size, cmap="viridis", norm=Normalize(*robust_limits(epi_raw)), linewidths=0, rasterized=True)
    set_axes(axes[1]); axes[1].set_title(f"EPIfrac: {epi}")
    fig.colorbar(b, ax=axes[1], fraction=0.046, pad=0.02).set_label("EPIfrac")
    fig.suptitle(f"{sample}: {cm} | {epi}", y=1.04, fontsize=9, fontweight="bold")
    save_pdf_svg(fig, out_dir / f"{sample}__{cm}__{epi}__raw_cmact_epifrac")


def pct_plot(score, pct, sample, cm, epi, out_dir, size):
    cm_pct, epi_pct = pct[f"CMact__{cm}"], pct[f"EPIfrac__{epi}"]
    cm_high, epi_high = cm_pct.ge(0.5), epi_pct.ge(0.5)
    category = pd.Series("both_low", index=score.index)
    category.loc[epi_high & ~cm_high] = "epi_high_only"
    category.loc[cm_high & ~epi_high] = "cm_high_only"
    category.loc[cm_high & epi_high] = "both_high"
    x, y = score["array_col"].to_numpy(), score["array_row"].to_numpy()
    fig, axes = plt.subplots(1, 3, figsize=(14.8, 4.8), constrained_layout=True)
    a = axes[0].scatter(x, y, c=cm_pct, cmap=CM_CMAP, vmin=0, vmax=1, s=size, linewidths=0, rasterized=True)
    set_axes(axes[0]); axes[0].set_title(f"CM percentile: {cm}")
    fig.colorbar(a, ax=axes[0], fraction=0.046, pad=0.02).set_label("CMact percentile")
    b = axes[1].scatter(x, y, c=epi_pct, cmap=EPI_CMAP, vmin=0, vmax=1, s=size, linewidths=0, rasterized=True)
    set_axes(axes[1]); axes[1].set_title(f"Epi percentile: {epi}")
    fig.colorbar(b, ax=axes[1], fraction=0.046, pad=0.02).set_label("EPIfrac percentile")
    for cat in CATEGORY_ORDER:
        mask = category.eq(cat).to_numpy()
        if mask.any():
            axes[2].scatter(x[mask], y[mask], c=CATEGORY_COLORS[cat], s=size, linewidths=0, rasterized=True)
    set_axes(axes[2]); axes[2].set_title("Percentile quadrant")
    handles = [Line2D([0], [0], marker="o", color="w", label=CATEGORY_LABELS[c], markerfacecolor=CATEGORY_COLORS[c], markersize=6) for c in CATEGORY_ORDER]
    axes[2].legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=True, fontsize=6)
    fig.suptitle(f"{sample}: {cm} | {epi}", y=1.04, fontsize=9, fontweight="bold")
    save_pdf_svg(fig, out_dir / f"{sample}__{cm}__{epi}__percentile_quadrant_0p5")


def main():
    sample = args().sample
    started = time.time()
    score_path = SCORE_DIR / f"{sample}_tangram_pseudobulk_epi_cm_spot_scores.csv"
    score = pd.read_csv(score_path)
    pairs = pd.read_csv(PAIR_MANIFEST).query("sample == @sample")
    if len(pairs) != 165:
        raise ValueError(f"{sample}: expected 165 pairs, found {len(pairs)}")
    feature_cols = [c for c in score if c.startswith("CMact__") or c.startswith("EPIfrac__")]
    pct = score[feature_cols].rank(method="average", pct=True)
    size = max(0.15, min(4.0, 25000.0 / len(score)))
    raw_dir = FIG_ROOT / "raw_cmact_epifrac_by_sample" / sample
    pct_dir = FIG_ROOT / "percentile_quadrant_by_sample" / sample
    raw_dir.mkdir(parents=True, exist_ok=True)
    pct_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for i, pair in enumerate(pairs.itertuples(index=False), 1):
        cm, epi = pair.CM, pair.epi_subtype
        raw_stem = raw_dir / f"{sample}__{cm}__{epi}__raw_cmact_epifrac"
        pct_stem = pct_dir / f"{sample}__{cm}__{epi}__percentile_quadrant_0p5"
        if not (raw_stem.with_suffix(".pdf").exists() and raw_stem.with_suffix(".svg").exists()):
            raw_plot(score, sample, cm, epi, raw_dir, size)
        if not (pct_stem.with_suffix(".pdf").exists() and pct_stem.with_suffix(".svg").exists()):
            pct_plot(score, pct, sample, cm, epi, pct_dir, size)
        rows.append({"sample": sample, "CM": cm, "epi_subtype": epi, "raw_pdf": str(raw_stem.with_suffix('.pdf')), "raw_svg": str(raw_stem.with_suffix('.svg')), "percentile_pdf": str(pct_stem.with_suffix('.pdf')), "percentile_svg": str(pct_stem.with_suffix('.svg')), "completed": True})
        if i % 10 == 0:
            print(f"{sample} {i}/165", flush=True)
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(TABLE_DIR / f"{sample}_plot_manifest.csv", index=False)
    summary = {"sample": sample, "n_pairs": len(rows), "n_spatial_cells": len(score), "point_size": size, "elapsed_seconds": round(time.time()-started, 3)}
    (TABLE_DIR / f"{sample}_plot_summary.json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
