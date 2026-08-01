"""Block 04 step 4: per-sample worker for figures + statistics."""
from __future__ import annotations

import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from scipy.stats import fisher_exact, spearmanr

SEED = 42
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
TAB = ROOT / "epi-cm-core-workflow/tables/04-spatial-validation-optional/01-tangram-spatial-mapping"
FIG = ROOT / "epi-cm-core-workflow/figures/04-spatial-validation-optional/01-tangram-spatial-mapping"
FIG_RAW = FIG / "raw_cmact_epifrac_by_sample"
FIG_PCT = FIG / "percentile_quadrant_by_sample"
TMP_STATS = TAB / "tmp_stats_per_sample"
for d in [FIG_RAW, FIG_PCT, TMP_STATS]:
    d.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "font.family": "DejaVu Sans",
})

CATEGORY_COLORS = {
    "both_low": "#d9d9d9", "epi_high_only": "#f59e0b",
    "cm_high_only": "#2b6cb0", "both_high": "#b91c1c",
}
CM_PERCENTILE_CMAP = LinearSegmentedColormap.from_list(
    "cm_pct", ["#f2f2f2", "#b7d7ea", "#4c78a8", "#2b6cb0"])
EPI_PERCENTILE_CMAP = LinearSegmentedColormap.from_list(
    "epi_pct", ["#f2f2f2", "#fde6b3", "#fbbf24", "#f59e0b"])
CATEGORY_ORDER = ["both_low", "epi_high_only", "cm_high_only", "both_high"]
CATEGORY_LABELS = {
    "both_low": "CM<0.5, Epi<0.5",
    "epi_high_only": "Epi>=0.5 only",
    "cm_high_only": "CM>=0.5 only",
    "both_high": "Both>=0.5",
}


def robust_limits(values, q=(0.01, 0.99)):
    arr = (pd.to_numeric(values, errors="coerce")
           .replace([np.inf, -np.inf], np.nan).dropna().to_numpy())
    if arr.size == 0: return 0.0, 1.0
    lo, hi = np.nanquantile(arr, q)
    if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
        lo, hi = float(np.nanmin(arr)), float(np.nanmax(arr))
    if lo == hi: hi = lo + 1e-9
    return float(lo), float(hi)


def percentile_rank(values):
    s = pd.to_numeric(values, errors="coerce")
    return s.rank(method="average", pct=True).fillna(0.0)


def classify_percentiles(cm_pct, epi_pct, cutoff=0.5):
    cm_high = pd.to_numeric(cm_pct, errors="coerce") >= cutoff
    epi_high = pd.to_numeric(epi_pct, errors="coerce") >= cutoff
    category = pd.Series("both_low", index=cm_pct.index, dtype="object")
    category.loc[epi_high & ~cm_high] = "epi_high_only"
    category.loc[cm_high & ~epi_high] = "cm_high_only"
    category.loc[cm_high & epi_high] = "both_high"
    return category


def set_spatial_axes(ax, x, y):
    ax.invert_yaxis(); ax.set_box_aspect(1)
    ax.set_xlabel("array_col"); ax.set_ylabel("array_row")


def save_pdf_svg(fig, stem):
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", dpi=300)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", dpi=300)
    plt.close(fig)


def plot_raw_cmact_epifrac(score, sample, cm, epi, out_dir, size):
    cm_raw = pd.to_numeric(score[f"CMact__{cm}"], errors="coerce")
    epi_raw = pd.to_numeric(score[f"EPIfrac__{epi}"], errors="coerce")
    x = score["array_col"].to_numpy(dtype=float)
    y = score["array_row"].to_numpy(dtype=float)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6), constrained_layout=True)
    sc_cm = axes[0].scatter(x, y, c=cm_raw, s=size, cmap="magma",
        norm=Normalize(*robust_limits(cm_raw)),
        marker="o", linewidths=0, rasterized=True)
    set_spatial_axes(axes[0], x, y); axes[0].set_title(f"CMact: {cm}")
    fig.colorbar(sc_cm, ax=axes[0], fraction=0.046, pad=0.02).set_label("CMact")
    sc_epi = axes[1].scatter(x, y, c=epi_raw, s=size, cmap="viridis",
        norm=Normalize(*robust_limits(epi_raw)),
        marker="o", linewidths=0, rasterized=True)
    set_spatial_axes(axes[1], x, y); axes[1].set_title(f"EPIfrac: {epi}")
    fig.colorbar(sc_epi, ax=axes[1], fraction=0.046, pad=0.02).set_label("EPIfrac")
    fig.suptitle(f"{sample}: {cm} | {epi}", y=1.04, fontsize=9, fontweight="bold")
    out_dir.mkdir(parents=True, exist_ok=True)
    save_pdf_svg(fig, out_dir / f"{sample}__{cm}__{epi}__raw_cmact_epifrac")


def plot_percentile_quadrant(score, sample, cm, epi, out_dir, size):
    cm_pct = percentile_rank(score[f"CMact__{cm}"])
    epi_pct = percentile_rank(score[f"EPIfrac__{epi}"])
    category = classify_percentiles(cm_pct, epi_pct, cutoff=0.5)
    x = score["array_col"].to_numpy(dtype=float)
    y = score["array_row"].to_numpy(dtype=float)
    fig, axes = plt.subplots(1, 3, figsize=(14.8, 4.8), constrained_layout=True)
    sc_cm = axes[0].scatter(x, y, c=cm_pct, cmap=CM_PERCENTILE_CMAP, vmin=0, vmax=1,
        s=size, marker="o", linewidths=0, rasterized=True)
    set_spatial_axes(axes[0], x, y); axes[0].set_title(f"CM percentile: {cm}")
    fig.colorbar(sc_cm, ax=axes[0], fraction=0.046, pad=0.02).set_label("CMact percentile")
    sc_epi = axes[1].scatter(x, y, c=epi_pct, cmap=EPI_PERCENTILE_CMAP, vmin=0, vmax=1,
        s=size, marker="o", linewidths=0, rasterized=True)
    set_spatial_axes(axes[1], x, y); axes[1].set_title(f"Epi percentile: {epi}")
    fig.colorbar(sc_epi, ax=axes[1], fraction=0.046, pad=0.02).set_label("EPIfrac percentile")
    for cat in CATEGORY_ORDER:
        mask = category.eq(cat).to_numpy()
        if mask.any():
            axes[2].scatter(x[mask], y[mask], c=CATEGORY_COLORS[cat], s=size,
                marker="o", linewidths=0, rasterized=True, label=CATEGORY_LABELS[cat])
    set_spatial_axes(axes[2], x, y); axes[2].set_title("Percentile quadrant")
    handles = [Line2D([0], [0], marker="o", color="w", label=CATEGORY_LABELS[cat],
        markerfacecolor=CATEGORY_COLORS[cat], markersize=6) for cat in CATEGORY_ORDER]
    axes[2].legend(handles=handles, loc="upper left",
        bbox_to_anchor=(1.02, 1.0), frameon=True, fontsize=6,
        borderpad=0.35, borderaxespad=0.0)
    fig.suptitle(f"{sample}: {cm} | {epi}", y=1.04, fontsize=9, fontweight="bold")
    out_dir.mkdir(parents=True, exist_ok=True)
    save_pdf_svg(fig, out_dir / f"{sample}__{cm}__{epi}__percentile_quadrant_0p5")


def compute_pair_statistics(score, sample, cm, epi, status):
    cm_raw = pd.to_numeric(score[f"CMact__{cm}"], errors="coerce").to_numpy()
    epi_raw = pd.to_numeric(score[f"EPIfrac__{epi}"], errors="coerce").to_numpy()
    ok = ~np.isnan(cm_raw) & ~np.isnan(epi_raw)
    n_spots = int(ok.sum())
    spearman_r = spearman_p = fisher_p = fisher_p_signed = np.nan
    a = b = c_ = d = 0
    if n_spots >= 3 and len(set(cm_raw[ok])) > 1 and len(set(epi_raw[ok])) > 1:
        spearman_r, spearman_p = spearmanr(cm_raw[ok], epi_raw[ok])
        cm_pct = percentile_rank(pd.Series(cm_raw)).to_numpy()
        epi_pct = percentile_rank(pd.Series(epi_raw)).to_numpy()
        cmh = (cm_pct >= 0.5) & ok; eph = (epi_pct >= 0.5) & ok
        a = int((cmh & eph).sum()); b = int((cmh & ~eph).sum())
        c_ = int((~cmh & eph).sum()); d = int((~cmh & ~eph).sum())
        if a + b + c_ + d > 0:
            _, fisher_p = fisher_exact([[a, b], [c_, d]])
        if np.isfinite(fisher_p):
            fisher_p_signed = float(np.sign((a * d) - (b * c_)) * fisher_p)
    return {
        "sample": sample, "CM": cm, "epi_subtype": epi, "status": status,
        "n_spots": n_spots,
        "spearman_r": float(spearman_r) if np.isfinite(spearman_r) else np.nan,
        "spearman_p": float(spearman_p) if np.isfinite(spearman_p) else np.nan,
        "n_both_high": a, "n_cm_high_only": b,
        "n_epi_high_only": c_, "n_both_low": d,
        "fisher_p": float(fisher_p) if np.isfinite(fisher_p) else np.nan,
        "fisher_p_signed": fisher_p_signed,
    }


def discover_cm_epi_pairs(score):
    cm_names = [c.split("__", 1)[1] for c in score.columns if c.startswith("CMact__")]
    epi_names = [c.split("__", 1)[1] for c in score.columns if c.startswith("EPIfrac__")]
    cm_names = sorted(set(cm_names)); epi_names = sorted(set(epi_names))
    pairs = [{"CM": cm, "epi_subtype": epi} for cm in cm_names for epi in epi_names]
    if not pairs: raise ValueError("No CMact__/EPIfrac__ pairs found")
    return cm_names, epi_names, pairs


def infer_sample_from_path(score_path):
    return score_path.stem.replace("_tangram_pseudobulk_epi_cm_spot_scores", "")


def process_sample(sample, status, size=4):
    score_path = TAB / f"{sample}_tangram_pseudobulk_epi_cm_spot_scores.csv"
    if not score_path.exists():
        raise FileNotFoundError(f"{sample}: spot-score CSV not found at {score_path}")
    score = pd.read_csv(score_path)
    if "sample" not in score.columns: score["sample"] = sample
    if "array_col" not in score.columns and "spatial_x" in score.columns:
        score["array_col"] = score["spatial_x"]
    if "array_row" not in score.columns and "spatial_y" in score.columns:
        score["array_row"] = score["spatial_y"]
    cm_names, epi_names, pairs = discover_cm_epi_pairs(score)
    print(f"[{sample}] {len(cm_names)} CMs x {len(epi_names)} epi = {len(pairs)} pairs "
          f"(status={status})", flush=True)
    raw_sample_dir = FIG_RAW / sample; pct_sample_dir = FIG_PCT / sample
    raw_sample_dir.mkdir(parents=True, exist_ok=True)
    pct_sample_dir.mkdir(parents=True, exist_ok=True)
    spearman_rows, fisher_rows = [], []
    t1 = time.time()
    for pair in pairs:
        cm = pair["CM"]; epi = pair["epi_subtype"]
        plot_raw_cmact_epifrac(score, sample, cm, epi, raw_sample_dir, size=size)
        plot_percentile_quadrant(score, sample, cm, epi, pct_sample_dir, size=size)
        stats = compute_pair_statistics(score, sample, cm, epi, status)
        spearman_rows.append({
            "sample": stats["sample"], "CM": stats["CM"], "epi_subtype": stats["epi_subtype"],
            "status": stats["status"], "n_spots": stats["n_spots"],
            "spearman_r": stats["spearman_r"], "spearman_p": stats["spearman_p"],
        })
        fisher_rows.append({
            "sample": stats["sample"], "CM": stats["CM"], "epi_subtype": stats["epi_subtype"],
            "status": stats["status"], "n_spots": stats["n_spots"],
            "n_both_high": stats["n_both_high"], "n_cm_high_only": stats["n_cm_high_only"],
            "n_epi_high_only": stats["n_epi_high_only"], "n_both_low": stats["n_both_low"],
            "fisher_p": stats["fisher_p"], "fisher_p_signed": stats["fisher_p_signed"],
        })
    elapsed = time.time() - t1
    print(f"  {sample}: {len(pairs)*2} figures + {len(pairs)} stats in {elapsed:.1f}s", flush=True)
    pd.DataFrame(spearman_rows).to_csv(TMP_STATS / f"{sample}_spearman.csv", index=False)
    pd.DataFrame(fisher_rows).to_csv(TMP_STATS / f"{sample}_fisher.csv", index=False)


def main():
    t0 = time.time()
    score_paths = sorted(TAB.glob("*_tangram_pseudobulk_epi_cm_spot_scores.csv"))
    if not score_paths:
        raise FileNotFoundError(f"No per-sample Tangram spot-score CSVs found in {TAB}")
    sample_scope = pd.read_csv(TAB / "sample_scope.csv")
    scope_lookup = dict(zip(sample_scope["sample"].astype(str), sample_scope["status"].astype(str)))
    if len(sys.argv) >= 2 and sys.argv[1] != "all":
        target = sys.argv[1]; status = scope_lookup.get(target, "unknown")
        process_sample(target, status)
    else:
        for sp in score_paths:
            sample = infer_sample_from_path(sp)
            status = scope_lookup.get(sample, "unknown")
            process_sample(sample, status)
    print(f"all done total={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()