#!/usr/bin/env python3
"""Reload one canonical spot-score CSV, plot all 100 pairs, and compute stats."""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import random
import re
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from scipy.stats import fisher_exact, norm, spearmanr

SEED = 42
random.seed(SEED)
np.random.seed(SEED)
ROOT = Path("/mnt/disk18t/lr_xcy/riku/brca_cm-epi_codex")
WORKFLOW = ROOT / "epi-cm-core-workflow"
MAPPING_TABLE_DIR = WORKFLOW / "tables/04-spatial-validation-optional-sig_genes/02-tangram-mapping"
STAT_TABLE_DIR = WORKFLOW / "tables/04-spatial-validation-optional-sig_genes/03-spatial-statistics-and-plotting/per_sample"
FIGURE_ROOT = WORKFLOW / "figures/04-spatial-validation-optional-sig_genes/03-spatial-statistics-and-plotting"
VALID_SAMPLES = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]
CUTOFF = 0.5

plt.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "font.family": "DejaVu Sans",
})

CM_PERCENTILE_CMAP = LinearSegmentedColormap.from_list("cm_percentile_light_to_blue", ["#f2f2f2", "#b7d7ea", "#4c78a8", "#2b6cb0"])
EPI_PERCENTILE_CMAP = LinearSegmentedColormap.from_list("epi_percentile_light_to_orange", ["#f2f2f2", "#fde6b3", "#fbbf24", "#f59e0b"])
CATEGORY_ORDER = ["both_low", "epi_high_only", "cm_high_only", "both_high"]
CATEGORY_LABELS = {"both_low": "CM<0.5, Epi<0.5", "epi_high_only": "Epi>=0.5 only", "cm_high_only": "CM>=0.5 only", "both_high": "Both>=0.5"}
CATEGORY_COLORS = {"both_low": "#d9d9d9", "epi_high_only": "#f59e0b", "cm_high_only": "#2b6cb0", "both_high": "#b91c1c"}


def pkg(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True, choices=VALID_SAMPLES)
    return parser.parse_args()


def cm_key(value: str) -> tuple[int, str]:
    match = re.search(r"CM(\d+)$", str(value))
    return (int(match.group(1)) if match else 10**9, str(value))


def robust_limits(values: pd.Series, q=(0.01, 0.99)) -> tuple[float, float]:
    arr = pd.to_numeric(values, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna().to_numpy()
    if arr.size == 0:
        return 0.0, 1.0
    lo, hi = np.nanquantile(arr, q)
    if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
        lo, hi = float(np.nanmin(arr)), float(np.nanmax(arr))
    if lo == hi:
        hi = lo + 1e-9
    return float(lo), float(hi)


def percentile_rank(values: pd.Series) -> pd.Series:
    return pd.to_numeric(values, errors="coerce").rank(method="average", pct=True).fillna(0.0)


def classify_percentiles(cm_pct: pd.Series, epi_pct: pd.Series) -> pd.Series:
    cm_high = cm_pct.ge(CUTOFF)
    epi_high = epi_pct.ge(CUTOFF)
    category = pd.Series("both_low", index=cm_pct.index, dtype="object")
    category.loc[epi_high & ~cm_high] = "epi_high_only"
    category.loc[cm_high & ~epi_high] = "cm_high_only"
    category.loc[cm_high & epi_high] = "both_high"
    return category


def set_spatial_axes(ax) -> None:
    ax.invert_yaxis()
    ax.set_box_aspect(1)
    ax.set_xlabel("array_col")
    ax.set_ylabel("array_row")


def save_pdf_svg(fig, stem: Path) -> str:
    pdf, svg = stem.with_suffix(".pdf"), stem.with_suffix(".svg")
    if pdf.exists() and svg.exists():
        plt.close(fig)
        return "reused_existing_complete_pair"
    if pdf.exists() or svg.exists():
        plt.close(fig)
        raise FileExistsError(f"Partial figure pair exists; refusing overwrite: {stem}")
    fig.savefig(pdf, bbox_inches="tight", dpi=300)
    fig.savefig(svg, bbox_inches="tight", dpi=300)
    plt.close(fig)
    return "generated"


def plot_raw(score, sample: str, cm: str, epi: str, out_dir: Path, size: float) -> str:
    cm_raw = pd.to_numeric(score[f"CMact__{cm}"], errors="coerce")
    epi_raw = pd.to_numeric(score[f"EPIfrac__{epi}"], errors="coerce")
    x, y = score["array_col"].to_numpy(float), score["array_row"].to_numpy(float)
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 3.2), constrained_layout=True)
    sc0 = axes[0].scatter(x, y, c=cm_raw, s=size, cmap="magma", norm=Normalize(*robust_limits(cm_raw)), linewidths=0, rasterized=False)
    set_spatial_axes(axes[0]); axes[0].set_title(f"CMact: {cm}"); fig.colorbar(sc0, ax=axes[0], fraction=0.046, pad=0.02).set_label("CMact")
    sc1 = axes[1].scatter(x, y, c=epi_raw, s=size, cmap="viridis", norm=Normalize(*robust_limits(epi_raw)), linewidths=0, rasterized=False)
    set_spatial_axes(axes[1]); axes[1].set_title(f"EPIfrac: {epi}"); fig.colorbar(sc1, ax=axes[1], fraction=0.046, pad=0.02).set_label("EPIfrac")
    fig.suptitle(f"{sample}: {cm} | {epi}", y=1.04, fontsize=9, fontweight="bold")
    return save_pdf_svg(fig, out_dir / f"{sample}__{cm}__{epi}__raw_cmact_epifrac")


def plot_percentile(score, sample: str, cm: str, epi: str, out_dir: Path, size: float) -> tuple[str, pd.Series, pd.Series, pd.Series]:
    cm_pct = percentile_rank(score[f"CMact__{cm}"])
    epi_pct = percentile_rank(score[f"EPIfrac__{epi}"])
    category = classify_percentiles(cm_pct, epi_pct)
    x, y = score["array_col"].to_numpy(float), score["array_row"].to_numpy(float)
    fig, axes = plt.subplots(1, 3, figsize=(14.8, 4.8), constrained_layout=True)
    sc0 = axes[0].scatter(x, y, c=cm_pct, cmap=CM_PERCENTILE_CMAP, vmin=0, vmax=1, s=size, marker="o", linewidths=0, rasterized=True)
    set_spatial_axes(axes[0]); axes[0].set_title(f"CM percentile: {cm}"); fig.colorbar(sc0, ax=axes[0], fraction=0.046, pad=0.02).set_label("CMact percentile")
    sc1 = axes[1].scatter(x, y, c=epi_pct, cmap=EPI_PERCENTILE_CMAP, vmin=0, vmax=1, s=size, marker="o", linewidths=0, rasterized=True)
    set_spatial_axes(axes[1]); axes[1].set_title(f"Epi percentile: {epi}"); fig.colorbar(sc1, ax=axes[1], fraction=0.046, pad=0.02).set_label("EPIfrac percentile")
    for cat in CATEGORY_ORDER:
        mask = category.eq(cat).to_numpy()
        if mask.any():
            axes[2].scatter(x[mask], y[mask], c=CATEGORY_COLORS[cat], s=size, marker="o", linewidths=0, rasterized=True)
    set_spatial_axes(axes[2]); axes[2].set_title("Percentile quadrant")
    handles = [Line2D([0], [0], marker="o", color="w", label=CATEGORY_LABELS[cat], markerfacecolor=CATEGORY_COLORS[cat], markersize=6) for cat in CATEGORY_ORDER]
    axes[2].legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=True, fontsize=6, borderpad=0.35, borderaxespad=0.0)
    fig.suptitle(f"{sample}: {cm} | {epi}", y=1.04, fontsize=9, fontweight="bold")
    status = save_pdf_svg(fig, out_dir / f"{sample}__{cm}__{epi}__percentile_quadrant_0p5")
    return status, cm_pct, epi_pct, category


def pair_statistics(score: pd.DataFrame, sample: str, cm: str, epi: str, cm_pct: pd.Series, epi_pct: pd.Series, category: pd.Series) -> dict[str, object]:
    x = pd.to_numeric(score[f"CMact__{cm}"], errors="coerce")
    y = pd.to_numeric(score[f"EPIfrac__{epi}"], errors="coerce")
    valid = x.notna() & y.notna() & np.isfinite(x) & np.isfinite(y)
    if valid.sum() >= 3 and x.loc[valid].nunique() > 1 and y.loc[valid].nunique() > 1:
        rho, spearman_p = spearmanr(x.loc[valid], y.loc[valid])
    else:
        rho, spearman_p = np.nan, np.nan
    counts = category.value_counts().reindex(CATEGORY_ORDER, fill_value=0)
    table = [[int(counts["both_high"]), int(counts["cm_high_only"])], [int(counts["epi_high_only"]), int(counts["both_low"])]]
    odds_ratio, fisher_p = fisher_exact(table, alternative="two-sided")
    cross = table[0][0] * table[1][1] - table[0][1] * table[1][0]
    direction = 1.0 if cross > 0 else (-1.0 if cross < 0 else 0.0)
    p_for_z = max(float(fisher_p), np.finfo(float).tiny)
    signed_z = direction * float(norm.isf(p_for_z / 2.0)) if direction else 0.0
    return {
        "sample": sample, "CM": cm, "epi_subtype": epi, "n_spots": int(valid.sum()),
        "spearman_rho": float(rho) if np.isfinite(rho) else np.nan, "spearman_p_value": float(spearman_p) if np.isfinite(spearman_p) else np.nan,
        "percentile_cutoff": CUTOFF, "both_low": int(counts["both_low"]), "epi_high_only": int(counts["epi_high_only"]), "cm_high_only": int(counts["cm_high_only"]), "both_high": int(counts["both_high"]),
        "fisher_odds_ratio": float(odds_ratio), "fisher_p_value": float(fisher_p), "fisher_direction": int(direction), "fisher_signed_z": signed_z,
        "raw_cm_min": float(x.min()), "raw_cm_max": float(x.max()), "raw_epi_min": float(y.min()), "raw_epi_max": float(y.max()),
        "cm_percentile_min": float(cm_pct.min()), "cm_percentile_max": float(cm_pct.max()), "epi_percentile_min": float(epi_pct.min()), "epi_percentile_max": float(epi_pct.max()),
    }


def main() -> None:
    started = time.time()
    sample = parse_args().sample
    score_path = MAPPING_TABLE_DIR / f"{sample}_tangram_pseudobulk_epi_cm_spot_scores.csv"
    completion_path = STAT_TABLE_DIR / f"{sample}_plot_stats_completion.json"
    stats_path = STAT_TABLE_DIR / f"{sample}_pair_statistics.csv"
    manifest_path = STAT_TABLE_DIR / f"{sample}_pair_manifest.csv"
    for path in [completion_path, stats_path, manifest_path]:
        if path.exists():
            raise FileExistsError(f"Refusing to overwrite existing table: {path}")
    score = pd.read_csv(score_path)
    if score["sample"].astype(str).nunique() != 1 or score["sample"].astype(str).iloc[0] != sample:
        raise ValueError(f"{sample}: score sample column mismatch")
    cm_names = sorted({column.split("__", 1)[1] for column in score.columns if column.startswith("CMact__")}, key=cm_key)
    epi_names = sorted({column.split("__", 1)[1] for column in score.columns if column.startswith("EPIfrac__")})
    if len(cm_names) != 10 or len(epi_names) != 10:
        raise ValueError(f"{sample}: expected 10 CM x 10 Epi, found {len(cm_names)} x {len(epi_names)}")
    STAT_TABLE_DIR.mkdir(parents=True, exist_ok=True)
    raw_dir = FIGURE_ROOT / "raw_cmact_epifrac_by_sample" / sample
    pct_dir = FIGURE_ROOT / "percentile_quadrant_by_sample" / sample
    raw_dir.mkdir(parents=True, exist_ok=True)
    pct_dir.mkdir(parents=True, exist_ok=True)
    size = float(max(4.0, min(14.0, 12000.0 / len(score))))
    stat_rows, manifest_rows = [], []
    generated, reused = 0, 0
    for cm in cm_names:
        for epi in epi_names:
            raw_status = plot_raw(score, sample, cm, epi, raw_dir, size)
            pct_status, cm_pct, epi_pct, category = plot_percentile(score, sample, cm, epi, pct_dir, size)
            generated += int(raw_status == "generated") + int(pct_status == "generated")
            reused += int(raw_status != "generated") + int(pct_status != "generated")
            stat_rows.append(pair_statistics(score, sample, cm, epi, cm_pct, epi_pct, category))
            manifest_rows.append({"sample": sample, "CM": cm, "epi_subtype": epi, "raw_pdf": str((raw_dir / f"{sample}__{cm}__{epi}__raw_cmact_epifrac.pdf").resolve()), "raw_svg": str((raw_dir / f"{sample}__{cm}__{epi}__raw_cmact_epifrac.svg").resolve()), "percentile_pdf": str((pct_dir / f"{sample}__{cm}__{epi}__percentile_quadrant_0p5.pdf").resolve()), "percentile_svg": str((pct_dir / f"{sample}__{cm}__{epi}__percentile_quadrant_0p5.svg").resolve())})
    pd.DataFrame(stat_rows).to_csv(stats_path, index=False)
    pd.DataFrame(manifest_rows).to_csv(manifest_path, index=False)
    pd.DataFrame([{"sample": sample, "code_file": str(Path(__file__).resolve()), "score_csv": str(score_path.resolve()), "n_spots": len(score), "cm_count": len(cm_names), "epi_count": len(epi_names), "pair_count": len(stat_rows), "raw_spearman_values": "original CMact/EPIfrac", "percentile_method": "pandas rank(method=average,pct=True)", "percentile_cutoff": CUTOFF, "fisher_table": "[[both_high,cm_high_only],[epi_high_only,both_low]]", "point_size": size, "raw_plot_rasterized": False, "percentile_plot_rasterized": True, "formats": "PDF,SVG", "seed": SEED}]).to_csv(STAT_TABLE_DIR / f"{sample}_plot_stats_parameters.csv", index=False)
    (STAT_TABLE_DIR / f"{sample}_package_versions.txt").write_text(f"python={sys.version.split()[0]}\nnumpy={pkg('numpy')}\npandas={pkg('pandas')}\nscipy={pkg('scipy')}\nmatplotlib={pkg('matplotlib')}\ncode={Path(__file__).resolve()}\nseed={SEED}\n", encoding="utf-8")
    completion = {"status": "completed", "sample": sample, "n_spots": len(score), "n_pairs": len(stat_rows), "expected_pairs": 100, "n_pair_figure_families_generated": generated, "n_pair_figure_families_reused": reused, "n_pdf": len(list(raw_dir.glob("*.pdf"))) + len(list(pct_dir.glob("*.pdf"))), "n_svg": len(list(raw_dir.glob("*.svg"))) + len(list(pct_dir.glob("*.svg"))), "elapsed_seconds": time.time() - started, "seed": SEED}
    completion_path.write_text(json.dumps(completion, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(completion, indent=2), flush=True)


if __name__ == "__main__":
    main()
