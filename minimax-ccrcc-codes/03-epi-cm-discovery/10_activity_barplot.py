"""Block 03 final: activity_df_tumor_vs_normal_mean_sd_barplot per SKILL.md.

Reads joint_CM_activity_tumor_vs_normal_mean_sd_summary.csv and produces a
side-by-side barplot (mean ± SD) for each CM, tumor vs normal-like, with
canonical CM order preserved.

Per updated SKILL.md:
  - "draw grouped bars for normal-like and tumor mean activity per CM with SD
    error bars"
  - "mark/annotate MWU/BH q significance only from the summary table"

This script:
  1. Loads activity_df_sample_by_CM.csv.
  2. Computes per-CM Mann-Whitney U test (tumor vs normal-like) + BH q.
  3. Updates joint_CM_activity_tumor_vs_normal_mean_sd_summary.csv with
     `mwu_p_value` and `q_value_bh` columns.
  4. Draws the barplot, annotating each CM with ns / * / ** / *** from q_value_bh.
"""
from __future__ import annotations

import time
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.family"] = "DejaVu Sans"
mpl.rcParams["axes.linewidth"] = 0.8

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
TAB = ROOT / "epi-cm-core-workflow/tables/03-epi-cm-discovery"
FIG = ROOT / "epi-cm-core-workflow/figures/03-epi-cm-discovery"
ACT = TAB / "activity_df_sample_by_CM.csv"
SUM = TAB / "joint_CM_activity_tumor_vs_normal_mean_sd_summary.csv"


def save_pdf_svg(fig, stem: Path):
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", dpi=300)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", dpi=300)
    plt.close(fig)


def q_label(q):
    if not np.isfinite(q):
        return "ns"
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return "ns"


def compute_mwu_q(activity_path: Path, summary_path: Path):
    """Compute per-CM Mann-Whitney U (tumor vs normal-like) and BH q.

    Updates the summary CSV with `mwu_p_value` and `q_value_bh` columns.
    """
    activity = pd.read_csv(activity_path)
    summary = pd.read_csv(summary_path)
    cm_cols = [c for c in activity.columns
               if c not in {"status", "non_epi_cells", "sample_id"}]
    tumor_vals = activity.loc[activity["status"] == "tumor", cm_cols]
    normal_vals = activity.loc[activity["status"] == "normal-like", cm_cols]

    rows = []
    for cm in cm_cols:
        t = pd.to_numeric(tumor_vals[cm], errors="coerce").dropna().to_numpy()
        n = pd.to_numeric(normal_vals[cm], errors="coerce").dropna().to_numpy()
        if len(t) < 1 or len(n) < 1:
            rows.append({"CM": cm, "mwu_p_value": np.nan})
            continue
        try:
            _, p = mannwhitneyu(t, n, alternative="two-sided")
        except ValueError:
            p = np.nan
        rows.append({"CM": cm, "mwu_p_value": float(p) if np.isfinite(p) else np.nan})
    p_df = pd.DataFrame(rows)
    valid = p_df["mwu_p_value"].notna()
    if valid.any():
        _, qvals, _, _ = multipletests(
            p_df.loc[valid, "mwu_p_value"].values, method="fdr_bh"
        )
        p_df.loc[valid, "q_value_bh"] = qvals
    else:
        p_df["q_value_bh"] = np.nan

    # Merge with summary (join on CM only; each CM has one tumor row + one normal row).
    summary_aug = summary.drop(columns=["mwu_p_value", "q_value_bh"], errors="ignore").merge(
        p_df, on="CM", how="left"
    )
    summary_aug.to_csv(summary_path, index=False)
    print(f"[barplot] updated summary with mwu_p_value + q_value_bh columns", flush=True)
    return summary_aug


def main():
    t0 = time.time()
    summary = compute_mwu_q(ACT, SUM)

    cls = pd.read_csv(TAB / "joint_module_classification.csv")[["CM", "class", "global_order"]]
    cm_order = cls.sort_values(["class", "global_order"])["CM"].astype(str).tolist()
    print(f"[barplot] CMs ({len(cm_order)}): {cm_order}", flush=True)

    pivot_mean = summary.pivot(index="CM", columns="status", values="mean").reindex(cm_order)
    pivot_sd = summary.pivot(index="CM", columns="status", values="sd").reindex(cm_order)
    qvals = summary.drop_duplicates("CM").set_index("CM").reindex(cm_order)["q_value_bh"]

    n = len(cm_order)
    fig, ax = plt.subplots(figsize=(max(9, n * 0.6), 5.6))
    x = np.arange(n)
    width = 0.38
    bars_t = ax.bar(x - width / 2, pivot_mean["tumor"].values, width,
                    yerr=pivot_sd["tumor"].values, capsize=3,
                    color="#d62728", edgecolor="black", linewidth=0.6, label="tumor")
    bars_n = ax.bar(x + width / 2, pivot_mean["normal-like"].values, width,
                    yerr=pivot_sd["normal-like"].values, capsize=3,
                    color="#1f77b4", edgecolor="black", linewidth=0.6, label="normal-like")

    # Significance stars above each CM pair (annotate MWU/BH q from summary table).
    ymax = np.nanmax(pivot_mean.values + pivot_sd.fillna(0).values)
    pad = ymax * 0.03
    star_y = ymax + pad
    bar_top_y = ymax + pad * 0.5
    for xi, cm in zip(x, cm_order):
        q = qvals.loc[cm]
        label = q_label(q)
        # Bracket between the two bars
        ax.plot([xi - width / 2, xi - width / 2, xi + width / 2, xi + width / 2],
                [bar_top_y, star_y, star_y, bar_top_y],
                color="black", linewidth=0.7)
        ax.text(xi, star_y + pad * 0.5, label, ha="center", va="bottom", fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(cm_order, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("CM activity (W)")
    ax.set_title("CM activity: tumor vs normal-like (mean ± SD; MWU + BH q stars)")
    ax.legend(loc="upper right", frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#d9d9d9", linewidth=0.5, alpha=0.7)
    ax.set_ylim(top=ymax * 1.25)
    plt.tight_layout()
    save_pdf_svg(fig, FIG / "activity_df_tumor_vs_normal_mean_sd_barplot")
    print(f"[barplot] activity_df_tumor_vs_normal_mean_sd_barplot saved", flush=True)
    print(f"[barplot] done. total={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()