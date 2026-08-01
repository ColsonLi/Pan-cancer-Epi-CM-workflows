"""Block 04 step 5: Stouffer aggregation + final heatmaps."""
from __future__ import annotations

import re
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

SEED = 42
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
TAB = ROOT / "epi-cm-core-workflow/tables/04-spatial-validation-optional/01-tangram-spatial-mapping"
FIG = ROOT / "epi-cm-core-workflow/figures/04-spatial-validation-optional/01-tangram-spatial-mapping"
FIG_STOUFFER = FIG / "stouffer_heatmaps"
TAB_STATS_ALL = ROOT / "epi-cm-core-workflow/tables/04-spatial-validation-optional/statistics/all-samples"
TAB_STATS_TUMOR = ROOT / "epi-cm-core-workflow/tables/04-spatial-validation-optional/statistics/tumor-only"
for d in [FIG_STOUFFER, TAB_STATS_ALL, TAB_STATS_TUMOR]:
    d.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "text.usetex": False, "font.family": "DejaVu Sans",
    "font.size": 8, "axes.titlesize": 9,
})


def q_label(q):
    if not np.isfinite(q): return "ns"
    if q < 0.001: return "***"
    if q < 0.01: return "**"
    if q < 0.05: return "*"
    return "ns"


def cm_sort_key(value):
    m = re.search(r"joint_(\d+)", str(value))
    return (int(m.group(1)) if m else 999, str(value))


def ordered_stouffer_matrix(df):
    z_value = df.pivot(index="epi_subtype", columns="CM", values="combined_signed_z")
    q_value = df.pivot(index="epi_subtype", columns="CM", values="q_value_bh")
    rows = sorted(z_value.index); cols = sorted(z_value.columns, key=cm_sort_key)
    return z_value.reindex(index=rows, columns=cols), q_value.reindex(index=rows, columns=cols)


def finite_symmetric_limit(values):
    arr = values.to_numpy(dtype=float)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0: return 1.0
    vmax = float(np.nanmax(np.abs(finite)))
    return vmax if vmax > 0 else 1.0


def q_annot(value_matrix, q_matrix):
    annot = pd.DataFrame("", index=value_matrix.index, columns=value_matrix.columns)
    for row in annot.index:
        for col in annot.columns:
            annot.loc[row, col] = q_label(q_matrix.loc[row, col])
    return annot


def plot_sample_stouffer_signed_z(input_csv, out_dir, output_stem, title):
    df = pd.read_csv(input_csv)
    z_matrix, q_matrix = ordered_stouffer_matrix(df)
    vmax = finite_symmetric_limit(z_matrix)
    fig_w = max(7.8, 0.52 * z_matrix.shape[1] + 2.6)
    fig_h = max(4.8, 0.42 * z_matrix.shape[0] + 1.9)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(z_matrix, cmap="coolwarm", center=0, vmin=-vmax, vmax=vmax,
        linewidths=0.25, linecolor="white",
        annot=q_annot(z_matrix, q_matrix), fmt="",
        annot_kws={"fontsize": 6.5, "color": "black", "linespacing": 0.9},
        cbar_kws={"label": "combined signed Z"}, square=True, ax=ax)
    ax.set_title(title); ax.set_xlabel("CM"); ax.set_ylabel("Epithelial subtype")
    ax.tick_params(axis="x", labelrotation=90); ax.tick_params(axis="y", labelrotation=0)
    fig.tight_layout()
    fig.savefig(out_dir / f"{output_stem}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(out_dir / f"{output_stem}.svg", bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  heatmap: {out_dir}/{output_stem}.pdf + .svg", flush=True)


def fisher_p_to_z(p_signed):
    if not np.isfinite(p_signed): return np.nan
    s = np.sign(p_signed)
    if s == 0: return 0.0
    ap = abs(p_signed)
    if ap >= 1.0: return 0.0
    if ap <= 0.0: return float(s * 40.0)
    return float(s * norm.isf(ap / 2.0))


def _stouffer_one_pair(group):
    s = group["fisher_p_signed"].to_numpy(dtype=float)
    z = np.array([fisher_p_to_z(p) for p in s], dtype=float)
    n_valid = int(np.sum(np.isfinite(z)))
    if n_valid == 0:
        return pd.Series({"n_samples": 0, "combined_signed_z": np.nan, "combined_p": np.nan})
    Z = float(np.nansum(z) / np.sqrt(n_valid))
    p_two = float(2.0 * (1.0 - norm.cdf(abs(Z)))) if np.isfinite(Z) else np.nan
    return pd.Series({"n_samples": n_valid, "combined_signed_z": Z, "combined_p": p_two})


def stouffer_signed_from_fisher(per_pair_signed):
    g = per_pair_signed.groupby(["CM", "epi_subtype"], as_index=False).apply(
        _stouffer_one_pair, include_groups=False).reset_index(drop=True)
    valid = g["combined_p"].notna()
    if valid.any():
        _, qvals, _, _ = multipletests(g.loc[valid, "combined_p"].values, method="fdr_bh")
        g.loc[valid, "q_value_bh"] = qvals
    else:
        g["q_value_bh"] = np.nan
    g = g.sort_values("q_value_bh", key=lambda x: x.fillna(1.0))
    return g


def merge_per_sample_stats():
    fisher_paths = sorted((TAB / "tmp_stats_per_sample").glob("*_fisher.csv"))
    spearman_paths = sorted((TAB / "tmp_stats_per_sample").glob("*_spearman.csv"))
    if not fisher_paths or not spearman_paths:
        raise FileNotFoundError(f"No per-sample temp stats in {TAB / 'tmp_stats_per_sample'}")
    fisher_df = pd.concat([pd.read_csv(p) for p in fisher_paths], ignore_index=True)
    spearman_df = pd.concat([pd.read_csv(p) for p in spearman_paths], ignore_index=True)
    if "fisher_p_signed" not in fisher_df.columns:
        a = fisher_df["n_both_high"].to_numpy()
        b = fisher_df["n_cm_high_only"].to_numpy()
        c_ = fisher_df["n_epi_high_only"].to_numpy()
        d = fisher_df["n_both_low"].to_numpy()
        p = fisher_df["fisher_p"].to_numpy()
        sign = np.sign((a * d) - (b * c_))
        fisher_df["fisher_p_signed"] = np.where(np.isnan(p), np.nan, sign * p)
    spearman_df.to_csv(TAB / "per_sample_spearman.csv", index=False)
    fisher_df.to_csv(TAB / "percentile_quadrant_fisher_per_sample.csv", index=False)
    print(f"per_sample_spearman: {len(spearman_df)} rows; "
          f"percentile_quadrant_fisher_per_sample: {len(fisher_df)} rows", flush=True)
    manifest_df = spearman_df[["sample", "CM", "epi_subtype", "status"]].copy()
    manifest_df.to_csv(TAB / "all_sample_cm_epi_pair_manifest.csv", index=False)
    print(f"all_sample_cm_epi_pair_manifest: {len(manifest_df)} rows", flush=True)
    return spearman_df, fisher_df


def main():
    t0 = time.time()
    spearman_df, fisher_df = merge_per_sample_stats()
    g_all = stouffer_signed_from_fisher(fisher_df)
    n_sig_all = int((g_all["q_value_bh"] < 0.05).sum()) if "q_value_bh" in g_all.columns else 0
    print(f"all_samples: stouffer {len(g_all)} pairs, {n_sig_all} q<0.05", flush=True)
    g_all.to_csv(TAB / "percentile_quadrant_fisher_sample_stouffer_all_samples.csv", index=False)
    g_all.to_csv(TAB_STATS_ALL / "percentile_quadrant_fisher_sample_stouffer_all_samples.csv", index=False)
    sub_tumor = fisher_df[fisher_df["status"] == "tumor"]
    g_tumor = None
    if len(sub_tumor) > 0:
        g_tumor = stouffer_signed_from_fisher(sub_tumor)
        n_sig_t = int((g_tumor["q_value_bh"] < 0.05).sum()) if "q_value_bh" in g_tumor.columns else 0
        print(f"tumor_only (12 samples, no R_cor/R_med): stouffer {len(g_tumor)} pairs, {n_sig_t} q<0.05", flush=True)
        g_tumor.to_csv(TAB / "percentile_quadrant_fisher_sample_stouffer_tumor_only.csv", index=False)
        g_tumor.to_csv(TAB_STATS_TUMOR / "percentile_quadrant_fisher_sample_stouffer_tumor_only.csv", index=False)
    plot_sample_stouffer_signed_z(
        TAB / "percentile_quadrant_fisher_sample_stouffer_all_samples.csv",
        FIG_STOUFFER, "all_samples_stouffer_signedZ_qstars",
        "All samples percentile quadrant Fisher Stouffer (signed Z, BH q)",
    )
    if g_tumor is not None:
        plot_sample_stouffer_signed_z(
            TAB / "percentile_quadrant_fisher_sample_stouffer_tumor_only.csv",
            FIG_STOUFFER, "tumor_only_stouffer_signedZ_qstars",
            "Tumor-only percentile quadrant Fisher Stouffer (signed Z, BH q)",
        )
    print(f"merge done total={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()