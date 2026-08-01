"""Module 03 step 13: generate display-scale variant tables + with_series_clustermap tables + heatmaps.

Per updated SKILL.md Module 03 (canonical plotting rules):
  Display scale definitions:
    raw               = original values, no display normalization
    standard_scale_col= per-column min-max to [0, 1]
    zscore            = per-column center by mean, divide by std (ddof=0)
    robust            = per-column center by median, divide by IQR

  CM activity heatmaps (sample x CM, status hidden labels):
    w_df_activity_sample_by_CM_raw.{pdf,svg}
    w_df_activity_sample_by_CM_zscore.{pdf,svg}
    w_df_activity_sample_by_CM_robust.{pdf,svg}
    w_df_activity_sample_by_CM_standard_scale_col.{pdf,svg}
    w_df_activity_sample_activity_per_CM_standard_scale_col_clustermap.{pdf,svg}   (no row annotations)
    w_df_activity_sample_activity_per_CM_standard_scale_col_with_series_clustermap.{pdf,svg}  (Series + Status tracks)
  Plus matrix + row-annotation tables for the with_series clustermap.

  CM loading heatmaps (cell_subtype x CM, after orientation validation):
    h_df_loading_cell_subtype_by_CM_raw.{pdf,svg}
    h_df_loading_cell_subtype_by_CM_zscore.{pdf,svg}
    h_df_loading_cell_subtype_by_CM_robust.{pdf,svg}
    h_df_loading_cell_subtype_by_CM_standard_scale_col.{pdf,svg}

Display scales are display-only — do NOT use them for statistical tests.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
TAB = ROOT / "epi-cm-core-workflow/tables/03-epi-cm-discovery"
FIG = ROOT / "epi-cm-core-workflow/figures/03-epi-cm-discovery"
TAB.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "font.family": "DejaVu Sans",
})


def save_pdf_svg(fig, stem):
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", dpi=300)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", dpi=300)
    plt.close(fig)


def standardize_column_scale(M):
    """Each CM column independently min-max scaled to [0, 1]."""
    M = M.astype(float).copy()
    for col in M.columns:
        v = M[col].to_numpy()
        lo, hi = np.nanmin(v), np.nanmax(v)
        if hi > lo:
            M[col] = (v - lo) / (hi - lo)
        else:
            M[col] = 0.0
    return M


def zscore_column(M, ddof=0):
    """Per-column center by mean, divide by std (ddof=0)."""
    M = M.astype(float).copy()
    for col in M.columns:
        v = M[col].to_numpy()
        m = np.nanmean(v)
        s = np.nanstd(v, ddof=ddof)
        if s > 0:
            M[col] = (v - m) / s
        else:
            M[col] = 0.0
    return M


def robust_column(M):
    """Per-column center by median, divide by IQR."""
    M = M.astype(float).copy()
    for col in M.columns:
        v = M[col].to_numpy()
        med = np.nanmedian(v)
        q1 = np.nanquantile(v, 0.25); q3 = np.nanquantile(v, 0.75)
        iqr = q3 - q1
        if iqr > 0:
            M[col] = (v - med) / iqr
        else:
            M[col] = 0.0
    return M


def detect_orientation(df, expected_index_name=None):
    """Detect orientation: rows=sample vs rows=cell_subtype.

    Returns: ('sample_x_cm', df) | ('subtype_x_cm', df) | ('cm_x_subtype', transposed df)
    H_df.csv is conventionally CM x cell_subtype unless filename says otherwise.
    """
    n = len(df)
    n_cols = len(df.columns)
    if expected_index_name in df.columns:
        if expected_index_name in {"cell_subtype"}:
            return "subtype_x_cm", df.set_index(expected_index_name)
        if expected_index_name in {"sample_id", "sample"}:
            return "sample_x_cm", df.set_index(expected_index_name)
    if n_cols > n:
        return "cm_x_subtype", df.set_index(df.columns[0]).T
    return "subtype_x_cm", df.set_index(df.columns[0])


def plot_heatmap(M, stem, title, cmap, cbar_label, hide_row_labels=False):
    fig, ax = plt.subplots(figsize=(max(6.0, 0.32 * M.shape[1] + 2.0),
                                     max(4.0, 0.32 * M.shape[0] + 1.5)))
    sns.heatmap(M, cmap=cmap, ax=ax, cbar_kws={"label": cbar_label},
                linewidths=0, linecolor=None)
    ax.set_title(title)
    ax.set_xlabel("CM"); ax.set_ylabel(M.index.name or "")
    if hide_row_labels:
        ax.set_yticklabels([])
    fig.tight_layout()
    save_pdf_svg(fig, Path(stem))
    print(f"  {stem}.pdf + .svg saved")


def main():
    # ---- Activity tables (sample x CM) ----
    act = pd.read_csv(TAB / "activity_df_sample_by_CM.csv")
    series_col = None
    for cand in ["series", "Series", "dataset", "Dataset"]:
        if cand in act.columns:
            series_col = cand; break
    if series_col is None:
        series_col = "series"
        act[series_col] = "all"

    cm_cols = [c for c in act.columns if c not in {"sample_id", "status", "non_epi_cells", series_col}]
    act_num = act.set_index("sample_id")[cm_cols].astype(float)
    raw = act_num.copy()
    raw.to_csv(TAB / "w_df_activity_sample_by_CM_raw.csv")
    zscore_column(act_num).to_csv(TAB / "w_df_activity_sample_by_CM_zscore.csv")
    robust_column(act_num).to_csv(TAB / "w_df_activity_sample_by_CM_robust.csv")
    standardize_column_scale(act_num).to_csv(TAB / "w_df_activity_sample_by_CM_standard_scale_col.csv")

    # Status / series row annotations for the with_series clustermap
    row_meta = act.set_index("sample_id")[["status"]].copy()
    row_meta[series_col] = act[series_col].values
    row_meta.to_csv(TAB / "w_df_activity_sample_by_CM_standard_scale_col_with_series_row_annotations.csv")

    # Heatmaps
    plot_heatmap(raw, FIG / "w_df_activity_sample_by_CM_raw",
                 "CM activity (raw) per sample", "viridis", "W activity",
                 hide_row_labels=True)
    plot_heatmap(pd.read_csv(TAB / "w_df_activity_sample_by_CM_zscore.csv", index_col=0),
                 FIG / "w_df_activity_sample_by_CM_zscore",
                 "CM activity (z-score per column)", "RdBu_r", "z-score",
                 hide_row_labels=True)
    plot_heatmap(pd.read_csv(TAB / "w_df_activity_sample_by_CM_robust.csv", index_col=0),
                 FIG / "w_df_activity_sample_by_CM_robust",
                 "CM activity (robust per column)", "RdBu_r", "robust Z",
                 hide_row_labels=True)
    plot_heatmap(pd.read_csv(TAB / "w_df_activity_sample_by_CM_standard_scale_col.csv", index_col=0),
                 FIG / "w_df_activity_sample_by_CM_standard_scale_col",
                 "CM activity (column min-max)", "viridis", "min-max [0,1]",
                 hide_row_labels=True)

    # Clustermap variants (seaborn.clustermap)
    std = pd.read_csv(TAB / "w_df_activity_sample_by_CM_standard_scale_col.csv", index_col=0)
    # Unannotated clustermap
    g = sns.clustermap(std, cmap="viridis", figsize=(max(7.0, 0.32 * std.shape[1] + 2.5),
                                                       max(6.0, 0.32 * std.shape[0] + 2.0)),
                       cbar_kws={"label": "min-max [0,1]"},
                       xticklabels=True, yticklabels=False)
    g.ax_heatmap.set_xlabel("CM"); g.ax_heatmap.set_ylabel("sample")
    g.savefig(FIG / "w_df_activity_sample_activity_per_CM_standard_scale_col_clustermap.pdf",
              bbox_inches="tight", dpi=300)
    g.savefig(FIG / "w_df_activity_sample_activity_per_CM_standard_scale_col_clustermap.svg",
              bbox_inches="tight", dpi=300)
    plt.close("all")
    print(f"  w_df_activity_sample_activity_per_CM_standard_scale_col_clustermap saved")

    # Save matrix for with_series clustermap (used to drive clustermap below)
    std.to_csv(TAB / "w_df_activity_sample_by_CM_standard_scale_col_with_series_plot_matrix.csv")

    # Annotated clustermap with Series + Status tracks
    status_palette = {"tumor": "#d62728", "normal-like": "#1f77b4"}
    series_palette = {"all": "#7f7f7f"}
    if isinstance(act[series_col].astype(str).unique().tolist(), list):
        for i, s in enumerate(sorted(set(act[series_col].astype(str)))):
            series_palette[s] = plt.get_cmap("tab10")(i)
    row_color_series = row_meta[series_col].astype(str).map(series_palette)
    row_color_status = row_meta["status"].astype(str).map(status_palette)
    row_colors_df = pd.DataFrame({
        "status": row_color_status,
        series_col: row_color_series,
    }, index=std.index)
    g2 = sns.clustermap(std, cmap="viridis",
                        figsize=(max(8.0, 0.32 * std.shape[1] + 4.0),
                                 max(7.0, 0.32 * std.shape[0] + 3.0)),
                        row_colors=row_colors_df,
                        cbar_kws={"label": "min-max [0,1]"},
                        xticklabels=True, yticklabels=False)
    g2.ax_heatmap.set_xlabel("CM"); g2.ax_heatmap.set_ylabel("sample")
    g2.savefig(FIG / "w_df_activity_sample_activity_per_CM_standard_scale_col_with_series_clustermap.pdf",
               bbox_inches="tight", dpi=300)
    g2.savefig(FIG / "w_df_activity_sample_activity_per_CM_standard_scale_col_with_series_clustermap.svg",
               bbox_inches="tight", dpi=300)
    plt.close("all")
    print(f"  w_df_activity_sample_activity_per_CM_standard_scale_col_with_series_clustermap saved")

    # ---- Loading tables (cell_subtype x CM) ----
    load = pd.read_csv(TAB / "loading_df_cell_subtype_by_CM.csv")
    ori, load_M = detect_orientation(load, expected_index_name="cell_subtype")
    print(f"[load] detected orientation: {ori}; shape={load_M.shape}", flush=True)
    if ori == "cm_x_subtype":
        load_M = load_M.T
        print(f"[load] transposed to cell_subtype x CM: {load_M.shape}", flush=True)
    raw_L = load_M.astype(float).copy()
    raw_L.to_csv(TAB / "h_df_loading_cell_subtype_by_CM_raw.csv")
    zscore_column(raw_L).to_csv(TAB / "h_df_loading_cell_subtype_by_CM_zscore.csv")
    robust_column(raw_L).to_csv(TAB / "h_df_loading_cell_subtype_by_CM_robust.csv")
    standardize_column_scale(raw_L).to_csv(TAB / "h_df_loading_cell_subtype_by_CM_standard_scale_col.csv")

    # Loading heatmaps: never cluster CM axis
    plot_heatmap(raw_L, FIG / "h_df_loading_cell_subtype_by_CM_raw",
                 "CM loading (H) per cell_subtype", "viridis", "H loading")
    plot_heatmap(pd.read_csv(TAB / "h_df_loading_cell_subtype_by_CM_zscore.csv", index_col=0),
                 FIG / "h_df_loading_cell_subtype_by_CM_zscore",
                 "CM loading (H, z-score per CM column)", "RdBu_r", "z-score")
    plot_heatmap(pd.read_csv(TAB / "h_df_loading_cell_subtype_by_CM_robust.csv", index_col=0),
                 FIG / "h_df_loading_cell_subtype_by_CM_robust",
                 "CM loading (H, robust per CM column)", "RdBu_r", "robust Z")
    plot_heatmap(pd.read_csv(TAB / "h_df_loading_cell_subtype_by_CM_standard_scale_col.csv", index_col=0),
                 FIG / "h_df_loading_cell_subtype_by_CM_standard_scale_col",
                 "CM loading (H, column min-max)", "viridis", "min-max [0,1]")

    # Also update the legacy cm_loading_heatmap to standard_scale_col variant for consistency
    # (existing cm_loading_heatmap may be raw; keep the new one canonical and rename the old if needed)
    print("[main] done.")


if __name__ == "__main__":
    main()