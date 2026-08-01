"""Block 03 final plotting per SKILL.md (canonical).

Produces:
  - joint_nmf_k_selection.{pdf,svg}
  - CM activity heatmap
  - CM loading heatmap
  - tumor_all_CM_nodeplot.{pdf,svg} + per-CM
  - normal_like_all_CM_nodeplot.{pdf,svg} + per-CM
  - top-node correlation heatmaps
  - epi-cm-association-spearman/* heatmaps and scatter
  - epi-cm-association-pearson/* heatmaps and scatter
"""
from __future__ import annotations

import json
import re
import time
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import networkx as nx
import seaborn as sns
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from scipy.stats import spearmanr, pearsonr

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.family"] = "DejaVu Sans"
mpl.rcParams["axes.linewidth"] = 0.8
sns.set_style("ticks")

SEED = 42
EDGE_THRESHOLD = 0.25
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
TAB = ROOT / "epi-cm-core-workflow/tables/03-epi-cm-discovery"
FIG = ROOT / "epi-cm-core-workflow/figures/03-epi-cm-discovery"
FIG.mkdir(parents=True, exist_ok=True)
(FIG / "tumor_nodeplots_by_cm").mkdir(parents=True, exist_ok=True)
(FIG / "normal_like_nodeplots_by_cm").mkdir(parents=True, exist_ok=True)
(FIG / "epi-cm-association-spearman").mkdir(parents=True, exist_ok=True)
(FIG / "epi-cm-association-spearman/scatter_tumor").mkdir(parents=True, exist_ok=True)
(FIG / "epi-cm-association-spearman/scatter_normal_like").mkdir(parents=True, exist_ok=True)
(FIG / "epi-cm-association-pearson").mkdir(parents=True, exist_ok=True)
(FIG / "epi-cm-association-pearson/scatter_tumor").mkdir(parents=True, exist_ok=True)
(FIG / "epi-cm-association-pearson/scatter_normal_like").mkdir(parents=True, exist_ok=True)


def save_pdf_svg(fig, stem: Path):
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", dpi=300)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", dpi=300)
    plt.close(fig)


def safe_filename(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s))


def closest_square_grid(n: int) -> tuple:
    if n <= 0:
        return 0, 0
    best = None
    for n_rows in range(1, int(np.ceil(np.sqrt(n))) + 1):
        n_cols = int(np.ceil(n / n_rows))
        empty = n_rows * n_cols - n
        aspect_gap = abs(n_cols - n_rows)
        candidate = (aspect_gap, empty, n_rows, n_cols)
        if best is None or candidate < best:
            best = candidate
    return best[2], best[3]


def circular_rank_layout(nodes: list[str]) -> dict:
    n = len(nodes)
    return {
        node: (float(np.cos(2 * np.pi * i / max(n, 1))),
               float(np.sin(2 * np.pi * i / max(n, 1))))
        for i, node in enumerate(nodes)
    }


def prefix_color_map(nodes: list[str]) -> dict:
    prefixes = sorted({str(node).split("_")[0] for node in nodes})
    palette = plt.get_cmap("tab20").colors
    return {prefix: palette[i % len(palette)] for i, prefix in enumerate(prefixes)}


def q_stars(q):
    if not np.isfinite(q):
        return "ns"
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return "ns"


def _fmt_p(p):
    """FIXED per SKILL scatter rule: format p/q for the statistic box."""
    if not np.isfinite(p):
        return "NA"
    return f"{p:.1e}" if p < 0.001 else f"{p:.3f}"


def draw_cm_nodeplot(ax, cm, cm_nodes, cm_edges_pass, node_colors, edge_norm, cmap, mode):
    g = nx.Graph()
    g.add_nodes_from(cm_nodes)
    pos = circular_rank_layout(cm_nodes)
    for _, e in cm_edges_pass.iterrows():
        g.add_edge(e["node_a"], e["node_b"], weight=float(e["pearson_r"]))
    for u, v, d in g.edges(data=True):
        w = d["weight"]
        c = cmap(edge_norm(w))
        lw = 0.5 + 3.0 * abs(w)
        nx.draw_networkx_edges(g, pos, edgelist=[(u, v)], edge_color=[c], width=lw, ax=ax, alpha=0.8)
    nc = [node_colors.get(str(n).split("_")[0], (0.5, 0.5, 0.5, 1.0)) for n in cm_nodes]
    ns = [200 + 100 * (i + 1) for i in range(len(cm_nodes))]
    nx.draw_networkx_nodes(g, pos, nodelist=cm_nodes, node_color=nc, node_size=ns, ax=ax, edgecolors="black", linewidths=0.5)
    nx.draw_networkx_labels(g, pos, labels={n: n for n in cm_nodes}, font_size=7, ax=ax)
    ax.set_title(f"{cm}\n({mode}, n={len(cm_edges_pass)} passing edges)", fontsize=9)
    ax.set_axis_off()


def heatmap_with_qstars(ax, rho, qval, title, stat_label="Spearman rho"):
    vmax = float(np.nanmax(np.abs(rho.values))) if rho.size else 1.0
    im = ax.imshow(rho.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(rho.shape[1]))
    ax.set_xticklabels(rho.columns, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(rho.shape[0]))
    ax.set_yticklabels(rho.index, fontsize=7)
    ax.set_title(title, fontsize=9)
    for i in range(rho.shape[0]):
        for j in range(rho.shape[1]):
            q = qval.iloc[i, j]
            s = q_stars(q)
            if s != "ns":
                ax.text(j, i, s, ha="center", va="center", color="black", fontsize=7, fontweight="bold")
    return im


def main() -> None:
    t0 = time.time()
    H = pd.read_csv(TAB / "H_df.csv", index_col=0)
    W = pd.read_csv(TAB / "W_df.csv", index_col=0)
    activity = pd.read_csv(TAB / "activity_df_sample_by_CM.csv", index_col=0)
    loading = pd.read_csv(TAB / "loading_df_cell_subtype_by_CM.csv", index_col=0)
    classification = pd.read_csv(TAB / "joint_module_classification.csv")
    ref_nodes = pd.read_csv(TAB / "balanced_joint_cm_reference_node_sets_after_edge_threshold.csv")
    edges = pd.read_csv(TAB / "status_specific_nodeplot_edges.csv")
    metrics = pd.read_csv(TAB / "joint_nmf_k_selection_metrics.csv")
    canonical_cms = H.index.astype(str).tolist()
    print(f"[plot] loaded inputs, time={time.time()-t0:.1f}s", flush=True)

    # ---- 1. K-selection figure ----
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.2), constrained_layout=False)
    panels = [
        ("best_balanced_explained_fraction", "Best explained", "#1f77b4"),
        ("stability_matched_cosine", "Stability (matched cosine)", "#2ca02c"),
        ("selection_score", "Selection score", "#d62728"),
    ]
    for ax, (col, title, color) in zip(axes, panels):
        x = metrics["k"].astype(int).values
        y = metrics[col].astype(float).values
        ax.plot(x, y, marker="o", markersize=4, linewidth=1.6, color=color)
        ax.set_title(title, fontweight="bold", fontsize=10)
        ax.set_xlabel("K", fontsize=9)
        ax.set_xticks(x)
        ax.grid(axis="y", color="#d9d9d9", linewidth=0.6, alpha=0.8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    best_k = int(metrics.loc[metrics["selected"], "k"].iloc[0])
    for ax, _ in zip(axes, panels):
        ax.axvline(best_k, color="black", linestyle="--", linewidth=1)
        ax.text(best_k, ax.get_ylim()[1] * 0.95, f"K={best_k}", ha="center", va="top", fontsize=9)
    plt.tight_layout()
    save_pdf_svg(fig, FIG / "joint_nmf_k_selection")
    print(f"[plot] K_selection saved (K={best_k})", flush=True)

    # ---- 2. CM activity heatmap (sample x CM, by status) ----
    act_for_plot = activity[["status"] + canonical_cms].copy()
    sample_order = sorted(act_for_plot.index, key=lambda s: (act_for_plot.loc[s, "status"], s))
    act_for_plot = act_for_plot.loc[sample_order]
    fig, ax = plt.subplots(figsize=(max(7, len(canonical_cms) * 0.5), 8))
    sns.heatmap(act_for_plot[canonical_cms].astype(float), cmap="viridis",
                cbar_kws={"label": "CM activity (W)"}, ax=ax)
    # Add a row-color bar on the left
    for i, s in enumerate(act_for_plot.index):
        c = "#d62728" if act_for_plot.loc[s, "status"] == "tumor" else "#1f77b4"
        ax.add_patch(plt.Rectangle((-1, i), 0.3, 1, facecolor=c, clip_on=False))
    ax.set_title("CM activity (W) by sample\n(red=tumor, blue=normal-like)", fontsize=10)
    plt.tight_layout()
    save_pdf_svg(fig, FIG / "cm_activity_heatmap")
    print(f"[plot] cm_activity_heatmap saved", flush=True)

    # ---- 3. CM loading heatmap (CM x subtype) ----
    fig, ax = plt.subplots(figsize=(12, max(4, len(canonical_cms) * 0.45)))
    sns.heatmap(H.astype(float), cmap="viridis", cbar_kws={"label": "Loading (H)"}, ax=ax)
    ax.set_title("CM loading (H): CM x non-epithelial subtype", fontsize=10)
    plt.tight_layout()
    save_pdf_svg(fig, FIG / "cm_loading_heatmap")
    print(f"[plot] cm_loading_heatmap saved", flush=True)

    # ---- 4. Per-CM nodeplots: tumor + normal-like ----
    edge_norm = Normalize(vmin=EDGE_THRESHOLD, vmax=1.0)
    cmap = plt.get_cmap("viridis")
    all_nodes = ref_nodes["node"].astype(str).tolist()
    node_colors = prefix_color_map(all_nodes)
    for fctx, label in [("tumor", "tumor"), ("normal-like", "normal_like")]:
        # all-CM overview grid
        n_rows, n_cols = closest_square_grid(len(canonical_cms))
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 5, n_rows * 5))
        axes_flat = np.atleast_1d(axes).flatten()
        for idx, cm in enumerate(canonical_cms):
            rnodes = ref_nodes.loc[ref_nodes["CM"] == cm, "node"].astype(str).tolist()
            sub_e = edges[(edges["context"] == fctx) & (edges["CM"] == cm) & (edges["edge_pass_r_ge_0.25"])]
            draw_cm_nodeplot(axes_flat[idx], cm, rnodes, sub_e, node_colors, edge_norm, cmap, label)
        for j in range(len(canonical_cms), len(axes_flat)):
            axes_flat[j].set_axis_off()
        plt.tight_layout()
        save_pdf_svg(fig, FIG / f"{label}_all_CM_nodeplot")
        # per-CM single
        for cm in canonical_cms:
            fig, ax = plt.subplots(figsize=(6, 6))
            rnodes = ref_nodes.loc[ref_nodes["CM"] == cm, "node"].astype(str).tolist()
            sub_e = edges[(edges["context"] == fctx) & (edges["CM"] == cm) & (edges["edge_pass_r_ge_0.25"])]
            draw_cm_nodeplot(ax, cm, rnodes, sub_e, node_colors, edge_norm, cmap, label)
            plt.tight_layout()
            save_pdf_svg(fig, FIG / f"{label}_nodeplots_by_cm" / f"{cm}_{label}_nodeplot")
    print(f"[plot] per-CM nodeplots saved", flush=True)

    # FIXED: top10 node correlation heatmaps are produced by 11_top10_node_heatmaps.py
    # (canonical SKILL implementation). The legacy per-CM heatmaps that wrote to
    # node_correlation_heatmaps/ are removed because SKILL only requires the
    # top10_*-no_edge_filter outputs (multi-panel + per-CM).

    # ---- 6. Epi-CM heatmap (per branch) and scatter ----
    epi_freq = pd.read_csv(TAB / "epi_subtype_frequency.csv", index_col=0)
    sample_status = pd.read_csv(TAB / "sample_status.csv", index_col=0)["status"]
    incl = pd.read_csv(TAB / "sample_inclusion_exclusion.csv", index_col=0)
    keep_samples = incl.index[incl["keep_for_epi_cm"]].astype(str).tolist()
    keep_samples = [s for s in keep_samples if s in epi_freq.index and s in W.index and s in sample_status.index]
    epi_color_map = prefix_color_map(list(epi_freq.columns))

    for method, stat_label, dir_name in [("spearman", "rho", "epi-cm-association-spearman"),
                                        ("pearson", "r", "epi-cm-association-pearson")]:
        assoc_dir = TAB / f"association-{method}"
        rho_suffix = "rho" if method == "spearman" else "pearson_r"
        q_suffix = "q" if method == "spearman" else "pearson_q"
        # heatmap (overall)
        rho_mat = pd.read_csv(assoc_dir / f"epi_cm_association_overall_{rho_suffix}_matrix.csv", index_col=0)
        q_mat = pd.read_csv(assoc_dir / f"epi_cm_association_overall_{q_suffix}_matrix.csv", index_col=0)
        fig, ax = plt.subplots(figsize=(max(8, rho_mat.shape[1] * 0.5), max(6, rho_mat.shape[0] * 0.5)))
        im = heatmap_with_qstars(ax, rho_mat, q_mat,
                                  f"Epi subtype freq x CM activity: {method} {stat_label}\n(n=82 samples)",
                                  stat_label=f"{method} {stat_label}")
        plt.colorbar(im, ax=ax, fraction=0.05, pad=0.04, label=f"{method} {stat_label}")
        plt.tight_layout()
        save_pdf_svg(fig, FIG / dir_name / f"epi_cm_association_overall_{method}_heatmap_qstars")
        # per-context heatmaps
        for ctx in ["tumor", "normal-like"]:
            rho_c = pd.read_csv(assoc_dir / f"epi_cm_association_{ctx}_{rho_suffix}_matrix.csv", index_col=0)
            q_c = pd.read_csv(assoc_dir / f"epi_cm_association_{ctx}_{q_suffix}_matrix.csv", index_col=0)
            fig, ax = plt.subplots(figsize=(max(8, rho_c.shape[1] * 0.5), max(6, rho_c.shape[0] * 0.5)))
            im = heatmap_with_qstars(ax, rho_c, q_c,
                                      f"{ctx}: Epi x CM {method} {stat_label}", stat_label=f"{method} {stat_label}")
            plt.colorbar(im, ax=ax, fraction=0.05, pad=0.04, label=f"{method} {stat_label}")
            plt.tight_layout()
            save_pdf_svg(fig, FIG / dir_name / f"epi_cm_association_{ctx}_{method}_heatmap_qstars")
        # FIXED per updated SKILL: full Cartesian product of Epi x CM pairs.
        # Do not restrict to selected / representative / significant / top-ranked.
        long_file = "epi_cm_association_overall_long.csv" if method == "spearman" else "epi_cm_association_overall_pearson_long.csv"
        long_df = pd.read_csv(assoc_dir / long_file)
        qmat_for_box = pd.read_csv(assoc_dir / ("epi_cm_association_overall_q_matrix.csv" if method == "spearman" else "epi_cm_association_overall_pearson_q_matrix.csv"), index_col=0)
        for _, row in long_df.iterrows():
            epi, cm = row["epi_subtype"], row["CM"]
            for status_value in ["tumor", "normal-like"]:
                samples_ctx = [s for s in keep_samples
                               if s in W.index and s in epi_freq.index
                               and sample_status.loc[s] == status_value]
                if len(samples_ctx) < 5:
                    continue
                x = W.loc[samples_ctx, cm].astype(float).values
                y = epi_freq.loc[samples_ctx, epi].astype(float).values
                if method == "spearman":
                    r_, p_ = spearmanr(x, y)
                else:
                    r_, p_ = pearsonr(x, y)
                # FIXED per updated SKILL: canonical scatter style
                # - sns.set_style("ticks")
                # - canvas 6 x 5 inches
                # - y-axis fixed to 0-1
                # - point_size=60, alpha=0.85, no marker edge
                # - line width=2, ci=95 (default seaborn confidence band)
                # - white rounded statistic box at upper-left:
                #   "Spearman rho=..., q=..., n=..." (use r/ρ label per method)
                # - per-pair stem: scatter_<CM>_vs_<epi_subtype>.{pdf,svg}
                # - color points by epithelial subtype palette
                sns.set_style("ticks")
                fig, ax = plt.subplots(figsize=(6, 5))
                sns.regplot(
                    x=x, y=y, ax=ax,
                    scatter_kws={"s": 60, "alpha": 0.85,
                                 "color": epi_color_map.get(epi, "#666666"),
                                 "edgecolors": "none"},
                    line_kws={"lw": 2, "color": "black"},
                    ci=95,
                )
                ax.set_ylim(0, 1)
                ax.set_xlabel(f"{cm} score (W)")
                ax.set_ylabel(f"{epi} fraction")
                test_label = "Spearman" if method == "spearman" else "Pearson"
                # FIXED: use q from the per-status summary (already computed above).
                q_val = qmat_for_box.loc[epi, cm] if (epi in qmat_for_box.index and cm in qmat_for_box.columns) else np.nan
                stat_text = (f"{test_label} ρ={r_:.3f}\nq={_fmt_p(q_val)}\nn={len(samples_ctx)}")
                ax.text(0.03, 0.97, stat_text, transform=ax.transAxes,
                        ha="left", va="top", fontsize=9,
                        bbox=dict(boxstyle="round,pad=0.4",
                                  facecolor="white", edgecolor="gray", alpha=0.9))
                fig.tight_layout()
                save_pdf_svg(fig, FIG / dir_name / f"scatter_{status_value.replace('-','_')}" / f"scatter_{safe_filename(cm)}_vs_{safe_filename(epi)}")
    print(f"[plot] Epi-CM heatmaps and scatters saved", flush=True)

    print(f"[plot] done. total={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()