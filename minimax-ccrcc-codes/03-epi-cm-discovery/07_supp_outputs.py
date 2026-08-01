"""Block 03 supplementary outputs per SKILL.md.

  1. joint_CM_activity_tumor_vs_normal_mean_sd_summary.csv
  2. Rename Spearman/Pearson association files to canonical names
  3. Tumor-centric edge-origin nodeplots (per-CM + all-CM)
  4. All-pair scatter plots (154 epi x CM) per method per status
"""
from __future__ import annotations

import re
import shutil
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


def safe_filename(s):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s))


def save_pdf_svg(fig, stem):
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", dpi=300)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", dpi=300)
    plt.close(fig)


def circular_rank_layout(nodes):
    n = len(nodes)
    return {node: (float(np.cos(2 * np.pi * i / max(n, 1))),
                   float(np.sin(2 * np.pi * i / max(n, 1))))
            for i, node in enumerate(nodes)}


def prefix_color_map(nodes):
    prefixes = sorted({str(node).split("_")[0] for node in nodes})
    palette = plt.get_cmap("tab20").colors
    return {prefix: palette[i % len(palette)] for i, prefix in enumerate(prefixes)}


def main():
    t0 = time.time()
    # ---- 1. CM activity mean/SD summary ----
    activity = pd.read_csv(TAB / "activity_df_sample_by_CM.csv", index_col=0)
    rows = []
    for cm in [c for c in activity.columns if c not in ("status", "non_epi_cells")]:
        for status_value in ["tumor", "normal-like"]:
            sub = activity[activity["status"] == status_value]
            vals = sub[cm].astype(float).values
            rows.append({
                "CM": cm, "status": status_value, "n_samples": int(len(vals)),
                "mean": float(np.mean(vals)),
                "sd": float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0,
                "median": float(np.median(vals)),
                "min": float(vals.min()),
                "max": float(vals.max()),
            })
    pd.DataFrame(rows).to_csv(TAB / "joint_CM_activity_tumor_vs_normal_mean_sd_summary.csv", index=False)
    print(f"[supp] joint_CM_activity_tumor_vs_normal_mean_sd_summary.csv saved", flush=True)

    # ---- 2. Rename Spearman/Pearson files to canonical names ----
    canonical_map = {
        # Spearman
        "epi_cm_association_overall_rho_matrix.csv": "balanced_joint_cm_epi_cm_association_overall_rho_matrix.csv",
        "epi_cm_association_overall_q_matrix.csv": "balanced_joint_cm_epi_cm_association_overall_q_matrix.csv",
        "epi_cm_association_overall_long.csv": "balanced_joint_cm_epi_cm_association_overall_spearman_long.csv",
        "epi_cm_association_tumor_rho_matrix.csv": "balanced_joint_cm_epi_cm_association_tumor_rho_matrix.csv",
        "epi_cm_association_tumor_q_matrix.csv": "balanced_joint_cm_epi_cm_association_tumor_q_matrix.csv",
        "epi_cm_association_tumor_long.csv": "balanced_joint_cm_epi_cm_association_tumor_spearman_long.csv",
        "epi_cm_association_tumor_significant_q0.05.csv": "balanced_joint_cm_epi_cm_association_tumor_spearman_significant_q0.05.csv",
        "epi_cm_association_normal-like_rho_matrix.csv": "balanced_joint_cm_epi_cm_association_normal-like_rho_matrix.csv",
        "epi_cm_association_normal-like_q_matrix.csv": "balanced_joint_cm_epi_cm_association_normal-like_q_matrix.csv",
        "epi_cm_association_normal-like_long.csv": "balanced_joint_cm_epi_cm_association_normal-like_spearman_long.csv",
        "epi_cm_association_normal-like_significant_q0.05.csv": "balanced_joint_cm_epi_cm_association_normal-like_spearman_significant_q0.05.csv",
        "epi_cm_association_overall_significant_q0.05.csv": "balanced_joint_cm_epi_cm_association_overall_spearman_significant_q0.05.csv",
        # Pearson
        "epi_cm_association_overall_pearson_r_matrix.csv": "balanced_joint_cm_epi_cm_association_overall_pearson_r_matrix.csv",
        "epi_cm_association_overall_pearson_q_matrix.csv": "balanced_joint_cm_epi_cm_association_overall_pearson_q_matrix.csv",
        "epi_cm_association_overall_pearson_long.csv": "balanced_joint_cm_epi_cm_association_overall_pearson_long.csv",
        "epi_cm_association_tumor_pearson_r_matrix.csv": "balanced_joint_cm_epi_cm_association_tumor_pearson_r_matrix.csv",
        "epi_cm_association_tumor_pearson_q_matrix.csv": "balanced_joint_cm_epi_cm_association_tumor_pearson_q_matrix.csv",
        "epi_cm_association_tumor_pearson_long.csv": "balanced_joint_cm_epi_cm_association_tumor_pearson_long.csv",
        "epi_cm_association_normal-like_pearson_r_matrix.csv": "balanced_joint_cm_epi_cm_association_normal-like_pearson_r_matrix.csv",
        "epi_cm_association_normal-like_pearson_q_matrix.csv": "balanced_joint_cm_epi_cm_association_normal-like_pearson_q_matrix.csv",
        "epi_cm_association_normal-like_pearson_long.csv": "balanced_joint_cm_epi_cm_association_normal-like_pearson_long.csv",
    }
    for old, new in canonical_map.items():
        src_sp = TAB / "association-spearman" / old
        src_pe = TAB / "association-pearson" / old
        if src_sp.exists():
            shutil.copy2(src_sp, TAB / "association-spearman" / new)
        if src_pe.exists():
            shutil.copy2(src_pe, TAB / "association-pearson" / new)
    print(f"[supp] canonical files copied", flush=True)

    # ---- 3. Tumor-centric edge-origin nodeplots ----
    ref_nodes = pd.read_csv(TAB / "balanced_joint_cm_reference_node_sets_after_edge_threshold.csv")
    edges = pd.read_csv(TAB / "status_specific_nodeplot_edges.csv")
    H = pd.read_csv(TAB / "H_df.csv", index_col=0)
    canonical_cms = H.index.astype(str).tolist()
    all_nodes = ref_nodes["node"].astype(str).tolist()
    node_colors = prefix_color_map(all_nodes)
    edge_norm = Normalize(vmin=EDGE_THRESHOLD, vmax=1.0)

    tumor_passing = edges[(edges["context"] == "tumor") & edges["edge_pass_r_ge_0.25"]].copy()
    # mark edge_origin: shared if same unordered pair also passes in normal-like
    normal_passing_pairs = set()
    for _, e in edges[(edges["context"] == "normal-like") & edges["edge_pass_r_ge_0.25"]].iterrows():
        normal_passing_pairs.add((e["CM"], tuple(sorted((e["node_a"], e["node_b"])))))
    tumor_passing["edge_origin"] = tumor_passing.apply(
        lambda r: "shared" if (r["CM"], tuple(sorted((r["node_a"], r["node_b"])))) in normal_passing_pairs
                     else "tumor_only", axis=1
    )
    tumor_passing.to_csv(TAB / "tumor_centric_edges_with_origin.csv", index=False)
    print(f"[supp] tumor-centric edges with origin saved: {len(tumor_passing)}", flush=True)

    tumor_dir = FIG / "tumor_centric_nodeplots_by_cm"
    tumor_dir.mkdir(parents=True, exist_ok=True)
    # all-CM overview grid
    def closest_square(n):
        if n <= 0:
            return 0, 0
        best = None
        for nr in range(1, int(np.ceil(np.sqrt(n))) + 1):
            nc = int(np.ceil(n / nr))
            cand = (abs(nc - nr), nr * nc - n, nr, nc)
            if best is None or cand < best:
                best = cand
        return best[2], best[3]
    n_rows, n_cols = closest_square(len(canonical_cms))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 5, n_rows * 5))
    axes_flat = np.atleast_1d(axes).flatten()
    FIXED_LW = 2.0  # FIXED edge width for tumor-centric nodeplots per skill
    for idx, cm in enumerate(canonical_cms):
        ax = axes_flat[idx]
        rnodes = ref_nodes.loc[ref_nodes["CM"] == cm, "node"].astype(str).tolist()
        sub_e = tumor_passing[tumor_passing["CM"] == cm]
        g = nx.Graph()
        g.add_nodes_from(rnodes)
        pos = circular_rank_layout(rnodes)
        for _, e in sub_e.iterrows():
            g.add_edge(e["node_a"], e["node_b"], weight=float(e["pearson_r"]), origin=e["edge_origin"])
        for u, v, d in g.edges(data=True):
            w = d["weight"]
            c = plt.cm.Reds(0.4 + 0.5 * w) if d["origin"] == "tumor_only" else plt.cm.Purples(0.4 + 0.5 * w)
            nx.draw_networkx_edges(g, pos, edgelist=[(u, v)], edge_color=[c], width=FIXED_LW, ax=ax, alpha=0.85)
        nc = [node_colors.get(str(n).split("_")[0], (0.5, 0.5, 0.5, 1.0)) for n in rnodes]
        ns = [200 + 100 * (i + 1) for i in range(len(rnodes))]
        nx.draw_networkx_nodes(g, pos, nodelist=rnodes, node_color=nc, node_size=ns, ax=ax, edgecolors="black", linewidths=0.5)
        nx.draw_networkx_labels(g, pos, labels={n: n for n in rnodes}, font_size=7, ax=ax)
        ax.set_title(f"{cm} (tumor-centric edge-origin)\n(n={len(sub_e)} passing)", fontsize=9)
        ax.set_axis_off()
    for j in range(len(canonical_cms), len(axes_flat)):
        axes_flat[j].set_axis_off()
    plt.tight_layout()
    save_pdf_svg(fig, FIG / "tumor_centric_nodeplot_edge_origin")

    # per-CM single panel
    for cm in canonical_cms:
        fig, ax = plt.subplots(figsize=(6, 6))
        rnodes = ref_nodes.loc[ref_nodes["CM"] == cm, "node"].astype(str).tolist()
        sub_e = tumor_passing[tumor_passing["CM"] == cm]
        g = nx.Graph()
        g.add_nodes_from(rnodes)
        pos = circular_rank_layout(rnodes)
        for _, e in sub_e.iterrows():
            g.add_edge(e["node_a"], e["node_b"], weight=float(e["pearson_r"]), origin=e["edge_origin"])
        for u, v, d in g.edges(data=True):
            w = d["weight"]
            c = plt.cm.Reds(0.4 + 0.5 * w) if d["origin"] == "tumor_only" else plt.cm.Purples(0.4 + 0.5 * w)
            nx.draw_networkx_edges(g, pos, edgelist=[(u, v)], edge_color=[c], width=FIXED_LW, ax=ax, alpha=0.85)
        nc = [node_colors.get(str(n).split("_")[0], (0.5, 0.5, 0.5, 1.0)) for n in rnodes]
        ns = [200 + 100 * (i + 1) for i in range(len(rnodes))]
        nx.draw_networkx_nodes(g, pos, nodelist=rnodes, node_color=nc, node_size=ns, ax=ax, edgecolors="black", linewidths=0.5)
        nx.draw_networkx_labels(g, pos, labels={n: n for n in rnodes}, font_size=7, ax=ax)
        ax.set_title(f"{cm} (tumor-centric)\n(n={len(sub_e)} passing edges)", fontsize=9)
        ax.set_axis_off()
        plt.tight_layout()
        save_pdf_svg(fig, tumor_dir / f"{cm}_tumor_centric_nodeplot_edge_origin")

    # legend figure
    fig, ax = plt.subplots(figsize=(5, 1.5))
    handles = [Line2D([0], [0], color=plt.cm.Reds(0.65), lw=2, label="tumor_only"),
               Line2D([0], [0], color=plt.cm.Purples(0.65), lw=2, label="shared (tumor ∩ normal)")]
    ax.legend(handles=handles, loc="center", ncol=2, frameon=False)
    ax.set_axis_off()
    save_pdf_svg(fig, FIG / "tumor_centric_nodeplot_edge_origin_legend")
    print(f"[supp] tumor-centric nodeplots saved", flush=True)

    # ---- 4. All-pair Epi-CM scatter per method per status ----
    epi_freq = pd.read_csv(TAB / "epi_subtype_frequency.csv", index_col=0)
    W = pd.read_csv(TAB / "W_df.csv", index_col=0)
    sample_status = pd.read_csv(TAB / "sample_status.csv", index_col=0)["status"]
    incl = pd.read_csv(TAB / "sample_inclusion_exclusion.csv", index_col=0)
    keep_samples = incl.index[incl["keep_for_epi_cm"]].astype(str).tolist()
    keep_samples = [s for s in keep_samples if s in epi_freq.index and s in W.index and s in sample_status.index]
    epi_color_map = prefix_color_map(list(epi_freq.columns))
    # validate palette coverage
    missing_palette = [c for c in epi_freq.columns if c not in epi_color_map]
    if missing_palette:
        for c in missing_palette:
            epi_color_map[c] = plt.get_cmap("tab20").colors[len(epi_color_map) % 20]

    for method, stat_label, dir_name in [("spearman", "rho", "epi-cm-association-spearman"),
                                        ("pearson", "r", "epi-cm-association-pearson")]:
        scatter_root = FIG / dir_name
        for status_value in ["tumor", "normal-like"]:
            sc_dir = scatter_root / f"scatter_{status_value.replace('-', '_')}"
            sc_dir.mkdir(parents=True, exist_ok=True)
            samples_ctx = [s for s in keep_samples
                           if s in W.index and s in epi_freq.index
                           and sample_status.loc[s] == status_value]
            for epi in epi_freq.columns:
                for cm in W.columns:
                    x = W.loc[samples_ctx, cm].astype(float).values
                    y = epi_freq.loc[samples_ctx, epi].astype(float).values
                    if method == "spearman":
                        r_, p_ = spearmanr(x, y)
                    else:
                        r_, p_ = pearsonr(x, y)
                    fig, ax = plt.subplots(figsize=(5, 4))
                    ax.scatter(x, y, s=20, color=epi_color_map.get(epi, (0.5, 0.5, 0.5)),
                               edgecolors="none", alpha=0.65)
                    sns.regplot(x=x, y=y, ax=ax, scatter=False, color="black", line_kws={"lw": 1})
                    ax.set_xlabel(f"{cm} score (W)")
                    ax.set_ylabel(f"{epi} fraction")
                    ax.set_title(f"{status_value}: {cm} vs {epi}\n{method} {stat_label}={r_:.3f}, p={p_:.2e}, n={len(samples_ctx)}", fontsize=8)
                    plt.tight_layout()
                    save_pdf_svg(fig, sc_dir / f"scatter_{safe_filename(cm)}_vs_{safe_filename(epi)}")
            plt.close("all")
        print(f"[supp] {method}: all-pair scatter saved (308 per method)", flush=True)

    print(f"[supp] done. total={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()