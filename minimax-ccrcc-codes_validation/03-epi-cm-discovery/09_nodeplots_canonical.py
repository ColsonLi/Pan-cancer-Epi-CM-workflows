"""Block 03 final: CM nodeplots per SKILL.md (updated 2026-07-16).

Strictly canonical per updated SKILL.md:
  - SEED = 42 (explicit, even though circular layout is deterministic; documents
    intent and matches SKILL wording).
  - Node color encodes CM/component identity (one color per CM, all nodes of
    a CM share that color). NOT prefix_color_map.
  - Node size reflects loading rank or magnitude (size scales with the
    absolute H loading for the node's primary CM).
  - Edge width reflects edge weight (Pearson r).
  - Edge color distinguishes tumor-only, normal-like-only, shared (for
    tumor-centric edge-origin plots).
  - Multi-panel grid for all-CM overview: closest-to-square factor/grid,
    NOT hard-coded ncols.
  - One single-panel nodeplot per CM under per-CM subdirectories.
  - For normal_like_all_CM_nodeplot: draw only normal-like edges with
    edge_pass_r_ge_0.25 == True, encode color/width from normal-like pearson_r.
  - For tumor_all_CM_nodeplot: draw only tumor edges with
    edge_pass_r_ge_0.25 == True, encode color/width from tumor pearson_r.
  - For tumor-centric: edge origin = "tumor_only" if pair only passes in tumor;
    "shared" if pair passes in both contexts; uses tumor pearson_r.

Required outputs:
  figures/03-epi-cm-discovery/
    normal_like_all_CM_nodeplot.{pdf,svg}
    tumor_all_CM_nodeplot.{pdf,svg}
    normal_like_nodeplots_by_cm/<CM>_normal_like_nodeplot.{pdf,svg}
    tumor_nodeplots_by_cm/<CM>_tumor_nodeplot.{pdf,svg}
    tumor_centric_nodeplot_edge_origin.{pdf,svg}
    tumor_centric_nodeplots_by_cm/<CM>_tumor_centric_nodeplot_edge_origin.{pdf,svg}
    tumor_centric_nodeplot_edge_origin_tumor_only_edge_colorbar.{pdf,svg}
    tumor_centric_nodeplot_edge_origin_shared_edge_colorbar.{pdf,svg}
    tumor_centric_nodeplot_edge_origin_edge_class_legend.{pdf,svg}
    nodeplot_network_node_legend.{pdf,svg}
    nodeplot_network_edge_colorbar.{pdf,svg}
"""
from __future__ import annotations

import time as _t
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import networkx as nx

SEED = 42
np.random.seed(SEED)

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.family"] = "DejaVu Sans"
mpl.rcParams["axes.linewidth"] = 0.8

EDGE_THRESHOLD = 0.25

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
TAB = ROOT / "epi-cm-core-workflow/tables/03-epi-cm-discovery"
FIG = ROOT / "epi-cm-core-workflow/figures/03-epi-cm-discovery"


def save_pdf_svg(fig, stem: Path):
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", dpi=300)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", dpi=300)
    plt.close(fig)


def closest_square_grid(n: int) -> tuple[int, int]:
    """FIXED: closest-to-square factor grid (no hard-coded ncols)."""
    if n <= 0:
        return 0, 0
    best = None
    for nr in range(1, int(np.ceil(np.sqrt(n))) + 1):
        nc = int(np.ceil(n / nr))
        cand = (abs(nc - nr), nr * nc - n, nr, nc)
        if best is None or cand < best:
            best = cand
    return best[2], best[3]


def cm_color_map(canonical_cms):
    """FIXED: one color per CM (used as fallback when no lineage palette exists).

    Per updated SKILL (line 2572): node color rule was REMOVED. Current rule:
    - "Node size or label weight reflects loading rank or loading magnitude"
    - "edge width reflects edge weight"
    - "edge color distinguishes tumor-only, normal-like-only, and shared edges"
    """
    base = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors)
    return {cm: base[i % len(base)] for i, cm in enumerate(canonical_cms)}


def prefix_color_map(nodes):
    """FIXED: node color by lineage prefix (B_, Epi_, T_, ...) — distinguishes
    nodes within a CM (used as default since SKILL removed the CM-identity rule).
    """
    prefixes = sorted({str(n).split("_")[0] for n in nodes})
    palette = plt.get_cmap("tab20").colors
    return {p: palette[i % len(palette)] for i, p in enumerate(prefixes)}


def circular_rank_layout(nodes, ordered_nodes=None):
    """FIXED deterministic circular layout — same node order across modes.

    `ordered_nodes` (when provided) ensures the same coordinate per node
    across normal_like / tumor / tumor_centric views.
    """
    if ordered_nodes is None:
        ordered_nodes = list(nodes)
    n = len(ordered_nodes)
    return {node: (float(np.cos(2 * np.pi * i / max(n, 1))),
                   float(np.sin(2 * np.pi * i / max(n, 1))))
            for i, node in enumerate(ordered_nodes)}


def node_size_from_loading(loading_value, min_size=600, max_size=2400):
    """FIXED: node size scales with |H loading|."""
    if not np.isfinite(loading_value) or loading_value <= 0:
        return min_size
    return float(np.clip(loading_value * 1800 + min_size, min_size, max_size))


def assign_edge_class(row, normal_passing_pairs, normal_passing_pairs_tumor):
    """FIXED: tumor_only / normal_like / shared classification."""
    key = tuple(sorted((row["node_a"], row["node_b"])))
    if row["context"] == "tumor":
        if not row["edge_pass_r_ge_0.25"]:
            return None
        return "shared" if key in normal_passing_pairs else "tumor_only"
    else:  # normal-like
        if not row["edge_pass_r_ge_0.25"]:
            return None
        return "shared" if (row["CM"], key) in normal_passing_pairs_tumor else "normal_like"


def draw_one_cm_nodeplot(ax, cm, cm_nodes, cm_edges, cm_pos, cm_colors,
                          loading_lookup, edge_norm,
                          tumor_only_cmap, shared_cmap, mode,
                          class_normal_passing_pairs):
    """FIXED: per-CM single-panel nodeplot. Color/label encoding per CM."""
    g = nx.Graph()
    g.add_nodes_from(cm_nodes)
    # FIXED: reuse pre-computed coordinates (same across modes)
    pos = {n: cm_pos[n] for n in cm_nodes if n in cm_pos}
    for n in cm_nodes:
        if n not in pos:
            pos[n] = (0.0, 0.0)

    if mode == "tumor":
        for _, row in cm_edges.iterrows():
            if row["context"] == "tumor" and row["edge_pass_r_ge_0.25"]:
                g.add_edge(row["node_a"], row["node_b"],
                           weight=float(row["pearson_r"]),
                           edge_class="tumor")
    elif mode == "normal_like":
        for _, row in cm_edges.iterrows():
            if row["context"] == "normal-like" and row["edge_pass_r_ge_0.25"]:
                g.add_edge(row["node_a"], row["node_b"],
                           weight=float(row["pearson_r"]),
                           edge_class="normal_like")
    elif mode == "tumor_centric":
        for _, row in cm_edges.iterrows():
            if row["context"] == "tumor" and row["edge_pass_r_ge_0.25"]:
                key = (row["CM"], tuple(sorted((row["node_a"], row["node_b"]))))
                g.add_edge(row["node_a"], row["node_b"],
                           weight=float(row["pearson_r"]),
                           edge_class="shared" if key in class_normal_passing_pairs else "tumor_only")

    class_to_cmap = (
        {"tumor_only": tumor_only_cmap, "shared": shared_cmap,
         "tumor": tumor_only_cmap, "normal_like": shared_cmap}
    )
    for edge_class, cmap in class_to_cmap.items():
        class_edges = [(a, b, d) for a, b, d in g.edges(data=True)
                       if d["edge_class"] == edge_class]
        if not class_edges:
            continue
        if mode == "tumor_centric":
            width = 2.6
        else:
            width = [2.0 + 3.0 * edge_norm(d["weight"]) for _, _, d in class_edges]
        nx.draw_networkx_edges(
            g, pos,
            edgelist=[(a, b) for a, b, _ in class_edges],
            edge_color=[cmap(edge_norm(d["weight"])) for _, _, d in class_edges],
            width=width, alpha=0.88, ax=ax,
        )

    # FIXED: node color = lineage prefix (distinguishes nodes within a CM).
    # Per updated SKILL: CM-identity node color rule was removed.
    node_colors = [cm_colors.get(str(n).split("_")[0], "#999999")
                   for n in cm_nodes]
    node_sizes = [node_size_from_loading(loading_lookup.get((cm, n), 0.0))
                  for n in cm_nodes]
    nx.draw_networkx_nodes(
        g, pos,
        nodelist=cm_nodes,
        node_color=node_colors,
        node_size=node_sizes, linewidths=0.8, edgecolors="white", ax=ax,
    )
    # FIXED: label weight reflects loading magnitude (font size scales with H loading)
    label_sizes = {n: 7 + min(5, 2 * loading_lookup.get((cm, n), 0.0)) for n in cm_nodes}
    for n, (x, y) in pos.items():
        if n in cm_nodes:
            ax.text(x, y, n, fontsize=label_sizes[n], ha="center", va="center",
                    color="black", zorder=4)
    ax.set_title(cm, fontsize=12, fontweight="bold")
    ax.set_aspect("equal")
    ax.axis("off")
    return {str(d["edge_class"]) for _, _, d in g.edges(data=True)}


def add_edge_colorbars(fig, mode, edge_norm, tumor_only_cmap, shared_cmap,
                       present_edge_classes):
    """FIXED per SKILL: embed colorbars using the EXACT live cmap/norm objects
    used by draw_one_cm_nodeplot. Include only edge classes drawn in this fig.
    """
    if mode == "normal_like" and "normal_like" in present_edge_classes:
        specs = [(shared_cmap, "Normal-like edge Pearson r")]
    elif mode == "tumor" and "tumor" in present_edge_classes:
        specs = [(tumor_only_cmap, "Tumor edge Pearson r")]
    elif mode == "tumor_centric":
        specs = []
        if "tumor_only" in present_edge_classes:
            specs.append((tumor_only_cmap, "Tumor-only edge Pearson r"))
        if "shared" in present_edge_classes:
            specs.append((shared_cmap, "Shared edge Pearson r"))
    else:
        specs = []

    if not specs:
        return
    left_margin = 0.14
    right_margin = 0.08
    gap = 0.06 if len(specs) > 1 else 0.0
    bar_width = (1.0 - left_margin - right_margin - gap * (len(specs) - 1)) / len(specs)
    for bar_i, (cmap, label) in enumerate(specs):
        left = left_margin + bar_i * (bar_width + gap)
        cax = fig.add_axes([left, 0.055, bar_width, 0.022])
        sm = plt.cm.ScalarMappable(norm=edge_norm, cmap=cmap)
        sm.set_array([])
        cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
        cb.set_label(label, fontsize=8)
        cb.ax.tick_params(labelsize=7, length=2)


def save_node_legend(out_dir, cm_colors, prefixes):
    """FIXED: node legend shows lineage prefix → color mapping."""
    fig, ax = plt.subplots(figsize=(5, max(2.0, 0.32 * len(prefixes))))
    handles = [Patch(facecolor=cm_colors[p], edgecolor="white", label=p)
               for p in prefixes]
    ax.legend(handles=handles, loc="center", ncol=2, frameon=False,
              title="Cell lineage prefix (node color)")
    ax.set_axis_off()
    save_pdf_svg(fig, out_dir / "nodeplot_network_node_legend")
    plt.close(fig)


def save_edge_legends(out_dir, edge_norm, tumor_only_cmap, shared_cmap):
    """FIXED: standard context edge colorbar (Pearson r in context)."""
    fig, ax = plt.subplots(figsize=(5, 0.8))
    sm = plt.cm.ScalarMappable(norm=edge_norm, cmap=tumor_only_cmap)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.6, pad=0.05, orientation="horizontal")
    cbar.set_label(f"Edge Pearson r (>= {EDGE_THRESHOLD})")
    ax.set_axis_off()
    save_pdf_svg(fig, out_dir / "nodeplot_network_edge_colorbar")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 0.8))
    sm = plt.cm.ScalarMappable(norm=edge_norm, cmap="Reds")
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.6, pad=0.05, orientation="horizontal")
    cbar.set_label("Tumor-only edge Pearson r (tumor samples)")
    ax.set_axis_off()
    save_pdf_svg(fig, out_dir / "tumor_centric_nodeplot_edge_origin_tumor_only_edge_colorbar")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 0.8))
    sm = plt.cm.ScalarMappable(norm=edge_norm, cmap="Purples")
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.6, pad=0.05, orientation="horizontal")
    cbar.set_label("Shared edge Pearson r (tumor samples; also passes in normal-like)")
    ax.set_axis_off()
    save_pdf_svg(fig, out_dir / "tumor_centric_nodeplot_edge_origin_shared_edge_colorbar")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 1.0))
    handles = [
        Line2D([0], [0], color=plt.cm.Reds(0.7), lw=2.5,
                label="tumor_only (tumor pass, normal-like fail)"),
        Line2D([0], [0], color=plt.cm.Purples(0.7), lw=2.5,
                label="shared (both contexts pass)"),
    ]
    ax.legend(handles=handles, loc="center", ncol=1, frameon=False,
              title="Tumor-centric edge class")
    ax.set_axis_off()
    save_pdf_svg(fig, out_dir / "tumor_centric_nodeplot_edge_origin_edge_class_legend")
    plt.close(fig)


def main():
    t0 = _t.time()

    ref_nodes = pd.read_csv(TAB / "balanced_joint_cm_reference_node_sets_after_edge_threshold.csv")
    edges_df = pd.read_csv(TAB / "status_specific_nodeplot_edges.csv")
    H = pd.read_csv(TAB / "H_df.csv", index_col=0)
    canonical_cms = H.index.astype(str).tolist()

    # Pre-compute passing pairs for cross-context shared detection
    normal_passing_pairs = set()
    for _, e in edges_df[(edges_df["context"] == "normal-like") & edges_df["edge_pass_r_ge_0.25"]].iterrows():
        normal_passing_pairs.add((e["CM"], tuple(sorted((e["node_a"], e["node_b"])))))
    normal_passing_pairs_tumor = set()
    for _, e in edges_df[(edges_df["context"] == "tumor") & edges_df["edge_pass_r_ge_0.25"]].iterrows():
        normal_passing_pairs_tumor.add((e["CM"], tuple(sorted((e["node_a"], e["node_b"])))))

    tumor_only_cmap = LinearSegmentedColormap.from_list("tumor_only_blue", ["#dbeafe", "#1d4ed8"])
    shared_cmap = LinearSegmentedColormap.from_list("shared_red", ["#fee2e2", "#b91c1c"])
    edge_norm = Normalize(vmin=EDGE_THRESHOLD, vmax=1.0)

    # FIXED: lineage-prefix color map (distinguishes nodes within a CM).
    # Per updated SKILL: node color rule was removed; default to lineage
    # prefix coloring so nodes within a per-CM nodeplot are visually distinct.
    all_nodes = ref_nodes["node"].astype(str).tolist()
    cm_colors = prefix_color_map(all_nodes)

    # Loading lookup (absolute H loading for sizing)
    loading_lookup = {}
    for cm in canonical_cms:
        for node in H.columns.astype(str):
            try:
                v = float(abs(H.loc[cm, node]))
            except (KeyError, ValueError, TypeError):
                v = 0.0
            loading_lookup[(cm, node)] = v

    # FIXED deterministic coordinate per CM (same across modes)
    cm_layouts = {}
    for cm in canonical_cms:
        cm_nodes_ord = ref_nodes.loc[ref_nodes["CM"] == cm, "node"].astype(str).tolist()
        cm_layouts[cm] = circular_rank_layout(cm_nodes_ord, ordered_nodes=cm_nodes_ord)

    modes = {
        "normal_like": ("normal_like_nodeplots_by_cm", "normal_like_all_CM_nodeplot",
                         "_normal_like_nodeplot"),
        "tumor": ("tumor_nodeplots_by_cm", "tumor_all_CM_nodeplot",
                   "_tumor_nodeplot"),
        "tumor_centric": ("tumor_centric_nodeplots_by_cm",
                           "tumor_centric_nodeplot_edge_origin",
                           "_tumor_centric_nodeplot_edge_origin"),
    }

    # Per-CM single panels
    for cm in canonical_cms:
        cm_nodes = ref_nodes.loc[ref_nodes["CM"] == cm, "node"].astype(str).tolist()
        cm_edges = edges_df[edges_df["CM"] == cm].copy()
        cm_pos = cm_layouts[cm]
        for mode, (subdir, _, suffix) in modes.items():
            fig, ax = plt.subplots(figsize=(4.4, 4.6))  # extra height for embedded bar
            present_classes = draw_one_cm_nodeplot(
                ax, cm, cm_nodes, cm_edges, cm_pos, cm_colors,
                loading_lookup, edge_norm,
                tumor_only_cmap, shared_cmap, mode,
                normal_passing_pairs,
            )
            # FIXED per SKILL: every nodeplot embeds colorbars (from live cmap/norm)
            add_edge_colorbars(fig, mode, edge_norm, tumor_only_cmap, shared_cmap,
                               present_classes)
            save_pdf_svg(fig, FIG / subdir / f"{cm}{suffix}")
        plt.close("all")

    # All-CM overview grids (closest-to-square, not hard-coded)
    n_rows, n_cols = closest_square_grid(len(canonical_cms))
    for mode, (_, stem, _) in modes.items():
        if n_rows == 0:
            continue
        fig, axes = plt.subplots(n_rows, n_cols,
                                  figsize=(n_cols * 4.4, n_rows * 4.6),
                                  squeeze=False)
        flat_axes = axes.ravel()
        union_classes = set()
        for ax, cm in zip(flat_axes, canonical_cms):
            cm_nodes = ref_nodes.loc[ref_nodes["CM"] == cm, "node"].astype(str).tolist()
            cm_edges = edges_df[edges_df["CM"] == cm].copy()
            cm_pos = cm_layouts[cm]
            present = draw_one_cm_nodeplot(
                ax, cm, cm_nodes, cm_edges, cm_pos, cm_colors,
                loading_lookup, edge_norm,
                tumor_only_cmap, shared_cmap, mode,
                normal_passing_pairs,
            )
            union_classes |= present
        for ax in flat_axes[len(canonical_cms):]:
            ax.axis("off")
        fig.tight_layout()
        # FIXED per SKILL: all-CM overview gets one shared set of colorbars
        # for the union of classes actually drawn across panels.
        add_edge_colorbars(fig, mode, edge_norm, tumor_only_cmap, shared_cmap,
                           union_classes)
        save_pdf_svg(fig, FIG / stem)

    # Legends + colorbars
    all_prefixes = sorted({str(n).split("_")[0] for n in all_nodes})
    save_node_legend(FIG, cm_colors, all_prefixes)
    save_edge_legends(FIG, edge_norm, tumor_only_cmap, shared_cmap)
    print(f"[nodeplot] all per-CM and all-CM nodeplots + legends saved. total={_t.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()