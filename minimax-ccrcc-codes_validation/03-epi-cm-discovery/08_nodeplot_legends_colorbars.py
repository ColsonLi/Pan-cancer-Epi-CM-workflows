"""Block 03 final: nodeplot legends and edge colorbars (canonical SKILL.md).

Per skill, write:
  - nodeplot_network_node_legend.{pdf,svg} (node color legend)
  - nodeplot_network_edge_colorbar.{pdf,svg} (edge color bar)
  - tumor_centric_nodeplot_edge_origin_tumor_only_edge_colorbar.{pdf,svg}
  - tumor_centric_nodeplot_edge_origin_shared_edge_colorbar.{pdf,svg}
  - tumor_centric_nodeplot_edge_origin_edge_class_legend.{pdf,svg} (re-drawn with correct name)
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.cm import ScalarMappable

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.family"] = "DejaVu Sans"
mpl.rcParams["axes.linewidth"] = 0.8

EDGE_THRESHOLD = 0.25

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
FIG = ROOT / "epi-cm-core-workflow/figures/03-epi-cm-discovery"
TAB = ROOT / "epi-cm-core-workflow/tables/03-epi-cm-discovery"


def save_pdf_svg(fig, stem: Path):
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", dpi=300)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", dpi=300)
    plt.close(fig)


def main():
    # ---- 1. nodeplot_network_node_legend ----
    # Node color: by subtype prefix/lineage (stable palette)
    ref_nodes = __import__("pandas").read_csv(TAB / "balanced_joint_cm_reference_node_sets_after_edge_threshold.csv")
    nodes = ref_nodes["node"].astype(str).tolist()
    prefixes = sorted({n.split("_")[0] for n in nodes})
    palette = plt.get_cmap("tab20").colors
    prefix_colors = {p: palette[i % len(palette)] for i, p in enumerate(prefixes)}

    fig, ax = plt.subplots(figsize=(max(4, 0.4 * len(prefixes) + 1), 1.6))
    handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=prefix_colors[p],
                        markersize=12, label=p)
               for p in prefixes]
    ax.legend(handles=handles, loc="center", ncol=min(len(prefixes), 6),
              frameon=False, title="Subtype prefix (lineage)")
    ax.set_axis_off()
    save_pdf_svg(fig, FIG / "nodeplot_network_node_legend")
    print(f"[legend] nodeplot_network_node_legend saved", flush=True)

    # ---- 2. nodeplot_network_edge_colorbar ----
    # Standard context nodeplot edges: r in [EDGE_THRESHOLD, 1.0], colormap viridis
    fig, ax = plt.subplots(figsize=(5, 0.8))
    norm = Normalize(vmin=EDGE_THRESHOLD, vmax=1.0)
    sm = ScalarMappable(norm=norm, cmap="viridis")
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.6, pad=0.05, orientation="horizontal")
    cbar.set_label(f"Edge Pearson r (>= {EDGE_THRESHOLD})")
    ax.set_axis_off()
    save_pdf_svg(fig, FIG / "nodeplot_network_edge_colorbar")
    print(f"[legend] nodeplot_network_edge_colorbar saved", flush=True)

    # ---- 3. tumor_centric_nodeplot_edge_origin_tumor_only_edge_colorbar ----
    # tumor_only edges: red colormap
    fig, ax = plt.subplots(figsize=(5, 0.8))
    sm = ScalarMappable(norm=norm, cmap="Reds")
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.6, pad=0.05, orientation="horizontal")
    cbar.set_label("Tumor-only edge Pearson r (tumor samples)")
    ax.set_axis_off()
    save_pdf_svg(fig, FIG / "tumor_centric_nodeplot_edge_origin_tumor_only_edge_colorbar")
    print(f"[legend] tumor_centric_nodeplot_edge_origin_tumor_only_edge_colorbar saved", flush=True)

    # ---- 4. tumor_centric_nodeplot_edge_origin_shared_edge_colorbar ----
    # shared edges: purple colormap
    fig, ax = plt.subplots(figsize=(5, 0.8))
    sm = ScalarMappable(norm=norm, cmap="Purples")
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.6, pad=0.05, orientation="horizontal")
    cbar.set_label("Shared edge Pearson r (tumor samples; also passes in normal-like)")
    ax.set_axis_off()
    save_pdf_svg(fig, FIG / "tumor_centric_nodeplot_edge_origin_shared_edge_colorbar")
    print(f"[legend] tumor_centric_nodeplot_edge_origin_shared_edge_colorbar saved", flush=True)

    # ---- 5. tumor_centric_nodeplot_edge_origin_edge_class_legend ----
    fig, ax = plt.subplots(figsize=(5, 1.0))
    handles = [
        Line2D([0], [0], color=plt.cm.Reds(0.7), lw=2.5, label="tumor_only (tumor pass, normal-like fail)"),
        Line2D([0], [0], color=plt.cm.Purples(0.7), lw=2.5, label="shared (both contexts pass)"),
    ]
    ax.legend(handles=handles, loc="center", ncol=1, frameon=False)
    ax.set_axis_off()
    save_pdf_svg(fig, FIG / "tumor_centric_nodeplot_edge_origin_edge_class_legend")
    print(f"[legend] tumor_centric_nodeplot_edge_origin_edge_class_legend saved", flush=True)

    print(f"[legend] all 5 legends/colorbars done", flush=True)


if __name__ == "__main__":
    main()