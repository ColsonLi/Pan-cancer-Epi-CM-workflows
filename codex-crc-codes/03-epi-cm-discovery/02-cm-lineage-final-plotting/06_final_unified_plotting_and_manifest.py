#!/usr/bin/env python3
"""Unified canonical plotting for CRC Epi-CM discovery.

Reads only canonical tables from submodule 01 and writes PDF/SVG figures plus
one figure manifest. Statistical results are never recomputed here.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import random
import re
import sys
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D


SEED = 42
random.seed(SEED)
np.random.seed(SEED)

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.family"] = "DejaVu Sans"
mpl.rcParams["axes.linewidth"] = 0.8
mpl.rcParams["figure.facecolor"] = "white"
sns.set_style("ticks")

WORKFLOW_ROOT = Path(__file__).resolve().parents[3]
TABLE_ROOT = WORKFLOW_ROOT / "tables/03-epi-cm-discovery"
ANALYSIS_ROOT = TABLE_ROOT / "01-cm-lineage-analysis"
PREP_DIR = ANALYSIS_ROOT / "01_prepare_inputs_and_frequency_tables"
NMF_DIR = ANALYSIS_ROOT / "02_balanced_joint_nmf"
NODE_DIR = ANALYSIS_ROOT / "03_cm_classification_nodes_and_edges"
SPEARMAN_DIR = ANALYSIS_ROOT / "04_epi_cm_association_spearman"
PEARSON_DIR = ANALYSIS_ROOT / "05_epi_cm_association_pearson"
FIG_ROOT = WORKFLOW_ROOT / "figures/03-epi-cm-discovery/02-cm-lineage-final-plotting"
PLOT_TABLE_DIR = TABLE_ROOT / "02-cm-lineage-final-plotting"
EPI_H5AD = (
    WORKFLOW_ROOT
    / "h5ad/02-cell_subtype_integration_clustering/05-subtype-deg-annotation/epithelial/adata_epi.h5ad"
)

MANIFEST: list[dict] = []


def cm_sort_key(cm: str) -> tuple[int, str]:
    match = re.search(r"_CM(\d+)$", str(cm))
    return (int(match.group(1)) if match else 10**9, str(cm))


def safe_filename_part(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")


def save_pdf_svg(
    fig,
    stem: Path,
    family: str,
    function: str,
    inputs: list[Path] | tuple[Path, ...],
    method: str = "",
    notes: str = "",
) -> tuple[Path, Path]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    pdf = stem.with_suffix(".pdf")
    svg = stem.with_suffix(".svg")
    fig.savefig(pdf, bbox_inches="tight", dpi=300)
    fig.savefig(svg, bbox_inches="tight", dpi=300)
    plt.close(fig)
    for path in (pdf, svg):
        MANIFEST.append(
            {
                "figure_file": str(path),
                "figure_family": family,
                "plotting_function": function,
                "direct_input_tables": ";".join(str(x) for x in inputs),
                "method": method,
                "output_format": path.suffix.lstrip("."),
                "notes": notes,
            }
        )
    return pdf, svg


def read_matrix(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path, index_col=0)
    if df.empty:
        raise ValueError(f"Empty matrix: {path}")
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    return df


def q_stars(value: float) -> str:
    if not np.isfinite(value):
        return "ns"
    if value < 0.001:
        return "***"
    if value < 0.01:
        return "**"
    if value < 0.05:
        return "*"
    return "ns"


def plot_joint_nmf_k_selection() -> None:
    path = NMF_DIR / "joint_nmf_k_selection_metrics.csv"
    metrics = pd.read_csv(path).sort_values("k")
    k = metrics["k"].astype(int).to_numpy()
    selected_k = int(metrics.loc[metrics["selected"].astype(bool), "k"].iloc[0])
    fig, axes = plt.subplots(1, 3, figsize=(10, 3), constrained_layout=True)
    panels = [
        ("best_balanced_explained_fraction", "Balanced fit", "Explained fraction", "#1f77b4"),
        ("stability_matched_cosine", "Stability", "Matched cosine", "#2ca02c"),
        ("selection_score", "Selection score", "Score", "#d62728"),
    ]
    for ax, (column, title, ylabel, color) in zip(axes, panels):
        ax.plot(k, metrics[column].astype(float), marker="o", ms=3.8, lw=1.8, color=color)
        ax.axvline(selected_k, color="black", ls="--", lw=1)
        ax.set_title(title, fontweight="bold")
        ax.set_xlabel("K")
        ax.set_ylabel(ylabel)
        ax.set_xticks(k)
        ax.tick_params(axis="x", labelsize=7)
        ax.grid(axis="y", color="#d9d9d9", lw=0.6, alpha=0.8)
        ax.spines[["top", "right"]].set_visible(False)
    save_pdf_svg(
        fig,
        FIG_ROOT / "nmf-k-selection/joint_nmf_k_selection",
        "nmf_k_selection",
        "plot_joint_nmf_k_selection",
        [path],
        notes=f"selected_K={selected_k}",
    )


def plot_activity_mean_sd_barplot() -> None:
    path = NODE_DIR / "joint_CM_activity_tumor_vs_normal_mean_sd_summary.csv"
    df = pd.read_csv(path)
    cms = df["CM"].astype(str).tolist()
    statuses = ["normal", "tumor"]
    palette = {"normal": "#377EB8", "tumor": "#E41A1C"}
    x = np.arange(len(cms))
    width = 0.36
    offsets = {"normal": -width / 2, "tumor": width / 2}
    fig, ax = plt.subplots(figsize=(max(7, len(cms) * 0.72), 3.8))
    for status in statuses:
        mean = df[f"{status}_mean"].astype(float).to_numpy()
        sd = df[f"{status}_sd"].astype(float).fillna(0).to_numpy()
        ax.bar(x + offsets[status], mean, width=width, color=palette[status], label=status, edgecolor="white")
        ax.errorbar(x + offsets[status], mean, yerr=sd, fmt="none", ecolor="black", lw=0.8, capsize=2)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(cms, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Mean CM activity +/- SD")
    ax.set_title("Joint CM activity in tumor vs normal samples", fontweight="bold")
    ax.grid(axis="y", color="#d9d9d9", lw=0.6, alpha=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False)
    fig.tight_layout()
    save_pdf_svg(
        fig,
        FIG_ROOT / "cm-activity/activity_df_tumor_vs_normal_mean_sd_barplot",
        "cm_activity_mean_sd",
        "plot_activity_mean_sd_barplot",
        [path],
    )


def status_series_colors(sample_meta: pd.DataFrame) -> tuple[pd.DataFrame, dict, dict]:
    status_palette = {"normal": "#377EB8", "tumor": "#E41A1C"}
    series_values = sorted(sample_meta["series"].astype(str).unique())
    colors = sns.color_palette("husl", n_colors=len(series_values)).as_hex()
    series_palette = dict(zip(series_values, colors))
    row_colors = pd.DataFrame(
        {
            "Series": sample_meta["series"].astype(str).map(series_palette),
            "Status": sample_meta["status"].astype(str).map(status_palette),
        },
        index=sample_meta.index,
    )
    return row_colors, series_palette, status_palette


def plot_clustermap(
    df: pd.DataFrame,
    stem: Path,
    title: str,
    cbar_label: str,
    cmap: str,
    center: float | None,
    vmin: float | None,
    vmax: float | None,
    family: str,
    inputs: list[Path],
    row_colors: pd.DataFrame | None = None,
    hide_y: bool = False,
    col_cluster: bool = False,
    row_cluster: bool = True,
) -> None:
    figsize = (max(6.2, df.shape[1] * 0.38 + 2.5), max(4.0, min(18.0, df.shape[0] * 0.08 + 2.8)))
    grid = sns.clustermap(
        df,
        cmap=cmap,
        center=center,
        vmin=vmin,
        vmax=vmax,
        cbar_kws={"label": cbar_label},
        col_cluster=col_cluster,
        row_cluster=row_cluster,
        row_colors=row_colors,
        colors_ratio=(0.06, 0.01),
        figsize=figsize,
        linewidths=0,
    )
    grid.fig.suptitle(title, y=1.02)
    if hide_y:
        grid.ax_heatmap.set_yticks([])
        grid.ax_heatmap.set_yticklabels([])
        grid.ax_heatmap.tick_params(left=False)
    grid.ax_heatmap.tick_params(axis="x", labelrotation=90)
    save_pdf_svg(grid.fig, stem, family, "plot_clustermap", inputs)


def plot_activity_heatmaps() -> None:
    status_path = PREP_DIR / "sample_status.csv"
    sample_meta = pd.read_csv(status_path, index_col=0)
    sample_meta.index = sample_meta.index.astype(str)
    row_colors, series_palette, status_palette = status_series_colors(sample_meta)
    methods = ("raw", "zscore", "robust", "standard_scale_col")
    for method in methods:
        path = NMF_DIR / f"w_df_activity_sample_by_CM_{method}.csv"
        df = read_matrix(path)
        df = df.loc[:, sorted(df.columns, key=cm_sort_key)]
        rc = row_colors.reindex(df.index)
        if method in {"zscore", "robust"}:
            finite = np.abs(df.to_numpy(float)[np.isfinite(df.to_numpy(float))])
            vmax = float(np.percentile(finite, 98)) if finite.size else 1.0
            cmap, center, vmin, vmax_value = "vlag", 0, -vmax, vmax
        elif method == "standard_scale_col":
            cmap, center, vmin, vmax_value = "viridis", None, 0, 1
        else:
            cmap, center, vmin, vmax_value = "viridis", None, 0, None
        plot_clustermap(
            df,
            FIG_ROOT / f"cm-activity/w_df_activity_sample_activity_per_CM_{method}_clustermap",
            f"W_df activity: sample activity per CM ({method})",
            "CM activity",
            cmap,
            center,
            vmin,
            vmax_value,
            "cm_activity_heatmap",
            [path, status_path],
            row_colors=rc,
            hide_y=True,
        )

    minmax_path = NMF_DIR / "w_df_activity_sample_by_CM_standard_scale_col.csv"
    df = read_matrix(minmax_path)
    df = df.loc[:, sorted(df.columns, key=cm_sort_key)]
    rc = row_colors.reindex(df.index)
    plot_clustermap(
        df,
        FIG_ROOT / "cm-activity/w_df_activity_sample_activity_per_CM_standard_scale_col_with_series_clustermap",
        "W_df activity: sample activity per CM (column min-max, with series)",
        "Column min-max activity",
        "viridis",
        None,
        0,
        1,
        "cm_activity_heatmap_annotated",
        [minmax_path, status_path],
        row_colors=rc,
        hide_y=True,
    )
    legend_rows = []
    for label, color in series_palette.items():
        legend_rows.append({"track": "Series", "label": label, "color": color})
    for label, color in status_palette.items():
        legend_rows.append({"track": "Status", "label": label, "color": color})
    pd.DataFrame(legend_rows).to_csv(PLOT_TABLE_DIR / "activity_heatmap_annotation_palette.csv", index=False)


def plot_loading_heatmaps() -> None:
    methods = ("raw", "zscore", "robust", "standard_scale_col")
    for method in methods:
        path = NMF_DIR / f"h_df_loading_cell_subtype_by_CM_{method}.csv"
        # Historical h_df orientation is CM x subtype; final plotting is subtype x CM.
        df = read_matrix(path).T
        df = df.loc[:, sorted(df.columns, key=cm_sort_key)]
        if method in {"zscore", "robust"}:
            finite = np.abs(df.to_numpy(float)[np.isfinite(df.to_numpy(float))])
            vmax = float(np.percentile(finite, 98)) if finite.size else 1.0
            cmap, center, vmin, vmax_value = "vlag", 0, -vmax, vmax
        elif method == "standard_scale_col":
            cmap, center, vmin, vmax_value = "viridis", None, 0, 1
        else:
            cmap, center, vmin, vmax_value = "viridis", None, 0, None
        plot_clustermap(
            df,
            FIG_ROOT / f"cm-loading/h_df_loading_cell_subtype_weights_per_CM_{method}_clustermap",
            f"H_df loading: cell subtype weights per CM ({method})",
            "CM loading",
            cmap,
            center,
            vmin,
            vmax_value,
            "cm_loading_heatmap",
            [path],
            hide_y=False,
        )


def plot_joint_module_top_subtype_heatmap() -> None:
    h_path = NMF_DIR / "H_df.csv"
    node_path = NODE_DIR / "joint_cm_cell_subtype_nodes_top20_from_H_df.csv"
    H = read_matrix(h_path)
    nodes = pd.read_csv(node_path)
    selected = nodes.loc[nodes["rank"].le(12), "cell_subtype"].drop_duplicates().tolist()
    frac = H.div(H.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    frac = frac.loc[sorted(frac.index, key=cm_sort_key), [x for x in selected if x in frac.columns]]
    fig, ax = plt.subplots(figsize=(max(9, frac.shape[1] * 0.18), max(4, frac.shape[0] * 0.34)))
    sns.heatmap(frac, cmap="Reds", linewidths=0.2, linecolor="white", cbar_kws={"label": "Loading fraction"}, ax=ax)
    ax.set_xlabel("Cell subtype")
    ax.set_ylabel("CM")
    ax.set_title("Top subtype loadings", fontweight="bold")
    ax.tick_params(axis="x", rotation=90)
    fig.tight_layout()
    save_pdf_svg(
        fig,
        FIG_ROOT / "cm-loading/joint_module_top_subtype_heatmap",
        "top_subtype_heatmap",
        "plot_joint_module_top_subtype_heatmap",
        [h_path, node_path],
    )


def closest_square_grid(n: int) -> tuple[int, int]:
    if n <= 0:
        return 0, 0
    best = None
    for rows in range(1, int(np.ceil(np.sqrt(n))) + 1):
        cols = int(np.ceil(n / rows))
        candidate = (abs(cols - rows), rows * cols - n, rows, cols)
        if best is None or candidate < best:
            best = candidate
    return best[2], best[3]


def circular_layout(nodes: list[str]) -> dict[str, tuple[float, float]]:
    return {
        node: (float(np.cos(2 * np.pi * i / max(len(nodes), 1))), float(np.sin(2 * np.pi * i / max(len(nodes), 1))))
        for i, node in enumerate(nodes)
    }


def draw_one_cm(
    ax,
    cm: str,
    nodes: list[str],
    edges: pd.DataFrame,
    mode: str,
    prefix_colors: dict[str, str],
    normal_cmap,
    tumor_cmap,
    edge_norm: Normalize,
) -> None:
    graph = nx.Graph()
    graph.add_nodes_from(nodes)
    pos = circular_layout(nodes)
    normal_pass = {
        tuple(sorted((str(r.node_a), str(r.node_b))))
        for r in edges.loc[edges["context"].eq("normal") & edges["edge_pass_r_ge_0.25"].astype(bool)].itertuples()
    }
    if mode == "normal":
        plot_edges = edges.loc[edges["context"].eq("normal") & edges["edge_pass_r_ge_0.25"].astype(bool)]
    else:
        plot_edges = edges.loc[edges["context"].eq("tumor") & edges["edge_pass_r_ge_0.25"].astype(bool)]
    for row in plot_edges.itertuples():
        a, b = str(row.node_a), str(row.node_b)
        if a not in pos or b not in pos:
            continue
        edge_class = mode
        if mode == "tumor_centric":
            edge_class = "shared" if tuple(sorted((a, b))) in normal_pass else "tumor_only"
        graph.add_edge(a, b, pearson_r=float(row.pearson_r), edge_class=edge_class)

    classes = ["normal"] if mode == "normal" else (["tumor"] if mode == "tumor" else ["tumor_only", "shared"])
    for edge_class in classes:
        selected = [(a, b, d) for a, b, d in graph.edges(data=True) if d["edge_class"] == edge_class]
        if not selected:
            continue
        cmap = normal_cmap if edge_class in {"normal", "shared"} else tumor_cmap
        colors = [cmap(edge_norm(d["pearson_r"])) for _, _, d in selected]
        widths = 2.6 if mode == "tumor_centric" else [2.0 + 3.0 * edge_norm(d["pearson_r"]) for _, _, d in selected]
        nx.draw_networkx_edges(graph, pos, edgelist=[(a, b) for a, b, _ in selected], edge_color=colors, width=widths, alpha=0.88, ax=ax)
    nx.draw_networkx_nodes(
        graph,
        pos,
        nodelist=nodes,
        node_color=[prefix_colors.get(n.split("_", 1)[0], "#999999") for n in nodes],
        node_size=1500,
        linewidths=0.8,
        edgecolors="white",
        ax=ax,
    )
    nx.draw_networkx_labels(graph, pos, labels={n: n for n in nodes}, font_size=7, ax=ax)
    ax.set_title(cm, fontsize=11, fontweight="bold")
    ax.set_aspect("equal")
    ax.axis("off")
    return {str(data["edge_class"]) for _, _, data in graph.edges(data=True)}


def add_edge_colorbars(
    fig,
    mode: str,
    edge_norm: Normalize,
    normal_cmap,
    tumor_cmap,
    present_edge_classes: set[str],
) -> None:
    """Embed only the edge-color scales actually used by a saved nodeplot."""
    if mode == "normal" and "normal" in present_edge_classes:
        specs = [(normal_cmap, "Normal edge Pearson r")]
    elif mode == "tumor" and "tumor" in present_edge_classes:
        specs = [(tumor_cmap, "Tumor edge Pearson r")]
    elif mode == "tumor_centric":
        specs = []
        if "tumor_only" in present_edge_classes:
            specs.append((tumor_cmap, "Tumor-only edge Pearson r"))
        if "shared" in present_edge_classes:
            specs.append((normal_cmap, "Shared edge Pearson r"))
    elif mode not in {"normal", "tumor", "tumor_centric"}:
        raise ValueError(f"Unsupported nodeplot mode: {mode}")
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
        sm = ScalarMappable(norm=edge_norm, cmap=cmap)
        sm.set_array([])
        cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
        cb.set_label(label, fontsize=8)
        cb.ax.tick_params(labelsize=7, length=2)


def nodeplot_legends(prefix_colors: dict[str, str], normal_cmap, tumor_cmap, edge_norm: Normalize) -> None:
    out = FIG_ROOT / "nodeplots"
    fig, ax = plt.subplots(figsize=(max(4, len(prefix_colors) * 0.8), 1.4))
    ax.axis("off")
    handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=c, label=p, markersize=8) for p, c in prefix_colors.items()]
    ax.legend(handles=handles, loc="center", ncol=min(5, len(handles)), frameon=False, title="Cell lineage prefix")
    save_pdf_svg(fig, out / "nodeplot_network_node_legend", "nodeplot_legend", "nodeplot_legends", [NODE_DIR / "balanced_joint_cm_reference_node_sets_after_edge_threshold.csv"])

    # The standalone standard-nodeplot reference contains both live palettes.
    fig, axes = plt.subplots(2, 1, figsize=(4.4, 1.35))
    for ax, cmap, label in (
        (axes[0], normal_cmap, "Normal edge Pearson r"),
        (axes[1], tumor_cmap, "Tumor edge Pearson r"),
    ):
        sm = ScalarMappable(norm=edge_norm, cmap=cmap)
        sm.set_array([])
        cb = fig.colorbar(sm, cax=ax, orientation="horizontal")
        cb.set_label(label, fontsize=8)
        cb.ax.tick_params(labelsize=7, length=2)
    fig.subplots_adjust(hspace=1.45, left=0.12, right=0.96, top=0.94, bottom=0.13)
    save_pdf_svg(fig, out / "nodeplot_network_edge_colorbar", "nodeplot_legend", "nodeplot_legends", [NODE_DIR / "status_specific_nodeplot_edges.csv"])

    for stem, cmap, label in [
        ("tumor_centric_nodeplot_edge_origin_tumor_only_edge_colorbar", tumor_cmap, "Tumor-only edge correlation in tumor samples (r)"),
        ("tumor_centric_nodeplot_edge_origin_shared_edge_colorbar", normal_cmap, "Shared edge correlation in tumor samples (r)"),
    ]:
        fig, ax = plt.subplots(figsize=(4.4, 0.55))
        sm = ScalarMappable(norm=edge_norm, cmap=cmap)
        sm.set_array([])
        cb = fig.colorbar(sm, cax=ax, orientation="horizontal")
        cb.set_label(label)
        save_pdf_svg(fig, out / stem, "nodeplot_legend", "nodeplot_legends", [NODE_DIR / "status_specific_nodeplot_edges.csv"])
    legend_r = float(edge_norm.vmax)
    fig, ax = plt.subplots(figsize=(3.8, 1.1))
    ax.axis("off")
    handles = [
        Line2D([0], [0], color=tumor_cmap(edge_norm(legend_r)), lw=3, alpha=0.88, label="Tumor only"),
        Line2D([0], [0], color=normal_cmap(edge_norm(legend_r)), lw=3, alpha=0.88, label="Shared with normal-like"),
    ]
    ax.legend(handles=handles, loc="center", frameon=False, ncol=2)
    save_pdf_svg(fig, out / "tumor_centric_nodeplot_edge_origin_edge_class_legend", "nodeplot_legend", "nodeplot_legends", [NODE_DIR / "status_specific_nodeplot_edges.csv"])


def plot_nodeplots() -> None:
    membership_path = NODE_DIR / "balanced_joint_cm_reference_node_sets_after_edge_threshold.csv"
    edge_path = NODE_DIR / "status_specific_nodeplot_edges.csv"
    nodes = pd.read_csv(membership_path).sort_values(["CM", "reference_node_rank"])
    edges = pd.read_csv(edge_path)
    cms = sorted(nodes["CM"].astype(str).unique(), key=cm_sort_key)
    prefixes = sorted(nodes["node"].astype(str).str.split("_", n=1).str[0].unique())
    colors = sns.color_palette("tab20", n_colors=max(3, len(prefixes))).as_hex()
    prefix_colors = dict(zip(prefixes, colors))
    normal_cmap = LinearSegmentedColormap.from_list("normal_shared", ["#fee2e2", "#b91c1c"])
    tumor_cmap = LinearSegmentedColormap.from_list("tumor_only", ["#dbeafe", "#1d4ed8"])
    edge_norm = Normalize(vmin=0.25, vmax=1.0)
    modes = {
        "normal": ("normal_like_nodeplots_by_cm", "_normal_like_nodeplot", "normal_like_all_CM_nodeplot"),
        "tumor": ("tumor_nodeplots_by_cm", "_tumor_nodeplot", "tumor_all_CM_nodeplot"),
        "tumor_centric": ("tumor_centric_nodeplots_by_cm", "_tumor_centric_nodeplot_edge_origin", "tumor_centric_nodeplot_edge_origin"),
    }
    for cm in cms:
        cm_nodes = nodes.loc[nodes["CM"].eq(cm), "node"].astype(str).tolist()
        cm_edges = edges.loc[edges["CM"].eq(cm)].copy()
        for mode, (subdir, suffix, _) in modes.items():
            canvas = (6.8, 4.8) if mode == "tumor_centric" else (4.8, 4.8)
            fig, ax = plt.subplots(figsize=canvas)
            present_edge_classes = draw_one_cm(
                ax, cm, cm_nodes, cm_edges, mode,
                prefix_colors, normal_cmap, tumor_cmap, edge_norm,
            )
            fig.tight_layout(rect=[0, 0.16, 1, 1])
            add_edge_colorbars(
                fig,
                mode=mode,
                edge_norm=edge_norm,
                normal_cmap=normal_cmap,
                tumor_cmap=tumor_cmap,
                present_edge_classes=present_edge_classes,
            )
            save_pdf_svg(
                fig,
                FIG_ROOT / "nodeplots" / subdir / f"{cm}{suffix}",
                "cm_nodeplot",
                "plot_nodeplots",
                [membership_path, edge_path],
                notes=f"context={mode}; status values in table are normal/tumor",
            )
    rows, cols = closest_square_grid(len(cms))
    for mode, (_, _, overview_stem) in modes.items():
        fig, axes = plt.subplots(rows, cols, figsize=(cols * 4.4, rows * 4.2 + 0.9), squeeze=False)
        overview_edge_classes = set()
        for ax, cm in zip(axes.ravel(), cms):
            cm_nodes = nodes.loc[nodes["CM"].eq(cm), "node"].astype(str).tolist()
            cm_edges = edges.loc[edges["CM"].eq(cm)].copy()
            panel_edge_classes = draw_one_cm(
                ax, cm, cm_nodes, cm_edges, mode,
                prefix_colors, normal_cmap, tumor_cmap, edge_norm,
            )
            overview_edge_classes.update(panel_edge_classes)
        for ax in axes.ravel()[len(cms):]:
            ax.axis("off")
        fig.tight_layout(rect=[0, 0.11, 1, 1])
        add_edge_colorbars(
            fig,
            mode=mode,
            edge_norm=edge_norm,
            normal_cmap=normal_cmap,
            tumor_cmap=tumor_cmap,
            present_edge_classes=overview_edge_classes,
        )
        save_pdf_svg(
            fig,
            FIG_ROOT / "nodeplots" / overview_stem,
            "cm_nodeplot_overview",
            "plot_nodeplots",
            [membership_path, edge_path],
            notes=f"context={mode}; closest-to-square grid",
        )
    nodeplot_legends(prefix_colors, normal_cmap, tumor_cmap, edge_norm)
    pd.DataFrame([{"prefix": k, "color": v} for k, v in prefix_colors.items()]).to_csv(PLOT_TABLE_DIR / "node_prefix_palette.csv", index=False)


def top10_nodes_by_cm(top_nodes: pd.DataFrame, cm: str) -> list[str]:
    """Return exactly the H/loading-ranked diagnostic top10 nodes for one CM."""
    rank_col = "rank" if "rank" in top_nodes.columns else "loading_rank"
    sub = top_nodes.loc[top_nodes["CM"].astype(str).eq(str(cm))].copy()
    sub = sub.sort_values(rank_col, kind="stable").head(10)
    node_col = "cell_subtype" if "cell_subtype" in sub.columns else "node"
    nodes = sub[node_col].astype(str).tolist()
    if len(nodes) != 10:
        raise ValueError(f"{cm}: expected exactly 10 diagnostic top nodes, found {len(nodes)}")
    return nodes


def plot_top10_node_correlation_heatmap(
    corr: pd.DataFrame,
    top_nodes: pd.DataFrame,
    cm: str,
    stem: Path,
    title: str,
    inputs: list[Path],
) -> None:
    nodes = top10_nodes_by_cm(top_nodes, cm)
    nodes = [node for node in nodes if node in corr.index and node in corr.columns]
    if len(nodes) < 2:
        raise ValueError(f"{cm}: fewer than 2 top10 nodes present in correlation matrix")
    plot_df = corr.loc[nodes, nodes].astype(float)
    fig, ax = plt.subplots(figsize=(max(4.2, 0.48 * len(nodes)), max(4.0, 0.48 * len(nodes))))
    sns.heatmap(
        plot_df,
        cmap="RdBu_r",
        center=0,
        vmin=-1,
        vmax=1,
        square=True,
        linewidths=0.2,
        linecolor="white",
        annot=False,
        cbar_kws={"label": "Pearson r"},
        ax=ax,
    )
    ax.set_title(title)
    ax.set_xticklabels(nodes, rotation=90, ha="center", fontsize=7)
    ax.set_yticklabels(nodes, rotation=0, fontsize=7)
    ax.tick_params(axis="both", length=0)
    fig.tight_layout()
    save_pdf_svg(
        fig,
        stem,
        "top10_node_correlation_heatmap_by_cm",
        "plot_top10_node_correlation_heatmap",
        inputs,
    )


def plot_top10_node_correlation_summary(
    corr: pd.DataFrame,
    top_nodes: pd.DataFrame,
    cm_order: list[str],
    stem: Path,
    context_label: str,
    inputs: list[Path],
    ncols: int = 4,
) -> None:
    """Draw one separate 10x10 top-node matrix per CM in a shared summary."""
    panel_data = []
    for cm in cm_order:
        nodes = top10_nodes_by_cm(top_nodes, cm)
        nodes = [node for node in nodes if node in corr.index and node in corr.columns]
        if len(nodes) < 2:
            raise ValueError(f"{cm}: fewer than 2 top10 nodes present in correlation matrix")
        panel_data.append((cm, nodes, corr.loc[nodes, nodes].astype(float)))

    if not panel_data:
        raise ValueError(f"No CM panels available for {context_label}")

    ncols = min(ncols, len(panel_data))
    nrows = int(np.ceil(len(panel_data) / ncols))
    panel_size = 4.0
    fig = plt.figure(figsize=(panel_size * ncols + 0.8, panel_size * nrows))
    grid = fig.add_gridspec(
        nrows,
        ncols + 1,
        width_ratios=[1.0] * ncols + [0.05],
        wspace=0.65,
        hspace=0.75,
    )
    cbar_ax = fig.add_subplot(grid[:, -1])

    for panel_i, (cm, nodes, plot_df) in enumerate(panel_data):
        row, col = divmod(panel_i, ncols)
        ax = fig.add_subplot(grid[row, col])
        sns.heatmap(
            plot_df,
            cmap="RdBu_r",
            center=0,
            vmin=-1,
            vmax=1,
            square=True,
            linewidths=0.2,
            linecolor="white",
            annot=False,
            cbar=(panel_i == 0),
            cbar_ax=cbar_ax if panel_i == 0 else None,
            cbar_kws={"label": "Pearson r"} if panel_i == 0 else None,
            ax=ax,
        )
        ax.set_title(str(cm))
        ax.set_xticklabels(nodes, rotation=90, ha="center", fontsize=6.5)
        ax.set_yticklabels(nodes, rotation=0, fontsize=6.5)
        ax.tick_params(axis="both", length=0)

    for empty_i in range(len(panel_data), nrows * ncols):
        row, col = divmod(empty_i, ncols)
        ax = fig.add_subplot(grid[row, col])
        ax.set_axis_off()

    fig.suptitle(f"{context_label}: top10 node correlations by CM", y=0.995)
    fig.subplots_adjust(top=0.94, bottom=0.08, left=0.06, right=0.95)
    save_pdf_svg(
        fig,
        stem,
        "top10_node_correlation_summary",
        "plot_top10_node_correlation_summary",
        inputs,
        notes="15 independent 10x10 CM panels; no cross-CM union matrix",
    )


def plot_top10_correlations() -> None:
    top_path = NODE_DIR / "joint_cm_cell_subtype_nodes_top10_from_H_df.csv"
    classification_path = NODE_DIR / "joint_module_classification.csv"
    top = pd.read_csv(top_path)
    classification = pd.read_csv(classification_path).sort_values("global_order", kind="stable")
    cm_order = classification["CM"].astype(str).tolist()
    if set(cm_order) != set(top["CM"].astype(str)):
        raise ValueError("Top10 CM IDs do not match joint_module_classification.csv")
    normal_path = NODE_DIR / "normal_node_node_correlation_matrix.csv"
    tumor_path = NODE_DIR / "tumor_node_node_correlation_matrix.csv"
    normal_corr = read_matrix(normal_path)
    tumor_corr = read_matrix(tumor_path)
    common_inputs = [top_path, classification_path]
    out_root = FIG_ROOT / "node-correlation-heatmaps"

    plot_top10_node_correlation_summary(
        normal_corr,
        top,
        cm_order,
        FIG_ROOT / "node-correlation-heatmaps/normal_like_top10_node_correlation_heatmap_no_edge_filter",
        "Normal",
        [normal_path, *common_inputs],
    )
    plot_top10_node_correlation_summary(
        tumor_corr,
        top,
        cm_order,
        FIG_ROOT / "node-correlation-heatmaps/tumor_top10_node_correlation_heatmap_no_edge_filter",
        "Tumor",
        [tumor_path, *common_inputs],
    )

    contexts = {
        "normal_like": (normal_corr, normal_path, "Normal"),
        "tumor": (tumor_corr, tumor_path, "Tumor"),
    }
    for context, (corr, corr_path, context_label) in contexts.items():
        for cm in cm_order:
            plot_top10_node_correlation_heatmap(
                corr,
                top,
                cm,
                out_root
                / "top10_node_correlation_heatmaps_no_edge_filter_by_cm"
                / context
                / f"{cm}_{context}_top10_node_correlation_heatmap_no_edge_filter",
                f"{cm}: top10 H-loading node correlations ({context_label})",
                [corr_path, *common_inputs],
            )


def plot_qstar_heatmap(value: pd.DataFrame, q: pd.DataFrame, stem: Path, title: str, label: str, inputs: list[Path], method: str) -> None:
    value = value.loc[:, sorted(value.columns, key=cm_sort_key)]
    q = q.reindex(index=value.index, columns=value.columns)
    annot = q.apply(lambda col: col.map(q_stars))
    vmax = float(np.nanmax(np.abs(value.to_numpy(float))))
    if not np.isfinite(vmax) or vmax == 0:
        vmax = 1.0
    fig, ax = plt.subplots(figsize=(max(9, value.shape[1] * 0.55 + 2), max(5, value.shape[0] * 0.42 + 1.8)))
    sns.heatmap(
        value.astype(float),
        cmap="coolwarm",
        center=0,
        vmin=-vmax,
        vmax=vmax,
        annot=annot,
        fmt="",
        annot_kws={"fontsize": 7},
        linewidths=0.25,
        linecolor="white",
        cbar_kws={"label": label},
        ax=ax,
    )
    ax.set_title(title)
    ax.set_xlabel("CM")
    ax.set_ylabel("Epithelial subtype")
    ax.tick_params(axis="x", labelrotation=90)
    ax.tick_params(axis="y", labelrotation=0)
    fig.tight_layout()
    save_pdf_svg(fig, stem, "epi_cm_qstar_heatmap", "plot_qstar_heatmap", inputs, method=method)


def plot_epi_cm_heatmaps(method: str) -> None:
    if method == "spearman":
        root = SPEARMAN_DIR
        stat = "rho"
        normal_value = root / "balanced_joint_cm_epi_cm_association_normal-like_rho_matrix.csv"
        normal_q = root / "balanced_joint_cm_epi_cm_association_normal-like_q_matrix.csv"
        tumor_value = root / "balanced_joint_cm_epi_cm_association_tumor_rho_matrix.csv"
        tumor_q = root / "balanced_joint_cm_epi_cm_association_tumor_q_matrix.csv"
    else:
        root = PEARSON_DIR
        stat = "r"
        normal_value = root / "balanced_joint_cm_epi_cm_association_normal-like_pearson_r_matrix.csv"
        normal_q = root / "balanced_joint_cm_epi_cm_association_normal-like_pearson_q_matrix.csv"
        tumor_value = root / "balanced_joint_cm_epi_cm_association_tumor_pearson_r_matrix.csv"
        tumor_q = root / "balanced_joint_cm_epi_cm_association_tumor_pearson_q_matrix.csv"
    out = FIG_ROOT / f"epi-cm-heatmaps-{method}"
    plot_qstar_heatmap(
        read_matrix(tumor_value), read_matrix(tumor_q),
        out / "balanced_joint_cm_epi_cm_association_tumor_heatmap_qstars",
        f"Balanced joint CM-Epi association (tumor; {method})", f"{method.title()} {stat}",
        [tumor_value, tumor_q], method,
    )
    plot_qstar_heatmap(
        read_matrix(normal_value), read_matrix(normal_q),
        out / "balanced_joint_cm_epi_cm_association_normal-like_heatmap_qstars",
        f"Balanced joint CM-Epi association (normal; {method})", f"{method.title()} {stat}",
        [normal_value, normal_q], method,
    )


def epithelial_palette() -> dict[str, str]:
    if not EPI_H5AD.exists():
        raise FileNotFoundError(EPI_H5AD)
    adata = ad.read_h5ad(EPI_H5AD, backed="r")
    labels = adata.obs["cell_subtype"]
    if not hasattr(labels, "cat"):
        raise ValueError("Epithelial cell_subtype is not categorical; cannot align palette")
    categories = labels.cat.categories.astype(str).tolist()
    colors = [str(x) for x in adata.uns.get("cell_subtype_colors", [])]
    adata.file.close()
    if len(categories) != len(colors):
        raise ValueError(f"Epithelial subtype palette mismatch: {len(categories)} labels, {len(colors)} colors")
    palette = dict(zip(categories, colors))
    pd.DataFrame([{"cell_subtype": k, "color": v} for k, v in palette.items()]).to_csv(
        PLOT_TABLE_DIR / "epithelial_cell_subtype_colors.csv", index=False
    )
    return palette


def fmt_num(value: float) -> str:
    return "NA" if not np.isfinite(value) else f"{value:.3f}"


def fmt_q(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    return f"{value:.1e}" if value < 0.001 else f"{value:.3f}"


def plot_all_scatter(method: str, palette: dict[str, str]) -> None:
    epi_path = PREP_DIR / "epi_subtype_frequency.csv"
    activity_path = NMF_DIR / "W_df.csv"
    status_path = PREP_DIR / "sample_status.csv"
    summary_path = (SPEARMAN_DIR if method == "spearman" else PEARSON_DIR) / f"epi_cm_{method}_all_pairs_long.csv"
    E = read_matrix(epi_path)
    C = read_matrix(activity_path)
    status = pd.read_csv(status_path, index_col=0)["status"].astype(str)
    status.index = status.index.astype(str)
    summary = pd.read_csv(summary_path)
    stat_col = "rho" if method == "spearman" else "r"
    expected_pairs = len(E.columns) * len(C.columns) * 2
    if len(summary) != expected_pairs:
        raise ValueError(f"{method}: expected {expected_pairs} rows, found {len(summary)}")
    missing_colors = sorted(set(E.columns) - set(palette))
    if missing_colors:
        raise ValueError(f"Missing epithelial colors: {missing_colors}")
    figure_rows = []
    for row in summary.itertuples(index=False):
        context = str(row.status)
        cm = str(row.CM)
        epi = str(row.epi_subtype)
        samples = status.index[status.eq(context)].intersection(E.index).intersection(C.index)
        plot_df = pd.DataFrame(
            {
                "CM_score": C.loc[samples, cm].astype(float),
                "Epi_fraction": E.loc[samples, epi].astype(float),
            },
            index=samples,
        ).dropna()
        fig, ax = plt.subplots(figsize=(6, 5))
        sns.regplot(
            data=plot_df,
            x="CM_score",
            y="Epi_fraction",
            color=palette[epi],
            scatter_kws={"s": 60, "alpha": 0.85, "edgecolor": "none"},
            line_kws={"lw": 2},
            ci=95,
            ax=ax,
        )
        ax.set_xlabel(f"{cm} score")
        ax.set_ylabel(f"{epi} fraction")
        ax.set_title(f"{context.title()} samples: {cm} vs {epi}")
        ax.set_ylim(0, 1)
        statistic = float(getattr(row, stat_col))
        ax.text(
            0.02,
            0.98,
            f"{method.title()} {stat_col}={fmt_num(statistic)}, q={fmt_q(float(row.q_value))}\nn={int(row.n_samples)}",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=10,
            bbox={"boxstyle": "round,pad=0.3", "fc": "white", "ec": "gray", "alpha": 0.8},
        )
        sns.despine(ax=ax)
        fig.tight_layout()
        stem = (
            FIG_ROOT
            / f"epi-cm-scatterplots-{method}/{context}"
            / f"scatter_{safe_filename_part(cm)}_vs_{safe_filename_part(epi)}"
        )
        pdf, svg = save_pdf_svg(
            fig,
            stem,
            "epi_cm_scatter_all_pairs",
            "plot_all_scatter",
            [epi_path, activity_path, status_path, summary_path, PLOT_TABLE_DIR / "epithelial_cell_subtype_colors.csv"],
            method=method,
            notes=f"status={context}; y-axis fixed 0..1; seaborn 95% CI",
        )
        for file in (pdf, svg):
            figure_rows.append(
                {
                    "method": method,
                    "status": context,
                    "CM": cm,
                    "epi_subtype": epi,
                    "figure_file": str(file),
                }
            )
    pd.DataFrame(figure_rows).to_csv(PLOT_TABLE_DIR / f"epi_cm_{method}_scatter_all_pairs_manifest.csv", index=False)


def write_reports() -> None:
    manifest = pd.DataFrame(MANIFEST)
    manifest.to_csv(PLOT_TABLE_DIR / "figure_manifest.csv", index=False)
    required_columns = {
        "figure_file", "figure_family", "plotting_function", "direct_input_tables", "method", "output_format", "notes"
    }
    if set(manifest.columns) != required_columns:
        raise ValueError(f"Figure manifest columns mismatch: {manifest.columns.tolist()}")
    files = [Path(x) for x in manifest["figure_file"]]
    missing = [str(x) for x in files if not x.exists() or x.stat().st_size == 0]
    pngs = list((WORKFLOW_ROOT / "figures/03-epi-cm-discovery").rglob("*.png"))
    validation = {
        "figure_count": len(files),
        "pdf_count": int(manifest["output_format"].eq("pdf").sum()),
        "svg_count": int(manifest["output_format"].eq("svg").sum()),
        "missing_or_empty_files": missing,
        "png_files": [str(x) for x in pngs],
        "all_manifest_files_exist": not missing,
        "pdf_svg_only": not pngs and set(manifest["output_format"]) == {"pdf", "svg"},
        "spearman_scatter_panels": int((manifest["figure_family"].eq("epi_cm_scatter_all_pairs") & manifest["method"].eq("spearman") & manifest["output_format"].eq("pdf")).sum()),
        "pearson_scatter_panels": int((manifest["figure_family"].eq("epi_cm_scatter_all_pairs") & manifest["method"].eq("pearson") & manifest["output_format"].eq("pdf")).sum()),
        "status_contract": ["normal", "tumor"],
        "spatial_validation": "skipped_no_spatial_input",
    }
    write_path = PLOT_TABLE_DIR / "validation_report.json"
    write_path.write_text(json.dumps(validation, indent=2, ensure_ascii=False), encoding="utf-8")
    if missing or pngs:
        raise RuntimeError(f"Final figure validation failed: {validation}")
    params = {
        "code_file": str(Path(__file__).resolve()),
        "analysis_table_root": str(ANALYSIS_ROOT),
        "figure_root": str(FIG_ROOT),
        "seed": SEED,
        "methods": ["spearman", "pearson"],
        "status_values": ["normal", "tumor"],
        "normal_like_filename_note": "Canonical filenames retain normal-like stems; plotted/data status is normal",
        "figure_formats": ["pdf", "svg"],
    }
    (PLOT_TABLE_DIR / "run_parameters.json").write_text(json.dumps(params, indent=2), encoding="utf-8")
    versions = (
        f"python={sys.version.split()[0]}\n"
        f"platform={platform.platform()}\n"
        f"environment={sys.executable}\n"
        f"numpy={np.__version__}\n"
        f"pandas={pd.__version__}\n"
        f"scipy={scipy.__version__}\n"
        f"matplotlib={mpl.__version__}\n"
        f"seaborn={sns.__version__}\n"
        f"networkx={nx.__version__}\n"
        f"seed={SEED}\n"
        "backend=CPU; canonical plotting stack has no GPU backend\n"
        f"code={Path(__file__).resolve()}\n"
    )
    (PLOT_TABLE_DIR / "package_versions.txt").write_text(versions, encoding="utf-8")
    (PLOT_TABLE_DIR / "readme.txt").write_text(
        "This plotting submodule reads only canonical tables from 01-cm-lineage-analysis.\n"
        "It does not recompute NMF, classification, correlations, p values, or q values.\n"
        "Figures are PDF/SVG only. All Spearman and Pearson Epi x CM pairs are plotted for normal and tumor.\n"
        "The canonical normal-like filename stem denotes the dataset's status=normal group; no AT/HD grouping is used.\n",
        encoding="utf-8",
    )
    pd.DataFrame(
        [
            {
                "step": "canonical_cm_plotting",
                "planned_backend": "CPU",
                "attempted_backend": "CPU",
                "status": "completed",
                "error_summary": "",
                "fallback_backend": "",
                "clean_input_reloaded": True,
                "final_backend_for_rerun": "CPU",
            }
        ]
    ).to_csv(PLOT_TABLE_DIR / "gpu_backend_capability_summary.csv", index=False)


def write_top10_redraw_reports() -> None:
    """Replace only the canonical top10-correlation rows in the full manifest."""
    redraw_manifest = pd.DataFrame(MANIFEST)
    expected = 64
    if len(redraw_manifest) != expected:
        raise ValueError(f"Expected {expected} redrawn PDF/SVG manifest rows, found {len(redraw_manifest)}")

    manifest_path = PLOT_TABLE_DIR / "figure_manifest.csv"
    previous = pd.read_csv(manifest_path)
    keep = ~previous["figure_file"].astype(str).str.contains("/node-correlation-heatmaps/", regex=False)
    combined = pd.concat([previous.loc[keep], redraw_manifest], ignore_index=True)
    combined.to_csv(manifest_path, index=False)
    redraw_manifest.to_csv(PLOT_TABLE_DIR / "top10_node_correlation_figure_manifest.csv", index=False)

    files = [Path(path) for path in combined["figure_file"].astype(str)]
    missing = [str(path) for path in files if not path.is_file() or path.stat().st_size == 0]
    pngs = list((WORKFLOW_ROOT / "figures/03-epi-cm-discovery").rglob("*.png"))
    validation_path = PLOT_TABLE_DIR / "validation_report.json"
    validation = json.loads(validation_path.read_text(encoding="utf-8"))
    validation.update(
        {
            "figure_count": len(combined),
            "pdf_count": int(combined["output_format"].eq("pdf").sum()),
            "svg_count": int(combined["output_format"].eq("svg").sum()),
            "missing_or_empty_files": missing,
            "png_files": [str(path) for path in pngs],
            "all_manifest_files_exist": not missing,
            "pdf_svg_only": not pngs and set(combined["output_format"]) == {"pdf", "svg"},
            "top10_node_correlation_contract": {
                "summary_layout": "15 independent CM-specific 10x10 panels per status",
                "per_cm_statuses": ["normal_like", "tumor"],
                "cm_count": 15,
                "figure_stems": 32,
                "pdf_svg_files": expected,
            },
        }
    )
    validation_path.write_text(json.dumps(validation, indent=2, ensure_ascii=False), encoding="utf-8")
    summary = {
        "status": "completed",
        "source_tables_reused": [
            str(NODE_DIR / "joint_cm_cell_subtype_nodes_top10_from_H_df.csv"),
            str(NODE_DIR / "joint_module_classification.csv"),
            str(NODE_DIR / "normal_node_node_correlation_matrix.csv"),
            str(NODE_DIR / "tumor_node_node_correlation_matrix.csv"),
        ],
        "statistics_recomputed": False,
        "cm_count": 15,
        "per_cm_matrices_per_status": 15,
        "summary_panels_per_status": 15,
        "figure_stems": 32,
        "pdf_svg_files": expected,
        "normal_status_value": "normal",
        "normal_like_filename_note": "normal_like output stems represent status=normal",
    }
    (PLOT_TABLE_DIR / "top10_node_correlation_redraw_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    if missing or pngs:
        raise RuntimeError(f"Top10 correlation redraw validation failed: missing={missing}, pngs={pngs}")


def write_nodeplot_redraw_reports() -> None:
    """Replace only canonical nodeplot rows after embedded-colorbar redraw."""
    redraw_manifest = pd.DataFrame(MANIFEST)
    expected = 106
    if len(redraw_manifest) != expected:
        raise ValueError(f"Expected {expected} nodeplot PDF/SVG manifest rows, found {len(redraw_manifest)}")

    manifest_path = PLOT_TABLE_DIR / "figure_manifest.csv"
    previous = pd.read_csv(manifest_path)
    keep = ~previous["figure_file"].astype(str).str.contains("/nodeplots/", regex=False)
    combined = pd.concat([previous.loc[keep], redraw_manifest], ignore_index=True)
    combined.to_csv(manifest_path, index=False)
    redraw_manifest.to_csv(PLOT_TABLE_DIR / "nodeplot_embedded_colorbar_figure_manifest.csv", index=False)

    files = [Path(path) for path in combined["figure_file"].astype(str)]
    missing = [str(path) for path in files if not path.is_file() or path.stat().st_size == 0]
    pngs = list((WORKFLOW_ROOT / "figures/03-epi-cm-discovery").rglob("*.png"))
    validation_path = PLOT_TABLE_DIR / "validation_report.json"
    validation = json.loads(validation_path.read_text(encoding="utf-8"))
    validation.update(
        {
            "figure_count": len(combined),
            "pdf_count": int(combined["output_format"].eq("pdf").sum()),
            "svg_count": int(combined["output_format"].eq("svg").sum()),
            "missing_or_empty_files": missing,
            "png_files": [str(path) for path in pngs],
            "all_manifest_files_exist": not missing,
            "pdf_svg_only": not pngs and set(combined["output_format"]) == {"pdf", "svg"},
            "nodeplot_embedded_colorbar_contract": {
                "normal_status_value": "normal",
                "single_cm_figures": 45,
                "all_cm_overviews": 3,
                "standalone_legend_colorbar_stems": 5,
                "figure_stems": 53,
                "pdf_svg_files": expected,
                "colorbar_source": "same live colormap and Normalize(vmin=0.25, vmax=1.0) as plotted edges",
            },
        }
    )
    validation_path.write_text(json.dumps(validation, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    summary = {
        "latest_skill_sha256": "6d2dcca7db984346b8e6042597ee2bba54d46335e63e0ce30ba4cf0483e93d13",
        "status": "completed",
        "statistics_recomputed": False,
        "source_tables_reused": [
            str(NODE_DIR / "balanced_joint_cm_reference_node_sets_after_edge_threshold.csv"),
            str(NODE_DIR / "status_specific_nodeplot_edges.csv"),
        ],
        "cm_count": 15,
        "normal_per_cm_with_embedded_bar": 15,
        "tumor_per_cm_with_embedded_bar": 15,
        "tumor_centric_per_cm_with_embedded_bars": 15,
        "all_cm_overviews_with_embedded_bars": 3,
        "standalone_legend_colorbar_stems": 5,
        "figure_stems": 53,
        "pdf_svg_files": expected,
        "png_files": [],
    }
    (PLOT_TABLE_DIR / "nodeplot_embedded_colorbar_redraw_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    if missing or pngs:
        raise RuntimeError(f"Nodeplot redraw validation failed: missing={missing}, pngs={pngs}")


def parse_cli_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--only-top10-node-correlations",
        action="store_true",
        help="Redraw only the status-specific top10 node-correlation figure family.",
    )
    parser.add_argument(
        "--only-nodeplots",
        action="store_true",
        help="Redraw only nodeplots and their embedded/standalone edge colorbars.",
    )
    return parser.parse_args()


def main() -> None:
    PLOT_TABLE_DIR.mkdir(parents=True, exist_ok=True)
    cli = parse_cli_args()
    if cli.only_nodeplots:
        plot_nodeplots()
        write_nodeplot_redraw_reports()
        print(json.dumps({"status": "completed", "figure_files": len(MANIFEST)}, ensure_ascii=False))
        return
    if cli.only_top10_node_correlations:
        plot_top10_correlations()
        write_top10_redraw_reports()
        print(json.dumps({"status": "completed", "figure_files": len(MANIFEST)}, ensure_ascii=False))
        return
    plot_joint_nmf_k_selection()
    plot_activity_mean_sd_barplot()
    plot_activity_heatmaps()
    plot_loading_heatmaps()
    plot_joint_module_top_subtype_heatmap()
    plot_nodeplots()
    plot_top10_correlations()
    plot_epi_cm_heatmaps("spearman")
    plot_epi_cm_heatmaps("pearson")
    palette = epithelial_palette()
    plot_all_scatter("spearman", palette)
    plot_all_scatter("pearson", palette)
    write_reports()
    print(json.dumps({"status": "completed", "figure_files": len(MANIFEST)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
