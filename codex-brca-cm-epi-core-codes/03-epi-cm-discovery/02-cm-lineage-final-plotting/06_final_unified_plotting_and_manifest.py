#!/usr/bin/env python3
"""Canonical unified final plotting for BRCA tumor-only Epi-CM discovery."""

from __future__ import annotations

import importlib.metadata
import json
import math
import re
import sys
import time
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.ticker import FormatStrFormatter

CODE_DIR = Path(__file__).resolve().parents[1] / "01-cm-lineage-analysis"
sys.path.insert(0, str(CODE_DIR))
from cm_analysis_common import (  # noqa: E402
    EDGE_R_THRESHOLD, FIGURE_ROOT, NMF_DIR, NODE_DIR, PLOT_TABLE_ROOT,
    PREP_DIR, SEED, SPEARMAN_DIR, cm_sort_key, write_json,
)


CODE_PATH = Path(__file__).resolve()
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.family"] = "DejaVu Sans"
mpl.rcParams["figure.facecolor"] = "white"
mpl.rcParams["axes.facecolor"] = "white"
sns.set_style("ticks")

MANIFEST_ROWS: list[dict[str, object]] = []


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def safe_part(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")


def save_pdf_svg(
    fig: mpl.figure.Figure,
    stem: Path,
    *,
    family: str,
    function: str,
    inputs: list[Path],
    method: str,
    notes: str = "",
) -> tuple[Path, Path]:
    stem.parent.mkdir(parents=True, exist_ok=True)
    pdf, svg = stem.with_suffix(".pdf"), stem.with_suffix(".svg")
    fig.savefig(pdf, bbox_inches="tight", dpi=300)
    fig.savefig(svg, bbox_inches="tight", dpi=300)
    plt.close(fig)
    for output in [pdf, svg]:
        MANIFEST_ROWS.append(
            {
                "figure_file": str(output), "figure_family": family,
                "plotting_function": function,
                "direct_input_tables": ";".join(map(str, inputs)),
                "method": method, "output_format": output.suffix.lstrip("."),
                "notes": notes,
            }
        )
    return pdf, svg


def q_star(value: float) -> str:
    if pd.isna(value):
        return ""
    if value < 0.001:
        return "***"
    if value < 0.01:
        return "**"
    if value < 0.05:
        return "*"
    return "ns"


def plot_k_selection(output_root: Path = FIGURE_ROOT) -> None:
    source = NMF_DIR / "joint_nmf_k_selection_metrics.csv"
    frame = pd.read_csv(source).sort_values("k")
    k = frame["k"].astype(int).to_numpy()
    selected = int(frame.loc[frame["selected"].astype(bool), "k"].iloc[0])
    fig, axes = plt.subplots(1, 3, figsize=(10, 3), constrained_layout=True)
    panels = [
        ("best_balanced_explained_fraction", "Balanced fit", "Explained fraction", "#1f77b4", np.arange(0.50, 0.96, 0.10)),
        ("stability_matched_cosine", "Stability", "Matched cosine", "#2ca02c", np.arange(0.50, 1.01, 0.10)),
        ("selection_score", "Selection score", "Score", "#d62728", np.arange(0.50, 0.91, 0.10)),
    ]
    for ax, (column, title, ylabel, color, yticks) in zip(axes, panels):
        ax.plot(k, frame[column].astype(float), marker="o", markersize=3.8, linewidth=1.8, color=color)
        ax.set_title(title, fontweight="bold")
        ax.set_xlabel("k")
        ax.set_ylabel(ylabel)
        ax.set_xticks(k)
        ax.set_xlim(k.min() - 0.5, k.max() + 0.5)
        ax.set_yticks(yticks)
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
        ax.tick_params(axis="x", labelsize=7, length=3)
        ax.tick_params(axis="y", labelsize=8, length=3)
        ax.grid(axis="y", color="#d9d9d9", linewidth=0.6, alpha=0.8)
        ax.axvline(selected, color="black", linestyle="--", linewidth=1)
        sns.despine(ax=ax)
    save_pdf_svg(
        fig,
        output_root / "joint_nmf_k_selection",
        family="K-selection",
        function="plot_k_selection",
        inputs=[source],
        method="balanced_joint_nmf",
        notes=(
            "Tumor-only equal row weights; K 2..20, seeds 0..4. "
            "Squared reconstruction error remains in the source table and enters selection through explained fraction; "
            "it is not repeated as an independent plotted selection metric."
        ),
    )


def heatmap_scale(frame: pd.DataFrame, method: str) -> tuple[str, float | None, float | None, float | None, str]:
    if method in {"zscore", "robust"}:
        values = frame.to_numpy(float)
        values = values[np.isfinite(values)]
        vmax = float(np.nanpercentile(np.abs(values), 98)) if values.size else 1.0
        vmax = vmax if np.isfinite(vmax) and vmax > 0 else 1.0
        return "vlag", -vmax, vmax, 0.0, method
    if method == "standard_scale_col":
        return "viridis", 0.0, 1.0, None, "Column min-max"
    return "viridis", 0.0, None, None, "Raw"


def save_clustermap(
    frame: pd.DataFrame,
    stem: Path,
    *,
    method: str,
    value_label: str,
    family: str,
    source: Path,
    row_colors: pd.DataFrame | None = None,
    hide_row_labels: bool = False,
) -> None:
    frame = frame.loc[:, sorted(frame.columns.astype(str), key=cm_sort_key)]
    cmap, vmin, vmax, center, scale_label = heatmap_scale(frame, method)
    height = max(4.0, min(16.0, frame.shape[0] * 0.20 + 2.5))
    width = max(5.0, frame.shape[1] * 0.45 + 3.0)
    graph = sns.clustermap(
        frame,
        cmap=cmap, vmin=vmin, vmax=vmax, center=center,
        cbar_kws={"label": f"{scale_label} {value_label}"},
        figsize=(width, height), row_colors=row_colors,
        col_cluster=False, row_cluster=True, linewidths=0,
        colors_ratio=(0.06, 0.01),
    )
    if hide_row_labels:
        graph.ax_heatmap.set_yticklabels([])
        graph.ax_heatmap.set_yticks([])
        graph.ax_heatmap.tick_params(left=False)
        graph.ax_heatmap.set_ylabel("")
    graph.ax_heatmap.tick_params(axis="x", labelrotation=90)
    graph.fig.tight_layout()
    save_pdf_svg(graph.fig, stem, family=family, function="save_clustermap", inputs=[source], method=method, notes="CM axis not clustered")


def plot_activity_heatmaps() -> None:
    for method in ["raw", "zscore", "robust", "standard_scale_col"]:
        source = NMF_DIR / f"w_df_activity_sample_by_CM_{method}.csv"
        frame = pd.read_csv(source, index_col=0)
        save_clustermap(
            frame,
            FIGURE_ROOT / f"w_df_activity_sample_activity_per_CM_{method}_clustermap",
            method=method, value_label="activity", family="CM activity heatmap",
            source=source, hide_row_labels=True,
        )
    matrix_source = NMF_DIR / "w_df_activity_sample_by_CM_standard_scale_col_with_series_plot_matrix.csv"
    annotation_source = NMF_DIR / "w_df_activity_sample_by_CM_standard_scale_col_with_series_row_annotations.csv"
    frame = pd.read_csv(matrix_source, index_col=0)
    annotations = pd.read_csv(annotation_source, index_col=0).loc[frame.index]
    series_values = sorted(annotations["Series"].astype(str).unique())
    series_palette = dict(zip(series_values, sns.color_palette("Set2", len(series_values))))
    status_palette = {"tumor": "#E41A1C", "normal-like": "#377EB8"}
    row_colors = pd.DataFrame(
        {
            "Series": annotations["Series"].astype(str).map(series_palette),
            "Status": annotations["Status"].astype(str).map(status_palette),
        }, index=frame.index,
    )
    source = matrix_source
    frame = frame.loc[:, sorted(frame.columns.astype(str), key=cm_sort_key)]
    graph = sns.clustermap(
        frame, cmap="viridis", vmin=0, vmax=1,
        cbar_kws={"label": "Column min-max activity"},
        figsize=(max(6, 0.5 * frame.shape[1] + 3), 7),
        row_colors=row_colors, col_cluster=False, row_cluster=True,
        colors_ratio=(0.08, 0.01),
    )
    graph.ax_heatmap.set_yticklabels([]); graph.ax_heatmap.set_yticks([])
    graph.ax_heatmap.tick_params(left=False, axis="x", labelrotation=90)
    graph.ax_heatmap.set_ylabel("")
    handles = [mpl.patches.Patch(facecolor=color, label=label) for label, color in series_palette.items()]
    handles += [mpl.patches.Patch(facecolor=color, label=label) for label, color in status_palette.items() if label in set(annotations["Status"])]
    graph.ax_col_dendrogram.legend(handles=handles, title="Series / Status", frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    graph.fig.tight_layout()
    save_pdf_svg(
        graph.fig,
        FIGURE_ROOT / "w_df_activity_sample_activity_per_CM_standard_scale_col_with_series_clustermap",
        family="CM activity heatmap", function="plot_activity_heatmaps",
        inputs=[matrix_source, annotation_source], method="standard_scale_col",
        notes="Exactly Series and Status row annotation tracks; sample labels hidden; CM axis not clustered",
    )


def plot_loading_heatmaps() -> None:
    for method in ["raw", "zscore", "robust", "standard_scale_col"]:
        source = NMF_DIR / f"h_df_loading_cell_subtype_by_CM_{method}.csv"
        raw = pd.read_csv(source, index_col=0)
        cm_names = [value for value in raw.index.astype(str) if "CM" in value]
        if len(cm_names) == raw.shape[0]:
            frame = raw.T
            orientation = "CM x cell_subtype detected and transposed"
        else:
            frame = raw
            orientation = "cell_subtype x CM detected"
        save_clustermap(
            frame,
            FIGURE_ROOT / f"h_df_loading_cell_subtype_weights_per_CM_{method}_clustermap",
            method=method, value_label="loading", family="CM loading heatmap",
            source=source, hide_row_labels=False,
        )
        MANIFEST_ROWS[-1]["notes"] = orientation + "; CM axis not clustered"
        MANIFEST_ROWS[-2]["notes"] = orientation + "; CM axis not clustered"


def plot_top_subtype_heatmap() -> None:
    h_source = NMF_DIR / "H_df.csv"
    nodes_source = NODE_DIR / "joint_cm_cell_subtype_nodes_top20_from_H_df.csv"
    h_df = pd.read_csv(h_source, index_col=0)
    nodes = pd.read_csv(nodes_source)
    selected = nodes.loc[nodes["rank"] <= 12, "cell_subtype"].drop_duplicates().tolist()
    fraction = h_df.div(h_df.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    plot_df = fraction.loc[sorted(fraction.index, key=cm_sort_key), [x for x in selected if x in fraction.columns]]
    fig, ax = plt.subplots(figsize=(max(9, 0.22 * plot_df.shape[1] + 2), 5.5))
    sns.heatmap(plot_df, cmap="Reds", linewidths=0.2, linecolor="white", cbar_kws={"label": "Loading fraction"}, ax=ax)
    ax.set_xlabel("Cell subtype"); ax.set_ylabel("CM"); ax.tick_params(axis="x", rotation=90)
    fig.tight_layout()
    save_pdf_svg(fig, FIGURE_ROOT / "joint_module_top_subtype_heatmap", family="Top CM loadings", function="plot_top_subtype_heatmap", inputs=[h_source, nodes_source], method="loading_fraction", notes="Top 12 H-loading subtypes per CM; canonical CM order")


def closest_square_grid(n: int) -> tuple[int, int]:
    best = None
    for rows in range(1, int(math.ceil(math.sqrt(n))) + 1):
        cols = int(math.ceil(n / rows)); candidate = (abs(cols - rows), rows * cols - n, rows, cols)
        if best is None or candidate < best: best = candidate
    return best[2], best[3]


def circular_layout(nodes: list[str]) -> dict[str, tuple[float, float]]:
    return {node: (float(np.cos(2 * np.pi * i / max(1, len(nodes)))), float(np.sin(2 * np.pi * i / max(1, len(nodes))))) for i, node in enumerate(nodes)}


def draw_tumor_nodeplot(
    ax: mpl.axes.Axes,
    cm: str,
    nodes: list[str],
    edges: pd.DataFrame,
    colors: dict[str, object],
    edge_norm: Normalize,
    edge_cmap,
) -> set[str]:
    graph = nx.Graph(); graph.add_nodes_from(nodes); pos = circular_layout(nodes)
    passing = edges.loc[edges["edge_pass_r_ge_0.25"].astype(bool)]
    for row in passing.itertuples(index=False):
        if row.node_a in pos and row.node_b in pos:
            graph.add_edge(row.node_a, row.node_b, weight=float(row.pearson_r))
    nx.draw_networkx_edges(graph, pos, edge_color=[edge_cmap(edge_norm(d["weight"])) for _, _, d in graph.edges(data=True)], width=[2 + 3 * edge_norm(d["weight"]) for _, _, d in graph.edges(data=True)], alpha=0.88, ax=ax)
    nx.draw_networkx_nodes(graph, pos, node_color=[colors[n.split("_", 1)[0]] for n in nodes], node_size=1350, edgecolors="white", linewidths=0.8, ax=ax)
    nx.draw_networkx_labels(graph, pos, font_size=7, ax=ax)
    ax.set_title(cm, fontweight="bold"); ax.set_aspect("equal"); ax.axis("off")
    return {"tumor"} if graph.number_of_edges() else set()


def add_embedded_tumor_edge_colorbar(
    fig: mpl.figure.Figure,
    edge_norm: Normalize,
    edge_cmap,
    present_edge_classes: set[str],
) -> None:
    """Embed the live tumor edge mapping below a nodeplot when edges are drawn."""
    if "tumor" not in present_edge_classes:
        return
    cax = fig.add_axes([0.16, 0.055, 0.68, 0.022])
    sm = ScalarMappable(norm=edge_norm, cmap=edge_cmap)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label("Tumor edge Pearson r", fontsize=8)
    cb.ax.tick_params(labelsize=7, length=2)


def plot_nodeplots() -> None:
    node_source = NODE_DIR / "tumor_network_nodes_from_H_df.csv"
    edge_source = NODE_DIR / "status_specific_nodeplot_edges.csv"
    nodes = pd.read_csv(node_source).sort_values(["CM", "rank"])
    edges = pd.read_csv(edge_source)
    cm_order = sorted(nodes["CM"].unique(), key=cm_sort_key)
    prefixes = sorted({value.split("_", 1)[0] for value in nodes["cell_subtype"].astype(str)})
    palette = plt.get_cmap("tab20").colors
    colors = {prefix: palette[i % len(palette)] for i, prefix in enumerate(prefixes)}
    edge_norm = Normalize(vmin=EDGE_R_THRESHOLD, vmax=1.0)
    edge_cmap = mpl.colormaps["Reds"]
    for cm in cm_order:
        cm_nodes = nodes.loc[nodes["CM"].eq(cm), "cell_subtype"].astype(str).tolist()
        cm_edges = edges.loc[(edges["CM"].eq(cm)) & (edges["context"].eq("tumor"))]
        fig, ax = plt.subplots(figsize=(5.2, 5.7))
        present_edge_classes = draw_tumor_nodeplot(ax, cm, cm_nodes, cm_edges, colors, edge_norm, edge_cmap)
        fig.tight_layout(rect=[0, 0.16, 1, 1])
        add_embedded_tumor_edge_colorbar(fig, edge_norm, edge_cmap, present_edge_classes)
        save_pdf_svg(fig, FIGURE_ROOT / "tumor_nodeplots_by_cm" / f"{cm}_tumor_nodeplot", family="CM nodeplot", function="plot_nodeplots", inputs=[node_source, edge_source], method="tumor Pearson subtype correlation", notes="Tumor passing edges r>=0.25; deterministic circular loading-rank layout; embedded tumor Pearson-r colorbar uses the exact live edge colormap and Normalize")
    rows, cols = closest_square_grid(len(cm_order))
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 5.0, rows * 4.8 + 0.9), squeeze=False)
    overview_edge_classes: set[str] = set()
    for ax, cm in zip(axes.ravel(), cm_order):
        cm_nodes = nodes.loc[nodes["CM"].eq(cm), "cell_subtype"].astype(str).tolist()
        cm_edges = edges.loc[(edges["CM"].eq(cm)) & (edges["context"].eq("tumor"))]
        overview_edge_classes.update(draw_tumor_nodeplot(ax, cm, cm_nodes, cm_edges, colors, edge_norm, edge_cmap))
    for ax in axes.ravel()[len(cm_order):]: ax.axis("off")
    fig.tight_layout(rect=[0, 0.11, 1, 1])
    add_embedded_tumor_edge_colorbar(fig, edge_norm, edge_cmap, overview_edge_classes)
    save_pdf_svg(fig, FIGURE_ROOT / "tumor_all_CM_nodeplot", family="CM nodeplot overview", function="plot_nodeplots", inputs=[node_source, edge_source], method="tumor Pearson subtype correlation", notes="Closest-to-square grid; one subtype-subtype network panel per CM; one shared embedded tumor Pearson-r colorbar uses the exact live edge colormap and Normalize")
    handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=colors[p], markersize=10, label=p) for p in prefixes]
    fig, ax = plt.subplots(figsize=(max(4, len(prefixes) * 0.9), 1.4)); ax.axis("off"); ax.legend(handles=handles, ncol=min(4, len(prefixes)), frameon=False, loc="center", title="Subtype prefix")
    save_pdf_svg(fig, FIGURE_ROOT / "nodeplot_network_node_legend", family="CM nodeplot legend", function="plot_nodeplots", inputs=[node_source], method="prefix palette")
    fig, ax = plt.subplots(figsize=(4.5, 0.55)); sm = ScalarMappable(norm=edge_norm, cmap=edge_cmap); sm.set_array([]); cb = fig.colorbar(sm, cax=ax, orientation="horizontal"); cb.set_label("Tumor subtype-subtype Pearson r")
    save_pdf_svg(fig, FIGURE_ROOT / "nodeplot_network_edge_colorbar", family="CM nodeplot legend", function="plot_nodeplots", inputs=[edge_source], method="tumor Pearson subtype correlation")


def plot_corr_heatmap(frame: pd.DataFrame, stem: Path, title: str, inputs: list[Path]) -> None:
    fig, ax = plt.subplots(figsize=(max(5, 0.38 * len(frame)), max(4.5, 0.38 * len(frame))))
    sns.heatmap(frame, cmap="RdBu_r", center=0, vmin=-1, vmax=1, square=True, annot=False, linewidths=0.1, cbar_kws={"label": "Pearson r"}, ax=ax)
    ax.set_title(title); ax.tick_params(axis="x", rotation=90, labelsize=6); ax.tick_params(axis="y", rotation=0, labelsize=6)
    fig.tight_layout()
    save_pdf_svg(fig, stem, family="Top-node correlation heatmap", function="plot_corr_heatmap", inputs=inputs, method="Pearson", notes="Top10 H-loading diagnostic nodes; no edge filter; no cell text")


def plot_node_correlation_heatmaps(output_root: Path = FIGURE_ROOT) -> None:
    """Plot tumor top10 correlations as one CM-specific matrix per panel.

    The former union-node tumor matrix remains a provenance table but is not a
    valid input for this figure family. Every panel is built from the per-CM
    top10 H/loading diagnostic matrix computed within the tumor samples.
    """
    source_dir = NODE_DIR / "top10_node_correlation_matrices_no_edge_filter_by_cm"
    sources = sorted(
        source_dir.glob("*_top10_node_correlation_matrix_no_edge_filter.csv"),
        key=lambda path: cm_sort_key(path.name.split("_top10")[0]),
    )
    if len(sources) != 10:
        raise ValueError(f"Expected 10 tumor CM top10 correlation matrices, found {len(sources)}")

    panel_data: list[tuple[str, list[str], pd.DataFrame, Path]] = []
    for source in sources:
        cm = source.name.split("_top10")[0]
        matrix = pd.read_csv(source, index_col=0).astype(float)
        if matrix.shape != (10, 10) or matrix.index.astype(str).tolist() != matrix.columns.astype(str).tolist():
            raise ValueError(f"{cm}: expected aligned 10 x 10 top-node correlation matrix, found {matrix.shape}")
        nodes = matrix.index.astype(str).tolist()
        panel_data.append((cm, nodes, matrix, source))
        plot_corr_heatmap(
            matrix,
            output_root / "top10_node_correlation_heatmaps_no_edge_filter_by_cm" / "tumor" / f"{cm}_tumor_top10_node_correlation_heatmap_no_edge_filter",
            f"{cm}: tumor top10-node correlation",
            [source],
        )

    ncols = min(4, len(panel_data))
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
    for panel_i, (cm, nodes, matrix, _) in enumerate(panel_data):
        row, col = divmod(panel_i, ncols)
        ax = fig.add_subplot(grid[row, col])
        sns.heatmap(
            matrix,
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
        ax.set_title(cm)
        ax.set_xticklabels(nodes, rotation=90, ha="center", fontsize=6.5)
        ax.set_yticklabels(nodes, rotation=0, fontsize=6.5)
        ax.tick_params(axis="both", length=0)
    for empty_i in range(len(panel_data), nrows * ncols):
        row, col = divmod(empty_i, ncols)
        ax = fig.add_subplot(grid[row, col])
        ax.set_axis_off()
    fig.suptitle("Tumor: top10 node correlations by CM", y=0.995)
    fig.subplots_adjust(top=0.94, bottom=0.08, left=0.06, right=0.95)
    save_pdf_svg(
        fig,
        output_root / "tumor_top10_node_correlation_heatmap_no_edge_filter",
        family="Top-node correlation heatmap summary",
        function="plot_node_correlation_heatmaps",
        inputs=[item[3] for item in panel_data],
        method="tumor Pearson",
        notes="One separate 10 x 10 top-H-loading node matrix per CM panel; no union-node matrix; no edge filter; no cell text.",
    )


def run_skill_update_20260716() -> None:
    """Write the updated K/top10 figure contract without replacing older figures."""
    version = "skill-update-20260716"
    output_root = FIGURE_ROOT / version
    table_root = PLOT_TABLE_ROOT / version
    completion_path = table_root / "skill_update_completion.json"
    if completion_path.exists():
        completion = json.loads(completion_path.read_text())
        if completion.get("status") == "completed":
            print(json.dumps({"status": "valid_existing_skill_update_reused", "version": version}, indent=2))
            return
        raise FileExistsError(f"Existing incomplete skill update: {completion_path}")
    if output_root.exists() or table_root.exists():
        raise FileExistsError(f"Versioned skill-update output already exists without valid completion: {version}")

    started = time.time()
    output_root.mkdir(parents=True, exist_ok=False)
    table_root.mkdir(parents=True, exist_ok=False)
    plot_k_selection(output_root)
    plot_node_correlation_heatmaps(output_root)
    manifest = pd.DataFrame(MANIFEST_ROWS)
    manifest.to_csv(table_root / "figure_manifest.csv", index=False)
    expected_stems = 12  # one K diagnostic, one tumor summary, and ten per-CM panels
    pdf_count = int(manifest["output_format"].eq("pdf").sum())
    svg_count = int(manifest["output_format"].eq("svg").sum())
    all_files_ok = all(Path(path).exists() and Path(path).stat().st_size > 0 for path in manifest["figure_file"])
    validation = pd.DataFrame(
        [
            {"check": "expected_pdf_count", "observed": pdf_count, "expected": expected_stems, "passed": pdf_count == expected_stems},
            {"check": "expected_svg_count", "observed": svg_count, "expected": expected_stems, "passed": svg_count == expected_stems},
            {"check": "all_manifest_files_exist_nonempty", "observed": all_files_ok, "expected": True, "passed": all_files_ok},
            {"check": "tumor_per_cm_status_directory", "observed": len(list((output_root / "top10_node_correlation_heatmaps_no_edge_filter_by_cm" / "tumor").glob("*.pdf"))), "expected": 10, "passed": len(list((output_root / "top10_node_correlation_heatmaps_no_edge_filter_by_cm" / "tumor").glob("*.pdf"))) == 10},
            {"check": "normal_like_outputs_absent_for_tumor_only", "observed": len(list(output_root.rglob("*normal_like*"))), "expected": 0, "passed": len(list(output_root.rglob("*normal_like*"))) == 0},
        ]
    )
    validation.to_csv(table_root / "validation_report.csv", index=False)
    (table_root / "readme.txt").write_text(
        "Versioned 2026-07-16 skill-compliance plotting update. The original final figures are preserved. "
        "This directory makes the fixed K-selection reference plot and tumor-only top10 node-correlation figures canonical for the updated skill: "
        "one 10x10 matrix per CM panel, ten status-scoped per-CM PDF/SVG pairs, and no union-node summary matrix.\n",
        encoding="utf-8",
    )
    (table_root / "package_versions.txt").write_text(
        f"python={sys.version.split()[0]}\nmatplotlib={package_version('matplotlib')}\nseaborn={package_version('seaborn')}\npandas={package_version('pandas')}\ncode={CODE_PATH}\nseed={SEED}\n",
        encoding="utf-8",
    )
    completion = {
        "status": "completed" if bool(validation["passed"].all()) else "failed",
        "version": version,
        "figure_root": str(output_root),
        "manifest_rows": len(manifest),
        "pdf_count": pdf_count,
        "svg_count": svg_count,
        "n_failed_checks": int((~validation["passed"]).sum()),
        "old_outputs_preserved": True,
        "code_file": str(CODE_PATH),
        "elapsed_seconds": time.time() - started,
    }
    write_json(completion, completion_path)
    print(json.dumps(completion, indent=2), flush=True)
    if completion["status"] != "completed":
        raise SystemExit(1)


def plot_qstar_heatmap(value: pd.DataFrame, q: pd.DataFrame, stem: Path, title: str, cbar: str, inputs: list[Path]) -> None:
    columns = sorted(value.columns.astype(str), key=cm_sort_key)
    value = value.loc[:, columns]; q = q.reindex(index=value.index, columns=columns)
    annot = q.map(q_star)
    fig, ax = plt.subplots(figsize=(max(7, 0.65 * value.shape[1] + 2), max(5, 0.5 * value.shape[0] + 1.5)))
    sns.heatmap(value.astype(float), cmap="coolwarm", center=0, vmin=-1, vmax=1, annot=annot, fmt="", annot_kws={"fontsize": 7}, linewidths=0.25, linecolor="white", cbar_kws={"label": cbar}, ax=ax)
    ax.set_title(title); ax.set_xlabel("CM"); ax.set_ylabel("Epithelial subtype"); ax.tick_params(axis="x", rotation=90); ax.tick_params(axis="y", rotation=0)
    fig.tight_layout()
    save_pdf_svg(fig, stem, family="Epi-CM association heatmap", function="plot_qstar_heatmap", inputs=inputs, method="Spearman", notes="Cell text is BH-FDR q significance only: ns/*/**/***")


def plot_epi_cm_heatmaps() -> None:
    rho_source = SPEARMAN_DIR / "balanced_joint_cm_epi_cm_association_tumor_rho_matrix.csv"
    q_source = SPEARMAN_DIR / "balanced_joint_cm_epi_cm_association_tumor_q_matrix.csv"
    plot_qstar_heatmap(pd.read_csv(rho_source, index_col=0), pd.read_csv(q_source, index_col=0), FIGURE_ROOT / "balanced_joint_cm_epi_cm_association_tumor_heatmap_qstars", "Balanced joint CM-Epi association (tumor)", "Spearman rho", [rho_source, q_source])


def plot_all_spearman_scatters() -> None:
    association_source = SPEARMAN_DIR / "balanced_joint_cm_epi_cm_association_tumor_spearman_long.csv"
    scatter_source = SPEARMAN_DIR / "balanced_joint_cm_epi_cm_association_tumor_spearman_scatter_source.csv"
    palette_source = Path("/mnt/disk18t/lr_xcy/riku/brca_cm-epi_codex/epi-cm-core-workflow/tables/02-cell_subtype_integration_clustering/06-project-subtypes-to-full-adata/cell_subtype_color_mapping.csv")
    association = pd.read_csv(association_source)
    source = pd.read_csv(scatter_source)
    palette_table = pd.read_csv(palette_source)
    palette = dict(zip(palette_table["cell_subtype"].astype(str), palette_table["hex_color"].astype(str)))
    missing = sorted(set(association["epi_subtype"].astype(str)) - set(palette))
    if missing: raise ValueError(f"Epithelial subtype palette missing: {missing}")
    out_dir = FIGURE_ROOT / "epi_cm_association_spearman" / "tumor" / "scatterplots"
    for row in association.itertuples(index=False):
        pair = source.loc[(source["epi_subtype"].eq(row.epi_subtype)) & (source["CM"].eq(row.CM))]
        fig, ax = plt.subplots(figsize=(6, 5))
        sns.regplot(data=pair, x="cm_activity", y="epi_frequency", ci=95, scatter_kws={"s": 60, "alpha": 0.85, "edgecolor": "none", "color": palette[row.epi_subtype]}, line_kws={"linewidth": 2, "color": palette[row.epi_subtype]}, ax=ax)
        ax.set_ylim(0, 1); ax.set_xlabel(f"{row.CM} activity"); ax.set_ylabel(f"{row.epi_subtype} frequency")
        q_text = "NA" if pd.isna(row.q_value) else f"{row.q_value:.3g}"
        rho_text = "NA" if pd.isna(row.rho) else f"{row.rho:.3f}"
        ax.text(0.03, 0.97, f"Spearman rho={rho_text}\nq={q_text}\nn={int(row.n_samples)}", transform=ax.transAxes, va="top", ha="left", bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.9, "edgecolor": "#cccccc"})
        sns.despine(ax=ax); fig.tight_layout()
        stem = out_dir / f"scatter_{safe_part(row.CM)}_vs_{safe_part(row.epi_subtype)}"
        save_pdf_svg(fig, stem, family="Epi-CM all-pair scatter", function="plot_all_spearman_scatters", inputs=[association_source, scatter_source, palette_source], method="Spearman", notes="All tumor epithelial subtype x CM pairs; seaborn regplot ci=95; y fixed 0..1; q annotation")


def main() -> None:
    started = time.time()
    completion_path = PLOT_TABLE_ROOT / "final_plotting_completion.json"
    if completion_path.exists():
        completion = json.loads(completion_path.read_text())
        if completion.get("status") == "completed":
            print(json.dumps({"status": "valid_existing_final_plots_reused"}, indent=2))
            return
        raise FileExistsError("Existing final plotting completion is invalid.")
    FIGURE_ROOT.mkdir(parents=True, exist_ok=True); PLOT_TABLE_ROOT.mkdir(parents=True, exist_ok=True)
    required = [
        NMF_DIR / "joint_nmf_k_selection_metrics.csv", NMF_DIR / "W_df.csv",
        NMF_DIR / "H_df.csv", NODE_DIR / "tumor_network_nodes_from_H_df.csv",
        NODE_DIR / "status_specific_nodeplot_edges.csv",
        SPEARMAN_DIR / "balanced_joint_cm_epi_cm_association_tumor_spearman_long.csv",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing: raise FileNotFoundError(f"Missing canonical plotting inputs: {missing}")
    plot_k_selection(); plot_activity_heatmaps(); plot_loading_heatmaps()
    plot_top_subtype_heatmap(); plot_nodeplots(); plot_node_correlation_heatmaps()
    plot_epi_cm_heatmaps(); plot_all_spearman_scatters()
    manifest = pd.DataFrame(MANIFEST_ROWS)
    manifest.to_csv(PLOT_TABLE_ROOT / "figure_manifest.csv", index=False)
    pdf_count = int(manifest["output_format"].eq("pdf").sum())
    svg_count = int(manifest["output_format"].eq("svg").sum())
    scatter_pdf = int((manifest["figure_family"].eq("Epi-CM all-pair scatter") & manifest["output_format"].eq("pdf")).sum())
    if pdf_count != svg_count or scatter_pdf != 100:
        raise ValueError(f"Figure inventory incomplete: pdf={pdf_count}, svg={svg_count}, scatter_pdf={scatter_pdf}")
    validation = pd.DataFrame(
        [
            {"check": "all_manifest_files_exist_nonempty", "passed": all(Path(p).exists() and Path(p).stat().st_size > 0 for p in manifest["figure_file"])},
            {"check": "pdf_svg_counts_match", "passed": pdf_count == svg_count},
            {"check": "all_100_spearman_scatter_pairs_plotted", "passed": scatter_pdf == 100},
            {"check": "normal_like_figures_skipped_for_tumor_only", "passed": not any("normal_like" in Path(p).name for p in manifest["figure_file"])},
            {"check": "pearson_figures_absent_by_default", "passed": not any("pearson" in str(p).lower() for p in manifest["figure_file"])},
        ]
    )
    validation.to_csv(PLOT_TABLE_ROOT / "validation_report.csv", index=False)
    (PLOT_TABLE_ROOT / "package_versions.txt").write_text(
        f"python={sys.version.split()[0]}\nmatplotlib={package_version('matplotlib')}\nseaborn={package_version('seaborn')}\nnetworkx={package_version('networkx')}\npandas={package_version('pandas')}\ncode={CODE_PATH}\nseed={SEED}\n",
        encoding="utf-8",
    )
    pd.DataFrame([{"detected_mode": "tumor_only", "final_plotting_script": str(CODE_PATH), "figure_root": str(FIGURE_ROOT), "pdf_count": pdf_count, "svg_count": svg_count, "spearman_scatter_pair_count": scatter_pdf, "normal_like_outputs_skipped": True, "tumor_normal_comparison_outputs_skipped": True, "pearson_optional_branch_skipped": True, "seed": SEED}]).to_csv(PLOT_TABLE_ROOT / "final_plotting_parameters.csv", index=False)
    completion = {"status": "completed", "detected_mode": "tumor_only", "manifest_rows": len(manifest), "pdf_count": pdf_count, "svg_count": svg_count, "spearman_scatter_pair_count": scatter_pdf, "all_manifest_files_exist_nonempty": bool(validation["passed"].all()), "normal_like_and_comparison_outputs_skipped_as_structurally_inapplicable": True, "pearson_optional_branch_skipped": True, "elapsed_seconds": time.time() - started}
    write_json(completion, completion_path)
    print(json.dumps(completion, indent=2), flush=True)


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "--skill-update-20260716":
        run_skill_update_20260716()
    elif len(sys.argv) == 1:
        main()
    else:
        raise SystemExit("Usage: 06_final_unified_plotting_and_manifest.py [--skill-update-20260716]")
