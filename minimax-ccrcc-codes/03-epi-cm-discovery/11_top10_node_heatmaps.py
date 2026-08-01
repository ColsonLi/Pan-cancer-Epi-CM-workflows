"""Block 03 final: top10_node_correlation_heatmaps (no edge filter) per SKILL.md.

Per updated SKILL.md (module 03):
  - Node-node correlation heatmap: use the top10 H-loading diagnostic nodes for
    each CM, regardless of whether those nodes survive edge filtering or appear
    in the nodeplot. Keep row and column order identical for square matrices,
    use a centered diverging palette fixed at -1 to 1, show colorbar ticks,
    and do not annotate cells with numeric correlation values. By default leave
    heatmap cells text-free; if the user explicitly requests significance
    labels, overlay q-value symbols only from the matching q table.
  - Source tables:
      `node_node_correlation_matrix.csv`  (Pearson r, pre-computed)
      `node_node_correlation_q_matrix.csv` (BH q values, pre-computed)
  - Record the selected correlation method (here: Pearson).
"""
from __future__ import annotations

import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests

SEED = 42
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
TAB = ROOT / "epi-cm-core-workflow/tables/03-epi-cm-discovery"
FIG = ROOT / "epi-cm-core-workflow/figures/03-epi-cm-discovery"
H_CANON = TAB / "H_df.csv"
NON_EPI_FREQ = TAB / "non_epi_subtype_frequency.csv"
SAMPLE_INCL = TAB / "sample_inclusion_exclusion.csv"
SAMPLE_STATUS = TAB / "sample_status.csv"

CORRELATION_METHOD = "pearson"  # documented per SKILL "Record the selected correlation method"


def build_top10_table(H: pd.DataFrame, n: int = 10) -> pd.DataFrame:
    """FIXED: top10 H-loading diagnostic nodes per CM (context-free diagnostic)."""
    rows = []
    for cm in H.index.astype(str):
        h = H.loc[cm].sort_values(ascending=False).head(n)
        for rank, (node, loading_val) in enumerate(h.items(), 1):
            rows.append({
                "CM": cm,
                "cell_subtype": str(node),
                "loading": float(loading_val),
                "rank": rank,
            })
    return pd.DataFrame(rows)


def compute_node_node_pearson(non_epi_freq: pd.DataFrame, samples_ctx, top10_nodes):
    """FIXED: pre-compute node x node Pearson matrix for the top10 nodes."""
    sub = non_epi_freq.loc[samples_ctx, top10_nodes]
    n = len(top10_nodes)
    rmat = np.eye(n)
    pmat = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            a = sub.iloc[:, i].values
            b = sub.iloc[:, j].values
            if a.std() == 0 or b.std() == 0:
                continue
            r, p = pearsonr(a, b)
            rmat[i, j] = r
            pmat[i, j] = p
    return pd.DataFrame(rmat, index=top10_nodes, columns=top10_nodes), \
           pd.DataFrame(pmat, index=top10_nodes, columns=top10_nodes)


def write_precomputed_matrices(corr_by_status):
    """FIXED: write consolidated node_node_correlation_matrix.csv and q_matrix.csv."""
    rows = []
    for status_ctx, cm_dict in corr_by_status.items():
        for cm, (rmat, pmat) in cm_dict.items():
            for i, a in enumerate(rmat.index):
                for j, b in enumerate(rmat.columns):
                    r = float(rmat.iloc[i, j])
                    p = float(pmat.iloc[i, j])
                    rows.append({
                        "CM": cm, "status_context": status_ctx,
                        "node_a": a, "node_b": b,
                        "pearson_r": r, "p_value": p,
                    })
    long_df = pd.DataFrame(rows)

    # BH q across all pairs (status_context × CM × pairs).
    long_df["q_value"] = np.nan
    for (ctx, cm), idx in long_df.groupby(["status_context", "CM"], sort=False).groups.items():
        p = long_df.loc[idx, "p_value"].to_numpy(float)
        valid = np.isfinite(p)
        q = np.full(p.shape, np.nan, dtype=float)
        if valid.any():
            q[valid] = multipletests(p[valid], method="fdr_bh")[1]
        long_df.loc[idx, "q_value"] = q

    long_df.to_csv(TAB / "node_node_correlation_matrix.csv", index=False)

    wide = long_df.pivot_table(
        index=["status_context", "CM", "node_a"],
        columns="node_b",
        values="pearson_r",
    )
    wide.to_csv(TAB / "node_node_correlation_matrix_wide.csv")
    wide_q = long_df.pivot_table(
        index=["status_context", "CM", "node_a"],
        columns="node_b",
        values="q_value",
    )
    wide_q.to_csv(TAB / "node_node_correlation_q_matrix.csv")
    return long_df


def plot_one(rmat, stem, title):
    fig, ax = plt.subplots(figsize=(max(4.2, 0.48 * len(rmat)), max(4.0, 0.48 * len(rmat))))
    sns.heatmap(rmat, cmap="RdBu_r", center=0, vmin=-1, vmax=1, square=True,
                linewidths=0.2, linecolor="white", annot=False,
                cbar_kws={"label": f"{CORRELATION_METHOD.title()} r"}, ax=ax)
    ax.set_title(title)
    ax.set_xticklabels(rmat.columns, rotation=90, ha="center", fontsize=7)
    ax.set_yticklabels(rmat.index, rotation=0, fontsize=7)
    ax.tick_params(axis="both", length=0)
    fig.tight_layout()
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", dpi=300)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", dpi=300)
    plt.close(fig)


def plot_top10_node_correlation_summary(
    cm_matrices, top_nodes, cm_order, stem, context_label, ncols=4,
):
    """FIXED per SKILL: multi-panel summary, one CM-specific matrix per panel.

    `cm_matrices` is {cm: (rmat, pmat)} for a single status_context.
    Per SKILL:
      - "The summary heatmap files are fixed multi-panel figures, not one merged
        all-node matrix."
      - "draw one separate square panel per CM"
      - "Each panel must contain only that CM's own top10 H/loading diagnostic
        nodes and their within-status 10 x 10 Pearson correlation matrix."
      - "Never take the union of nodes from different CMs"
    """
    panel_data = []
    for cm in cm_order:
        nodes = (
            top_nodes.loc[top_nodes["CM"].astype(str).eq(str(cm))]
            .sort_values("rank", kind="stable").head(10)["cell_subtype"].astype(str).tolist()
        )
        if not nodes:
            continue
        rmat_pcm, _ = cm_matrices[cm]
        nodes_present = [n for n in nodes if n in rmat_pcm.index and n in rmat_pcm.columns]
        if len(nodes_present) < 2:
            continue
        panel_data.append((cm, nodes_present, rmat_pcm.loc[nodes_present, nodes_present].astype(float)))

    if not panel_data:
        raise ValueError(f"No CM panels available for {context_label}")

    ncols = min(ncols, len(panel_data))
    nrows = int(np.ceil(len(panel_data) / ncols))
    panel_size = 4.0
    fig = plt.figure(figsize=(panel_size * ncols + 0.8, panel_size * nrows + 0.8))
    grid = fig.add_gridspec(
        nrows, ncols + 1,
        width_ratios=[1.0] * ncols + [0.05],
        wspace=0.65, hspace=0.75,
    )
    cbar_ax = fig.add_subplot(grid[:, -1])

    for panel_i, (cm, nodes, plot_df) in enumerate(panel_data):
        row, col = divmod(panel_i, ncols)
        ax = fig.add_subplot(grid[row, col])
        sns.heatmap(
            plot_df, cmap="RdBu_r", center=0, vmin=-1, vmax=1, square=True,
            linewidths=0.2, linecolor="white", annot=False,
            cbar=(panel_i == 0), cbar_ax=cbar_ax if panel_i == 0 else None,
            cbar_kws={"label": f"{CORRELATION_METHOD.title()} r"} if panel_i == 0 else None,
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

    fig.suptitle(f"{context_label}", y=0.995, fontsize=12)
    fig.subplots_adjust(top=0.94, bottom=0.08, left=0.06, right=0.95)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", dpi=300)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", dpi=300)
    plt.close(fig)


def main():
    t0 = time.time()
    H = pd.read_csv(H_CANON, index_col=0)
    non_epi_freq = pd.read_csv(NON_EPI_FREQ, index_col=0)
    incl = pd.read_csv(SAMPLE_INCL, index_col=0)
    status = pd.read_csv(SAMPLE_STATUS, index_col=0)["status"]
    canonical_cms = H.index.astype(str).tolist()

    keep_samples = incl.index[incl["keep_for_cm"]].astype(str).tolist()
    keep_samples = [s for s in keep_samples if s in non_epi_freq.index and s in status.index]
    is_tumor = [s for s in keep_samples if status.loc[s] == "tumor"]
    is_normal = [s for s in keep_samples if status.loc[s] == "normal-like"]
    print(f"[heat] tumor={len(is_tumor)}, normal={len(is_normal)}", flush=True)

    # Top10 diagnostic table.
    top10_table = build_top10_table(H, n=10)
    top10_table.to_csv(TAB / "joint_cm_cell_subtype_nodes_top10_from_H_df.csv", index=False)
    print(f"[heat] top10 table: {len(top10_table)} rows", flush=True)
    cm_top10 = {cm: top10_table.loc[top10_table["CM"] == cm, "cell_subtype"].astype(str).tolist()
                for cm in canonical_cms}

    # Pre-compute correlations.
    corr_by_status = {
        "tumor": {cm: compute_node_node_pearson(non_epi_freq, is_tumor, cm_top10[cm])
                  for cm in canonical_cms},
        "normal-like": {cm: compute_node_node_pearson(non_epi_freq, is_normal, cm_top10[cm])
                        for cm in canonical_cms},
    }
    write_precomputed_matrices(corr_by_status)
    print(f"[heat] wrote node_node_correlation_matrix.csv + q_matrix.csv", flush=True)

    # Per-CM heatmaps.
    out_dir = FIG / "top10_node_correlation_heatmaps_no_edge_filter_by_cm"
    (out_dir / "tumor").mkdir(parents=True, exist_ok=True)
    (out_dir / "normal-like").mkdir(parents=True, exist_ok=True)
    for cm in canonical_cms:
        for status_ctx in ["tumor", "normal-like"]:
            rmat, _ = corr_by_status[status_ctx][cm]
            if len(rmat) < 2:
                continue
            sub_dir = out_dir / status_ctx
            plot_one(
                rmat,
                sub_dir / f"{cm}_{status_ctx.replace('-', '_')}_top10_node_correlation_heatmap_no_edge_filter",
                f"{cm} {status_ctx} top10 node {CORRELATION_METHOD.title()} r",
            )

    # Multi-panel summary heatmaps: one panel per CM (NOT a merged all-CMs matrix).
    # Per SKILL: "Never take the union of nodes from different CMs to draw one
    # giant correlation matrix". Each panel shows that CM's own top10 nodes.
    for status_ctx in ["tumor", "normal-like"]:
        stem = FIG / f"{status_ctx.replace('-', '_')}_top10_node_correlation_heatmap_no_edge_filter"
        cm_matrices = corr_by_status[status_ctx]
        plot_top10_node_correlation_summary(
            cm_matrices, top10_table, canonical_cms, stem,
            context_label=f"{status_ctx} top10 node {CORRELATION_METHOD.title()} r (no edge filter)",
            ncols=4,
        )

    print(f"[heat] done. total={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()