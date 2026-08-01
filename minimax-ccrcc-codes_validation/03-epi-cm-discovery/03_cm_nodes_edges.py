"""Block 03 step 3: CM Nodes and Edges (per SKILL.md).

Canonical:
  - top_n_subtypes = 20, top_n_nodes = 10, edge_r_threshold = 0.25
  - Per-CM node candidate from H_df.csv (top loading)
  - Pearson node x node correlation in tumor and normal-like
  - Retain nodes with at least one edge r>=0.25 in either context
  - canonical node table columns: CM, reference_node_rank, node (no extra cols)
  - status_specific_nodeplot_edges.csv: context, CM, node_a, node_b, pearson_r, edge_pass_r_ge_0.25
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

SEED = 42
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
TAB = ROOT / "epi-cm-core-workflow/tables/03-epi-cm-discovery"
H_CANON = TAB / "H_df.csv"
LOADING = TAB / "loading_df_cell_subtype_by_CM.csv"
NON_EPI_FREQ = TAB / "non_epi_subtype_frequency.csv"
SAMPLE_INCL = TAB / "sample_inclusion_exclusion.csv"
SAMPLE_STATUS = TAB / "sample_status.csv"

TOP_N_SUBTYPES = 20
TOP_N_NODES = 10
EDGE_R_THRESHOLD = 0.25


def main() -> None:
    t0 = time.time()
    H = pd.read_csv(H_CANON, index_col=0)
    loading = pd.read_csv(LOADING, index_col=0)
    non_epi_freq = pd.read_csv(NON_EPI_FREQ, index_col=0)
    incl = pd.read_csv(SAMPLE_INCL, index_col=0)
    status = pd.read_csv(SAMPLE_STATUS, index_col=0)["status"]
    print(f"[node] H: {H.shape}, loading: {loading.shape}", flush=True)

    canonical_cms = H.index.astype(str).tolist()
    candidate_nodes = {}
    for cm in canonical_cms:
        h = H.loc[cm].sort_values(ascending=False)
        top = h.head(TOP_N_SUBTYPES).index.astype(str).tolist()
        candidate_nodes[cm] = top[:TOP_N_NODES]

    keep_samples = incl.index[incl["keep_for_epi_cm"]].astype(str).tolist()
    keep_samples = [s for s in keep_samples if s in non_epi_freq.index and s in status.index]
    is_tumor = set([s for s in keep_samples if status.loc[s] == "tumor"])
    is_normal = set([s for s in keep_samples if status.loc[s] == "normal-like"])
    print(f"[node] tumor={len(is_tumor)}, normal={len(is_normal)}", flush=True)

    status_edges = []
    node_node_rho = {}
    node_node_q = {}
    for cm in canonical_cms:
        cand = candidate_nodes[cm]
        for context_name, sample_set in [("tumor", is_tumor), ("normal-like", is_normal)]:
            samples_ctx = sorted([s for s in keep_samples if s in sample_set])
            sub_freq = non_epi_freq.loc[samples_ctx, cand]
            n = len(cand)
            rmat = np.full((n, n), np.nan)
            pmat = np.full((n, n), np.nan)
            for i in range(n):
                for j in range(n):
                    if i == j:
                        rmat[i, j] = 1.0
                        pmat[i, j] = 0.0
                        continue
                    a = sub_freq.iloc[:, i].values
                    b = sub_freq.iloc[:, j].values
                    if a.std() == 0 or b.std() == 0:
                        continue
                    r, p = pearsonr(a, b)
                    rmat[i, j] = r
                    pmat[i, j] = p
            rdf = pd.DataFrame(rmat, index=cand, columns=cand)
            pdf = pd.DataFrame(pmat, index=cand, columns=cand)
            node_node_rho[(cm, context_name)] = rdf
            node_node_q[(cm, context_name)] = pdf
            for i in range(n):
                for j in range(i + 1, n):
                    r = rmat[i, j]
                    if np.isnan(r):
                        continue
                    p = pmat[i, j] if not np.isnan(pmat[i, j]) else 1.0
                    status_edges.append({
                        "context": context_name,
                        "CM": cm,
                        "node_a": cand[i],
                        "node_b": cand[j],
                        "pearson_r": float(r),
                        "pval": float(p),
                        "edge_pass_r_ge_0.25": bool(r >= EDGE_R_THRESHOLD),
                    })

    retained_nodes = {}
    for cm in canonical_cms:
        cand = candidate_nodes[cm]
        sub = pd.DataFrame(status_edges)
        sub = sub[(sub["CM"] == cm) & (sub["edge_pass_r_ge_0.25"])]
        passed = set(sub["node_a"].tolist()) | set(sub["node_b"].tolist())
        retained = [n for n in cand if n in passed]
        retained_nodes[cm] = retained

    # Canonical node membership table (only 3 columns)
    rows = []
    for cm in canonical_cms:
        for rank, n in enumerate(retained_nodes[cm], 1):
            rows.append({"CM": cm, "reference_node_rank": rank, "node": n})
    pd.DataFrame(rows).to_csv(TAB / "balanced_joint_cm_reference_node_sets_after_edge_threshold.csv", index=False)

    # Status-specific node tables
    tumor_nodes = []
    normal_nodes = []
    for cm in canonical_cms:
        for n in retained_nodes[cm]:
            r = float(H.loc[cm, n])
            tumor_nodes.append({"CM": cm, "node": n, "loading_in_H": r, "context": "tumor"})
            normal_nodes.append({"CM": cm, "node": n, "loading_in_H": r, "context": "normal-like"})
    pd.DataFrame(tumor_nodes).to_csv(TAB / "tumor_network_nodes_from_H_df.csv", index=False)
    pd.DataFrame(normal_nodes).to_csv(TAB / "normal_like_network_nodes_from_H_df.csv", index=False)

    # status_specific_nodeplot_edges.csv
    kept_pairs = set()
    for cm in canonical_cms:
        rnodes = retained_nodes[cm]
        for a in rnodes:
            for b in rnodes:
                if a < b:
                    kept_pairs.add((cm, (a, b)))
    edges_kept = [r for r in status_edges if (r["CM"], tuple(sorted((r["node_a"], r["node_b"])))) in kept_pairs]
    edges_df = pd.DataFrame(edges_kept)
    edges_df.to_csv(TAB / "status_specific_nodeplot_edges.csv", index=False)

    # Per-CM context node-node correlation matrices (on retained nodes)
    for cm in canonical_cms:
        rnodes = retained_nodes[cm]
        for context_name in ["tumor", "normal-like"]:
            rdf = node_node_rho[(cm, context_name)].loc[rnodes, rnodes]
            qdf = node_node_q[(cm, context_name)].loc[rnodes, rnodes]
            tag = context_name.replace("-", "_")
            rdf.to_csv(TAB / f"node_node_correlation_{tag}_{cm}.csv")
            qdf.to_csv(TAB / f"node_node_q_{tag}_{cm}.csv")

    summary = {
        "n_CMs": len(canonical_cms),
        "top_n_subtypes": TOP_N_SUBTYPES,
        "top_n_nodes": TOP_N_NODES,
        "edge_r_threshold": EDGE_R_THRESHOLD,
        "n_total_edges": len(edges_kept),
        "n_retained_nodes_total": len(rows),
        "retained_per_cm": {cm: len(retained_nodes[cm]) for cm in canonical_cms},
    }
    with open(TAB / "cm_nodes_edges_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[node] retained={len(rows)} nodes, edges={len(edges_kept)}", flush=True)
    print(f"[node] done. total={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()