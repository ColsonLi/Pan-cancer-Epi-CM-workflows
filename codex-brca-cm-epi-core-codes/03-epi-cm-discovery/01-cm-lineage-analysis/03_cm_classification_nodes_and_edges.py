#!/usr/bin/env python3
"""Build canonical tumor-only CM nodes, edges, and correlation tables."""

from __future__ import annotations

import importlib.metadata
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from cm_analysis_common import (
    EDGE_R_THRESHOLD, NMF_DIR, NODE_DIR, PREP_DIR, SEED, TOP_N_NODES,
    TOP_N_SUBTYPES, bh_fdr, cm_sort_key, write_json,
)


CODE_PATH = Path(__file__).resolve()


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def subtype_lineage(value: str) -> str:
    return str(value).split("_", 1)[0]


def corr_and_q(frame: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    columns = frame.columns.astype(str).tolist()
    corr = pd.DataFrame(np.eye(len(columns)), index=columns, columns=columns, dtype=float)
    pmat = pd.DataFrame(np.zeros((len(columns), len(columns))), index=columns, columns=columns, dtype=float)
    pairs: list[tuple[str, str, float]] = []
    for a, b in itertools.combinations(columns, 2):
        if frame[a].nunique() < 2 or frame[b].nunique() < 2:
            r, p = np.nan, np.nan
        else:
            r, p = pearsonr(frame[a], frame[b])
        corr.loc[a, b] = corr.loc[b, a] = r
        pmat.loc[a, b] = pmat.loc[b, a] = p
        pairs.append((a, b, p))
    pair_frame = pd.DataFrame(pairs, columns=["node_a", "node_b", "p_value"])
    pair_frame["q_value"] = bh_fdr(pair_frame["p_value"])
    qmat = pd.DataFrame(np.zeros((len(columns), len(columns))), index=columns, columns=columns, dtype=float)
    for row in pair_frame.itertuples(index=False):
        qmat.loc[row.node_a, row.node_b] = qmat.loc[row.node_b, row.node_a] = row.q_value
    return corr, pmat, qmat


def main() -> None:
    started = time.time()
    completion_path = NODE_DIR / "nodes_edges_completion.json"
    if completion_path.exists():
        completion = json.loads(completion_path.read_text())
        if completion.get("status") == "completed":
            print(json.dumps({"status": "valid_existing_nodes_edges_reused"}, indent=2))
            return
        raise FileExistsError("Existing nodes/edges completion is invalid.")
    NODE_DIR.mkdir(parents=True, exist_ok=True)
    mode = json.loads((PREP_DIR / "cm_status_mode.json").read_text())["detected_mode"]
    if mode != "tumor_only":
        raise ValueError(f"Expected tumor_only route, found {mode}")
    h_df = pd.read_csv(NMF_DIR / "H_df.csv", index_col=0)
    h_df.index = h_df.index.astype(str)
    frequency = pd.read_csv(NMF_DIR / "non_epi_subtype_frequency_column_minmax.csv", index_col=0)
    sample_status = pd.read_csv(PREP_DIR / "sample_status.csv", index_col=0)
    if set(sample_status["status"].astype(str)) != {"tumor"}:
        raise ValueError("Tumor node analysis received non-tumor samples.")
    if set(h_df.columns.astype(str)) != set(frequency.columns.astype(str)):
        raise ValueError("H subtype axis differs from column-minmax frequency columns.")
    cm_order = sorted(h_df.index.tolist(), key=cm_sort_key)
    h_df = h_df.loc[cm_order]

    all_rows: list[dict[str, object]] = []
    for cm in cm_order:
        ranked = h_df.loc[cm].sort_values(ascending=False)
        for rank, (subtype, loading) in enumerate(ranked.items(), start=1):
            all_rows.append(
                {
                    "CM": cm, "cell_subtype": subtype,
                    "loading": float(loading), "rank": rank,
                    "cell_lineage": subtype_lineage(subtype),
                    "status_context": "tumor",
                }
            )
    nodes_all = pd.DataFrame(all_rows)
    nodes_top20 = nodes_all.loc[nodes_all["rank"] <= TOP_N_SUBTYPES].copy()
    nodes_top10 = nodes_all.loc[nodes_all["rank"] <= TOP_N_NODES].copy()
    nodes_all.to_csv(NODE_DIR / "joint_cm_cell_subtype_nodes_all_from_H_df.csv", index=False)
    nodes_top20.to_csv(NODE_DIR / "joint_cm_cell_subtype_nodes_top20_from_H_df.csv", index=False)
    nodes_top10.to_csv(NODE_DIR / "joint_cm_cell_subtype_nodes_top10_from_H_df.csv", index=False)

    membership_rows: list[dict[str, object]] = []
    membership_diagnostic: list[dict[str, object]] = []
    edge_rows: list[dict[str, object]] = []
    cm_audit: list[dict[str, object]] = []
    per_cm_corr_dir = NODE_DIR / "top10_node_correlation_matrices_no_edge_filter_by_cm"
    per_cm_corr_dir.mkdir(parents=True, exist_ok=True)
    for cm in cm_order:
        top = nodes_top10.loc[nodes_top10["CM"].eq(cm)].sort_values("rank")
        candidates = top["cell_subtype"].astype(str).tolist()
        corr, pmat, qmat = corr_and_q(frequency[candidates])
        corr.to_csv(per_cm_corr_dir / f"{cm}_top10_node_correlation_matrix_no_edge_filter.csv")
        pmat.to_csv(per_cm_corr_dir / f"{cm}_top10_node_correlation_p_matrix_no_edge_filter.csv")
        qmat.to_csv(per_cm_corr_dir / f"{cm}_top10_node_correlation_q_matrix_no_edge_filter.csv")
        passing_pairs = {
            tuple(sorted((a, b)))
            for a, b in itertools.combinations(candidates, 2)
            if np.isfinite(corr.loc[a, b]) and float(corr.loc[a, b]) >= EDGE_R_THRESHOLD
        }
        retained = [
            node for node in candidates
            if any(node in pair for pair in passing_pairs)
        ]
        for rank, node in enumerate(retained, start=1):
            membership_rows.append({"CM": cm, "reference_node_rank": rank, "node": node})
            h_rank = int(top.loc[top["cell_subtype"].eq(node), "rank"].iloc[0])
            loading = float(top.loc[top["cell_subtype"].eq(node), "loading"].iloc[0])
            membership_diagnostic.append(
                {"CM": cm, "reference_node_rank": rank, "node": node,
                 "loading_rank_in_H_top10": h_rank, "loading": loading,
                 "cell_lineage": subtype_lineage(node)}
            )
        for a, b in itertools.combinations(retained, 2):
            r = float(corr.loc[a, b])
            edge_rows.append(
                {"context": "tumor", "CM": cm, "node_a": a, "node_b": b,
                 "pearson_r": r,
                 "edge_pass_r_ge_0.25": bool(np.isfinite(r) and r >= EDGE_R_THRESHOLD)}
            )
        cm_audit.append(
            {"CM": cm, "n_top10_candidates": len(candidates),
             "n_passing_edges": len(passing_pairs), "n_retained_nodes": len(retained),
             "panel_drawable": len(retained) >= 2 and len(passing_pairs) >= 1,
             "status_context": "tumor"}
        )

    membership = pd.DataFrame(membership_rows, columns=["CM", "reference_node_rank", "node"])
    if list(membership.columns) != ["CM", "reference_node_rank", "node"]:
        raise AssertionError("Canonical membership schema changed.")
    membership.to_csv(NODE_DIR / "balanced_joint_cm_reference_node_sets_after_edge_threshold.csv", index=False)
    pd.DataFrame(membership_diagnostic).to_csv(
        NODE_DIR / "balanced_joint_cm_reference_node_sets_after_edge_threshold_diagnostic.csv",
        index=False,
    )
    edges = pd.DataFrame(
        edge_rows,
        columns=["context", "CM", "node_a", "node_b", "pearson_r", "edge_pass_r_ge_0.25"],
    )
    edges.to_csv(NODE_DIR / "status_specific_nodeplot_edges.csv", index=False)
    node_table = pd.DataFrame(membership_diagnostic).rename(columns={"node": "cell_subtype", "reference_node_rank": "rank"})
    node_table["status_context"] = "tumor"
    node_table[["CM", "cell_subtype", "loading", "rank", "cell_lineage", "status_context"]].to_csv(
        NODE_DIR / "tumor_network_nodes_from_H_df.csv", index=False
    )
    pd.DataFrame(cm_audit).to_csv(NODE_DIR / "nodeplot_panel_audit.csv", index=False)

    union_top10 = nodes_top10.sort_values(["CM", "rank"])["cell_subtype"].drop_duplicates().tolist()
    global_corr, global_p, global_q = corr_and_q(frequency[union_top10])
    global_corr.to_csv(NODE_DIR / "node_node_correlation_matrix.csv")
    global_p.to_csv(NODE_DIR / "node_node_correlation_p_matrix.csv")
    global_q.to_csv(NODE_DIR / "node_node_correlation_q_matrix.csv")
    global_corr.to_csv(NODE_DIR / "tumor_top10_node_correlation_matrix_no_edge_filter.csv")
    global_q.to_csv(NODE_DIR / "tumor_top10_node_correlation_q_matrix_no_edge_filter.csv")
    pd.DataFrame(
        [
            {"detected_mode": mode, "status_context": "tumor",
             "top_n_subtypes": TOP_N_SUBTYPES, "plot_top_n_subtypes": 12,
             "top_n_nodes": TOP_N_NODES, "edge_r_threshold": EDGE_R_THRESHOLD,
             "node_correlation_method": "Pearson",
             "edge_pass_rule": "pearson_r >= 0.25 in tumor context",
             "normal_like_context_available": False,
             "code_file": str(CODE_PATH), "seed": SEED}
        ]
    ).to_csv(NODE_DIR / "nodes_edges_parameters.csv", index=False)
    (NODE_DIR / "package_versions.txt").write_text(
        f"python={sys.version.split()[0]}\npandas={package_version('pandas')}\nscipy={package_version('scipy')}\ncode={CODE_PATH}\nseed={SEED}\n",
        encoding="utf-8",
    )
    completion = {
        "status": "completed", "detected_mode": mode,
        "n_cms": len(cm_order), "n_all_loading_nodes": len(nodes_all),
        "n_top10_diagnostic_nodes": len(nodes_top10),
        "n_reference_nodes_after_edge_threshold": len(membership),
        "n_status_specific_edge_rows": len(edges),
        "n_passing_tumor_edges": int(edges["edge_pass_r_ge_0.25"].sum()) if len(edges) else 0,
        "all_cm_panels_drawable": bool(pd.DataFrame(cm_audit)["panel_drawable"].all()),
        "normal_like_outputs_skipped": True,
        "elapsed_seconds": time.time() - started,
    }
    write_json(completion, completion_path)
    print(json.dumps(completion, indent=2), flush=True)


if __name__ == "__main__":
    main()
