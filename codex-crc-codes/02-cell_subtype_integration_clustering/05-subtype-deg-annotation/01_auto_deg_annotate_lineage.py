#!/usr/bin/env python3
"""Run subtype DEGs, agent-authorized DEG-based naming, and final annotation."""

from __future__ import annotations

import argparse
import gc
import platform
import random
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc

matplotlib.use("Agg")

SEED = 42
random.seed(SEED)
np.random.seed(SEED)

WORKFLOW_ROOT = Path("/mnt/disk18t/lr_xcy/riku/crc_val/epi-cm-core-workflow")
MODULE = "02-cell_subtype_integration_clustering"
STEP_SELECTED = "04-subtype-selected-clustering"
STEP_ANNOTATION = "05-subtype-deg-annotation"
CODE_FILE = Path(__file__)

LINEAGE_CONFIG = {
    "b_cells": {"abbr": "b", "prefix": "B", "selected_file": "adata_b_selected_clustered.h5ad"},
    "cycling_immune": {"abbr": "cycling", "prefix": "Cycling", "selected_file": "adata_cycling_selected_clustered.h5ad"},
    "endothelial": {"abbr": "endo", "prefix": "Endo", "selected_file": "adata_endo_selected_clustered.h5ad"},
    "epithelial": {"abbr": "epi", "prefix": "Epi", "selected_file": "adata_epi_selected_clustered.h5ad"},
    "mast": {"abbr": "mast", "prefix": "Mast", "selected_file": "adata_mast_selected_clustered.h5ad"},
    "myeloid": {"abbr": "mye", "prefix": "Mye", "selected_file": "adata_mye_selected_clustered.h5ad"},
    "plasma_cells": {"abbr": "plasma", "prefix": "Plasma", "selected_file": "adata_plasma_selected_clustered.h5ad"},
    "schwann": {"abbr": "schwann", "prefix": "Schwann", "selected_file": "adata_schwann_selected_clustered.h5ad"},
    "stromal": {"abbr": "stromal", "prefix": "Stromal", "selected_file": "adata_stromal_selected_clustered.h5ad"},
    "t_cells": {"abbr": "t", "prefix": "T", "selected_file": "adata_t_selected_clustered.h5ad"},
}

STATE_MARKERS = {
    "proliferative": {"MKI67", "TOP2A", "PCNA", "UBE2C", "STMN1", "HMGB2"},
    "stress_response": {"FOS", "JUN", "JUNB", "DUSP1", "HSPA1A", "HSPA1B", "HSPH1"},
    "interferon_response": {"ISG15", "IFIT1", "IFIT2", "IFIT3", "IFI6", "MX1", "OAS1"},
    "cytokine_chemokine": {"CCL2", "CCL3", "CCL4", "CCL5", "CXCL8", "CXCL9", "CXCL10"},
    "metabolic": {"APOE", "APOC1", "FABP4", "FABP5", "MT1G", "MT1X"},
    "hypoxia": {"VEGFA", "CA9", "ENO1", "LDHA", "NDRG1"},
    "antigen_presentation": {"HLA-DRA", "HLA-DRB1", "HLA-DPA1", "CD74", "B2M"},
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lineage", choices=sorted(LINEAGE_CONFIG), required=True)
    parser.add_argument("--allow-existing", action="store_true")
    return parser.parse_args()


def res_token(resolution: float) -> str:
    return str(float(resolution)).replace(".", "p")


def safe_name(value: str) -> str:
    return str(value).replace("/", "_").replace(" ", "_")


def export_rank_genes(adata: ad.AnnData, groupby: str, out_dir: Path, filename_token: str) -> dict[str, pd.DataFrame]:
    out_dir.mkdir(parents=True, exist_ok=True)
    sc.tl.rank_genes_groups(adata, groupby=groupby, method="t-test", use_raw=True)
    result = adata.uns["rank_genes_groups"]
    groups = [str(g) for g in result["names"].dtype.names]
    exported: dict[str, pd.DataFrame] = {}
    for group in groups:
        df = pd.DataFrame(
            {
                "gene": result["names"][group],
                "score": result["scores"][group],
                "logfoldchanges": result["logfoldchanges"][group],
                "pvals": result["pvals"][group],
                "pvals_adj": result["pvals_adj"][group],
            }
        )
        df.to_csv(out_dir / f"{safe_name(group)}_degs_{filename_token}.csv", index=False)
        exported[group] = df
    return exported


def eligible_positive_rows(df: pd.DataFrame) -> pd.DataFrame:
    pos = df[np.isfinite(df["logfoldchanges"]) & (df["logfoldchanges"] > 0)].copy()
    if pos.empty:
        pos = df.copy()
    pos["gene"] = pos["gene"].astype(str)
    return pos


def evidence_tuple(row: pd.Series) -> tuple[float, float, float, float]:
    p_adj = float(row["pvals_adj"]) if pd.notna(row["pvals_adj"]) else np.inf
    p_val = float(row["pvals"]) if pd.notna(row["pvals"]) else np.inf
    score = float(row["score"]) if pd.notna(row["score"]) else -np.inf
    lfc = float(row["logfoldchanges"]) if pd.notna(row["logfoldchanges"]) else -np.inf
    return p_adj, p_val, -score, -lfc


def infer_state(gene: str, top_genes: list[str], lineage: str) -> str:
    top_set = set(top_genes[:20])
    for state, markers in STATE_MARKERS.items():
        if gene in markers or top_set.intersection(markers):
            return state
    return f"{lineage}_marker_defined_state"


def build_mapping(raw_degs: dict[str, pd.DataFrame], prefix: str, lineage: str) -> pd.DataFrame:
    candidates: dict[str, list[pd.Series]] = {}
    for cluster, df in raw_degs.items():
        pos = eligible_positive_rows(df)
        candidates[cluster] = [row for _, row in pos.iterrows()]

    chosen: dict[str, pd.Series] = {}
    duplicate_notes: dict[str, str] = {}
    remaining = set(candidates)
    while remaining:
        first_choice: dict[str, str] = {}
        for cluster in remaining:
            rows = candidates[cluster]
            if not rows:
                raise RuntimeError(f"No DEG rows left for cluster {cluster}")
            first_choice[cluster] = str(rows[0]["gene"])
        gene_to_clusters: dict[str, list[str]] = {}
        for cluster, gene in first_choice.items():
            gene_to_clusters.setdefault(gene, []).append(cluster)
        for gene, clusters in gene_to_clusters.items():
            if len(clusters) == 1:
                cluster = clusters[0]
                chosen[cluster] = candidates[cluster].pop(0)
                remaining.remove(cluster)
                continue
            winner = sorted(clusters, key=lambda c: evidence_tuple(candidates[c][0]))[0]
            chosen[winner] = candidates[winner].pop(0)
            remaining.remove(winner)
            for loser in clusters:
                if loser == winner:
                    continue
                candidates[loser].pop(0)
                duplicate_notes[loser] = (
                    f"Top positive DEG {gene} duplicated; {winner} kept {prefix}_{gene} by stronger evidence, "
                    "this cluster used the next-ranked eligible positive DEG from its own DEG table."
                )

    rows = []
    for cluster in sorted(chosen, key=lambda x: (len(str(x)), str(x))):
        row = chosen[cluster]
        gene = str(row["gene"])
        df = raw_degs[cluster]
        top_genes = eligible_positive_rows(df)["gene"].astype(str).head(10).tolist()
        state = infer_state(gene, top_genes, lineage)
        rationale = (
            f"{gene} is the most significant eligible positive DEG for raw cluster {cluster} "
            "after ordering by adjusted P value, raw P value, Scanpy score, and positive logfoldchange."
        )
        if cluster in duplicate_notes:
            rationale += " " + duplicate_notes[cluster]
        rows.append(
            {
                "cluster": cluster,
                "cell_subtype": f"{prefix}_{gene}",
                "functional_state": state,
                "selected_prefix": prefix,
                "most_significant_deg_gene": gene,
                "gene_selection_rationale": rationale,
                "top_supporting_markers": ";".join(top_genes[:10]),
                "marker_evidence": (
                    f"score={row['score']}; logfoldchanges={row['logfoldchanges']}; "
                    f"pvals={row['pvals']}; pvals_adj={row['pvals_adj']}"
                ),
                "annotation_confidence": "agent_authorized_deg_based",
                "annotation_note": "Automatic DEG-based subtype name following updated project skill.",
            }
        )
    mapping = pd.DataFrame(rows)
    if mapping["cell_subtype"].duplicated().any():
        dup = mapping.loc[mapping["cell_subtype"].duplicated(), "cell_subtype"].tolist()
        raise RuntimeError(f"Duplicate cell_subtype labels remain: {dup}")
    return mapping


def main() -> None:
    args = parse_args()
    cfg = LINEAGE_CONFIG[args.lineage]
    selected_params = (
        WORKFLOW_ROOT
        / "tables"
        / MODULE
        / STEP_SELECTED
        / args.lineage
        / "selected_clustering_parameters.csv"
    )
    params = pd.read_csv(selected_params).iloc[0]
    cluster_key = str(params["cluster_key"])
    n_pcs = int(params["n_pcs"])
    n_neighbors = int(params["n_neighbors"])
    resolution = float(params["resolution"])
    res_tok = res_token(resolution)
    param_token = f"pcs{n_pcs}_nn{n_neighbors}_res{res_tok}"

    selected_h5ad = (
        WORKFLOW_ROOT
        / "h5ad"
        / MODULE
        / STEP_SELECTED
        / args.lineage
        / cfg["selected_file"]
    )
    output_h5ad = (
        WORKFLOW_ROOT
        / "h5ad"
        / MODULE
        / STEP_ANNOTATION
        / args.lineage
        / f"adata_{cfg['abbr']}.h5ad"
    )
    table_dir = WORKFLOW_ROOT / "tables" / MODULE / STEP_ANNOTATION / args.lineage
    figure_dir = WORKFLOW_ROOT / "figures" / MODULE / STEP_ANNOTATION / args.lineage
    raw_deg_dir = table_dir / f"degs_{cluster_key}_{param_token}"
    subtype_deg_dir = table_dir / f"degs_cell_subtype_{param_token}"
    candidate_table = table_dir / "candidate_subtype_annotation_table.csv"
    final_mapping = table_dir / "final_subtype_mapping.csv"
    subtype_counts = table_dir / f"{args.lineage}_subtype_counts.csv"
    final_umap = figure_dir / "umap_cell_subtype.pdf"
    marker_dotplot = figure_dir / "dotplot_cell_subtype_markers.pdf"

    outputs = [output_h5ad, candidate_table, final_mapping, subtype_counts, final_umap, marker_dotplot]
    existing = [str(path) for path in outputs if path.exists()]
    if existing and not args.allow_existing:
        raise FileExistsError("Subtype annotation output already exists; refusing to overwrite:\n" + "\n".join(existing))
    if args.allow_existing and output_h5ad.exists() and final_mapping.exists():
        print(f"{args.lineage}: final annotation already exists")
        return

    output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    sc.settings.autoshow = False
    sc.settings.figdir = str(figure_dir)
    sc.settings.set_figure_params(figsize=(3, 3), dpi=150)

    adata = ad.read_h5ad(selected_h5ad)
    if adata.raw is None:
        raise RuntimeError(f"Missing raw in {selected_h5ad}")
    if cluster_key not in adata.obs.columns:
        raise KeyError(f"Missing selected raw cluster key {cluster_key}")

    adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")
    raw_degs = export_rank_genes(adata, cluster_key, raw_deg_dir, f"{cluster_key}_{param_token}")
    mapping = build_mapping(raw_degs, cfg["prefix"], args.lineage)
    mapping["lineage"] = args.lineage
    mapping["cluster_key"] = cluster_key
    mapping["n_pcs"] = n_pcs
    mapping["n_neighbors"] = n_neighbors
    mapping["resolution"] = resolution
    mapping["agent_authorization"] = "User requested continuing to complete subtype naming on 2026-07-11."
    mapping.to_csv(candidate_table, index=False)
    mapping.to_csv(final_mapping, index=False)

    subtype_map = mapping.set_index("cluster")["cell_subtype"].to_dict()
    state_map = mapping.set_index("cluster")["functional_state"].to_dict()
    rationale_map = mapping.set_index("cluster")["gene_selection_rationale"].to_dict()
    evidence_map = mapping.set_index("cluster")["marker_evidence"].to_dict()
    adata.obs["cell_subtype"] = adata.obs[cluster_key].astype(str).map(subtype_map).astype("category")
    adata.obs["functional_state"] = adata.obs[cluster_key].astype(str).map(state_map).astype("category")
    adata.obs["annotation_source"] = f"{args.lineage}:{cluster_key}:{param_token}"
    adata.obs["annotation_confidence"] = "agent_authorized_deg_based"
    adata.obs["gene_selection_rationale"] = adata.obs[cluster_key].astype(str).map(rationale_map)
    adata.obs["marker_evidence"] = adata.obs[cluster_key].astype(str).map(evidence_map)
    adata.obs["cell_type_previous"] = adata.obs["cell_type"].astype(str).values if "cell_type" in adata.obs else ""
    adata.obs["cell_type"] = cfg["prefix"]

    counts = adata.obs["cell_subtype"].astype(str).value_counts().rename_axis("cell_subtype").reset_index(name="n_cells")
    counts.to_csv(subtype_counts, index=False)

    subtype_degs = export_rank_genes(adata, "cell_subtype", subtype_deg_dir, f"cell_subtype_{param_token}")
    del subtype_degs

    marker_genes = mapping["most_significant_deg_gene"].astype(str).drop_duplicates().tolist()
    sc.pl.umap(adata, color="cell_subtype", save="_cell_subtype.pdf", show=False)
    sc.pl.dotplot(
        adata,
        marker_genes,
        groupby="cell_subtype",
        use_raw=True,
        standard_scale="var",
        save="cell_subtype_markers.pdf",
        show=False,
    )

    adata.write_h5ad(output_h5ad)

    summary = pd.DataFrame(
        [
            {
                "lineage": args.lineage,
                "selected_h5ad": str(selected_h5ad),
                "output_h5ad": str(output_h5ad),
                "cluster_key": cluster_key,
                "n_clusters": int(mapping.shape[0]),
                "n_cell_subtypes": int(mapping["cell_subtype"].nunique()),
                "raw_deg_dir": str(raw_deg_dir),
                "post_annotation_deg_dir": str(subtype_deg_dir),
                "candidate_annotation_table": str(candidate_table),
                "final_mapping_csv": str(final_mapping),
                "final_umap_pdf": str(final_umap),
                "marker_dotplot_pdf": str(marker_dotplot),
                "agent_authorized_final_labeling": True,
                "random_seed": SEED,
                "code_file": str(CODE_FILE),
            }
        ]
    )
    summary.to_csv(table_dir / "subtype_annotation_summary.csv", index=False)

    with (table_dir / "package_versions.txt").open("w") as fh:
        fh.write(f"python: {platform.python_version()}\n")
        fh.write(f"anndata: {ad.__version__}\n")
        fh.write(f"scanpy: {sc.__version__}\n")
        fh.write(f"numpy: {np.__version__}\n")
        fh.write(f"pandas: {pd.__version__}\n")
        fh.write(f"code_file: {CODE_FILE}\n")

    with (table_dir / "readme.txt").open("w") as fh:
        fh.write(f"{args.lineage} subtype DEG annotation completed.\n")
        fh.write(f"Selected clustered input: {selected_h5ad}\n")
        fh.write(f"Output annotated h5ad: {output_h5ad}\n")
        fh.write(f"Raw cluster key: {cluster_key}\n")
        fh.write("Subtype names were assigned from the most significant eligible positive DEG per raw cluster.\n")
        fh.write("Duplicate top positive DEG labels were resolved using the next-ranked eligible DEG as required by the updated skill.\n")

    print(summary.to_string(index=False))
    del adata
    gc.collect()


if __name__ == "__main__":
    main()
