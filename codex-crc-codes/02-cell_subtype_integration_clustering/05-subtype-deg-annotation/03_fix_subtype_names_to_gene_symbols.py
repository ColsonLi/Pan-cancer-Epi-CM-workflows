#!/usr/bin/env python3
"""Replace Ensembl-ID subtype suffixes with gene symbols and refresh final outputs."""

from __future__ import annotations

import gc
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
GTF_CSV = Path(
    "/mnt/disk18t/lr_xcy/riku/external_atlas_for_cm-epi_validation/CRC_S-BIAD2208/"
    "input_datasets_extracted/load_datasets/harmonize_datasets/artifacts/adata_var_gtf.csv"
)
MODULE = "02-cell_subtype_integration_clustering"
STEP = "05-subtype-deg-annotation"
CODE_FILE = Path(__file__)

LINEAGE_CONFIG = {
    "b_cells": {"abbr": "b", "prefix": "B"},
    "cycling_immune": {"abbr": "cycling", "prefix": "Cycling"},
    "endothelial": {"abbr": "endo", "prefix": "Endo"},
    "epithelial": {"abbr": "epi", "prefix": "Epi"},
    "mast": {"abbr": "mast", "prefix": "Mast"},
    "myeloid": {"abbr": "mye", "prefix": "Mye"},
    "plasma_cells": {"abbr": "plasma", "prefix": "Plasma"},
    "schwann": {"abbr": "schwann", "prefix": "Schwann"},
    "stromal": {"abbr": "stromal", "prefix": "Stromal"},
    "t_cells": {"abbr": "t", "prefix": "T"},
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


def safe_name(value: str) -> str:
    return str(value).replace("/", "_").replace(" ", "_")


def load_gene_map() -> dict[str, str]:
    gtf = pd.read_csv(GTF_CSV, usecols=["ensembl", "GeneSymbol"])
    gtf = gtf.dropna().drop_duplicates("ensembl")
    return dict(zip(gtf["ensembl"].astype(str), gtf["GeneSymbol"].astype(str)))


def symbol_for(gene_id: str, gene_map: dict[str, str]) -> str:
    gene_id = str(gene_id)
    return gene_map.get(gene_id.split(".")[0], gene_id)


def evidence_tuple(row: pd.Series) -> tuple[float, float, float, float]:
    p_adj = float(row["pvals_adj"]) if pd.notna(row["pvals_adj"]) else np.inf
    p_val = float(row["pvals"]) if pd.notna(row["pvals"]) else np.inf
    score = float(row["score"]) if pd.notna(row["score"]) else -np.inf
    lfc = float(row["logfoldchanges"]) if pd.notna(row["logfoldchanges"]) else -np.inf
    return p_adj, p_val, -score, -lfc


def eligible_positive_rows(df: pd.DataFrame, gene_map: dict[str, str]) -> list[pd.Series]:
    pos = df[np.isfinite(df["logfoldchanges"]) & (df["logfoldchanges"] > 0)].copy()
    if pos.empty:
        pos = df.copy()
    pos["gene"] = pos["gene"].astype(str)
    pos["gene_symbol"] = pos["gene"].map(lambda g: symbol_for(g, gene_map))
    return [row for _, row in pos.iterrows()]


def infer_state(gene_symbol: str, top_symbols: list[str], lineage: str) -> str:
    top_set = set(top_symbols[:20])
    for state, markers in STATE_MARKERS.items():
        if gene_symbol in markers or top_set.intersection(markers):
            return state
    return f"{lineage}_marker_defined_state"


def build_symbol_mapping(raw_degs: dict[str, pd.DataFrame], prefix: str, lineage: str, gene_map: dict[str, str]) -> pd.DataFrame:
    candidates = {cluster: eligible_positive_rows(df, gene_map) for cluster, df in raw_degs.items()}
    chosen: dict[str, pd.Series] = {}
    duplicate_notes: dict[str, str] = {}
    assigned_symbols: set[str] = set()
    remaining = set(candidates)
    while remaining:
        current: dict[str, pd.Series] = {}
        for cluster in list(remaining):
            while candidates[cluster] and str(candidates[cluster][0]["gene_symbol"]) in assigned_symbols:
                skipped = str(candidates[cluster].pop(0)["gene_symbol"])
                duplicate_notes[cluster] = (
                    f"Gene symbol {skipped} was already assigned to a stronger cluster; "
                    "this cluster used the next-ranked eligible positive DEG from its own DEG table."
                )
            if not candidates[cluster]:
                raise RuntimeError(f"No eligible unique gene symbols left for {lineage} cluster {cluster}")
            current[cluster] = candidates[cluster][0]
        symbol_to_clusters: dict[str, list[str]] = {}
        for cluster, row in current.items():
            symbol_to_clusters.setdefault(str(row["gene_symbol"]), []).append(cluster)
        progressed = False
        for symbol, clusters in symbol_to_clusters.items():
            if len(clusters) == 1:
                cluster = clusters[0]
                chosen[cluster] = candidates[cluster].pop(0)
                assigned_symbols.add(symbol)
                remaining.remove(cluster)
                progressed = True
                continue
            winner = sorted(clusters, key=lambda c: evidence_tuple(candidates[c][0]))[0]
            chosen[winner] = candidates[winner].pop(0)
            assigned_symbols.add(symbol)
            remaining.remove(winner)
            progressed = True
            for loser in clusters:
                if loser == winner:
                    continue
                candidates[loser].pop(0)
                duplicate_notes[loser] = (
                    f"Top positive gene symbol {symbol} duplicated; {winner} kept {prefix}_{symbol} by stronger evidence, "
                    "this cluster used the next-ranked eligible positive DEG from its own DEG table."
                )
        if not progressed:
            raise RuntimeError(f"Could not resolve subtype gene symbols for {lineage}")

    rows = []
    for cluster in sorted(chosen, key=lambda x: (len(str(x)), str(x))):
        row = chosen[cluster]
        gene_id = str(row["gene"])
        gene_symbol = str(row["gene_symbol"])
        top = eligible_positive_rows(raw_degs[cluster], gene_map)
        top_symbols = []
        for top_row in top:
            symbol = str(top_row["gene_symbol"])
            if symbol not in top_symbols:
                top_symbols.append(symbol)
            if len(top_symbols) >= 10:
                break
        state = infer_state(gene_symbol, top_symbols, lineage)
        rationale = (
            f"{gene_symbol} ({gene_id}) is the most significant eligible positive DEG for raw cluster {cluster} "
            "after mapping Ensembl IDs to GeneSymbol using adata_var_gtf.csv and ordering by adjusted P value, "
            "raw P value, Scanpy score, and positive logfoldchange."
        )
        if cluster in duplicate_notes:
            rationale += " " + duplicate_notes[cluster]
        rows.append(
            {
                "cluster": cluster,
                "cell_subtype": f"{prefix}_{gene_symbol}",
                "functional_state": state,
                "selected_prefix": prefix,
                "most_significant_deg_gene": gene_symbol,
                "most_significant_deg_gene_id": gene_id,
                "gene_selection_rationale": rationale,
                "top_supporting_markers": ";".join(top_symbols),
                "marker_evidence": (
                    f"score={row['score']}; logfoldchanges={row['logfoldchanges']}; "
                    f"pvals={row['pvals']}; pvals_adj={row['pvals_adj']}"
                ),
                "annotation_confidence": "agent_authorized_deg_based",
                "annotation_note": "Corrected subtype suffix from Ensembl ID to GeneSymbol following updated project skill.",
            }
        )
    mapping = pd.DataFrame(rows)
    if mapping["cell_subtype"].duplicated().any():
        raise RuntimeError(f"Duplicate subtype labels remain for {lineage}")
    return mapping


def export_rank_genes(adata: ad.AnnData, groupby: str, out_dir: Path, filename_token: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for old in out_dir.glob("*.csv"):
        old.unlink()
    sc.tl.rank_genes_groups(adata, groupby=groupby, method="t-test", use_raw=True)
    result = adata.uns["rank_genes_groups"]
    for group in result["names"].dtype.names:
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


def main() -> None:
    gene_map = load_gene_map()
    for lineage, cfg in LINEAGE_CONFIG.items():
        table_dir = WORKFLOW_ROOT / "tables" / MODULE / STEP / lineage
        figure_dir = WORKFLOW_ROOT / "figures" / MODULE / STEP / lineage
        h5ad_path = WORKFLOW_ROOT / "h5ad" / MODULE / STEP / lineage / f"adata_{cfg['abbr']}.h5ad"
        summary_path = table_dir / "subtype_annotation_summary.csv"
        summary = pd.read_csv(summary_path).iloc[0]
        raw_deg_dir = Path(summary["raw_deg_dir"])
        post_deg_dir = Path(summary["post_annotation_deg_dir"])
        cluster_key = str(summary["cluster_key"])
        raw_degs = {}
        for csv_path in sorted(raw_deg_dir.glob(f"*_degs_{cluster_key}_*.csv")):
            cluster = csv_path.name.split("_degs_")[0]
            raw_degs[cluster] = pd.read_csv(csv_path)
        mapping = build_symbol_mapping(raw_degs, cfg["prefix"], lineage, gene_map)
        mapping["lineage"] = lineage
        mapping["cluster_key"] = cluster_key
        mapping["agent_authorization"] = "User requested deleting/redoing subtype naming; corrected to GeneSymbol on 2026-07-11."
        mapping.to_csv(table_dir / "candidate_subtype_annotation_table.csv", index=False)
        mapping.to_csv(table_dir / "final_subtype_mapping.csv", index=False)

        adata = ad.read_h5ad(h5ad_path)
        subtype_map = mapping.set_index("cluster")["cell_subtype"].to_dict()
        state_map = mapping.set_index("cluster")["functional_state"].to_dict()
        rationale_map = mapping.set_index("cluster")["gene_selection_rationale"].to_dict()
        evidence_map = mapping.set_index("cluster")["marker_evidence"].to_dict()
        adata.obs["cell_subtype_previous"] = adata.obs["cell_subtype"].astype(str).values
        adata.obs["cell_subtype"] = adata.obs[cluster_key].astype(str).map(subtype_map).astype("category")
        adata.obs["functional_state"] = adata.obs[cluster_key].astype(str).map(state_map).astype("category")
        adata.obs["gene_selection_rationale"] = adata.obs[cluster_key].astype(str).map(rationale_map)
        adata.obs["marker_evidence"] = adata.obs[cluster_key].astype(str).map(evidence_map)
        adata.obs["annotation_note"] = "Subtype labels corrected from Ensembl IDs to GeneSymbol."

        counts = adata.obs["cell_subtype"].astype(str).value_counts().rename_axis("cell_subtype").reset_index(name="n_cells")
        counts.to_csv(table_dir / f"{lineage}_subtype_counts.csv", index=False)
        export_rank_genes(adata, "cell_subtype", post_deg_dir, post_deg_dir.name.replace("degs_", ""))

        sc.settings.autoshow = False
        sc.settings.figdir = str(figure_dir)
        sc.settings.set_figure_params(figsize=(3, 3), dpi=150)
        markers = mapping["most_significant_deg_gene_id"].astype(str).drop_duplicates().tolist()
        sc.pl.umap(adata, color="cell_subtype", save="_cell_subtype.pdf", show=False)
        sc.pl.dotplot(adata, markers, groupby="cell_subtype", use_raw=True, standard_scale="var", save="cell_subtype_markers.pdf", show=False)
        adata.write_h5ad(h5ad_path)

        summary_df = pd.DataFrame([summary.to_dict()])
        summary_df["gene_symbol_correction"] = True
        summary_df["gene_symbol_source"] = str(GTF_CSV)
        summary_df.to_csv(summary_path, index=False)
        with (table_dir / "gene_symbol_correction_readme.txt").open("w") as fh:
            fh.write("Subtype names were corrected from Ensembl IDs to GeneSymbol using adata_var_gtf.csv.\n")
            fh.write("The old Ensembl-based subtype names are preserved in obs['cell_subtype_previous'].\n")
            fh.write("Post-annotation DEG CSVs and final UMAP/dotplot were regenerated after correction.\n")
        print(lineage, "corrected", mapping["cell_subtype"].tolist())
        del adata
        gc.collect()


if __name__ == "__main__":
    main()
