#!/usr/bin/env python3
"""Build a marker-reviewed candidate subtype table for one BRCA lineage."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path("/mnt/disk18t/lr_xcy/riku/brca_cm-epi_codex")
WORKFLOW = ROOT / "epi-cm-core-workflow"
BLOCK = "02-cell_subtype_integration_clustering"
REVIEW_DEPTH = 50

LINEAGE_CONFIG = {
    "epithelial": {"prefix": "Epi", "label": "Epithelial Cells"},
    "t_cells": {"prefix": "T", "label": "T Cells"},
    "myeloid": {"prefix": "Mye", "label": "Myeloid Cells"},
    "b_cells": {"prefix": "B", "label": "B Cells"},
    "plasma": {"prefix": "Plasma", "label": "Plasma Cells"},
    "endothelial": {"prefix": "Endo", "label": "Endothelial Cells"},
    "stromal": {"prefix": "Stromal", "label": "Stromal Cells"},
    "perivascular": {"prefix": "PVL", "label": "Perivascular Cells"},
}

REVIEW = {
    "epithelial": {
        "0": ("interferon-responsive epithelial", "epithelial", "high"),
        "1": ("luminal-secretory oxidative epithelial", "epithelial", "moderate"),
        "2": ("luminal immediate-early epithelial", "epithelial", "high"),
        "3": ("proliferative cycling epithelial", "epithelial", "high"),
        "4": ("stress/low-information epithelial", "epithelial", "low"),
        "5": ("immune-like epithelial", "immune-like", "high"),
        "6": ("basal/progenitor-like epithelial", "epithelial", "high"),
        "7": ("immediate-early activated epithelial", "epithelial", "high"),
        "8": ("mesenchymal/neural-like epithelial", "epithelial", "moderate"),
        "9": ("XBP1-high secretory luminal epithelial", "epithelial", "low"),
    },
    "t_cells": {
        "0": ("CCR7/IL7R CD4 naive-like T", "T", "moderate"),
        "1": ("activated HLA-II cytotoxic CD8 T", "T", "high"),
        "2": ("IL7R CD4 memory-like T", "T", "high"),
        "3": ("interferon-responsive CD8 T", "T", "high"),
        "4": ("ZNF683 tissue-resident effector CD8 T", "T", "high"),
        "5": ("immediate-early activated CD4 T", "T", "high"),
        "6": ("NK-like cytotoxic T", "NK-like", "high"),
        "7": ("CXCL13 exhausted/helper-like T", "T", "high"),
        "8": ("FOXP3/CTLA4 regulatory-exhausted T", "T", "high"),
        "9": ("NKT-like cytotoxic T", "NKT-like", "high"),
    },
    "myeloid": {
        "0": ("SPP1/GPNMB lipid-associated macrophage", "myeloid", "high"),
        "1": ("proliferative cycling myeloid", "myeloid", "high"),
        "2": ("CXCL10 interferon macrophage", "myeloid", "high"),
        "3": ("epithelial-like inflammatory myeloid", "epithelial-like", "moderate"),
        "4": ("FCN1/S100A8 classical monocyte", "myeloid", "high"),
        "5": ("APOE/C1Q lipid-associated macrophage", "myeloid", "high"),
        "6": ("C1QC resident-like macrophage", "myeloid", "high"),
        "7": ("cycling/vascular-like myeloid", "off-lineage-like", "low"),
        "8": ("LAMP3/CCR7 dendritic cell", "myeloid", "high"),
    },
    "b_cells": {
        "0": ("interferon-responsive memory B", "B", "high"),
        "1": ("NR4A/CD83 activated memory B", "B", "high"),
        "2": ("CD27/IGHA memory B", "B", "high"),
        "3": ("ribosomal low-information B", "B", "low"),
        "4": ("DC/T-like B", "off-lineage-like", "low"),
        "5": ("TCL1A/IGHD naive B", "B", "high"),
        "6": ("GZMB/LILRA4 plasmacytoid-DC-like B", "pDC-like", "high"),
    },
    "plasma": {
        "0": ("kappa-light-chain antibody-secreting", "plasma", "high"),
        "1": ("lambda/ECM-associated antibody-secreting", "plasma", "moderate"),
        "2": ("antigen-presentation/translational plasma", "plasma", "moderate"),
        "3": ("immune-like activated plasma", "immune-like", "low"),
        "4": ("stress-response plasma", "plasma", "moderate"),
        "5": ("lambda-light-chain antibody-secreting", "plasma", "high"),
        "6": ("stress/epithelial-like plasma", "epithelial-like", "low"),
    },
    "endothelial": {
        "0": ("ACKR1/CCL14 venous HEV-like endothelial", "endothelial", "high"),
        "1": ("ACKR1 HLA-II activated venous endothelial", "endothelial", "high"),
        "2": ("VWF vascular endothelial", "endothelial", "high"),
        "3": ("RGS5/PDGFRB perivascular-like endothelial", "perivascular-like", "high"),
        "4": ("interferon/immune-like endothelial", "immune-like", "moderate"),
        "5": ("PLVAP/KDR angiogenic capillary endothelial", "endothelial", "high"),
        "6": ("FABP4/CA4 capillary endothelial", "endothelial", "high"),
        "7": ("CXCL12/GJA4 arterial-like endothelial", "endothelial", "high"),
        "8": ("inflamed perivascular-like endothelial", "perivascular-like", "moderate"),
        "9": ("CCL21/PROX1 lymphatic endothelial", "endothelial", "high"),
    },
    "stromal": {
        "0": ("COL1A1/POSTN matrix myCAF", "stromal", "high"),
        "1": ("HMGA proliferative epithelial-like CAF", "epithelial-like", "moderate"),
        "2": ("immediate-early myCAF", "stromal", "high"),
        "3": ("COL3A1/THY1 matrix myCAF", "stromal", "high"),
        "4": ("APOD/CXCL12 iCAF", "stromal", "high"),
        "5": ("MGP/CFD matrix-adipogenic iCAF", "stromal", "high"),
        "6": ("SOD2/CCL2 inflammatory CAF", "stromal", "high"),
        "7": ("immediate-early CXCL12 iCAF", "stromal", "high"),
    },
    "perivascular": {
        "0": ("CALD1/ACTA2 differentiated contractile PVL", "perivascular", "high"),
        "1": ("immediate-early differentiated PVL", "perivascular", "high"),
        "2": ("oxidative antigen-presenting PVL", "perivascular", "moderate"),
        "3": ("C1S/CCL19 inflammatory immature PVL", "perivascular", "high"),
        "4": ("ADIRF/MYH11 smooth-muscle-like PVL", "perivascular", "high"),
        "5": ("TAGLN/MYH11 contractile PVL", "perivascular", "high"),
        "6": ("THY1/RGS5 matrix immature PVL", "perivascular", "high"),
    },
}

LINEAGE_SLUG = "epithelial"
CODE_PATH = Path(__file__).resolve()


def configure_lineage(lineage_slug: str, code_path: Path | None = None) -> None:
    global LINEAGE_SLUG, CODE_PATH
    if lineage_slug not in LINEAGE_CONFIG:
        raise ValueError(f"Unsupported lineage slug: {lineage_slug}")
    LINEAGE_SLUG = lineage_slug
    if code_path is not None:
        CODE_PATH = code_path


def evidence_tuple(row: pd.Series) -> tuple[float, float, float, float]:
    return (
        float(row["pvals_adj"]) if np.isfinite(row["pvals_adj"]) else np.inf,
        float(row["pvals"]) if np.isfinite(row["pvals"]) else np.inf,
        -float(row["score"]) if np.isfinite(row["score"]) else np.inf,
        -abs(float(row["logfoldchanges"]))
        if np.isfinite(row["logfoldchanges"])
        else np.inf,
    )


def main() -> None:
    table_dir = WORKFLOW / "tables" / BLOCK / "05-subtype-deg-annotation" / LINEAGE_SLUG
    output = table_dir / "candidate_subtype_annotation.csv"
    parameters = table_dir / "candidate_annotation_review_parameters.csv"
    if output.exists() or parameters.exists():
        if output.exists() and parameters.exists():
            candidate = pd.read_csv(output, dtype={"cluster": str})
            if (
                len(candidate) == len(REVIEW[LINEAGE_SLUG])
                and candidate["cluster"].is_unique
                and candidate["proposed_cell_subtype"].is_unique
                and candidate["needs_user_confirmation"].all()
            ):
                print(json.dumps({"lineage": LINEAGE_SLUG, "status": "valid_existing_candidate_reused", "rows": len(candidate)}, indent=2))
                return
        raise FileExistsError("Partial or invalid candidate annotation outputs already exist.")

    audit = pd.read_csv(table_dir / "raw_cluster_deg_audit.csv", dtype={"raw_cluster": str})
    if set(audit["raw_cluster"]) != set(REVIEW[LINEAGE_SLUG]):
        raise ValueError("Curated review rows do not exactly cover raw clusters.")
    minor = pd.read_csv(
        table_dir / "raw_cluster_by_published_minor_fractions.csv", index_col=0
    )
    minor.index = minor.index.astype(str)

    frames: dict[str, pd.DataFrame] = {}
    top_gene: dict[str, str] = {}
    for row in audit.itertuples(index=False):
        frame = pd.read_csv(row.deg_csv)
        positive = frame[frame["logfoldchanges"] > 0].copy()
        if positive.empty:
            raise ValueError(f"No positive DEGs for cluster {row.raw_cluster}")
        frames[str(row.raw_cluster)] = positive
        top_gene[str(row.raw_cluster)] = str(positive.iloc[0]["gene"])

    gene_to_clusters: dict[str, list[str]] = {}
    for cluster, gene in top_gene.items():
        gene_to_clusters.setdefault(gene, []).append(cluster)
    selected_gene: dict[str, str] = {}
    duplicate_note: dict[str, str] = {}
    reserved: set[str] = set()
    for gene, clusters in gene_to_clusters.items():
        if len(clusters) == 1:
            selected_gene[clusters[0]] = gene
            reserved.add(gene)
            duplicate_note[clusters[0]] = ""
            continue
        ranked = sorted(
            clusters,
            key=lambda cluster: evidence_tuple(frames[cluster].iloc[0]),
        )
        winner = ranked[0]
        selected_gene[winner] = gene
        reserved.add(gene)
        duplicate_note[winner] = (
            f"Top positive DEG {gene} was shared by clusters {','.join(clusters)}; "
            f"cluster {winner} retained it based on stronger adjusted-P/raw-P/score evidence."
        )
        for loser in ranked[1:]:
            duplicate_note[loser] = (
                f"Top positive DEG {gene} was assigned to cluster {winner} based on stronger "
                "adjusted-P/raw-P/score evidence; next-ranked eligible unique positive DEG used."
            )
    for cluster in sorted(frames, key=lambda value: int(value)):
        if cluster in selected_gene:
            continue
        for gene in frames[cluster]["gene"].astype(str):
            if gene not in reserved:
                selected_gene[cluster] = gene
                reserved.add(gene)
                break
        if cluster not in selected_gene:
            raise ValueError(f"Cannot resolve a unique positive DEG for cluster {cluster}")

    prefix = LINEAGE_CONFIG[LINEAGE_SLUG]["prefix"]
    total_cells = int(audit["n_cells"].sum())
    rows: list[dict[str, object]] = []
    for row in audit.sort_values("raw_cluster", key=lambda x: x.astype(int)).itertuples(index=False):
        cluster = str(row.raw_cluster)
        gene = selected_gene[cluster]
        chosen = frames[cluster].loc[frames[cluster]["gene"].astype(str).eq(gene)].iloc[0]
        state, marker_lineage, confidence = REVIEW[LINEAGE_SLUG][cluster]
        supporting = frames[cluster]["gene"].astype(str).head(15).tolist()
        if cluster in minor.index and len(minor.columns):
            author = minor.loc[cluster].sort_values(ascending=False)
            author_label = str(author.index[0])
            author_fraction = float(author.iloc[0])
        else:
            author_label = "unavailable"
            author_fraction = np.nan
        if duplicate_note[cluster]:
            rationale = (
                f"{gene} is the next-ranked eligible unique positive DEG in the full t-test table. "
                + duplicate_note[cluster]
            )
        else:
            rationale = (
                f"{gene} is the first positive-logfoldchange gene in the full t-test DEG table "
                "sorted by Scanpy rank_genes_groups evidence."
            )
        prefix_rationale = (
            f"Default current-lineage prefix {prefix} derived from "
            f"{LINEAGE_CONFIG[LINEAGE_SLUG]['label']}; marker evidence does not move cells "
            "between lineage h5ad objects."
        )
        rows.append(
            {
                "cluster": cluster,
                "n_cells": int(row.n_cells),
                "fraction_of_lineage": int(row.n_cells) / total_cells,
                "proposed_cell_subtype": f"{prefix}_{gene}",
                "functional_state": state,
                "selected_prefix": prefix,
                "most_significant_deg_gene": gene,
                "top_positive_deg_gene_before_duplicate_resolution": top_gene[cluster],
                "gene_score": float(chosen["score"]),
                "gene_logfoldchanges": float(chosen["logfoldchanges"]),
                "gene_pvals": float(chosen["pvals"]),
                "gene_pvals_adj": float(chosen["pvals_adj"]),
                "gene_selection_rationale": rationale,
                "top supporting markers": ";".join(supporting),
                "marker-supported lineage": marker_lineage,
                "prefix rationale": prefix_rationale,
                "published_minor_dominant_label_validation_only": author_label,
                "published_minor_dominant_fraction_validation_only": author_fraction,
                "review_depth": REVIEW_DEPTH,
                "user_modified_marker": "",
                "confidence": confidence,
                "needs_user_confirmation": True,
            }
        )
    candidate = pd.DataFrame(rows)
    if not candidate["cluster"].is_unique or not candidate["proposed_cell_subtype"].is_unique:
        raise ValueError("Candidate mapping is not one-to-one and unique.")
    candidate.to_csv(output, index=False)
    pd.DataFrame(
        [
            {
                "lineage": LINEAGE_SLUG,
                "raw_cluster_deg_audit": str(table_dir / "raw_cluster_deg_audit.csv"),
                "candidate_annotation_table": str(output),
                "review_depth": REVIEW_DEPTH,
                "review_scope": "same first 50 rows for every raw cluster",
                "published_labels_role": "independent_validation_only",
                "prefix_source": "conservative default from broad lineage",
                "duplicate_gene_rule": "padj then pval then score then absolute positive logfoldchange",
                "agent_led_final_labeling_authorized": False,
                "final_mapping_written": False,
                "code_file": str(CODE_PATH),
            }
        ]
    ).to_csv(parameters, index=False)
    print(
        json.dumps(
            {
                "lineage": LINEAGE_SLUG,
                "status": "candidate_annotation_completed",
                "rows": len(candidate),
                "unique_proposed_cell_subtypes": int(candidate["proposed_cell_subtype"].nunique()),
                "needs_user_confirmation": True,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
