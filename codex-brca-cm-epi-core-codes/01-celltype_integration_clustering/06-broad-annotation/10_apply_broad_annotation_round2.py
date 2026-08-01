#!/usr/bin/env python3
"""Apply broad labels, export round-2 DEGs, and save the annotated BRCA atlas."""

from __future__ import annotations

import importlib.metadata
import json
import re
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


ROOT = Path("/mnt/disk18t/lr_xcy/riku/brca_cm-epi_codex")
WORKFLOW = ROOT / "epi-cm-core-workflow"
BLOCK = "01-celltype_integration_clustering"
INPUT_H5AD = (
    WORKFLOW
    / "h5ad"
    / BLOCK
    / "05-clustering-parameter-search"
    / "selected"
    / "adata_inte.h5ad"
)
OUT_H5AD = WORKFLOW / "h5ad" / BLOCK / "06-broad-annotation" / "adata_anno.h5ad"
TABLE_DIR = WORKFLOW / "tables" / BLOCK / "06-broad-annotation"
FIGURE_DIR = WORKFLOW / "figures" / BLOCK / "06-broad-annotation"
MAPPING = TABLE_DIR / "broad_annotation_mapping.csv"
ROUND1_DEG_DIR = TABLE_DIR / "degs_leiden_res0p8_pcs20_nn30_res0p8"
ROUND2_DEG_DIR = (
    TABLE_DIR / "degs_leiden_coarse_pcs20_nn30_res0p8_myo_merged"
)
CODE_PATH = (
    WORKFLOW
    / "codes"
    / BLOCK
    / "06-broad-annotation"
    / "10_apply_broad_annotation_round2.py"
)

RAW_KEY = "leiden_res0p8"
COARSE_KEY = "leiden_coarse"
METHOD = "t-test"
SEED = 42
CATEGORY_ORDER = [
    "Epithelial Cells",
    "T Cells",
    "Myeloid Cells",
    "B Cells",
    "Plasma Cells",
    "Endothelial Cells",
    "Stromal Cells",
    "Perivascular Cells",
]
MARKERS = {
    "Epithelial Cells": ["EPCAM", "KRT8", "KRT14"],
    "T Cells": ["CD3D", "CD3E", "IL7R"],
    "Myeloid Cells": ["LST1", "TYROBP", "C1QC"],
    "B Cells": ["MS4A1", "CD79A", "CD74"],
    "Plasma Cells": ["MZB1", "JCHAIN", "TNFRSF17"],
    "Endothelial Cells": ["PECAM1", "VWF", "EMCN"],
    "Stromal Cells": ["COL1A1", "COL1A2", "DCN"],
    "Perivascular Cells": ["RGS5", "PDGFRB", "MCAM"],
}


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")


def main() -> None:
    started = time.time()
    OUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
    ROUND2_DEG_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    mapping = pd.read_csv(MAPPING, dtype={"raw_cluster": str})
    if mapping["raw_cluster"].duplicated().any() or len(mapping) != 15:
        raise ValueError("Broad annotation mapping must contain 15 unique raw clusters.")
    if set(mapping["leiden_coarse"]) != set(CATEGORY_ORDER):
        raise ValueError("Broad annotation mapping labels do not match the declared order.")
    if not mapping["review_depth"].eq(200).all():
        raise ValueError("All raw clusters must use the same 200-gene review depth.")

    round1_files = sorted(ROUND1_DEG_DIR.glob("*_degs_*.csv"))
    if len(round1_files) != 15:
        raise ValueError("Round-1 full-length DEG evidence is incomplete.")
    for path in round1_files:
        if sum(1 for _ in path.open(encoding="utf-8")) - 1 != 27716:
            raise ValueError(f"Round-1 DEG table is not full length: {path}")

    adata = sc.read_h5ad(INPUT_H5AD)
    if adata.raw is None or adata.raw.n_vars != 27716:
        raise ValueError("Broad annotation requires preserved raw normalized/log expression.")
    if RAW_KEY not in adata.obs or adata.obs[RAW_KEY].isna().any():
        raise ValueError(f"Missing selected raw clustering column: {RAW_KEY}")
    raw_groups = set(adata.obs[RAW_KEY].astype(str).unique())
    if raw_groups != set(mapping["raw_cluster"]):
        raise ValueError("Broad mapping does not cover the selected raw clusters exactly.")

    label_map = mapping.set_index("raw_cluster")["leiden_coarse"].to_dict()
    coarse = adata.obs[RAW_KEY].astype(str).map(label_map)
    if coarse.isna().any():
        raise ValueError("Broad mapping produced missing labels.")
    adata.obs[COARSE_KEY] = pd.Categorical(
        coarse, categories=CATEGORY_ORDER, ordered=True
    )
    adata.obs["cell_type"] = pd.Categorical(
        coarse, categories=CATEGORY_ORDER, ordered=True
    )
    if not adata.obs[COARSE_KEY].astype(str).equals(adata.obs["cell_type"].astype(str)):
        raise ValueError("Initial cell_type is not identical to leiden_coarse.")

    print(
        f"[Broad annotation] labels applied: cells={adata.n_obs}, "
        f"raw_clusters={len(raw_groups)}, broad_labels={adata.obs[COARSE_KEY].nunique()}",
        flush=True,
    )
    sc.tl.rank_genes_groups(
        adata,
        groupby=COARSE_KEY,
        method=METHOD,
        use_raw=True,
        n_genes=adata.raw.n_vars,
        key_added="rank_genes_groups_leiden_coarse",
    )
    result = adata.uns["rank_genes_groups_leiden_coarse"]
    audit_rows: list[dict[str, object]] = []
    for group in CATEGORY_ORDER:
        frame = pd.DataFrame(
            {
                "gene": np.asarray(result["names"][group]).astype(str),
                "score": np.asarray(result["scores"][group], dtype=float),
                "logfoldchanges": np.asarray(
                    result["logfoldchanges"][group], dtype=float
                ),
                "pvals": np.asarray(result["pvals"][group], dtype=float),
                "pvals_adj": np.asarray(result["pvals_adj"][group], dtype=float),
            }
        )
        if len(frame) != adata.raw.n_vars or frame["gene"].duplicated().any():
            raise ValueError(f"Round-2 DEG table is not full length and unique: {group}")
        out = (
            ROUND2_DEG_DIR
            / f"{safe_name(group)}_degs_leiden_coarse_pcs20_nn30_res0p8.csv"
        )
        frame.to_csv(out, index=False)
        audit_rows.append(
            {
                "leiden_coarse": group,
                "n_cells": int((adata.obs[COARSE_KEY].astype(str) == group).sum()),
                "n_deg_rows": len(frame),
                "full_length": True,
                "use_raw": True,
                "method": METHOD,
                "deg_csv": str(out),
            }
        )
        print(f"[Broad DEG round 2] wrote group={group} rows={len(frame)}", flush=True)
    pd.DataFrame(audit_rows).to_csv(
        TABLE_DIR / "round2_leiden_coarse_deg_audit.csv", index=False
    )

    marker_genes = [gene for genes in MARKERS.values() for gene in genes]
    missing_markers = sorted(set(marker_genes) - set(adata.raw.var_names.astype(str)))
    if missing_markers:
        raise ValueError(f"Broad dotplot markers are absent from adata.raw: {missing_markers}")
    pd.DataFrame(
        [
            {"leiden_coarse": group, "marker_order": order, "gene": gene}
            for group, genes in MARKERS.items()
            for order, gene in enumerate(genes, start=1)
        ]
    ).to_csv(TABLE_DIR / "broad_marker_dotplot_genes.csv", index=False)

    counts = adata.obs[COARSE_KEY].value_counts(sort=False).reindex(CATEGORY_ORDER)
    composition = pd.DataFrame(
        {
            "leiden_coarse": CATEGORY_ORDER,
            "n_cells": counts.to_numpy(dtype=int),
            "fraction": (counts / counts.sum()).to_numpy(dtype=float),
        }
    )
    composition.to_csv(TABLE_DIR / "broad_cell_type_composition.csv", index=False)
    sample_counts = pd.crosstab(adata.obs["sample"], adata.obs[COARSE_KEY]).reindex(
        columns=CATEGORY_ORDER, fill_value=0
    )
    sample_fractions = sample_counts.div(sample_counts.sum(axis=1), axis=0)
    sample_counts.to_csv(TABLE_DIR / "broad_cell_type_by_sample_counts.csv")
    sample_fractions.to_csv(TABLE_DIR / "broad_cell_type_by_sample_fractions.csv")
    # A stacked composition plot has no Scanpy plotting equivalent; use the
    # official pandas/Matplotlib interface and record that special case here.
    ax = sample_fractions.plot(kind="bar", stacked=True, figsize=(12, 4), width=0.85)
    ax.set_xlabel("sample")
    ax.set_ylabel("cell fraction")
    ax.legend(title=COARSE_KEY, bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=7)
    ax.figure.tight_layout()
    ax.figure.savefig(FIGURE_DIR / "broad_cell_type_composition_by_sample.pdf")
    plt.close(ax.figure)

    sc.settings.figdir = FIGURE_DIR
    sc.settings.file_format_figs = "pdf"
    sc.set_figure_params(figsize=(3, 3), dpi=150, fontsize=8)
    sc.pl.umap(
        adata,
        color=RAW_KEY,
        show=False,
        save="_raw_leiden_res0p8.pdf",
    )
    sc.pl.umap(
        adata,
        color=COARSE_KEY,
        show=False,
        save="_leiden_coarse.pdf",
    )
    sc.pl.dotplot(
        adata,
        var_names=MARKERS,
        groupby=COARSE_KEY,
        standard_scale="var",
        use_raw=True,
        show=False,
        save="_all_leiden_coarse.pdf",
    )

    adata.uns["broad_annotation_parameters"] = {
        "source_h5ad": str(INPUT_H5AD),
        "raw_cluster_column": RAW_KEY,
        "mapping_csv": str(MAPPING),
        "review_depth": 200,
        "review_depth_reason": (
            "Cluster 5 was cell-cycle dominated in the first 50/100 rows; "
            "EPCAM appeared at rank 131, so all clusters were reviewed to 200."
        ),
        "round1_deg_dir": str(ROUND1_DEG_DIR),
        "round2_deg_dir": str(ROUND2_DEG_DIR),
        "deg_method": METHOD,
        "deg_use_raw": True,
        "category_order": CATEGORY_ORDER,
        "annotation_revision": (
            "User-requested broad-level merge of raw cluster 11 "
            "Myoepithelial Cells into Epithelial Cells; raw cluster labels retained."
        ),
        "marker_dotplot_use_raw": True,
        "marker_dotplot_standard_scale": "var",
        "seed": SEED,
    }

    versions = {
        "python": sys.version.split()[0],
        "scanpy": package_version("scanpy"),
        "anndata": package_version("anndata"),
        "matplotlib": package_version("matplotlib"),
        "code": str(CODE_PATH),
        "seed": str(SEED),
    }
    (TABLE_DIR / "round2_package_versions.txt").write_text(
        "\n".join(f"{key}={value}" for key, value in versions.items()) + "\n",
        encoding="utf-8",
    )
    (TABLE_DIR / "round2_readme.txt").write_text(
        f"""BRCA broad annotation and round-2 DEGs

Input: {INPUT_H5AD}
Output: {OUT_H5AD}
Raw cluster column retained: {RAW_KEY}
Broad annotation columns: {COARSE_KEY}, cell_type
Mapping evidence: {MAPPING}
Round-1 full DEGs: {ROUND1_DEG_DIR}
Round-2 full DEGs: {ROUND2_DEG_DIR}

All clusters were reviewed to a uniform depth of 200 rows from the saved
full-length round-1 DEG tables. The expansion from 50/100 was required because
cluster 5 was cycling-dominant and EPCAM appeared at rank 131. Round-2 DEGs use
adata.raw normalized/log expression with method=t-test and use_raw=True. The
default marker dotplot uses exactly three ordered canonical genes per observed
broad category, use_raw=True, and standard_scale=var.
Raw cluster 11 retains its technical cluster identity but is merged into
Epithelial Cells at the broad-label level by user request.
""",
        encoding="utf-8",
    )

    print(f"[Broad annotation] writing {OUT_H5AD}", flush=True)
    adata.write_h5ad(OUT_H5AD, compression="gzip")
    figures = sorted(path.name for path in FIGURE_DIR.glob("*.pdf"))
    report = {
        "n_cells": int(adata.n_obs),
        "n_raw_clusters": int(adata.obs[RAW_KEY].nunique()),
        "n_broad_labels": int(adata.obs[COARSE_KEY].nunique()),
        "broad_label_counts": {
            key: int(value) for key, value in counts.to_dict().items()
        },
        "n_round1_full_deg_csvs": len(round1_files),
        "n_round2_full_deg_csvs": len(audit_rows),
        "review_depth": 200,
        "figures": figures,
        "output_h5ad": str(OUT_H5AD),
        "elapsed_seconds": time.time() - started,
    }
    (TABLE_DIR / "round2_completion.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
