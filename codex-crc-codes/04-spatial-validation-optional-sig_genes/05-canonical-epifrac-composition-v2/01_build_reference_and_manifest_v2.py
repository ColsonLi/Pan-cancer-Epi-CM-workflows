#!/usr/bin/env python3
"""Build DEG-derived subtype pseudobulk and the 37-region spatial manifest."""

from __future__ import annotations

import hashlib
import json
import platform
import random
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy
from scipy import sparse

SEED = 42
random.seed(SEED)
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/crc_codex")
WF = ROOT / "epi-cm-core-workflow"
TASK = "05-canonical-epifrac-composition-v2"
BLOCK = "04-spatial-validation-optional-sig_genes"
CODE_DIR = WF / f"codes/{BLOCK}/{TASK}"
TABLE_DIR = WF / f"tables/{BLOCK}/{TASK}/01-reference-and-manifest"
H5AD_DIR = WF / f"h5ad/{BLOCK}/{TASK}/01-reference-and-manifest"
UNPROCESSED_SPATIAL = WF / "h5ad/04-spatial-validation-optional/01-input-audit/adata_xenium_unprocessed.h5ad"
FINAL_SC = WF / "h5ad/02-cell_subtype_integration_clustering/06-project-subtypes-to-full-adata/adata_anno_cellsubtype.h5ad"
DEG_ROOT = WF / "tables/02-cell_subtype_integration_clustering/05-subtype-deg-annotation"
CM_LOADING = WF / "tables/03-epi-cm-discovery/01-cm-lineage-analysis/02_balanced_joint_nmf/loading_df_cell_subtype_by_CM_fraction.csv"
REFERENCE_OUT = H5AD_DIR / "adata_cell_subtype_pseudobulk_deg_panel.h5ad"


def dec(values) -> list[str]:
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values]


def categorical(group: h5py.Group, key: str) -> tuple[list[str], np.ndarray]:
    node = group[key]
    if isinstance(node, h5py.Group) and "categories" in node:
        return dec(node["categories"][()]), node["codes"][()].astype(np.int32)
    values = np.asarray(dec(node[()]), dtype=object)
    cats, codes = np.unique(values, return_inverse=True)
    return cats.tolist(), codes.astype(np.int32)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def build_sample_manifest() -> pd.DataFrame:
    with h5py.File(UNPROCESSED_SPATIAL, "r") as f:
        obs = f["obs"]
        sample_cats, sample_codes = categorical(obs, "name")
        tissue_cats, tissue_codes = categorical(obs, "tissue_region")
        patient_cats, patient_codes = categorical(obs, "patient_id")
    rows = []
    for sample_i, sample in enumerate(sample_cats):
        mask = sample_codes == sample_i
        tissues = np.unique(tissue_codes[mask])
        patients = np.unique(patient_codes[mask])
        if len(tissues) != 1 or len(patients) != 1:
            raise ValueError(f"{sample}: non-unique tissue/patient mapping")
        tissue = tissue_cats[int(tissues[0])]
        patient = patient_cats[int(patients[0])]
        status = "tumor" if tissue in {"core", "margin"} else "normal"
        rows.append(
            {
                "sample": sample,
                "patient_id": patient,
                "tissue_region": tissue,
                "status": status,
                "n_spatial_cells": int(mask.sum()),
                "include_all_samples": True,
                "include_tumor_only": status == "tumor",
                "exclusion_reason_all_samples": "",
                "exclusion_reason_tumor_only": "" if status == "tumor" else "status=normal",
                "status_source": "crca_xenium.h5ad:obs[tissue_region]; core/margin=tumor, normal=normal",
            }
        )
    result = pd.DataFrame(rows).sort_values("sample", kind="stable").reset_index(drop=True)
    if len(result) != 37:
        raise ValueError(f"Expected 37 spatial samples, found {len(result)}")
    return result


def build_deg_marker_table(panel: set[str], reference_genes: set[str]) -> pd.DataFrame:
    rows = []
    paths = sorted(DEG_ROOT.glob("*/degs_cell_subtype*/*.csv"))
    subtype_to_paths: dict[str, list[Path]] = {}
    for path in paths:
        subtype = path.stem.split("_degs_", 1)[0]
        subtype_to_paths.setdefault(subtype, []).append(path)
    duplicated_sources = {key: value for key, value in subtype_to_paths.items() if len(value) != 1}
    if len(paths) != 87 or len(subtype_to_paths) != 87 or duplicated_sources:
        raise ValueError(
            f"Require exactly one post-annotation DEG CSV for each of 87 final subtypes; "
            f"files={len(paths)}, subtypes={len(subtype_to_paths)}, duplicates={duplicated_sources}"
        )
    for path in paths:
        subtype = path.stem.split("_degs_", 1)[0]
        deg = pd.read_csv(path)
        required = {"gene", "score", "logfoldchanges", "pvals_adj"}
        missing = sorted(required - set(deg.columns))
        if missing:
            raise KeyError(f"{path}: missing required post-annotation DEG columns {missing}")
        original_score = deg["score"].copy()
        original_lfc = deg["logfoldchanges"].copy()
        original_padj = deg["pvals_adj"].copy()
        deg["logfoldchanges"] = pd.to_numeric(deg["logfoldchanges"], errors="coerce")
        deg["pvals_adj"] = pd.to_numeric(deg["pvals_adj"], errors="coerce")
        deg["score"] = pd.to_numeric(deg["score"], errors="coerce")
        if (
            deg.loc[original_score.notna(), "score"].isna().any()
            or deg.loc[original_lfc.notna(), "logfoldchanges"].isna().any()
            or deg.loc[original_padj.notna(), "pvals_adj"].isna().any()
        ):
            raise ValueError(f"{path}: score/logfoldchanges/pvals_adj contains non-numeric values")
        deg["gene"] = deg["gene"].astype(str)

        # FIXED ORDER: complete single-cell marker selection before any Xenium
        # or reference-gene intersection is applied.
        deg = deg[(deg["pvals_adj"] < 0.05) & (deg["logfoldchanges"] > 0)].copy()
        deg = deg.sort_values("score", ascending=False, kind="stable")
        deg = deg.drop_duplicates("gene", keep="first").head(100)
        if deg.empty:
            raise ValueError(f"{path}: no positive unique DEG remains for {subtype}")
        deg["deg_rank"] = np.arange(1, len(deg) + 1)
        for row in deg.itertuples(index=False):
            rows.append(
                {
                    "cell_subtype": subtype,
                    "gene": str(row.gene),
                    "deg_rank": int(row.deg_rank),
                    "score": float(row.score),
                    "logfoldchanges": float(row.logfoldchanges),
                    "pvals_adj": float(row.pvals_adj),
                    "selection_stage": "post_annotation_padj_lt_0p05_positive_score_desc_deduplicated_top100_before_xenium_intersection",
                    "source_deg_csv": str(path),
                }
            )
    selected = pd.DataFrame(rows)
    if selected["cell_subtype"].nunique() != 87:
        raise ValueError("Positive DEG selection is missing one or more final subtypes")

    # FIXED INTERSECTION POSITION: only after every subtype has completed its
    # independent top-100 marker selection may genes be intersected with the
    # single-cell raw feature set and Xenium panel.
    selected["in_spatial_panel"] = selected["gene"].isin(panel)
    selected["in_reference_raw"] = selected["gene"].isin(reference_genes)
    selected["used_for_tangram"] = selected["in_spatial_panel"] & selected["in_reference_raw"]
    selected.to_csv(TABLE_DIR / "marker_genes_used.csv", index=False)
    return selected.loc[selected["used_for_tangram"]].copy()


def build_pseudobulk(marker_genes: list[str]) -> ad.AnnData:
    with h5py.File(FINAL_SC, "r") as f:
        subtype_categories, subtype_codes = categorical(f["obs"], "cell_subtype")
        raw_var_names = dec(f["raw"]["var"]["_index"][()])
        gene_to_i = {gene: i for i, gene in enumerate(raw_var_names)}
        selected_i = np.asarray([gene_to_i[g] for g in marker_genes], dtype=np.int64)
        x = f["raw"]["X"]
        if not isinstance(x, h5py.Group) or "indptr" not in x:
            raise TypeError("Expected sparse CSR adata.raw.X")
        data = x["data"][()]
        indices = x["indices"][()]
        indptr = x["indptr"][()]
        raw_x = sparse.csr_matrix(
            (data, indices, indptr), shape=(len(subtype_codes), len(raw_var_names))
        )

    selected = raw_x[:, selected_i].tocsr()
    valid = subtype_codes >= 0
    rows = np.arange(len(subtype_codes), dtype=np.int64)[valid]
    groups = subtype_codes[valid]
    group_indicator = sparse.csr_matrix(
        (np.ones(len(rows), dtype=np.float32), (groups, rows)),
        shape=(len(subtype_categories), len(subtype_codes)),
    )
    sums = (group_indicator @ selected).toarray().astype(np.float32)
    counts = np.bincount(groups, minlength=len(subtype_categories)).astype(np.float32)
    means = sums / counts[:, None]
    obs = pd.DataFrame(index=pd.Index(subtype_categories, name="cell_subtype"))
    obs["cell_subtype"] = obs.index.astype(str)
    obs["cell_class"] = np.where(obs.index.str.startswith("Epi_"), "epithelial", "non_epithelial")
    obs["n_reference_cells"] = counts.astype(np.int64)
    var = pd.DataFrame(index=pd.Index(marker_genes, name="gene"))
    result = ad.AnnData(X=means, obs=obs, var=var)
    result.uns["expression_source"] = "adata_anno_cellsubtype.raw normalized/log expression"
    result.uns["marker_source"] = "all significant positive cell_subtype DEG genes intersected with Xenium panel"
    result.uns["random_seed"] = SEED
    return result


def main() -> None:
    started = time.time()
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    H5AD_DIR.mkdir(parents=True, exist_ok=True)
    manifest = build_sample_manifest()
    manifest.to_csv(TABLE_DIR / "spatial_sample_scope_manifest.csv", index=False)

    with h5py.File(UNPROCESSED_SPATIAL, "r") as f:
        panel = set(dec(f["var"]["_index"][()]))
    with h5py.File(FINAL_SC, "r") as f:
        reference_genes = set(dec(f["raw"]["var"]["_index"][()]))
    if len(panel) != 380:
        raise ValueError(f"CRC Xenium panel contract failed: expected 380 genes, found {len(panel)}")
    if len(panel & reference_genes) != 356:
        raise ValueError(
            f"CRC Xenium/single-cell raw overlap contract failed: expected 356 genes, found {len(panel & reference_genes)}"
        )

    markers = build_deg_marker_table(panel, reference_genes)
    markers.to_csv(TABLE_DIR / "cell_subtype_deg_markers_in_spatial_panel.csv", index=False)
    selected_markers = pd.read_csv(TABLE_DIR / "marker_genes_used.csv")
    marker_union = (
        selected_markers.groupby("gene", as_index=False)
        .agg(
            n_subtypes_top100=("cell_subtype", "nunique"),
            best_score=("score", "max"),
            in_spatial_panel=("in_spatial_panel", "max"),
            in_reference_raw=("in_reference_raw", "max"),
            used_for_tangram=("used_for_tangram", "max"),
        )
        .sort_values(["used_for_tangram", "best_score", "gene"], ascending=[False, False, True], kind="stable")
    )
    marker_union.to_csv(TABLE_DIR / "marker_gene_union_intersection.csv", index=False)
    marker_genes = sorted(markers["gene"].unique().tolist())
    if len(marker_genes) != 247:
        raise ValueError(
            f"Post-top100 Xenium/reference three-way intersection contract failed: expected 247 genes, found {len(marker_genes)}"
        )
    pd.DataFrame({"gene": marker_genes}).to_csv(TABLE_DIR / "tangram_common_deg_genes.csv", index=False)
    pd.DataFrame(
        [
            {"metric": "subtype_top100_marker_rows", "value": len(selected_markers)},
            {"metric": "subtype_top100_unique_gene_union", "value": selected_markers["gene"].nunique()},
            {"metric": "xenium_panel_genes", "value": len(panel)},
            {"metric": "xenium_panel_genes_in_single_cell_raw", "value": len(panel & reference_genes)},
            {"metric": "top100_union_genes_in_xenium_panel", "value": len(set(selected_markers["gene"]) & panel)},
            {"metric": "final_tangram_three_way_intersection", "value": len(marker_genes)},
            {"metric": "top100_union_genes_not_in_xenium_panel", "value": len(set(selected_markers["gene"]) - panel)},
            {"metric": "panel_and_raw_genes_not_selected_top100", "value": len((panel & reference_genes) - set(selected_markers["gene"]))},
            {"metric": "panel_genes_absent_from_single_cell_raw", "value": len(panel - reference_genes)},
        ]
    ).to_csv(TABLE_DIR / "marker_gene_intersection_summary.csv", index=False)
    reference = build_pseudobulk(marker_genes)
    reference.obs[["n_reference_cells"]].rename_axis("subtype").assign(
        n_cells_total=lambda x: x["n_reference_cells"],
        n_cells_used=lambda x: x["n_reference_cells"],
        all_cells_used=True,
    )[["n_cells_total", "n_cells_used", "all_cells_used"]].to_csv(TABLE_DIR / "reference_cells_used.csv")

    loadings = pd.read_csv(CM_LOADING).set_index("cell_subtype")
    missing_nodes = sorted(set(loadings.index) - set(reference.obs_names))
    epi_names = sorted(reference.obs_names[reference.obs["cell_class"].eq("epithelial")])
    if missing_nodes or len(epi_names) != 11 or loadings.shape != (76, 15):
        raise ValueError(
            f"Reference contract failed: missing_nodes={missing_nodes}, epi={len(epi_names)}, H={loadings.shape}"
        )
    loadings.to_csv(TABLE_DIR / "cm_loading_fraction_non_epi_subtype_by_CM.csv")
    pd.DataFrame({"cell_subtype": epi_names}).to_csv(TABLE_DIR / "epithelial_subtypes.csv", index=False)
    reference.write_h5ad(REFERENCE_OUT, compression="lzf")

    audit = {
        "completed_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "unprocessed_spatial_h5ad": str(UNPROCESSED_SPATIAL),
        "unprocessed_spatial_sha256": sha256(UNPROCESSED_SPATIAL),
        "unprocessed_spatial_was_modified": False,
        "single_cell_reference": str(FINAL_SC),
        "single_cell_expression_source": "raw",
        "n_spatial_samples": len(manifest),
        "n_spatial_patients": int(manifest["patient_id"].nunique()),
        "n_reference_subtypes": int(reference.n_obs),
        "n_epithelial_subtypes": len(epi_names),
        "n_non_epithelial_cm_nodes": int(loadings.shape[0]),
        "n_cm": int(loadings.shape[1]),
        "n_xenium_panel_genes": len(panel),
        "n_xenium_panel_genes_in_single_cell_raw": len(panel & reference_genes),
        "n_top100_marker_rows": len(selected_markers),
        "n_top100_unique_gene_union": int(selected_markers["gene"].nunique()),
        "n_tangram_deg_panel_genes": int(reference.n_vars),
        "elapsed_seconds": round(time.time() - started, 3),
        "seed": SEED,
    }
    (TABLE_DIR / "reference_and_manifest_summary.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    (TABLE_DIR / "readme.txt").write_text(
        "Block 04 reference and manifest.\n"
        "The unprocessed Xenium audit H5AD is byte-identical to the source and was not normalized.\n"
        "Pseudobulk rows are all 87 final cell_subtypes; values are means from single-cell adata.raw.\n"
        "Markers are saved post-annotation subtype DEGs filtered to pvals_adj < 0.05 and positive logFC, stable score-descending, deduplicated, top 100 per subtype.\n"
        "Tangram genes are the three-way intersection of that marker union, single-cell adata.raw, and the 380-gene Xenium panel.\n"
        "CM activity downstream is the Tangram-projected non-epithelial subtype abundance multiplied by canonical H loading fractions.\n",
        encoding="utf-8",
    )
    (TABLE_DIR / "package_versions.txt").write_text(
        "\n".join(
            [
                f"python={platform.python_version()}",
                f"anndata={ad.__version__}",
                f"numpy={np.__version__}",
                f"pandas={pd.__version__}",
                f"scipy={scipy.__version__}",
                f"h5py={h5py.__version__}",
                f"code_file={Path(__file__)}",
                f"seed={SEED}",
            ]
        ) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
