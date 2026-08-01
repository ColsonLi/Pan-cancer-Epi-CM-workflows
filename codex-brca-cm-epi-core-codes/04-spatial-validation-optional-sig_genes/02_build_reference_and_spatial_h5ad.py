#!/usr/bin/env python3
"""Build a significant-positive-DEG reference and six clean spatial h5ads."""

from __future__ import annotations

import importlib.metadata
import importlib.util
import json
import random
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

SEED = 42
random.seed(SEED)
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/brca_cm-epi_codex")
WORKFLOW = ROOT / "epi-cm-core-workflow"
DATA = Path("/mnt/disk18t/lr_xcy/riku/external_atlas_for_cm-epi_validation/Breast_Wu2021_Zenodo4739739/spatial")
TABLE_DIR = WORKFLOW / "tables/04-spatial-validation-optional-sig_genes/01-input-audit-and-reference-v4"
H5AD_DIR = WORKFLOW / "h5ad/04-spatial-validation-optional-sig_genes/01-input-audit-and-reference"
H5AD_DIR.mkdir(parents=True, exist_ok=False)

DEG_ROOT = WORKFLOW / "tables/02-cell_subtype_integration_clustering/05-subtype-deg-annotation"
LINEAGE_H5AD_ROOT = WORKFLOW / "h5ad/02-cell_subtype_integration_clustering/05-subtype-deg-annotation"
PROJECTED_ATLAS = WORKFLOW / "h5ad/02-cell_subtype_integration_clustering/06-project-subtypes-to-full-adata/adata_anno_cellsubtype.h5ad"
SAMPLES = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]
TOP_MARKERS = 100

HELPER_PATH = Path(__file__).with_name("01_audit_spatial_inputs_v3.py")
SPEC = importlib.util.spec_from_file_location("spatial_audit_helpers", HELPER_PATH)
if SPEC is None or SPEC.loader is None:
    raise ImportError(HELPER_PATH)
HELPER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(HELPER)


def pkg(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def discover_deg_tables() -> dict[str, Path]:
    paths = sorted(DEG_ROOT.glob("*/degs_cell_subtype_*/*_degs_cell_subtype_*.csv"))
    result: dict[str, Path] = {}
    for path in paths:
        subtype = path.name.split("_degs_cell_subtype_", 1)[0]
        if subtype in result:
            raise ValueError(f"Duplicate post-annotation DEG table for {subtype}")
        result[subtype] = path
    if len(result) != 68:
        raise ValueError(f"Expected 68 post-annotation subtype DEG tables, found {len(result)}")
    return result


def select_markers(deg_paths: dict[str, Path]) -> tuple[pd.DataFrame, list[str], pd.DataFrame]:
    rows: list[dict[str, object]] = []
    selection_rows: list[dict[str, object]] = []
    union: list[str] = []
    seen_union: set[str] = set()
    for subtype, path in sorted(deg_paths.items()):
        frame = pd.read_csv(path)
        required = {"gene", "pvals_adj", "score", "logfoldchanges"}
        if not required.issubset(frame.columns):
            raise ValueError(f"{path}: missing {sorted(required.difference(frame.columns))}")
        for column in ["pvals_adj", "score", "logfoldchanges"]:
            frame[column] = pd.to_numeric(frame[column], errors="raise")
            if not np.isfinite(frame[column].to_numpy(dtype=float)).all():
                raise ValueError(f"{path}: non-finite {column} values")
        frame["source_order"] = np.arange(len(frame))
        # FIXED: canonical significant-positive DEG contract from the current skill.
        significant = frame.loc[frame["pvals_adj"].lt(0.05) & frame["logfoldchanges"].gt(0)].copy()
        significant = significant.sort_values(["score", "source_order"], ascending=[False, True], kind="stable")
        significant_unique = significant.drop_duplicates("gene", keep="first")
        selected = significant_unique.head(TOP_MARKERS)
        if selected.empty:
            raise ValueError(f"No significant positive DEG marker for {subtype}: {path}")
        selection_rows.append({
            "subtype": subtype,
            "source_deg_csv": str(path.resolve()),
            "n_source_rows": len(frame),
            "n_significant_positive_rows": len(significant),
            "n_significant_positive_unique_genes": len(significant_unique),
            "n_selected": len(selected),
            "fewer_than_100_selected": len(selected) < TOP_MARKERS,
            "pvals_adj_threshold": 0.05,
            "logfoldchanges_rule": "> 0",
            "marker_order": "score descending; stable source order for exact ties",
        })
        for rank, record in enumerate(selected.itertuples(index=False), start=1):
            gene = str(record.gene)
            rows.append({
                "subtype": subtype,
                "gene": gene,
                "marker_rank": rank,
                "pvals_adj": float(record.pvals_adj),
                "score": float(record.score),
                "logfoldchanges": float(record.logfoldchanges),
                "source_deg_csv": str(path.resolve()),
            })
            if gene not in seen_union:
                seen_union.add(gene)
                union.append(gene)
    return pd.DataFrame(rows), union, pd.DataFrame(selection_rows)


def build_reference(marker_table: pd.DataFrame, marker_union: list[str], deg_subtypes: set[str]) -> tuple[ad.AnnData, pd.DataFrame]:
    lineage_paths = sorted(LINEAGE_H5AD_ROOT.glob("*/*.h5ad"))
    if len(lineage_paths) != 8:
        raise ValueError(f"Expected 8 final lineage h5ads, found {len(lineage_paths)}")
    means: dict[str, np.ndarray] = {}
    cell_rows: list[dict[str, object]] = []
    raw_gene_reference: list[str] | None = None
    for path in lineage_paths:
        a = ad.read_h5ad(path, backed="r")
        if a.raw is None:
            raise ValueError(f"Missing adata.raw: {path}")
        if "cell_subtype" not in a.obs.columns:
            raise ValueError(f"Missing cell_subtype: {path}")
        if raw_gene_reference is None:
            raw_gene_reference = a.raw.var_names.astype(str).tolist()
        elif raw_gene_reference != a.raw.var_names.astype(str).tolist():
            raise ValueError(f"adata.raw gene order differs: {path}")
        missing = [gene for gene in marker_union if gene not in a.raw.var_names]
        if missing:
            raise ValueError(f"{path}: marker genes absent from raw: {missing[:20]}")
        gene_index = a.raw.var_names.get_indexer(marker_union)
        expression = a.raw.X[:, gene_index]
        subtype_series = a.obs["cell_subtype"].astype(str)
        for subtype in sorted(subtype_series.unique()):
            mask = subtype_series.eq(subtype).to_numpy()
            n_cells = int(mask.sum())
            if n_cells == 0:
                raise ValueError(f"Zero-cell subtype: {subtype}")
            means[subtype] = np.asarray(expression[mask, :].mean(axis=0)).ravel().astype(np.float32)
            cell_rows.append({
                "subtype": subtype,
                "lineage_h5ad": str(path.resolve()),
                "n_cells_total": n_cells,
                "n_cells_used": n_cells,
                "all_cells_used": True,
                "expression_source": "adata.raw normalized/log expression",
            })
        a.file.close()
    if set(means) != deg_subtypes:
        raise ValueError(f"DEG/reference subtype mismatch: DEG-only={sorted(deg_subtypes-set(means))}, reference-only={sorted(set(means)-deg_subtypes)}")
    order = sorted(means)
    obs = pd.DataFrame({"cell_subtype": order, "n_cells": [next(row["n_cells_total"] for row in cell_rows if row["subtype"] == subtype) for subtype in order]}, index=pd.Index(order, name="pseudocell"))
    obs["cell_subtype"] = pd.Categorical(obs["cell_subtype"], categories=order)
    var = pd.DataFrame(index=pd.Index(marker_union, name="gene"))
    var["n_subtypes_marker"] = marker_table.groupby("gene")["subtype"].nunique().reindex(marker_union).fillna(0).astype(int).to_numpy()
    asc = ad.AnnData(X=np.vstack([means[subtype] for subtype in order]), obs=obs, var=var)
    asc.uns["reference_method"] = "subtype mean from every annotated cell using adata.raw normalized/log expression"
    asc.uns["marker_method"] = "pvals_adj < 0.05 and logfoldchanges > 0; top 100 unique genes per subtype by score descending"
    asc.uns["source_projected_atlas"] = str(PROJECTED_ATLAS.resolve())
    asc.uns["seed"] = SEED
    return asc, pd.DataFrame(cell_rows).sort_values("subtype")


def build_spatial_h5ad(sample: str) -> tuple[ad.AnnData, dict[str, object]]:
    fdir = DATA / "filtered_count_matrices" / f"{sample}_filtered_count_matrix"
    meta_path = DATA / "metadata" / f"{sample}_metadata.csv"
    spatial_dir = DATA / "spatial" / f"{sample}_spatial"
    barcodes = HELPER.lines(fdir / "barcodes.tsv.gz")
    genes = HELPER.unique_matrix_symbols(HELPER.symbols(fdir / "features.tsv.gz"))
    counts = HELPER.matrix(fdir / "matrix.mtx.gz").T.tocsr().astype(np.float32)
    if counts.shape != (len(barcodes), len(genes)):
        raise ValueError(f"{sample}: matrix/feature/barcode shape mismatch")
    meta = pd.read_csv(meta_path, index_col=0)
    meta.index = meta.index.astype(str)
    if barcodes != meta.index.tolist():
        raise ValueError(f"{sample}: metadata order differs from filtered barcode order")
    pos = pd.read_csv(spatial_dir / "tissue_positions_list.csv", header=None, names=["barcode", "source_in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]).set_index("barcode")
    pos.index = pos.index.astype(str)
    coordinates = pos.loc[barcodes]
    obs = meta.copy()
    obs["barcode"] = barcodes
    obs["sample"] = sample
    for column in coordinates.columns:
        obs[column] = coordinates[column].to_numpy()
    obs["spatial_x"] = obs["pxl_col"].astype(float)
    obs["spatial_y"] = obs["pxl_row"].astype(float)
    obs.index = pd.Index([f"{barcode}_{sample}" for barcode in barcodes], name="spot_id")
    var = pd.DataFrame(index=pd.Index(genes, name="gene"))
    var["source_feature_name"] = genes
    a = ad.AnnData(X=counts, obs=obs, var=var)
    a.layers["counts"] = a.X.copy()
    a.obsm["spatial"] = a.obs[["spatial_x", "spatial_y"]].to_numpy(dtype=float)
    a.uns["sample"] = sample
    a.uns["expression_state"] = "raw integer-like counts from curated filtered count matrix"
    a.uns["source_count_matrix"] = str((fdir / "matrix.mtx.gz").resolve())
    a.uns["source_metadata"] = str(meta_path.resolve())
    a.uns["source_positions"] = str((spatial_dir / "tissue_positions_list.csv").resolve())
    if not a.obs_names.is_unique or not a.var_names.is_unique:
        raise ValueError(f"{sample}: non-unique h5ad names")
    total_counts = np.asarray(a.X.sum(axis=1)).ravel()
    if not np.allclose(total_counts, pd.to_numeric(meta["nCount_RNA"]).to_numpy(), rtol=0, atol=0):
        raise ValueError(f"{sample}: spatial h5ad count sums differ from metadata")
    return a, {
        "sample": sample,
        "n_spots": a.n_obs,
        "n_genes": a.n_vars,
        "matrix_nnz": int(a.X.nnz),
        "n_source_in_tissue_0_retained": int(a.obs["source_in_tissue"].astype(int).eq(0).sum()),
        "all_coordinates_finite": bool(np.isfinite(a.obsm["spatial"]).all()),
        "count_sum_matches_metadata": True,
    }


def main() -> None:
    audit = json.loads((TABLE_DIR / "input_audit_completion.json").read_text())
    if audit.get("status") != "completed":
        raise RuntimeError("Final v4 spatial input audit is not complete")
    deg_paths = discover_deg_tables()
    marker_table, marker_union, selection_summary = select_markers(deg_paths)
    reference, cells_used = build_reference(marker_table, marker_union, set(deg_paths))
    reference_path = H5AD_DIR / "subtype_pseudobulk_reference_top100_significant_positive_deg.h5ad"
    reference.write_h5ad(reference_path, compression="gzip")
    marker_table.to_csv(TABLE_DIR / "marker_genes_used.csv", index=False)
    selection_summary.to_csv(TABLE_DIR / "marker_selection_summary.csv", index=False)
    cells_used.to_csv(TABLE_DIR / "reference_cells_used.csv", index=False)
    marker_union_table = pd.DataFrame({"gene": marker_union})
    marker_union_table["n_subtypes_marker"] = marker_table.groupby("gene")["subtype"].nunique().reindex(marker_union).to_numpy()
    marker_union_table.to_csv(TABLE_DIR / "marker_gene_union.csv", index=False)

    spatial_rows = []
    intersection_rows = []
    reference_genes = set(reference.var_names.astype(str))
    for sample in SAMPLES:
        spatial, row = build_spatial_h5ad(sample)
        output = H5AD_DIR / f"{sample}_spatial_raw_counts.h5ad"
        spatial.write_h5ad(output, compression="gzip")
        row["output_h5ad"] = str(output.resolve())
        spatial_rows.append(row)
        present = reference_genes.intersection(spatial.var_names.astype(str))
        for gene in marker_union:
            intersection_rows.append({"sample": sample, "gene": gene, "present_in_reference": gene in reference_genes, "present_in_spatial": gene in spatial.var_names, "used_for_tangram": gene in present})
    pd.DataFrame(spatial_rows).to_csv(TABLE_DIR / "spatial_h5ad_manifest.csv", index=False)
    pd.DataFrame(intersection_rows).to_csv(TABLE_DIR / "marker_gene_intersection_by_sample.csv", index=False)

    completion = {
        "status": "completed",
        "n_subtypes": reference.n_obs,
        "n_epithelial_subtypes": int(reference.obs["cell_subtype"].astype(str).str.startswith("Epi_").sum()),
        "n_marker_rows": len(marker_table),
        "n_marker_union_genes": len(marker_union),
        "n_subtypes_with_fewer_than_100_markers": int(selection_summary["fewer_than_100_selected"].sum()),
        "n_reference_cells_total": int(cells_used["n_cells_used"].sum()),
        "all_reference_cells_used": bool(cells_used["all_cells_used"].all() and (cells_used["n_cells_used"] == cells_used["n_cells_total"]).all()),
        "n_spatial_samples": len(spatial_rows),
        "n_spatial_spots": int(sum(row["n_spots"] for row in spatial_rows)),
        "reference_h5ad": str(reference_path.resolve()),
        "code_file": str(Path(__file__).resolve()),
        "seed": SEED,
    }
    (TABLE_DIR / "reference_and_spatial_h5ad_completion.json").write_text(json.dumps(completion, indent=2) + "\n", encoding="utf-8")
    pd.DataFrame([{
        "code_file": str(Path(__file__).resolve()),
        "post_annotation_deg_root": str(DEG_ROOT.resolve()),
        "lineage_h5ad_root": str(LINEAGE_H5AD_ROOT.resolve()),
        "marker_count_per_subtype": TOP_MARKERS,
        "marker_filter": "pvals_adj < 0.05 and logfoldchanges > 0",
        "marker_order": "score descending; stable source order for ties",
        "reference_expression": "subtype mean of annotated single-cell adata.raw",
        "reference_cell_cap": "none",
        "spatial_expression": "curated filtered raw counts; normalize/log deferred to mapping",
        "seed": SEED,
    }]).to_csv(TABLE_DIR / "reference_build_parameters.csv", index=False)
    with (TABLE_DIR / "readme.txt").open("a", encoding="utf-8") as handle:
        handle.write("Reference build: 68 subtype pseudocells use all annotated cells and up to 100 significant positive post-annotation DEG markers per subtype (pvals_adj < 0.05, logfoldchanges > 0, score descending). Subtypes with fewer than 100 eligible unique genes retain all eligible genes and are recorded in marker_selection_summary.csv. Six per-sample h5ads contain curated raw counts plus barcode, array/pixel coordinates, sample, and Classification metadata.\n")
    with (TABLE_DIR / "package_versions.txt").open("a", encoding="utf-8") as handle:
        handle.write(f"reference_builder_python={sys.version.split()[0]}\nreference_builder_anndata={pkg('anndata')}\nreference_builder_code={Path(__file__).resolve()}\n")
    print(json.dumps(completion, indent=2), flush=True)


if __name__ == "__main__":
    main()
