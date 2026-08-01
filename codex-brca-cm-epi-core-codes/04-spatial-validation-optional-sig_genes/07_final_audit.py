#!/usr/bin/env python3
"""Final numerical, inventory, format, and provenance audit for Module 04."""

from __future__ import annotations

import importlib.metadata
import json
import random
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

SEED = 42
random.seed(SEED)
np.random.seed(SEED)
ROOT = Path("/mnt/disk18t/lr_xcy/riku/brca_cm-epi_codex")
WORKFLOW = ROOT / "epi-cm-core-workflow"
MODULE = "04-spatial-validation-optional-sig_genes"
TABLE_MODULE = WORKFLOW / "tables" / MODULE
FIGURE_MODULE = WORKFLOW / "figures" / MODULE
H5AD_MODULE = WORKFLOW / "h5ad" / MODULE
CODE_MODULE = WORKFLOW / "codes" / MODULE
MAPPING_DIR = TABLE_MODULE / "02-tangram-mapping"
STATS_DIR = TABLE_MODULE / "03-spatial-statistics-and-plotting"
FIGURE_STATS = FIGURE_MODULE / "03-spatial-statistics-and-plotting"
OUT = TABLE_MODULE / "04-final-audit"
OUT.mkdir(parents=True, exist_ok=False)
SAMPLES = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]
RASTER_SUFFIXES = {".png", ".jpg", ".jpeg", ".tif", ".tiff", ".webp", ".bmp"}


def pkg(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def signature_ok(path: Path) -> bool:
    size = path.stat().st_size
    if size <= 0:
        return False
    with path.open("rb") as handle:
        head = handle.read(1024)
        handle.seek(max(0, size - 2048))
        tail = handle.read()
    if path.suffix.lower() == ".pdf":
        return head.startswith(b"%PDF") and b"%%EOF" in tail
    if path.suffix.lower() == ".svg":
        return b"<svg" in head and b"</svg>" in tail
    return True


def add(rows: list[dict[str, object]], criterion: str, observed, expected, passed: bool, evidence: str) -> None:
    rows.append({"criterion": criterion, "observed": observed, "expected": expected, "passed": bool(passed), "evidence": evidence})


def audit_marker_contract(markers: pd.DataFrame) -> tuple[pd.DataFrame, bool]:
    required = {"subtype", "gene", "marker_rank", "pvals_adj", "score", "logfoldchanges", "source_deg_csv"}
    if not required.issubset(markers.columns):
        missing = sorted(required.difference(markers.columns))
        return pd.DataFrame([{"subtype": "ALL", "n_observed": len(markers), "n_expected": np.nan, "passed": False, "reason": f"missing columns: {missing}"}]), False
    audit_rows = []
    for subtype, observed in markers.groupby("subtype", sort=True):
        observed = observed.sort_values("marker_rank", kind="stable").copy()
        sources = observed["source_deg_csv"].astype(str).unique().tolist()
        if len(sources) != 1 or not Path(sources[0]).exists():
            audit_rows.append({"subtype": subtype, "n_observed": len(observed), "n_expected": np.nan, "passed": False, "reason": "source DEG path is missing or non-unique"})
            continue
        source = pd.read_csv(sources[0])
        source["_source_order"] = np.arange(len(source))
        source["pvals_adj"] = pd.to_numeric(source["pvals_adj"], errors="raise")
        source["score"] = pd.to_numeric(source["score"], errors="raise")
        source["logfoldchanges"] = pd.to_numeric(source["logfoldchanges"], errors="raise")
        expected = source.loc[source["pvals_adj"].lt(0.05) & source["logfoldchanges"].gt(0)].copy()
        expected = expected.sort_values(["score", "_source_order"], ascending=[False, True], kind="stable").drop_duplicates("gene", keep="first").head(100)
        gene_match = observed["gene"].astype(str).tolist() == expected["gene"].astype(str).tolist()
        rank_match = observed["marker_rank"].astype(int).tolist() == list(range(1, len(observed) + 1))
        values_match = len(observed) == len(expected)
        for column in ["pvals_adj", "score", "logfoldchanges"]:
            values_match = values_match and np.allclose(observed[column].to_numpy(float), expected[column].to_numpy(float), rtol=1e-12, atol=1e-12)
        passed = bool(gene_match and rank_match and values_match)
        audit_rows.append({
            "subtype": subtype,
            "n_observed": len(observed),
            "n_expected": len(expected),
            "max_pvals_adj": float(observed["pvals_adj"].max()),
            "min_logfoldchanges": float(observed["logfoldchanges"].min()),
            "gene_order_exact": gene_match,
            "rank_exact": rank_match,
            "numeric_values_exact": values_match,
            "passed": passed,
            "reason": "" if passed else "selected markers differ from independent source-table recomputation",
        })
    report = pd.DataFrame(audit_rows)
    return report, len(report) == 68 and report["passed"].astype(bool).all()


def main() -> None:
    rows: list[dict[str, object]] = []
    inventory_rows: list[dict[str, object]] = []

    input_completion = json.loads((TABLE_MODULE / "01-input-audit-and-reference-v4/input_audit_completion.json").read_text())
    reference_completion = json.loads((TABLE_MODULE / "01-input-audit-and-reference-v4/reference_and_spatial_h5ad_completion.json").read_text())
    mapping_completion = json.loads((MAPPING_DIR / "mapping_all_samples_completion.json").read_text())
    stats_completion = json.loads((STATS_DIR / "statistics_and_heatmaps_completion.json").read_text())
    for name, value in [("input audit", input_completion), ("reference build", reference_completion), ("Tangram mapping", mapping_completion), ("statistics/heatmaps", stats_completion)]:
        add(rows, f"{name} completion status", value.get("status"), "completed", value.get("status") == "completed", "completion JSON")

    reference_path = H5AD_MODULE / "01-input-audit-and-reference/subtype_pseudobulk_reference_top100_significant_positive_deg.h5ad"
    reference = ad.read_h5ad(reference_path, backed="r")
    add(rows, "reference subtype count", reference.n_obs, 68, reference.n_obs == 68, str(reference_path))
    add(rows, "reference epithelial subtype count", int(reference.obs["cell_subtype"].astype(str).str.startswith("Epi_").sum()), 10, int(reference.obs["cell_subtype"].astype(str).str.startswith("Epi_").sum()) == 10, str(reference_path))
    reference.file.close()
    cells_used = pd.read_csv(TABLE_MODULE / "01-input-audit-and-reference-v4/reference_cells_used.csv")
    add(rows, "reference cells used", int(cells_used["n_cells_used"].sum()), 80406, int(cells_used["n_cells_used"].sum()) == 80406 and cells_used["all_cells_used"].astype(bool).all() and (cells_used["n_cells_used"] == cells_used["n_cells_total"]).all(), "reference_cells_used.csv")
    markers = pd.read_csv(TABLE_MODULE / "01-input-audit-and-reference-v4/marker_genes_used.csv")
    marker_contract, marker_contract_ok = audit_marker_contract(markers)
    marker_contract.to_csv(OUT / "marker_deg_contract_reaudit.csv", index=False)
    expected_marker_rows = int(pd.to_numeric(marker_contract["n_expected"], errors="coerce").sum())
    add(rows, "post-annotation significant marker rows", len(markers), expected_marker_rows, len(markers) == expected_marker_rows, "marker_genes_used.csv and independent source-table recomputation")
    add(rows, "marker DEG contract", marker_contract_ok, True, marker_contract_ok, "marker_deg_contract_reaudit.csv; pvals_adj < 0.05, logfoldchanges > 0, score descending, unique top 100 or all available")
    selection_summary = pd.read_csv(TABLE_MODULE / "01-input-audit-and-reference-v4/marker_selection_summary.csv")
    summary_ok = len(selection_summary) == 68 and selection_summary["n_selected"].between(1, 100).all() and selection_summary["n_selected"].sum() == len(markers)
    add(rows, "marker selection summary", f"{len(selection_summary)} subtypes; {int(selection_summary['fewer_than_100_selected'].astype(bool).sum())} below 100", "68 subtypes; fewer-than-100 explicitly recorded", summary_ok, "marker_selection_summary.csv")

    total_spots = 0
    score_paths = []
    for sample in SAMPLES:
        score_path = MAPPING_DIR / f"{sample}_tangram_pseudobulk_epi_cm_spot_scores.csv"
        score_paths.append(score_path)
        score = pd.read_csv(score_path)
        epi_cols = [column for column in score.columns if column.startswith("EPIfrac__")]
        cm_cols = [column for column in score.columns if column.startswith("CMact__")]
        total_spots += len(score)
        add(rows, f"{sample} subtype/CM score dimensions", f"{len(epi_cols)} Epi, {len(cm_cols)} CM", "10 Epi, 10 CM", len(epi_cols) == 10 and len(cm_cols) == 10, str(score_path))
        epi_values = score[epi_cols].to_numpy(float)
        cm_values = score[cm_cols].to_numpy(float)
        passed = np.isfinite(epi_values).all() and np.isfinite(cm_values).all() and (epi_values >= -1e-8).all() and (cm_values >= -1e-8).all() and np.allclose(epi_values.sum(axis=1), 1.0, rtol=1e-6, atol=1e-8)
        add(rows, f"{sample} score numerical invariants", passed, True, bool(passed), "finite/nonnegative CMact and EPIfrac; EPIfrac rows sum to one")
        coordinate_ok = np.isfinite(score[["array_row", "array_col", "spatial_x", "spatial_y"]].to_numpy(float)).all()
        add(rows, f"{sample} spot ID and coordinate invariants", f"unique={score['spot_id'].astype(str).is_unique}; finite={coordinate_ok}", "unique=True; finite=True", score["spot_id"].astype(str).is_unique and coordinate_ok, str(score_path))
    add(rows, "total spot-score rows", total_spots, 15611, total_spots == 15611, "six canonical spot-score CSVs")
    add(rows, "spot-score CSV count", len(score_paths), 6, len(score_paths) == 6 and all(path.exists() for path in score_paths), str(MAPPING_DIR))
    mapper_paths = sorted((H5AD_MODULE / "02-tangram-mapping").glob("*_tangram_mapper.h5ad"))
    add(rows, "Tangram mapper h5ad count", len(mapper_paths), 6, len(mapper_paths) == 6, str(H5AD_MODULE / "02-tangram-mapping"))

    pair_stats = pd.read_csv(STATS_DIR / "per_sample_all_pair_statistics.csv")
    pair_manifest = pd.read_csv(STATS_DIR / "all_sample_cm_epi_pair_manifest.csv")
    add(rows, "sample x CM x Epi statistics", len(pair_stats), 600, len(pair_stats) == 600 and pair_stats.groupby("sample").size().eq(100).all(), "per_sample_all_pair_statistics.csv")
    add(rows, "all-pair figure manifest", len(pair_manifest), 600, len(pair_manifest) == 600 and pair_manifest[["sample", "CM", "epi_subtype"]].drop_duplicates().shape[0] == 600, "all_sample_cm_epi_pair_manifest.csv")
    for column in ["raw_pdf", "raw_svg", "percentile_pdf", "percentile_svg"]:
        paths = [Path(value) for value in pair_manifest[column].astype(str)]
        ok = len(paths) == 600 and len(set(paths)) == 600 and all(path.exists() and signature_ok(path) for path in paths)
        add(rows, f"manifest {column} files", len(paths), 600, ok, column)

    all_stouffer = pd.read_csv(STATS_DIR / "statistics/all-samples/percentile_quadrant_fisher_sample_stouffer_all_samples.csv")
    tumor_stouffer = pd.read_csv(STATS_DIR / "statistics/tumor-only/percentile_quadrant_fisher_sample_stouffer_tumor_only.csv")
    add(rows, "all-samples Stouffer pair count", len(all_stouffer), 100, len(all_stouffer) == 100 and all_stouffer["n_samples"].eq(6).all(), "all-samples Stouffer CSV")
    add(rows, "tumor-only Stouffer pair count", len(tumor_stouffer), 100, len(tumor_stouffer) == 100 and tumor_stouffer["n_samples"].eq(6).all(), "tumor-only Stouffer CSV")
    scope_equal = all_stouffer[["CM", "epi_subtype", "combined_signed_z", "p_value", "q_value_bh"]].reset_index(drop=True).equals(tumor_stouffer[["CM", "epi_subtype", "combined_signed_z", "p_value", "q_value_bh"]].reset_index(drop=True))
    add(rows, "scope results identical for six tumor specimens", scope_equal, True, scope_equal, "explicit sample_scope.csv includes all six in both scopes")

    figure_files = sorted(path for path in FIGURE_MODULE.rglob("*") if path.is_file())
    pdf_files = [path for path in figure_files if path.suffix.lower() == ".pdf"]
    svg_files = [path for path in figure_files if path.suffix.lower() == ".svg"]
    raster_files = [path for path in figure_files if path.suffix.lower() in RASTER_SUFFIXES]
    add(rows, "sig_genes Module04 PDF count", len(pdf_files), 1202, len(pdf_files) == 1202, str(FIGURE_MODULE))
    add(rows, "sig_genes Module04 SVG count", len(svg_files), 1202, len(svg_files) == 1202, str(FIGURE_MODULE))
    add(rows, "sig_genes Module04 raster result files", len(raster_files), 0, len(raster_files) == 0, str(FIGURE_MODULE))
    signatures = all(signature_ok(path) for path in pdf_files + svg_files)
    add(rows, "all PDF/SVG signatures and EOF markers", signatures, True, signatures, f"{len(pdf_files)+len(svg_files)} files")

    all_module_files = [path for parent in [TABLE_MODULE, FIGURE_MODULE, H5AD_MODULE, CODE_MODULE] for path in parent.rglob("*") if path.is_file()]
    forbidden = [path for path in all_module_files if path.suffix.lower() in RASTER_SUFFIXES]
    add(rows, "no raster outputs in any sig_genes Module04 category", len(forbidden), 0, len(forbidden) == 0, "tables/figures/h5ad/codes")
    for path in all_module_files:
        category = path.relative_to(WORKFLOW).parts[0]
        inventory_rows.append({"category": category, "relative_path": str(path.relative_to(WORKFLOW)), "bytes": path.stat().st_size, "suffix": path.suffix.lower()})
    pd.DataFrame(inventory_rows).to_csv(OUT / "module_file_inventory.csv", index=False)
    audit = pd.DataFrame(rows)
    audit.to_csv(OUT / "validation_report.csv", index=False)
    failures = audit.loc[~audit["passed"].astype(bool)].to_dict("records")
    completion = {"status": "completed" if not failures else "failed", "n_validation_criteria": len(audit), "n_failed_criteria": len(failures), "failures": failures, "n_module_files": len(inventory_rows), "n_pdf": len(pdf_files), "n_svg": len(svg_files), "n_raster": len(raster_files), "n_spots": total_spots, "n_sample_pair_rows": len(pair_stats), "n_pairs_per_scope": len(all_stouffer), "code_file": str(Path(__file__).resolve()), "seed": SEED}
    (OUT / "final_audit_completion.json").write_text(json.dumps(completion, indent=2) + "\n", encoding="utf-8")
    pd.DataFrame([{"code_file": str(Path(__file__).resolve()), "module": MODULE, "marker_contract": "pvals_adj < 0.05 and logfoldchanges > 0; score descending; unique top 100 or all available", "expected_samples": 6, "expected_spots": 15611, "expected_sample_pair_rows": 600, "expected_pdf": 1202, "expected_svg": 1202, "expected_raster": 0, "signature_checks": "PDF header+EOF; SVG open+close tags", "seed": SEED}]).to_csv(OUT / "run_parameters.csv", index=False)
    (OUT / "readme.txt").write_text("Final sig_genes Module04 audit independently reconstructs every subtype marker list from its saved post-annotation DEG CSV under pvals_adj < 0.05 and logfoldchanges > 0, score-descending unique top-100/all-available rules. It also covers all-cell pseudobulk provenance, six score CSVs and mapper h5ads, EPIfrac/CMact numerical invariants, the complete 600-row sample-pair family, both 100-pair Stouffer scopes, all manifest figure paths, PDF/SVG signatures, and zero raster result files. module_file_inventory.csv is the checkable file inventory.\n", encoding="utf-8")
    (OUT / "package_versions.txt").write_text(f"python={sys.version.split()[0]}\nanndata={pkg('anndata')}\nnumpy={pkg('numpy')}\npandas={pkg('pandas')}\ncode={Path(__file__).resolve()}\nseed={SEED}\n", encoding="utf-8")
    print(json.dumps(completion, indent=2), flush=True)
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
