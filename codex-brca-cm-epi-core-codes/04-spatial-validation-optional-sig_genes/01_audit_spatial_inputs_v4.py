#!/usr/bin/env python3
"""Final audit using the curated filtered raw-count matrices as canonical.

The three CID raw 10x matrices use a different gene-symbol annotation release
from their curated filtered matrices, so raw-vs-filtered count equality is only
tested when every filtered symbol can be aligned losslessly.
"""

from __future__ import annotations

import importlib.metadata
import importlib.util
import json
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 42
random.seed(SEED)
np.random.seed(SEED)
ROOT = Path("/mnt/disk18t/lr_xcy/riku/brca_cm-epi_codex")
WORKFLOW = ROOT / "epi-cm-core-workflow"
DATA = Path("/mnt/disk18t/lr_xcy/riku/external_atlas_for_cm-epi_validation/Breast_Wu2021_Zenodo4739739/spatial")
OUT = WORKFLOW / "tables/04-spatial-validation-optional-sig_genes/01-input-audit-and-reference-v4"
OUT.mkdir(parents=True, exist_ok=False)
SAMPLES = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]

HELPER_PATH = Path(__file__).with_name("01_audit_spatial_inputs_v3.py")
SPEC = importlib.util.spec_from_file_location("spatial_audit_v3_helpers", HELPER_PATH)
if SPEC is None or SPEC.loader is None:
    raise ImportError(HELPER_PATH)
HELPER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(HELPER)


def pkg(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def main() -> None:
    rows, classes, overrides, annotation_rows = [], [], [], []
    for sample in SAMPLES:
        fdir = DATA / "filtered_count_matrices" / f"{sample}_filtered_count_matrix"
        rdir = DATA / "raw_count_matrices" / f"{sample}_raw_feature_bc_matrix"
        meta_path = DATA / "metadata" / f"{sample}_metadata.csv"
        spatial_dir = DATA / "spatial" / f"{sample}_spatial"
        fm_path, rm_path = fdir / "matrix.mtx.gz", rdir / "matrix.mtx.gz"
        f_bars, r_bars = HELPER.lines(fdir / "barcodes.tsv.gz"), HELPER.lines(rdir / "barcodes.tsv.gz")
        f_names = HELPER.unique_matrix_symbols(HELPER.symbols(fdir / "features.tsv.gz"))
        r_names = HELPER.symbols(rdir / "features.tsv.gz")
        meta = pd.read_csv(meta_path, index_col=0)
        meta.index = meta.index.astype(str)
        pos = pd.read_csv(spatial_dir / "tissue_positions_list.csv", header=None, names=["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"])
        pos["barcode"] = pos["barcode"].astype(str)
        pos_idx = pos.set_index("barcode")
        f_m = HELPER.matrix(fm_path)
        filtered_totals = np.asarray(f_m.sum(axis=0)).ravel()
        metadata_totals = pd.to_numeric(meta.loc[f_bars, "nCount_RNA"], errors="coerce").to_numpy()

        raw_originals = set(r_names)
        f_keys = [HELPER.aggregate_key(name, raw_originals) for name in f_names]
        r_keys = [HELPER.aggregate_key(name, raw_originals) for name in r_names]
        missing_symbols = sorted(set(f_keys).difference(r_keys))
        exact_status = "not_applicable_different_gene_annotation_release"
        exact_nnz = np.nan
        exact_max = np.nan
        if not missing_symbols:
            r_m = HELPER.matrix(rm_path)
            r_bar_idx = {barcode: i for i, barcode in enumerate(r_bars)}
            r_m = r_m[:, [r_bar_idx[barcode] for barcode in f_bars]]
            _, diff_nnz, diff_max = HELPER.compare_by_symbol(f_m, f_names, r_m, r_names, sample)
            exact_status = "passed" if diff_nnz == 0 else "failed"
            exact_nnz = diff_nnz
            exact_max = diff_max
        annotation_rows.append({
            "sample": sample,
            "filtered_gene_rows": len(f_names),
            "raw_gene_rows": len(r_names),
            "filtered_symbols_absent_from_raw_annotation": len(missing_symbols),
            "example_absent_symbols": ";".join(missing_symbols[:20]),
            "raw_filtered_exact_overlap_status": exact_status,
            "raw_filtered_delta_nnz": exact_nnz,
            "raw_filtered_delta_max_abs": exact_max,
        })

        flags = pos_idx.loc[f_bars, "in_tissue"].astype(int).to_numpy()
        off = [barcode for barcode, flag in zip(f_bars, flags) if flag == 0]
        overrides.extend({"sample": sample, "barcode": barcode, "source_in_tissue": 0, "retained_because": "curated filtered+metadata spot with finite real coordinates"} for barcode in off)
        f_set, r_set, meta_set, pos_set = set(f_bars), set(r_bars), set(meta.index), set(pos["barcode"])
        h = HELPER.header(fm_path)
        rows.append({
            "sample": sample,
            "canonical_expression_input": "curated filtered raw-count matrix",
            "filtered_matrix_format": "gzip" if fm_path.open("rb").read(2) == b"\x1f\x8b" else "plain_text_despite_gz_suffix",
            "filtered_n_gene_rows": len(f_names), "filtered_n_spots": len(f_bars), "filtered_matrix_rows": h[0], "filtered_matrix_cols": h[1], "filtered_matrix_nnz": h[2],
            "raw_n_spots": len(r_bars), "metadata_rows": len(meta), "positions_rows": len(pos), "source_in_tissue_spots": int(pos["in_tissue"].astype(int).sum()), "curated_spots_with_source_in_tissue_0": len(off),
            "metadata_barcode_exact_filtered_order": f_bars == meta.index.tolist(), "filtered_barcode_set_eq_metadata": f_set == meta_set, "filtered_barcode_subset_raw": f_set.issubset(r_set), "filtered_barcode_subset_positions": f_set.issubset(pos_set),
            "all_curated_spots_have_finite_coordinates": bool(np.isfinite(pos_idx.loc[f_bars, ["array_row", "array_col", "pxl_row", "pxl_col"]].to_numpy(dtype=float)).all()),
            "metadata_nCount_matches_filtered": bool(np.allclose(filtered_totals, metadata_totals, rtol=0, atol=0)),
            "raw_filtered_exact_overlap_status": exact_status,
        })
        for label, count in meta["Classification"].fillna("NA").replace("", "NA").astype(str).value_counts().items():
            classes.append({"sample": sample, "Classification": label, "n_spots": int(count)})

    audit = pd.DataFrame(rows)
    annotation = pd.DataFrame(annotation_rows)
    audit.to_csv(OUT / "spatial_input_audit.csv", index=False)
    annotation.to_csv(OUT / "raw_filtered_gene_annotation_compatibility.csv", index=False)
    pd.DataFrame(classes).to_csv(OUT / "spot_classification_counts.csv", index=False)
    pd.DataFrame(overrides).to_csv(OUT / "curated_spots_with_source_in_tissue_0.csv", index=False)
    pd.DataFrame([{"sample": s, "status": "tumor", "status_label": "Primary tumor specimen", "include_all_samples": True, "include_tumor_only": True, "exclusion_reason_all_samples": "", "exclusion_reason_tumor_only": "", "status_source": "user_confirmed_primary_tumor_study_design; Classification is spot-level"} for s in SAMPLES]).to_csv(OUT / "sample_scope.csv", index=False)
    required = ["metadata_barcode_exact_filtered_order", "filtered_barcode_set_eq_metadata", "filtered_barcode_subset_raw", "filtered_barcode_subset_positions", "all_curated_spots_have_finite_coordinates", "metadata_nCount_matches_filtered"]
    failures = [{"sample": row["sample"], "criterion": key} for row in audit.to_dict("records") for key in required if not bool(row[key])]
    failures.extend({"sample": row["sample"], "criterion": "raw_filtered_exact_overlap_status"} for row in annotation.to_dict("records") if row["raw_filtered_exact_overlap_status"] == "failed")
    completion = {"status": "completed" if not failures else "failed", "n_samples": len(SAMPLES), "n_curated_spots": int(audit["filtered_n_spots"].sum()), "n_curated_spots_with_source_in_tissue_0": int(audit["curated_spots_with_source_in_tissue_0"].sum()), "n_exact_raw_overlap_samples": int(annotation["raw_filtered_exact_overlap_status"].eq("passed").sum()), "n_gene_annotation_mismatch_samples": int(annotation["raw_filtered_exact_overlap_status"].str.startswith("not_applicable").sum()), "n_failed_invariants": len(failures), "failures": failures, "canonical_input_is_filtered_raw_counts": True, "all_six_primary_tumor_specimens": True, "code_file": str(Path(__file__).resolve()), "seed": SEED}
    (OUT / "input_audit_completion.json").write_text(json.dumps(completion, indent=2) + "\n", encoding="utf-8")
    pd.DataFrame([{"code_file": str(Path(__file__).resolve()), "data_root": str(DATA), "samples": ",".join(SAMPLES), "canonical_expression_input": "filtered_count_matrices raw integer counts", "raw_matrix_role": "barcode/source audit; exact count comparison only when gene annotation is losslessly compatible", "spot_inclusion": "all curated filtered+metadata spots with finite real coordinates", "status": "all six primary tumor specimens", "seed": SEED}]).to_csv(OUT / "run_parameters.csv", index=False)
    (OUT / "readme.txt").write_text("Canonical expression input is the curated filtered raw-count matrix because it matches metadata spot order and nCount_RNA exactly. The source files are plain text despite .gz suffixes and are read by magic-byte detection without rewriting. The 114/116 raw matrices permit exact gene-symbol count overlap; the CID raw matrices use another gene-annotation release, documented separately, and are not mixed into expression. Curated source in_tissue=0 spots are retained only when they have metadata and finite real coordinates. All six sections are primary-tumor specimens; Classification is spot-level.\n", encoding="utf-8")
    (OUT / "package_versions.txt").write_text(f"python={sys.version.split()[0]}\nnumpy={pkg('numpy')}\npandas={pkg('pandas')}\nscipy={pkg('scipy')}\ncode={Path(__file__).resolve()}\nhelper={HELPER_PATH}\nseed={SEED}\n", encoding="utf-8")
    print(json.dumps(completion, indent=2), flush=True)
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
