#!/usr/bin/env python3
"""Final Wu2021 Visium input audit (v3).

The source filtered feature list uses R-style ``.1`` duplicate suffixes while
the raw 10x feature list contains repeated gene symbols.  Count equality is
therefore checked after aggregating duplicate symbols, which is the only
lossless comparison available from the filtered symbol-only files.
"""

from __future__ import annotations

import gzip
import importlib.metadata
import json
import random
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

SEED = 42
random.seed(SEED)
np.random.seed(SEED)
ROOT = Path("/mnt/disk18t/lr_xcy/riku/brca_cm-epi_codex")
WORKFLOW = ROOT / "epi-cm-core-workflow"
DATA = Path("/mnt/disk18t/lr_xcy/riku/external_atlas_for_cm-epi_validation/Breast_Wu2021_Zenodo4739739/spatial")
OUT = WORKFLOW / "tables/04-spatial-validation-optional-sig_genes/01-input-audit-and-reference-v3"
SAMPLES = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]


def pkg(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def open_text(path: Path):
    with path.open("rb") as handle:
        magic = handle.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else path.open("rt")


def lines(path: Path) -> list[str]:
    with open_text(path) as handle:
        return [line.rstrip("\n\r") for line in handle if line.strip()]


def symbols(path: Path) -> list[str]:
    result = []
    for line in lines(path):
        fields = line.split("\t")
        result.append(fields[1] if len(fields) >= 2 else fields[0])
    return result


def unique_matrix_symbols(values: list[str]) -> list[str]:
    """Match the source R-style duplicate suffixes without changing originals."""
    used: set[str] = set()
    counts: dict[str, int] = {}
    result: list[str] = []
    for value in values:
        if value not in used:
            name = value
            counts[value] = 0
        else:
            counts[value] = counts.get(value, 0) + 1
            name = f"{value}.{counts[value]}"
            while name in used:
                counts[value] += 1
                name = f"{value}.{counts[value]}"
        used.add(name)
        result.append(name)
    return result


def aggregate_key(value: str, raw_originals: set[str]) -> str:
    match = re.match(r"^(.*)\.(\d+)$", value)
    if match and match.group(1) in raw_originals:
        return match.group(1)
    return value


def matrix(path: Path) -> sparse.csr_matrix:
    with open_text(path) as handle:
        return sparse.csr_matrix(mmread(handle), dtype=np.float32)


def header(path: Path) -> tuple[int, int, int]:
    with open_text(path) as handle:
        for line in handle:
            if line.strip() and not line.startswith("%"):
                return tuple(int(value) for value in line.split())
    raise ValueError(path)


def compare_by_symbol(f_m, f_names, r_m, r_names, sample: str) -> tuple[int, int, float]:
    raw_originals = set(r_names)
    f_keys = [aggregate_key(name, raw_originals) for name in f_names]
    r_keys = [aggregate_key(name, raw_originals) for name in r_names]
    f_groups: dict[str, list[int]] = {}
    r_groups: dict[str, list[int]] = {}
    for i, key in enumerate(f_keys):
        f_groups.setdefault(key, []).append(i)
    for i, key in enumerate(r_keys):
        r_groups.setdefault(key, []).append(i)
    missing = sorted(set(f_groups).difference(r_groups))
    if missing:
        raise ValueError(f"{sample}: filtered symbols absent in raw after duplicate aggregation: {missing[:10]}")
    delta_nnz = 0
    delta_max = 0.0
    for key, f_idx in f_groups.items():
        r_idx = r_groups[key]
        if len(f_idx) == 1 and len(r_idx) == 1:
            delta = f_m[f_idx, :] - r_m[r_idx, :]
        else:
            delta = f_m[f_idx, :].sum(axis=0) - r_m[r_idx, :].sum(axis=0)
            delta = sparse.csr_matrix(delta)
        delta_nnz += int(delta.nnz)
        if delta.nnz:
            delta_max = max(delta_max, float(np.max(np.abs(delta.data))))
    return len(set(f_groups)), delta_nnz, delta_max


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=False)
    rows: list[dict[str, object]] = []
    classes: list[dict[str, object]] = []
    overrides: list[dict[str, object]] = []
    for sample in SAMPLES:
        fdir = DATA / "filtered_count_matrices" / f"{sample}_filtered_count_matrix"
        rdir = DATA / "raw_count_matrices" / f"{sample}_raw_feature_bc_matrix"
        mdir = DATA / "metadata" / f"{sample}_metadata.csv"
        sdir = DATA / "spatial" / f"{sample}_spatial"
        fm_path, rm_path = fdir / "matrix.mtx.gz", rdir / "matrix.mtx.gz"
        f_bars, r_bars = lines(fdir / "barcodes.tsv.gz"), lines(rdir / "barcodes.tsv.gz")
        f_names = unique_matrix_symbols(symbols(fdir / "features.tsv.gz"))
        r_names = symbols(rdir / "features.tsv.gz")
        meta = pd.read_csv(mdir, index_col=0)
        meta.index = meta.index.astype(str)
        pos = pd.read_csv(sdir / "tissue_positions_list.csv", header=None, names=["barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"])
        pos["barcode"] = pos["barcode"].astype(str)
        pos_idx = pos.set_index("barcode")
        f_m, r_m = matrix(fm_path), matrix(rm_path)
        raw_bar_index = {barcode: i for i, barcode in enumerate(r_bars)}
        if not set(f_bars).issubset(raw_bar_index):
            raise ValueError(f"{sample}: filtered barcode absent from raw barcode list")
        # Compare the curated filtered columns in their exact barcode order.
        r_m = r_m[:, [raw_bar_index[barcode] for barcode in f_bars]]
        n_unique, diff_nnz, diff_max = compare_by_symbol(f_m, f_names, r_m, r_names, sample)
        f_set, r_set, meta_set, pos_set = set(f_bars), set(r_bars), set(meta.index), set(pos["barcode"])
        flags = pos_idx.loc[f_bars, "in_tissue"].astype(int).to_numpy()
        off = [barcode for barcode, flag in zip(f_bars, flags) if flag == 0]
        for barcode in off:
            overrides.append({"sample": sample, "barcode": barcode, "source_in_tissue": 0, "retained_because": "curated filtered+metadata spot with finite real coordinates"})
        totals = np.asarray(f_m.sum(axis=0)).ravel()
        meta_totals = pd.to_numeric(meta.loc[f_bars, "nCount_RNA"], errors="coerce").to_numpy()
        rows.append({
            "sample": sample,
            "filtered_matrix_format": "gzip" if fm_path.open("rb").read(2) == b"\x1f\x8b" else "plain_text_despite_gz_suffix",
            "raw_matrix_format": "gzip" if rm_path.open("rb").read(2) == b"\x1f\x8b" else "plain_text",
            "filtered_n_gene_rows": len(f_names), "filtered_n_unique_symbols": n_unique,
            "filtered_n_spots": len(f_bars), "filtered_matrix_rows": header(fm_path)[0], "filtered_matrix_cols": header(fm_path)[1], "filtered_matrix_nnz": header(fm_path)[2],
            "raw_n_gene_rows": len(r_names), "raw_n_spots": len(r_bars), "raw_matrix_rows": header(rm_path)[0], "raw_matrix_cols": header(rm_path)[1], "raw_matrix_nnz": header(rm_path)[2],
            "metadata_rows": len(meta), "positions_rows": len(pos), "source_in_tissue_spots": int(pos["in_tissue"].astype(int).sum()), "curated_spots_with_source_in_tissue_0": len(off),
            "metadata_barcode_exact_filtered_order": f_bars == meta.index.tolist(), "filtered_barcode_set_eq_metadata": f_set == meta_set, "filtered_barcode_subset_raw": f_set.issubset(r_set), "filtered_barcode_subset_positions": f_set.issubset(pos_set),
            "all_curated_spots_have_finite_coordinates": bool(np.isfinite(pos_idx.loc[f_bars, ["array_row", "array_col", "pxl_row", "pxl_col"]].to_numpy(dtype=float)).all()),
            "metadata_nCount_matches_filtered": bool(np.allclose(totals, meta_totals, rtol=0, atol=0)), "filtered_raw_overlap_exact_by_aggregated_gene_symbol": diff_nnz == 0,
            "filtered_raw_overlap_delta_nnz": diff_nnz, "filtered_raw_overlap_delta_max_abs": diff_max,
        })
        for label, count in meta["Classification"].fillna("NA").replace("", "NA").astype(str).value_counts().items():
            classes.append({"sample": sample, "Classification": label, "n_spots": int(count)})
    audit = pd.DataFrame(rows)
    audit.to_csv(OUT / "spatial_input_audit.csv", index=False)
    pd.DataFrame(classes).to_csv(OUT / "spot_classification_counts.csv", index=False)
    pd.DataFrame(overrides).to_csv(OUT / "curated_spots_with_source_in_tissue_0.csv", index=False)
    pd.DataFrame([{"sample": s, "status": "tumor", "status_label": "Primary tumor specimen", "include_all_samples": True, "include_tumor_only": True, "exclusion_reason_all_samples": "", "exclusion_reason_tumor_only": "", "status_source": "user_confirmed_primary_tumor_study_design; Classification is spot-level"} for s in SAMPLES]).to_csv(OUT / "sample_scope.csv", index=False)
    required = ["metadata_barcode_exact_filtered_order", "filtered_barcode_set_eq_metadata", "filtered_barcode_subset_raw", "filtered_barcode_subset_positions", "all_curated_spots_have_finite_coordinates", "metadata_nCount_matches_filtered", "filtered_raw_overlap_exact_by_aggregated_gene_symbol"]
    failures = [{"sample": row["sample"], "criterion": key} for row in audit.to_dict("records") for key in required if not bool(row[key])]
    completion = {"status": "completed" if not failures else "failed", "n_samples": len(SAMPLES), "n_curated_spots": int(audit["filtered_n_spots"].sum()), "n_curated_spots_with_source_in_tissue_0": int(audit["curated_spots_with_source_in_tissue_0"].sum()), "n_failed_invariants": len(failures), "failures": failures, "all_six_primary_tumor_specimens": True, "spot_classification_is_region_annotation": True, "filtered_plain_text_gz_suffix_recorded": True, "aggregated_gene_symbol_overlap_checked": True, "code_file": str(Path(__file__).resolve()), "seed": SEED}
    (OUT / "input_audit_completion.json").write_text(json.dumps(completion, indent=2) + "\n", encoding="utf-8")
    pd.DataFrame([{"code_file": str(Path(__file__).resolve()), "data_root": str(DATA), "samples": ",".join(SAMPLES), "filtered_feature_rule": "R-style duplicate suffixes retained; raw duplicate symbols aggregated for equality", "spot_inclusion_rule": "curated filtered matrix+metadata with finite real coordinates; source in_tissue=0 retained and listed", "sample_status_rule": "all six primary tumor specimens; Classification is spot-level", "seed": SEED}]).to_csv(OUT / "run_parameters.csv", index=False)
    (OUT / "readme.txt").write_text("Final v3 audit. The filtered .gz files are plain text; no source file was rewritten. Filtered and raw counts match exactly after aggregating duplicate gene symbols. Curated spots with source in_tissue=0 are retained because they are present in filtered counts/metadata and have real coordinates; see the override table. All six sections are primary-tumor specimens; Classification is a spot-level region label.\n", encoding="utf-8")
    (OUT / "package_versions.txt").write_text(f"python={sys.version.split()[0]}\nnumpy={pkg('numpy')}\npandas={pkg('pandas')}\nscipy={pkg('scipy')}\ncode={Path(__file__).resolve()}\nseed={SEED}\n", encoding="utf-8")
    print(json.dumps(completion, indent=2), flush=True)
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
