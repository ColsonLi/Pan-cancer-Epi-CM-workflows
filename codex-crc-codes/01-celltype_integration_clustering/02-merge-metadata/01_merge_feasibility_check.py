#!/usr/bin/env python3
"""Pre-merge feasibility checks for selected CRC samples.

This script does not load expression matrices into AnnData. It inspects h5ad
metadata and CSR pointers directly, then writes the official merge feasibility
tables required before deciding whether a full unfiltered merge is practical.
"""

from __future__ import annotations

import platform
import random
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

SEED = 42
random.seed(SEED)
np.random.seed(SEED)

ANALYSIS_ROOT = Path("/mnt/disk18t/lr_xcy/riku/crc_val")
WORKFLOW_ROOT = ANALYSIS_ROOT / "epi-cm-core-workflow"
SELECTED_SAMPLES = (
    WORKFLOW_ROOT
    / "tables/01-celltype_integration_clustering/01-select-samples/02_eligible_samples_for_merge.csv"
)
TABLE_DIR = (
    WORKFLOW_ROOT / "tables/01-celltype_integration_clustering/02-merge-metadata"
)
CODE_FILE = Path(__file__)

INITIAL_MIN_GENE_CANDIDATES = [50, 100, 200]
CHUNK_ROWS = 5_000_000


def _decode_array(values) -> list[str]:
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values]


def _categorical_counts(obs_group, column: str) -> tuple[list[str], np.ndarray]:
    node = obs_group[column]
    if isinstance(node, h5py.Group) and "categories" in node and "codes" in node:
        categories = _decode_array(node["categories"][()])
        codes = node["codes"][()]
        return categories, codes
    values = np.asarray(_decode_array(node[()]))
    categories, codes = np.unique(values, return_inverse=True)
    return categories.tolist(), codes.astype(np.int64)


def _var_names(var_group) -> list[str]:
    if "_index" in var_group:
        return _decode_array(var_group["_index"][()])
    if "symbol" in var_group:
        return _decode_array(var_group["symbol"][()])
    raise KeyError("Could not find var/_index or var/symbol in h5ad")


def _csr_shape_and_nnz(x_node) -> tuple[int, int, int]:
    if isinstance(x_node, h5py.Group):
        shape = tuple(int(v) for v in x_node.attrs["shape"])
        nnz = int(x_node["data"].shape[0])
        return shape[0], shape[1], nnz
    shape = tuple(int(v) for v in x_node.shape)
    return shape[0], shape[1], int(np.prod(shape))


def _count_selected_rows_by_thresholds(
    h5ad_path: str,
    selected_sample_ids: set[str],
    thresholds: list[int],
) -> tuple[dict[str, int], dict[int, int], dict[int, int]]:
    with h5py.File(h5ad_path, "r") as f:
        categories, codes = _categorical_counts(f["obs"], "sample_id")
        selected_codes = np.array(
            [i for i, sample_id in enumerate(categories) if sample_id in selected_sample_ids],
            dtype=codes.dtype,
        )
        if selected_codes.size == 0:
            return {sid: 0 for sid in selected_sample_ids}, {t: 0 for t in thresholds}, {
                t: 0 for t in thresholds
            }

        per_sample = {}
        code_counts = np.bincount(codes[codes >= 0], minlength=len(categories))
        for sid in selected_sample_ids:
            per_sample[sid] = int(code_counts[categories.index(sid)]) if sid in categories else 0

        keep_by_threshold = {t: 0 for t in thresholds}
        nnz_by_threshold = {t: 0 for t in thresholds}
        n_obs, _, _ = _csr_shape_and_nnz(f["X"])

        if isinstance(f["X"], h5py.Group) and "indptr" in f["X"]:
            indptr_ds = f["X"]["indptr"]
            for start in range(0, n_obs, CHUNK_ROWS):
                end = min(start + CHUNK_ROWS, n_obs)
                code_chunk = codes[start:end]
                selected_mask = np.isin(code_chunk, selected_codes)
                if not selected_mask.any():
                    continue
                ptr = indptr_ds[start : end + 1]
                nnz_per_row = np.diff(ptr)
                for threshold in thresholds:
                    mask = selected_mask & (nnz_per_row >= threshold)
                    keep_by_threshold[threshold] += int(mask.sum())
                    nnz_by_threshold[threshold] += int(nnz_per_row[mask].sum())
        else:
            # Dense matrices are not expected here. Count selected rows only.
            for threshold in thresholds:
                keep_by_threshold[threshold] = int(np.isin(codes, selected_codes).sum())
                nnz_by_threshold[threshold] = np.nan

    return per_sample, keep_by_threshold, nnz_by_threshold


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    selected = pd.read_csv(SELECTED_SAMPLES)

    file_rows = []
    sample_rows = []
    common_genes: set[str] | None = None

    for h5ad_path, sub in selected.groupby("h5ad", sort=False):
        path = Path(h5ad_path)
        sample_ids = set(sub["sample_id"].astype(str))
        with h5py.File(path, "r") as f:
            n_obs, n_vars, nnz = _csr_shape_and_nnz(f["X"])
            genes = set(_var_names(f["var"]))
            common_genes = genes if common_genes is None else common_genes & genes
            matrix_encoding = f["X"].attrs.get("encoding-type", "dense")
            matrix_encoding = (
                matrix_encoding.decode("utf-8")
                if isinstance(matrix_encoding, bytes)
                else str(matrix_encoding)
            )

        per_sample, keep_by_threshold, nnz_by_threshold = _count_selected_rows_by_thresholds(
            h5ad_path, sample_ids, INITIAL_MIN_GENE_CANDIDATES
        )

        selected_rows = int(sum(per_sample.values()))
        file_row = {
            "h5ad": h5ad_path,
            "file": path.name,
            "dataset": ";".join(sorted(sub["dataset"].astype(str).unique())),
            "n_selected_samples": int(sub["sample_id"].nunique()),
            "n_obs_file": int(n_obs),
            "n_vars_file": int(n_vars),
            "nnz_file": int(nnz),
            "mean_nnz_per_obs_file": float(nnz / n_obs) if n_obs else np.nan,
            "n_selected_obs_unfiltered": selected_rows,
            "file_size_gb": path.stat().st_size / 1e9,
            "matrix_encoding": matrix_encoding,
        }
        for threshold in INITIAL_MIN_GENE_CANDIDATES:
            file_row[f"n_selected_obs_n_genes_ge_{threshold}"] = keep_by_threshold[threshold]
            file_row[f"selected_nnz_n_genes_ge_{threshold}"] = nnz_by_threshold[threshold]
        file_rows.append(file_row)

        for _, sample in sub.iterrows():
            row = sample.to_dict()
            sid = str(sample["sample_id"])
            row["source_h5ad"] = h5ad_path
            row["n_obs_unfiltered_by_h5ad_sample_id"] = per_sample.get(sid, 0)
            sample_rows.append(row)

    file_df = pd.DataFrame(file_rows).sort_values(
        "n_selected_obs_unfiltered", ascending=False
    )
    sample_df = pd.DataFrame(sample_rows).sort_values(
        ["dataset", "sample_id"], kind="mergesort"
    )

    common = sorted(common_genes or [])
    pd.DataFrame({"gene": common}).to_csv(
        TABLE_DIR / "inner_join_common_genes.csv", index=False
    )
    file_df.to_csv(TABLE_DIR / "merge_feasibility_by_h5ad.csv", index=False)
    sample_df.to_csv(TABLE_DIR / "selected_sample_manifest.csv", index=False)

    summary = pd.DataFrame(
        [
            {
                "step": "merge_feasibility_check",
                "selected_sample_csv": str(SELECTED_SAMPLES),
                "n_selected_samples": int(selected["sample_id"].nunique()),
                "n_selected_datasets": int(selected["dataset"].nunique()),
                "n_source_h5ad_files": int(selected["h5ad"].nunique()),
                "join_strategy": "inner",
                "n_common_genes": len(common),
                "n_selected_obs_unfiltered": int(
                    file_df["n_selected_obs_unfiltered"].sum()
                ),
                "largest_file": str(file_df.iloc[0]["file"]),
                "largest_file_selected_obs_unfiltered": int(
                    file_df.iloc[0]["n_selected_obs_unfiltered"]
                ),
                "initial_min_gene_candidates": ";".join(
                    map(str, INITIAL_MIN_GENE_CANDIDATES)
                ),
                "code_file": str(CODE_FILE),
                "random_seed": SEED,
                "notes": (
                    "This check uses h5ad obs/sample_id and CSR indptr only; "
                    "it does not load X into AnnData."
                ),
            }
        ]
    )
    for threshold in INITIAL_MIN_GENE_CANDIDATES:
        summary.loc[0, f"n_selected_obs_n_genes_ge_{threshold}"] = int(
            file_df[f"n_selected_obs_n_genes_ge_{threshold}"].sum()
        )
        summary.loc[0, f"selected_nnz_n_genes_ge_{threshold}"] = int(
            file_df[f"selected_nnz_n_genes_ge_{threshold}"].sum()
        )
    summary.to_csv(TABLE_DIR / "merge_feasibility_summary.csv", index=False)

    (TABLE_DIR / "readme.txt").write_text(
        "\n".join(
            [
                "02-merge-metadata feasibility check",
                "",
                f"Input selected samples: {SELECTED_SAMPLES}",
                f"Code: {CODE_FILE}",
                "Join strategy planned by skill default: inner gene join.",
                "Outputs:",
                "- merge_feasibility_by_h5ad.csv",
                "- selected_sample_manifest.csv",
                "- inner_join_common_genes.csv",
                "- merge_feasibility_summary.csv",
                "- package_versions.txt",
                "",
                "No expression matrix was loaded or merged in this feasibility pass.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    versions = [
        f"python={platform.python_version()}",
        f"platform={platform.platform()}",
        f"numpy={np.__version__}",
        f"pandas={pd.__version__}",
        f"h5py={h5py.__version__}",
        f"environment=/mnt/disk18t/lr_xcy/riku/crc_val/uv_envs/main/.venv",
        f"code_file={CODE_FILE}",
        f"random_seed={SEED}",
    ]
    (TABLE_DIR / "package_versions.txt").write_text(
        "\n".join(versions) + "\n", encoding="utf-8"
    )

    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
