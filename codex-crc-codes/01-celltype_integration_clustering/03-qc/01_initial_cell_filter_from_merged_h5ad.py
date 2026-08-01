#!/usr/bin/env python3
"""Initial cell filter for the huge pre-QC merged h5ad.

The skill requires saving the unfiltered merge before QC. That object is too
large for AnnData/Scanpy to load directly, so this script applies only the
first QC step on disk: keep cells with at least MIN_GENES_INITIAL detected
genes. All later QC steps operate on the resulting AnnData-readable h5ad.
"""

from __future__ import annotations

import os
import platform
import random
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

SEED = 42
random.seed(SEED)
np.random.seed(SEED)

ANALYSIS_ROOT = Path("/mnt/disk18t/lr_xcy/riku/crc_val")
WORKFLOW_ROOT = ANALYSIS_ROOT / "epi-cm-core-workflow"
INPUT_H5AD = (
    WORKFLOW_ROOT
    / "h5ad/01-celltype_integration_clustering/02-merge-metadata/adata_merge.h5ad"
)
TABLE_DIR = WORKFLOW_ROOT / "tables/01-celltype_integration_clustering/03-qc"
MERGE_TABLE_DIR = WORKFLOW_ROOT / "tables/01-celltype_integration_clustering/02-merge-metadata"
H5AD_DIR = WORKFLOW_ROOT / "h5ad/01-celltype_integration_clustering/03-qc"
OUTPUT_H5AD = H5AD_DIR / "adata_initial_cell_filtered.h5ad"
PART_H5AD = H5AD_DIR / "adata_initial_cell_filtered.h5ad.part"
CODE_FILE = Path(__file__)

MIN_GENES_INITIAL = 200
SCAN_ROWS = 100_000
SPARSE_CHUNK = 8_000_000


def _now() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def _copy_root_schema(src: h5py.File, out_path: Path, n_obs: int) -> h5py.File:
    if OUTPUT_H5AD.exists() or PART_H5AD.exists():
        raise FileExistsError(
            f"Output already exists: {OUTPUT_H5AD} or {PART_H5AD}. "
            "Cleanup must be explicit before rerun."
        )
    out = h5py.File(PART_H5AD, "w")
    out.attrs.update(src.attrs)

    obs = out.create_group("obs")
    obs.attrs.update(src["obs"].attrs)
    old_cols = list(src["obs"].attrs["column-order"])
    old_cols = [c.decode("utf-8") if isinstance(c, bytes) else str(c) for c in old_cols]
    new_cols = old_cols + ["n_genes_by_counts", "total_counts"]
    obs.attrs["column-order"] = np.asarray(new_cols, dtype=object)

    for key, node in src["obs"].items():
        if isinstance(node, h5py.Group):
            group = obs.create_group(key)
            group.attrs.update(node.attrs)
            group.create_dataset(
                "categories",
                data=node["categories"][()],
                dtype=h5py.string_dtype(encoding="utf-8"),
            )
            group["categories"].attrs.update(node["categories"].attrs)
            group.create_dataset(
                "codes",
                shape=(n_obs,),
                maxshape=(n_obs,),
                chunks=(min(100_000, max(n_obs, 1)),),
                dtype=node["codes"].dtype,
                compression="lzf",
            )
            group["codes"].attrs.update(node["codes"].attrs)
        else:
            dtype = node.dtype
            if dtype.kind == "O":
                dtype = h5py.string_dtype(encoding="utf-8")
            obs.create_dataset(
                key,
                shape=(n_obs,),
                maxshape=(n_obs,),
                chunks=(min(100_000, max(n_obs, 1)),),
                dtype=dtype,
                compression="lzf" if dtype != h5py.string_dtype(encoding="utf-8") else None,
            )
            obs[key].attrs.update(node.attrs)

    for metric in ["n_genes_by_counts", "total_counts"]:
        obs.create_dataset(
            metric,
            shape=(n_obs,),
            maxshape=(n_obs,),
            chunks=(min(100_000, max(n_obs, 1)),),
            dtype=np.float64 if metric == "total_counts" else np.int32,
            compression="lzf",
        )
        obs[metric].attrs["encoding-type"] = "array"
        obs[metric].attrs["encoding-version"] = "0.2.0"

    src.copy("var", out)
    for name in ["layers", "obsm", "obsp", "uns", "varm", "varp"]:
        group = out.create_group(name)
        group.attrs["encoding-type"] = "dict"
        group.attrs["encoding-version"] = "0.1.0"

    x = out.create_group("X")
    x.attrs.update(src["X"].attrs)
    x.attrs["shape"] = np.asarray([n_obs, int(src["X"].attrs["shape"][1])], dtype=np.int64)
    x.create_dataset(
        "data",
        shape=(0,),
        maxshape=(None,),
        chunks=(SPARSE_CHUNK,),
        dtype=np.int32,
        compression="lzf",
    )
    x.create_dataset(
        "indices",
        shape=(0,),
        maxshape=(None,),
        chunks=(SPARSE_CHUNK,),
        dtype=np.int32,
        compression="lzf",
    )
    x.create_dataset(
        "indptr",
        shape=(n_obs + 1,),
        maxshape=(n_obs + 1,),
        chunks=(min(100_000, n_obs + 1),),
        dtype=np.int64,
        compression="lzf",
    )
    x["indptr"][0] = 0
    return out


def _expected_keep_count(src: h5py.File) -> int:
    """Count cells passing the threshold in the actual inner-joined matrix."""
    indptr = src["X/indptr"]
    n_obs = int(src["X"].attrs["shape"][0])
    keep_total = 0
    for start in range(0, n_obs, 5_000_000):
        end = min(start + 5_000_000, n_obs)
        keep_total += int((np.diff(indptr[start : end + 1]) >= MIN_GENES_INITIAL).sum())
    return keep_total


def _sample_categories(src: h5py.File) -> list[str]:
    n_obs = int(src["X"].attrs["shape"][0])
    del n_obs
    return [
        x.decode("utf-8") if isinstance(x, bytes) else str(x)
        for x in src["obs/sample/categories"][()]
    ]


def _append_sparse(data_ds, indices_ds, data, indices) -> int:
    old = int(data_ds.shape[0])
    new = old + int(len(data))
    if new > old:
        data_ds.resize((new,))
        indices_ds.resize((new,))
        data_ds[old:new] = data.astype(np.int32, copy=False)
        indices_ds[old:new] = indices.astype(np.int32, copy=False)
    return new


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    H5AD_DIR.mkdir(parents=True, exist_ok=True)

    with h5py.File(INPUT_H5AD, "r") as src:
        n_obs_before = int(src["X"].attrs["shape"][0])
        n_vars_before = int(src["X"].attrs["shape"][1])
        n_keep = _expected_keep_count(src)
        sample_categories = _sample_categories(src)
        before_counts = np.zeros(len(sample_categories), dtype=np.int64)
        after_counts = np.zeros(len(sample_categories), dtype=np.int64)
        out = _copy_root_schema(src, PART_H5AD, n_keep)
        out_data = out["X/data"]
        out_indices = out["X/indices"]
        out_indptr = out["X/indptr"]
        in_indptr = src["X/indptr"]
        in_data = src["X/data"]
        in_indices = src["X/indices"]
        sample_codes = src["obs/sample/codes"]

        out_row = 0
        out_nnz = 0
        progress_path = TABLE_DIR / "01_initial_cell_filter_progress.tsv"
        with progress_path.open("w", encoding="utf-8") as progress:
            progress.write("time\tinput_start\tinput_end\toutput_rows_total\tnnz_total\n")

        for start in range(0, n_obs_before, SCAN_ROWS):
            end = min(start + SCAN_ROWS, n_obs_before)
            ptr = in_indptr[start : end + 1]
            n_genes = np.diff(ptr)
            keep = n_genes >= MIN_GENES_INITIAL
            codes = sample_codes[start:end]
            before_counts += np.bincount(codes[codes >= 0], minlength=len(sample_categories))
            after_counts += np.bincount(codes[keep], minlength=len(sample_categories))
            if not keep.any():
                continue

            local_rows = np.nonzero(keep)[0].astype(np.int64)
            source_rows = local_rows + start
            counts = n_genes[local_rows].astype(np.int64, copy=False)
            data_start = int(ptr[0])
            data_end = int(ptr[-1])
            all_data = in_data[data_start:data_end]
            all_indices = in_indices[data_start:data_end]
            local_row_ids = np.repeat(
                np.arange(end - start, dtype=np.int32), n_genes.astype(np.int64)
            )
            keep_nnz = keep[local_row_ids]
            batch_data = all_data[keep_nnz]
            batch_indices = all_indices[keep_nnz]

            totals_all = np.bincount(
                local_row_ids[keep_nnz],
                weights=np.asarray(batch_data, dtype=np.float64),
                minlength=end - start,
            )
            totals = totals_all[local_rows]

            row_slice = slice(out_row, out_row + len(local_rows))
            for key, node in src["obs"].items():
                if isinstance(node, h5py.Group):
                    out["obs"][key]["codes"][row_slice] = node["codes"][source_rows]
                else:
                    out["obs"][key][row_slice] = node[source_rows]
            out["obs/n_genes_by_counts"][row_slice] = counts.astype(np.int32)
            out["obs/total_counts"][row_slice] = totals

            out_nnz = _append_sparse(out_data, out_indices, batch_data, batch_indices)
            out_indptr[out_row + 1 : out_row + len(local_rows) + 1] = (
                np.cumsum(counts, dtype=np.int64) + out_indptr[out_row]
            )
            out_row += len(local_rows)

            with progress_path.open("a", encoding="utf-8") as progress:
                progress.write(f"{_now()}\t{start}\t{end}\t{out_row}\t{out_nnz}\n")

        out.flush()
        out.close()

    if out_row != n_keep:
        raise RuntimeError(
            f"Initial filter wrote {out_row} rows, expected {n_keep}"
        )
    os.rename(PART_H5AD, OUTPUT_H5AD)

    sample_report = pd.DataFrame(
        {
            "sample": sample_categories,
            "n_cells_before_initial_filter": before_counts,
            "n_cells_after_initial_filter": after_counts,
        }
    )
    sample_report["n_cells_removed_initial_filter"] = (
        sample_report["n_cells_before_initial_filter"]
        - sample_report["n_cells_after_initial_filter"]
    )

    params = pd.DataFrame(
        [
            {
                "step": "initial_cell_filter",
                "input_h5ad": str(INPUT_H5AD),
                "output_h5ad_or_object": str(OUTPUT_H5AD),
                "n_obs_before": int(n_obs_before),
                "n_obs_after": int(n_keep),
                "n_vars_before": n_vars_before,
                "n_vars_after": n_vars_before,
                "min_genes_initial": MIN_GENES_INITIAL,
                "backend_package": "h5py_streaming",
                "code_file": str(CODE_FILE),
                "random_seed": SEED,
                "notes": "Initial cell filtering is done on disk because the saved pre-QC merge has 117.8M rows.",
            }
        ]
    )
    params.to_csv(TABLE_DIR / "01_initial_cell_filter_parameters.csv", index=False)
    sample_report = sample_report[sample_report["sample"] != "NA"].copy()
    sample_report.to_csv(TABLE_DIR / "initial_cell_filter_by_sample.csv", index=False)

    with (TABLE_DIR / "readme.txt").open("a", encoding="utf-8") as fh:
        fh.write("\n03-qc initial cell filter completed.\n")
        fh.write(f"Input: {INPUT_H5AD}\n")
        fh.write(f"Output: {OUTPUT_H5AD}\n")
        fh.write(f"min_genes_initial={MIN_GENES_INITIAL}\n")

    versions = [
        f"python={platform.python_version()}",
        f"platform={platform.platform()}",
        f"numpy={np.__version__}",
        f"pandas={pd.__version__}",
        f"h5py={h5py.__version__}",
        "environment=/mnt/disk18t/lr_xcy/riku/crc_val/uv_envs/main/.venv",
        f"code_file={CODE_FILE}",
        f"random_seed={SEED}",
    ]
    (TABLE_DIR / "package_versions_initial_filter.txt").write_text(
        "\n".join(versions) + "\n", encoding="utf-8"
    )
    print(params.to_string(index=False))


if __name__ == "__main__":
    main()
