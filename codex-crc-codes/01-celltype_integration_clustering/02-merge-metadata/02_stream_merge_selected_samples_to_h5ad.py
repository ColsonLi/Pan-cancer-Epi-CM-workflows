#!/usr/bin/env python3
"""Stream selected samples into a merged raw-count h5ad.

The selected CRC inputs include very large raw droplet matrices. Loading them
through AnnData before the initial QC filter is not practical, so this script
uses the h5ad on-disk schema directly while preserving the requested pre-QC
object: selected rows are written first, QC is intentionally not applied here.
"""

from __future__ import annotations

import json
import os
import platform
import random
import shutil
import subprocess
import sys
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
TABLE_DIR = WORKFLOW_ROOT / "tables/01-celltype_integration_clustering/02-merge-metadata"
H5AD_DIR = WORKFLOW_ROOT / "h5ad/01-celltype_integration_clustering/02-merge-metadata"
SELECTED_SAMPLES = (
    WORKFLOW_ROOT
    / "tables/01-celltype_integration_clustering/01-select-samples/02_eligible_samples_for_merge.csv"
)
FEASIBILITY_BY_H5AD = TABLE_DIR / "merge_feasibility_by_h5ad.csv"
COMMON_GENES = TABLE_DIR / "inner_join_common_genes.csv"
OUT_H5AD = H5AD_DIR / "adata_merge.h5ad"
PART_H5AD = H5AD_DIR / "adata_merge.h5ad.part"
POST_MERGE_GENE_ID_CONVERTER = (
    ANALYSIS_ROOT
    / "h5ad_gene_id_conversion_after_merge/convert_all_h5ad_gene_ids.py"
)
CODE_FILE = Path(__file__)

CHUNK_NNZ_TARGET = 12_000_000
MIN_ROWS_PER_CHUNK = 20_000
MAX_ROWS_PER_CHUNK = 2_000_000
STRING_CHUNK = 100_000
SPARSE_CHUNK = 8_000_000
TEXT_NA = "NA"


def _now() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def _decode_array(values) -> list[str]:
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values]


def _categorical_codes(obs_group, column: str) -> tuple[list[str], np.ndarray]:
    node = obs_group[column]
    if isinstance(node, h5py.Group) and "categories" in node and "codes" in node:
        return _decode_array(node["categories"][()]), node["codes"][()]
    values = np.asarray(_decode_array(node[()]))
    categories, codes = np.unique(values, return_inverse=True)
    return categories.tolist(), codes.astype(np.int64)


def _var_names(var_group) -> list[str]:
    if "_index" in var_group:
        return _decode_array(var_group["_index"][()])
    if "symbol" in var_group:
        return _decode_array(var_group["symbol"][()])
    raise KeyError("Could not find var/_index or var/symbol")


def _string_ds(group, name: str, shape: int, chunks: int = STRING_CHUNK):
    return group.create_dataset(
        name,
        shape=(shape,),
        maxshape=(shape,),
        chunks=(min(chunks, max(shape, 1)),),
        dtype=h5py.string_dtype(encoding="utf-8"),
    )


def _array_ds(group, name: str, shape: int, dtype, chunks: int = STRING_CHUNK):
    return group.create_dataset(
        name,
        shape=(shape,),
        maxshape=(shape,),
        chunks=(min(chunks, max(shape, 1)),),
        dtype=dtype,
        compression="lzf",
    )


def _code_dtype(n_categories: int):
    if n_categories <= np.iinfo(np.int8).max:
        return np.int8
    if n_categories <= np.iinfo(np.int16).max:
        return np.int16
    return np.int32


def _clean_value(value) -> str:
    if pd.isna(value):
        return TEXT_NA
    text = str(value)
    if text == "" or text.lower() == "nan":
        return TEXT_NA
    return text


def _create_empty_dict_group(root, name: str):
    group = root.create_group(name)
    group.attrs["encoding-type"] = "dict"
    group.attrs["encoding-version"] = "0.1.0"
    return group


def _init_h5ad(
    out_path: Path,
    n_obs: int,
    genes: list[str],
    selected: pd.DataFrame,
) -> tuple[h5py.File, dict[str, h5py.Dataset], dict[str, dict[str, int]], h5py.Dataset, h5py.Dataset, h5py.Dataset]:
    H5AD_DIR.mkdir(parents=True, exist_ok=True)
    if out_path.exists() or PART_H5AD.exists():
        raise FileExistsError(
            f"Output already exists: {out_path} or {PART_H5AD}. "
            "The skill forbids overwriting existing outputs without an explicit cleanup request."
        )

    f = h5py.File(PART_H5AD, "w")
    f.attrs["encoding-type"] = "anndata"
    f.attrs["encoding-version"] = "0.1.0"

    obs = f.create_group("obs")
    obs.attrs["encoding-type"] = "dataframe"
    obs.attrs["encoding-version"] = "0.2.0"
    obs.attrs["_index"] = "_index"

    base_cols = ["sample", "series", "status", "original_barcode", "source_h5ad", "source_row_index"]
    metadata_cols = [c for c in selected.columns if c not in {"h5ad"}]
    obs_columns = []
    for col in base_cols + metadata_cols:
        if col not in obs_columns:
            obs_columns.append(col)
    obs.attrs["column-order"] = np.asarray(obs_columns, dtype=object)

    obs_datasets: dict[str, h5py.Dataset] = {}
    obs_datasets["_index"] = _string_ds(obs, "_index", n_obs)
    obs_datasets["original_barcode"] = _string_ds(obs, "original_barcode", n_obs)
    obs_datasets["source_row_index"] = _array_ds(obs, "source_row_index", n_obs, np.int64)

    categorical_maps: dict[str, dict[str, int]] = {}
    categorical_columns = [c for c in obs_columns if c not in {"original_barcode", "source_row_index"}]
    selected_for_codes = selected.copy()
    selected_for_codes["sample"] = selected_for_codes["sample_id"].map(_clean_value)
    selected_for_codes["series"] = selected_for_codes["dataset"].map(_clean_value)
    selected_for_codes["source_h5ad"] = selected_for_codes["h5ad"].map(_clean_value)
    for col in categorical_columns:
        group = obs.create_group(col)
        group.attrs["encoding-type"] = "categorical"
        group.attrs["encoding-version"] = "0.2.0"
        group.attrs["ordered"] = False
        values = sorted({_clean_value(v) for v in selected_for_codes[col].tolist()})
        if TEXT_NA not in values:
            values.append(TEXT_NA)
        group.create_dataset(
            "categories",
            data=np.asarray(values, dtype=object),
            dtype=h5py.string_dtype(encoding="utf-8"),
        )
        group["categories"].attrs["encoding-type"] = "string-array"
        group["categories"].attrs["encoding-version"] = "0.2.0"
        codes = _array_ds(group, "codes", n_obs, _code_dtype(len(values)))
        group["codes"].attrs["encoding-type"] = "array"
        group["codes"].attrs["encoding-version"] = "0.2.0"
        obs_datasets[col] = codes
        categorical_maps[col] = {value: i for i, value in enumerate(values)}

    obs["_index"].attrs["encoding-type"] = "string-array"
    obs["_index"].attrs["encoding-version"] = "0.2.0"
    obs["original_barcode"].attrs["encoding-type"] = "string-array"
    obs["original_barcode"].attrs["encoding-version"] = "0.2.0"
    obs["source_row_index"].attrs["encoding-type"] = "array"
    obs["source_row_index"].attrs["encoding-version"] = "0.2.0"

    var = f.create_group("var")
    var.attrs["encoding-type"] = "dataframe"
    var.attrs["encoding-version"] = "0.2.0"
    var.attrs["_index"] = "_index"
    var.attrs["column-order"] = np.asarray([], dtype="S1")
    var.create_dataset(
        "_index",
        data=np.asarray(genes, dtype=object),
        dtype=h5py.string_dtype(encoding="utf-8"),
    )
    var["_index"].attrs["encoding-type"] = "string-array"
    var["_index"].attrs["encoding-version"] = "0.2.0"

    x = f.create_group("X")
    x.attrs["encoding-type"] = "csr_matrix"
    x.attrs["encoding-version"] = "0.1.0"
    x.attrs["shape"] = np.asarray([n_obs, len(genes)], dtype=np.int64)
    data = x.create_dataset(
        "data",
        shape=(0,),
        maxshape=(None,),
        chunks=(SPARSE_CHUNK,),
        dtype=np.int64,
        compression="lzf",
    )
    indices = x.create_dataset(
        "indices",
        shape=(0,),
        maxshape=(None,),
        chunks=(SPARSE_CHUNK,),
        dtype=np.int32,
        compression="lzf",
    )
    indptr = x.create_dataset(
        "indptr",
        shape=(n_obs + 1,),
        maxshape=(n_obs + 1,),
        chunks=(min(STRING_CHUNK, n_obs + 1),),
        dtype=np.int64,
        compression="lzf",
    )
    indptr[0] = 0

    for name in ["layers", "obsm", "obsp", "uns", "varm", "varp"]:
        _create_empty_dict_group(f, name)

    return f, obs_datasets, categorical_maps, data, indices, indptr


def _append_sparse(data_ds, indices_ds, data_chunk: np.ndarray, index_chunk: np.ndarray) -> int:
    old = int(data_ds.shape[0])
    new = old + int(data_chunk.shape[0])
    if new > old:
        data_ds.resize((new,))
        indices_ds.resize((new,))
        data_ds[old:new] = data_chunk.astype(np.int64, copy=False)
        indices_ds[old:new] = index_chunk.astype(np.int32, copy=False)
    return new


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    selected = pd.read_csv(SELECTED_SAMPLES)
    feasibility = pd.read_csv(FEASIBILITY_BY_H5AD)
    common_genes = pd.read_csv(COMMON_GENES)["gene"].astype(str).tolist()
    gene_to_out = {gene: i for i, gene in enumerate(common_genes)}
    n_obs_total = int(feasibility["n_selected_obs_unfiltered"].sum())

    f, obs_datasets, categorical_maps, data_ds, indices_ds, indptr_ds = _init_h5ad(
        OUT_H5AD, n_obs_total, common_genes, selected
    )

    progress_path = TABLE_DIR / "merge_progress.tsv"
    with progress_path.open("w", encoding="utf-8") as progress:
        progress.write(
            "time\th5ad\tdataset\tstart_row\tend_row\toutput_rows_total\tnnz_total\n"
        )

    selected_by_sample = selected.set_index("sample_id", drop=False)
    output_row = 0
    nnz_total = 0
    manifest_rows = []
    sample_manifest = []

    try:
        for _, file_info in feasibility.iterrows():
            h5ad_path = str(file_info["h5ad"])
            sub = selected[selected["h5ad"].astype(str) == h5ad_path].copy()
            selected_sample_ids = set(sub["sample_id"].astype(str))
            if not selected_sample_ids:
                continue
            with h5py.File(h5ad_path, "r") as src:
                source_genes = _var_names(src["var"])
                source_to_out = np.full(len(source_genes), -1, dtype=np.int32)
                for source_i, gene in enumerate(source_genes):
                    target_i = gene_to_out.get(gene)
                    if target_i is not None:
                        source_to_out[source_i] = target_i

                sample_categories, sample_codes = _categorical_codes(src["obs"], "sample_id")
                selected_source_codes = np.asarray(
                    [i for i, sid in enumerate(sample_categories) if sid in selected_sample_ids],
                    dtype=sample_codes.dtype,
                )
                if selected_source_codes.size == 0:
                    raise ValueError(f"No selected sample_id matched in {h5ad_path}")

                source_code_to_sample = {i: sid for i, sid in enumerate(sample_categories)}
                source_code_to_meta_code: dict[str, np.ndarray] = {}
                for col, cmap in categorical_maps.items():
                    arr = np.full(len(sample_categories), cmap[TEXT_NA], dtype=np.int32)
                    for source_code, sid in source_code_to_sample.items():
                        if sid in selected_by_sample.index:
                            row = selected_by_sample.loc[sid]
                            if isinstance(row, pd.DataFrame):
                                row = row.iloc[0]
                            if col == "sample":
                                value = sid
                            elif col == "series":
                                value = row["dataset"]
                            elif col == "source_h5ad":
                                value = row["h5ad"]
                            else:
                                value = row[col] if col in row.index else TEXT_NA
                            arr[source_code] = cmap.get(_clean_value(value), cmap[TEXT_NA])
                    source_code_to_meta_code[col] = arr

                mean_nnz = max(float(file_info["mean_nnz_per_obs_file"]), 1.0)
                rows_per_chunk = int(CHUNK_NNZ_TARGET / mean_nnz)
                rows_per_chunk = max(MIN_ROWS_PER_CHUNK, min(MAX_ROWS_PER_CHUNK, rows_per_chunk))
                n_source_obs = int(file_info["n_obs_file"])
                x = src["X"]
                indptr_src = x["indptr"]
                indices_src = x["indices"]
                data_src = x["data"]
                index_src = src["obs"]["_index"]
                file_output_start = output_row
                file_nnz_start = nnz_total

                for start in range(0, n_source_obs, rows_per_chunk):
                    end = min(start + rows_per_chunk, n_source_obs)
                    code_chunk = sample_codes[start:end]
                    selected_mask = np.isin(code_chunk, selected_source_codes)
                    n_selected = int(selected_mask.sum())
                    if n_selected == 0:
                        continue

                    selected_codes_chunk = code_chunk[selected_mask]
                    source_row_indices = np.nonzero(selected_mask)[0].astype(np.int64) + start
                    barcodes_all = index_src[start:end]
                    barcodes = np.asarray(barcodes_all, dtype=object)[selected_mask]
                    barcodes = np.asarray(
                        [
                            x.decode("utf-8") if isinstance(x, bytes) else str(x)
                            for x in barcodes
                        ],
                        dtype=object,
                    )
                    sample_ids = np.asarray(
                        [source_code_to_sample[int(c)] for c in selected_codes_chunk],
                        dtype=object,
                    )
                    merged_ids = np.asarray(
                        [f"{sample}|{barcode}" for sample, barcode in zip(sample_ids, barcodes)],
                        dtype=object,
                    )

                    out_slice = slice(output_row, output_row + n_selected)
                    obs_datasets["_index"][out_slice] = merged_ids
                    obs_datasets["original_barcode"][out_slice] = barcodes
                    obs_datasets["source_row_index"][out_slice] = source_row_indices
                    for col, mapper in source_code_to_meta_code.items():
                        obs_datasets[col][out_slice] = mapper[selected_codes_chunk].astype(
                            obs_datasets[col].dtype, copy=False
                        )

                    ptr = indptr_src[start : end + 1].astype(np.int64, copy=False)
                    row_nnz_all = np.diff(ptr)
                    local_row_ids_all = np.repeat(
                        np.arange(end - start, dtype=np.int32), row_nnz_all
                    )
                    data_start = int(ptr[0])
                    data_end = int(ptr[-1])
                    src_indices_all = indices_src[data_start:data_end]
                    src_data_all = data_src[data_start:data_end]
                    keep_selected_nnz = selected_mask[local_row_ids_all]
                    if keep_selected_nnz.any():
                        selected_local_rows = local_row_ids_all[keep_selected_nnz]
                        out_local_rows = np.cumsum(selected_mask, dtype=np.int64)[
                            selected_local_rows
                        ] - 1
                        selected_src_indices = src_indices_all[keep_selected_nnz]
                        out_indices_all = source_to_out[selected_src_indices]
                        keep_common = out_indices_all >= 0
                        out_local_rows = out_local_rows[keep_common]
                        out_indices = out_indices_all[keep_common]
                        out_data = src_data_all[keep_selected_nnz][keep_common]
                        row_counts = np.bincount(out_local_rows, minlength=n_selected)
                    else:
                        out_indices = np.asarray([], dtype=np.int32)
                        out_data = np.asarray([], dtype=np.int64)
                        row_counts = np.zeros(n_selected, dtype=np.int64)

                    nnz_total = _append_sparse(data_ds, indices_ds, out_data, out_indices)
                    indptr_ds[output_row + 1 : output_row + n_selected + 1] = (
                        np.cumsum(row_counts, dtype=np.int64) + indptr_ds[output_row]
                    )
                    output_row += n_selected

                    with progress_path.open("a", encoding="utf-8") as progress:
                        progress.write(
                            f"{_now()}\t{h5ad_path}\t{file_info['dataset']}\t{start}\t{end}\t{output_row}\t{nnz_total}\n"
                        )

                manifest_rows.append(
                    {
                        "h5ad": h5ad_path,
                        "dataset": file_info["dataset"],
                        "n_selected_samples": int(file_info["n_selected_samples"]),
                        "n_obs_written": int(output_row - file_output_start),
                        "nnz_written_after_inner_gene_join": int(nnz_total - file_nnz_start),
                        "n_common_genes": len(common_genes),
                        "join_strategy": "inner",
                        "source_rows_per_chunk": rows_per_chunk,
                    }
                )

                for sid in sorted(selected_sample_ids):
                    sample_manifest.append(
                        {
                            "sample": sid,
                            "series": str(sub[sub["sample_id"].astype(str) == sid]["dataset"].iloc[0]),
                            "source_h5ad": h5ad_path,
                            "included_in_merge": True,
                        }
                    )

        if output_row != n_obs_total:
            raise RuntimeError(f"Wrote {output_row} rows, expected {n_obs_total}")

        f.flush()
    except Exception:
        f.close()
        raise
    else:
        f.close()
        os.rename(PART_H5AD, OUT_H5AD)

    # The merge is complete and safely finalized before gene identifiers are
    # converted. This is a separate post-merge metadata step; it does not rerun
    # the merge or rewrite the expression matrix.
    subprocess.run(
        [sys.executable, str(POST_MERGE_GENE_ID_CONVERTER), "--only", str(OUT_H5AD)],
        check=True,
    )

    pd.DataFrame(manifest_rows).to_csv(TABLE_DIR / "merge_manifest.csv", index=False)
    pd.DataFrame(sample_manifest).drop_duplicates().to_csv(
        TABLE_DIR / "merged_sample_manifest.csv", index=False
    )
    summary = pd.DataFrame(
        [
            {
                "step": "stream_merge_selected_samples",
                "input_selected_sample_csv": str(SELECTED_SAMPLES),
                "output_h5ad": str(OUT_H5AD),
                "n_obs_after": int(output_row),
                "n_vars_after": len(common_genes),
                "nnz_after_inner_gene_join": int(nnz_total),
                "join_strategy": "inner",
                "qc_applied": False,
                "raw_counts_preserved_in_X": True,
                "layers_copied_from_source": False,
                "raw_copied_from_source": False,
                "obsm_copied_from_source": False,
                "required_obs_columns": "sample;series;status;original_barcode",
                "post_merge_gene_id_conversion": True,
                "post_merge_gene_id_converter": str(POST_MERGE_GENE_ID_CONVERTER),
                "code_file": str(CODE_FILE),
                "random_seed": SEED,
            }
        ]
    )
    summary.to_csv(TABLE_DIR / "merge_summary.csv", index=False)

    readme = TABLE_DIR / "readme.txt"
    with readme.open("a", encoding="utf-8") as fh:
        fh.write("\nRaw selected-sample merge completed.\n")
        fh.write(f"Output h5ad: {OUT_H5AD}\n")
        fh.write("The object is pre-QC; no initial cell/gene/QC filtering was applied.\n")
        fh.write("Source layers/raw/obsm were intentionally not copied.\n")
        fh.write(f"After merge finalization, gene IDs were converted by: {POST_MERGE_GENE_ID_CONVERTER}\n")
        fh.write("The post-merge conversion preserves Ensembl IDs in var['ensembl_id'] and does not rewrite X.\n")

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
    (TABLE_DIR / "package_versions_merge.txt").write_text(
        "\n".join(versions) + "\n", encoding="utf-8"
    )

    print(json.dumps(summary.iloc[0].to_dict(), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
