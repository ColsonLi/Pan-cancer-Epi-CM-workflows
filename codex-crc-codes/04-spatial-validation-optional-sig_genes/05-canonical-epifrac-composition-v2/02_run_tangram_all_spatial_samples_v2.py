#!/usr/bin/env python3
"""Run canonical Tangram for all 37 Xenium regions with resumable CSV handoff."""

from __future__ import annotations

import argparse
import gc
import json
import platform
import random
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import tangram as tg
import torch

SEED = 42
random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)
torch.cuda.manual_seed_all(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/crc_codex")
WF = ROOT / "epi-cm-core-workflow"
TASK = "05-canonical-epifrac-composition-v2"
BLOCK = "04-spatial-validation-optional-sig_genes"
BASE_TABLE = WF / f"tables/{BLOCK}/{TASK}"
TABLE_DIR = BASE_TABLE / "02-tangram-mapping"
SCORE_DIR = TABLE_DIR / "spot_scores"
SPATIAL = WF / "h5ad/04-spatial-validation-optional/01-input-audit/adata_xenium_unprocessed.h5ad"
REFERENCE = WF / f"h5ad/{BLOCK}/{TASK}/01-reference-and-manifest/adata_cell_subtype_pseudobulk_deg_panel.h5ad"
COMMON_GENES = BASE_TABLE / "01-reference-and-manifest/tangram_common_deg_genes.csv"
SAMPLE_MANIFEST = BASE_TABLE / "01-reference-and-manifest/spatial_sample_scope_manifest.csv"
CM_LOADING = BASE_TABLE / "01-reference-and-manifest/cm_loading_fraction_non_epi_subtype_by_CM.csv"
EPI_SUBTYPES = BASE_TABLE / "01-reference-and-manifest/epithelial_subtypes.csv"
RUN_MANIFEST = TABLE_DIR / "tangram_sample_run_manifest.csv"

MODE = "cells"
DEVICE = "cuda:0"
NUM_EPOCHS = 350
LEARNING_RATE = 0.05
TARGET_SUM = 1e4


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--only-sample")
    return p.parse_args()


def atomic_csv(df: pd.DataFrame, path: Path) -> None:
    temp = path.with_name(path.name + ".tmp")
    df.to_csv(temp, index=False)
    temp.replace(path)


def save_manifest(df: pd.DataFrame) -> None:
    atomic_csv(df, RUN_MANIFEST)


def valid_existing_score(path: Path, n_expected: int, cm_names: list[str], epi_names: list[str]) -> bool:
    if not path.exists():
        return False
    head = pd.read_csv(path, nrows=2)
    required = {
        "sample", "spot_id", "barcode", "array_row", "array_col", "spatial_x", "spatial_y",
        *[f"CMact__{x}" for x in cm_names],
        *[f"EPIfrac__{x}" for x in epi_names],
    }
    if not required.issubset(head.columns):
        return False
    with path.open("rb") as fh:
        n_rows = sum(1 for _ in fh) - 1
    return n_rows == n_expected


def decode(values) -> list[str]:
    return [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in values]


def load_clean_spatial_counts():
    """Read only valid raw counts/coordinates/IDs; skip the broken unused scvi layer."""
    with h5py.File(SPATIAL, "r") as f:
        var_names = decode(f["var"]["_index"][()])
        cell_ids = np.asarray(decode(f["obs"]["cell_id"][()]), dtype=object)
        name_node = f["obs"]["name"]
        sample_categories = decode(name_node["categories"][()])
        sample_codes = name_node["codes"][()].astype(np.int16)
        coordinates = f["obsm"]["spatial"][()]
        x = f["layers"]["counts"]
        counts = sparse.csr_matrix(
            (x["data"][()], x["indices"][()], x["indptr"][()]),
            shape=tuple(int(v) for v in x.attrs["shape"]),
        )
    return counts, var_names, cell_ids, sample_categories, sample_codes, coordinates


def prepare_spatial_sample(
    counts,
    var_names,
    cell_ids,
    sample_categories,
    sample_codes,
    coordinates,
    sample: str,
    common: list[str],
) -> ad.AnnData:
    sample_i = sample_categories.index(sample)
    rows = np.flatnonzero(sample_codes == sample_i)
    barcode = cell_ids[rows].astype(str)
    obs = pd.DataFrame(index=pd.Index([f"{x}_{sample}" for x in barcode], dtype=str))
    obs["barcode"] = barcode
    obs["sample"] = sample
    asp = ad.AnnData(
        X=counts[rows, :].copy(),
        obs=obs,
        var=pd.DataFrame(index=pd.Index(var_names, dtype=str)),
    )
    asp.obsm["spatial"] = coordinates[rows, :].copy()
    if not asp.obs_names.is_unique:
        raise ValueError(f"{sample}: non-unique standardized obs_names")
    asp.layers["counts"] = asp.X.copy()
    sc.pp.normalize_total(asp, target_sum=TARGET_SUM)
    sc.pp.log1p(asp)
    asp = asp[:, common].copy()
    return asp


def map_one(
    spatial_parts,
    reference_template: ad.AnnData,
    sample_row: pd.Series,
    common: list[str],
    loadings: pd.DataFrame,
    epi_names: list[str],
) -> tuple[Path, dict[str, object]]:
    sample = str(sample_row["sample"])
    started = time.time()
    asp = prepare_spatial_sample(*spatial_parts, sample, common)
    asc = reference_template.copy()
    tg.pp_adatas(asc, asp, genes=common)
    mapper = tg.map_cells_to_space(
        asc,
        asp,
        mode=MODE,
        device=DEVICE,
        num_epochs=NUM_EPOCHS,
        learning_rate=LEARNING_RATE,
        random_state=SEED,
        verbose=False,
    )
    tg.project_cell_annotations(mapper, asp, annotation="cell_subtype")
    abundance = asp.obsm["tangram_ct_pred"].copy()
    if not isinstance(abundance, pd.DataFrame):
        abundance = pd.DataFrame(abundance, index=asp.obs_names)
    missing = sorted(set(loadings.index) - set(abundance.columns))
    missing_epi = sorted(set(epi_names) - set(abundance.columns))
    if missing or missing_epi:
        raise ValueError(f"{sample}: missing Tangram subtype columns: CM={missing}, Epi={missing_epi}")

    epi_abundance = abundance.loc[:, epi_names].astype(float)
    epi_total = epi_abundance.sum(axis=1)
    epi = epi_abundance.div(epi_total.replace(0, np.nan), axis=0).fillna(0.0)
    positive_epi = epi_total.gt(0)
    if positive_epi.any() and not np.allclose(
        epi.loc[positive_epi].sum(axis=1).to_numpy(), 1.0, rtol=1e-6, atol=1e-8
    ):
        raise ValueError(f"{sample}: EPIfrac row-normalization failed")
    cm = abundance.loc[:, loadings.index].to_numpy(dtype=float) @ loadings.to_numpy(dtype=float)
    cm = pd.DataFrame(cm, index=abundance.index, columns=loadings.columns)
    coords = np.asarray(asp.obsm["spatial"], dtype=float)
    score = pd.DataFrame(index=asp.obs_names)
    score["sample"] = sample
    score["spot_id"] = score.index.astype(str)
    score["barcode"] = asp.obs["barcode"].astype(str).to_numpy()
    score["patient_id"] = str(sample_row["patient_id"])
    score["tissue_region"] = str(sample_row["tissue_region"])
    score["status"] = str(sample_row["status"])
    score["array_row"] = coords[:, 1]
    score["array_col"] = coords[:, 0]
    score["spatial_x"] = coords[:, 0]
    score["spatial_y"] = coords[:, 1]
    for name in cm.columns:
        score[f"CMact__{name}"] = cm[name].to_numpy()
    for name in epi.columns:
        score[f"EPIfrac__{name}"] = epi[name].to_numpy()
    out = SCORE_DIR / f"{sample}_tangram_pseudobulk_epi_cm_spot_scores.csv"
    atomic_csv(score.reset_index(drop=True), out)
    metrics = {
        "sample": sample,
        "n_spatial_cells": int(asp.n_obs),
        "n_common_genes": len(common),
        "n_reference_subtypes": int(asc.n_obs),
        "n_cm": len(loadings.columns),
        "n_epi_subtypes": len(epi_names),
        "n_positive_epi_rows": int(positive_epi.sum()),
        "n_zero_epi_rows": int((~positive_epi).sum()),
        "epifrac_definition": "tangram_ct_pred epithelial columns row-normalized within spatial observation",
        "mode": MODE,
        "device": DEVICE,
        "num_epochs": NUM_EPOCHS,
        "learning_rate": LEARNING_RATE,
        "random_state": SEED,
        "score_csv": str(out),
        "elapsed_seconds": round(time.time() - started, 3),
    }
    del mapper, abundance, cm, epi, asc, asp, score
    gc.collect()
    torch.cuda.empty_cache()
    return out, metrics


def main() -> None:
    cli = args()
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA unavailable")
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    SCORE_DIR.mkdir(parents=True, exist_ok=True)
    sample_manifest = pd.read_csv(SAMPLE_MANIFEST)
    if cli.only_sample:
        if cli.only_sample not in set(sample_manifest["sample"]):
            raise KeyError(cli.only_sample)
        run_samples = sample_manifest[sample_manifest["sample"].eq(cli.only_sample)].copy()
    else:
        run_samples = sample_manifest.copy()
    common = pd.read_csv(COMMON_GENES)["gene"].astype(str).tolist()
    loadings = pd.read_csv(CM_LOADING).set_index("cell_subtype")
    epi_names = pd.read_csv(EPI_SUBTYPES)["cell_subtype"].astype(str).tolist()
    reference = ad.read_h5ad(REFERENCE)
    reference.obs["cell_subtype"] = reference.obs_names.astype(str)

    run_manifest = sample_manifest.copy()
    run_manifest["status_run"] = "pending"
    run_manifest["score_csv"] = ""
    run_manifest["elapsed_seconds"] = np.nan
    run_manifest["error"] = ""
    for i, row in run_manifest.iterrows():
        out = SCORE_DIR / f"{row['sample']}_tangram_pseudobulk_epi_cm_spot_scores.csv"
        if valid_existing_score(out, int(row["n_spatial_cells"]), list(loadings.columns), epi_names):
            run_manifest.loc[i, ["status_run", "score_csv"]] = ["completed", str(out)]
    save_manifest(run_manifest)

    spatial_parts = load_clean_spatial_counts()
    for _, sample_row in run_samples.iterrows():
        sample = str(sample_row["sample"])
        i = run_manifest.index[run_manifest["sample"].eq(sample)][0]
        if run_manifest.loc[i, "status_run"] == "completed":
            print(f"SKIP completed {sample}", flush=True)
            continue
        run_manifest.loc[i, "status_run"] = "running"
        save_manifest(run_manifest)
        print(f"START {sample} n={int(sample_row['n_spatial_cells'])}", flush=True)
        try:
            out, metrics = map_one(spatial_parts, reference, sample_row, common, loadings, epi_names)
            run_manifest.loc[i, "status_run"] = "completed"
            run_manifest.loc[i, "score_csv"] = str(out)
            run_manifest.loc[i, "elapsed_seconds"] = metrics["elapsed_seconds"]
            (TABLE_DIR / f"{sample}_mapping_metrics.json").write_text(
                json.dumps(metrics, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
            )
            print(f"DONE {sample} elapsed={metrics['elapsed_seconds']}", flush=True)
        except Exception as exc:
            run_manifest.loc[i, "status_run"] = "failed"
            run_manifest.loc[i, "error"] = repr(exc)
            save_manifest(run_manifest)
            raise
        save_manifest(run_manifest)

    cm_names = list(loadings.columns)
    pair_manifest = pd.MultiIndex.from_product(
        [sample_manifest["sample"].astype(str), cm_names, epi_names],
        names=["sample", "CM", "epi_subtype"],
    ).to_frame(index=False)
    pair_manifest.to_csv(TABLE_DIR / "all_sample_cm_epi_pair_manifest.csv", index=False)
    summary = {
        "n_expected_samples": 37,
        "n_completed_samples": int(run_manifest["status_run"].eq("completed").sum()),
        "n_expected_pairs": int(len(pair_manifest)),
        "n_cm": len(cm_names),
        "n_epi_subtypes": len(epi_names),
        "canonical_parameters": {
            "mode": MODE,
            "device": DEVICE,
            "num_epochs": NUM_EPOCHS,
            "learning_rate": LEARNING_RATE,
            "random_state": SEED,
        },
        "python": platform.python_version(),
        "torch": torch.__version__,
    }
    (TABLE_DIR / "tangram_mapping_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
