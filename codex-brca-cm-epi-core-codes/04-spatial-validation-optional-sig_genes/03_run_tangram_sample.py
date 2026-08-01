#!/usr/bin/env python3
"""Run canonical Tangram pseudobulk mapping for exactly one spatial sample."""

from __future__ import annotations

import argparse
import importlib.metadata
import importlib.util
import json
import random
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.optimize import nnls

SEED = 42
random.seed(SEED)
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/brca_cm-epi_codex")
WORKFLOW = ROOT / "epi-cm-core-workflow"
REFERENCE = WORKFLOW / "h5ad/04-spatial-validation-optional-sig_genes/01-input-audit-and-reference/subtype_pseudobulk_reference_top100_significant_positive_deg.h5ad"
SPATIAL_INPUT_DIR = WORKFLOW / "h5ad/04-spatial-validation-optional-sig_genes/01-input-audit-and-reference"
TABLE_DIR = WORKFLOW / "tables/04-spatial-validation-optional-sig_genes/02-tangram-mapping"
H5AD_DIR = WORKFLOW / "h5ad/04-spatial-validation-optional-sig_genes/02-tangram-mapping"
NMF_DIR = WORKFLOW / "tables/03-epi-cm-discovery/01-cm-lineage-analysis/02_balanced_joint_nmf"
H_DF_PATH = NMF_DIR / "H_df.csv"
MINMAX_PATH = NMF_DIR / "column_minmax_params.csv"
VALID_SAMPLES = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]
EPOCHS = 350
LEARNING_RATE = 0.05
DEVICE = "cuda:0"


def pkg(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True, choices=VALID_SAMPLES)
    return parser.parse_args()


def normalize_spatial_sample(adata: ad.AnnData, target_sum: float = 1e4) -> ad.AnnData:
    adata.var_names_make_unique()
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    else:
        adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    return adata


def abundance_frame(asp: ad.AnnData, asc: ad.AnnData) -> pd.DataFrame:
    value = asp.obsm["tangram_ct_pred"]
    if isinstance(value, pd.DataFrame):
        frame = value.copy()
    else:
        categories = asc.obs["cell_subtype"].astype("category").cat.categories.astype(str).tolist()
        frame = pd.DataFrame(np.asarray(value), index=asp.obs_names.astype(str), columns=categories)
    frame.index = frame.index.astype(str)
    if set(frame.index) != set(asp.obs_names.astype(str)):
        raise ValueError("Tangram abundance spot IDs do not match spatial h5ad")
    return frame.reindex(asp.obs_names.astype(str)).astype(float)


def compute_cm_activity(abundance: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Refit canonical Module03 H activities from Tangram non-Epi composition."""
    h_df = pd.read_csv(H_DF_PATH, index_col=0)
    params = pd.read_csv(MINMAX_PATH, index_col=0)
    non_epi = h_df.columns.astype(str).tolist()
    missing = sorted(set(non_epi).difference(abundance.columns))
    if missing:
        raise ValueError(f"Tangram abundance missing Module03 non-Epi subtype columns: {missing}")
    non_epi_abundance = abundance.loc[:, non_epi].clip(lower=0.0)
    non_epi_total = non_epi_abundance.sum(axis=1)
    non_epi_composition = non_epi_abundance.div(non_epi_total.replace(0, np.nan), axis=0).fillna(0.0)
    params = params.reindex(non_epi)
    if params[["min", "range"]].isna().any().any():
        raise ValueError("Module03 column-minmax parameters do not cover all non-Epi subtypes")
    ranges = params["range"].replace(0, np.nan)
    v = non_epi_composition.sub(params["min"], axis=1).div(ranges, axis=1).fillna(0.0).clip(lower=0.0)
    design = h_df.loc[:, non_epi].to_numpy(dtype=float).T
    activities = np.zeros((len(v), h_df.shape[0]), dtype=float)
    for row_index, values in enumerate(v.to_numpy(dtype=float)):
        activities[row_index, :] = nnls(design, values, maxiter=design.shape[1] * 10)[0]
    cm = pd.DataFrame(activities, index=abundance.index, columns=h_df.index.astype(str))
    return cm, non_epi_composition, v


def main() -> None:
    started = time.time()
    sample = parse_args().sample
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    H5AD_DIR.mkdir(parents=True, exist_ok=True)
    score_path = TABLE_DIR / f"{sample}_tangram_pseudobulk_epi_cm_spot_scores.csv"
    mapper_path = H5AD_DIR / f"{sample}_tangram_mapper.h5ad"
    completion_path = TABLE_DIR / f"{sample}_mapping_completion.json"
    for path in [score_path, mapper_path, completion_path]:
        if path.exists():
            raise FileExistsError(f"Refusing to overwrite existing output: {path}")

    # FIXED Tangram preflight gate: no score/NNLS/marker projection substitute.
    if importlib.util.find_spec("tangram") is None:
        raise RuntimeError("Tangram is required for canonical spatial validation")
    if importlib.util.find_spec("torch") is None:
        raise RuntimeError("PyTorch is required for Tangram CUDA mapping")
    import torch
    import tangram as tg

    torch.manual_seed(SEED)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(SEED)
    if not torch.cuda.is_available():
        raise RuntimeError("Tangram canonical GPU route requires a usable cuda:0 device")
    torch.cuda.set_device(0)

    asc = ad.read_h5ad(REFERENCE)
    asp_path = SPATIAL_INPUT_DIR / f"{sample}_spatial_raw_counts.h5ad"
    asp = normalize_spatial_sample(ad.read_h5ad(asp_path))
    if "cell_subtype" not in asc.obs.columns:
        raise ValueError("Reference missing cell_subtype")
    common = [gene for gene in asc.var_names.astype(str) if gene in asp.var_names]
    if not common:
        raise ValueError(f"{sample}: no marker/intersected Tangram genes")

    # FIXED canonical Tangram pseudobulk call and parameters.
    tg.pp_adatas(asc, asp, genes=common)
    mapper = tg.map_cells_to_space(
        asc,
        asp,
        mode="cells",
        device=DEVICE,
        num_epochs=EPOCHS,
        learning_rate=LEARNING_RATE,
        random_state=SEED,
        verbose=False,
    )
    tg.project_cell_annotations(mapper, asp, annotation="cell_subtype")
    abundance = abundance_frame(asp, asc)
    if (abundance.to_numpy() < -1e-8).any():
        raise ValueError("Tangram subtype abundance contains negative values")
    abundance = abundance.clip(lower=0.0)

    epi_subtypes = [column for column in abundance.columns if str(column).startswith("Epi_")]
    if len(epi_subtypes) != 10:
        raise ValueError(f"Expected 10 epithelial subtype columns, found {len(epi_subtypes)}")
    epi_abundance = abundance.loc[:, epi_subtypes]
    epi_total = epi_abundance.sum(axis=1)
    epi_frac = epi_abundance.div(epi_total.replace(0, np.nan), axis=0).fillna(0.0)
    positive_epi = epi_total.gt(0)
    if positive_epi.any() and not np.allclose(epi_frac.loc[positive_epi].sum(axis=1), 1.0, rtol=1e-6, atol=1e-8):
        raise ValueError("EPIfrac row-normalization failed")

    cm_df, non_epi_composition, minmax_v = compute_cm_activity(abundance)
    obs_columns = ["barcode", "array_row", "array_col", "spatial_x", "spatial_y", "source_in_tissue", "Classification", "subtype", "patientid"]
    missing_obs = [column for column in obs_columns if column not in asp.obs.columns]
    if missing_obs:
        raise KeyError(f"{sample}: spatial obs missing {missing_obs}")
    score = asp.obs.loc[:, obs_columns].copy()
    score.insert(0, "sample", sample)
    score.insert(1, "spot_id", score.index.astype(str))
    score = score.join(cm_df.add_prefix("CMact__"), how="left")
    score = score.join(epi_frac.add_prefix("EPIfrac__"), how="left")
    if score.filter(like="CMact__").isna().any().any() or score.filter(like="EPIfrac__").isna().any().any():
        raise ValueError(f"{sample}: missing values after spot-score joins")

    # Mandatory canonical handoff: write CSV before mapper/audit h5ad or statistics.
    score.to_csv(score_path, index=False)
    mapper.uns["sample"] = sample
    mapper.uns["spot_score_csv"] = str(score_path.resolve())
    mapper.uns["tangram_parameters"] = {"mode": "cells", "device": DEVICE, "num_epochs": EPOCHS, "learning_rate": LEARNING_RATE, "random_state": SEED}
    mapper.write_h5ad(mapper_path, compression="gzip")

    pd.DataFrame({
        "subtype": abundance.columns,
        "min": abundance.min(axis=0).to_numpy(),
        "max": abundance.max(axis=0).to_numpy(),
        "mean": abundance.mean(axis=0).to_numpy(),
        "sum": abundance.sum(axis=0).to_numpy(),
    }).to_csv(TABLE_DIR / f"{sample}_tangram_subtype_abundance_summary.csv", index=False)
    pd.DataFrame([{
        "sample": sample,
        "code_file": str(Path(__file__).resolve()),
        "reference_h5ad": str(REFERENCE.resolve()),
        "spatial_h5ad": str(asp_path.resolve()),
        "n_spots": asp.n_obs,
        "n_reference_pseudocells": asc.n_obs,
        "n_common_genes": len(common),
        "spatial_normalization": "normalize_total target_sum=1e4 then log1p from layers[counts]",
        "tangram_mode": "cells",
        "device": DEVICE,
        "num_epochs": EPOCHS,
        "learning_rate": LEARNING_RATE,
        "random_state": SEED,
        "epi_fraction_denominator": "sum of all 10 projected epithelial subtype abundances within spot",
        "cm_activity": "NNLS refit of Module03 H on Tangram non-Epi composition after Module03 saved column-minmax transform",
        "h_df": str(H_DF_PATH.resolve()),
        "column_minmax_params": str(MINMAX_PATH.resolve()),
    }]).to_csv(TABLE_DIR / f"{sample}_mapping_parameters.csv", index=False)
    (TABLE_DIR / f"{sample}_package_versions.txt").write_text(
        f"python={sys.version.split()[0]}\nanndata={pkg('anndata')}\nscanpy={pkg('scanpy')}\ntangram-sc={pkg('tangram-sc')}\ntorch={pkg('torch')}\nnumpy={pkg('numpy')}\npandas={pkg('pandas')}\nscipy={pkg('scipy')}\ntorch_cuda={torch.version.cuda}\ngpu={torch.cuda.get_device_name(0)}\ncode={Path(__file__).resolve()}\nseed={SEED}\n",
        encoding="utf-8",
    )
    completion = {
        "status": "completed",
        "sample": sample,
        "n_spots": asp.n_obs,
        "n_common_genes": len(common),
        "n_subtypes": abundance.shape[1],
        "n_epi_subtypes": len(epi_subtypes),
        "n_cms": cm_df.shape[1],
        "n_positive_epi_spots": int(positive_epi.sum()),
        "epi_fraction_positive_rows_sum_one": True,
        "spot_score_written_before_mapper": True,
        "tangram_available": True,
        "torch_cuda_available": True,
        "device": DEVICE,
        "elapsed_seconds": time.time() - started,
        "score_csv": str(score_path.resolve()),
        "mapper_h5ad": str(mapper_path.resolve()),
        "seed": SEED,
    }
    completion_path.write_text(json.dumps(completion, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(completion, indent=2), flush=True)
    del mapper, asp, asc
    torch.cuda.empty_cache()


if __name__ == "__main__":
    main()
