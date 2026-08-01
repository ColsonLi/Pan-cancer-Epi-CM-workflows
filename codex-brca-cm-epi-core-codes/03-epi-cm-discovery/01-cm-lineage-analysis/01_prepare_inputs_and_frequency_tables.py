#!/usr/bin/env python3
"""Prepare BRCA tumor-only subtype frequency matrices and status preflight."""

from __future__ import annotations

import importlib.metadata
import json
import re
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from cm_analysis_common import PREP_DIR, ROOT, SEED, WORKFLOW, write_json


INPUT_H5AD = (
    WORKFLOW / "h5ad" / "02-cell_subtype_integration_clustering"
    / "06-project-subtypes-to-full-adata" / "adata_anno_cellsubtype.h5ad"
)
CODE_PATH = Path(__file__).resolve()
MIN_NON_EPI_CELLS = 50
MIN_EPI_CELLS = 1
EPI_CELL_TYPE = "Epithelial Cells"

PREFIX_TO_CELL_TYPE = {
    "Epi": "Epithelial Cells", "T": "T Cells", "Mye": "Myeloid Cells",
    "B": "B Cells", "Plasma": "Plasma Cells", "Endo": "Endothelial Cells",
    "Stromal": "Stromal Cells", "PVL": "Perivascular Cells",
}


def status_token(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(value).strip().lower()).strip()


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def main() -> None:
    started = time.time()
    completion_path = PREP_DIR / "input_preparation_completion.json"
    if completion_path.exists():
        completion = json.loads(completion_path.read_text())
        if completion.get("status") == "completed":
            print(json.dumps({"status": "valid_existing_input_preparation_reused"}, indent=2))
            return
        raise FileExistsError("Existing input-preparation completion is invalid.")
    PREP_DIR.mkdir(parents=True, exist_ok=True)
    atlas = ad.read_h5ad(INPUT_H5AD, backed="r")
    try:
        required = {"sample", "series", "status", "cell_type", "cell_subtype"}
        if not required.issubset(atlas.obs.columns):
            raise ValueError(f"Projected atlas lacks obs columns: {sorted(required - set(atlas.obs))}")
        obs = atlas.obs[list(required)].copy()
    finally:
        atlas.file.close()
    obs = obs.dropna().copy()
    obs["sample"] = obs["sample"].astype(str)
    obs["series"] = obs["series"].astype(str)
    obs["raw_status"] = obs["status"].astype(str)
    obs["cell_type"] = obs["cell_type"].astype(str)
    obs["cell_subtype"] = obs["cell_subtype"].astype(str)

    obs["subtype_prefix"] = obs["cell_subtype"].str.split("_", n=1).str[0]
    obs["prefix_derived_cell_type"] = obs["subtype_prefix"].map(PREFIX_TO_CELL_TYPE)
    conflicts = obs.loc[
        obs["prefix_derived_cell_type"].isna()
        | obs["prefix_derived_cell_type"].ne(obs["cell_type"]),
        ["cell_subtype", "subtype_prefix", "prefix_derived_cell_type", "cell_type"],
    ].drop_duplicates()
    conflicts.to_csv(PREP_DIR / "cell_subtype_prefix_cell_type_conflicts.csv", index=False)
    if len(conflicts):
        raise ValueError("Subtype-prefix-derived cell types conflict with projected cell_type.")

    raw_sample_meta = obs[["sample", "series", "raw_status"]].drop_duplicates()
    if raw_sample_meta["sample"].duplicated().any():
        raise ValueError("Conflicting series/status metadata within biological samples.")
    aliases = {
        "tumor": "tumor", "tumour": "tumor", "primary tumor": "tumor",
        "primary tumour": "tumor", "metastatic tumor": "tumor",
        "metastatic tumour": "tumor", "metastasis": "tumor",
        "normal": "normal-like", "normal like": "normal-like",
        "adjacent normal": "normal-like", "normal adjacent": "normal-like",
    }
    raw_sample_meta["status"] = raw_sample_meta["raw_status"].map(
        lambda value: aliases.get(status_token(value))
    )
    if raw_sample_meta["status"].isna().any():
        bad = raw_sample_meta.loc[raw_sample_meta["status"].isna(), "raw_status"].unique()
        raise ValueError(f"Unresolved statuses: {bad}")

    epi_mask = obs["cell_type"].eq(EPI_CELL_TYPE)
    compartment = pd.DataFrame(index=sorted(obs["sample"].unique()))
    compartment.index.name = "sample"
    compartment["epithelial_cells"] = obs.loc[epi_mask].groupby("sample").size()
    compartment["non_epi_cells"] = obs.loc[~epi_mask].groupby("sample").size()
    compartment = compartment.fillna(0).astype(int)
    compartment.to_csv(PREP_DIR / "sample_compartment_cell_counts.csv")

    inclusion = raw_sample_meta.set_index("sample").join(compartment)
    inclusion["keep_for_cm"] = inclusion["non_epi_cells"] >= MIN_NON_EPI_CELLS
    inclusion["keep_for_epi_cm"] = inclusion["keep_for_cm"] & (
        inclusion["epithelial_cells"] >= MIN_EPI_CELLS
    )
    inclusion["cm_exclusion_reason"] = np.where(
        inclusion["keep_for_cm"], "", f"non_epi_cells < {MIN_NON_EPI_CELLS}"
    )
    inclusion["epi_cm_exclusion_reason"] = np.where(
        inclusion["keep_for_epi_cm"], "",
        np.where(~inclusion["keep_for_cm"], inclusion["cm_exclusion_reason"],
                 f"epithelial_cells < {MIN_EPI_CELLS}"),
    )
    inclusion.reset_index().to_csv(PREP_DIR / "sample_inclusion_exclusion.csv", index=False)
    cm_samples = inclusion.index[inclusion["keep_for_cm"]].astype(str)
    if len(cm_samples) < 2:
        raise ValueError("Fewer than two samples pass CM eligibility.")

    non_epi_obs = obs.loc[~epi_mask & obs["sample"].isin(cm_samples)].copy()
    epi_obs = obs.loc[epi_mask & obs["sample"].isin(cm_samples)].copy()
    non_epi_counts = pd.crosstab(non_epi_obs["sample"], non_epi_obs["cell_subtype"]).astype(int)
    non_epi_counts = non_epi_counts.reindex(index=cm_samples, fill_value=0).sort_index(axis=1)
    non_epi_frequency = non_epi_counts.div(non_epi_counts.sum(axis=1), axis=0).fillna(0.0)
    epi_counts = pd.crosstab(epi_obs["sample"], epi_obs["cell_subtype"]).astype(int)
    epi_counts = epi_counts.reindex(index=cm_samples, fill_value=0).sort_index(axis=1)
    epi_frequency = epi_counts.div(epi_counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    if not np.allclose(non_epi_frequency.sum(axis=1), 1.0):
        raise ValueError("Non-epithelial frequencies do not sum to one per sample.")
    epi_positive = epi_counts.sum(axis=1) > 0
    if not np.allclose(epi_frequency.loc[epi_positive].sum(axis=1), 1.0):
        raise ValueError("Positive epithelial frequency rows do not sum to one.")
    non_epi_counts.to_csv(PREP_DIR / "non_epi_subtype_counts.csv")
    non_epi_frequency.to_csv(PREP_DIR / "non_epi_subtype_frequency.csv")
    epi_counts.to_csv(PREP_DIR / "epi_subtype_counts.csv")
    epi_frequency.to_csv(PREP_DIR / "epi_subtype_frequency.csv")

    sample_status = inclusion.loc[cm_samples, ["status", "raw_status", "series", "non_epi_cells", "epithelial_cells"]].copy()
    sample_status.index.name = "sample"
    sample_status.to_csv(PREP_DIR / "sample_status.csv")
    observed = set(sample_status["status"])
    if observed == {"tumor", "normal-like"}:
        mode, balanced, classification = "tumor_normal", True, True
    elif observed == {"tumor"}:
        mode, balanced, classification = "tumor_only", False, False
    elif observed == {"normal-like"}:
        mode, balanced, classification = "normal_only", False, False
    else:
        raise ValueError(f"Unsupported status set: {observed}")
    detection = (
        sample_status.groupby("status").size().rename("n_samples").reset_index()
        .rename(columns={"status": "canonical_status"})
    )
    detection.insert(0, "detected_mode", mode)
    detection["n_total_samples"] = sample_status.index.nunique()
    detection["status_source"] = "projected_atlas_obs_status"
    detection["status_balance_applied"] = balanced
    detection["cm_classification_available"] = classification
    detection.to_csv(PREP_DIR / "cm_status_mode_detection.csv", index=False)
    write_json(
        {
            "detected_mode": mode,
            "n_total_samples": int(sample_status.index.nunique()),
            "status_source": "projected_atlas_obs_status",
            "status_balance_applied": balanced,
            "cm_classification_available": classification,
            "raw_to_canonical_status": dict(
                sample_status[["raw_status", "status"]].drop_duplicates().itertuples(index=False, name=None)
            ),
            "included_samples": sample_status.index.astype(str).tolist(),
        },
        PREP_DIR / "cm_status_mode.json",
    )
    pd.DataFrame(
        [
            {
                "input_h5ad": str(INPUT_H5AD), "sample_col": "sample",
                "status_col": "status", "series_col": "series",
                "cell_type_col": "cell_type", "cell_subtype_col": "cell_subtype",
                "epithelial_cell_type": EPI_CELL_TYPE,
                "min_non_epi_cells_per_sample": MIN_NON_EPI_CELLS,
                "min_epi_cells_per_sample": MIN_EPI_CELLS,
                "detected_mode": mode, "n_candidate_samples": len(inclusion),
                "n_keep_for_cm": int(inclusion["keep_for_cm"].sum()),
                "n_keep_for_epi_cm": int(inclusion["keep_for_epi_cm"].sum()),
                "n_epi_subtypes": epi_frequency.shape[1],
                "n_non_epi_subtypes": non_epi_frequency.shape[1],
                "code_file": str(CODE_PATH), "seed": SEED,
            }
        ]
    ).to_csv(PREP_DIR / "input_preparation_parameters.csv", index=False)
    (PREP_DIR / "package_versions.txt").write_text(
        f"python={sys.version.split()[0]}\npandas={package_version('pandas')}\nanndata={package_version('anndata')}\ncode={CODE_PATH}\nseed={SEED}\n",
        encoding="utf-8",
    )
    completion = {
        "status": "completed", "detected_mode": mode,
        "n_candidate_samples": len(inclusion),
        "n_keep_for_cm": int(inclusion["keep_for_cm"].sum()),
        "n_keep_for_epi_cm": int(inclusion["keep_for_epi_cm"].sum()),
        "n_epi_subtypes": epi_frequency.shape[1],
        "n_non_epi_subtypes": non_epi_frequency.shape[1],
        "non_epi_frequency_rows_sum_to_one": True,
        "prefix_cell_type_conflicts": 0,
        "elapsed_seconds": time.time() - started,
    }
    write_json(completion, completion_path)
    print(json.dumps(completion, indent=2), flush=True)


if __name__ == "__main__":
    main()
