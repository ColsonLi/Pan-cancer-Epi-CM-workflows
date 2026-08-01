"""Block 03 step 1: Derive Frequency Matrices (canonical per SKILL.md).

Uses the exact canonical helpers from the skill:
  - build_non_epi_and_epi_frequencies
  - column_minmax_normalize
  - sample_status with status_majority_fraction and non_epi_cells

Outputs:
  non_epi_subtype_counts.csv
  non_epi_subtype_frequency.csv
  non_epi_subtype_frequency_column_minmax.csv
  column_minmax_params.csv
  epi_subtype_counts.csv
  epi_subtype_frequency.csv
  sample_status.csv
  sample_compartment_cell_counts.csv
  sample_inclusion_exclusion.csv
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

SEED = 42
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
ATLAS = ROOT / "epi-cm-core-workflow/h5ad/02-cell_subtype_integration_clustering/06-ccrcc-subtypes-to-full-adata/adata_anno_cellsubtype.h5ad"
TAB = ROOT / "epi-cm-core-workflow/tables/03-epi-cm-discovery"
TAB.mkdir(parents=True, exist_ok=True)

MIN_NON_EPI_CELLS = 50
MIN_EPI_CELLS_ASSOC = 1


def build_non_epi_and_epi_frequencies(obs: pd.DataFrame, sample_col: str, status_col: str,
                                    celltype_col: str, subtype_col: str, epi_label: str,
                                    min_non_epi_cells_per_sample: int):
    """Canonical helper: split obs into epi / non_epi, compute counts + frequency."""
    required = {sample_col, status_col, celltype_col, subtype_col}
    missing = required - set(obs.columns)
    if missing:
        raise ValueError(f"Missing required obs columns: {sorted(missing)}")

    clean = obs.dropna(subset=list(required)).copy()
    is_epi = clean[celltype_col].astype(str).eq(epi_label)
    non_epi = clean.loc[~is_epi].copy()
    epi = clean.loc[is_epi].copy()
    if non_epi.empty:
        raise ValueError("No non-epithelial cells remain for CM NMF.")

    non_epi_cells = non_epi.groupby(sample_col, observed=True).size().rename("non_epi_cells")
    keep_samples = non_epi_cells.index[non_epi_cells >= min_non_epi_cells_per_sample]
    non_epi = non_epi.loc[non_epi[sample_col].isin(keep_samples)].copy()
    if non_epi.empty:
        raise ValueError("No samples pass min_non_epi_cells_per_sample.")

    non_epi_counts = pd.crosstab(non_epi[sample_col], non_epi[subtype_col]).astype(float)
    non_epi_counts = non_epi_counts.sort_index(axis=0).sort_index(axis=1)
    non_epi_frequency = non_epi_counts.div(non_epi_counts.sum(axis=1), axis=0).fillna(0.0)

    status_counts = pd.crosstab(non_epi[sample_col], non_epi[status_col]).reindex(non_epi_frequency.index).fillna(0)
    sample_status = pd.DataFrame(index=non_epi_frequency.index)
    sample_status["status"] = status_counts.idxmax(axis=1)
    sample_status["status_majority_fraction"] = status_counts.max(axis=1) / status_counts.sum(axis=1)
    sample_status["non_epi_cells"] = non_epi_counts.sum(axis=1).astype(int)

    epi_counts = pd.crosstab(epi[sample_col], epi[subtype_col]).astype(float)
    epi_counts = epi_counts.reindex(non_epi_frequency.index).fillna(0.0)
    epi_frequency = epi_counts.div(epi_counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    epi_frequency = epi_frequency.sort_index(axis=0).sort_index(axis=1)
    return non_epi_counts, non_epi_frequency, epi_counts, epi_frequency, sample_status


def main() -> None:
    t0 = time.time()
    print(f"[freq] reading {ATLAS.name}…", flush=True)
    adata = sc.read_h5ad(ATLAS)
    obs = adata.obs.copy()
    obs["sample_id"] = obs["sample"].astype(str)
    obs["status"] = obs["status"].astype(str)
    obs["cell_type"] = obs["cell_type"].astype(str)
    obs["cell_subtype"] = obs["cell_subtype"].astype(str)
    print(f"[freq] loaded {adata.shape}, time={time.time()-t0:.1f}s", flush=True)

    non_epi_counts, non_epi_frequency, epi_counts, epi_frequency, sample_status = (
        build_non_epi_and_epi_frequencies(
            obs, "sample_id", "status", "cell_type", "cell_subtype", "Epi",
            min_non_epi_cells_per_sample=MIN_NON_EPI_CELLS,
        )
    )
    print(f"[freq] samples: {len(sample_status)}, "
          f"non_epi sub: {non_epi_counts.shape[1]}, epi sub: {epi_counts.shape[1]}", flush=True)
    print(f"[freq] status: {sample_status['status'].value_counts().to_dict()}", flush=True)

    # Compartment counts
    compartment = pd.DataFrame(index=sample_status.index)
    compartment["non_epi_cells"] = sample_status["non_epi_cells"].astype(int)
    compartment["epi_cells"] = epi_counts.sum(axis=1).reindex(compartment.index).fillna(0).astype(int)
    compartment.to_csv(TAB / "sample_compartment_cell_counts.csv")
    print(f"[freq] compartment counts (median): {compartment.median().to_dict()}", flush=True)

    # Inclusion / exclusion
    inclusion = compartment.copy()
    inclusion["has_status"] = inclusion.index.isin(sample_status.index)
    inclusion["has_non_epi"] = inclusion["non_epi_cells"] >= MIN_NON_EPI_CELLS
    inclusion["has_epi_for_association"] = inclusion["epi_cells"] >= MIN_EPI_CELLS_ASSOC
    inclusion["keep_for_cm"] = inclusion["has_status"] & inclusion["has_non_epi"]
    inclusion["keep_for_epi_cm"] = inclusion["keep_for_cm"] & inclusion["has_epi_for_association"]

    def exclusion_reason(row):
        rs = []
        if not row["has_status"]:
            rs.append("missing_status")
        if not row["has_non_epi"]:
            rs.append(f"non_epi_cells<{MIN_NON_EPI_CELLS}")
        if not row["has_epi_for_association"]:
            rs.append(f"epi_cells<{MIN_EPI_CELLS_ASSOC}")
        return ";".join(rs) if rs else "kept"

    inclusion["exclusion_reason"] = inclusion.apply(exclusion_reason, axis=1)
    inclusion.to_csv(TAB / "sample_inclusion_exclusion.csv")
    print(f"[freq] keep_for_cm={inclusion['keep_for_cm'].sum()}, keep_for_epi_cm={inclusion['keep_for_epi_cm'].sum()}", flush=True)

    # Save raw counts and frequency
    non_epi_counts.astype(int).to_csv(TAB / "non_epi_subtype_counts.csv")
    non_epi_frequency.to_csv(TAB / "non_epi_subtype_frequency.csv")
    epi_counts.astype(int).to_csv(TAB / "epi_subtype_counts.csv")
    epi_frequency.to_csv(TAB / "epi_subtype_frequency.csv")
    sample_status.to_csv(TAB / "sample_status.csv")

    # Column-minmax normalize (on samples passing CM filter)
    keep_cm_samples = inclusion.index[inclusion["keep_for_cm"]]
    non_epi_freq_cm = non_epi_frequency.loc[keep_cm_samples]
    column_min = non_epi_freq_cm.min(axis=0)
    column_range = (non_epi_freq_cm.max(axis=0) - column_min).replace(0, np.nan)
    V_column_minmax = ((non_epi_freq_cm - column_min) / column_range).fillna(0.0).clip(lower=0.0)
    if (V_column_minmax.var(axis=0) > 0).sum() == 0:
        raise ValueError("Column-minmax produced zero valid non-epithelial subtype columns.")
    V_column_minmax.to_csv(TAB / "non_epi_subtype_frequency_column_minmax.csv")
    pd.DataFrame({"min": column_min, "range": column_range.fillna(0.0)}).rename_axis("cell_subtype").to_csv(
        TAB / "column_minmax_params.csv"
    )
    print(f"[freq] V_column_minmax: {V_column_minmax.shape}", flush=True)

    summary = {
        "n_samples_total": int(compartment.shape[0]),
        "n_samples_keep_for_cm": int(inclusion["keep_for_cm"].sum()),
        "n_samples_keep_for_epi_cm": int(inclusion["keep_for_epi_cm"].sum()),
        "n_epi_subtypes": int(epi_counts.shape[1]),
        "n_non_epi_subtypes": int(non_epi_counts.shape[1]),
        "min_non_epi_per_sample": MIN_NON_EPI_CELLS,
        "min_epi_per_sample_association": MIN_EPI_CELLS_ASSOC,
        "atlas_cells": int(adata.n_obs),
    }
    with open(TAB / "freq_matrices_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[freq] done. total={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()