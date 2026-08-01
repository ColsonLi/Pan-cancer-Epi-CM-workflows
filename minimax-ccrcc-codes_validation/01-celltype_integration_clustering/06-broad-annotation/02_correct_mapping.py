"""Correct cell_type mapping based on top DEG evidence (not argmax).

The auto-argmax in 01_broad_annotation.py mis-mapped:
  - C12: argmax -> T, but top DEGs are all ribosomal -> LowQuality
  - C13: argmax -> Endo, but only 15 cells with stress markers -> Doublets

This script:
  1. Reads adata_anno.h5ad
  2. Applies a manually curated cluster -> cell_type mapping
  3. Writes adata_anno.h5ad with corrected cell_type
  4. Also flags leiden_coarse 13 as low quality (low-quality / doublets)

Manual mapping rationale (top DEG evidence):
  C0  -> Endo    (EPAS1, FLT1, PLVAP, TCF4)
  C1  -> Epi     (PDZK1IP1, FXYD2, CRYAB)
  C2  -> S       (CALD1, RGS5, MGP, BGN, ADIRF)
  C3  -> Mast    (TPSB2, TPSAB1, CPA3, KIT, MS4A2)
  C4  -> T       (STMN1, MKI67 + CD3D/E) - cycling T
  C5  -> B       (CD79A, MS4A1, BANK1, IGHM)
  C6  -> Mye     (HLA-DRA, CD74, TYROBP, CST3, AIF1) - DC/mono
  C7  -> NK      (NKG7, KLRD1, GNLY, GZMB)
  C8  -> Epi     (CD24, FXYD2, SPP1) - proximal tubule
  C9  -> T       (CCL5, NKG7, CD3D, GZMK, TRAC) - effector T
  C10 -> T       (IL7R, CD3D, CD2, TRAC) - naive T
  C11 -> Mye     (TYROBP, LST1, S100A4, FCER1G) - monocyte
  C12 -> LowQuality (RPS27, RPL41, RPL13A, RPS15A) - ribosomal/dying
  C13 -> Doublets  (15 cells, NDUFA4L2/HSPA1A stress signature)
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

SEED = 42
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
IN_H5AD = ROOT / "epi-cm-core-workflow/h5ad/01-celltype_integration_clustering/06-broad-annotation/adata_anno.h5ad"
OUT_H5AD = IN_H5AD  # in-place
TAB = ROOT / "epi-cm-core-workflow/tables/01-celltype_integration_clustering/06-broad-annotation"
FIG = ROOT / "epi-cm-core-workflow/figures/01-celltype_integration_clustering/06-broad-annotation"

MANUAL_MAP = {
    "0": "Endo",
    "1": "Epi",
    "2": "S",
    "3": "Mast",
    "4": "T",
    "5": "B",
    "6": "Mye",
    "7": "NK",
    "8": "Epi",
    "9": "T",
    "10": "T",
    "11": "Mye",
    "12": "LowQuality",
    "13": "Doublets",
}

# Subtype labels (for B-cell cluster C5 - check if it has plasma cells; for now keep B)
# Add cell_subtype init empty
CELLTYPE_RATIONALE = {
    "Endo":      "C0: EPAS1, FLT1, PLVAP, TCF4",
    "Epi":       "C1, C8: PDZK1IP1, FXYD2, CD24, SPP1 (kidney tubule epithelium)",
    "S":         "C2: CALD1, RGS5, MGP, BGN (pericyte/stromal)",
    "Mast":      "C3: TPSB2, TPSAB1, CPA3, KIT",
    "T":         "C4, C9, C10: CD3D/E + various states (cycling C4; effector C9; naive C10)",
    "B":         "C5: CD79A, MS4A1, BANK1, IGHM (B cells; pDC not seen at this resolution)",
    "Mye":       "C6, C11: HLA-DRA, CD74, TYROBP, LST1 (DC + monocyte)",
    "NK":        "C7: NKG7, KLRD1, GNLY, GZMB",
    "LowQuality": "C12: top DEGs all ribosomal / low complexity",
    "Doublets":  "C13: 15 cells, stress signature NDUFA4L2/HSPA1A",
}


def main() -> None:
    t0 = time.time()
    print(f"[corr] reading {IN_H5AD.name}…", flush=True)
    adata = sc.read_h5ad(IN_H5AD)
    print(f"[corr] loaded {adata.shape} in {time.time()-t0:.1f}s", flush=True)

    # Apply manual mapping
    auto_ct = adata.obs["cell_type"].astype(str).copy()
    adata.obs["cell_type_auto"] = auto_ct.values  # preserve
    adata.obs["cell_type"] = adata.obs["leiden_coarse"].astype(str).map(MANUAL_MAP).astype("category")

    print(f"[corr] corrected cell_type counts:")
    print(adata.obs["cell_type"].value_counts(), flush=True)

    # Quality summary
    n_total = adata.n_obs
    n_keep = int((~adata.obs["cell_type"].isin(["LowQuality", "Doublets"])).sum())
    print(f"[corr] retained: {n_keep}/{n_total} ({100*n_keep/n_total:.2f}%)", flush=True)

    # Write corrected h5ad (overwrite)
    print(f"[corr] writing corrected {OUT_H5AD.name}…", flush=True)
    adata.write_h5ad(OUT_H5AD, compression="gzip")

    # Save updated map
    pd.DataFrame([
        {"cluster": cl, "cell_type": ct, "rationale": CELLTYPE_RATIONALE.get(ct, "")}
        for cl, ct in MANUAL_MAP.items()
    ]).to_csv(TAB / "cluster_to_celltype_map.csv", index=False)
    print(f"[corr] cluster_to_celltype_map.csv updated", flush=True)

    # Re-plot UMAP with corrected cell_type
    sc.settings.figdir = str(FIG)
    sc.set_figure_params(figsize=(5, 5), dpi=150)
    sc.pl.umap(adata, color="cell_type", save="_cell_type_corrected.pdf", show=False)

    # Per-celltype cell counts
    grp = adata.obs.groupby(["cell_type", "status"], observed=True).size().unstack(fill_value=0)
    grp.to_csv(TAB / "celltype_status_counts.csv")

    grp2 = adata.obs.groupby(["cell_type", "series"], observed=True).size().unstack(fill_value=0)
    grp2.to_csv(TAB / "celltype_series_counts.csv")

    # Per-cluster cell counts (final)
    pd.DataFrame(
        adata.obs.groupby(["leiden_coarse", "cell_type"], observed=True).size().rename("n_cells")
    ).reset_index().to_csv(TAB / "cluster_celltype_counts.csv", index=False)

    print(f"[corr] done. total={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()