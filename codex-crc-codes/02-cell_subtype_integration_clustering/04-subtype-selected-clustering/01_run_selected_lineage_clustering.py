#!/usr/bin/env python3
"""Rerun user-selected subtype clustering from saved lineage Harmony h5ad."""

from __future__ import annotations

import argparse
import gc
import platform
import random
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import scanpy as sc

matplotlib.use("Agg")

SEED = 42
random.seed(SEED)
np.random.seed(SEED)

WORKFLOW_ROOT = Path("/mnt/disk18t/lr_xcy/riku/crc_val/epi-cm-core-workflow")
MODULE = "02-cell_subtype_integration_clustering"
STEP = "04-subtype-selected-clustering"
CODE_FILE = Path(__file__)
BASIS = "X_pca_inte"
NEIGHBORS_KEY = "neighbors"
UMAP_KEY = "X_umap"

LINEAGE_CONFIG = {
    "b_cells": {"abbr": "b", "label": "B Cells", "pcs": 40, "nn": 40, "res": 0.5},
    "cycling_immune": {"abbr": "cycling", "label": "Cycling Immune Cells", "pcs": 35, "nn": 40, "res": 0.1},
    "endothelial": {"abbr": "endo", "label": "Endothelial Cells", "pcs": 20, "nn": 30, "res": 0.1},
    "epithelial": {"abbr": "epi", "label": "Epithelial Cells", "pcs": 30, "nn": 30, "res": 0.3},
    "mast": {"abbr": "mast", "label": "Mast Cells", "pcs": 20, "nn": 20, "res": 0.3},
    "myeloid": {"abbr": "mye", "label": "Myeloid Cells", "pcs": 40, "nn": 40, "res": 0.2},
    "plasma_cells": {"abbr": "plasma", "label": "Plasma Cells", "pcs": 25, "nn": 40, "res": 0.4},
    "schwann": {"abbr": "schwann", "label": "Schwann Cells", "pcs": 10, "nn": 30, "res": 0.3},
    "stromal": {"abbr": "stromal", "label": "Stromal Cells", "pcs": 20, "nn": 40, "res": 0.3},
    "t_cells": {"abbr": "t", "label": "T Cells", "pcs": 20, "nn": 40, "res": 0.4},
}

STATUS_MAP = {
    "primary_crc_no_metastasis_no_lymph_node": "tumor",
    "normal_colon_no_metastasis_no_lymph_node": "normal",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lineage", choices=sorted(LINEAGE_CONFIG), required=True)
    parser.add_argument("--allow-existing", action="store_true")
    return parser.parse_args()


def res_token(resolution: float) -> str:
    return str(resolution).replace(".", "p")


def main() -> None:
    args = parse_args()
    cfg = LINEAGE_CONFIG[args.lineage]
    n_pcs = int(cfg["pcs"])
    n_neighbors = int(cfg["nn"])
    resolution = float(cfg["res"])
    cluster_key = f"leiden_res{res_token(resolution)}"

    input_h5ad = (
        WORKFLOW_ROOT
        / "h5ad"
        / MODULE
        / "02-subtype-harmony"
        / args.lineage
        / f"adata_{args.lineage}_harmony.h5ad"
    )
    output_h5ad = (
        WORKFLOW_ROOT
        / "h5ad"
        / MODULE
        / STEP
        / args.lineage
        / f"adata_{cfg['abbr']}_selected_clustered.h5ad"
    )
    table_dir = WORKFLOW_ROOT / "tables" / MODULE / STEP / args.lineage
    figure_dir = WORKFLOW_ROOT / "figures" / MODULE / STEP / args.lineage
    selected_table = table_dir / "selected_clustering_parameters.csv"
    counts_table = table_dir / f"{cluster_key}_cluster_counts.csv"
    figure_pdf = figure_dir / f"umap_{cluster_key}_selected.pdf"

    outputs = [output_h5ad, selected_table, counts_table, figure_pdf]
    existing = [str(path) for path in outputs if path.exists()]
    if existing and not args.allow_existing:
        raise FileExistsError("Selected clustering output already exists; refusing to overwrite:\n" + "\n".join(existing))
    if args.allow_existing and all(path.exists() for path in outputs):
        print(f"{args.lineage}: selected clustering already exists")
        return

    output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    sc.settings.autoshow = False
    sc.settings.figdir = str(figure_dir)
    sc.settings.set_figure_params(figsize=(3, 3), dpi=150)

    import cupy as cp
    import rapids_singlecell as rsc

    adata = ad.read_h5ad(input_h5ad)
    if BASIS not in adata.obsm:
        raise KeyError(f"Missing obsm[{BASIS!r}] in {input_h5ad}")
    if adata.raw is None:
        raise RuntimeError(f"Missing raw in {input_h5ad}")

    if "status" in adata.obs.columns:
        adata.obs["status_original"] = adata.obs["status"].astype(str).values
        mapped = adata.obs["status_original"].map(STATUS_MAP)
        adata.obs["status"] = mapped.fillna(adata.obs["status_original"]).astype("category")

    rsc.get.anndata_to_GPU(adata)
    rsc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        use_rep=BASIS,
        random_state=SEED,
        key_added=NEIGHBORS_KEY,
    )
    rsc.tl.umap(adata, random_state=SEED, neighbors_key=NEIGHBORS_KEY, key_added=UMAP_KEY)
    rsc.tl.leiden(
        adata,
        resolution=resolution,
        key_added=cluster_key,
        random_state=SEED,
        neighbors_key=NEIGHBORS_KEY,
    )
    rsc.get.anndata_to_CPU(adata)
    cp.get_default_memory_pool().free_all_blocks()

    counts = adata.obs[cluster_key].astype(str).value_counts().rename_axis("cluster").reset_index(name="n_cells")
    counts.to_csv(counts_table, index=False)
    sc.pl.umap(adata, color=cluster_key, save=f"_{cluster_key}_selected.pdf", show=False)

    adata.write_h5ad(output_h5ad)
    params = pd.DataFrame(
        [
            {
                "lineage": args.lineage,
                "target_leiden_coarse": cfg["label"],
                "cell_abbrev": cfg["abbr"],
                "source_harmony_h5ad": str(input_h5ad),
                "selected_h5ad": str(output_h5ad),
                "n_obs": int(adata.n_obs),
                "n_vars": int(adata.n_vars),
                "n_pcs": n_pcs,
                "n_neighbors": n_neighbors,
                "resolution": resolution,
                "basis": BASIS,
                "neighbors_key": NEIGHBORS_KEY,
                "umap_key": UMAP_KEY,
                "cluster_key": cluster_key,
                "n_clusters": int(counts.shape[0]),
                "selected_by_user": True,
                "status_values_standardized": True,
                "status_mapping": ";".join(f"{k}->{v}" for k, v in STATUS_MAP.items()),
                "backend_package": "rapids_singlecell",
                "random_seed": SEED,
                "code_file": str(CODE_FILE),
                "selected_umap_pdf": str(figure_pdf),
                "cluster_count_table": str(counts_table),
            }
        ]
    )
    params.to_csv(selected_table, index=False)

    with (table_dir / "package_versions.txt").open("w") as fh:
        fh.write(f"python: {platform.python_version()}\n")
        fh.write(f"anndata: {ad.__version__}\n")
        fh.write(f"scanpy: {sc.__version__}\n")
        fh.write(f"numpy: {np.__version__}\n")
        fh.write(f"pandas: {pd.__version__}\n")
        try:
            import rapids_singlecell as rsc_version

            fh.write(f"rapids_singlecell: {rsc_version.__version__}\n")
        except Exception:
            fh.write("rapids_singlecell: unavailable_version\n")
        fh.write(f"code_file: {CODE_FILE}\n")

    with (table_dir / "readme.txt").open("w") as fh:
        fh.write(f"{args.lineage} selected subtype clustering completed.\n")
        fh.write(f"Input Harmony h5ad: {input_h5ad}\n")
        fh.write(f"Selected clustering h5ad: {output_h5ad}\n")
        fh.write(f"User-selected parameters: pcs={n_pcs}, n_neighbors={n_neighbors}, resolution={resolution}\n")
        fh.write(f"Raw cluster key: {cluster_key}\n")
        fh.write("Status values standardized to tumor/normal in this selected object.\n")

    print(params.to_string(index=False))
    del adata
    gc.collect()


if __name__ == "__main__":
    main()
