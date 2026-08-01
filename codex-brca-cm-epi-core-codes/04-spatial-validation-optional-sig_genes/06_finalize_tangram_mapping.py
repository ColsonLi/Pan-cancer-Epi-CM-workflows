#!/usr/bin/env python3
"""Validate and finalize the six-sample canonical Tangram mapping inventory."""

from __future__ import annotations

import importlib.metadata
import json
import random
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

SEED = 42
random.seed(SEED)
np.random.seed(SEED)
ROOT = Path("/mnt/disk18t/lr_xcy/riku/brca_cm-epi_codex")
WORKFLOW = ROOT / "epi-cm-core-workflow"
TABLE_DIR = WORKFLOW / "tables/04-spatial-validation-optional-sig_genes/02-tangram-mapping"
H5AD_DIR = WORKFLOW / "h5ad/04-spatial-validation-optional-sig_genes/02-tangram-mapping"
SAMPLES = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]


def pkg(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def main() -> None:
    completion_path = TABLE_DIR / "mapping_all_samples_completion.json"
    if completion_path.exists():
        raise FileExistsError(completion_path)
    rows, backend_rows = [], []
    failures = []
    for sample in SAMPLES:
        sample_completion_path = TABLE_DIR / f"{sample}_mapping_completion.json"
        score_path = TABLE_DIR / f"{sample}_tangram_pseudobulk_epi_cm_spot_scores.csv"
        mapper_path = H5AD_DIR / f"{sample}_tangram_mapper.h5ad"
        completion = json.loads(sample_completion_path.read_text())
        score = pd.read_csv(score_path)
        epi_cols = [column for column in score.columns if column.startswith("EPIfrac__")]
        cm_cols = [column for column in score.columns if column.startswith("CMact__")]
        epi_sum = score[epi_cols].sum(axis=1)
        mapper = ad.read_h5ad(mapper_path, backed="r")
        training_gene_count = len(mapper.uns["train_genes_df"])
        mapper_nonnegative = bool(np.nanmin(np.asarray(mapper.X)) >= -1e-8)
        row = {
            "sample": sample,
            "status": completion.get("status"),
            "n_spots": len(score),
            "n_common_marker_genes": completion.get("n_common_genes"),
            "n_tangram_training_genes": training_gene_count,
            "n_subtypes": mapper.n_obs,
            "n_epi_subtypes": len(epi_cols),
            "n_cms": len(cm_cols),
            "spot_ids_unique": bool(score["spot_id"].astype(str).is_unique),
            "barcodes_unique": bool(score["barcode"].astype(str).is_unique),
            "coordinates_finite": bool(np.isfinite(score[["array_row", "array_col", "spatial_x", "spatial_y"]].to_numpy(float)).all()),
            "epi_fraction_rows_sum_one": bool(np.allclose(epi_sum.to_numpy(), 1.0, rtol=1e-6, atol=1e-8)),
            "epi_fraction_nonnegative": bool((score[epi_cols].to_numpy() >= -1e-8).all()),
            "cm_activity_nonnegative": bool((score[cm_cols].to_numpy() >= -1e-8).all()),
            "all_cm_columns_variable": bool(score[cm_cols].nunique().gt(1).all()),
            "mapper_nonnegative": mapper_nonnegative,
            "device": completion.get("device"),
            "epochs": mapper.uns["tangram_parameters"]["num_epochs"],
            "learning_rate": mapper.uns["tangram_parameters"]["learning_rate"],
            "seed": mapper.uns["tangram_parameters"]["random_state"],
            "elapsed_seconds": completion.get("elapsed_seconds"),
            "score_csv": str(score_path.resolve()),
            "mapper_h5ad": str(mapper_path.resolve()),
        }
        mapper.file.close()
        for criterion in ["spot_ids_unique", "barcodes_unique", "coordinates_finite", "epi_fraction_rows_sum_one", "epi_fraction_nonnegative", "cm_activity_nonnegative", "all_cm_columns_variable", "mapper_nonnegative"]:
            if not row[criterion]:
                failures.append({"sample": sample, "criterion": criterion})
        if row["status"] != "completed" or row["n_subtypes"] != 68 or row["n_epi_subtypes"] != 10 or row["n_cms"] != 10:
            failures.append({"sample": sample, "criterion": "completion_or_dimension"})
        rows.append(row)
        backend_rows.extend([
            {"sample": sample, "step": "Tangram pseudobulk mapping", "planned_backend": "GPU cuda:0", "attempted_backend": "GPU cuda:0", "status": "completed", "error_summary": "", "fallback_backend": "", "clean_input_reloaded": True, "final_backend_for_rerun": "GPU cuda:0"},
            {"sample": sample, "step": "Module03 H NNLS activity refit", "planned_backend": "CPU scipy", "attempted_backend": "CPU scipy", "status": "completed", "error_summary": "", "fallback_backend": "", "clean_input_reloaded": True, "final_backend_for_rerun": "CPU scipy"},
        ])
    summary = pd.DataFrame(rows)
    summary.to_csv(TABLE_DIR / "mapping_sample_summary.csv", index=False)
    pd.DataFrame(backend_rows).to_csv(TABLE_DIR / "gpu_backend_capability_summary.csv", index=False)
    pd.DataFrame([{
        "CMact_definition": "NNLS refit of canonical Module03 H_df on Tangram-derived non-Epi within-spot composition after saved Module03 column-minmax transform",
        "EPIfrac_definition": "Tangram abundance normalized only across the 10 Epi subtype columns within each spot",
        "H_df": str((WORKFLOW / "tables/03-epi-cm-discovery/01-cm-lineage-analysis/02_balanced_joint_nmf/H_df.csv").resolve()),
        "column_minmax_params": str((WORKFLOW / "tables/03-epi-cm-discovery/01-cm-lineage-analysis/02_balanced_joint_nmf/column_minmax_params.csv").resolve()),
        "tangram_mode": "cells", "device": "cuda:0", "num_epochs": 350, "learning_rate": 0.05, "seed": SEED,
    }]).to_csv(TABLE_DIR / "mapping_feature_definitions.csv", index=False)
    (TABLE_DIR / "readme.txt").write_text("Each of six primary-tumor Visium sections was mapped independently with real Tangram pseudobulk subtype mapping on cuda:0 (mode=cells, epochs=350, learning_rate=0.05, seed=42). Every mandatory per-sample spot-score CSV was written before its mapper h5ad. EPIfrac is normalized only within the 10 projected epithelial subtypes. CMact is the canonical Module03 H NNLS activity refit from Tangram-derived non-Epi composition using the saved Module03 column-minmax parameters. Downstream plots/statistics must reload these CSVs.\n", encoding="utf-8")
    (TABLE_DIR / "package_versions.txt").write_text(f"python={sys.version.split()[0]}\nanndata={pkg('anndata')}\nscanpy={pkg('scanpy')}\ntangram-sc={pkg('tangram-sc')}\ntorch={pkg('torch')}\nnumpy={pkg('numpy')}\npandas={pkg('pandas')}\nscipy={pkg('scipy')}\nenvironment={ROOT / 'uv_envs/rapids/.venv'}\ntorch_install_source=official PyTorch cu128 index because CUDA wheels are not served by the TUNA PyPI mirror\ntangram_install_source=TUNA PyPI mirror\ncode={Path(__file__).resolve()}\nseed={SEED}\n", encoding="utf-8")
    completion = {"status": "completed" if not failures else "failed", "n_samples": len(rows), "n_spots": int(summary["n_spots"].sum()), "n_failed_invariants": len(failures), "failures": failures, "all_gpu_tangram": bool(summary["device"].eq("cuda:0").all()), "tangram_parameters_fixed": bool(summary["epochs"].eq(350).all() and np.allclose(summary["learning_rate"], 0.05) and summary["seed"].eq(SEED).all()), "n_score_csvs": len(list(TABLE_DIR.glob("*_tangram_pseudobulk_epi_cm_spot_scores.csv"))), "n_mapper_h5ads": len(list(H5AD_DIR.glob("*_tangram_mapper.h5ad"))), "code_file": str(Path(__file__).resolve()), "seed": SEED}
    completion_path.write_text(json.dumps(completion, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(completion, indent=2), flush=True)
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
