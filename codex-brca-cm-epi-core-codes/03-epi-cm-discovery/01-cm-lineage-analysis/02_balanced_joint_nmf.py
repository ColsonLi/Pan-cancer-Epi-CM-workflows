#!/usr/bin/env python3
"""Run canonical tumor-only K selection, NMF, and NNLS CM activity refit."""

from __future__ import annotations

import importlib.metadata
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from cm_analysis_common import (
    CM_K_RANGE, NMF_ALPHA_H, NMF_ALPHA_W, NMF_L1_RATIO, NMF_MAX_ITER,
    NMF_TOL, NMF_DIR, PREP_DIR, RANK_SELECTION_SEEDS, SEED,
    column_minmax_normalize, evaluate_cm_rank, fit_nmf,
    refit_activity_by_nnls, transform_cm_columns, write_json,
)


CODE_PATH = Path(__file__).resolve()


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def main() -> None:
    started = time.time()
    completion_path = NMF_DIR / "nmf_completion.json"
    if completion_path.exists():
        completion = json.loads(completion_path.read_text())
        if completion.get("status") == "completed":
            print(json.dumps({"status": "valid_existing_nmf_reused"}, indent=2))
            return
        raise FileExistsError("Existing NMF completion is invalid.")
    NMF_DIR.mkdir(parents=True, exist_ok=True)
    mode = json.loads((PREP_DIR / "cm_status_mode.json").read_text())["detected_mode"]
    if mode != "tumor_only":
        raise ValueError(f"This BRCA run expected tumor_only mode, found {mode}")
    frequency = pd.read_csv(PREP_DIR / "non_epi_subtype_frequency.csv", index_col=0)
    sample_status = pd.read_csv(PREP_DIR / "sample_status.csv", index_col=0)
    sample_status.index = sample_status.index.astype(str)
    frequency.index = frequency.index.astype(str)
    if not frequency.index.equals(sample_status.index):
        raise ValueError("NMF frequency and sample_status IDs/order differ.")
    if set(sample_status["status"].astype(str)) != {"tumor"}:
        raise ValueError("Tumor-only NMF received a non-tumor status.")

    normalized, minmax_params = column_minmax_normalize(frequency)
    normalized.to_csv(NMF_DIR / "non_epi_subtype_frequency_column_minmax.csv")
    minmax_params.to_csv(NMF_DIR / "column_minmax_params.csv")
    weights = pd.Series(1.0, index=normalized.index, name="joint_nmf_row_weight")
    weight_audit = weights.to_frame().join(sample_status[["status", "non_epi_cells"]])
    weight_audit["status_balance_applied"] = False
    weight_audit["weight_reason"] = "equal weights for tumor_only cohort"
    weight_audit.to_csv(NMF_DIR / "group_balanced_sample_weights.csv")
    v = normalized.to_numpy(dtype=float)
    metrics, selected_k = evaluate_cm_rank(v)
    metrics.to_csv(NMF_DIR / "joint_nmf_k_selection_metrics.csv", index=False)
    best_seed = int(metrics.loc[metrics["selected"], "best_seed"].iloc[0])
    _, h_raw, final_error, final_n_iter = fit_nmf(v, selected_k, best_seed)
    w_raw = refit_activity_by_nnls(h_raw, v)

    raw_cms = [f"CM{i + 1}" for i in range(selected_k)]
    final_cms = [f"t_CM{i + 1}" for i in range(selected_k)]
    w_raw_df = pd.DataFrame(w_raw, index=normalized.index, columns=raw_cms)
    h_raw_df = pd.DataFrame(h_raw, index=raw_cms, columns=normalized.columns)
    mapping = pd.DataFrame(
        {
            "raw_component": raw_cms,
            "raw_component_order": range(1, selected_k + 1),
            "CM": final_cms,
            "class": "tCM", "class_prefix": "t",
            "global_order": range(1, selected_k + 1),
            "classification_available": False,
            "classification_basis": "single_status_presence_fallback",
            "classification_reason": "normal-like samples absent",
        }
    )
    mapping.to_csv(NMF_DIR / "raw_to_canonical_CM_mapping.csv", index=False)
    w_df = w_raw_df.rename(columns=dict(zip(raw_cms, final_cms))).loc[:, final_cms]
    h_df = h_raw_df.rename(index=dict(zip(raw_cms, final_cms))).loc[final_cms]
    loading = h_df.T.copy()
    loading_fraction = loading.div(loading.sum(axis=0).replace(0, np.nan), axis=1).fillna(0.0)
    classification = mapping.copy()
    classification["tumor_mean"] = [float(w_df[cm].mean()) for cm in final_cms]
    classification["normal_like_mean"] = np.nan
    classification["delta"] = np.nan
    classification["p"] = np.nan
    classification["q"] = np.nan
    classification.to_csv(NMF_DIR / "joint_module_classification.csv", index=False)
    w_df.to_csv(NMF_DIR / "W_df.csv")
    h_df.to_csv(NMF_DIR / "H_df.csv")
    w_df.to_csv(NMF_DIR / "activity_df_sample_by_CM.csv")
    w_df.T.to_csv(NMF_DIR / "activity_df_CM_by_sample.csv")
    loading.to_csv(NMF_DIR / "loading_df_cell_subtype_by_CM.csv")
    loading_fraction.to_csv(NMF_DIR / "loading_df_cell_subtype_by_CM_fraction.csv")
    loading.to_csv(NMF_DIR / "balanced_joint_cm_subtype_loadings_raw_from_H_df.csv")
    loading_fraction.to_csv(NMF_DIR / "balanced_joint_cm_subtype_loadings_fraction_from_H_df.csv")
    for method in ["raw", "zscore", "robust", "standard_scale_col"]:
        transform_cm_columns(w_df, method).to_csv(
            NMF_DIR / f"w_df_activity_sample_by_CM_{method}.csv"
        )
        transform_cm_columns(loading, method).T.to_csv(
            NMF_DIR / f"h_df_loading_cell_subtype_by_CM_{method}.csv"
        )
    transform_cm_columns(w_df, "standard_scale_col").to_csv(
        NMF_DIR / "w_df_activity_sample_by_CM_standard_scale_col_with_series_plot_matrix.csv"
    )
    sample_status[["series", "status"]].rename(
        columns={"series": "Series", "status": "Status"}
    ).to_csv(NMF_DIR / "w_df_activity_sample_by_CM_standard_scale_col_with_series_row_annotations.csv")
    skipped = pd.DataFrame(
        {
            "output_or_step": [
                "status-balanced weighting",
                "status-derived sharedCM/normalCM/tCM classification",
                "activity_df_tumor_vs_normal_mean_sd_barplot",
                "normal-like nodeplots and correlation heatmaps",
                "normal-like Epi-CM associations",
                "tumor_centric_nodeplot_edge_origin",
                "Pearson optional association branch",
            ],
            "status": "skipped",
            "reason": [
                "tumor_only route uses equal weights",
                "normal-like samples absent; tumor presence fallback naming",
                "normal-like comparison status absent",
                "normal-like samples absent",
                "normal-like samples absent",
                "normal-like edge origin cannot be defined",
                "Pearson was not requested; Spearman is the default branch",
            ],
        }
    )
    skipped.to_csv(NMF_DIR / "single_status_skipped_outputs.csv", index=False)
    write_json(
        {
            "selected_module_k": selected_k, "best_seed": best_seed,
            "candidate_Ks": metrics["k"].astype(int).tolist(),
            "rank_selection_seeds": list(RANK_SELECTION_SEEDS),
            "detected_mode": mode, "present_status": "tumor",
            "absent_status": "normal-like", "status_balance_applied": False,
            "cm_classification_available": False,
            "final_cm_naming": "t_CM1..t_CMK in raw component order",
            "classification_basis": "single_status_presence_fallback",
        },
        NMF_DIR / "selected_module_k.json",
    )
    pd.DataFrame(
        [
            {
                "input_frequency": str(PREP_DIR / "non_epi_subtype_frequency.csv"),
                "detected_mode": mode, "n_samples": normalized.shape[0],
                "n_non_epi_subtypes": normalized.shape[1],
                "column_minmax": True, "equal_row_weights": True,
                "candidate_k_min": CM_K_RANGE[0], "candidate_k_max": CM_K_RANGE[1],
                "rank_selection_seeds": ";".join(map(str, RANK_SELECTION_SEEDS)),
                "selected_k": selected_k, "best_seed": best_seed,
                "nmf_max_iter": NMF_MAX_ITER, "nmf_tol": NMF_TOL,
                "nmf_alpha_w": NMF_ALPHA_W, "nmf_alpha_h": NMF_ALPHA_H,
                "nmf_l1_ratio": NMF_L1_RATIO, "final_nmf_error": final_error,
                "final_nmf_n_iter": final_n_iter,
                "activity_refit": "scipy.optimize.nnls on unweighted column-minmax V",
                "code_file": str(CODE_PATH), "seed": SEED,
            }
        ]
    ).to_csv(NMF_DIR / "nmf_parameters.csv", index=False)
    (NMF_DIR / "package_versions.txt").write_text(
        f"python={sys.version.split()[0]}\npandas={package_version('pandas')}\nscikit-learn={package_version('scikit-learn')}\nscipy={package_version('scipy')}\ncode={CODE_PATH}\nseed={SEED}\n",
        encoding="utf-8",
    )
    pd.DataFrame(
        [
            {"step": "K-selection NMF", "planned_backend": "CPU sklearn", "attempted_backend": "CPU sklearn", "status": "completed", "error_summary": "", "fallback_backend": "", "clean_input_reloaded": True, "final_backend_for_rerun": "CPU sklearn"},
            {"step": "NNLS activity refit", "planned_backend": "CPU scipy", "attempted_backend": "CPU scipy", "status": "completed", "error_summary": "", "fallback_backend": "", "clean_input_reloaded": True, "final_backend_for_rerun": "CPU scipy"},
        ]
    ).to_csv(NMF_DIR / "gpu_backend_capability_summary.csv", index=False)
    completion = {
        "status": "completed", "detected_mode": mode,
        "n_samples": normalized.shape[0], "n_non_epi_subtypes": normalized.shape[1],
        "candidate_k_count": len(metrics), "selected_k": selected_k,
        "best_seed": best_seed, "w_shape": list(w_df.shape),
        "h_shape": list(h_df.shape), "column_minmax_nonzero_columns": int((normalized.var(axis=0) > 0).sum()),
        "canonical_cm_names_applied": True, "nnls_activity_refit_complete": True,
        "elapsed_seconds": time.time() - started,
    }
    write_json(completion, completion_path)
    print(json.dumps(completion, indent=2), flush=True)


if __name__ == "__main__":
    main()
