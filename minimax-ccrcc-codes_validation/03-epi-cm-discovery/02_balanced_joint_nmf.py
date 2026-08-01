"""Block 03 step 2: Balanced Joint NMF (canonical per SKILL.md).

Strictly follows skill:
  - candidate K range = 2 to 20 inclusive
  - random seeds = 0, 1, 2, 3, 4
  - input = group-balanced weighted column-minmax non-epithelial frequency
  - V_weighted = sqrt(weight)[:, None] * V_column_minmax
  - selection_score = best_balanced_explained_fraction + 0.05*stability - 0.01*K
  - final fit with selected K + best_seed; NNLS refit on unweighted V to fixed H
  - raw CMs classified as sharedCM/normalCM/tCM; canonical names s_CM/n_CM/t_CM
"""
from __future__ import annotations

import json
import math
import time
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.optimize import nnls, linear_sum_assignment
from sklearn.decomposition import NMF
from sklearn.metrics.pairwise import cosine_similarity

SEED = 42
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
TAB = ROOT / "epi-cm-core-workflow/tables/03-epi-cm-discovery"
V_CSV = TAB / "non_epi_subtype_frequency_column_minmax.csv"
SAMPLE_STATUS_CSV = TAB / "sample_status.csv"
SAMPLE_INCL_CSV = TAB / "sample_inclusion_exclusion.csv"


@dataclass(frozen=True)
class CMNMFConfig:
    sample_col: str = "sample_id"
    status_col: str = "status"
    celltype_col: str = "cell_type"
    subtype_col: str = "cell_subtype"
    epi_label: str = "Epi"
    normal_status: str = "normal-like"
    tumor_status: str = "tumor"
    min_non_epi_cells_per_sample: int = 50
    forced_cm_k: int | None = None
    cm_k_range: tuple = (2, 20)
    rank_selection_seeds: tuple = (0, 1, 2, 3, 4)
    nmf_max_iter: int = 3000
    nmf_tol: float = 1e-5
    nmf_alpha_w: float = 0.0
    nmf_alpha_h: float = 1e-3
    nmf_l1_ratio: float = 0.1
    normal_group_total_weight: float = 0.5
    tumor_group_total_weight: float = 0.5
    normal_specific_max_ratio: float = 0.5
    tumor_specific_min_ratio: float = 2.0
    min_active_fraction_for_specific: float = 0.05


def write_json(obj: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2, ensure_ascii=False)


def column_minmax_normalize(freq: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    column_min = freq.min(axis=0)
    column_range = (freq.max(axis=0) - column_min).replace(0, np.nan)
    V = ((freq - column_min) / column_range).fillna(0.0).clip(lower=0.0)
    params = pd.DataFrame({"min": column_min, "range": column_range.fillna(0.0)})
    params.index.name = "cell_subtype"
    if (V.var(axis=0) > 0).sum() == 0:
        raise ValueError("Column-minmax matrix has no non-zero-variance subtype columns.")
    return V.astype(float), params


def make_status_balanced_weights(sample_status: pd.DataFrame, cfg: CMNMFConfig) -> pd.Series:
    status = sample_status["status"].astype(str)
    normal_mask = status.eq(cfg.normal_status)
    tumor_mask = status.eq(cfg.tumor_status)
    if normal_mask.sum() == 0 or tumor_mask.sum() == 0:
        raise ValueError("Need both normal-like and tumor samples for status-balanced CM NMF.")
    weights = pd.Series(0.0, index=sample_status.index, name="joint_nmf_row_weight")
    weights.loc[normal_mask] = cfg.normal_group_total_weight / normal_mask.sum()
    weights.loc[tumor_mask] = cfg.tumor_group_total_weight / tumor_mask.sum()
    return weights / weights.mean()


def make_nmf(n_components: int, seed: int, cfg: CMNMFConfig) -> NMF:
    return NMF(
        n_components=n_components,
        init="nndsvda",
        solver="cd",
        beta_loss="frobenius",
        random_state=seed,
        max_iter=cfg.nmf_max_iter,
        tol=cfg.nmf_tol,
        alpha_W=cfg.nmf_alpha_w,
        alpha_H=cfg.nmf_alpha_h,
        l1_ratio=cfg.nmf_l1_ratio,
    )


def row_l1_normalize(matrix: np.ndarray) -> np.ndarray:
    denom = matrix.sum(axis=1, keepdims=True)
    return np.divide(matrix, denom, out=np.zeros_like(matrix, dtype=float), where=denom != 0)


def matched_cosine_mean(H_a: np.ndarray, H_b: np.ndarray) -> float:
    sim = cosine_similarity(row_l1_normalize(H_a), row_l1_normalize(H_b))
    r, c = linear_sum_assignment(-sim)
    return float(sim[r, c].mean())


def stability_matched_cosine(H_by_seed: list[np.ndarray]) -> float:
    if len(H_by_seed) < 2:
        return float("nan")
    values = [
        matched_cosine_mean(H_by_seed[i], H_by_seed[j])
        for i in range(len(H_by_seed)) for j in range(i + 1, len(H_by_seed))
    ]
    return float(np.mean(values))


def evaluate_cm_rank(V_weighted: np.ndarray, cfg: CMNMFConfig) -> tuple[pd.DataFrame, int]:
    max_rank = min(cfg.cm_k_range[1], V_weighted.shape[0], V_weighted.shape[1])
    min_rank = max(2, cfg.cm_k_range[0])
    if max_rank < min_rank:
        raise ValueError(f"Invalid CM K range after matrix size check: {min_rank}..{max_rank}")
    candidate_ks = list(range(min_rank, max_rank + 1))
    denom = float(np.square(V_weighted).sum())
    if denom == 0:
        raise ValueError("Weighted NMF input is all zero.")

    rows = []
    for k in candidate_ks:
        H_by_seed, seed_errors, seed_explained = [], [], []
        for seed in cfg.rank_selection_seeds:
            model = make_nmf(n_components=k, seed=seed, cfg=cfg)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=UserWarning)
                W = model.fit_transform(V_weighted)
            H = model.components_
            err = float(np.square(V_weighted - W @ H).sum())
            H_by_seed.append(H)
            seed_errors.append(err)
            seed_explained.append(1.0 - err / denom)
        best_i = int(np.argmin(seed_errors))
        rows.append({
            "k": k,
            "n_seeds": len(cfg.rank_selection_seeds),
            "mean_balanced_explained_fraction": float(np.nanmean(seed_explained)),
            "best_balanced_explained_fraction": float(np.nanmax(seed_explained)),
            "mean_reconstruction_error": float(np.mean(seed_errors)),
            "best_reconstruction_error": float(np.min(seed_errors)),
            "stability_matched_cosine": stability_matched_cosine(H_by_seed),
            "best_seed": int(cfg.rank_selection_seeds[best_i]),
        })
    metrics = pd.DataFrame(rows)
    metrics["selection_score"] = (
        metrics["best_balanced_explained_fraction"].fillna(0.0)
        + 0.05 * metrics["stability_matched_cosine"].fillna(0.0)
        - 0.01 * metrics["k"]
    )
    if cfg.forced_cm_k is not None:
        selected_k = int(cfg.forced_cm_k)
        if selected_k not in metrics["k"].tolist():
            raise ValueError(f"forced_cm_k={selected_k} outside evaluated ranks: {metrics['k'].tolist()}")
    else:
        selected_k = int(
            metrics.sort_values(
                ["selection_score", "best_balanced_explained_fraction"],
                ascending=False,
            ).iloc[0]["k"]
        )
    metrics["selected"] = metrics["k"].eq(selected_k)
    return metrics, selected_k


def fit_best_seed_nmf(V_weighted: np.ndarray, k: int, seeds: Iterable[int], cfg: CMNMFConfig):
    best_error = math.inf
    best_W = best_H = None
    for seed in seeds:
        model = make_nmf(n_components=k, seed=seed, cfg=cfg)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            W = model.fit_transform(V_weighted)
        H = model.components_
        err = float(np.square(V_weighted - W @ H).sum())
        if err < best_error:
            best_error, best_W, best_H = err, W, H
    if best_W is None or best_H is None:
        raise RuntimeError("Final NMF failed for all seeds.")
    return best_W, best_H


def refit_activity_by_nnls(H: np.ndarray, V_unweighted: np.ndarray) -> np.ndarray:
    design = H.T
    W = np.zeros((V_unweighted.shape[0], H.shape[0]), dtype=float)
    for i in range(V_unweighted.shape[0]):
        W[i, :] = nnls(design, V_unweighted[i, :], maxiter=design.shape[1] * 10)[0]
    return W


def classify_raw_cms(W: pd.DataFrame, sample_status: pd.DataFrame, cfg: CMNMFConfig) -> pd.DataFrame:
    status = sample_status["status"].astype(str)
    normal = W.loc[status.eq(cfg.normal_status)]
    tumor = W.loc[status.eq(cfg.tumor_status)]
    if normal.empty or tumor.empty:
        raise ValueError("Need both normal-like and tumor samples to classify CMs.")
    rows = []
    for raw_cm in W.columns:
        normal_mean = float(normal[raw_cm].mean())
        tumor_mean = float(tumor[raw_cm].mean())
        normal_active = float((normal[raw_cm] > 1e-8).mean())
        tumor_active = float((tumor[raw_cm] > 1e-8).mean())
        ratio = tumor_mean / normal_mean if normal_mean > 0 else np.inf
        if ratio >= cfg.tumor_specific_min_ratio and tumor_active >= cfg.min_active_fraction_for_specific:
            cm_class = "tCM"
        elif ratio <= cfg.normal_specific_max_ratio and normal_active >= cfg.min_active_fraction_for_specific:
            cm_class = "normalCM"
        else:
            cm_class = "sharedCM"
        rows.append({
            "raw_component": raw_cm,
            "raw_component_order": int(str(raw_cm).replace("CM", "")),
            "class": cm_class,
            "normal_mean_usage": normal_mean,
            "tumor_mean_usage": tumor_mean,
            "tumor_to_normal_mean_ratio": ratio,
            "normal_active_fraction": normal_active,
            "tumor_active_fraction": tumor_active,
        })
    return pd.DataFrame(rows).sort_values("raw_component_order", kind="stable")


def assign_canonical_cm_names(classification: pd.DataFrame) -> pd.DataFrame:
    prefix = {"sharedCM": "s", "normalCM": "n", "tCM": "t"}
    out = classification.sort_values("raw_component_order", kind="stable").copy()
    out["global_order"] = range(1, len(out) + 1)
    out["class_prefix"] = out["class"].map(prefix)
    if out["class_prefix"].isna().any():
        raise ValueError("Unknown CM class in classification table.")
    out["CM"] = out["class_prefix"] + "_CM" + out["global_order"].astype(str)
    return out


def run_balanced_joint_cm_nmf(non_epi_frequency: pd.DataFrame, sample_status: pd.DataFrame,
                              out_dir: Path, cfg: CMNMFConfig):
    V_column_minmax, minmax_params = column_minmax_normalize(non_epi_frequency)
    sample_weights = make_status_balanced_weights(
        sample_status.loc[V_column_minmax.index], cfg
    )
    V = V_column_minmax.to_numpy(dtype=float)
    V_weighted = V * np.sqrt(sample_weights.to_numpy(dtype=float))[:, None]

    metrics, selected_k = evaluate_cm_rank(V_weighted, cfg)
    best_seed = int(metrics.loc[metrics["k"].eq(selected_k), "best_seed"].iloc[0])
    _, H_raw = fit_best_seed_nmf(V_weighted, selected_k, seeds=(best_seed,), cfg=cfg)
    W_raw = refit_activity_by_nnls(H_raw, V)

    raw_cms = [f"CM{i + 1}" for i in range(selected_k)]
    W_raw_df = pd.DataFrame(W_raw, index=V_column_minmax.index, columns=raw_cms)
    H_raw_df = pd.DataFrame(H_raw, index=raw_cms, columns=V_column_minmax.columns)
    classification = assign_canonical_cm_names(
        classify_raw_cms(W_raw_df, sample_status.loc[W_raw_df.index], cfg)
    )
    raw_to_cm = dict(zip(classification["raw_component"], classification["CM"]))
    cm_order = classification["CM"].tolist()

    W_df = W_raw_df.rename(columns=raw_to_cm).loc[:, cm_order]
    H_df = H_raw_df.rename(index=raw_to_cm).loc[cm_order, :]
    loading_df = H_df.T.copy()
    loading_fraction = loading_df.div(loading_df.sum(axis=0).replace(0, np.nan), axis=1).fillna(0.0)
    activity_df = W_df.join(sample_status[["status", "non_epi_cells"]])

    out_dir.mkdir(parents=True, exist_ok=True)
    V_column_minmax.to_csv(out_dir / "non_epi_subtype_frequency_column_minmax.csv")
    minmax_params.to_csv(out_dir / "column_minmax_params.csv")
    sample_weights.to_frame().join(sample_status[["status", "non_epi_cells"]]).to_csv(
        out_dir / "group_balanced_sample_weights.csv"
    )
    metrics.to_csv(out_dir / "joint_nmf_k_selection_metrics.csv", index=False)
    write_json({
        "selected_module_k": int(selected_k),
        "best_seed": best_seed,
        "candidate_Ks": metrics["k"].astype(int).tolist(),
        "rank_selection_seeds": list(cfg.rank_selection_seeds),
    }, out_dir / "selected_module_k.json")
    classification.to_csv(out_dir / "joint_module_classification.csv", index=False)
    classification[["raw_component", "raw_component_order", "CM", "class", "class_prefix", "global_order"]].to_csv(
        out_dir / "raw_to_canonical_CM_mapping.csv", index=False
    )
    W_df.to_csv(out_dir / "W_df.csv")
    H_df.to_csv(out_dir / "H_df.csv")
    loading_df.to_csv(out_dir / "loading_df_cell_subtype_by_CM.csv")
    loading_fraction.to_csv(out_dir / "loading_df_cell_subtype_by_CM_fraction.csv")
    activity_df.to_csv(out_dir / "activity_df_sample_by_CM.csv")
    return W_df, H_df, loading_df, activity_df, metrics, classification


def main() -> None:
    t0 = time.time()
    cfg = CMNMFConfig()
    V = pd.read_csv(V_CSV, index_col=0)
    sample_status = pd.read_csv(SAMPLE_STATUS_CSV, index_col=0)
    incl = pd.read_csv(SAMPLE_INCL_CSV, index_col=0)
    # keep only samples passing keep_for_cm
    keep_samples = incl.index[incl["keep_for_cm"]].astype(str).tolist()
    keep_samples = [s for s in keep_samples if s in V.index and s in sample_status.index]
    V = V.loc[keep_samples]
    sample_status = sample_status.loc[keep_samples]
    print(f"[nmf] V: {V.shape}, samples={len(keep_samples)}, time={time.time()-t0:.1f}s", flush=True)

    n_samples_valid, n_subs_valid = V.shape
    actual_max_K = min(20, n_samples_valid, n_subs_valid)
    print(f"[nmf] n_samples={n_samples_valid}, n_subtypes_valid={n_subs_valid}, actual_max_K={actual_max_K}", flush=True)

    W_df, H_df, loading_df, activity_df, metrics, classification = run_balanced_joint_cm_nmf(
        V, sample_status, TAB, cfg
    )
    best_k = int(metrics.loc[metrics["selected"], "k"].iloc[0])
    best_seed = int(metrics.loc[metrics["selected"], "best_seed"].iloc[0])
    n_shared = int((classification["class"] == "sharedCM").sum())
    n_t = int((classification["class"] == "tCM").sum())
    n_n = int((classification["class"] == "normalCM").sum())
    print(f"[nmf] SELECTED K={best_k}, best_seed={best_seed}, classes: shared={n_shared}, t={n_t}, normal={n_n}", flush=True)
    print(classification[["raw_component", "CM", "class"]].to_string(), flush=True)
    print(f"[nmf] done. total={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()