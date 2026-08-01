#!/usr/bin/env python3
"""Shared fixed helpers for BRCA Block 03 CM-lineage analysis."""

from __future__ import annotations

import json
import math
import random
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment, nnls
from scipy.stats import pearsonr, spearmanr
from sklearn.decomposition import NMF
from sklearn.metrics.pairwise import cosine_similarity
from statsmodels.stats.multitest import multipletests


SEED = 42
random.seed(SEED)
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/brca_cm-epi_codex")
WORKFLOW = ROOT / "epi-cm-core-workflow"
BLOCK = "03-epi-cm-discovery"
ANALYSIS_TABLE_ROOT = WORKFLOW / "tables" / BLOCK / "01-cm-lineage-analysis"
PLOT_TABLE_ROOT = WORKFLOW / "tables" / BLOCK / "02-cm-lineage-final-plotting"
FIGURE_ROOT = WORKFLOW / "figures" / BLOCK / "02-cm-lineage-final-plotting"

PREP_DIR = ANALYSIS_TABLE_ROOT / "01_prepare_inputs_and_frequency_tables"
NMF_DIR = ANALYSIS_TABLE_ROOT / "02_balanced_joint_nmf"
NODE_DIR = ANALYSIS_TABLE_ROOT / "03_cm_classification_nodes_and_edges"
SPEARMAN_DIR = ANALYSIS_TABLE_ROOT / "04_epi_cm_association_spearman"

CM_K_RANGE = (2, 20)
RANK_SELECTION_SEEDS = (0, 1, 2, 3, 4)
NMF_MAX_ITER = 3000
NMF_TOL = 1e-5
NMF_ALPHA_W = 0.0
NMF_ALPHA_H = 1e-3
NMF_L1_RATIO = 0.1
EDGE_R_THRESHOLD = 0.25
TOP_N_SUBTYPES = 20
PLOT_TOP_N_SUBTYPES = 12
TOP_N_NODES = 10


def write_json(value: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def column_minmax_normalize(freq: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    column_min = freq.min(axis=0)
    column_range = (freq.max(axis=0) - column_min).replace(0, np.nan)
    normalized = ((freq - column_min) / column_range).fillna(0.0).clip(lower=0.0)
    params = pd.DataFrame({"min": column_min, "range": column_range.fillna(0.0)})
    params.index.name = "cell_subtype"
    if (normalized.var(axis=0) > 0).sum() == 0:
        raise ValueError("Column-minmax matrix has no non-zero-variance subtype columns.")
    return normalized.astype(float), params


def make_nmf(n_components: int, seed: int) -> NMF:
    return NMF(
        n_components=n_components,
        init="nndsvda",
        solver="cd",
        beta_loss="frobenius",
        random_state=seed,
        max_iter=NMF_MAX_ITER,
        tol=NMF_TOL,
        alpha_W=NMF_ALPHA_W,
        alpha_H=NMF_ALPHA_H,
        l1_ratio=NMF_L1_RATIO,
    )


def row_l1_normalize(matrix: np.ndarray) -> np.ndarray:
    denom = matrix.sum(axis=1, keepdims=True)
    return np.divide(matrix, denom, out=np.zeros_like(matrix, dtype=float), where=denom != 0)


def matched_cosine_mean(h_a: np.ndarray, h_b: np.ndarray) -> float:
    similarity = cosine_similarity(row_l1_normalize(h_a), row_l1_normalize(h_b))
    row_ind, col_ind = linear_sum_assignment(-similarity)
    return float(similarity[row_ind, col_ind].mean())


def stability_matched_cosine(h_by_seed: list[np.ndarray]) -> float:
    values = [
        matched_cosine_mean(h_by_seed[i], h_by_seed[j])
        for i in range(len(h_by_seed))
        for j in range(i + 1, len(h_by_seed))
    ]
    return float(np.mean(values)) if values else float("nan")


def evaluate_cm_rank(v_weighted: np.ndarray) -> tuple[pd.DataFrame, int]:
    max_rank = min(CM_K_RANGE[1], v_weighted.shape[0], v_weighted.shape[1])
    min_rank = max(2, CM_K_RANGE[0])
    if max_rank < min_rank:
        raise ValueError(f"Invalid CM K range: {min_rank}..{max_rank}")
    candidate_ks = list(range(min_rank, max_rank + 1))
    denom = float(np.square(v_weighted).sum())
    if denom == 0:
        raise ValueError("Weighted NMF input is all zero.")
    rows: list[dict[str, object]] = []
    for k in candidate_ks:
        h_by_seed: list[np.ndarray] = []
        seed_errors: list[float] = []
        seed_explained: list[float] = []
        converged: list[bool] = []
        iterations: list[int] = []
        for seed in RANK_SELECTION_SEEDS:
            model = make_nmf(k, seed)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=UserWarning)
                w = model.fit_transform(v_weighted)
            h = model.components_
            error = float(np.square(v_weighted - w @ h).sum())
            h_by_seed.append(h)
            seed_errors.append(error)
            seed_explained.append(1.0 - error / denom)
            iterations.append(int(model.n_iter_))
            converged.append(int(model.n_iter_) < NMF_MAX_ITER)
        best_i = int(np.argmin(seed_errors))
        rows.append(
            {
                "k": k,
                "n_seeds": len(RANK_SELECTION_SEEDS),
                "mean_balanced_explained_fraction": float(np.nanmean(seed_explained)),
                "best_balanced_explained_fraction": float(np.nanmax(seed_explained)),
                "mean_reconstruction_error": float(np.mean(seed_errors)),
                "best_reconstruction_error": float(np.min(seed_errors)),
                "stability_matched_cosine": stability_matched_cosine(h_by_seed),
                "best_seed": int(RANK_SELECTION_SEEDS[best_i]),
                "all_seeds_converged": bool(all(converged)),
                "max_n_iter_observed": int(max(iterations)),
            }
        )
    metrics = pd.DataFrame(rows)
    metrics["selection_score"] = (
        metrics["best_balanced_explained_fraction"].fillna(0.0)
        + 0.05 * metrics["stability_matched_cosine"].fillna(0.0)
        - 0.01 * metrics["k"]
    )
    selected_k = int(
        metrics.sort_values(
            ["selection_score", "best_balanced_explained_fraction"],
            ascending=False,
        ).iloc[0]["k"]
    )
    metrics["selected"] = metrics["k"].eq(selected_k)
    return metrics, selected_k


def fit_nmf(v_weighted: np.ndarray, k: int, seed: int) -> tuple[np.ndarray, np.ndarray, float, int]:
    model = make_nmf(k, seed)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        w = model.fit_transform(v_weighted)
    h = model.components_
    error = float(np.square(v_weighted - w @ h).sum())
    return w, h, error, int(model.n_iter_)


def refit_activity_by_nnls(h: np.ndarray, v_unweighted: np.ndarray) -> np.ndarray:
    design = h.T
    w = np.zeros((v_unweighted.shape[0], h.shape[0]), dtype=float)
    for index in range(v_unweighted.shape[0]):
        w[index, :] = nnls(design, v_unweighted[index, :], maxiter=design.shape[1] * 10)[0]
    return w


def transform_cm_columns(frame: pd.DataFrame, method: str) -> pd.DataFrame:
    x = frame.astype(float).copy()
    if method == "raw":
        return x
    if method == "standard_scale_col":
        minimum = x.min(axis=0)
        value_range = (x.max(axis=0) - minimum).replace(0, np.nan)
        return x.sub(minimum, axis=1).div(value_range, axis=1).fillna(0.0).clip(0, 1)
    if method == "zscore":
        return x.sub(x.mean(axis=0), axis=1).div(
            x.std(axis=0, ddof=0).replace(0, np.nan), axis=1
        ).fillna(0.0)
    if method == "robust":
        median = x.median(axis=0)
        iqr = (x.quantile(0.75) - x.quantile(0.25)).replace(0, np.nan)
        return x.sub(median, axis=1).div(iqr, axis=1).fillna(0.0)
    raise ValueError(f"Unknown display transform: {method}")


def bh_fdr(pvalues: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=pvalues.index, dtype=float)
    valid = pvalues.notna() & np.isfinite(pvalues.astype(float))
    if valid.any():
        result.loc[valid] = multipletests(pvalues.loc[valid].astype(float), method="fdr_bh")[1]
    return result


def safe_corr(x: pd.Series, y: pd.Series, method: str) -> tuple[float, float, int]:
    valid = x.notna() & y.notna()
    n = int(valid.sum())
    if n < 3 or x.loc[valid].nunique() < 2 or y.loc[valid].nunique() < 2:
        return float("nan"), float("nan"), n
    if method == "spearman":
        stat, p = spearmanr(x.loc[valid], y.loc[valid])
    elif method == "pearson":
        stat, p = pearsonr(x.loc[valid], y.loc[valid])
    else:
        raise ValueError(f"Unsupported correlation method: {method}")
    return float(stat), float(p), n


def cm_sort_key(value: str) -> tuple[int, str]:
    text = str(value)
    try:
        return int(text.split("CM")[-1]), text
    except ValueError:
        return 10**9, text

