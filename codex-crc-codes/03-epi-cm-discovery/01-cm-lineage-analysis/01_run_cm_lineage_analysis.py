#!/usr/bin/env python3
"""Canonical CRC Epi-CM discovery analysis (tables only; no final plots).

This script implements Block 03 / submodule 01 from the project workflow skill:
frequency-table preparation, status-balanced joint NMF, canonical CM naming,
node/edge derivation, and exhaustive Spearman plus Pearson Epi-CM association.
"""

from __future__ import annotations

import json
import platform
import random
import re
import sys
import warnings
from dataclasses import asdict, dataclass
from itertools import combinations
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy
import sklearn
from scipy.optimize import linear_sum_assignment, nnls
from scipy.stats import mannwhitneyu, pearsonr, spearmanr
from sklearn.decomposition import NMF
from sklearn.metrics.pairwise import cosine_similarity
from statsmodels.stats.multitest import multipletests


SEED = 42
random.seed(SEED)
np.random.seed(SEED)

WORKFLOW_ROOT = Path(__file__).resolve().parents[3]
INPUT_H5AD = (
    WORKFLOW_ROOT
    / "h5ad/02-cell_subtype_integration_clustering/06-project-subtypes-to-full-adata/adata_anno_cellsubtype.h5ad"
)
TABLE_ROOT = WORKFLOW_ROOT / "tables/03-epi-cm-discovery"
ANALYSIS_ROOT = TABLE_ROOT / "01-cm-lineage-analysis"
PREP_DIR = ANALYSIS_ROOT / "01_prepare_inputs_and_frequency_tables"
NMF_DIR = ANALYSIS_ROOT / "02_balanced_joint_nmf"
NODE_DIR = ANALYSIS_ROOT / "03_cm_classification_nodes_and_edges"
SPEARMAN_DIR = ANALYSIS_ROOT / "04_epi_cm_association_spearman"
PEARSON_DIR = ANALYSIS_ROOT / "05_epi_cm_association_pearson"


@dataclass(frozen=True)
class Config:
    sample_col: str = "sample"
    series_col: str = "series"
    status_col: str = "status"
    celltype_col: str = "cell_type"
    subtype_col: str = "cell_subtype"
    epithelial_cell_type: str = "Epithelial Cells"
    normal_status: str = "normal"
    tumor_status: str = "tumor"
    min_non_epi_cells_per_sample: int = 50
    min_epi_cells_per_sample: int = 1
    forced_cm_k: int | None = None
    cm_k_min: int = 2
    cm_k_max: int = 20
    rank_selection_seeds: tuple[int, ...] = (0, 1, 2, 3, 4)
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
    top_n_subtypes: int = 20
    plot_top_n_subtypes: int = 12
    top_n_nodes: int = 10
    edge_r_threshold: float = 0.25


CFG = Config()


def ensure_dirs() -> None:
    for path in (PREP_DIR, NMF_DIR, NODE_DIR, SPEARMAN_DIR, PEARSON_DIR):
        path.mkdir(parents=True, exist_ok=True)


def bh_fdr(values: np.ndarray | pd.Series) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    q = np.full(p.shape, np.nan, dtype=float)
    valid = np.isfinite(p)
    if valid.any():
        q[valid] = multipletests(p[valid], method="fdr_bh")[1]
    return q


def write_json(payload: dict, path: Path) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")


def prefix_cell_type(cell_subtype: pd.Series) -> pd.Series:
    mapping = {
        "Epi": "Epithelial Cells",
        "T": "T Cells",
        "B": "B Cells",
        "Mye": "Myeloid Cells",
        "Endo": "Endothelial Cells",
        "Stromal": "Stromal Cells",
        "Mast": "Mast Cells",
        "Plasma": "Plasma Cells",
        "Cycling": "Cycling Immune Cells",
        "Schwann": "Schwann Cells",
    }
    prefix = cell_subtype.astype(str).str.split("_", n=1).str[0]
    derived = prefix.map(mapping)
    if derived.isna().any():
        missing = sorted(prefix.loc[derived.isna()].unique().tolist())
        raise ValueError(f"No broad-cell mapping for subtype prefixes: {missing}")
    return derived.astype(str)


def load_obs() -> pd.DataFrame:
    if not INPUT_H5AD.exists():
        raise FileNotFoundError(INPUT_H5AD)
    adata = ad.read_h5ad(INPUT_H5AD, backed="r")
    required = [CFG.sample_col, CFG.series_col, CFG.status_col, CFG.celltype_col, CFG.subtype_col]
    missing = [c for c in required if c not in adata.obs.columns]
    if missing:
        adata.file.close()
        raise ValueError(f"Projected atlas is missing obs columns: {missing}")
    obs = adata.obs[required].copy()
    adata.file.close()
    obs.index = obs.index.astype(str)
    if not obs.index.is_unique:
        raise ValueError("Projected atlas obs_names are not unique")
    if obs[required].isna().any().any():
        counts = obs[required].isna().sum()
        raise ValueError(f"Missing required obs values: {counts[counts > 0].to_dict()}")
    for col in required:
        obs[col] = obs[col].astype(str)
    observed_status = set(obs[CFG.status_col].unique())
    expected_status = {CFG.normal_status, CFG.tumor_status}
    if observed_status != expected_status:
        raise ValueError(f"Expected status={sorted(expected_status)}, observed={sorted(observed_status)}")
    for col in (CFG.status_col, CFG.series_col):
        conflicts = obs[[CFG.sample_col, col]].drop_duplicates()[CFG.sample_col].duplicated().sum()
        if conflicts:
            raise ValueError(f"{conflicts} samples have conflicting {col} values")
    obs["cell_type_from_subtype_prefix"] = prefix_cell_type(obs[CFG.subtype_col])
    conflicts = obs.loc[
        obs[CFG.celltype_col].ne(obs["cell_type_from_subtype_prefix"]),
        [CFG.sample_col, CFG.celltype_col, CFG.subtype_col, "cell_type_from_subtype_prefix"],
    ]
    conflicts.to_csv(PREP_DIR / "cell_type_subtype_prefix_conflicts.csv", index=True, index_label="cell_id")
    # FIXED: subtype prefix is authoritative for CM matrix compartment assignment.
    obs[CFG.celltype_col] = obs["cell_type_from_subtype_prefix"]
    return obs


def prepare_frequency_tables(obs: pd.DataFrame):
    metadata = obs[[CFG.sample_col, CFG.status_col, CFG.series_col]].drop_duplicates()
    sample_status = metadata.set_index(CFG.sample_col).sort_index()
    sample_status.index.name = "sample"

    epi_mask = obs[CFG.celltype_col].eq(CFG.epithelial_cell_type)
    epi = obs.loc[epi_mask]
    non_epi = obs.loc[~epi_mask]
    if epi.empty or non_epi.empty:
        raise ValueError(f"Empty compartment: epi={len(epi)}, non_epi={len(non_epi)}")

    non_epi_counts_all = pd.crosstab(non_epi[CFG.sample_col], non_epi[CFG.subtype_col]).sort_index().astype(int)
    epi_counts_all = pd.crosstab(epi[CFG.sample_col], epi[CFG.subtype_col]).sort_index().astype(int)

    compartment = pd.DataFrame(index=sample_status.index)
    compartment["non_epi_cells"] = non_epi_counts_all.sum(axis=1).reindex(compartment.index).fillna(0).astype(int)
    compartment["epi_cells"] = epi_counts_all.sum(axis=1).reindex(compartment.index).fillna(0).astype(int)
    compartment.index.name = "sample"

    inclusion = compartment.copy()
    inclusion["has_status"] = True
    inclusion["has_non_epi"] = inclusion["non_epi_cells"].ge(CFG.min_non_epi_cells_per_sample)
    inclusion["has_epi_for_association"] = inclusion["epi_cells"].ge(CFG.min_epi_cells_per_sample)
    inclusion["keep_for_cm"] = inclusion["has_status"] & inclusion["has_non_epi"]
    inclusion["keep_for_epi_cm"] = inclusion["keep_for_cm"] & inclusion["has_epi_for_association"]

    def reason(row: pd.Series) -> str:
        reasons = []
        if not row["has_status"]:
            reasons.append("missing_status")
        if not row["has_non_epi"]:
            reasons.append("non_epi_cells_below_threshold")
        if not row["has_epi_for_association"]:
            reasons.append("epi_cells_below_association_threshold")
        return ";".join(reasons) if reasons else "kept"

    inclusion["exclusion_reason"] = inclusion.apply(reason, axis=1)
    keep_cm = inclusion.index[inclusion["keep_for_cm"]]
    keep_epi = inclusion.index[inclusion["keep_for_epi_cm"]]
    if len(keep_cm) < 2:
        raise ValueError("Fewer than two samples pass CM eligibility")

    non_epi_counts = non_epi_counts_all.reindex(keep_cm).fillna(0).astype(int)
    non_epi_freq = non_epi_counts.div(non_epi_counts.sum(axis=1), axis=0).fillna(0.0)
    epi_counts = epi_counts_all.reindex(keep_epi).fillna(0).astype(int)
    epi_freq = epi_counts.div(epi_counts.sum(axis=1), axis=0).fillna(0.0)
    sample_status_cm = sample_status.loc[keep_cm].copy()
    sample_status_cm["non_epi_cells"] = compartment.loc[keep_cm, "non_epi_cells"]
    sample_status_cm["epi_cells"] = compartment.loc[keep_cm, "epi_cells"]

    non_epi_counts.to_csv(PREP_DIR / "non_epi_subtype_counts.csv")
    non_epi_freq.to_csv(PREP_DIR / "non_epi_subtype_frequency.csv")
    epi_counts.to_csv(PREP_DIR / "epi_subtype_counts.csv")
    epi_freq.to_csv(PREP_DIR / "epi_subtype_frequency.csv")
    sample_status_cm.to_csv(PREP_DIR / "sample_status.csv")
    compartment.to_csv(PREP_DIR / "sample_compartment_cell_counts.csv")
    inclusion.to_csv(PREP_DIR / "sample_inclusion_exclusion.csv")
    return non_epi_freq, epi_freq, sample_status_cm, inclusion


def column_minmax_normalize(freq: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    column_min = freq.min(axis=0)
    column_range = (freq.max(axis=0) - column_min).replace(0, np.nan)
    normalized = ((freq - column_min) / column_range).fillna(0.0).clip(lower=0.0)
    params = pd.DataFrame({"min": column_min, "range": column_range.fillna(0.0)})
    params.index.name = "cell_subtype"
    valid = normalized.var(axis=0).gt(0)
    if valid.sum() == 0:
        raise ValueError("Column-minmax matrix has no non-zero-variance subtype columns")
    return normalized.astype(float), params


def balanced_weights(sample_status: pd.DataFrame) -> pd.Series:
    status = sample_status[CFG.status_col].astype(str)
    normal_mask = status.eq(CFG.normal_status)
    tumor_mask = status.eq(CFG.tumor_status)
    if normal_mask.sum() == 0 or tumor_mask.sum() == 0:
        raise ValueError("Both normal and tumor samples are required")
    weights = pd.Series(0.0, index=sample_status.index, name="joint_nmf_row_weight")
    weights.loc[normal_mask] = CFG.normal_group_total_weight / normal_mask.sum()
    weights.loc[tumor_mask] = CFG.tumor_group_total_weight / tumor_mask.sum()
    return weights / weights.mean()


def make_nmf(k: int, seed: int) -> NMF:
    return NMF(
        n_components=k,
        init="nndsvda",
        solver="cd",
        beta_loss="frobenius",
        random_state=seed,
        max_iter=CFG.nmf_max_iter,
        tol=CFG.nmf_tol,
        alpha_W=CFG.nmf_alpha_w,
        alpha_H=CFG.nmf_alpha_h,
        l1_ratio=CFG.nmf_l1_ratio,
    )


def row_l1(matrix: np.ndarray) -> np.ndarray:
    denom = matrix.sum(axis=1, keepdims=True)
    return np.divide(matrix, denom, out=np.zeros_like(matrix, dtype=float), where=denom != 0)


def matched_cosine(H_a: np.ndarray, H_b: np.ndarray) -> float:
    sim = cosine_similarity(row_l1(H_a), row_l1(H_b))
    rows, cols = linear_sum_assignment(-sim)
    return float(sim[rows, cols].mean())


def evaluate_rank(V_weighted: np.ndarray) -> tuple[pd.DataFrame, int]:
    max_rank = min(CFG.cm_k_max, V_weighted.shape[0], V_weighted.shape[1])
    min_rank = max(2, CFG.cm_k_min)
    if max_rank < min_rank:
        raise ValueError(f"Invalid K range: {min_rank}..{max_rank}")
    denom = float(np.square(V_weighted).sum())
    if denom == 0:
        raise ValueError("Weighted NMF input is all zero")
    rows = []
    for k in range(min_rank, max_rank + 1):
        Hs, errors, explained = [], [], []
        for seed in CFG.rank_selection_seeds:
            model = make_nmf(k, seed)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                W = model.fit_transform(V_weighted)
            H = model.components_
            err = float(np.square(V_weighted - W @ H).sum())
            Hs.append(H)
            errors.append(err)
            explained.append(1.0 - err / denom)
        stabilities = [matched_cosine(Hs[i], Hs[j]) for i in range(len(Hs)) for j in range(i + 1, len(Hs))]
        best_i = int(np.argmin(errors))
        rows.append(
            {
                "k": k,
                "n_seeds": len(CFG.rank_selection_seeds),
                "mean_balanced_explained_fraction": float(np.mean(explained)),
                "best_balanced_explained_fraction": float(np.max(explained)),
                "mean_reconstruction_error": float(np.mean(errors)),
                "best_reconstruction_error": float(np.min(errors)),
                "stability_matched_cosine": float(np.mean(stabilities)),
                "best_seed": int(CFG.rank_selection_seeds[best_i]),
            }
        )
    metrics = pd.DataFrame(rows)
    metrics["selection_score"] = (
        metrics["best_balanced_explained_fraction"]
        + 0.05 * metrics["stability_matched_cosine"]
        - 0.01 * metrics["k"]
    )
    if CFG.forced_cm_k is None:
        selected_k = int(
            metrics.sort_values(
                ["selection_score", "best_balanced_explained_fraction"], ascending=False
            ).iloc[0]["k"]
        )
    else:
        selected_k = int(CFG.forced_cm_k)
        if selected_k not in metrics["k"].tolist():
            raise ValueError(f"Forced K {selected_k} was not evaluated")
    metrics["selected"] = metrics["k"].eq(selected_k)
    return metrics, selected_k


def fit_selected(V_weighted: np.ndarray, k: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    model = make_nmf(k, seed)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        W = model.fit_transform(V_weighted)
    return W, model.components_


def refit_nnls(H: np.ndarray, V: np.ndarray) -> np.ndarray:
    design = H.T
    W = np.zeros((V.shape[0], H.shape[0]), dtype=float)
    for i in range(V.shape[0]):
        W[i, :] = nnls(design, V[i, :], maxiter=design.shape[1] * 10)[0]
    return W


def classify_cms(W_raw: pd.DataFrame, sample_status: pd.DataFrame) -> pd.DataFrame:
    status = sample_status[CFG.status_col].astype(str)
    normal_W = W_raw.loc[status.eq(CFG.normal_status)]
    tumor_W = W_raw.loc[status.eq(CFG.tumor_status)]
    rows = []
    for raw_cm in W_raw.columns:
        normal_values = normal_W[raw_cm].astype(float)
        tumor_values = tumor_W[raw_cm].astype(float)
        normal_mean = float(normal_values.mean())
        tumor_mean = float(tumor_values.mean())
        ratio = tumor_mean / normal_mean if normal_mean > 0 else np.inf
        normal_active = float(normal_values.gt(1e-8).mean())
        tumor_active = float(tumor_values.gt(1e-8).mean())
        if ratio >= CFG.tumor_specific_min_ratio and tumor_active >= CFG.min_active_fraction_for_specific:
            cm_class = "tCM"
        elif ratio <= CFG.normal_specific_max_ratio and normal_active >= CFG.min_active_fraction_for_specific:
            cm_class = "normalCM"
        else:
            cm_class = "sharedCM"
        _, p = mannwhitneyu(tumor_values, normal_values, alternative="two-sided")
        rows.append(
            {
                "raw_component": raw_cm,
                "raw_component_order": int(re.sub(r"\D", "", raw_cm)),
                "class": cm_class,
                "normal_mean": normal_mean,
                "tumor_mean": tumor_mean,
                "delta": tumor_mean - normal_mean,
                "tumor_to_normal_mean_ratio": ratio,
                "normal_active_fraction": normal_active,
                "tumor_active_fraction": tumor_active,
                "p": float(p),
            }
        )
    out = pd.DataFrame(rows).sort_values("raw_component_order", kind="stable")
    out["q"] = bh_fdr(out["p"])
    prefix = {"sharedCM": "s", "normalCM": "n", "tCM": "t"}
    out["global_order"] = np.arange(1, len(out) + 1)
    out["class_prefix"] = out["class"].map(prefix)
    out["CM"] = out["class_prefix"] + "_CM" + out["global_order"].astype(str)
    return out


def transform_columns(df: pd.DataFrame, method: str) -> pd.DataFrame:
    x = df.astype(float).copy()
    if method == "raw":
        return x
    if method == "standard_scale_col":
        mn = x.min(axis=0)
        rng = (x.max(axis=0) - mn).replace(0, np.nan)
        return x.sub(mn, axis=1).div(rng, axis=1).fillna(0.0).clip(0, 1)
    if method == "zscore":
        return x.sub(x.mean(axis=0), axis=1).div(x.std(axis=0, ddof=0).replace(0, np.nan), axis=1).fillna(0.0)
    if method == "robust":
        med = x.median(axis=0)
        iqr = (x.quantile(0.75) - x.quantile(0.25)).replace(0, np.nan)
        return x.sub(med, axis=1).div(iqr, axis=1).fillna(0.0)
    raise ValueError(method)


def fit_balanced_joint_nmf(non_epi_freq: pd.DataFrame, sample_status: pd.DataFrame):
    V_df, minmax = column_minmax_normalize(non_epi_freq)
    weights = balanced_weights(sample_status.loc[V_df.index])
    V = V_df.to_numpy(float)
    V_weighted = V * np.sqrt(weights.to_numpy(float))[:, None]
    metrics, selected_k = evaluate_rank(V_weighted)
    best_seed = int(metrics.loc[metrics["selected"], "best_seed"].iloc[0])
    _, H_raw = fit_selected(V_weighted, selected_k, best_seed)
    W_raw = refit_nnls(H_raw, V)

    raw_cms = [f"CM{i + 1}" for i in range(selected_k)]
    W_raw_df = pd.DataFrame(W_raw, index=V_df.index, columns=raw_cms)
    H_raw_df = pd.DataFrame(H_raw, index=raw_cms, columns=V_df.columns)
    classification = classify_cms(W_raw_df, sample_status.loc[W_raw_df.index])
    raw_to_cm = dict(zip(classification["raw_component"], classification["CM"]))
    cm_order = classification["CM"].tolist()
    W_df = W_raw_df.rename(columns=raw_to_cm).loc[:, cm_order]
    H_df = H_raw_df.rename(index=raw_to_cm).loc[cm_order]
    loading = H_df.T
    loading_fraction = loading.div(loading.sum(axis=0).replace(0, np.nan), axis=1).fillna(0.0)

    V_df.to_csv(NMF_DIR / "non_epi_subtype_frequency_column_minmax.csv")
    minmax.to_csv(NMF_DIR / "column_minmax_params.csv")
    weights.to_frame().join(sample_status[[CFG.status_col, CFG.series_col, "non_epi_cells"]]).to_csv(
        NMF_DIR / "group_balanced_sample_weights.csv"
    )
    metrics.to_csv(NMF_DIR / "joint_nmf_k_selection_metrics.csv", index=False)
    write_json(
        {
            "selected_module_k": selected_k,
            "best_seed": best_seed,
            "candidate_Ks": metrics["k"].astype(int).tolist(),
            "rank_selection_seeds": list(CFG.rank_selection_seeds),
            "forced": CFG.forced_cm_k is not None,
        },
        NMF_DIR / "selected_module_k.json",
    )
    classification.to_csv(NODE_DIR / "joint_module_classification.csv", index=False)
    classification[
        ["raw_component", "raw_component_order", "CM", "class", "class_prefix", "global_order"]
    ].to_csv(NODE_DIR / "raw_to_canonical_CM_mapping.csv", index=False)
    W_df.to_csv(NMF_DIR / "W_df.csv")
    H_df.to_csv(NMF_DIR / "H_df.csv")
    W_df.T.to_csv(NMF_DIR / "activity_df_CM_by_sample.csv")
    loading.to_csv(NMF_DIR / "loading_df_cell_subtype_by_CM.csv")
    loading_fraction.to_csv(NMF_DIR / "loading_df_cell_subtype_by_CM_fraction.csv")
    loading.to_csv(NMF_DIR / "balanced_joint_cm_subtype_loadings_raw_from_H_df.csv")
    loading_fraction.to_csv(NMF_DIR / "balanced_joint_cm_subtype_loadings_fraction_from_H_df.csv")

    activity = W_df.join(sample_status[[CFG.status_col, CFG.series_col, "non_epi_cells", "epi_cells"]])
    activity.to_csv(NMF_DIR / "activity_df_sample_by_CM.csv")
    for method in ("raw", "zscore", "robust", "standard_scale_col"):
        transform_columns(W_df, method).to_csv(NMF_DIR / f"w_df_activity_sample_by_CM_{method}.csv")
        # Historical h_df files are CM x subtype; transform in subtype x CM orientation first.
        transform_columns(loading, method).T.to_csv(NMF_DIR / f"h_df_loading_cell_subtype_by_CM_{method}.csv")
    w_minmax = transform_columns(W_df, "standard_scale_col")
    w_minmax.to_csv(NMF_DIR / "w_df_activity_sample_by_CM_standard_scale_col_with_series_plot_matrix.csv")
    sample_status.loc[W_df.index, [CFG.series_col, CFG.status_col]].rename(
        columns={CFG.series_col: "Series", CFG.status_col: "Status"}
    ).to_csv(NMF_DIR / "w_df_activity_sample_by_CM_standard_scale_col_with_series_row_annotations.csv")

    summary_rows = []
    for cm in cm_order:
        normal_values = W_df.loc[sample_status[CFG.status_col].eq(CFG.normal_status), cm]
        tumor_values = W_df.loc[sample_status[CFG.status_col].eq(CFG.tumor_status), cm]
        _, p = mannwhitneyu(tumor_values, normal_values, alternative="two-sided")
        summary_rows.append(
            {
                "CM": cm,
                "normal_mean": float(normal_values.mean()),
                "normal_sd": float(normal_values.std(ddof=1)),
                "tumor_mean": float(tumor_values.mean()),
                "tumor_sd": float(tumor_values.std(ddof=1)),
                "delta_tumor_minus_normal": float(tumor_values.mean() - normal_values.mean()),
                "mwu_p": float(p),
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary["mwu_q"] = bh_fdr(summary["mwu_p"])
    summary.to_csv(NODE_DIR / "joint_CM_activity_tumor_vs_normal_mean_sd_summary.csv", index=False)

    orientation = pd.DataFrame(
        [
            ("W_df.csv", "sample x canonical CM"),
            ("H_df.csv", "canonical CM x non-epithelial cell_subtype"),
            ("activity_df_sample_by_CM.csv", "sample x canonical CM plus metadata"),
            ("loading_df_cell_subtype_by_CM.csv", "cell_subtype x canonical CM"),
            ("h_df_loading_cell_subtype_by_CM_*.csv", "canonical CM x cell_subtype"),
        ],
        columns=["table", "orientation"],
    )
    orientation.to_csv(NMF_DIR / "matrix_orientation.csv", index=False)
    return V_df, W_df, H_df, loading, classification


def lineage_from_subtype(name: str) -> str:
    return str(name).split("_", 1)[0]


def corr_and_q_matrix(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    cols = list(df.columns)
    corr = pd.DataFrame(np.eye(len(cols)), index=cols, columns=cols, dtype=float)
    pmat = pd.DataFrame(np.zeros((len(cols), len(cols))), index=cols, columns=cols, dtype=float)
    pair_p, pairs = [], []
    for a, b in combinations(cols, 2):
        x, y = df[a].astype(float), df[b].astype(float)
        ok = x.notna() & y.notna()
        if ok.sum() < 3 or x[ok].nunique() < 2 or y[ok].nunique() < 2:
            r, p = np.nan, np.nan
        else:
            r, p = pearsonr(x[ok], y[ok])
        corr.loc[a, b] = corr.loc[b, a] = r
        pmat.loc[a, b] = pmat.loc[b, a] = p
        pair_p.append(p)
        pairs.append((a, b))
    q_values = bh_fdr(np.array(pair_p))
    qmat = pd.DataFrame(np.zeros((len(cols), len(cols))), index=cols, columns=cols, dtype=float)
    for (a, b), q in zip(pairs, q_values):
        qmat.loc[a, b] = qmat.loc[b, a] = q
    return corr, qmat


def build_nodes_and_edges(
    V_df: pd.DataFrame,
    H_df: pd.DataFrame,
    sample_status: pd.DataFrame,
) -> None:
    rows_all = []
    for cm in H_df.index:
        ranked = H_df.loc[cm].sort_values(ascending=False, kind="stable")
        for rank, (subtype, loading) in enumerate(ranked.items(), start=1):
            rows_all.append(
                {
                    "CM": cm,
                    "cell_subtype": subtype,
                    "loading": float(loading),
                    "rank": rank,
                    "cell_lineage": lineage_from_subtype(subtype),
                }
            )
    all_nodes = pd.DataFrame(rows_all)
    top20 = all_nodes.loc[all_nodes["rank"].le(CFG.top_n_subtypes)].copy()
    top10 = all_nodes.loc[all_nodes["rank"].le(CFG.top_n_nodes)].copy()
    all_nodes.to_csv(NODE_DIR / "joint_cm_cell_subtype_nodes_all_from_H_df.csv", index=False)
    top20.to_csv(NODE_DIR / "joint_cm_cell_subtype_nodes_top20_from_H_df.csv", index=False)
    top10.to_csv(NODE_DIR / "joint_cm_cell_subtype_nodes_top10_from_H_df.csv", index=False)

    status = sample_status.loc[V_df.index, CFG.status_col].astype(str)
    retained_rows, edge_rows = [], []
    for cm in H_df.index:
        candidates = top10.loc[top10["CM"].eq(cm)].sort_values("rank")["cell_subtype"].tolist()
        normal_corr = V_df.loc[status.eq(CFG.normal_status), candidates].corr(method="pearson")
        tumor_corr = V_df.loc[status.eq(CFG.tumor_status), candidates].corr(method="pearson")
        retained = set()
        for a, b in combinations(candidates, 2):
            for matrix in (normal_corr, tumor_corr):
                r = matrix.loc[a, b]
                if np.isfinite(r) and r >= CFG.edge_r_threshold:
                    retained.update((a, b))
        ordered_retained = [node for node in candidates if node in retained]
        for rank, node in enumerate(ordered_retained, start=1):
            retained_rows.append({"CM": cm, "reference_node_rank": rank, "node": node})
        for context, matrix in ((CFG.normal_status, normal_corr), (CFG.tumor_status, tumor_corr)):
            for a, b in combinations(ordered_retained, 2):
                r = float(matrix.loc[a, b]) if np.isfinite(matrix.loc[a, b]) else np.nan
                edge_rows.append(
                    {
                        "context": context,
                        "CM": cm,
                        "node_a": a,
                        "node_b": b,
                        "pearson_r": r,
                        "edge_pass_r_ge_0.25": bool(np.isfinite(r) and r >= CFG.edge_r_threshold),
                    }
                )
    membership = pd.DataFrame(retained_rows, columns=["CM", "reference_node_rank", "node"])
    missing_cms = sorted(set(H_df.index) - set(membership["CM"]))
    if missing_cms:
        raise ValueError(f"No edge-threshold reference nodes retained for CMs: {missing_cms}")
    edges = pd.DataFrame(edge_rows)
    membership.to_csv(NODE_DIR / "balanced_joint_cm_reference_node_sets_after_edge_threshold.csv", index=False)
    edges.to_csv(NODE_DIR / "status_specific_nodeplot_edges.csv", index=False)

    loading_lookup = H_df.stack().rename("loading")
    node_records = []
    for row in membership.itertuples(index=False):
        rank = int(top10.loc[(top10["CM"].eq(row.CM)) & (top10["cell_subtype"].eq(row.node)), "rank"].iloc[0])
        for context in (CFG.normal_status, CFG.tumor_status):
            node_records.append(
                {
                    "CM": row.CM,
                    "cell_subtype": row.node,
                    "loading": float(loading_lookup.loc[(row.CM, row.node)]),
                    "rank": rank,
                    "cell_lineage": lineage_from_subtype(row.node),
                    "status_context": context,
                }
            )
    node_context = pd.DataFrame(node_records)
    node_context.loc[node_context["status_context"].eq(CFG.tumor_status)].to_csv(
        NODE_DIR / "tumor_network_nodes_from_H_df.csv", index=False
    )
    node_context.loc[node_context["status_context"].eq(CFG.normal_status)].to_csv(
        NODE_DIR / "normal_like_network_nodes_from_H_df.csv", index=False
    )

    union_top10 = top10.sort_values(["CM", "rank"])["cell_subtype"].drop_duplicates().tolist()
    for context in (CFG.normal_status, CFG.tumor_status):
        corr, qmat = corr_and_q_matrix(V_df.loc[status.eq(context), union_top10])
        corr.to_csv(NODE_DIR / f"{context}_node_node_correlation_matrix.csv")
        qmat.to_csv(NODE_DIR / f"{context}_node_node_correlation_q_matrix.csv")
        if context == CFG.tumor_status:
            corr.to_csv(NODE_DIR / "node_node_correlation_matrix.csv")
            qmat.to_csv(NODE_DIR / "node_node_correlation_q_matrix.csv")


def safe_corr(x: pd.Series, y: pd.Series, method: str) -> tuple[float, float, int]:
    ok = x.notna() & y.notna()
    n = int(ok.sum())
    if n < 3 or x[ok].nunique() < 2 or y[ok].nunique() < 2:
        return np.nan, np.nan, n
    if method == "spearman":
        r, p = spearmanr(x[ok], y[ok])
    elif method == "pearson":
        r, p = pearsonr(x[ok], y[ok])
    else:
        raise ValueError(method)
    return float(r), float(p), n


def run_association(
    epi_freq: pd.DataFrame,
    W_df: pd.DataFrame,
    sample_status: pd.DataFrame,
    inclusion: pd.DataFrame,
    method: str,
    out_dir: Path,
) -> pd.DataFrame:
    eligible = inclusion.index[inclusion["keep_for_epi_cm"]]
    common = epi_freq.index.intersection(W_df.index).intersection(sample_status.index).intersection(eligible)
    E = epi_freq.loc[common].astype(float)
    C = W_df.loc[common].astype(float)
    status = sample_status.loc[common, CFG.status_col].astype(str)
    rows = []
    stat_col = "rho" if method == "spearman" else "r"
    for context in (CFG.tumor_status, CFG.normal_status):
        samples = status.index[status.eq(context)]
        context_rows = []
        for epi in E.columns:
            for cm in C.columns:
                r, p, n = safe_corr(E.loc[samples, epi], C.loc[samples, cm], method)
                context_rows.append(
                    {
                        "method": method,
                        "status": context,
                        "CM": cm,
                        "epi_subtype": epi,
                        stat_col: r,
                        "p_value": p,
                        "n_samples": n,
                    }
                )
        context_df = pd.DataFrame(context_rows)
        context_df["q_value"] = bh_fdr(context_df["p_value"])
        rows.append(context_df)
        value = context_df.pivot(index="epi_subtype", columns="CM", values=stat_col).reindex(index=E.columns, columns=C.columns)
        pmat = context_df.pivot(index="epi_subtype", columns="CM", values="p_value").reindex(index=E.columns, columns=C.columns)
        qmat = context_df.pivot(index="epi_subtype", columns="CM", values="q_value").reindex(index=E.columns, columns=C.columns)
        context_stem = "normal-like" if context == CFG.normal_status else "tumor"
        if method == "spearman":
            prefix = f"balanced_joint_cm_epi_cm_association_{context_stem}"
            value.to_csv(out_dir / f"{prefix}_rho_matrix.csv")
            pmat.to_csv(out_dir / f"{prefix}_p_matrix.csv")
            qmat.to_csv(out_dir / f"{prefix}_q_matrix.csv")
        else:
            prefix = f"balanced_joint_cm_epi_cm_association_{context_stem}_pearson"
            value.to_csv(out_dir / f"{prefix}_r_matrix.csv")
            pmat.to_csv(out_dir / f"{prefix}_p_matrix.csv")
            qmat.to_csv(out_dir / f"{prefix}_q_matrix.csv")
    summary = pd.concat(rows, ignore_index=True)
    summary.to_csv(out_dir / f"epi_cm_{method}_all_pairs_long.csv", index=False)
    return summary


def write_provenance(obs: pd.DataFrame, non_epi_freq: pd.DataFrame, epi_freq: pd.DataFrame) -> None:
    params = {
        "input_h5ad": str(INPUT_H5AD),
        "code_file": str(Path(__file__).resolve()),
        "seed": SEED,
        "config": asdict(CFG),
        "n_cells": int(len(obs)),
        "n_samples_total": int(obs[CFG.sample_col].nunique()),
        "n_non_epi_subtypes": int(non_epi_freq.shape[1]),
        "n_epi_subtypes": int(epi_freq.shape[1]),
        "status_values": sorted(obs[CFG.status_col].unique().tolist()),
        "pearson_branch": "executed for compact-workflow completion checklist",
        "spatial_validation": "skipped: no spatial input supplied or requested",
        "gpu_note": "scikit-learn NMF, SciPy NNLS/correlations, and plotting have no canonical GPU backend in this skill; CPU used",
        "skill_resource_note": "Referenced inventory/source directories were absent; explicit canonical code in SKILL.md was used",
    }
    write_json(params, ANALYSIS_ROOT / "run_parameters.json")
    (ANALYSIS_ROOT / "readme.txt").write_text(
        "Input: projected full atlas h5ad from Block 02.\n"
        "Code: 01_run_cm_lineage_analysis.py.\n"
        "Outputs: frequency tables, balanced joint NMF K selection and canonical W/H, CM classification, node/edge tables, exhaustive Spearman and Pearson Epi-CM association tables.\n"
        "Status contract: obs['status'] contains exactly normal and tumor; no AT/HD grouping is used.\n"
        "Optional Block 04 spatial validation: skipped because no spatial input was supplied/requested.\n",
        encoding="utf-8",
    )
    versions = (
        f"python={sys.version.split()[0]}\n"
        f"platform={platform.platform()}\n"
        f"environment={sys.executable}\n"
        f"anndata={ad.__version__ if hasattr(ad, '__version__') else 'see importlib metadata'}\n"
        f"numpy={np.__version__}\n"
        f"pandas={pd.__version__}\n"
        f"scipy={scipy.__version__}\n"
        f"scikit-learn={sklearn.__version__}\n"
        f"seed={SEED}\n"
        "backend=CPU; canonical sklearn NMF and scipy NNLS/statistics have no GPU implementation in the skill\n"
        f"code={Path(__file__).resolve()}\n"
    )
    (ANALYSIS_ROOT / "package_versions.txt").write_text(versions, encoding="utf-8")
    pd.DataFrame(
        [
            {
                "step": "frequency_tables_and_status_balancing",
                "planned_backend": "CPU",
                "attempted_backend": "CPU",
                "status": "completed",
                "error_summary": "",
                "fallback_backend": "",
                "clean_input_reloaded": True,
                "final_backend_for_rerun": "CPU",
            },
            {
                "step": "sklearn_balanced_joint_nmf_and_scipy_nnls",
                "planned_backend": "CPU",
                "attempted_backend": "CPU",
                "status": "completed",
                "error_summary": "",
                "fallback_backend": "",
                "clean_input_reloaded": True,
                "final_backend_for_rerun": "CPU",
            },
            {
                "step": "correlations_and_mann_whitney_tests",
                "planned_backend": "CPU",
                "attempted_backend": "CPU",
                "status": "completed",
                "error_summary": "",
                "fallback_backend": "",
                "clean_input_reloaded": True,
                "final_backend_for_rerun": "CPU",
            },
        ]
    ).to_csv(ANALYSIS_ROOT / "gpu_backend_capability_summary.csv", index=False)


def main() -> None:
    ensure_dirs()
    obs = load_obs()
    non_epi_freq, epi_freq, sample_status, inclusion = prepare_frequency_tables(obs)
    V_df, W_df, H_df, _, _ = fit_balanced_joint_nmf(non_epi_freq, sample_status)
    build_nodes_and_edges(V_df, H_df, sample_status)
    run_association(epi_freq, W_df, sample_status, inclusion, "spearman", SPEARMAN_DIR)
    run_association(epi_freq, W_df, sample_status, inclusion, "pearson", PEARSON_DIR)
    write_provenance(obs, non_epi_freq, epi_freq)
    print(
        json.dumps(
            {
                "status": "completed",
                "n_cells": len(obs),
                "n_samples_cm": len(W_df),
                "n_epi_subtypes": epi_freq.shape[1],
                "n_non_epi_subtypes": non_epi_freq.shape[1],
                "selected_k": W_df.shape[1],
                "canonical_CMs": W_df.columns.tolist(),
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
