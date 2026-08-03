#!/usr/bin/env python3
"""Reaggregate and plot the non-metastatic primary-tumor spatial cohort.

Tangram and per-region tests are intentionally reused.  A spatial patient is
excluded when the author metadata for the linked patient contains evidence of
distant metastasis or lymph-node involvement, using the same fields and
patterns as the single-cell sample-selection script.
"""

from __future__ import annotations

import hashlib
import json
import math
import platform
import re
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import seaborn as sns
from scipy import stats


ROOT = Path("/mnt/disk18t/lr_xcy/riku/crc_codex/epi-cm-core-workflow")
BLOCK = "04-spatial-validation-optional-sig_genes"
TASK = "05-canonical-epifrac-composition-v2"
BASE_TABLE = ROOT / f"tables/{BLOCK}/{TASK}"
SOURCE_STATS = BASE_TABLE / "03-all-pair-statistics-and-plots"
SOURCE_SCOPE = BASE_TABLE / "01-reference-and-manifest/spatial_sample_scope_manifest.csv"
SPATIAL_H5AD = ROOT / "h5ad/04-spatial-validation-optional/01-input-audit/adata_xenium_unprocessed.h5ad"
AUTHOR_METADATA = Path(
    "/mnt/disk18t/lr_xcy/riku/external_atlas_for_cm-epi_validation/CRC_S-BIAD2208/"
    "input_datasets_extracted/load_datasets/harmonize_datasets/artifacts/merged_sample_metadata.csv"
)

OUT_NAME = "04-nonmetastatic-cohort-reaggregation"
TABLE_DIR = BASE_TABLE / OUT_NAME
FIG_DIR = ROOT / f"figures/{BLOCK}/{TASK}/03-all-pair-statistics-and-plots/stouffer_heatmaps"

plt.rcParams.update(
    {
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "font.family": "DejaVu Sans",
        "font.size": 8,
        "axes.titlesize": 9,
    }
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def lower(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip().str.lower()


def decode(values) -> list[str]:
    return [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values]


def read_obs_column(obs: h5py.Group, key: str) -> np.ndarray:
    node = obs[key]
    if isinstance(node, h5py.Group) and "categories" in node and "codes" in node:
        categories = decode(node["categories"][()])
        codes = node["codes"][()].astype(int)
        return np.asarray([categories[code] if code >= 0 else "" for code in codes], dtype=object)
    return np.asarray(decode(node[()]), dtype=object)


def spatial_patient_mapping() -> pd.DataFrame:
    with h5py.File(SPATIAL_H5AD, "r") as handle:
        obs = handle["obs"]
        frame = pd.DataFrame(
            {
                "patient_id": read_obs_column(obs, "patient_id"),
                "sample": read_obs_column(obs, "name"),
                "stage": read_obs_column(obs, "stage"),
                "crca_patient_id": read_obs_column(obs, "crca_patient_id"),
            }
        )
    frame = frame.drop_duplicates()
    conflicts = frame.groupby("patient_id").agg(
        n_crca_patient_ids=("crca_patient_id", lambda x: x[x != ""].nunique()),
        n_stages=("stage", lambda x: x[x != ""].nunique()),
    )
    if (conflicts[["n_crca_patient_ids", "n_stages"]] > 1).any().any():
        raise ValueError(f"Non-unique spatial patient mapping:\n{conflicts}")
    return (
        frame.groupby("patient_id", as_index=False)
        .agg(
            crca_patient_id=("crca_patient_id", lambda x: next((v for v in x if v), "")),
            stage=("stage", lambda x: next((v for v in x if v), "")),
        )
        .sort_values("patient_id", kind="stable")
    )


def author_patient_evidence() -> pd.DataFrame:
    metadata = pd.read_csv(AUTHOR_METADATA, low_memory=False)
    sample_type = lower(metadata["sample_type"])
    sample_tissue = lower(metadata["sample_tissue"])
    m_stage = lower(metadata["tumor_stage_TNM_M"])
    n_stage = lower(metadata["tumor_stage_TNM_N"])
    tnm = lower(metadata["tumor_stage_TNM"])
    metastasis = (
        sample_type.str.contains("metastasis", na=False)
        | sample_tissue.eq("liver")
        | m_stage.str.match(r"^m1", na=False)
        | tnm.str.match(r"^iv", na=False)
    )
    lymph = (
        sample_type.eq("lymph node")
        | sample_tissue.str.contains("lymph", na=False)
        | n_stage.str.match(r"^n[12]", na=False)
        | n_stage.str.contains(r"npos|n\+", regex=True, na=False)
    )
    status = lower(metadata["treatment_status_before_resection"])
    drug = lower(metadata["treatment_drug"])
    allowed = {"", "nan", "na", "none", "naive"}
    treatment_ok = status.isin(allowed) & drug.isin(allowed)
    evidence = pd.DataFrame(
        {
            "crca_patient_id": metadata["patient_id"].fillna("").astype(str),
            "metastasis_evidence": metastasis,
            "lymph_node_evidence": lymph,
            "treatment_naive_no_drug": treatment_ok,
            "tnm": metadata["tumor_stage_TNM"].fillna("").astype(str),
            "n_stage": metadata["tumor_stage_TNM_N"].fillna("").astype(str),
            "m_stage": metadata["tumor_stage_TNM_M"].fillna("").astype(str),
            "treatment_status": metadata["treatment_status_before_resection"].fillna("").astype(str),
            "treatment_drug": metadata["treatment_drug"].fillna("").astype(str),
        }
    )
    evidence = evidence[evidence["crca_patient_id"] != ""]

    def joined(values: pd.Series) -> str:
        cleaned = sorted({str(value).strip() for value in values if str(value).strip() and str(value) != "nan"})
        return ";".join(cleaned)

    return evidence.groupby("crca_patient_id", as_index=False).agg(
        author_metadata_rows=("crca_patient_id", "size"),
        metastasis_evidence=("metastasis_evidence", "any"),
        lymph_node_evidence=("lymph_node_evidence", "any"),
        treatment_naive_no_drug=("treatment_naive_no_drug", "all"),
        author_tnm=("tnm", joined),
        author_n_stage=("n_stage", joined),
        author_m_stage=("m_stage", joined),
        author_treatment_status=("treatment_status", joined),
        author_treatment_drug=("treatment_drug", joined),
    )


def build_matched_scope() -> pd.DataFrame:
    scope = pd.read_csv(SOURCE_SCOPE)
    mapping = spatial_patient_mapping()
    evidence = author_patient_evidence()
    out = scope.merge(mapping, on="patient_id", how="left", validate="many_to_one", suffixes=("", "_h5ad"))
    out = out.merge(evidence, on="crca_patient_id", how="left", validate="many_to_one")
    out["author_metadata_available"] = out["author_metadata_rows"].notna()
    for column in ["metastasis_evidence", "lymph_node_evidence"]:
        out[column] = out[column].map(lambda value: bool(value) if pd.notna(value) else False)
    out["include_matched_all_samples"] = ~(out["metastasis_evidence"] | out["lymph_node_evidence"])
    out["include_matched_tumor_only"] = out["include_matched_all_samples"] & out["status"].eq("tumor")
    out["exclusion_reason"] = np.select(
        [
            out["metastasis_evidence"] & out["lymph_node_evidence"],
            out["metastasis_evidence"],
            out["lymph_node_evidence"],
        ],
        ["patient_has_metastasis_and_lymph_node_evidence", "patient_has_metastasis_evidence", "patient_has_lymph_node_evidence"],
        default="",
    )
    out["treatment_evidence_note"] = np.where(
        out["author_metadata_available"],
        np.where(
            out["treatment_naive_no_drug"].map(lambda value: bool(value) if pd.notna(value) else False),
            "author_metadata=naive_no_drug",
            "author_metadata=treated_or_drug",
        ),
        "not_available_in_spatial_h5ad_or_linked_author_metadata",
    )
    return out.sort_values("sample", kind="stable").reset_index(drop=True)


def bh(pvalues: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvalues, errors="coerce").to_numpy(dtype=float)
    out = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    values = p[valid]
    if len(values):
        order = np.argsort(values)
        ranked = values[order]
        adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
        adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
        restored = np.empty_like(adjusted)
        restored[order] = np.minimum(adjusted, 1.0)
        out[valid] = restored
    return pd.Series(out, index=pvalues.index)


def combine_stouffer(fisher: pd.DataFrame, scope_name: str) -> pd.DataFrame:
    rows = []
    for (cm, epi), group in fisher.groupby(["CM", "epi_subtype"], sort=False):
        z = pd.to_numeric(group["signed_z"], errors="coerce").dropna().to_numpy(dtype=float)
        combined = float(z.sum() / math.sqrt(len(z))) if len(z) else np.nan
        p_value = float(2.0 * stats.norm.sf(abs(combined))) if len(z) else np.nan
        rows.append(
            {
                "scope": scope_name,
                "CM": cm,
                "epi_subtype": epi,
                "n_samples": len(z),
                "combined_signed_z": combined,
                "p_value": p_value,
            }
        )
    result = pd.DataFrame(rows)
    result["q_value_bh"] = bh(result["p_value"])
    return result


def cm_sort_key(value: str) -> int:
    match = re.search(r"CM(\d+)", str(value))
    return int(match.group(1)) if match else 999


def q_label(q_value: float) -> str:
    if not np.isfinite(q_value) or q_value >= 0.05:
        return "ns"
    if q_value < 0.001:
        return "***"
    if q_value < 0.01:
        return "**"
    return "*"


def plot_tumor_heatmap(table: pd.DataFrame) -> tuple[Path, Path]:
    z = table.pivot(index="epi_subtype", columns="CM", values="combined_signed_z")
    q = table.pivot(index="epi_subtype", columns="CM", values="q_value_bh")
    rows = sorted(z.index)
    columns = sorted(z.columns, key=cm_sort_key)
    z = z.reindex(index=rows, columns=columns)
    q = q.reindex(index=rows, columns=columns)
    annotation = q.map(q_label)
    finite = z.to_numpy()[np.isfinite(z.to_numpy())]
    vmax = max(float(np.max(np.abs(finite))) if len(finite) else 1.0, 1.0)

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig, axis = plt.subplots(figsize=(10.4, 6.6))
    sns.heatmap(
        z,
        cmap="coolwarm",
        center=0,
        vmin=-vmax,
        vmax=vmax,
        linewidths=0.25,
        linecolor="white",
        annot=annotation,
        fmt="",
        annot_kws={"fontsize": 6.5, "color": "black"},
        square=True,
        cbar_kws={"label": "Combined signed Stouffer Z"},
        ax=axis,
    )
    axis.set_title(
        "Tumor spatial validation in the non-metastatic cohort\n"
        "9 patients, 17 primary-tumor regions",
        fontweight="bold",
    )
    axis.set_xlabel("CM")
    axis.set_ylabel("Epithelial subtype")
    axis.tick_params(axis="x", rotation=45)
    axis.tick_params(axis="y", rotation=0)
    fig.tight_layout()
    # Replace the canonical tumor-only integrated heatmap in place.  The
    # all-samples heatmap remains untouched.
    stem = FIG_DIR / "tumor_only_stouffer_signedZ_qstars"
    pdf = stem.with_suffix(".pdf")
    svg = stem.with_suffix(".svg")
    fig.savefig(pdf, dpi=300, bbox_inches="tight")
    fig.savefig(svg, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return pdf, svg


def compare_with_original(strict: pd.DataFrame, original_path: Path, scope_name: str) -> pd.DataFrame:
    original = pd.read_csv(original_path)
    left = original.rename(
        columns={
            "n_samples": "original_n_samples",
            "combined_signed_z": "original_combined_signed_z",
            "p_value": "original_p_value",
            "q_value_bh": "original_q_value_bh",
        }
    ).drop(columns=["scope"])
    right = strict.rename(
        columns={
            "n_samples": "matched_n_samples",
            "combined_signed_z": "matched_combined_signed_z",
            "p_value": "matched_p_value",
            "q_value_bh": "matched_q_value_bh",
        }
    ).drop(columns=["scope"])
    comparison = left.merge(right, on=["CM", "epi_subtype"], validate="one_to_one")
    comparison.insert(0, "scope", scope_name)
    comparison["original_significant_q_lt_0p05"] = comparison["original_q_value_bh"] < 0.05
    comparison["matched_significant_q_lt_0p05"] = comparison["matched_q_value_bh"] < 0.05
    comparison["significance_changed"] = (
        comparison["original_significant_q_lt_0p05"] != comparison["matched_significant_q_lt_0p05"]
    )
    comparison["direction_flipped"] = (
        np.sign(comparison["original_combined_signed_z"]) != np.sign(comparison["matched_combined_signed_z"])
    )
    return comparison


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    scope = build_matched_scope()
    tumor_samples = set(scope.loc[scope["include_matched_tumor_only"], "sample"].astype(str))
    included_patients = set(scope.loc[scope["include_matched_all_samples"], "patient_id"].astype(str))
    excluded_patients = set(scope.loc[~scope["include_matched_all_samples"], "patient_id"].astype(str))
    if included_patients != set("abcfjlmno"):
        raise ValueError(f"Unexpected included patients: {sorted(included_patients)}")
    if excluded_patients != set("deghik"):
        raise ValueError(f"Unexpected excluded patients: {sorted(excluded_patients)}")
    if len(tumor_samples) != 17:
        raise ValueError(f"Expected 17 retained tumor regions, found {len(tumor_samples)}")

    fisher_all = pd.read_csv(SOURCE_STATS / "percentile_quadrant_fisher_per_sample.csv")
    spearman_all = pd.read_csv(SOURCE_STATS / "per_sample_spearman.csv")
    fisher = fisher_all[fisher_all["sample"].astype(str).isin(tumor_samples)].copy()
    spearman = spearman_all[spearman_all["sample"].astype(str).isin(tumor_samples)].copy()
    if fisher["sample"].nunique() != 17 or len(fisher) != 17 * 15 * 11:
        raise ValueError(f"Unexpected tumor-only Fisher shape: {fisher.shape}")
    if spearman["sample"].nunique() != 17 or len(spearman) != 17 * 15 * 11:
        raise ValueError(f"Unexpected tumor-only Spearman shape: {spearman.shape}")

    tumor_only = combine_stouffer(fisher, "nonmetastatic-primary-tumor-only")
    retained_scope = scope.loc[scope["include_matched_tumor_only"]].copy()
    if not retained_scope["status"].eq("tumor").all() or not retained_scope["tissue_region"].isin(["core", "margin"]).all():
        raise ValueError("Retained scope contains non-tumor or non-primary-tumor regions")

    retained_scope.to_csv(TABLE_DIR / "tumor_spatial_sample_scope_nonmetastatic.csv", index=False)
    fisher.to_csv(TABLE_DIR / "percentile_quadrant_fisher_per_sample_tumor_only.csv", index=False)
    spearman.to_csv(TABLE_DIR / "per_sample_spearman_tumor_only.csv", index=False)
    tumor_dir = TABLE_DIR / "statistics/tumor-only"
    tumor_dir.mkdir(parents=True, exist_ok=True)
    tumor_out = tumor_dir / "percentile_quadrant_fisher_sample_stouffer_nonmetastatic_tumor_only.csv"
    tumor_only.to_csv(tumor_out, index=False)

    comparison = compare_with_original(
        tumor_only,
        SOURCE_STATS / "statistics/tumor-only/percentile_quadrant_fisher_sample_stouffer_tumor_only.csv",
        "tumor-only",
    )
    comparison.to_csv(TABLE_DIR / "comparison_with_original_tumor_only_aggregation.csv", index=False)

    patient_table = scope[
        [
            "patient_id",
            "crca_patient_id",
            "stage",
            "author_metadata_available",
            "metastasis_evidence",
            "lymph_node_evidence",
            "author_tnm",
            "author_n_stage",
            "author_m_stage",
            "treatment_evidence_note",
            "include_matched_all_samples",
            "exclusion_reason",
        ]
    ].drop_duplicates("patient_id").sort_values("patient_id", kind="stable")
    patient_table.to_csv(TABLE_DIR / "patient_level_eligibility_audit.csv", index=False)

    pdf, svg = plot_tumor_heatmap(tumor_only)
    heatmap_manifest = {
        "input": str(tumor_out),
        "input_sha256": sha256(tumor_out),
        "scope": "nonmetastatic-primary-tumor-only",
        "n_patients": len(included_patients),
        "n_primary_tumor_regions": len(tumor_samples),
        "n_cm_epi_pairs": len(tumor_only),
        "n_significant_q_lt_0p05": int((tumor_only["q_value_bh"] < 0.05).sum()),
        "pdf": str(pdf),
        "pdf_sha256": sha256(pdf),
        "svg": str(svg),
        "svg_sha256": sha256(svg),
    }
    (TABLE_DIR / "heatmap_manifest.json").write_text(
        json.dumps(heatmap_manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    summary = {
        "purpose": "Tumor-only reaggregation after applying the single-cell patient-level metastasis/lymph-node evidence exclusion rule.",
        "tangram_rerun": False,
        "per_region_statistics_recomputed": False,
        "retained_tissue_regions": ["core", "margin"],
        "adjacent_normal_regions_retained": 0,
        "heatmaps_generated": 1,
        "heatmap_files_generated": 2,
        "n_original_patients": int(scope["patient_id"].nunique()),
        "n_included_patients": len(included_patients),
        "included_patients": sorted(included_patients),
        "n_excluded_patients": len(excluded_patients),
        "excluded_patients": sorted(excluded_patients),
        "n_included_primary_tumor_regions": len(tumor_samples),
        "n_pair_region_rows_tumor_only": len(fisher),
        "n_pairs": len(tumor_only),
        "n_significant_pairs_tumor_only_q_lt_0p05": int((tumor_only["q_value_bh"] < 0.05).sum()),
        "n_significance_changes_tumor_only_vs_original": int(comparison["significance_changed"].sum()),
        "n_direction_flips_tumor_only_vs_original": int(comparison["direction_flipped"].sum()),
        "unlinked_treatment_metadata_patients": sorted(
            patient_table.loc[~patient_table["author_metadata_available"], "patient_id"].astype(str).tolist()
        ),
        "input_sha256": {
            "source_scope": sha256(SOURCE_SCOPE),
            "source_fisher": sha256(SOURCE_STATS / "percentile_quadrant_fisher_per_sample.csv"),
            "source_spearman": sha256(SOURCE_STATS / "per_sample_spearman.csv"),
            "author_metadata": sha256(AUTHOR_METADATA),
        },
        "outputs": {
            "tumor_scope_csv": str(TABLE_DIR / "tumor_spatial_sample_scope_nonmetastatic.csv"),
            "tumor_stouffer_csv": str(tumor_out),
            "heatmap_pdf": str(pdf),
            "heatmap_svg": str(svg),
        },
        "python": platform.python_version(),
        "scipy": scipy.__version__,
    }
    (TABLE_DIR / "reaggregation_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    (TABLE_DIR / "README.txt").write_text(
        "Non-metastatic primary-tumor spatial reaggregation\n"
        "\n"
        "Only the original Xenium tumor-core and invasive-margin regions are retained.\n"
        "Adjacent-normal regions are excluded. Tangram and per-region tests were not\n"
        "rerun. Existing per-region\n"
        "Fisher signed Z values were subset and combined with the same unweighted signed\n"
        "Stouffer method. Patients d/e/g/h were excluded for lymph-node evidence; patients\n"
        "i/k were excluded for both metastatic and lymph-node evidence. Included patients\n"
        "are a/b/c/f/j/l/m/n/o, contributing 17 primary-tumor regions.\n"
        "The canonical tumor-only integrated heatmap is replaced in place; the canonical\n"
        "all-samples heatmap and original per-region/statistics tables remain unchanged.\n"
        "\n"
        "Treatment metadata are author-confirmed naive/no-drug for linked included patients\n"
        "f/j/l/m/n. Treatment fields are unavailable for spatial patients a/b/c/o; this is\n"
        "recorded in the cohort audit and is not silently treated as confirmed naive.\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
