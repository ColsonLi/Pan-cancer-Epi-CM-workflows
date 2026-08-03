#!/usr/bin/env python3
from __future__ import annotations

import json
import platform
import re
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from bs4 import BeautifulSoup


ANALYSIS_ROOT = Path("/mnt/disk18t/lr_xcy/riku/crc_val")
DATA_ROOT = Path("/mnt/disk18t/lr_xcy/riku/external_atlas_for_cm-epi_validation/CRC_S-BIAD2208")
LOAD_ROOT = DATA_ROOT / "input_datasets_extracted/load_datasets"
HARMONIZE_HTML = LOAD_ROOT / "harmonize_datasets/01-Harmonize_datasets.html"
SOURCE_METADATA = LOAD_ROOT / "harmonize_datasets/artifacts/merged_sample_metadata.csv"

OUT_DIR = ANALYSIS_ROOT / "epi-cm-core-workflow/tables/01-celltype_integration_clustering/01-select-samples"
CODE_FILE = ANALYSIS_ROOT / "epi-cm-core-workflow/codes/01-celltype_integration_clustering/01-select-samples/01_build_concat_dataset_and_sample_eligibility_tables.py"

DATASET_STATUS_CSV = OUT_DIR / "01_concat_dataset_used_and_not_used.csv"
ELIGIBLE_SAMPLES_CSV = OUT_DIR / "02_eligible_samples_for_merge.csv"
INELIGIBLE_SAMPLES_CSV = OUT_DIR / "03_ineligible_samples_for_merge.csv"
SUMMARY_CSV = OUT_DIR / "sample_selection_summary.csv"
README = OUT_DIR / "readme.txt"
VERSIONS = OUT_DIR / "package_versions.txt"


def norm_str(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def lower(series: pd.Series) -> pd.Series:
    return norm_str(series).str.lower()


def read_obs_column_h5py(path: Path, column: str, values: bool = False) -> list[str] | np.ndarray:
    with h5py.File(path, "r") as f:
        obs = f["obs"]
        if column not in obs:
            return [] if not values else np.array([], dtype=str)
        node = obs[column]
        if isinstance(node, h5py.Group) and "categories" in node:
            cats = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in node["categories"][:]], dtype=object)
            if not values:
                return cats.tolist()
            codes = node["codes"][:]
            out = np.empty(codes.shape[0], dtype=object)
            valid = codes >= 0
            out[valid] = cats[codes[valid]]
            out[~valid] = ""
            return out
        arr = node[:]
        if arr.dtype.kind == "S":
            arr = np.array([x.decode() for x in arr], dtype=object)
        else:
            arr = arr.astype(str)
        return arr if values else sorted(pd.unique(arr).tolist())


def read_concat_used_dataset_counts() -> dict[str, int]:
    soup = BeautifulSoup(HARMONIZE_HTML.read_text(errors="ignore"), "html.parser")
    used: dict[str, int] = {}
    for pre in soup.find_all("pre"):
        text = pre.get_text()
        if "Name: count, dtype: int64" not in text:
            continue
        if "Chen_2024" not in text or "Bian_2018" not in text:
            continue
        for line in text.splitlines():
            line = line.rstrip()
            if not line or line.startswith("dataset") or line.startswith("Name:"):
                continue
            match = re.match(r"^(.+?)\s+([0-9]+)$", line)
            if match:
                used[match.group(1).strip()] = int(match.group(2))
    if not used:
        raise RuntimeError(f"Could not parse concat-used dataset list from {HARMONIZE_HTML}")
    return used


def local_h5ad_index() -> pd.DataFrame:
    rows = []
    for path in sorted(LOAD_ROOT.rglob("*-adata.h5ad")):
        if "harmonize_datasets" in path.parts:
            continue
        try:
            datasets = read_obs_column_h5py(path, "dataset", values=False)
        except Exception as exc:
            datasets = [path.name.replace("-adata.h5ad", "")]
            rows.append({"dataset": datasets[0], "h5ad": str(path), "h5ad_file_stem": datasets[0], "read_note": repr(exc)})
            continue
        for dataset in datasets:
            rows.append(
                {
                    "dataset": dataset,
                    "h5ad": str(path),
                    "h5ad_file_stem": path.name.replace("-adata.h5ad", ""),
                    "read_note": "dataset_from_internal_obs",
                }
            )
    return pd.DataFrame(rows)


def sample_h5ad_mapping(datasets_needed: set[str]) -> tuple[pd.DataFrame, dict[tuple[str, str], str], dict[tuple[str, str], int]]:
    rows = []
    sample_to_h5ad: dict[tuple[str, str], str] = {}
    sample_to_ncells: dict[tuple[str, str], int] = {}
    for path in sorted(LOAD_ROOT.rglob("*-adata.h5ad")):
        if "harmonize_datasets" in path.parts:
            continue
        try:
            internal_datasets = set(read_obs_column_h5py(path, "dataset", values=False))
            if not (internal_datasets & datasets_needed):
                continue
            dataset_values = read_obs_column_h5py(path, "dataset", values=True)
            sample_values = read_obs_column_h5py(path, "sample_id", values=True)
            if len(dataset_values) == 0 or len(sample_values) == 0:
                rows.append({"h5ad": str(path), "status": "missing_dataset_or_sample_id_obs"})
                continue
            obs = pd.DataFrame({"dataset": dataset_values, "sample_id": sample_values})
            obs["dataset"] = norm_str(obs["dataset"])
            obs["sample_id"] = norm_str(obs["sample_id"])
            obs = obs[obs["dataset"].isin(datasets_needed)]
            if len(obs):
                vc = obs.value_counts(["dataset", "sample_id"])
                for (dataset, sample_id), count in vc.items():
                    sample_to_h5ad[(dataset, sample_id)] = str(path)
                    sample_to_ncells[(dataset, sample_id)] = int(count)
                for dataset, n in obs.groupby("dataset", observed=True).size().items():
                    rows.append({"dataset": dataset, "h5ad": str(path), "n_cells_dataset_in_h5ad": int(n), "status": "mapped"})
        except Exception as exc:
            rows.append({"h5ad": str(path), "status": "read_error", "error": repr(exc)})
    return pd.DataFrame(rows), sample_to_h5ad, sample_to_ncells


def add_patient_evidence(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    sample_type = lower(out.get("sample_type", pd.Series("", index=out.index)))
    sample_tissue = lower(out.get("sample_tissue", pd.Series("", index=out.index)))
    m_stage = lower(out.get("tumor_stage_TNM_M", pd.Series("", index=out.index)))
    n_stage = lower(out.get("tumor_stage_TNM_N", pd.Series("", index=out.index)))
    tnm = lower(out.get("tumor_stage_TNM", pd.Series("", index=out.index)))

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
        | n_stage.str.contains(r"npos|n\\+", regex=True, na=False)
    )
    out["sample_metastasis_evidence"] = metastasis
    out["sample_lymph_node_evidence"] = lymph
    out["sample_metastasis_or_lymph_node_evidence"] = metastasis | lymph

    patient_key = norm_str(out.get("patient_id", pd.Series("", index=out.index)))
    fallback = norm_str(out["dataset"]) + "." + norm_str(out["sample_id"])
    patient_key = patient_key.where(patient_key.ne(""), fallback)
    out["patient_key_for_exclusion"] = patient_key
    out["has_patient_metastasis_or_lymph_node_evidence"] = out.groupby(patient_key, observed=True)[
        "sample_metastasis_or_lymph_node_evidence"
    ].transform("any")
    return out


def reason_join(items: list[str]) -> str:
    return ";".join([x for x in items if x]) or "eligible"


def build_sample_tables(used_counts: dict[str, int]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(SOURCE_METADATA, low_memory=False)
    df = add_patient_evidence(df)
    df["concat_merge_used_dataset"] = norm_str(df["dataset"]).isin(set(used_counts))
    df["concat_merged_cell_count_by_dataset"] = norm_str(df["dataset"]).map(used_counts).fillna(0).astype(int)

    mapping, sample_to_h5ad, sample_to_ncells = sample_h5ad_mapping(set(norm_str(df.loc[df["concat_merge_used_dataset"], "dataset"])))
    df["h5ad"] = [
        sample_to_h5ad.get((str(dataset), str(sample_id)), "")
        for dataset, sample_id in zip(norm_str(df["dataset"]), norm_str(df["sample_id"]), strict=True)
    ]
    df["n_cells_in_source_h5ad"] = [
        sample_to_ncells.get((str(dataset), str(sample_id)), 0)
        for dataset, sample_id in zip(norm_str(df["dataset"]), norm_str(df["sample_id"]), strict=True)
    ]

    enrichment = lower(df.get("enrichment_cell_types", pd.Series("", index=df.index)))
    no_presort = enrichment.eq("naive") | enrichment.eq("") | enrichment.eq("nan")
    sample_type = lower(df.get("sample_type", pd.Series("", index=df.index)))
    sample_tissue = lower(df.get("sample_tissue", pd.Series("", index=df.index)))
    tumor_source = lower(df.get("tumor_source", pd.Series("", index=df.index)))
    medical_condition = lower(df.get("medical_condition", pd.Series("", index=df.index)))
    treatment_status = lower(df.get("treatment_status_before_resection", pd.Series("", index=df.index)))
    treatment_drug = lower(df.get("treatment_drug", pd.Series("", index=df.index)))

    is_crc_tumor = sample_tissue.eq("colon") & sample_type.eq("tumor") & medical_condition.eq("colorectal cancer")
    is_crc_adjacent_normal = (
        sample_tissue.eq("colon")
        & sample_type.eq("normal")
        & tumor_source.isin(["normal", "", "nan"])
        & medical_condition.eq("colorectal cancer")
    )
    is_healthy_normal_colon = (
        sample_tissue.eq("colon")
        & sample_type.eq("normal")
        & tumor_source.isin(["normal", "", "nan"])
        & medical_condition.eq("healthy")
    )
    is_requested_sample_type = is_crc_tumor | is_crc_adjacent_normal | is_healthy_normal_colon
    treatment_naive_no_drug = treatment_status.isin(["", "nan", "na", "none", "naive"]) & treatment_drug.isin(["", "nan", "na", "none", "naive"])
    is_crc_patient_source = medical_condition.eq("colorectal cancer")
    no_met_lymph = ~df["has_patient_metastasis_or_lymph_node_evidence"].fillna(False).astype(bool)
    has_h5ad_sample = norm_str(df["h5ad"]).ne("") & df["n_cells_in_source_h5ad"].gt(0)

    eligible = df["concat_merge_used_dataset"] & no_presort & no_met_lymph & is_requested_sample_type & treatment_naive_no_drug & has_h5ad_sample
    df["eligible_for_merge"] = eligible
    df["treatment_naive_no_drug_evidence"] = treatment_naive_no_drug
    df["sample_source_cn"] = np.select(
        [
            is_crc_tumor,
            is_crc_adjacent_normal,
            is_healthy_normal_colon,
        ],
        [
            "非分选非转移CRC肿瘤样本",
            "CRC患者癌旁/正常结肠样本",
            "健康人正常结肠样本",
        ],
        default="不符合主合并来源要求",
    )

    reasons = []
    for i, row in df.iterrows():
        rs = []
        if not bool(row["concat_merge_used_dataset"]):
            rs.append("dataset_not_used_in_concat_merge")
        if not bool(no_presort.loc[i]):
            rs.append(f"not_naive_or_celltype_presorted:enrichment_cell_types={row.get('enrichment_cell_types', '')}")
        if bool(row["has_patient_metastasis_or_lymph_node_evidence"]):
            parts = []
            if bool(row["sample_metastasis_evidence"]):
                parts.append("sample_metastasis_evidence")
            if bool(row["sample_lymph_node_evidence"]):
                parts.append("sample_lymph_node_evidence")
            rs.append("patient_has_metastasis_or_lymph_node_evidence:" + ",".join(parts))
        if not bool(is_requested_sample_type.loc[i]):
            rs.append(
                "not_requested_sample_type:"
                f"medical_condition={row.get('medical_condition', '')},sample_type={row.get('sample_type', '')},"
                f"sample_tissue={row.get('sample_tissue', '')},tumor_source={row.get('tumor_source', '')}"
            )
        if not bool(treatment_naive_no_drug.loc[i]):
            rs.append(
                "not_treatment_naive_or_has_drug_evidence:"
                f"treatment_status_before_resection={row.get('treatment_status_before_resection', '')},"
                f"treatment_drug={row.get('treatment_drug', '')}"
            )
        if bool(row["concat_merge_used_dataset"]) and not bool(has_h5ad_sample.loc[i]):
            rs.append("selected_dataset_but_sample_not_found_in_local_h5ad")
        reasons.append(reason_join(rs))
    df["eligibility_reason"] = reasons
    df.loc[df["eligible_for_merge"], "eligibility_reason"] = "eligible"
    df["status"] = np.where(
        lower(df["sample_type"]).eq("tumor"),
        "primary_crc_no_metastasis_no_lymph_node",
        np.where(lower(df["sample_type"]).eq("normal"), "normal_colon_no_metastasis_no_lymph_node", "not_selected"),
    )
    df["include_reason"] = np.where(
        df["eligible_for_merge"],
        "concat_used_no_presort_no_metastasis_no_lymph_node_tumor_or_normal_colon",
        "",
    )

    keep_cols = [
        "sample_id",
        "dataset",
        "study_id",
        "patient_id",
        "sample_type",
        "medical_condition",
        "sample_tissue",
        "anatomic_region",
        "anatomic_location",
        "tumor_source",
        "tumor_stage_TNM",
        "tumor_stage_TNM_T",
        "tumor_stage_TNM_N",
        "tumor_stage_TNM_M",
        "platform",
        "matrix_type",
        "enrichment_cell_types",
        "tissue_cell_state",
        "tissue_processing_lab",
        "hospital_location",
        "country",
        "treatment_status_before_resection",
        "treatment_drug",
        "treatment_response",
        "RECIST",
        "study_pmid",
        "study_doi",
        "NCBI_BioProject_accession",
        "SRA_sample_accession",
        "GEO_sample_accession",
        "ENA_sample_accession",
        "synapse_sample_accession",
        "concat_merge_used_dataset",
        "concat_merged_cell_count_by_dataset",
        "eligible_for_merge",
        "treatment_naive_no_drug_evidence",
        "eligibility_reason",
        "sample_source_cn",
        "status",
        "include_reason",
        "patient_key_for_exclusion",
        "sample_metastasis_evidence",
        "sample_lymph_node_evidence",
        "has_patient_metastasis_or_lymph_node_evidence",
        "n_cells_in_source_h5ad",
        "h5ad",
    ]
    keep_cols = [c for c in keep_cols if c in df.columns]
    eligible_df = df.loc[df["eligible_for_merge"], keep_cols].drop_duplicates(["dataset", "sample_id", "h5ad"])
    ineligible_df = df.loc[~df["eligible_for_merge"], keep_cols].drop_duplicates(["dataset", "sample_id", "h5ad", "eligibility_reason"])
    return eligible_df.sort_values(["dataset", "sample_type", "sample_id"]), ineligible_df.sort_values(["dataset", "sample_id"]), mapping


def build_dataset_status(used_counts: dict[str, int], mapping: pd.DataFrame) -> pd.DataFrame:
    local = local_h5ad_index()
    all_datasets = sorted(set(local["dataset"]) | set(used_counts))
    local_map = local.drop_duplicates("dataset").set_index("dataset")["h5ad"].to_dict()
    rows = []
    known_not_used_reason = {
        "HTAN_MSK": "not present in harmonize concat dataset value_counts; local later Moorman/HTAN_MSK input not used by concat merge",
        "Liu_2022_CD45Pos": "harmonize notebook note: only normalized counts available; not present in concat output",
        "Masuda_2022_CD3Pos": "harmonize notebook note: no sample mapping and some samples stimulated in culture; not present in concat output",
    }
    for dataset in all_datasets:
        used = dataset in used_counts
        local_h5ad = local_map.get(dataset, "")
        if used:
            reason = "used_in_harmonize_concat_merge"
            if not local_h5ad:
                reason = "used_in_harmonize_concat_merge_but_no_matching_local_input_h5ad_found"
        else:
            reason = known_not_used_reason.get(dataset, "not present in harmonize concat dataset value_counts")
        rows.append(
            {
                "dataset": dataset,
                "local_h5ad": local_h5ad,
                "used_in_concat_merge": used,
                "concat_merged_cell_count": int(used_counts.get(dataset, 0)),
                "reason": reason,
            }
        )
    return pd.DataFrame(rows).sort_values(["used_in_concat_merge", "dataset"], ascending=[False, True])


def write_reports(dataset_status: pd.DataFrame, eligible: pd.DataFrame, ineligible: pd.DataFrame) -> None:
    summary = pd.DataFrame(
        [
            {"metric": "datasets_total_in_status_table", "value": int(len(dataset_status))},
            {"metric": "datasets_used_in_concat_merge", "value": int(dataset_status["used_in_concat_merge"].sum())},
            {"metric": "datasets_not_used_in_concat_merge", "value": int((~dataset_status["used_in_concat_merge"]).sum())},
            {"metric": "eligible_samples_total", "value": int(len(eligible))},
            {"metric": "eligible_datasets_total", "value": int(eligible["dataset"].nunique())},
            {"metric": "eligible_tumor_samples", "value": int(lower(eligible["sample_type"]).eq("tumor").sum())},
            {"metric": "eligible_normal_samples", "value": int(lower(eligible["sample_type"]).eq("normal").sum())},
            {
                "metric": "eligible_crc_tumor_samples",
                "value": int((eligible["sample_source_cn"] == "非分选非转移CRC肿瘤样本").sum()),
            },
            {
                "metric": "eligible_crc_adjacent_or_normal_colon_samples",
                "value": int((eligible["sample_source_cn"] == "CRC患者癌旁/正常结肠样本").sum()),
            },
            {
                "metric": "eligible_healthy_normal_colon_samples",
                "value": int((eligible["sample_source_cn"] == "健康人正常结肠样本").sum()),
            },
            {
                "metric": "eligible_treatment_naive_no_drug_samples",
                "value": int(eligible["treatment_naive_no_drug_evidence"].fillna(False).astype(bool).sum()),
            },
            {"metric": "eligible_source_cells_total", "value": int(eligible["n_cells_in_source_h5ad"].sum())},
            {"metric": "ineligible_samples_total", "value": int(len(ineligible))},
        ]
    )
    summary.to_csv(SUMMARY_CSV, index=False)
    README.write_text(
        "\n".join(
            [
                "01-select-samples",
                "",
                "Three requested CSV files:",
                f"1. Dataset concat-use status: {DATASET_STATUS_CSV}",
                f"2. Eligible samples for merge: {ELIGIBLE_SAMPLES_CSV}",
                f"3. Ineligible samples: {INELIGIBLE_SAMPLES_CSV}",
                "",
                "Concat-used datasets are parsed from the source harmonize notebook HTML:",
                str(HARMONIZE_HTML),
                "Specifically from the post-concat adata.obs['dataset'].value_counts() output, not from final_crc_atlas-adata.h5ad.",
                "",
                "Eligible sample definition:",
                "- dataset used in the harmonize concat merge",
                "- no cell-type presort/enrichment: enrichment_cell_types is naive/blank",
                "- no patient-level metastasis or lymph-node evidence",
                "- no drug-treatment evidence: treatment_status_before_resection is blank/naive and treatment_drug is blank/naive/none",
                "- requested samples are exactly: CRC tumor, CRC patient adjacent/normal colon, healthy human normal colon",
                "- exclude colorectal polyp, inflammation/IBD, blood, lymph node, liver/metastasis, and sorted/enriched datasets",
                "- local h5ad contains the sample_id",
                "",
                "中文标准:",
                "- 只合并三类样本：非分选非转移CRC肿瘤、CRC患者癌旁/正常结肠、健康人正常结肠。",
                "- 所有合并样本不能有药物治疗证据：treatment_status_before_resection 为空或 naive，且 treatment_drug 为空/naive/none。",
                "- 不合并：息肉患者、炎症/IBD患者、血液、淋巴结、转移灶、分选/富集样本。",
                "",
                f"Summary: {SUMMARY_CSV}",
            ]
        )
        + "\n"
    )
    VERSIONS.write_text(
        "\n".join(
            [
                f"python={sys.version.split()[0]}",
                f"python_executable={sys.executable}",
                f"platform={platform.platform()}",
                f"anndata={ad.__version__}",
                f"numpy={np.__version__}",
                f"pandas={pd.__version__}",
                f"beautifulsoup4_available=True",
                f"code_file={CODE_FILE}",
            ]
        )
        + "\n"
    )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    used_counts = read_concat_used_dataset_counts()
    eligible, ineligible, mapping = build_sample_tables(used_counts)
    dataset_status = build_dataset_status(used_counts, mapping)
    dataset_status.to_csv(DATASET_STATUS_CSV, index=False)
    eligible.to_csv(ELIGIBLE_SAMPLES_CSV, index=False)
    ineligible.to_csv(INELIGIBLE_SAMPLES_CSV, index=False)
    write_reports(dataset_status, eligible, ineligible)
    print(
        json.dumps(
            {
                "dataset_status_csv": str(DATASET_STATUS_CSV),
                "eligible_samples_csv": str(ELIGIBLE_SAMPLES_CSV),
                "ineligible_samples_csv": str(INELIGIBLE_SAMPLES_CSV),
                "datasets_used": int(dataset_status["used_in_concat_merge"].sum()),
                "datasets_not_used": int((~dataset_status["used_in_concat_merge"]).sum()),
                "eligible_samples": int(len(eligible)),
                "eligible_datasets": int(eligible["dataset"].nunique()),
                "eligible_source_cells": int(eligible["n_cells_in_source_h5ad"].sum()),
            },
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
