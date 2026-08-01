#!/usr/bin/env python3
"""Run all-pair tumor-only epithelial subtype versus CM Spearman analysis."""

from __future__ import annotations

import importlib.metadata
import json
import sys
import time
from pathlib import Path

import pandas as pd

from cm_analysis_common import (
    NMF_DIR, PREP_DIR, SEED, SPEARMAN_DIR, bh_fdr, cm_sort_key,
    safe_corr, write_json,
)


CODE_PATH = Path(__file__).resolve()
METHOD = "spearman"


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def q_star(q: float) -> str:
    if pd.isna(q):
        return ""
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return "ns"


def main() -> None:
    started = time.time()
    completion_path = SPEARMAN_DIR / "spearman_completion.json"
    if completion_path.exists():
        completion = json.loads(completion_path.read_text())
        if completion.get("status") == "completed":
            print(json.dumps({"status": "valid_existing_spearman_reused"}, indent=2))
            return
        raise FileExistsError("Existing Spearman completion is invalid.")
    SPEARMAN_DIR.mkdir(parents=True, exist_ok=True)
    mode = json.loads((PREP_DIR / "cm_status_mode.json").read_text())["detected_mode"]
    if mode != "tumor_only":
        raise ValueError(f"Expected tumor_only route, found {mode}")
    epi = pd.read_csv(PREP_DIR / "epi_subtype_frequency.csv", index_col=0)
    activity = pd.read_csv(NMF_DIR / "activity_df_sample_by_CM.csv", index_col=0)
    inclusion = pd.read_csv(PREP_DIR / "sample_inclusion_exclusion.csv").set_index("sample")
    sample_status = pd.read_csv(PREP_DIR / "sample_status.csv", index_col=0)
    eligible = inclusion.index[inclusion["keep_for_epi_cm"].astype(bool)]
    common = activity.index.intersection(epi.index).intersection(eligible)
    if len(common) != int(inclusion["keep_for_epi_cm"].sum()):
        raise ValueError("Eligible Epi-CM sample IDs do not align across inputs.")
    if set(sample_status.loc[common, "status"].astype(str)) != {"tumor"}:
        raise ValueError("Tumor-only Spearman input contains non-tumor samples.")
    epi = epi.loc[common]
    activity = activity.loc[common]
    epi_order = epi.columns.astype(str).tolist()
    cm_order = sorted(activity.columns.astype(str), key=cm_sort_key)
    activity = activity[cm_order]

    rows: list[dict[str, object]] = []
    scatter_rows: list[dict[str, object]] = []
    for epi_subtype in epi_order:
        for cm in cm_order:
            rho, p, n = safe_corr(epi[epi_subtype], activity[cm], METHOD)
            rows.append(
                {"context": "tumor", "epi_subtype": epi_subtype, "CM": cm,
                 "rho": rho, "p_value": p, "n_samples": n}
            )
            for sample in common:
                scatter_rows.append(
                    {"context": "tumor", "method": METHOD, "sample": sample,
                     "epi_subtype": epi_subtype, "CM": cm,
                     "epi_frequency": float(epi.loc[sample, epi_subtype]),
                     "cm_activity": float(activity.loc[sample, cm])}
                )
    association = pd.DataFrame(rows)
    association["q_value"] = bh_fdr(association["p_value"])
    association["significance"] = association["q_value"].map(q_star)
    association.to_csv(
        SPEARMAN_DIR / "balanced_joint_cm_epi_cm_association_tumor_spearman_long.csv",
        index=False,
    )
    pd.DataFrame(scatter_rows).to_csv(
        SPEARMAN_DIR / "balanced_joint_cm_epi_cm_association_tumor_spearman_scatter_source.csv",
        index=False,
    )
    for value, filename in [
        ("rho", "balanced_joint_cm_epi_cm_association_tumor_rho_matrix.csv"),
        ("p_value", "balanced_joint_cm_epi_cm_association_tumor_p_matrix.csv"),
        ("q_value", "balanced_joint_cm_epi_cm_association_tumor_q_matrix.csv"),
        ("significance", "balanced_joint_cm_epi_cm_association_tumor_significance_matrix.csv"),
    ]:
        association.pivot(index="epi_subtype", columns="CM", values=value).reindex(
            index=epi_order, columns=cm_order
        ).to_csv(SPEARMAN_DIR / filename)

    pd.DataFrame(
        [
            {"method": METHOD, "detected_mode": mode, "status_context": "tumor",
             "epi_frequency": str(PREP_DIR / "epi_subtype_frequency.csv"),
             "cm_activity": str(NMF_DIR / "activity_df_sample_by_CM.csv"),
             "sample_inclusion": str(PREP_DIR / "sample_inclusion_exclusion.csv"),
             "n_eligible_samples": len(common), "n_epi_subtypes": len(epi_order),
             "n_cms": len(cm_order), "n_pairs": len(association),
             "fdr_family": "all epithelial subtype x CM pairs within tumor Spearman context",
             "q_star_thresholds": "ns q>=0.05; * q<0.05; ** q<0.01; *** q<0.001",
             "pearson_branch_run": False, "code_file": str(CODE_PATH), "seed": SEED}
        ]
    ).to_csv(SPEARMAN_DIR / "spearman_parameters.csv", index=False)
    (SPEARMAN_DIR / "package_versions.txt").write_text(
        f"python={sys.version.split()[0]}\npandas={package_version('pandas')}\nscipy={package_version('scipy')}\nstatsmodels={package_version('statsmodels')}\ncode={CODE_PATH}\nseed={SEED}\n",
        encoding="utf-8",
    )
    completion = {
        "status": "completed", "detected_mode": mode,
        "status_context": "tumor", "method": METHOD,
        "n_eligible_samples": len(common), "n_epi_subtypes": len(epi_order),
        "n_cms": len(cm_order), "expected_pairs": len(epi_order) * len(cm_order),
        "completed_pairs": len(association),
        "bh_fdr_across_complete_pair_family": True,
        "scatter_source_rows": len(scatter_rows),
        "pearson_optional_branch_skipped": True,
        "elapsed_seconds": time.time() - started,
    }
    write_json(completion, completion_path)
    print(json.dumps(completion, indent=2), flush=True)


if __name__ == "__main__":
    main()
