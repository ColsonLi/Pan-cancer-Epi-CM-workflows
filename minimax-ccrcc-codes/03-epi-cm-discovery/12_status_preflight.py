"""Module 03 step 12: preflight cm_status_mode outputs per updated SKILL.md.

Required outputs:
  tables/03-epi-cm-discovery/cm_status_mode_detection.csv
  tables/03-epi-cm-discovery/cm_status_mode.json

Detects the Block 03 mode (tumor_normal / tumor_only / normal_only /
unsupported) from the trusted sample_status source.

In this project the canonical sample status for module 03 comes from the
activity_df_sample_by_CM.csv table (sample_id, status, ...) which records
single-cell sample status (61 tumor + 21 normal-like -> detected mode = tumor_normal).
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
TAB = ROOT / "epi-cm-core-workflow/tables/03-epi-cm-discovery"
ACT = TAB / "activity_df_sample_by_CM.csv"
OUT_DIR = TAB
OUT_DIR.mkdir(parents=True, exist_ok=True)


def _status_token(value):
    return re.sub(r"[^a-z0-9]+", " ", str(value).strip().lower()).strip()


def detect_block03_status_mode(
    sample_status,
    out_dir,
    sample_col="sample_id",
    status_col="status",
    status_source="sample_metadata",
    user_declared_tumor_only=False,
    user_declared_normal_only=False,
    status_aliases=None,
):
    """FIXED per updated SKILL.md Module 03 status mode preflight."""
    frame = sample_status.copy()
    if sample_col not in frame.columns:
        if frame.index.name is None:
            raise ValueError(f"Missing {sample_col!r} column and named sample index.")
        frame = frame.reset_index().rename(columns={frame.index.name: sample_col})

    frame[sample_col] = frame[sample_col].astype(str)
    if frame[sample_col].duplicated().any():
        if status_col not in frame.columns:
            raise ValueError("Duplicate sample rows require an explicit status column.")
        conflicts = frame.groupby(sample_col)[status_col].nunique(dropna=False)
        conflicts = conflicts[conflicts > 1]
        if len(conflicts):
            raise ValueError(
                "Conflicting status labels within samples: "
                + ", ".join(conflicts.index.astype(str))
            )
        frame = frame.drop_duplicates(sample_col, keep="first")

    if status_col not in frame.columns or frame[status_col].isna().all():
        if user_declared_tumor_only == user_declared_normal_only:
            raise ValueError(
                "Missing sample status. Supply trusted metadata or explicitly declare "
                "exactly one of tumor-only or normal-only."
            )
        if user_declared_tumor_only:
            frame[status_col] = "tumor"
            status_source = "user_declared_tumor_only"
        else:
            frame[status_col] = "normal-like"
            status_source = "user_declared_normal_only"
    elif frame[status_col].isna().any():
        missing_samples = frame.loc[frame[status_col].isna(), sample_col].astype(str)
        raise ValueError("Missing status for samples: " + ", ".join(missing_samples))

    aliases = {
        "tumor": "tumor", "tumour": "tumor",
        "primary tumor": "tumor", "primary tumour": "tumor",
        "metastatic tumor": "tumor", "metastatic tumour": "tumor",
        "metastasis": "tumor",
        "normal": "normal-like", "normal like": "normal-like",
        "adjacent normal": "normal-like", "normal adjacent": "normal-like",
    }
    if status_aliases is not None:
        aliases.update({_status_token(k): v for k, v in status_aliases.items()})

    frame["raw_status"] = frame[status_col].astype(str)
    frame["canonical_status"] = frame["raw_status"].map(
        lambda value: aliases.get(_status_token(value))
    )
    unresolved = frame.loc[frame["canonical_status"].isna(), "raw_status"].unique()
    if len(unresolved):
        raise ValueError(
            "Unresolved sample status labels; provide a trusted status_aliases map: "
            + ", ".join(map(str, unresolved))
        )

    observed = set(frame["canonical_status"])
    if observed == {"tumor", "normal-like"}:
        mode = "tumor_normal"; status_balance_applied = True; cm_classification_available = True
    elif observed == {"tumor"}:
        mode = "tumor_only"; status_balance_applied = False; cm_classification_available = False
    elif observed == {"normal-like"}:
        mode = "normal_only"; status_balance_applied = False; cm_classification_available = False
    else:
        raise ValueError(f"Unsupported Block 03 status set: {sorted(observed)}")

    frame[status_col] = frame["canonical_status"]
    counts = (frame.groupby("canonical_status", observed=True)[sample_col]
              .nunique().reset_index(name="n_samples"))
    counts["n_total_samples"] = int(frame[sample_col].nunique())
    counts["detected_mode"] = mode
    counts["status_source"] = status_source
    counts["status_balance_applied"] = status_balance_applied
    counts["cm_classification_available"] = cm_classification_available

    out_dir.mkdir(parents=True, exist_ok=True)
    counts.to_csv(out_dir / "cm_status_mode_detection.csv", index=False)
    with (out_dir / "cm_status_mode.json").open("w", encoding="utf-8") as handle:
        json.dump({
            "detected_mode": mode,
            "n_total_samples": int(frame[sample_col].nunique()),
            "status_source": status_source,
            "status_balance_applied": status_balance_applied,
            "cm_classification_available": cm_classification_available,
            "raw_to_canonical_status": dict(
                frame[["raw_status", "canonical_status"]]
                .drop_duplicates().itertuples(index=False, name=None)),
            "included_samples": frame[sample_col].astype(str).tolist(),
        }, handle, indent=2, ensure_ascii=False)
    return frame, mode


def main():
    if not ACT.exists():
        raise FileNotFoundError(f"Required input missing: {ACT}")
    sample_status = pd.read_csv(ACT)[["sample_id", "status"]]
    print(f"[preflight] read {len(sample_status)} rows from {ACT.name}", flush=True)
    print(f"[preflight] status counts:\n{sample_status['status'].value_counts()}", flush=True)
    frame, mode = detect_block03_status_mode(
        sample_status, OUT_DIR, sample_col="sample_id", status_col="status",
        status_source="activity_df_sample_by_CM.csv")
    print(f"[preflight] detected_mode={mode}", flush=True)
    print(f"[preflight] wrote {OUT_DIR}/cm_status_mode_detection.csv", flush=True)
    print(f"[preflight] wrote {OUT_DIR}/cm_status_mode.json", flush=True)


if __name__ == "__main__":
    main()