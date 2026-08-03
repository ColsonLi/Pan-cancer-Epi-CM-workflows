#!/usr/bin/env python3
"""Prepare all observed broad lineages and queue subtype GPU grids one lineage at a time."""

from __future__ import annotations

import csv
import subprocess
import sys
import time
from pathlib import Path

WORKDIR = Path("/mnt/disk18t/lr_xcy/riku/crc_val")
WORKFLOW_ROOT = WORKDIR / "epi-cm-core-workflow"
MAIN_PY = WORKDIR / "uv_envs/main/.venv/bin/python"
RAPIDS_PY = WORKDIR / "uv_envs/rapids/.venv/bin/python"
MODULE = "02-cell_subtype_integration_clustering"

LINEAGES = [
    ("epithelial", "Epithelial Cells"),
    ("t_cells", "T Cells"),
    ("plasma_cells", "Plasma Cells"),
    ("myeloid", "Myeloid Cells"),
    ("b_cells", "B Cells"),
    ("stromal", "Stromal Cells"),
    ("endothelial", "Endothelial Cells"),
    ("mast", "Mast Cells"),
    ("cycling_immune", "Cycling Immune Cells"),
    ("schwann", "Schwann Cells"),
]

SELECT_SCRIPT = WORKFLOW_ROOT / "codes" / MODULE / "01-lineage-selection" / "01_extract_lineage_from_qc.py"
HARMONY_SCRIPT = WORKFLOW_ROOT / "codes" / MODULE / "02-subtype-harmony" / "01_run_lineage_harmony.py"
GRID_SCRIPT = WORKFLOW_ROOT / "codes" / MODULE / "03-subtype-clustering-grid" / "01_run_lineage_gpu_grid.py"
RUN_TABLE_DIR = WORKFLOW_ROOT / "tables" / MODULE / "00-run"
RUN_LOG_DIR = RUN_TABLE_DIR / "logs"
MANIFEST = RUN_TABLE_DIR / "all_lineage_subtype_grid_run_manifest.csv"


def path_exists(step: str, lineage: str) -> bool:
    if step == "selection":
        return (WORKFLOW_ROOT / "h5ad" / MODULE / "01-lineage-selection" / lineage / f"adata_{lineage}_qc.h5ad").exists()
    if step == "harmony":
        return (WORKFLOW_ROOT / "h5ad" / MODULE / "02-subtype-harmony" / lineage / f"adata_{lineage}_harmony.h5ad").exists()
    if step == "grid":
        return (WORKFLOW_ROOT / "tables" / MODULE / "03-subtype-clustering-grid" / lineage / "clustering_grid_completion_check.csv").exists()
    raise ValueError(step)


def append_manifest(row: dict) -> None:
    RUN_TABLE_DIR.mkdir(parents=True, exist_ok=True)
    exists = MANIFEST.exists()
    with MANIFEST.open("a", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["time", "lineage", "target_label", "step", "status", "returncode", "log"])
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def run_step(cmd: list[str], lineage: str, target_label: str, step: str) -> None:
    RUN_LOG_DIR.mkdir(parents=True, exist_ok=True)
    log = RUN_LOG_DIR / f"{lineage}_{step}.log"
    append_manifest(
        {
            "time": time.strftime("%Y-%m-%d %H:%M:%S"),
            "lineage": lineage,
            "target_label": target_label,
            "step": step,
            "status": "running",
            "returncode": "",
            "log": str(log),
        }
    )
    with log.open("a") as fh:
        fh.write(f"\n[{time.strftime('%Y-%m-%d %H:%M:%S')}] START {' '.join(cmd)}\n")
        fh.flush()
        proc = subprocess.run(cmd, cwd=WORKDIR, stdout=fh, stderr=subprocess.STDOUT, text=True)
        fh.write(f"\n[{time.strftime('%Y-%m-%d %H:%M:%S')}] END returncode={proc.returncode}\n")
    status = "completed" if proc.returncode == 0 else "failed"
    append_manifest(
        {
            "time": time.strftime("%Y-%m-%d %H:%M:%S"),
            "lineage": lineage,
            "target_label": target_label,
            "step": step,
            "status": status,
            "returncode": proc.returncode,
            "log": str(log),
        }
    )
    if proc.returncode != 0:
        raise RuntimeError(f"{step} failed for {lineage}; see {log}")


def main() -> None:
    RUN_TABLE_DIR.mkdir(parents=True, exist_ok=True)
    RUN_LOG_DIR.mkdir(parents=True, exist_ok=True)
    with (RUN_TABLE_DIR / "readme.txt").open("w") as fh:
        fh.write("All observed broad lineages are prepared for subtype clustering.\n")
        fh.write("Selection and Harmony run as needed; GPU clustering grids are run one lineage at a time to avoid VRAM contention.\n")
        fh.write("Grid range for every lineage: pcs 10-50 step 5, n_neighbors 10-50 step 5, resolution 0.1-1.5 step 0.1.\n")
        fh.write("Candidate h5ad files are not saved during grid search.\n")

    for lineage, target_label in LINEAGES:
        if not path_exists("selection", lineage):
            run_step(
                [str(MAIN_PY), str(SELECT_SCRIPT), "--lineage", lineage, "--target-label", target_label],
                lineage,
                target_label,
                "selection",
            )
        else:
            run_step(
                [str(MAIN_PY), str(SELECT_SCRIPT), "--lineage", lineage, "--target-label", target_label, "--allow-existing"],
                lineage,
                target_label,
                "selection_validate",
            )

        if not path_exists("harmony", lineage):
            run_step([str(MAIN_PY), str(HARMONY_SCRIPT), "--lineage", lineage], lineage, target_label, "harmony")
        else:
            run_step([str(MAIN_PY), str(HARMONY_SCRIPT), "--lineage", lineage, "--allow-existing"], lineage, target_label, "harmony_validate")

        if not path_exists("grid", lineage):
            run_step([str(RAPIDS_PY), str(GRID_SCRIPT), "--lineage", lineage], lineage, target_label, "gpu_grid")
        else:
            run_step([str(RAPIDS_PY), str(GRID_SCRIPT), "--lineage", lineage, "--allow-existing"], lineage, target_label, "gpu_grid_validate")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
