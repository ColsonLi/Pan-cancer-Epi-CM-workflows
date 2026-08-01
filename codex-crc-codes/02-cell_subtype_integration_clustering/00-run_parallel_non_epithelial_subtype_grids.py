#!/usr/bin/env python3
"""Run non-epithelial subtype preparation in parallel and GPU grids with bounded concurrency."""

from __future__ import annotations

import csv
import subprocess
import time
from collections import deque
from pathlib import Path

WORKDIR = Path("/mnt/disk18t/lr_xcy/riku/crc_val")
WORKFLOW_ROOT = WORKDIR / "epi-cm-core-workflow"
MAIN_PY = WORKDIR / "uv_envs/main/.venv/bin/python"
RAPIDS_PY = WORKDIR / "uv_envs/rapids/.venv/bin/python"
MODULE = "02-cell_subtype_integration_clustering"

CPU_MAX_JOBS = 4
GPU_MAX_JOBS = 2

# Smaller lineages first so the user gets usable review grids sooner.
LINEAGES = [
    ("schwann", "Schwann Cells"),
    ("cycling_immune", "Cycling Immune Cells"),
    ("mast", "Mast Cells"),
    ("endothelial", "Endothelial Cells"),
    ("stromal", "Stromal Cells"),
    ("b_cells", "B Cells"),
    ("myeloid", "Myeloid Cells"),
    ("plasma_cells", "Plasma Cells"),
    ("t_cells", "T Cells"),
]

SELECT_SCRIPT = WORKFLOW_ROOT / "codes" / MODULE / "01-lineage-selection" / "01_extract_lineage_from_qc.py"
HARMONY_SCRIPT = WORKFLOW_ROOT / "codes" / MODULE / "02-subtype-harmony" / "01_run_lineage_harmony.py"
GRID_SCRIPT = WORKFLOW_ROOT / "codes" / MODULE / "03-subtype-clustering-grid" / "01_run_lineage_gpu_grid.py"
RUN_TABLE_DIR = WORKFLOW_ROOT / "tables" / MODULE / "00-run-parallel"
RUN_LOG_DIR = RUN_TABLE_DIR / "logs"
MANIFEST = RUN_TABLE_DIR / "parallel_non_epithelial_run_manifest.csv"


def timestamp() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def append_manifest(row: dict) -> None:
    RUN_TABLE_DIR.mkdir(parents=True, exist_ok=True)
    exists = MANIFEST.exists()
    with MANIFEST.open("a", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["time", "lineage", "target_label", "stage", "pid", "status", "returncode", "log"],
        )
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def selection_h5ad(lineage: str) -> Path:
    return WORKFLOW_ROOT / "h5ad" / MODULE / "01-lineage-selection" / lineage / f"adata_{lineage}_qc.h5ad"


def harmony_h5ad(lineage: str) -> Path:
    return WORKFLOW_ROOT / "h5ad" / MODULE / "02-subtype-harmony" / lineage / f"adata_{lineage}_harmony.h5ad"


def grid_done(lineage: str) -> bool:
    path = WORKFLOW_ROOT / "tables" / MODULE / "03-subtype-clustering-grid" / lineage / "clustering_grid_completion_check.csv"
    return path.exists()


def prep_cmd(lineage: str, target_label: str) -> list[str]:
    commands = []
    if selection_h5ad(lineage).exists():
        commands.append(
            f"{MAIN_PY} {SELECT_SCRIPT} --lineage {lineage} --target-label '{target_label}' --allow-existing"
        )
    else:
        commands.append(f"{MAIN_PY} {SELECT_SCRIPT} --lineage {lineage} --target-label '{target_label}'")
    if harmony_h5ad(lineage).exists():
        commands.append(f"{MAIN_PY} {HARMONY_SCRIPT} --lineage {lineage} --allow-existing")
    else:
        commands.append(f"{MAIN_PY} {HARMONY_SCRIPT} --lineage {lineage}")
    return ["bash", "-lc", " && ".join(commands)]


def grid_cmd(lineage: str) -> list[str]:
    if grid_done(lineage):
        return [str(RAPIDS_PY), str(GRID_SCRIPT), "--lineage", lineage, "--allow-existing"]
    return [str(RAPIDS_PY), str(GRID_SCRIPT), "--lineage", lineage]


def launch(lineage: str, target_label: str, stage: str, cmd: list[str]) -> tuple[subprocess.Popen, object, Path]:
    RUN_LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = RUN_LOG_DIR / f"{lineage}_{stage}.log"
    fh = log_path.open("a")
    fh.write(f"\n[{timestamp()}] START {' '.join(cmd)}\n")
    fh.flush()
    proc = subprocess.Popen(cmd, cwd=WORKDIR, stdout=fh, stderr=subprocess.STDOUT, text=True)
    append_manifest(
        {
            "time": timestamp(),
            "lineage": lineage,
            "target_label": target_label,
            "stage": stage,
            "pid": proc.pid,
            "status": "running",
            "returncode": "",
            "log": str(log_path),
        }
    )
    return proc, fh, log_path


def finish(lineage: str, target_label: str, stage: str, proc: subprocess.Popen, fh: object, log_path: Path) -> bool:
    rc = proc.returncode
    fh.write(f"\n[{timestamp()}] END returncode={rc}\n")
    fh.close()
    append_manifest(
        {
            "time": timestamp(),
            "lineage": lineage,
            "target_label": target_label,
            "stage": stage,
            "pid": proc.pid,
            "status": "completed" if rc == 0 else "failed",
            "returncode": rc,
            "log": str(log_path),
        }
    )
    return rc == 0


def main() -> None:
    RUN_TABLE_DIR.mkdir(parents=True, exist_ok=True)
    RUN_LOG_DIR.mkdir(parents=True, exist_ok=True)
    with (RUN_TABLE_DIR / "readme.txt").open("w") as fh:
        fh.write("Non-epithelial subtype workflow launched in parallel.\n")
        fh.write(f"CPU prep concurrency: {CPU_MAX_JOBS}\n")
        fh.write(f"GPU grid concurrency: {GPU_MAX_JOBS}\n")
        fh.write("Grid range: pcs 10-50 step 5; n_neighbors 10-50 step 5; resolution 0.1-1.5 step 0.1.\n")
        fh.write("Candidate h5ad files are not saved during grid search.\n")

    prep_queue = deque(LINEAGES)
    grid_queue: deque[tuple[str, str]] = deque()
    prep_running: dict[str, tuple[str, subprocess.Popen, object, Path]] = {}
    grid_running: dict[str, tuple[str, subprocess.Popen, object, Path]] = {}
    failed: set[str] = set()

    while prep_queue or prep_running or grid_queue or grid_running:
        while prep_queue and len(prep_running) < CPU_MAX_JOBS:
            lineage, target_label = prep_queue.popleft()
            proc, fh, log_path = launch(lineage, target_label, "prep", prep_cmd(lineage, target_label))
            prep_running[lineage] = (target_label, proc, fh, log_path)

        for lineage, (target_label, proc, fh, log_path) in list(prep_running.items()):
            if proc.poll() is not None:
                ok = finish(lineage, target_label, "prep", proc, fh, log_path)
                prep_running.pop(lineage)
                if ok:
                    grid_queue.append((lineage, target_label))
                else:
                    failed.add(lineage)

        while grid_queue and len(grid_running) < GPU_MAX_JOBS:
            lineage, target_label = grid_queue.popleft()
            if lineage in failed:
                continue
            proc, fh, log_path = launch(lineage, target_label, "gpu_grid", grid_cmd(lineage))
            grid_running[lineage] = (target_label, proc, fh, log_path)

        for lineage, (target_label, proc, fh, log_path) in list(grid_running.items()):
            if proc.poll() is not None:
                ok = finish(lineage, target_label, "gpu_grid", proc, fh, log_path)
                grid_running.pop(lineage)
                if not ok:
                    failed.add(lineage)

        time.sleep(10)

    if failed:
        raise RuntimeError("Failed lineages: " + ", ".join(sorted(failed)))


if __name__ == "__main__":
    main()
