#!/usr/bin/env python3
"""Prove that a completed grid resumes without redrawing valid figures."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path("/mnt/disk18t/lr_xcy/riku/brca_cm-epi_codex")
WORKFLOW = ROOT / "epi-cm-core-workflow"
BLOCK = "01-celltype_integration_clustering"
GRID = "05-clustering-parameter-search"
FIGURE_ROOT = WORKFLOW / "figures" / BLOCK / GRID
TABLE_ROOT = WORKFLOW / "tables" / BLOCK / GRID
GRID_SCRIPT = WORKFLOW / "codes" / BLOCK / GRID / "04_run_full_grid.py"
OUTPUT = TABLE_ROOT / "resume_skip_audit.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def snapshot() -> dict[str, tuple[int, str]]:
    figures = sorted(FIGURE_ROOT.rglob("umap_leiden_grid.pdf"))
    if len(figures) != 81:
        raise AssertionError(f"Expected 81 figures before resume, found {len(figures)}")
    return {
        str(path.relative_to(FIGURE_ROOT)): (path.stat().st_mtime_ns, sha256(path))
        for path in figures
    }


def main() -> None:
    before = snapshot()
    completed = subprocess.run(
        [sys.executable, "-u", str(GRID_SCRIPT)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    after = snapshot()
    if before != after:
        changed = sorted(key for key in before if before.get(key) != after.get(key))
        raise AssertionError(f"Resume changed existing figures: {changed}")
    if "restored_complete_graphs=81" not in completed.stdout:
        raise AssertionError("Resume did not report 81 restored complete graphs.")
    skipped_messages = completed.stdout.count("existing valid outputs skipped:")
    if skipped_messages != 81:
        raise AssertionError(
            f"Resume reported {skipped_messages} skipped graphs rather than 81."
        )
    summary = {
        "pass": True,
        "restored_complete_graphs": 81,
        "existing_valid_graphs_skipped": skipped_messages,
        "figures_before": len(before),
        "figures_after": len(after),
        "figure_hashes_unchanged": True,
        "figure_mtimes_unchanged": True,
        "figures_redrawn": 0,
    }
    OUTPUT.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
