#!/usr/bin/env python3
"""Run DEG-based subtype annotation for all selected lineages."""

from __future__ import annotations

import subprocess
from pathlib import Path

WORKDIR = Path("/mnt/disk18t/lr_xcy/riku/crc_val")
MAIN_PY = WORKDIR / "uv_envs/main/.venv/bin/python"
SCRIPT = (
    WORKDIR
    / "epi-cm-core-workflow/codes/02-cell_subtype_integration_clustering/05-subtype-deg-annotation/01_auto_deg_annotate_lineage.py"
)
LINEAGES = [
    "b_cells",
    "cycling_immune",
    "endothelial",
    "epithelial",
    "mast",
    "myeloid",
    "plasma_cells",
    "schwann",
    "stromal",
    "t_cells",
]


def main() -> None:
    for lineage in LINEAGES:
        cmd = [str(MAIN_PY), str(SCRIPT), "--lineage", lineage, "--allow-existing"]
        print("RUN", " ".join(cmd), flush=True)
        subprocess.run(cmd, cwd=WORKDIR, check=True)


if __name__ == "__main__":
    main()
