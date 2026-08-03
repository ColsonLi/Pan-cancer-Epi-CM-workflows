#!/usr/bin/env python3
"""Run selected clustering for all user-selected lineages."""

from __future__ import annotations

import subprocess
from pathlib import Path

WORKDIR = Path("/mnt/disk18t/lr_xcy/riku/crc_val")
RAPIDS_PY = WORKDIR / "uv_envs/rapids/.venv/bin/python"
SCRIPT = (
    WORKDIR
    / "epi-cm-core-workflow/codes/02-cell_subtype_integration_clustering/04-subtype-selected-clustering/01_run_selected_lineage_clustering.py"
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
        cmd = [str(RAPIDS_PY), str(SCRIPT), "--lineage", lineage, "--allow-existing"]
        print("RUN", " ".join(cmd), flush=True)
        subprocess.run(cmd, cwd=WORKDIR, check=True)


if __name__ == "__main__":
    main()
