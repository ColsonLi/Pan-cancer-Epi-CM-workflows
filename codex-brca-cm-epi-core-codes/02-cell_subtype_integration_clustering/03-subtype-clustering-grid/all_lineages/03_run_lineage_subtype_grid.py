#!/usr/bin/env python3
"""CLI wrapper for the validated Module 02 complete subtype grid."""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path


ROOT = Path("/mnt/disk18t/lr_xcy/riku/brca_cm-epi_codex")
BASE_CODE = (
    ROOT
    / "epi-cm-core-workflow"
    / "codes"
    / "02-cell_subtype_integration_clustering"
    / "03-subtype-clustering-grid"
    / "epithelial"
    / "03_run_epithelial_subtype_grid.py"
)
CODE_PATH = Path(__file__).resolve()
CHOICES = [
    "t_cells",
    "myeloid",
    "b_cells",
    "plasma",
    "endothelial",
    "stromal",
    "perivascular",
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lineage", required=True, choices=CHOICES)
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location("module02_lineage_grid", BASE_CODE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load base grid implementation: {BASE_CODE}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.configure_lineage(args.lineage, code_path=CODE_PATH)
    module.main()


if __name__ == "__main__":
    main()
