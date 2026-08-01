#!/usr/bin/env python3
"""Run the user-selected B-cell subtype clustering."""

import importlib.util
from pathlib import Path

BASE = Path(__file__).parents[1] / "epithelial" / "04_run_epithelial_selected_clustering.py"
spec = importlib.util.spec_from_file_location("module02_selected_base", BASE)
if spec is None or spec.loader is None:
    raise RuntimeError(f"Cannot load selected-clustering implementation: {BASE}")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
module.configure_lineage("b_cells", code_path=Path(__file__).resolve())
module.main()
