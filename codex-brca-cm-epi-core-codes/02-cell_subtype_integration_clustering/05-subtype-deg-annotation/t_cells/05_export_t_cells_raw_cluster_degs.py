#!/usr/bin/env python3
"""Export full-length T-cell raw-cluster DEGs."""
import importlib.util
from pathlib import Path
BASE = Path(__file__).parents[1] / "epithelial" / "05_export_epithelial_raw_cluster_degs.py"
spec = importlib.util.spec_from_file_location("module02_raw_deg_base", BASE)
if spec is None or spec.loader is None: raise RuntimeError(f"Cannot load {BASE}")
module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
module.configure_lineage("t_cells", code_path=Path(__file__).resolve()); module.main()
