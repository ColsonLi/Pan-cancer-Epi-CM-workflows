#!/usr/bin/env python3
"""Build the endothelial candidate subtype table."""
import importlib.util
from pathlib import Path
BASE=Path(__file__).parents[1]/"epithelial"/"06_build_epithelial_candidate_annotation.py"
spec=importlib.util.spec_from_file_location("module02_candidate_base",BASE)
if spec is None or spec.loader is None: raise RuntimeError(f"Cannot load {BASE}")
module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
module.configure_lineage("endothelial",code_path=Path(__file__).resolve()); module.main()
