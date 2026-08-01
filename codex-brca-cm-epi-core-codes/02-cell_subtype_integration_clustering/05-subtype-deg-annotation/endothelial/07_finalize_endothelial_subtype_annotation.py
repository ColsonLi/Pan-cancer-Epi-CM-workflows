#!/usr/bin/env python3
"""Finalize endothelial subtype annotation."""
import importlib.util
from pathlib import Path
BASE=Path(__file__).parents[1]/"epithelial"/"07_finalize_epithelial_subtype_annotation.py"
spec=importlib.util.spec_from_file_location("module02_final_annotation_base",BASE)
if spec is None or spec.loader is None: raise RuntimeError(f"Cannot load {BASE}")
module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
module.configure_lineage("endothelial",code_path=Path(__file__).resolve()); module.main()
