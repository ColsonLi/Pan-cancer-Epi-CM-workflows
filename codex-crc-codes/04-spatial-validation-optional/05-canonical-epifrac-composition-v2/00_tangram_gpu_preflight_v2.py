#!/usr/bin/env python3
"""Canonical Tangram/CUDA preflight for CRC spatial validation."""

import importlib.metadata
import json
import platform
from pathlib import Path

import torch
import tangram as tg

ROOT = Path("/mnt/disk18t/lr_xcy/riku/crc_codex/epi-cm-core-workflow")
OUT = ROOT / "tables/04-spatial-validation-optional-sig_genes/05-canonical-epifrac-composition-v2/02-tangram-mapping/tangram_gpu_preflight.json"
OUT.parent.mkdir(parents=True, exist_ok=True)

if not torch.cuda.is_available():
    raise RuntimeError("CUDA is not available; canonical Tangram requires cuda:0 on this machine")
x = torch.tensor([1.0, 2.0], device="cuda:0")
result = {
    "python": platform.python_version(),
    "tangram_sc": importlib.metadata.version("tangram-sc"),
    "torch": torch.__version__,
    "cuda_available": True,
    "cuda_runtime": torch.version.cuda,
    "device": torch.cuda.get_device_name(0),
    "device_total_memory_bytes": torch.cuda.get_device_properties(0).total_memory,
    "minimal_cuda_tensor_sum": float(x.sum().item()),
    "canonical_device": "cuda:0",
    "canonical_mode": "cells",
    "canonical_num_epochs": 350,
    "canonical_learning_rate": 0.05,
    "canonical_random_state": 42,
    "status": "passed",
}
OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
print(json.dumps(result, indent=2))
