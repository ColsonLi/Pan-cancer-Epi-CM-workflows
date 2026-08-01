#!/usr/bin/env python3
"""Export full-length raw-cluster DEGs for one selected BRCA lineage."""

from __future__ import annotations

import importlib.metadata
import json
import re
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


ROOT = Path("/mnt/disk18t/lr_xcy/riku/brca_cm-epi_codex")
WORKFLOW = ROOT / "epi-cm-core-workflow"
BLOCK = "02-cell_subtype_integration_clustering"
METHOD = "t-test"
SEED = 42
REVIEW_DEPTH = 50

LINEAGE_CONFIG = {
    "epithelial": {"abbrev": "epi", "label": "Epithelial Cells"},
    "t_cells": {"abbrev": "t", "label": "T Cells"},
    "myeloid": {"abbrev": "mye", "label": "Myeloid Cells"},
    "b_cells": {"abbrev": "b", "label": "B Cells"},
    "plasma": {"abbrev": "plasma", "label": "Plasma Cells"},
    "endothelial": {"abbrev": "endo", "label": "Endothelial Cells"},
    "stromal": {"abbrev": "stromal", "label": "Stromal Cells"},
    "perivascular": {"abbrev": "pvl", "label": "Perivascular Cells"},
}

LINEAGE_SLUG = "epithelial"
CODE_PATH = Path(__file__).resolve()


def configure_lineage(lineage_slug: str, code_path: Path | None = None) -> None:
    global LINEAGE_SLUG, CODE_PATH
    if lineage_slug not in LINEAGE_CONFIG:
        raise ValueError(f"Unsupported lineage slug: {lineage_slug}")
    LINEAGE_SLUG = lineage_slug
    if code_path is not None:
        CODE_PATH = code_path


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", str(value)).strip("_")


def res_token(resolution: float) -> str:
    return f"{resolution:.1f}".replace(".", "p")


def paths_and_parameters() -> tuple[dict[str, Path], pd.Series]:
    config = LINEAGE_CONFIG[LINEAGE_SLUG]
    selected_table = (
        WORKFLOW
        / "tables"
        / BLOCK
        / "04-subtype-selected-clustering"
        / LINEAGE_SLUG
        / "selected_clustering.csv"
    )
    selected = pd.read_csv(selected_table).iloc[0]
    raw_key = str(selected["raw_cluster_key"])
    pcs = int(selected["n_pcs"])
    nn = int(selected["n_neighbors"])
    resolution = float(selected["resolution"])
    token = f"{raw_key}_pcs{pcs}_nn{nn}_res{res_token(resolution)}"
    table_dir = WORKFLOW / "tables" / BLOCK / "05-subtype-deg-annotation" / LINEAGE_SLUG
    paths = {
        "input_h5ad": (
            WORKFLOW
            / "h5ad"
            / BLOCK
            / "04-subtype-selected-clustering"
            / LINEAGE_SLUG
            / f"adata_{config['abbrev']}_selected_clustered.h5ad"
        ),
        "selected_table": selected_table,
        "table_dir": table_dir,
        "deg_dir": table_dir / f"degs_{token}",
        "audit": table_dir / "raw_cluster_deg_audit.csv",
        "major_counts": table_dir / "raw_cluster_by_published_major_counts.csv",
        "minor_counts": table_dir / "raw_cluster_by_published_minor_counts.csv",
        "minor_fractions": table_dir / "raw_cluster_by_published_minor_fractions.csv",
        "parameters": table_dir / "raw_cluster_deg_parameters.csv",
        "completion": table_dir / "raw_cluster_deg_completion.json",
        "readme": table_dir / "raw_cluster_deg_readme.txt",
        "versions": table_dir / "raw_cluster_deg_package_versions.txt",
    }
    return paths, selected


def validate_existing(paths: dict[str, Path], selected: pd.Series) -> bool:
    if not paths["completion"].exists():
        other = [
            paths["audit"], paths["major_counts"], paths["minor_counts"],
            paths["minor_fractions"], paths["parameters"], paths["readme"],
            paths["versions"],
        ]
        if paths["deg_dir"].exists() or any(path.exists() for path in other):
            raise FileExistsError(
                "Partial raw-cluster DEG outputs already exist; refusing overwrite."
            )
        return False
    completion = json.loads(paths["completion"].read_text())
    audit = pd.read_csv(paths["audit"])
    expected_groups = int(completion["n_raw_clusters"])
    files = sorted(paths["deg_dir"].glob("*_degs_*.csv"))
    valid = (
        completion.get("status") == "completed"
        and completion.get("method") == METHOD
        and completion.get("use_raw") is True
        and len(audit) == expected_groups
        and len(files) == expected_groups
        and audit["n_deg_rows"].eq(int(completion["rows_per_deg_csv"])).all()
        and audit["full_length"].all()
        and str(selected["raw_cluster_key"]) == completion.get("groupby")
    )
    if not valid:
        raise FileExistsError("Existing raw-cluster DEG outputs failed validation.")
    return True


def main() -> None:
    started = time.time()
    paths, selected = paths_and_parameters()
    if validate_existing(paths, selected):
        print(
            json.dumps(
                {
                    "lineage": LINEAGE_SLUG,
                    "status": "valid_existing_raw_cluster_degs_reused",
                    "deg_dir": str(paths["deg_dir"]),
                },
                indent=2,
            )
        )
        return

    raw_key = str(selected["raw_cluster_key"])
    adata = sc.read_h5ad(paths["input_h5ad"])
    if adata.raw is None:
        raise ValueError("Raw-cluster DEGs require adata.raw normalized/log expression.")
    if raw_key not in adata.obs or adata.obs[raw_key].isna().any():
        raise ValueError(f"Missing or incomplete selected raw cluster column: {raw_key}")
    if "cell_subtype" in adata.obs:
        raise ValueError("Raw-cluster DEG export must precede subtype annotation.")
    if not adata.obs["leiden_coarse"].astype(str).eq(
        str(LINEAGE_CONFIG[LINEAGE_SLUG]["label"])
    ).all():
        raise ValueError("Off-lineage cells found in selected clustered input.")
    groups = adata.obs[raw_key].cat.categories.astype(str).tolist()
    if len(groups) != int(selected["cluster_count"]):
        raise ValueError(
            f"Selected cluster count mismatch: groups={len(groups)}, "
            f"parameters={selected['cluster_count']}"
        )

    print(
        f"[{LINEAGE_SLUG} raw DEG] cells={adata.n_obs} raw_genes={adata.raw.n_vars} "
        f"groups={len(groups)} method={METHOD} use_raw=True",
        flush=True,
    )
    sc.tl.rank_genes_groups(
        adata,
        groupby=raw_key,
        method=METHOD,
        use_raw=True,
        n_genes=adata.raw.n_vars,
        key_added="rank_genes_groups_raw_cluster",
    )
    result = adata.uns["rank_genes_groups_raw_cluster"]
    required = ["names", "scores", "logfoldchanges", "pvals", "pvals_adj"]
    missing = [key for key in required if key not in result]
    if missing:
        raise ValueError(f"rank_genes_groups result lacks fields: {missing}")

    paths["deg_dir"].mkdir(parents=True, exist_ok=True)
    audit_rows: list[dict[str, object]] = []
    pcs = int(selected["n_pcs"])
    nn = int(selected["n_neighbors"])
    resolution = float(selected["resolution"])
    suffix = f"{raw_key}_pcs{pcs}_nn{nn}_res{res_token(resolution)}"
    for group in groups:
        frame = pd.DataFrame(
            {
                "gene": np.asarray(result["names"][group]).astype(str),
                "score": np.asarray(result["scores"][group], dtype=float),
                "logfoldchanges": np.asarray(
                    result["logfoldchanges"][group], dtype=float
                ),
                "pvals": np.asarray(result["pvals"][group], dtype=float),
                "pvals_adj": np.asarray(result["pvals_adj"][group], dtype=float),
            }
        )
        if len(frame) != adata.raw.n_vars:
            raise ValueError(f"DEG table for cluster {group} is not full length.")
        if frame["gene"].isna().any() or frame["gene"].duplicated().any():
            raise ValueError(f"DEG genes for cluster {group} are missing or duplicated.")
        out = paths["deg_dir"] / f"{safe_name(group)}_degs_{suffix}.csv"
        frame.to_csv(out, index=False)
        audit_rows.append(
            {
                "lineage": LINEAGE_SLUG,
                "raw_cluster": group,
                "n_cells": int((adata.obs[raw_key].astype(str) == group).sum()),
                "n_deg_rows": len(frame),
                "full_length": True,
                "use_raw": True,
                "method": METHOD,
                "review_depth_initial": REVIEW_DEPTH,
                "deg_csv": str(out),
            }
        )
        print(
            f"[{LINEAGE_SLUG} raw DEG] wrote cluster={group} rows={len(frame)}",
            flush=True,
        )

    paths["table_dir"].mkdir(parents=True, exist_ok=True)
    audit = pd.DataFrame(audit_rows)
    audit.to_csv(paths["audit"], index=False)
    if "celltype_major" in adata.obs:
        pd.crosstab(
            adata.obs[raw_key].astype(str), adata.obs["celltype_major"].astype(str)
        ).to_csv(paths["major_counts"])
    else:
        pd.DataFrame(index=groups).to_csv(paths["major_counts"])
    if "celltype_minor" in adata.obs:
        minor = pd.crosstab(
            adata.obs[raw_key].astype(str), adata.obs["celltype_minor"].astype(str)
        )
        minor.to_csv(paths["minor_counts"])
        minor.div(minor.sum(axis=1), axis=0).to_csv(paths["minor_fractions"])
    else:
        pd.DataFrame(index=groups).to_csv(paths["minor_counts"])
        pd.DataFrame(index=groups).to_csv(paths["minor_fractions"])

    parameters = {
        "lineage": LINEAGE_SLUG,
        "input_h5ad": str(paths["input_h5ad"]),
        "selected_parameters_table": str(paths["selected_table"]),
        "groupby": raw_key,
        "n_cells": int(adata.n_obs),
        "n_raw_genes": int(adata.raw.n_vars),
        "n_raw_clusters": len(groups),
        "n_pcs": pcs,
        "n_neighbors": nn,
        "resolution": resolution,
        "method": METHOD,
        "use_raw": True,
        "n_genes_exported_per_cluster": int(adata.raw.n_vars),
        "combined_deg_csv_saved": False,
        "top_only_deg_csv_saved": False,
        "review_depth_initial": REVIEW_DEPTH,
        "published_labels_role": "independent_validation_only",
        "code_file": str(CODE_PATH),
        "seed": SEED,
    }
    pd.DataFrame([parameters]).to_csv(paths["parameters"], index=False)
    versions = {
        "python": sys.version.split()[0],
        "scanpy": package_version("scanpy"),
        "anndata": package_version("anndata"),
        "code": str(CODE_PATH),
        "input_h5ad": str(paths["input_h5ad"]),
        "groupby": raw_key,
        "method": METHOD,
        "use_raw": "True",
        "review_depth_initial": str(REVIEW_DEPTH),
        "seed": str(SEED),
    }
    paths["versions"].write_text(
        "\n".join(f"{key}={value}" for key, value in versions.items()) + "\n",
        encoding="utf-8",
    )
    paths["readme"].write_text(
        f"""BRCA {LINEAGE_CONFIG[LINEAGE_SLUG]['label']} raw-cluster DEGs

Input: {paths['input_h5ad']}
Groupby: {raw_key}
Method: {METHOD}
Expression: adata.raw normalized/log expression via use_raw=True
Output: {paths['deg_dir']}

One full-length DEG CSV with {adata.raw.n_vars} rows is written per raw cluster.
No combined or top-only DEG CSV is produced. Annotation review starts from the
first {REVIEW_DEPTH} rows of each saved full table. Published author labels are
saved only as independent validation crosstabs and do not replace DEG/marker
review.
""",
        encoding="utf-8",
    )
    completion = {
        "lineage": LINEAGE_SLUG,
        "status": "completed",
        "groupby": raw_key,
        "n_cells": int(adata.n_obs),
        "n_raw_genes": int(adata.raw.n_vars),
        "n_raw_clusters": len(groups),
        "n_full_length_deg_csvs": len(audit_rows),
        "rows_per_deg_csv": int(adata.raw.n_vars),
        "method": METHOD,
        "use_raw": True,
        "candidate_annotation_ready": True,
        "final_mapping_confirmed": False,
        "elapsed_seconds": time.time() - started,
    }
    paths["completion"].write_text(
        json.dumps(completion, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(completion, indent=2), flush=True)


if __name__ == "__main__":
    main()
