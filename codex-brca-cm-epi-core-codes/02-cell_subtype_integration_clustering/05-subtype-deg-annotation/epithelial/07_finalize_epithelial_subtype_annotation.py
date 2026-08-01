#!/usr/bin/env python3
"""Finalize one BRCA lineage annotation from the agent-reviewed mapping."""

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
AUTHORIZATION_TEXT = "接下来完成所有内容"
AUTHORIZATION_INTERPRETATION = (
    "User authorized completion of all remaining workflow content; recorded as "
    "agent-led final subtype labeling authorization."
)

LINEAGE_CONFIG = {
    "epithelial": {"abbrev": "epi", "label": "Epithelial Cells", "prefix": "Epi"},
    "t_cells": {"abbrev": "t", "label": "T Cells", "prefix": "T"},
    "myeloid": {"abbrev": "mye", "label": "Myeloid Cells", "prefix": "Mye"},
    "b_cells": {"abbrev": "b", "label": "B Cells", "prefix": "B"},
    "plasma": {"abbrev": "plasma", "label": "Plasma Cells", "prefix": "Plasma"},
    "endothelial": {"abbrev": "endo", "label": "Endothelial Cells", "prefix": "Endo"},
    "stromal": {"abbrev": "stromal", "label": "Stromal Cells", "prefix": "Stromal"},
    "perivascular": {"abbrev": "pvl", "label": "Perivascular Cells", "prefix": "PVL"},
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


def res_token(resolution: float) -> str:
    return f"{resolution:.1f}".replace(".", "p")


def safe_group(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_")


def paths_and_parameters() -> tuple[dict[str, Path], pd.Series]:
    config = LINEAGE_CONFIG[LINEAGE_SLUG]
    selected_table = (
        WORKFLOW / "tables" / BLOCK / "04-subtype-selected-clustering"
        / LINEAGE_SLUG / "selected_clustering.csv"
    )
    selected = pd.read_csv(selected_table).iloc[0]
    pcs = int(selected["n_pcs"])
    nn = int(selected["n_neighbors"])
    resolution = float(selected["resolution"])
    token = f"cell_subtype_pcs{pcs}_nn{nn}_res{res_token(resolution)}"
    table_dir = WORKFLOW / "tables" / BLOCK / "05-subtype-deg-annotation" / LINEAGE_SLUG
    figure_dir = WORKFLOW / "figures" / BLOCK / "05-subtype-deg-annotation" / LINEAGE_SLUG
    h5ad_dir = WORKFLOW / "h5ad" / BLOCK / "05-subtype-deg-annotation" / LINEAGE_SLUG
    paths = {
        "input_h5ad": (
            WORKFLOW / "h5ad" / BLOCK / "04-subtype-selected-clustering"
            / LINEAGE_SLUG / f"adata_{config['abbrev']}_selected_clustered.h5ad"
        ),
        "selected_table": selected_table,
        "table_dir": table_dir,
        "figure_dir": figure_dir,
        "h5ad_dir": h5ad_dir,
        "candidate": table_dir / "candidate_subtype_annotation.csv",
        "agent_reviewed_candidate": table_dir / "candidate_subtype_annotation_agent_reviewed.csv",
        "mapping": table_dir / "final_raw_cluster_to_subtype_mapping.csv",
        "authorization": table_dir / "final_label_authorization.csv",
        "post_deg_dir": table_dir / f"degs_{token}",
        "post_deg_audit": table_dir / "post_annotation_cell_subtype_deg_audit.csv",
        "subtype_counts": table_dir / "subtype_counts.csv",
        "sample_counts": table_dir / "subtype_by_sample_counts.csv",
        "sample_fractions": table_dir / "subtype_by_sample_fractions.csv",
        "parameters": table_dir / "final_subtype_annotation_parameters.csv",
        "completion": table_dir / "final_subtype_annotation_completion.json",
        "readme": table_dir / "readme.txt",
        "versions": table_dir / "package_versions.txt",
        "umap": figure_dir / "umap_cell_subtype.pdf",
        "dotplot": figure_dir / "dotplot__cell_subtype_markers.pdf",
        "output_h5ad": h5ad_dir / f"adata_{config['abbrev']}.h5ad",
    }
    return paths, selected


def pdf_is_readable(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    import subprocess
    result = subprocess.run(["pdfinfo", str(path)], capture_output=True, text=True)
    return result.returncode == 0 and "Pages:" in result.stdout


def validate_existing(paths: dict[str, Path], selected: pd.Series) -> bool:
    required = [
        paths["agent_reviewed_candidate"], paths["mapping"], paths["authorization"],
        paths["post_deg_audit"], paths["subtype_counts"], paths["sample_counts"],
        paths["sample_fractions"], paths["parameters"], paths["completion"],
        paths["readme"], paths["versions"], paths["umap"], paths["dotplot"],
        paths["output_h5ad"],
    ]
    present = [path.exists() for path in required]
    post_deg_exists = paths["post_deg_dir"].exists()
    if not any(present) and not post_deg_exists:
        return False
    if not all(present) or not post_deg_exists:
        return False
    completion = json.loads(paths["completion"].read_text())
    mapping = pd.read_csv(paths["mapping"], dtype={"cluster": str})
    audit = pd.read_csv(paths["post_deg_audit"])
    files = list(paths["post_deg_dir"].glob("*_degs_*.csv"))
    raw_key = str(selected["raw_cluster_key"])
    saved = ad.read_h5ad(paths["output_h5ad"], backed="r")
    try:
        valid = (
            completion.get("status") == "completed"
            and completion.get("agent_led_final_labeling_authorized") is True
            and len(mapping) == int(selected["cluster_count"])
            and mapping["cluster"].is_unique
            and mapping["cell_subtype"].is_unique
            and len(audit) == len(mapping)
            and len(files) == len(mapping)
            and raw_key in saved.obs
            and "cell_subtype" in saved.obs
            and "functional_state" in saved.obs
            and int(saved.obs["cell_subtype"].nunique()) == len(mapping)
            and saved.raw is not None
            and pdf_is_readable(paths["umap"])
            and pdf_is_readable(paths["dotplot"])
        )
    finally:
        saved.file.close()
    if not valid:
        raise FileExistsError("Existing final subtype annotation outputs failed validation.")
    return True


def correct_duplicate_winner_rationale(candidate: pd.DataFrame) -> pd.DataFrame:
    candidate = candidate.copy()
    for idx, row in candidate.iterrows():
        same_gene = (
            str(row["most_significant_deg_gene"])
            == str(row["top_positive_deg_gene_before_duplicate_resolution"])
        )
        rationale = str(row["gene_selection_rationale"])
        if same_gene and "was shared by clusters" in rationale and "next-ranked" in rationale:
            gene = str(row["most_significant_deg_gene"])
            duplicate_clause = rationale.split("Top positive DEG", 1)[1]
            candidate.loc[idx, "gene_selection_rationale"] = (
                f"{gene} is the first positive-logfoldchange gene in the full t-test DEG "
                f"table sorted by Scanpy rank_genes_groups evidence. Top positive DEG"
                f"{duplicate_clause}"
            )
    return candidate


def main() -> None:
    started = time.time()
    config = LINEAGE_CONFIG[LINEAGE_SLUG]
    paths, selected = paths_and_parameters()
    if validate_existing(paths, selected):
        print(json.dumps({"lineage": LINEAGE_SLUG, "status": "valid_existing_final_annotation_reused", "output_h5ad": str(paths["output_h5ad"])}, indent=2))
        return

    candidate = pd.read_csv(paths["candidate"], dtype={"cluster": str})
    candidate = correct_duplicate_winner_rationale(candidate)
    if (
        not candidate["cluster"].is_unique
        or not candidate["proposed_cell_subtype"].is_unique
        or len(candidate) != int(selected["cluster_count"])
    ):
        raise ValueError("Candidate annotation is not a complete one-to-one mapping.")
    if not candidate["selected_prefix"].eq(config["prefix"]).all():
        raise ValueError("Candidate prefix does not match the recorded lineage prefix.")
    candidate["agent_led_final_labeling_authorized"] = True
    candidate["needs_user_confirmation"] = False
    if paths["agent_reviewed_candidate"].exists():
        existing_candidate = pd.read_csv(
            paths["agent_reviewed_candidate"], dtype={"cluster": str}
        )
        compare_columns = [
            "cluster", "proposed_cell_subtype", "functional_state",
            "most_significant_deg_gene", "gene_selection_rationale",
        ]
        if not existing_candidate[compare_columns].equals(candidate[compare_columns]):
            raise FileExistsError("Existing agent-reviewed candidate table differs from this run.")
    else:
        candidate.to_csv(paths["agent_reviewed_candidate"], index=False)

    mapping = candidate.rename(columns={"proposed_cell_subtype": "cell_subtype"}).copy()
    mapping["annotation_note"] = np.where(
        mapping["marker-supported lineage"].astype(str).str.lower().isin(
            [LINEAGE_SLUG.lower(), config["prefix"].lower(), config["label"].lower(),
             "epithelial", "t", "myeloid", "b", "plasma", "endothelial",
             "stromal", "perivascular"]
        ),
        "Retained in current broad lineage with standard prefix.",
        "Marker evidence is alternate-lineage-like; cells retained in current lineage as required.",
    )
    required_mapping = ["cluster", "cell_subtype", "functional_state", "gene_selection_rationale"]
    if mapping[required_mapping].isna().any().any():
        raise ValueError("Final mapping contains missing required values.")
    if paths["mapping"].exists():
        existing_mapping = pd.read_csv(paths["mapping"], dtype={"cluster": str})
        compare_columns = [
            "cluster", "cell_subtype", "functional_state",
            "most_significant_deg_gene", "gene_selection_rationale",
            "annotation_note",
        ]
        if not existing_mapping[compare_columns].equals(mapping[compare_columns]):
            raise FileExistsError("Existing final mapping differs from this run.")
    else:
        mapping.to_csv(paths["mapping"], index=False)
    authorization = pd.DataFrame(
        [
            {
                "lineage": LINEAGE_SLUG,
                "authorization_source": "user_message",
                "authorization_text": AUTHORIZATION_TEXT,
                "authorization_interpretation": AUTHORIZATION_INTERPRETATION,
                "agent_led_final_labeling_authorized": True,
                "selected_prefix": config["prefix"],
                "candidate_table": str(paths["agent_reviewed_candidate"]),
                "final_mapping": str(paths["mapping"]),
            }
        ]
    )
    if paths["authorization"].exists():
        existing_authorization = pd.read_csv(paths["authorization"])
        compare_columns = [
            "lineage", "authorization_source", "authorization_text",
            "authorization_interpretation", "selected_prefix",
        ]
        if not existing_authorization[compare_columns].equals(
            authorization[compare_columns]
        ):
            raise FileExistsError("Existing authorization record differs from this run.")
    else:
        authorization.to_csv(paths["authorization"], index=False)

    adata = sc.read_h5ad(paths["input_h5ad"])
    raw_key = str(selected["raw_cluster_key"])
    if adata.raw is None:
        raise ValueError("Final annotation requires adata.raw normalized/log expression.")
    if raw_key not in adata.obs or adata.obs[raw_key].isna().any():
        raise ValueError(f"Missing selected raw cluster key: {raw_key}")
    raw_groups = set(adata.obs[raw_key].astype(str).unique())
    if raw_groups != set(mapping["cluster"]):
        raise ValueError("Final mapping does not cover every selected raw cluster exactly once.")

    mapping_idx = mapping.set_index("cluster")
    subtype_order = mapping["cell_subtype"].astype(str).tolist()
    subtype = adata.obs[raw_key].astype(str).map(mapping_idx["cell_subtype"])
    state = adata.obs[raw_key].astype(str).map(mapping_idx["functional_state"])
    confidence = adata.obs[raw_key].astype(str).map(mapping_idx["confidence"])
    evidence_lineage = adata.obs[raw_key].astype(str).map(mapping_idx["marker-supported lineage"])
    prefix_rationale = adata.obs[raw_key].astype(str).map(mapping_idx["prefix rationale"])
    note = adata.obs[raw_key].astype(str).map(mapping_idx["annotation_note"])
    if any(series.isna().any() for series in [subtype, state, confidence, evidence_lineage, prefix_rationale, note]):
        raise ValueError("Mapping produced missing annotation values.")
    adata.obs["cell_subtype"] = pd.Categorical(
        subtype, categories=subtype_order, ordered=True
    )
    adata.obs["functional_state"] = pd.Categorical(state)
    adata.obs["annotation_confidence"] = pd.Categorical(confidence)
    adata.obs["marker_evidence_lineage"] = pd.Categorical(evidence_lineage)
    adata.obs["prefix_rationale"] = prefix_rationale.astype(str)
    adata.obs["annotation_note"] = note.astype(str)
    adata.obs["annotation_source"] = "Module02 agent-authorized DEG/marker review"
    if "cell_type" not in adata.obs or not adata.obs["cell_type"].astype(str).eq(config["label"]).all():
        raise ValueError("Trusted selected input cell_type is absent or inconsistent with broad lineage.")

    print(
        f"[{LINEAGE_SLUG} final annotation] cells={adata.n_obs} "
        f"subtypes={adata.obs['cell_subtype'].nunique()} post-DEG use_raw=True",
        flush=True,
    )
    sc.tl.rank_genes_groups(
        adata,
        groupby="cell_subtype",
        method=METHOD,
        use_raw=True,
        n_genes=adata.raw.n_vars,
        key_added="rank_genes_groups_cell_subtype",
    )
    result = adata.uns["rank_genes_groups_cell_subtype"]
    required = ["names", "scores", "logfoldchanges", "pvals", "pvals_adj"]
    missing = [key for key in required if key not in result]
    if missing:
        raise ValueError(f"Post-annotation DEG result lacks fields: {missing}")

    paths["post_deg_dir"].mkdir(parents=True, exist_ok=True)
    pcs = int(selected["n_pcs"])
    nn = int(selected["n_neighbors"])
    resolution = float(selected["resolution"])
    suffix = f"cell_subtype_pcs{pcs}_nn{nn}_res{res_token(resolution)}"
    audit_rows: list[dict[str, object]] = []
    for group in subtype_order:
        frame = pd.DataFrame(
            {
                "gene": np.asarray(result["names"][group]).astype(str),
                "score": np.asarray(result["scores"][group], dtype=float),
                "logfoldchanges": np.asarray(result["logfoldchanges"][group], dtype=float),
                "pvals": np.asarray(result["pvals"][group], dtype=float),
                "pvals_adj": np.asarray(result["pvals_adj"][group], dtype=float),
            }
        )
        if len(frame) != adata.raw.n_vars or frame["gene"].duplicated().any():
            raise ValueError(f"Post-annotation DEG table is not full length: {group}")
        out = paths["post_deg_dir"] / f"{safe_group(group)}_degs_{suffix}.csv"
        if out.exists():
            existing_rows = sum(1 for _ in out.open(encoding="utf-8")) - 1
            if existing_rows != adata.raw.n_vars:
                raise FileExistsError(
                    f"Existing post-annotation DEG is not full length: {out}"
                )
        else:
            frame.to_csv(out, index=False)
        audit_rows.append(
            {
                "cell_subtype": group,
                "n_cells": int((adata.obs["cell_subtype"].astype(str) == group).sum()),
                "n_deg_rows": len(frame),
                "full_length": True,
                "method": METHOD,
                "use_raw": True,
                "review_depth": REVIEW_DEPTH,
                "deg_csv": str(out),
            }
        )
    post_deg_audit = pd.DataFrame(audit_rows)
    if paths["post_deg_audit"].exists():
        existing_audit = pd.read_csv(paths["post_deg_audit"])
        if not existing_audit.equals(post_deg_audit):
            raise FileExistsError(
                "Existing post-annotation DEG audit differs from this run."
            )
    else:
        post_deg_audit.to_csv(paths["post_deg_audit"], index=False)

    counts = adata.obs["cell_subtype"].value_counts(sort=False).rename_axis("cell_subtype").reset_index(name="n_cells")
    counts["fraction_of_lineage"] = counts["n_cells"] / adata.n_obs
    if not paths["subtype_counts"].exists():
        counts.to_csv(paths["subtype_counts"], index=False)
    sample_counts = pd.crosstab(
        adata.obs["sample"].astype(str), adata.obs["cell_subtype"].astype(str)
    ).reindex(columns=subtype_order, fill_value=0)
    if not paths["sample_counts"].exists():
        sample_counts.to_csv(paths["sample_counts"])
    if not paths["sample_fractions"].exists():
        sample_counts.div(
            sample_counts.sum(axis=1).replace(0, np.nan), axis=0
        ).fillna(0).to_csv(paths["sample_fractions"])

    paths["figure_dir"].mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = paths["figure_dir"]
    sc.settings.file_format_figs = "pdf"
    sc.set_figure_params(figsize=(3, 3), dpi=150, fontsize=8)
    if not paths["umap"].exists():
        sc.pl.umap(adata, color="cell_subtype", show=False, save="_cell_subtype.pdf")
    if not pdf_is_readable(paths["umap"]):
        raise FileNotFoundError(f"Unreadable subtype UMAP: {paths['umap']}")
    marker_genes = mapping["most_significant_deg_gene"].astype(str).tolist()
    missing_markers = sorted(set(marker_genes) - set(adata.raw.var_names.astype(str)))
    if missing_markers:
        raise ValueError(f"Subtype naming markers absent from adata.raw: {missing_markers}")
    if not paths["dotplot"].exists():
        sc.pl.dotplot(
            adata,
            var_names=marker_genes,
            groupby="cell_subtype",
            use_raw=True,
            standard_scale="var",
            show=False,
            save="_cell_subtype_markers.pdf",
        )
    if not pdf_is_readable(paths["dotplot"]):
        raise FileNotFoundError(f"Unreadable subtype marker dotplot: {paths['dotplot']}")

    adata.uns["subtype_annotation"] = {
        "lineage": LINEAGE_SLUG,
        "raw_cluster_key": raw_key,
        "mapping_csv": str(paths["mapping"]),
        "agent_led_final_labeling_authorized": True,
        "selected_prefix": config["prefix"],
        "method": METHOD,
        "use_raw": True,
        "review_depth": REVIEW_DEPTH,
    }
    paths["h5ad_dir"].mkdir(parents=True, exist_ok=True)
    if not paths["output_h5ad"].exists():
        adata.write_h5ad(paths["output_h5ad"], compression="gzip")
    saved = ad.read_h5ad(paths["output_h5ad"], backed="r")
    try:
        if (
            "cell_subtype" not in saved.obs
            or "functional_state" not in saved.obs
            or int(saved.obs["cell_subtype"].nunique()) != len(mapping)
            or saved.raw is None
        ):
            raise ValueError(f"Saved annotated h5ad failed validation: {paths['output_h5ad']}")
    finally:
        saved.file.close()

    parameters = {
        "lineage": LINEAGE_SLUG,
        "input_selected_h5ad": str(paths["input_h5ad"]),
        "output_annotated_h5ad": str(paths["output_h5ad"]),
        "raw_cluster_key": raw_key,
        "n_cells": int(adata.n_obs),
        "n_raw_genes": int(adata.raw.n_vars),
        "n_raw_clusters": len(mapping),
        "n_cell_subtypes": int(adata.obs["cell_subtype"].nunique()),
        "n_pcs": pcs,
        "n_neighbors": nn,
        "resolution": resolution,
        "method": METHOD,
        "use_raw": True,
        "review_depth": REVIEW_DEPTH,
        "selected_prefix": config["prefix"],
        "agent_led_final_labeling_authorized": True,
        "candidate_table": str(paths["agent_reviewed_candidate"]),
        "final_mapping": str(paths["mapping"]),
        "post_annotation_deg_dir": str(paths["post_deg_dir"]),
        "umap": str(paths["umap"]),
        "dotplot": str(paths["dotplot"]),
        "code_file": str(CODE_PATH),
        "seed": SEED,
    }
    pd.DataFrame([parameters]).to_csv(paths["parameters"], index=False)
    versions = {
        "python": sys.version.split()[0],
        "scanpy": package_version("scanpy"),
        "anndata": package_version("anndata"),
        "code": str(CODE_PATH),
        "method": METHOD,
        "use_raw": "True",
        "seed": str(SEED),
    }
    paths["versions"].write_text(
        "\n".join(f"{key}={value}" for key, value in versions.items()) + "\n",
        encoding="utf-8",
    )
    paths["readme"].write_text(
        f"""BRCA {config['label']} final subtype DEG annotation

Input selected h5ad: {paths['input_h5ad']}
Raw cluster key: {raw_key}
Candidate table: {paths['agent_reviewed_candidate']}
Final mapping: {paths['mapping']}
Output annotated h5ad: {paths['output_h5ad']}
Post-annotation DEG directory: {paths['post_deg_dir']}

The mapping covers every raw cluster exactly once and uses unique
<lineage-prefix>_<most-significant-positive-DEG> subtype names. Full-length
post-annotation DEGs use t-test with use_raw=True. Off-lineage-like marker
patterns remain inside the current lineage object and are documented rather
than moved. The selected UMAP graph is preserved.
""",
        encoding="utf-8",
    )
    completion = {
        "lineage": LINEAGE_SLUG,
        "status": "completed",
        "n_cells": int(adata.n_obs),
        "n_raw_clusters": len(mapping),
        "n_cell_subtypes": int(adata.obs["cell_subtype"].nunique()),
        "mapping_complete_one_to_one": True,
        "cell_subtype_written": True,
        "functional_state_written": True,
        "post_annotation_full_length_degs_complete": True,
        "final_umap_complete": True,
        "marker_dotplot_complete": True,
        "annotated_h5ad": str(paths["output_h5ad"]),
        "agent_led_final_labeling_authorized": True,
        "elapsed_seconds": time.time() - started,
    }
    paths["completion"].write_text(json.dumps(completion, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(completion, indent=2), flush=True)


if __name__ == "__main__":
    main()
