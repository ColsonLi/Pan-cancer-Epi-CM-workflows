#!/usr/bin/env python3
"""Project all final BRCA lineage subtype annotations to the strict full atlas."""

from __future__ import annotations

import colorsys
import importlib.metadata
import json
import subprocess
import sys
import time
from pathlib import Path

import anndata as ad
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import scanpy as sc


ROOT = Path("/mnt/disk18t/lr_xcy/riku/brca_cm-epi_codex")
WORKFLOW = ROOT / "epi-cm-core-workflow"
BLOCK = "02-cell_subtype_integration_clustering"
SOURCE_H5AD = (
    WORKFLOW
    / "h5ad"
    / "01-celltype_integration_clustering"
    / "07-score-rank-qc"
    / "adata_anno_score_genes_rank_consistent.h5ad"
)
OUTPUT_H5AD = (
    WORKFLOW
    / "h5ad"
    / BLOCK
    / "06-project-subtypes-to-full-adata"
    / "adata_anno_cellsubtype.h5ad"
)
TABLE_ROOT = WORKFLOW / "tables" / BLOCK / "06-project-subtypes-to-full-adata"
FIGURE_ROOT = WORKFLOW / "figures" / BLOCK / "06-project-subtypes-to-full-adata"
CODE_PATH = Path(__file__).resolve()
SEED = 42

LINEAGES = {
    "epithelial": {"abbrev": "epi", "label": "Epithelial Cells"},
    "t_cells": {"abbrev": "t", "label": "T Cells"},
    "myeloid": {"abbrev": "mye", "label": "Myeloid Cells"},
    "b_cells": {"abbrev": "b", "label": "B Cells"},
    "plasma": {"abbrev": "plasma", "label": "Plasma Cells"},
    "endothelial": {"abbrev": "endo", "label": "Endothelial Cells"},
    "stromal": {"abbrev": "stromal", "label": "Stromal Cells"},
    "perivascular": {"abbrev": "pvl", "label": "Perivascular Cells"},
}
BROAD_ORDER = [config["label"] for config in LINEAGES.values()]


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def pdf_is_readable(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    result = subprocess.run(["pdfinfo", str(path)], capture_output=True, text=True)
    return result.returncode == 0 and "Pages:" in result.stdout


def subtype_shades(base_color: str, n: int) -> list[str]:
    base = mcolors.to_rgb(base_color)
    hue, saturation, value = colorsys.rgb_to_hsv(*base)
    colors: list[str] = []
    for index in range(n):
        fraction = 0 if n == 1 else index / (n - 1)
        shifted_hue = (hue + (fraction - 0.5) * 0.10) % 1.0
        shifted_saturation = min(0.95, max(0.45, saturation * (0.75 + 0.35 * fraction)))
        shifted_value = min(0.95, max(0.45, 0.92 - 0.38 * fraction))
        colors.append(mcolors.to_hex(colorsys.hsv_to_rgb(shifted_hue, shifted_saturation, shifted_value)))
    return colors


def plot_umap_formats(adata: ad.AnnData, color: str | list[str], suffix: str, **kwargs: object) -> tuple[Path, Path]:
    sc.pl.umap(adata, color=color, show=False, save=f"_{suffix}.pdf", **kwargs)
    sc.pl.umap(adata, color=color, show=False, save=f"_{suffix}.svg", **kwargs)
    pdf = Path(sc.settings.figdir) / f"umap_{suffix}.pdf"
    svg = Path(sc.settings.figdir) / f"umap_{suffix}.svg"
    if not pdf_is_readable(pdf) or not svg.exists() or svg.stat().st_size == 0:
        raise FileNotFoundError(f"UMAP PDF/SVG output validation failed: {pdf}, {svg}")
    return pdf, svg


def main() -> None:
    started = time.time()
    required_final = [
        OUTPUT_H5AD,
        TABLE_ROOT / "projection_summary.csv",
        TABLE_ROOT / "projection_completion.json",
        TABLE_ROOT / "cell_subtype_color_mapping.csv",
        FIGURE_ROOT / "umap_leiden_coarse_vs_projected_cell_type.pdf",
        FIGURE_ROOT / "umap_projected_cell_subtype.pdf",
    ]
    if any(path.exists() for path in required_final):
        if all(path.exists() for path in required_final):
            completion = json.loads((TABLE_ROOT / "projection_completion.json").read_text())
            if completion.get("status") == "completed" and completion.get("unmatched_full_cells") == 0:
                print(json.dumps({"status": "valid_existing_projection_reused", "output_h5ad": str(OUTPUT_H5AD)}, indent=2))
                return
        raise FileExistsError("Partial subtype projection outputs already exist; refusing overwrite.")

    full = sc.read_h5ad(SOURCE_H5AD)
    if full.n_obs != 80406 or not full.obs_names.is_unique:
        raise ValueError(f"Unexpected strict full atlas cells/IDs: {full.n_obs}")
    if "X_umap" not in full.obsm or full.obsm["X_umap"].shape != (full.n_obs, 2):
        raise ValueError("Strict full atlas lacks the validated full-atlas X_umap.")
    required_obs = {"leiden_coarse", "cell_type", "sample", "series", "status"}
    if not required_obs.issubset(full.obs.columns):
        raise ValueError(f"Full atlas lacks obs columns: {sorted(required_obs - set(full.obs))}")

    full.obs["leiden_coarse_before_subtype_projection"] = full.obs["leiden_coarse"].copy()
    full.obs["cell_type_before_subtype_projection"] = full.obs["cell_type"].copy()
    projection_parts: list[pd.DataFrame] = []
    reports: list[dict[str, object]] = []
    subtype_order: list[str] = []
    lineage_subtype_order: dict[str, list[str]] = {}
    all_projected_ids: list[str] = []

    for lineage, config in LINEAGES.items():
        annotated_h5ad = (
            WORKFLOW / "h5ad" / BLOCK / "05-subtype-deg-annotation"
            / lineage / f"adata_{config['abbrev']}.h5ad"
        )
        parameters = pd.read_csv(
            WORKFLOW / "tables" / BLOCK / "05-subtype-deg-annotation"
            / lineage / "final_subtype_annotation_parameters.csv"
        ).iloc[0]
        raw_key = str(parameters["raw_cluster_key"])
        lineage_adata = ad.read_h5ad(annotated_h5ad, backed="r")
        try:
            required_lineage_obs = {
                "cell_type", "cell_subtype", "functional_state",
                "annotation_confidence", "marker_evidence_lineage",
                "prefix_rationale", "annotation_note", "annotation_source", raw_key,
            }
            if not required_lineage_obs.issubset(lineage_adata.obs.columns):
                raise ValueError(
                    f"{lineage} annotated h5ad lacks columns: "
                    f"{sorted(required_lineage_obs - set(lineage_adata.obs))}"
                )
            ids = lineage_adata.obs_names.astype(str)
            missing = ids.difference(full.obs_names.astype(str))
            duplicated = ids[ids.duplicated()]
            if len(missing) or len(duplicated):
                raise ValueError(
                    f"{lineage} projection IDs invalid: missing={len(missing)}, duplicated={len(duplicated)}"
                )
            subtype_categories = lineage_adata.obs["cell_subtype"].cat.categories.astype(str).tolist()
            lineage_subtype_order[lineage] = subtype_categories
            subtype_order.extend(subtype_categories)
            part = pd.DataFrame(index=ids)
            part["cell_type"] = lineage_adata.obs["cell_type"].astype(str).to_numpy()
            part["cell_subtype"] = lineage_adata.obs["cell_subtype"].astype(str).to_numpy()
            part["functional_state"] = lineage_adata.obs["functional_state"].astype(str).to_numpy()
            part["annotation_confidence"] = lineage_adata.obs["annotation_confidence"].astype(str).to_numpy()
            part["marker_evidence_lineage"] = lineage_adata.obs["marker_evidence_lineage"].astype(str).to_numpy()
            part["prefix_rationale"] = lineage_adata.obs["prefix_rationale"].astype(str).to_numpy()
            part["annotation_note"] = lineage_adata.obs["annotation_note"].astype(str).to_numpy()
            part["annotation_source"] = lineage_adata.obs["annotation_source"].astype(str).to_numpy()
            part["lineage_raw_cluster"] = lineage_adata.obs[raw_key].astype(str).to_numpy()
            part["lineage_raw_cluster_key"] = raw_key
            part["subtype_projection_lineage"] = lineage
            part["subtype_projection_source_h5ad"] = str(annotated_h5ad)
            projection_parts.append(part)
            all_projected_ids.extend(ids.tolist())
            report = {
                "lineage": lineage,
                "target_leiden_coarse": config["label"],
                "source_h5ad": str(annotated_h5ad),
                "target_h5ad": str(SOURCE_H5AD),
                "n_source_cells": int(lineage_adata.n_obs),
                "n_exact_obs_name_matches": int(lineage_adata.n_obs),
                "n_unmatched_source_cells": 0,
                "n_duplicated_source_ids": 0,
                "n_cell_subtypes": len(subtype_categories),
                "raw_cluster_key": raw_key,
                "projection_rule": "exact obs_names match",
                "old_columns_backed_up": "leiden_coarse;cell_type",
                "projected_columns": "cell_type;cell_subtype;functional_state;annotation_confidence;marker_evidence_lineage;prefix_rationale;annotation_note;annotation_source;lineage_raw_cluster;lineage_raw_cluster_key",
            }
            reports.append(report)
        finally:
            lineage_adata.file.close()

    if len(all_projected_ids) != len(set(all_projected_ids)):
        raise ValueError("Annotated lineage h5ad objects overlap in cell IDs.")
    if set(all_projected_ids) != set(full.obs_names.astype(str)):
        raise ValueError(
            "Annotated lineage union does not exactly cover strict full atlas: "
            f"union={len(set(all_projected_ids))}, full={full.n_obs}"
        )
    projection = pd.concat(projection_parts, axis=0)
    projection = projection.loc[full.obs_names.astype(str)]
    if projection.isna().any().any():
        raise ValueError("Projection table contains missing values.")
    if len(subtype_order) != len(set(subtype_order)):
        raise ValueError("cell_subtype labels are not globally unique across lineages.")

    full.obs["cell_type"] = pd.Categorical(
        projection["cell_type"], categories=BROAD_ORDER, ordered=True
    )
    full.obs["cell_subtype"] = pd.Categorical(
        projection["cell_subtype"], categories=subtype_order, ordered=True
    )
    full.obs["functional_state"] = pd.Categorical(projection["functional_state"])
    full.obs["annotation_confidence"] = pd.Categorical(projection["annotation_confidence"])
    full.obs["marker_evidence_lineage"] = pd.Categorical(projection["marker_evidence_lineage"])
    for column in [
        "prefix_rationale", "annotation_note", "annotation_source",
        "lineage_raw_cluster", "lineage_raw_cluster_key",
        "subtype_projection_lineage", "subtype_projection_source_h5ad",
    ]:
        full.obs[column] = projection[column].astype(str).to_numpy()
    if not full.obs["cell_type"].astype(str).equals(
        full.obs["cell_type_before_subtype_projection"].astype(str)
    ):
        raise ValueError("Projected broad cell_type changed from the trusted strict atlas.")

    broad_colors = list(full.uns.get("leiden_coarse_colors", []))
    if len(broad_colors) != len(BROAD_ORDER):
        broad_colors = [mcolors.to_hex(color) for color in sc.pl.palettes.default_20[:len(BROAD_ORDER)]]
    full.obs["leiden_coarse"] = pd.Categorical(
        full.obs["leiden_coarse"].astype(str), categories=BROAD_ORDER, ordered=True
    )
    full.uns["leiden_coarse_colors"] = broad_colors
    full.uns["cell_type_colors"] = broad_colors

    subtype_colors: list[str] = []
    color_rows: list[dict[str, object]] = []
    for (lineage, config), base_color in zip(LINEAGES.items(), broad_colors):
        subtypes = lineage_subtype_order[lineage]
        colors = subtype_shades(base_color, len(subtypes))
        subtype_colors.extend(colors)
        for subtype, color in zip(subtypes, colors):
            color_rows.append(
                {
                    "lineage": lineage,
                    "cell_type": config["label"],
                    "cell_subtype": subtype,
                    "hex_color": color,
                    "broad_base_color": base_color,
                }
            )
    full.uns["cell_subtype_colors"] = subtype_colors
    TABLE_ROOT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(color_rows).to_csv(TABLE_ROOT / "cell_subtype_color_mapping.csv", index=False)
    pd.DataFrame(
        {"cell_type": BROAD_ORDER, "hex_color": broad_colors}
    ).to_csv(TABLE_ROOT / "broad_cell_type_color_mapping.csv", index=False)

    OUTPUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
    full.write_h5ad(OUTPUT_H5AD, compression="gzip")
    saved = ad.read_h5ad(OUTPUT_H5AD, backed="r")
    try:
        if (
            saved.n_obs != 80406
            or "cell_subtype" not in saved.obs
            or int(saved.obs["cell_subtype"].nunique()) != len(subtype_order)
            or saved.obs["cell_subtype"].isna().any()
            or saved.raw is None
        ):
            raise ValueError("Saved projected full atlas failed validation.")
    finally:
        saved.file.close()

    FIGURE_ROOT.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = FIGURE_ROOT
    sc.settings.file_format_figs = "pdf"
    sc.set_figure_params(figsize=(5, 5), dpi=150, fontsize=8)
    broad_pdf, broad_svg = plot_umap_formats(
        full,
        ["leiden_coarse", "cell_type"],
        "leiden_coarse_vs_projected_cell_type",
        ncols=2,
        wspace=0.4,
    )
    sc.set_figure_params(figsize=(10, 10), dpi=150, fontsize=7)
    subtype_pdf, subtype_svg = plot_umap_formats(
        full,
        "cell_subtype",
        "projected_cell_subtype",
        legend_fontsize=6,
    )

    report_by_lineage = {row["lineage"]: row for row in reports}
    color_lookup = dict(
        zip(
            pd.DataFrame(color_rows)["cell_subtype"],
            pd.DataFrame(color_rows)["hex_color"],
        )
    )
    for lineage, config in LINEAGES.items():
        annotated_h5ad = Path(str(report_by_lineage[lineage]["source_h5ad"]))
        lineage_adata = sc.read_h5ad(annotated_h5ad)
        subtypes = lineage_subtype_order[lineage]
        lineage_adata.obs["cell_subtype"] = pd.Categorical(
            lineage_adata.obs["cell_subtype"].astype(str),
            categories=subtypes,
            ordered=True,
        )
        lineage_adata.uns["cell_subtype_colors"] = [color_lookup[value] for value in subtypes]
        lineage_figdir = FIGURE_ROOT / lineage
        lineage_figdir.mkdir(parents=True, exist_ok=True)
        sc.settings.figdir = lineage_figdir
        sc.set_figure_params(figsize=(5, 5), dpi=150, fontsize=8)
        line_pdf, line_svg = plot_umap_formats(
            lineage_adata,
            "cell_subtype",
            f"{lineage}_cell_subtype_with_projected_subtype_palette",
            legend_fontsize=7,
        )
        report_by_lineage[lineage]["broad_palette_csv"] = str(
            TABLE_ROOT / "broad_cell_type_color_mapping.csv"
        )
        report_by_lineage[lineage]["cell_subtype_palette_csv"] = str(
            TABLE_ROOT / "cell_subtype_color_mapping.csv"
        )
        report_by_lineage[lineage]["umap_basis"] = "X_umap"
        report_by_lineage[lineage]["lineage_subtype_umap_pdf"] = str(line_pdf)
        report_by_lineage[lineage]["lineage_subtype_umap_svg"] = str(line_svg)
        lineage_table_dir = TABLE_ROOT / lineage
        lineage_table_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame([report_by_lineage[lineage]]).to_csv(
            lineage_table_dir / "projection_match_report.csv", index=False
        )

    summary = pd.DataFrame(list(report_by_lineage.values()))
    summary.to_csv(TABLE_ROOT / "projection_summary.csv", index=False)
    pd.DataFrame(
        [
            {
                "source_h5ad": str(SOURCE_H5AD),
                "output_h5ad": str(OUTPUT_H5AD),
                "n_full_cells": int(full.n_obs),
                "n_projected_cells": int(len(projection)),
                "n_cell_types": int(full.obs["cell_type"].nunique()),
                "n_cell_subtypes": int(full.obs["cell_subtype"].nunique()),
                "projection_rule": "exact obs_names match",
                "old_columns_backed_up": "leiden_coarse;cell_type",
                "broad_comparison_pdf": str(broad_pdf),
                "broad_comparison_svg": str(broad_svg),
                "subtype_umap_pdf": str(subtype_pdf),
                "subtype_umap_svg": str(subtype_svg),
                "code_file": str(CODE_PATH),
                "seed": SEED,
            }
        ]
    ).to_csv(TABLE_ROOT / "projection_parameters.csv", index=False)
    versions = {
        "python": sys.version.split()[0],
        "scanpy": package_version("scanpy"),
        "anndata": package_version("anndata"),
        "matplotlib": package_version("matplotlib"),
        "code": str(CODE_PATH),
        "seed": str(SEED),
    }
    (TABLE_ROOT / "package_versions.txt").write_text(
        "\n".join(f"{key}={value}" for key, value in versions.items()) + "\n",
        encoding="utf-8",
    )
    (TABLE_ROOT / "readme.txt").write_text(
        f"""BRCA subtype projection to strict full atlas

Source full atlas: {SOURCE_H5AD}
Output projected atlas: {OUTPUT_H5AD}
Code: {CODE_PATH}

All eight final annotated lineage h5ad objects were projected by exact obs_names
match. Their non-overlapping union exactly covers all {full.n_obs} cells. Original
leiden_coarse and cell_type columns are retained in backup columns before the
validated broad cell_type and 68 unique cell_subtype labels are written. Palette
mappings and PDF/SVG projection QC figures are saved with per-lineage match
reports.
""",
        encoding="utf-8",
    )
    completion = {
        "status": "completed",
        "source_full_h5ad": str(SOURCE_H5AD),
        "output_projected_h5ad": str(OUTPUT_H5AD),
        "n_full_cells": int(full.n_obs),
        "n_projected_cells": int(len(projection)),
        "unmatched_full_cells": 0,
        "duplicated_projection_ids": 0,
        "n_lineages": len(LINEAGES),
        "n_cell_subtypes": int(full.obs["cell_subtype"].nunique()),
        "old_columns_backed_up": True,
        "projection_qc_figures_complete": True,
        "elapsed_seconds": time.time() - started,
    }
    (TABLE_ROOT / "projection_completion.json").write_text(
        json.dumps(completion, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(completion, indent=2), flush=True)


if __name__ == "__main__":
    main()
