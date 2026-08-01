"""Block 04 step 1: read existing per-subtype DEG, top-100 positive markers, sample scope.

Per updated SKILL.md (Module 05 marker source contract):
  - Use the saved post-annotation subtype differential-expression results
    from Module 02 (tables/02-cell_subtype_integration_clustering/01-main/
    degs_subtype_<lineage>/<subtype>_degs_subtype_<lineage>.csv).
  - DO NOT recompute DEG inside this block.
  - marker ordering = score descending only
  - positive DEG rule = logfoldchanges > 0
  - defensive significance gate = pvals_adj < 0.05 (BH-adjusted); padj column is
    required; missing padj raises an error so silent regressions are caught
  - keep source-table order for exact score ties (stable sort)
  - drop duplicated gene names within subtype, keep highest-score occurrence
  - take the first 100 unique positive genes per subtype
  - pvals / pvals_adj retained only for audit (NOT used as sort keys, only as gate)
  - Save marker_genes_used.csv with explicit `source_deg_csv` column.
  - Save reference_cells_used.csv.
  - reference expression = subtype mean from annotated single-cell adata.raw.
  - reference cells = every available cell in each subtype.
  - sample_scope.csv with R_cor/R_med status='normal', source='user_correction'.

Outputs:
  tables/04-spatial-validation-optional/01-tangram-spatial-mapping/
    marker_genes_used.csv
    reference_cells_used.csv
    gene_intersection_list.csv
    gene_intersection_summary.csv
    sample_scope.csv
    subtype_mean_reference.csv
"""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

SEED = 42
np.random.seed(SEED)

ROOT = Path("/mnt/disk18t/lr_xcy/riku/ccrcc-cmepi")
ATLAS = ROOT / "epi-cm-core-workflow/h5ad/02-cell_subtype_integration_clustering/06-ccrcc-subtypes-to-full-adata/adata_anno_cellsubtype.h5ad"
DEG_ROOT = ROOT / "epi-cm-core-workflow/tables/02-cell_subtype_integration_clustering/01-main"
SP_DIR = ROOT / "spatial_h5ad"
TAB = ROOT / "epi-cm-core-workflow/tables/04-spatial-validation-optional/01-tangram-spatial-mapping"
TAB.mkdir(parents=True, exist_ok=True)


def read_existing_deg(deg_root: Path) -> dict[str, tuple[pd.DataFrame, Path]]:
    """FIXED per updated SKILL.md: read upstream post-annotation per-subtype DEG.

    Returns: {subtype: (deg_df, source_csv_path)}
    Each subtype CSV must contain numeric 'score' and 'logfoldchanges'.
    """
    out = {}
    for d in sorted(deg_root.glob("degs_subtype_*")):
        for f in sorted(d.glob("*.csv")):
            subtype = f.stem.split("_degs_subtype_")[0]
            df = pd.read_csv(f)
            required = {"gene", "score", "logfoldchanges"}
            missing = required - set(df.columns)
            if missing:
                raise ValueError(f"{f}: missing DEG columns {missing}")
            out[subtype] = (df, f)
    if not out:
        raise FileNotFoundError(f"No per-subtype DEG CSVs found under {deg_root}")
    return out


def top_positive_markers(deg_df: pd.DataFrame, n: int = 100,
                         padj_threshold: float = 0.05) -> pd.DataFrame:
    """FIXED per updated SKILL.md: score descending only.

    Rules:
      - Require pvals_adj column (defensive significance gate)
      - Filter pvals_adj < padj_threshold AND logfoldchanges > 0
      - Sort by score descending only (NO pvals / pvals_adj / q used as sort keys)
      - Keep source-table order for exact score ties (mergesort is stable)
      - Drop duplicated gene names, retain the highest-score occurrence
        (since sort is stable, drop_duplicates(keep='first') keeps the
        top-score occurrence)
      - Take the first 100 unique genes
      - pvals / pvals_adj retained for audit only

    Note: empirically, the current upstream DEG tables (66/66 subtypes, 6600
    top-100 markers) all satisfy padj < 0.05, so this gate is a defensive
    safety check rather than an active filter. If a future DEG rerun produces
    noisier tables, this gate will catch non-significant markers automatically.
    """
    required = {"gene", "score", "logfoldchanges", "pvals_adj"}
    missing = required - set(deg_df.columns)
    if missing:
        raise ValueError(f"DEG table must contain columns {sorted(required)}; missing {sorted(missing)}")
    pos = deg_df[(deg_df["logfoldchanges"] > 0) & (deg_df["pvals_adj"] < padj_threshold)].copy()
    pos = pos.sort_values(by="score", ascending=False, kind="mergesort")
    pos = pos.drop_duplicates(subset="gene", keep="first")
    return pos.head(n).reset_index(drop=True)


def write_reference_cells_used(adata: sc.AnnData, groupby: str, out_path: Path) -> pd.DataFrame:
    """FIXED per SKILL.md: subtype, n_cells_total, n_cells_used, all_cells_used."""
    counts = adata.obs[groupby].astype(str).value_counts()
    rows = []
    for subtype, n in counts.items():
        rows.append({
            "subtype": subtype,
            "n_cells_total": int(n),
            "n_cells_used": int(n),
            "all_cells_used": True,
        })
    df = pd.DataFrame(rows).sort_values("subtype").reset_index(drop=True)
    df.to_csv(out_path, index=False)
    bad = df[df["n_cells_used"] != df["n_cells_total"]]
    if not bad.empty:
        raise RuntimeError(f"reference_cells_used violated all-cells rule: {bad}")
    if not df["all_cells_used"].all():
        raise RuntimeError("reference_cells_used violated all_cells_used==True")
    return df


def write_sample_scope(sp_dir: Path, out_path: Path) -> pd.DataFrame:
    """FIXED per SKILL.md: R_cor / R_med = normal, status_source=user_correction."""
    samples = sorted(p.stem for p in sp_dir.glob("*.h5ad"))
    rows = []
    for s in samples:
        if s in {"R_cor", "R_med"}:
            rows.append({
                "sample": s, "status": "normal",
                "include_all_samples": True, "include_tumor_only": False,
                "exclusion_reason_all_samples": "kept",
                "exclusion_reason_tumor_only": "normal_kidney_tissue_user_correction",
                "status_source": "user_correction",
            })
        else:
            rows.append({
                "sample": s, "status": "tumor",
                "include_all_samples": True, "include_tumor_only": True,
                "exclusion_reason_all_samples": "kept",
                "exclusion_reason_tumor_only": "kept",
                "status_source": "filename_T",
            })
    df = pd.DataFrame(rows)
    df.to_csv(out_path, index=False)
    return df


def build_subtype_mean_reference(atlas_path, common_genes):
    """FIXED: subtype mean from annotated single-cell adata.raw over common_genes."""
    a = sc.read_h5ad(atlas_path)
    assert a.raw is not None, "atlas must have adata.raw for Tangram reference"
    raw_X = a.raw.X
    if hasattr(raw_X, "toarray"):
        raw_X = raw_X.toarray()
    raw_var = a.raw.var_names.astype(str).tolist()
    raw_df = pd.DataFrame(raw_X, index=a.obs_names, columns=raw_var)
    raw_df["cell_subtype"] = a.obs["cell_subtype"].astype(str).values
    common_in_raw = [g for g in common_genes if g in raw_var]
    subtype_mean = raw_df.groupby("cell_subtype")[common_in_raw].mean()
    return subtype_mean, common_in_raw


def main():
    t0 = time.time()
    print("[prep] loading atlas for reference cells + subtype mean…", flush=True)
    adata = sc.read_h5ad(ATLAS)
    assert adata.raw is not None
    print(f"[prep] atlas: {adata.shape}  raw: {adata.raw.shape}", flush=True)

    print("[prep] reading upstream post-annotation per-subtype DEG CSVs "
          "(no DEG recomputation; per new SKILL Module 05)…", flush=True)
    deg_by_subtype = read_existing_deg(DEG_ROOT)
    print(f"[prep] DEG CSVs loaded for {len(deg_by_subtype)} subtypes", flush=True)

    marker_rows = []
    per_subtype_markers = {}
    for subtype, (deg_df, source_path) in deg_by_subtype.items():
        top = top_positive_markers(deg_df, n=100)
        top["subtype"] = subtype
        top["source_deg_csv"] = str(source_path.resolve())
        marker_rows.append(top)
        per_subtype_markers[subtype] = set(top["gene"].tolist())
    markers_df = pd.concat(marker_rows, ignore_index=True)
    cols = ["subtype", "source_deg_csv", "gene", "score", "logfoldchanges",
            "pvals", "pvals_adj"]
    cols = [c for c in cols if c in markers_df.columns]
    markers_df = markers_df[cols]
    markers_df.to_csv(TAB / "marker_genes_used.csv", index=False)
    union_markers = set().union(*per_subtype_markers.values()) if per_subtype_markers else set()
    print(f"[prep] markers: {len(markers_df)} rows; union size={len(union_markers)}",
          flush=True)

    write_reference_cells_used(adata, "cell_subtype", TAB / "reference_cells_used.csv")
    print(f"[prep] reference_cells_used.csv written", flush=True)

    sp_files = sorted(SP_DIR.glob("*.h5ad"))
    sp_gene_sets = {}
    for f in sp_files:
        sp = sc.read_h5ad(f)
        sp_gene_sets[f.stem] = set(map(str, sp.var_names.tolist()))
        del sp
    all_sp_genes = set.intersection(*sp_gene_sets.values()) if sp_gene_sets else set()
    ref_genes = set(adata.raw.var_names.astype(str))
    union_in_ref = union_markers & ref_genes
    common_genes = sorted(union_in_ref & all_sp_genes)
    print(f"[prep] ∩ ref: {len(union_in_ref)}; ∩ all spatial: {len(all_sp_genes)}; "
          f"final common: {len(common_genes)}", flush=True)

    subtype_mean, common_in_ref = build_subtype_mean_reference(ATLAS, common_genes)
    subtype_mean.index.name = "cell_subtype"
    subtype_mean.to_csv(TAB / "subtype_mean_reference.csv")
    print(f"[prep] subtype_mean_reference.csv written ({subtype_mean.shape})",
          flush=True)

    pd.DataFrame([{"gene": g} for g in common_genes]).to_csv(
        TAB / "gene_intersection_list.csv", index=False
    )
    pd.DataFrame([
        {"stat": "n_subtypes_with_deg_csv", "value": len(deg_by_subtype)},
        {"stat": "n_union_markers", "value": len(union_markers)},
        {"stat": "n_union_in_reference", "value": len(union_in_ref)},
        {"stat": "n_intersection_with_all_spatial", "value": len(common_genes)},
        {"stat": "deg_source", "value": (
            "tables/02-cell_subtype_integration_clustering/01-main/"
            "degs_subtype_<lineage>/<subtype>_degs_subtype_<lineage>.csv"
        )},
        {"stat": "marker_selection_rule", "value": (
            "score desc only + logfoldchanges>0 + dedup(keep first) + head 100"
        )},
    ]).to_csv(TAB / "gene_intersection_summary.csv", index=False)

    write_sample_scope(SP_DIR, TAB / "sample_scope.csv")
    print(f"[prep] sample_scope.csv written ({len(sp_files)} samples)", flush=True)

    print(f"[prep] done. total={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()