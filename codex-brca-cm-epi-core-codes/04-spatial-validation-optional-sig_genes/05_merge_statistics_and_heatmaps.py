#!/usr/bin/env python3
"""Merge six per-sample all-pair results and draw the two canonical heatmaps."""

from __future__ import annotations

import importlib.metadata
import json
import random
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import norm

SEED = 42
random.seed(SEED)
np.random.seed(SEED)
ROOT = Path("/mnt/disk18t/lr_xcy/riku/brca_cm-epi_codex")
WORKFLOW = ROOT / "epi-cm-core-workflow"
TABLE_ROOT = WORKFLOW / "tables/04-spatial-validation-optional-sig_genes/03-spatial-statistics-and-plotting"
PER_SAMPLE_DIR = TABLE_ROOT / "per_sample"
SCOPE_SOURCE = WORKFLOW / "tables/04-spatial-validation-optional-sig_genes/01-input-audit-and-reference-v4/sample_scope.csv"
FIGURE_ROOT = WORKFLOW / "figures/04-spatial-validation-optional-sig_genes/03-spatial-statistics-and-plotting/stouffer_heatmaps"
SAMPLES = ["1142243F", "1160920F", "CID4290", "CID4465", "CID44971", "CID4535"]

plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "text.usetex": False, "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9})


def pkg(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def bh_fdr(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.notna() & np.isfinite(pd.to_numeric(values, errors="coerce"))
    if not valid.any():
        return result
    p = values.loc[valid].astype(float).to_numpy()
    order = np.argsort(p, kind="stable")
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1].clip(0, 1)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    result.loc[valid] = restored
    return result


def cm_key(value: str) -> tuple[int, str]:
    match = re.search(r"CM(\d+)$", str(value))
    return (int(match.group(1)) if match else 10**9, str(value))


def q_label(q: float) -> str:
    if not np.isfinite(q): return "ns"
    if q < 0.001: return "***"
    if q < 0.01: return "**"
    if q < 0.05: return "*"
    return "ns"


def combine_stouffer(frame: pd.DataFrame, scope: str) -> pd.DataFrame:
    rows = []
    for (cm, epi), group in frame.groupby(["CM", "epi_subtype"], sort=False):
        values = pd.to_numeric(group["fisher_signed_z"], errors="coerce")
        valid = values.notna() & np.isfinite(values)
        n = int(valid.sum())
        z = float(values.loc[valid].sum() / np.sqrt(n)) if n else np.nan
        p = float(2 * norm.sf(abs(z))) if np.isfinite(z) else np.nan
        rows.append({"scope": scope, "CM": cm, "epi_subtype": epi, "combined_signed_z": z, "p_value": p, "n_samples": n, "samples": ";".join(group.loc[valid, "sample"].astype(str)), "direction": "positive" if z > 0 else ("negative" if z < 0 else "neutral")})
    result = pd.DataFrame(rows)
    result["q_value_bh"] = bh_fdr(result["p_value"])
    result["significance"] = result["q_value_bh"].map(q_label)
    return result.sort_values(["CM", "epi_subtype"], key=lambda s: s.map(cm_key) if s.name == "CM" else s, kind="stable")


def plot_heatmap(frame: pd.DataFrame, output_dir: Path, stem: str, title: str) -> tuple[Path, Path]:
    z = frame.pivot(index="epi_subtype", columns="CM", values="combined_signed_z")
    q = frame.pivot(index="epi_subtype", columns="CM", values="q_value_bh")
    rows = sorted(z.index)
    cols = sorted(z.columns, key=cm_key)
    z, q = z.reindex(index=rows, columns=cols), q.reindex(index=rows, columns=cols)
    annot = q.applymap(q_label)
    finite = z.to_numpy()[np.isfinite(z.to_numpy())]
    vmax = float(np.max(np.abs(finite))) if finite.size else 1.0
    if vmax == 0: vmax = 1.0
    fig_w, fig_h = max(7.8, 0.52 * z.shape[1] + 2.6), max(4.8, 0.42 * z.shape[0] + 1.9)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(z, cmap="coolwarm", center=0, vmin=-vmax, vmax=vmax, linewidths=0.25, linecolor="white", annot=annot, fmt="", annot_kws={"fontsize": 6.5, "color": "black", "linespacing": 0.9}, cbar_kws={"label": "combined signed Z"}, square=True, ax=ax)
    ax.set_title(title); ax.set_xlabel("CM"); ax.set_ylabel("Epithelial subtype"); ax.tick_params(axis="x", labelrotation=90); ax.tick_params(axis="y", labelrotation=0)
    fig.tight_layout()
    output_dir.mkdir(parents=True, exist_ok=False)
    pdf, svg = output_dir / f"{stem}.pdf", output_dir / f"{stem}.svg"
    fig.savefig(pdf, bbox_inches="tight", dpi=300); fig.savefig(svg, bbox_inches="tight", dpi=300); plt.close(fig)
    return pdf, svg


def main() -> None:
    completion_path = TABLE_ROOT / "statistics_and_heatmaps_completion.json"
    if completion_path.exists():
        raise FileExistsError(completion_path)
    frames, manifests = [], []
    for sample in SAMPLES:
        completion = json.loads((PER_SAMPLE_DIR / f"{sample}_plot_stats_completion.json").read_text())
        if completion.get("status") != "completed" or completion.get("n_pairs") != 100:
            raise RuntimeError(f"Incomplete sample plotting/statistics: {sample}")
        frames.append(pd.read_csv(PER_SAMPLE_DIR / f"{sample}_pair_statistics.csv"))
        manifests.append(pd.read_csv(PER_SAMPLE_DIR / f"{sample}_pair_manifest.csv"))
    all_stats = pd.concat(frames, ignore_index=True)
    manifest = pd.concat(manifests, ignore_index=True)
    if len(all_stats) != 600 or len(manifest) != 600:
        raise ValueError("Expected 600 sample x CM x Epi rows")
    all_stats["spearman_q_value_bh_600"] = bh_fdr(all_stats["spearman_p_value"])
    all_stats["fisher_q_value_bh_600"] = bh_fdr(all_stats["fisher_p_value"])
    all_stats.to_csv(TABLE_ROOT / "per_sample_all_pair_statistics.csv", index=False)
    all_stats[["sample", "CM", "epi_subtype", "n_spots", "spearman_rho", "spearman_p_value", "spearman_q_value_bh_600"]].to_csv(TABLE_ROOT / "per_sample_spearman.csv", index=False)
    all_stats[["sample", "CM", "epi_subtype", "percentile_cutoff", "both_low", "epi_high_only", "cm_high_only", "both_high", "fisher_odds_ratio", "fisher_p_value", "fisher_q_value_bh_600", "fisher_direction", "fisher_signed_z"]].to_csv(TABLE_ROOT / "percentile_quadrant_fisher_per_sample.csv", index=False)
    manifest.to_csv(TABLE_ROOT / "all_sample_cm_epi_pair_manifest.csv", index=False)

    scope = pd.read_csv(SCOPE_SOURCE)
    if set(scope["sample"].astype(str)) != set(SAMPLES):
        raise ValueError("Sample scope does not cover the six spatial samples")
    scope.to_csv(TABLE_ROOT / "sample_scope.csv", index=False)
    all_members = scope.loc[scope["include_all_samples"].astype(bool), "sample"].astype(str).tolist()
    tumor_members = scope.loc[scope["include_tumor_only"].astype(bool), "sample"].astype(str).tolist()
    all_stouffer = combine_stouffer(all_stats.loc[all_stats["sample"].isin(all_members)], "all-samples")
    tumor_stouffer = combine_stouffer(all_stats.loc[all_stats["sample"].isin(tumor_members)], "tumor-only")
    all_dir = TABLE_ROOT / "statistics/all-samples"; tumor_dir = TABLE_ROOT / "statistics/tumor-only"
    all_dir.mkdir(parents=True, exist_ok=False); tumor_dir.mkdir(parents=True, exist_ok=False)
    all_csv = all_dir / "percentile_quadrant_fisher_sample_stouffer_all_samples.csv"
    tumor_csv = tumor_dir / "percentile_quadrant_fisher_sample_stouffer_tumor_only.csv"
    all_stouffer.to_csv(all_csv, index=False); tumor_stouffer.to_csv(tumor_csv, index=False)
    all_pdf, all_svg = plot_heatmap(all_stouffer, FIGURE_ROOT / "all-samples", "all_samples_stouffer_signedZ_qstars", "Spatial CM-Epi Fisher Stouffer (all samples)")
    tumor_pdf, tumor_svg = plot_heatmap(tumor_stouffer, FIGURE_ROOT / "tumor-only", "tumor_only_stouffer_signedZ_qstars", "Spatial CM-Epi Fisher Stouffer (tumor only)")

    pd.DataFrame([{"code_file": str(Path(__file__).resolve()), "sample_scope_source": str(SCOPE_SOURCE.resolve()), "all_samples": ";".join(all_members), "tumor_only_samples": ";".join(tumor_members), "same_membership_because_all_six_primary_tumor_specimens": set(all_members) == set(tumor_members), "per_sample_pairs": 600, "stouffer_pairs_per_scope": 100, "fisher_signed_z": "sign(cross-product difference) * two-sided normal quantile", "stouffer": "sum signed sample Z / sqrt(n)", "stouffer_bh_family": "100 CM x Epi pairs separately per scope", "heatmap_annotation": "BH q-value stars/ns only; no numeric cell labels", "seed": SEED}]).to_csv(TABLE_ROOT / "statistics_and_heatmap_parameters.csv", index=False)
    (TABLE_ROOT / "readme.txt").write_text("Inputs are the six mandatory Tangram spot-score CSVs reloaded by the per-sample plotting/statistics scripts. Every 10 CM x 10 Epi pair is tested in every sample. Raw CMact/EPIfrac values are used for Spearman; within-sample percentile ranks at 0.5 define Fisher quadrants. Signed Fisher Z values are combined by Stouffer across explicit all-samples and tumor-only scopes. Both scopes contain the same six primary-tumor specimens and are still saved separately.\n", encoding="utf-8")
    (TABLE_ROOT / "package_versions.txt").write_text(f"python={sys.version.split()[0]}\nnumpy={pkg('numpy')}\npandas={pkg('pandas')}\nscipy={pkg('scipy')}\nmatplotlib={pkg('matplotlib')}\nseaborn={pkg('seaborn')}\ncode={Path(__file__).resolve()}\nseed={SEED}\n", encoding="utf-8")
    completion = {"status": "completed", "n_samples": len(SAMPLES), "n_sample_pair_rows": len(all_stats), "n_pairs_per_scope": len(all_stouffer), "all_samples_members": all_members, "tumor_only_members": tumor_members, "scope_memberships_identical": set(all_members) == set(tumor_members), "all_samples_table": str(all_csv.resolve()), "tumor_only_table": str(tumor_csv.resolve()), "heatmaps": [str(all_pdf.resolve()), str(all_svg.resolve()), str(tumor_pdf.resolve()), str(tumor_svg.resolve())], "seed": SEED}
    completion_path.write_text(json.dumps(completion, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(completion, indent=2), flush=True)


if __name__ == "__main__":
    main()
