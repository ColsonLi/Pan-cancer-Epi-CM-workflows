#!/usr/bin/env python3
"""Final integrity audit for the latest-skill canonical spatial run."""

from __future__ import annotations

import importlib.metadata
import json
from pathlib import Path
import platform

import numpy as np
import pandas as pd


ROOT = Path("/mnt/disk18t/lr_xcy/riku/crc_codex/epi-cm-core-workflow")
TASK = "05-canonical-epifrac-composition-v2"
BLOCK = "04-spatial-validation-optional-sig_genes"
TABLE = ROOT / f"tables/{BLOCK}/{TASK}"
FIGURE = ROOT / f"figures/{BLOCK}/{TASK}"
OUT = TABLE / "04-completion-audit"
EXPECTED_SAMPLES = 37
EXPECTED_PAIRS = 6105
EXPECTED_FIGURES = 24420


def version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "not-installed"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    checks = []

    def check(name, observed, expected):
        checks.append({"check": name, "observed": observed, "expected": expected, "passed": observed == expected})

    ref = TABLE / "01-reference-and-manifest"
    mapping = TABLE / "02-tangram-mapping"
    stats = TABLE / "03-all-pair-statistics-and-plots"
    fig = FIGURE / "03-all-pair-statistics-and-plots"

    scope = pd.read_csv(ref / "spatial_sample_scope_manifest.csv")
    check("spatial_samples", len(scope), EXPECTED_SAMPLES)
    check("tumor_samples", int(scope.status.eq("tumor").sum()), 26)
    check("normal_samples", int(scope.status.eq("normal").sum()), 11)

    marker = pd.read_csv(ref / "marker_genes_used.csv")
    check("marker_subtypes", marker.cell_subtype.nunique(), 87)
    check("marker_rows", len(marker), 8617)
    marker_significant_positive = bool(
        pd.to_numeric(marker["pvals_adj"], errors="coerce").lt(0.05).all()
        and pd.to_numeric(marker["logfoldchanges"], errors="coerce").gt(0).all()
    )
    check("marker_rows_significant_positive", int(marker_significant_positive), 1)
    used_genes = pd.read_csv(ref / "tangram_common_deg_genes.csv").gene.astype(str)
    check("tangram_genes", used_genes.nunique(), 247)
    reference_cells = pd.read_csv(ref / "reference_cells_used.csv")
    all_cells_ok = bool(
        reference_cells["all_cells_used"].astype(bool).all()
        and reference_cells["n_cells_total"].eq(reference_cells["n_cells_used"]).all()
    )
    check("reference_subtypes_all_cells_used", int(all_cells_ok), 1)

    score_paths = sorted((mapping / "spot_scores").glob("*_tangram_pseudobulk_epi_cm_spot_scores.csv"))
    check("spot_score_tables", len(score_paths), EXPECTED_SAMPLES)
    run_manifest = pd.read_csv(mapping / "tangram_sample_run_manifest.csv")
    check("tangram_completed_samples", int(run_manifest.status_run.eq("completed").sum()), EXPECTED_SAMPLES)
    pair_manifest = pd.read_csv(mapping / "all_sample_cm_epi_pair_manifest.csv")
    check("sample_cm_epi_pairs", len(pair_manifest), EXPECTED_PAIRS)

    epi_rows = epi_zero = epi_bad = 0
    epi_max_error = 0.0
    for path in score_paths:
        columns = pd.read_csv(path, nrows=0).columns
        epi_cols = [column for column in columns if column.startswith("EPIfrac__")]
        check(f"{path.stem}_epi_columns", len(epi_cols), 11)
        for chunk in pd.read_csv(path, usecols=epi_cols, chunksize=100000):
            sums = chunk.sum(axis=1).to_numpy(float)
            epi_rows += len(sums)
            zero = np.isclose(sums, 0.0, atol=1e-8)
            epi_zero += int(zero.sum())
            errors = np.abs(sums[~zero] - 1.0)
            if errors.size:
                epi_bad += int((errors > 1e-6).sum())
                epi_max_error = max(epi_max_error, float(errors.max()))
    check("epifrac_rows", epi_rows, 3706544)
    check("epifrac_bad_positive_rows", epi_bad, 0)

    worker = pd.read_csv(stats / "plot_worker_status.csv")
    check("plot_workers_completed", int(worker.status.isin(["completed", "already_completed"]).sum()), EXPECTED_SAMPLES)
    check("plot_workers_failed", int(worker.status.eq("failed").sum()), 0)

    manifest_paths = sorted((stats / "plot_manifests_by_sample").glob("*_plot_manifest.csv"))
    check("sample_plot_manifests", len(manifest_paths), EXPECTED_SAMPLES)
    merged = pd.concat([pd.read_csv(path) for path in manifest_paths], ignore_index=True)
    check("merged_plot_manifest_rows", len(merged), EXPECTED_PAIRS)
    path_columns = ["raw_pdf", "raw_svg", "percentile_pdf", "percentile_svg"]
    figure_paths = [Path(path) for column in path_columns for path in merged[column].astype(str)]
    check("pair_figure_paths", len(figure_paths), EXPECTED_FIGURES)
    check("pair_figures_present", sum(path.is_file() for path in figure_paths), EXPECTED_FIGURES)
    check("pair_figures_nonempty", sum(path.is_file() and path.stat().st_size > 0 for path in figure_paths), EXPECTED_FIGURES)
    check("png_files", sum(1 for _ in fig.rglob("*.png")), 0)

    spearman = pd.read_csv(stats / "per_sample_spearman.csv")
    fisher = pd.read_csv(stats / "percentile_quadrant_fisher_per_sample.csv")
    stouffer_all = pd.read_csv(stats / "percentile_quadrant_fisher_sample_stouffer_all_samples.csv")
    stouffer_tumor = pd.read_csv(stats / "percentile_quadrant_fisher_sample_stouffer_tumor_only.csv")
    check("spearman_rows", len(spearman), EXPECTED_PAIRS)
    check("fisher_rows", len(fisher), EXPECTED_PAIRS)
    check("stouffer_all_rows", len(stouffer_all), 165)
    check("stouffer_tumor_rows", len(stouffer_tumor), 165)

    merged.to_csv(OUT / "all_sample_plot_manifest.csv", index=False)
    checks_df = pd.DataFrame(checks)
    checks_df.to_csv(OUT / "spatial_validation_checks.csv", index=False)
    passed = bool(checks_df.passed.all())
    summary = {
        "latest_skill_sha256": "d2416f0a721104872d00095ffb21c35d4ca44b33dff09758d8fc1bb493d9f0d0",
        "all_checks_passed": passed,
        "n_spatial_samples": len(scope),
        "n_tumor_samples": int(scope.status.eq("tumor").sum()),
        "n_tangram_genes": used_genes.nunique(),
        "n_pairs": len(merged),
        "n_pair_figure_files": len(figure_paths),
        "n_png": 0,
        "epifrac_definition": "within-observation epithelial-internal composition; positive rows sum to 1",
        "epifrac_rows": epi_rows,
        "epifrac_zero_rows": epi_zero,
        "epifrac_bad_positive_rows": epi_bad,
        "epifrac_max_abs_sum_error": epi_max_error,
        "marker_ordering": "pvals_adj < 0.05; logfoldchanges > 0; score descending stable; duplicate genes dropped; top 100",
    }
    (OUT / "spatial_validation_completion_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    versions = {name: version(name) for name in ["anndata", "h5py", "matplotlib", "numpy", "pandas", "scipy", "tangram-sc", "torch"]}
    versions["python"] = platform.python_version()
    (OUT / "package_versions.txt").write_text("\n".join(f"{k}={v}" for k, v in versions.items()) + "\n")
    (OUT / "readme.txt").write_text(
        "Latest-skill canonical CRC spatial validation.\n"
        "EPIfrac is normalized across the 11 epithelial subtype Tangram abundances within each spatial cell.\n"
        "Markers use saved post-annotation cell_subtype DEG tables, pvals_adj < 0.05, positive logFC, stable score-descending ordering, deduplication, and top 100 per subtype.\n"
        "Tangram uses 247 spatial-panel genes, mode=cells, cuda:0, 350 epochs, learning_rate=0.05, seed=42.\n"
        f"Completion audit: {'PASS' if passed else 'FAIL'}.\n"
    )
    if not passed:
        raise RuntimeError(checks_df.loc[~checks_df.passed].to_dict("records"))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
