---
name: project-epi-cm-core-workflow
description: Compact disease/project Epi-CM workflow covering single-cell integration, clustering/subtype annotation, epithelial-CM discovery, and optional spatial validation. Use when an agent should run or audit only these four core blocks instead of the larger numbered module suite.
---

# Project Epi-CM Core Workflow

Use this skill as a compact, self-contained workflow with four blocks only:

```text
01-celltype_integration_clustering
02-cell_subtype_integration_clustering
03-epi-cm-discovery
04-spatial-validation-optional
```

This skill is parallel to the larger modular skill set. It does not replace the
numbered module skills; it gives an agent a smaller canonical route when the
user asks for only integration, clustering, Epi-CM discovery, and optional
spatial validation.

## Non-Negotiable Rules

- Do not read or use `_to_download`, `old`, temporary, checkpoint, or example-only
  outputs unless the user explicitly names them as the current input.
- Do not invent data, random tables, pseudo-results, or substitute methods when
  an input is missing. Stop and ask for the missing input.
- Do not downsample or subsample cells unless the user explicitly requests that
  reduced-cell branch. If full data cannot run, record the blocker and ask.
- Do not delete, overwrite, move, or clear existing user outputs unless the user
  explicitly names the exact path and operation.
- Save figures as PDF and/or SVG only. Do not create PNG previews or temporary
  PNGs.
- Keep code, parameters, package versions, and readme/provenance tables with the
  outputs for every executed block.
- Use existing project labels and columns when present. Do not force a kidney,
  cancer, epithelial, CM, or status label into another dataset unless supported
  by its metadata and markers.
- For expression-based DEG and marker plotting, use `adata.raw` when present and
  the method supports it. Record the expression layer actually used.
- If a package/method has a GPU implementation and the local environment can use
  it, use it. If GPU fails from OOM, record OOM and use the CPU equivalent for
  that step; for non-OOM GPU failures, repair compatibility before falling back.

## Output Layout

Under the user-specified analysis root, create one compact workflow directory:

```text
<analysis_root>/epi-cm-core-workflow/
  codes/
  h5ad/
  tables/
  figures/
```

Use the four block names as secondary directories only where real outputs are
written:

```text
codes/01-celltype_integration_clustering/
h5ad/01-celltype_integration_clustering/
tables/01-celltype_integration_clustering/
figures/01-celltype_integration_clustering/

codes/02-cell_subtype_integration_clustering/
h5ad/02-cell_subtype_integration_clustering/
tables/02-cell_subtype_integration_clustering/
figures/02-cell_subtype_integration_clustering/

codes/03-epi-cm-discovery/
h5ad/03-epi-cm-discovery/
tables/03-epi-cm-discovery/
figures/03-epi-cm-discovery/

codes/04-spatial-validation-optional/
h5ad/04-spatial-validation-optional/
tables/04-spatial-validation-optional/
figures/04-spatial-validation-optional/
```

Do not pre-create empty subdirectories. A directory should exist only if it
receives a real output. Add `tables/<block>/readme.txt` for each executed block.

## Pre-Run Plan

Before long-running work, state a short plan:

```text
goal / expected result
main inputs
which of the four blocks will run or be skipped
key parameters or default choices
expected h5ad/tables/figures
validation checks
```

For inspection-only requests, one or two sentences are enough.

## 01-Celltype Integration Clustering

This block is the former Module 01, copied as the canonical instruction source
for the compact workflow. Inside this big skill, write this block's outputs under
`epi-cm-core-workflow/{codes,h5ad,tables,figures}/01-celltype_integration_clustering/`
unless the user explicitly asks to use the original numbered module output tree.
Do not use the shorter summary as a substitute for these copied rules.

# 01-project-singlecell-integration

Use this skill for sample merge, QC, normalization, dimensional reduction, batch correction, broad clustering, and major cell-type annotation. Detailed lineage-specific subtype annotation belongs in Module 02. Do not write lineage/subtype reclustering outputs into the Module 01 directory.

Project-context rule: this skill set is for disease CM-lineage analysis across replaceable project contexts. For any new disease, organ, species, cohort, or platform, set the project label, disease context, tissue/status vocabulary, broad-cell marker panel, subtype marker panel, sample metadata fields, and downstream biological focus from the user-provided or dataset-supported context. Do not force organ-specific, disease-specific, status-specific, lineage-specific, or CM-specific labels into a project unless the data and user goal support them.

## Pre-Execution Plan Requirement

Before executing code from this skill, write a concise method-and-result plan that the user can review and copy as the goal. Keep it result-oriented rather than overly procedural. Include only:

```text
analysis goal / expected result
method route to use
main inputs or provided intermediates
major code modules to run or skip
expected output figures/tables
key validation criterion
```

Do not start long-running analysis, dependency installation, or file-rewriting steps until this short plan has been stated. For simple inspection-only tasks, one or two sentences are enough.

If the user does not provide a manual choice for parameters, thresholds, method options, output naming, or optional branches, use the documented default settings in this skill and state that the default was used.

## Project Organization and Figure Output Contract

Treat the active working directory or user-specified analysis root as the output boundary for a project run, then create one module output directory inside it for each module. Module output directory names must use a replaceable project slug: `01-<project_slug>-singlecell-integration/`, `02-<project_slug>-cell-annotation/`, and so on. Derive `<project_slug>` from the user-provided project, disease, cancer type, organ, or cohort label using lowercase ASCII letters/numbers/underscores only; for example, a BRCA project uses `brca`, creating `01-brca-singlecell-integration/`. If the user does not provide a slug, ask once or use a sanitized slug from the active workdir/project name and record it. Do not create run-output directories inside the reusable skill source directory; that directory is only the reusable instruction/code source. Inside the module output directory, use one shared four-directory layout: `figures/`, `tables/`, `codes/`, and `h5ad/`. Use a category-first layout: create the four top-level category directories first. Create secondary-task subdirectories inside those category directories only when that category will immediately receive at least one real output file for that task. Do not pre-create empty secondary-task subdirectories under `figures/`, `tables/`, `codes/`, or `h5ad/` as placeholders, progress markers, or mirrored layout stubs. A directory creation command for a secondary task or candidate must be coupled to writing a real output file there; if the output is not generated, do not leave that task directory behind. Secondary tasks are logical analysis units inside that module run, not top-level output folders. This is a write-location rule, not a read restriction: a module run may read/reuse files and already generated outputs from other modules as inputs, but newly generated outputs for the current run must be written inside the active workdir's module output directory. Across the whole skill workflow, an agent must not delete, clear, overwrite, or move any existing output file or directory anywhere unless the user explicitly names the exact path and operation. During normal module execution, do not write, move, overwrite, clear, or delete files or directories inside any other numbered module output directory. Also do not delete, clear, overwrite, or move any existing output file or directory inside the current module unless the user explicitly names the exact path and operation. If a new result would conflict with an existing output, write to a new versioned path or stop and ask. Deleting any output directory is never part of a module run; it requires a separate explicit cleanup request naming the exact path. For example, when the workdir is `<analysis_root>` and the project slug is `brca`, Module 01 outputs go under `<analysis_root>/01-brca-singlecell-integration/{figures,tables,codes,h5ad}/<secondary-task>/`, not directly under `<analysis_root>/{figures,tables,codes,h5ad}/`, and not under the reusable skill source directory. Unless explicitly stated otherwise, relative output paths in this skill such as `h5ad/03-qc/adata_qc.h5ad` are relative to the active project-slugged module output directory, not to the skill source directory and not to the workdir root.

Use stable numbered secondary-task names that describe the analysis step, lineage, method, or figure group. Use the same secondary-task name under each shared `figures/`, `tables/`, `codes/`, and `h5ad/` category directory that receives an output from the same analysis, so files from that analysis stay aligned without creating empty placeholder directories. For example:

```text
<workdir>/01-<project_slug>-singlecell-integration/
  codes/
    03-qc/
      02_qc.py
    05-clustering-parameter-search/
      run_pcs10_nn10to40_res0p2to0p5.py
  h5ad/
    03-qc/
      adata_qc.h5ad
    05-clustering-parameter-search/
      selected/
        adata_inte.h5ad
  figures/
    03-qc/
      qc_violin.pdf
    05-clustering-parameter-search/
      pcs-10_nn-20_res-0p5/
        umap_leiden_res0p5.pdf
  tables/
    03-qc/
      qc_report.csv
    05-clustering-parameter-search/
      pcs-10_nn-20_res-0p5/
        clustering_parameters.csv
```

By default, save executable/reproducibility code under the module output directory's shared `codes/<secondary-task>/`, using ordered names such as `01_read_merge.ipynb`, `02_qc.py`, or `03_integrate.R`. Save AnnData-like objects under `h5ad/<secondary-task>/` as `.h5ad`, `.loom`, `.rds`, or equivalent files with stable names. Save corresponding figure files under `figures/<secondary-task>/` and use ordered names such as `01_umap.pdf` or `02_marker_dotplot.svg`. Save text-like and tabular outputs under `tables/<secondary-task>/`, such as CSV/TSV/XLSX/TXT/JSON/YAML logs, manifests, reports, mapping files, and parameter records. `figures/` should contain figure files only. `tables/` should contain text-like and tabular outputs only. `codes/` should contain executable/reproducibility code only. `h5ad/` should contain AnnData-like/intermediate object files only. Add `tables/<secondary-task>/readme.txt` documenting the input files, including any cross-module input/output files that were read, code order, h5ad/loom/rds objects, output figures/tables, and any skipped optional branches.

Do not write new h5ad, code, figures, or tables directly into the working-directory root, directly into the module output-directory root, or into the skill source root. Executable run outputs should live under the module output directory's shared four category directories. If a simple task has only one natural step, still use a small secondary-task subdirectory such as `01-main`, `01-qc`, or `01-epithelial-subclustering`, but create it only under category directories that receive real output files.

If one analysis step outputs multiple files or figures, put that output set in the same named secondary-task subdirectory under `figures/`, `tables/`, `codes/`, or `h5ad/`, using the same analysis prefix when possible.

If an output already exists, do not rerun only to recreate it in the new layout. Do not move or delete existing outputs for layout cleanup unless the user explicitly names the exact path and operation. Prefer to leave existing outputs in place, copy them into the organized location only when provenance is recorded, then update the corresponding code paths so future runs write to the same organized location.

When a task creates a run, lineage, candidate parameter set, or method variant, create matching candidate subdirectories inside the active secondary-task subdirectory only under category parents that receive output files for that candidate. For multi-candidate or multi-condition runs, use matching candidate names under each relevant parent when needed; for example, create `figures/05-clustering-parameter-search/pcs-20_nn-30_res-0p5/` only if figures will be saved and create `tables/05-clustering-parameter-search/pcs-20_nn-30_res-0p5/` only if parameter or cluster-count tables will be saved. A directory name containing a single resolution, such as `pcs-25_nn-30_res-0p3`, must contain only outputs for that exact resolution. If one directory intentionally stores multiple resolutions for the same graph, do not name it after one resolution; use an explicit aggregate name such as `pcs-25_nn-30_res-all` or an explicit range such as `pcs-25_nn-30_res-0p2-0p4`. Inside such aggregate directories, every output filename and every table row must still include the exact algorithm and resolution, such as `louvain_res0p2`, `louvain_res0p3`, or `leiden_res0p4`. For Leiden clustering parameter searches in this module, candidate outputs are intentionally lightweight: do not create h5ad candidate directories and do not save per-candidate AnnData objects. Keep candidate code files under the shared `codes/` with matching secondary-task and parameter-coded paths when code is emitted. Create `h5ad/05-clustering-parameter-search/selected/` only after the user selects final parameters and the selected graph/UMAP/Leiden run is rerun from the clean upstream Harmony object. After the user has reviewed the grid, selected final parameters, and the selected run has been rerun, redrawn, and saved, keep the candidate grid figures/tables/code in place. Rerunning or redrawing the selected/final clustering is not a cleanup request and must not delete grid-search outputs. Delete grid-search outputs only if the user gives a separate explicit cleanup request naming the exact path or output group to remove. Do not delete the final selected h5ad, selected-parameter record, broad-annotation inputs, candidate manifests, or any non-grid outputs during any cleanup.

Module 01 uses these secondary tasks by default:

```text
01-raw-count-generation        optional Cell Ranger count generation from FASTQ
02-merge-metadata              sample reading, per-dataset merge, cross-dataset merge, metadata recovery
03-qc                          QC metric calculation, Scrublet/doublet filtering, final cell/gene filtering
04-integration-harmony         normalize/log, HVG, raw preservation, regression, scaling, PCA, Harmony
05-clustering-parameter-search neighbor/PC/resolution candidate runs, UMAPs, manual broad Leiden parameter choice
06-broad-annotation            raw Leiden DEG export, broad marker review, leiden_coarse/cell_type annotation
07-score-rank-qc               broad score_genes/rank evidence and consistent-cell filtered handoff
08-unintegrated-diagnostics    optional unintegrated UMAP/clustering and batch-effect diagnostics
```

Optional backup integration task: `04-integration-scvi`. This task is not part of the default path. Create it only when the user explicitly requests scVI, supplies a valid scVI-integrated object, or approves scVI after the default Harmony route cannot complete. Keep it parallel to `04-integration-harmony`; never overwrite Harmony outputs with scVI outputs.

Broad clustering and broad annotation have separate column-name responsibilities. The clustering task must save only raw technical cluster labels such as `leiden_res0p5` or `louvain_res0p5`. The column `leiden_coarse` is created only after broad marker/DEG annotation is completed in `06-broad-annotation`; it is the completed broad cell-type annotation column, not a clustering output column and not a synonym for the selected Leiden/Louvain result.

Keep each secondary task aligned across the module output directory's relevant category directories that actually contain outputs for that task. For example, the broad annotation code goes under `codes/06-broad-annotation/`, its h5ad output under `h5ad/06-broad-annotation/`, its UMAP/dotplot output under `figures/06-broad-annotation/`, and its DEG/parameter tables under `tables/06-broad-annotation/`. If a task has no output in one of these categories, do not create an empty secondary-task directory there.

Required Module 01 h5ad outputs:

```text
h5ad/02-merge-metadata/adata_merge.h5ad
  merged pre-QC atlas; must preserve obs['sample'], obs['series'], obs['status'], obs['original_barcode'], unique obs_names, and raw/count expression.

h5ad/03-qc/adata_qc.h5ad
  full-cell, post-QC atlas; must contain QC metric columns and doublet columns when doublet detection ran.

h5ad/04-integration-harmony/adata_harmony.h5ad
  post-normalize/log/HVG/regress/scale/PCA/Harmony output; must contain adata.raw, obsm['X_pca'], and obsm['X_pca_inte'] or the documented adjusted basis. Do not add graph/UMAP/Leiden outputs in this object unless they are inherited non-final diagnostics.

h5ad/05-clustering-parameter-search/selected/adata_inte.h5ad
  do not create this during parameter search, but after the user manually chooses final broad-clustering parameters this selected Leiden-completed object is required. Rerun graph construction, UMAP, and Leiden from a clean read of h5ad/04-integration-harmony/adata_harmony.h5ad using the selected parameters, then save exactly this one selected Leiden-completed object. This selected clustering object must contain the selected raw cluster column such as leiden_res<RES> or louvain_res<RES>, but it must not create or rename that column to leiden_coarse. Do not save per-candidate h5ad files during parameter search.

tables/05-clustering-parameter-search/selected_clustering.csv
  create this only after the user manually chooses a completed candidate for broad annotation; it must contain selected_h5ad = h5ad/05-clustering-parameter-search/selected/adata_inte.h5ad, final raw Leiden/Louvain cluster column, n_pcs, n_neighbors, resolution, and manual-selection note. The final raw cluster column must not be leiden_coarse.

h5ad/06-broad-annotation/adata_anno.h5ad
  broad-annotated atlas; must contain obs['leiden_coarse'] and obs['cell_type'] initialized from it when required.

h5ad/07-score-rank-qc/adata_anno_score_genes_rank.h5ad
  unfiltered scored broad-annotation audit object.

h5ad/07-score-rank-qc/adata_anno_score_genes_rank_consistent.h5ad
  default downstream handoff after best-rank/label consistency filtering.
```

AnnData naming rule: h5ad filenames must use the current stable names. Use names such as `adata_anno.h5ad`, `adata_anno_score_genes_rank.h5ad`, and `adata_anno_score_genes_rank_consistent.h5ad`.

Optional h5ad outputs are allowed only when their branch is actually run, for example `h5ad/08-unintegrated-diagnostics/adata_unintegrated_diagnostics.h5ad` for unintegrated comparison. Do not save h5ad files directly in the Module 01 output-directory root.

When the optional scVI branch is approved and run, save its object as `h5ad/04-integration-scvi/adata_scvi.h5ad`. The object must contain `obsm["X_scVI"]`, the selected count layer used by scVI, and the original sample/series/status/original_barcode metadata. Save the trained model under `h5ad/04-integration-scvi/scvi_model/` or, if the model is external, write the exact model path to `tables/04-integration-scvi/scvi_model_path.txt`. Do not create `h5ad/04-integration-harmony/adata_harmony.h5ad` from scVI output and do not label `X_scVI` as Harmony.

QC-to-Harmony parameter-table contract:

From QC through Harmony, every executed step must write its own one-row-or-long-format CSV under the matching secondary-task `tables/` directory. Do not rely only on one merged notebook output or console log. Every table should include `step`, `input_h5ad`, `output_h5ad_or_object`, `n_obs_before`, `n_obs_after`, `n_vars_before`, `n_vars_after`, exact parameters, backend/package, code file, random seed when relevant, and notes/fallbacks. Required tables:

```text
tables/03-qc/01_initial_cell_filter_parameters.csv
tables/03-qc/02_scrublet_parameters.csv
tables/03-qc/02_scrublet_failed_sample_exclusions.csv
tables/03-qc/03_doublet_filter_parameters.csv
tables/03-qc/04_mt_ribo_qc_metric_parameters.csv
tables/03-qc/05_gene_filter_parameters.csv
tables/03-qc/06_final_cell_filter_parameters.csv
tables/03-qc/qc_report.csv

tables/04-integration-harmony/01_normalize_total_parameters.csv
tables/04-integration-harmony/02_log1p_parameters.csv
tables/04-integration-harmony/03_highly_variable_genes_parameters.csv
tables/04-integration-harmony/04_raw_assignment_and_hvg_subset_record.csv
tables/04-integration-harmony/05_regress_out_parameters.csv
tables/04-integration-harmony/06_scale_parameters.csv
tables/04-integration-harmony/07_pca_parameters.csv
tables/04-integration-harmony/08_harmony_integrate_parameters.csv
tables/04-integration-harmony/integration_harmony_summary.csv
```

If the optional scVI branch is executed, write separate scVI parameter/provenance tables under `tables/04-integration-scvi/`:

```text
tables/04-integration-scvi/01_scvi_input_counts_audit.csv
tables/04-integration-scvi/02_scvi_hvg_or_gene_selection_parameters.csv
tables/04-integration-scvi/03_scvi_setup_parameters.csv
tables/04-integration-scvi/04_scvi_train_parameters.csv
tables/04-integration-scvi/05_scvi_latent_and_model_paths.csv
tables/04-integration-scvi/integration_scvi_summary.csv
tables/04-integration-scvi/package_versions.txt
```

The QC parameter values remain dataset-selected, not hard-coded. The requirement is that each chosen value and each resulting cell/gene count change is recorded in the corresponding table.

Each analysis that produces an output should have corresponding source code under the current module's `codes/<secondary-task>/`. Acceptable code artifacts include `.ipynb`, `.py`, `.R`, and `.sh`, depending on the language actually used. Do not leave a figure, table, or exported result that can only be traced to manual GUI editing. If an analysis uses Python, keep the notebook and/or `.py` script that generates it; if it uses R, keep the `.R` script or R notebook; if both languages are used, keep both code artifacts under `codes/<secondary-task>/` with clear ordered prefixes. When converting notebooks to upload/download versions, keep the executable cells needed to reproduce the outputs and remove stale display output only when requested.

Each executed run should also create or update a parameter/provenance report under `tables/`, such as `tables/run_parameters.txt`, `tables/run_parameters.csv`, or a step-specific report in the same output subdirectory. The report should list the code file used, input files/objects, output files, exact parameters, random seeds, selected candidate/final settings, skipped steps, fallback decisions, and any user-approved method changes. Separate parameter changes made for biological/clustering quality from parameter changes made for resource limits such as GPU memory, CPU memory, runtime, or package failure. Record the tuning purpose, affected step, old/new setting, reason, and whether the final output is a biological final result or a resource-constrained fallback result.

Do not substitute another analysis method, algorithm, statistical test, visualization strategy, database, or input layer without explicit user permission. A missing software environment or missing package is not a reason to switch to an already installed alternative; first install or repair the environment for the target method. If the specified method cannot run after installation/repair attempts, stop that module, document the blocker in `tables/<secondary-task>/readme.txt`, and ask for confirmation before using any alternative. Any approved or documented method change should state why the original method was unsuitable or failed and why the replacement method is appropriate for the same analysis goal.

When a task, notebook run, script run, or long interactive kernel finishes, promptly close the process/kernel/session and release CPU memory and GPU memory. Do not leave idle Python, R, Jupyter, CellChat, RAPIDS, PyTorch, TensorFlow, or CUDA processes holding RAM/VRAM after the requested work is complete.

After each secondary task finishes, create or update `tables/<secondary-task>/package_versions.txt` describing the packages and tools used by that secondary task. Include Python packages, R packages, command-line tools, CUDA/GPU libraries when relevant, interpreter/R version, environment name or path, and the code files that used them. Optionally keep a module-level index at `tables/package_versions.txt` that points to the secondary-task reports, but do not replace the secondary-task reports with only a root-level file.

Install missing dependencies when they are required to execute the specified method or its approved acceleration path. This includes installing a compatible GPU-accelerated implementation when the method supports it and the machine has a usable GPU/CUDA driver, for example installing `rapids-singlecell`/RAPIDS to run Scanpy-style preprocessing through `rsc`. Dependency installation is allowed to make the requested method work; method substitution is not allowed without explicit user permission. Do not use an already installed alternative package or method just because the target environment is missing. For packages or methods that already provide GPU acceleration, enable and use the GPU-accelerated path after installing any missing compatible GPU packages and verifying imports/minimal execution. If the expected `rsc`/GPU path is installed but broken, incompatible, stalled, hung, or repeatedly fails for reasons other than GPU OOM, including a CUDA/CUDA-tag mismatch such as `cu11` vs `cu12` wheels, CuPy/RAPIDS/PyTorch wheels incompatible with the visible driver, missing CUDA runtime libraries, `libucx`/UCX errors, or `cuCtxGetDevice`/CUDA context errors, first diagnose, repair, reinstall, clean stale processes, release VRAM/RAM when safe, and retry the same `rsc`/GPU method without changing the requested method. Choose a compatible wheel, channel, or uv environment automatically from `nvidia-smi`, Python version, platform, and package compatibility information; do not ask the user to choose the CUDA tag. Ask the user only before system-driver changes, OS package changes that require elevated privileges, deleting an existing environment, or replacing a working environment used by other analyses. If GPU runs out of memory, inspect active GPU processes, close only stale or idle processes left by previous tasks/kernels when they can be safely identified, release VRAM, record the OOM, and switch that OOM step directly to the equivalent CPU implementation; do not repair, reinstall, or repeatedly retry GPU solely for OOM. Do not terminate unrelated active user processes unless the user explicitly approves. If no usable GPU is present, document the reason and use the equivalent `scanpy`/`sc` CPU implementation for the same method. If a usable GPU is present but the compatible `rsc`/GPU path still cannot run after non-OOM repair and cleanup attempts, stop that step, document the blocker in `tables/<secondary-task>/readme.txt` and `tables/<secondary-task>/package_versions.txt`, and ask the user before using a CPU fallback. GPU OOM after VRAM release is a pre-approved CPU fallback and does not require separate user confirmation. Cell sampling/downsampling is not a valid fallback for Module 01.

This GPU backend rule applies to all GPU-capable code in every module, not only Module 01. If the first full execution of the required task completes successfully with the planned backends, do not rerun only to validate the backend plan. If a GPU-accelerated step fails from GPU OOM, release VRAM, record the OOM, and treat CPU fallback for that step as pre-approved. If a GPU-accelerated step fails for a non-OOM reason and the user approves CPU fallback after documented repair/cleanup attempts, treat fallback at step granularity: run the equivalent `scanpy`/`sc` CPU operation only for the failed step when needed, then continue later steps with `rsc`/GPU whenever those later steps have a valid GPU path. Do not mark the entire remaining workflow CPU-only because one GPU step failed. Record a backend capability table under the relevant `tables/<secondary-task>/` directory with one row per executed step, including `step`, `planned_backend`, `attempted_backend`, `status`, `error_summary`, `fallback_backend`, `clean_input_reloaded`, and `final_backend_for_rerun`. Also export a human-readable backend summary as `tables/<secondary-task>/gpu_backend_capability_summary.csv` and, when useful, `tables/<secondary-task>/gpu_backend_capability_summary.txt`; these files must state which steps can use GPU/`rsc`, which steps must use CPU/`sc`, and the reason for each CPU step. When any GPU failure, fallback, or partial object mutation occurs during this exploratory/profiling pass, finish the required task only to learn which steps can use GPU, then start a fresh Python/R process, reload the nearest clean upstream input h5ad or original input files from disk, and rerun the whole required task once using the recorded final backend plan. In that final rerun, every step marked GPU-capable must use GPU/`rsc`, except steps with recorded GPU OOM or approved non-OOM CPU fallback, which should use the recorded CPU/`sc` fallback so the final run avoids repeating known mid-run GPU failures while still using GPU wherever it works. This clean rerun is required to guarantee object and variable purity. Do not reuse in-memory `adata`, views, layers, arrays, GPU buffers, fitted models, neighbor graphs, or partially written `.uns`/`.obsm`/`.var`/`.obs` fields from the profiling pass. Do not present outputs from a partial-failure/profiling pass as canonical final outputs. For mutating AnnData steps such as normalize/log/HVG/regress/scale/PCA/Harmony/neighbors/UMAP/Leiden, switching from `rsc` to `sc` after failure must start from the nearest saved clean upstream object, not from a possibly half-mutated in-memory object. Preserve existing outputs; if final rerun output paths conflict, write versioned paths or stop and ask.

Python package management rule: use `uv` by default for Python dependencies, and save uv-managed environments under a dedicated subdirectory in the total analysis directory for the dataset/project. Create the layout `uv_envs/<category>/.venv` under the top-level analysis root, where `<category>` is a stable dependency category such as `main`, `rapids`, `velocity_cellrank`, `cellchat_liana`, or `survival`. Use `uv_envs/main/.venv` for the default shared Python stack, and create another category only when dependency compatibility requires it. Install with `uv pip install --python uv_envs/<category>/.venv/bin/python --index-url https://pypi.tuna.tsinghua.edu.cn/simple ...` or activate that category `.venv` before `uv pip install --index-url https://pypi.tuna.tsinghua.edu.cn/simple ...`. Do not create a separate per-module environment unless it is also a documented dependency category, do not create or reuse a global environment, do not put the run environment inside the skill source directory, and do not run bare `pip install` into the system/user Python unless `uv` is unavailable or the user explicitly requests it. Keep any category-level `pyproject.toml`, `uv.lock`, or requirements export inside `uv_envs/<category>/`, and record all environment paths, categories, package versions, and mirror fallbacks under the relevant `tables/<secondary-task>/package_versions.txt`. Prefer TUNA mirrors for package downloads: use the TUNA PyPI mirror for `uv pip install`, and use TUNA CRAN/Bioconductor mirrors for R packages when practical. If the TUNA mirror is unavailable, stale, or missing a required package, fall back to the official source only for the affected dependency and record the mirror fallback in `tables/<secondary-task>/package_versions.txt`. If `uv` is unavailable, install `uv` first when possible; otherwise document the fallback package manager in `tables/<secondary-task>/package_versions.txt`. R packages, Cell Ranger, velocyto, CUDA drivers, and system libraries are outside `uv` and should be installed with their appropriate manager while still recording versions.

For this module, RAPIDS/`rapids-singlecell` is a hard preflight gate before CPU Scanpy fallback. Before running the QC or integration path with `scanpy/sc` CPU equivalents, create and run a small preflight check that records GPU visibility, `rapids_singlecell` import status, a minimal AnnData GPU transfer test, and CUDA/CuPy/RAPIDS versions in the relevant task table directory. If the preflight report shows a visible GPU but `rapids_singlecell` is not importable, install or repair a compatible RAPIDS/`rapids-singlecell` environment, rerun the preflight, and only then continue. CUDA tag/version problems, for example wrong `cu11`/`cu12` wheels, incompatible CuPy/RAPIDS packages, missing CUDA runtime libraries, `libucx`/UCX errors, or `cuCtxGetDevice`/CUDA context errors, are part of this automatic environment-repair step. Select the compatible package set from the visible driver and Python version without asking the user to choose the CUDA tag, unless the repair requires system-driver changes, elevated OS package changes, deleting an existing environment, or replacing a working shared environment. If `rsc` starts but OOMs, release VRAM, record the OOM, and switch that step directly to `scanpy/sc` CPU fallback without GPU repair or repeated GPU retries. If `rsc` stalls, hangs, or fails for a non-OOM reason, diagnose the GPU/process/environment issue, clean only safely identifiable stale processes, release memory, and retry `rsc`. Direct CPU `scanpy/sc` fallback is allowed without an install/repair attempt when no usable GPU is visible or when the attempted GPU step OOMs after VRAM release. When a usable GPU is visible, CPU fallback after non-OOM `rsc` failure requires documented repair attempts and explicit user approval.

If CPU fallback is approved for one `rsc` step, continue to test later GPU-capable steps with `rsc` during the same required workflow instead of converting all downstream steps to CPU. If the mixed backend workflow completes after any fallback, use it as backend profiling, close that process/kernel, reload the clean upstream input from disk in a fresh process, and rerun the final workflow with the recorded per-step backend plan. If the original backend workflow completes once without GPU failure or fallback, no additional rerun is required.

Do not create a symlink in the current directory that points to an older output just to make it look renamed. If an output needs a new name, either regenerate it in the correct output directory or make a real copied file with documented provenance.

Python random seed rule: every Python script and notebook must define and use the default fixed random seed near the top of the file: `SEED = 42`, unless the user explicitly specifies another seed. Set `random.seed(SEED)` and `numpy.random.seed(SEED)` when those libraries are used, and set framework-specific seeds for stochastic packages when relevant, such as PyTorch, TensorFlow, scvi-tools, scVelo, veloVI, CellRank, scikit-learn, Scanpy, or UMAP. Pass `random_state=SEED`, `seed=SEED`, or the package-equivalent argument to every function that supports it, including PCA, neighbors/UMAP, Leiden/clustering, train/test splits, model fitting, bootstrapping, permutations, and plotting layouts when applicable. If a function has no seed argument or still has nondeterministic GPU behavior, document that limitation. Record the fixed seed in the relevant secondary-task parameter tables, per-candidate parameter tables, and `tables/<secondary-task>/package_versions.txt`; do not change the seed between candidate groups unless the user explicitly requests seed sensitivity testing.

Expression-layer rule: for expression-based plotting, marker visualization, dotplots, violin plots, gene scoring, or other Scanpy-compatible calculations, if the function exposes a `use_raw` argument and `adata.raw` is present, set `use_raw=True` by default unless the user explicitly requests a different layer or the method requires counts. Record any exception and the expression layer actually used.

Default figure formats are PDF and/or SVG only. Do not create, save, convert, or request PNG files for final, intermediate, diagnostic, preview, thumbnail, or temporary figure outputs. If a tool defaults to PNG, override it to PDF/SVG or stop and ask; do not leave `.png` files under `figures/`, `tables/`, `codes/`, or `h5ad/`.

Global figure style for every module:
- Use the module-specific `Module Figure Style Contract` in this SKILL.md as the plotting-style source of truth. Do not rely on any material outside this SKILL.md to infer final figure style. Keep exactly one current canonical plotting route for each final figure, and label any other route as non-canonical or exploratory; do not duplicate final PDFs/SVGs for the same panel or file stem.
- Apply the style for the corresponding analysis type and module-specific section first; use these global rules only as baseline guardrails. Do not copy a figure style from an unrelated analysis type just because the file format or package is similar.
- Save figure outputs as PDF and/or SVG only; never create PNG previews, thumbnails, diagnostics, or temporary figure files.
- Keep text editable whenever the backend supports it. In Matplotlib set `pdf.fonttype = 42`, `ps.fonttype = 42`, `svg.fonttype = "none"`, and use a standard sans-serif font such as DejaVu Sans or Arial.
- Keep axis ticks visible on quantitative and categorical plots, including heatmaps, dotplots, barplots, forest plots, scatter plots, survival plots, and UMAP panels with axes. Do not call `axis("off")`, remove tick labels, or hide spines unless the plot type is a pure network/chord/graph layout where axes have no coordinate meaning.
- Use clean white backgrounds, black axis text, readable tick labels, and legends/colorbars outside or to the right when practical. Rotate dense x-axis labels rather than letting them overlap.
- Prefer the official plotting interface for the relevant package before manual low-level drawing. Use manual Matplotlib/ggplot2/Seaborn layout code only when the package interface cannot express the required final figure or when this skill gives explicit canonical manual code; record that reason in the code or parameter log.
- Size every figure element for the exported canvas: labels, legends, colorbars, risk tables, titles, p/q labels, arrows, node labels, and panel titles must be readable and must not collide. If any element overlaps or is clipped in the generated PDF/SVG, increase canvas size, margins, row/column spacing, legend placement, font size, label wrapping, or panel spacing and rerender before considering the figure complete.
- For UMAP-like embedding panels, keep panels square, use the official Scanpy/scVelo plotting interface when possible, keep framed axes, and use consistent palettes for the same labels across full-atlas and lineage-specific plots.
- For heatmaps, use a colorbar with visible ticks. Use centered diverging palettes for signed correlations/effects, sequential palettes for nonnegative scores or transition weights, and method-specific labels in titles/captions.
- For dotplots, keep gene and group orders explicit, preserve tick labels, and use standard scaling only when the method requires display scaling.
- For forest/survival plots, keep confidence intervals, hazard-ratio/reference lines, p/q labels, and axis ticks visible.
- For network/chord/directed-graph plots, use PDF/SVG, editable labels, clear legends/colorbars when edge weights are encoded, and document any intentional axis removal.


Scanpy plotting hard rule: For ordinary `sc.pl.*` outputs that support `save=...`, do not create Matplotlib axes, do not pass `ax=...`, and do not save with `fig.savefig` or `ax.figure.savefig`. Use `sc.settings.figdir` plus the Scanpy `save` argument. For multi-panel Scanpy plots, use the package interface such as `color=[...]`, `ncols`, `wspace`, or `standard_scale` instead of manual subplots. Manual `ax=...` is allowed only for documented special overlays or incompatible per-panel settings that Scanpy cannot express; record the reason in the run-parameter table or a code comment.


For both Python and R plotting, use the official plotting interface of the relevant package by default, such as Scanpy, Matplotlib, Seaborn, CellRank, scVelo, ggplot2, ComplexHeatmap, or CellChat plotting APIs. For Python-generated single-panel plots, use a default square canvas of 2.5 x 2.5 inches unless the user specifies another size or the plot type clearly requires more space, such as multi-panel layouts, heatmaps, wide dotplots, survival plots, or network/chord diagrams. For UMAP plots, use the default Scanpy save path and keep the call minimal. Set figure parameters with `sc.set_figure_params(figsize=(3, 3), dpi=150)` or `sc.settings.set_figure_params(figsize=(3, 3), dpi=150)`, set `sc.settings.figdir` to the target output directory, then call `sc.pl.umap(adata, color="leiden_coarse", save="_name.pdf")`. Keep the default Scanpy-style framed axes and outside legends; do not manually create `plt.subplots`, pass `ax=...`, or call `fig.savefig` for ordinary UMAPs. Use `ncols` only when `color` contains multiple objects, for example `sc.pl.umap(adata, color=["cell_subtype", "status", "cnv_score"], ncols=3, wspace=0.4, save="_celltype_status_cnvscore.pdf")`. Do not use `ncols` for a single-color UMAP. Keep `save` as a suffix/name handled by Scanpy rather than a full path. Manual axes and `fig.savefig` are reserved for special cases where different panels require incompatible per-panel parameters or post-processing that the official `color=[...]` interface cannot express, and the reason must be documented because manual saving can greatly increase PDF/SVG size. For ordinary UMAPs, specifically avoid `return_fig=True` followed by `ax.savefig(..., bbox_inches="tight")`; on large atlases this can inflate files from sub-MB Scanpy-saved outputs to multi-MB PDFs or very large SVGs. When the task is only to inspect, audit, or explain existing UMAP code/output size, do not modify the source code or rerun plotting unless the user explicitly asks for a fix or rerender. Keep all figure text editable as text whenever the plotting backend supports it; do not convert labels, legends, tick labels, titles, or annotations to outlines/paths unless the user explicitly requests it. For multi-panel figures, especially multi-panel UMAPs, verify that each UMAP panel remains square, legends/colorbars/titles do not overlap, and adjacent panels do not collide. It is acceptable to adjust the default figure width, height, `wspace`, `hspace`, legend font size, or margins to prevent overlap and preserve square UMAP panels. Do not add extra custom titles to UMAP panels; use the Scanpy default title derived from `color` unless the user explicitly asks for custom titles. Do not add automatic bitmap/raster conversion rules; let the user decide figure-size tradeoffs from the actual output files. Do not draw sample-colored UMAPs as default final figures; generate sample-colored UMAPs only when the user asks for them or when they are needed as integration/batch-mixing diagnostics, and label them as diagnostic outputs.

## Module Figure Style Contract

Use the following single-cell integration figure styles unless the user
explicitly asks for a different style. Do not mention implementation provenance from prior runs in generated reusable code, figure labels, captions, or readme files.

- QC plots: keep axis ticks and sample labels visible, use separate panels or
  files for n_genes, counts, mitochondrial fraction, doublet score, and retained
  cell counts, and save the threshold table used to generate the figures.
- Integration/clustering UMAPs: use official Scanpy plotting with square framed
  panels, one candidate figure per parameter set or per documented grid layout,
  and stable palettes for the same label column. Do not save candidate h5ad
  objects during grid search unless the user selects the final parameters.
- Broad annotation UMAPs: plot raw cluster labels and `leiden_coarse` labels
  separately, using the same embedding and figure dimensions, so manual
  annotation changes are visually auditable.
- Broad marker dotplots: use `use_raw=True` when available,
  `standard_scale='var'`, explicit marker order, explicit broad-label order,
  visible gene/group ticks, and PDF/SVG output only.
- Score/rank consistency QC UMAPs: show predicted best label and final
  `leiden_coarse` on comparable square panels, save consistent and inconsistent
  cell summaries, and keep diagnostic plots clearly labeled as QC.

## Cross-Module Contract

The modular skills are not fully isolated pipelines. Each module can be entered independently when valid intermediates are supplied, but all modules should share the same project contracts:

```text
sample_id / obs['sample']
series
status
original_barcode
cell_type
cell_subtype
score matrix sample IDs
spatial sample IDs
bulk/survival sample IDs
```

Do not let each module invent its own sample names, subtype names, broad cell-type mapping, or score-matrix orientation. When a module starts from supplied intermediate files, first validate that these identifiers match the upstream/downstream modules or write an explicit mapping table.

## Gene Identifier Audit Contract

Before any AnnData object is used for normalization, HVG selection, DEG, marker scoring, pathway analysis, Tangram, CellChat, pySCENIC, inferCNV, or cross-dataset gene matching, audit whether `adata.var_names` contains gene IDs or gene names. Repeat the audit for `adata.raw.var_names` when `adata.raw` exists; do not assume `.X` and `.raw` use the same identifier type.

Inspect all of the following:

```text
adata.var_names
adata.var columns: gene_id, gene_ids, ensembl_id, feature_id
adata.var columns: gene_name, gene_names, gene_symbol, symbol, feature_name
adata.raw.var_names and adata.raw.var columns when adata.raw exists
```

Classify identifiers using the complete vector, not one example. Strip an Ensembl version suffix only for detection, for example `ENSG00000141510.18 -> ENSG00000141510`. Treat identifiers matching `^ENS[A-Z]*G[0-9]+(?:\.[0-9]+)?$` as Ensembl-like gene IDs. Record `var_names_type` as `gene_id` when at least 90% are Ensembl-like, `gene_name` when at most 10% are Ensembl-like and the values otherwise behave as gene symbols/names, and `mixed_or_unknown` otherwise. Also record total features, Ensembl-like count/fraction, unique count, duplicate count, available ID/name columns, missing values, and the first 20 non-sensitive examples.

Write the audit before downstream gene-based analysis:

```text
tables/<secondary-task>/gene_identifier_audit.csv
tables/<secondary-task>/gene_identifier_examples.csv
```

Use one row per audited object layer (`var_names` and, when present, `raw.var_names`) in the audit table. Include `input_h5ad`, `object_layer`, `var_names_type`, `n_features`, `n_unique`, `n_duplicates`, `ensembl_like_n`, `ensembl_like_fraction`, `candidate_gene_id_columns`, `candidate_gene_name_columns`, `selected_gene_id_source`, `selected_gene_name_source`, and `action_taken`.

Do not silently rename genes. If `var_names` contains gene IDs but a trusted gene-name column exists, preserve the original IDs in `var['gene_id_original']`, preserve the selected names in `var['gene_name']`, and record which column was used. If no trusted gene-name mapping exists and a downstream step requires gene symbols, stop and request or build a documented species/build-matched annotation mapping before continuing. Never map human and mouse annotations interchangeably. When mapping creates duplicate gene names, report the duplicates and keep a reversible ID-to-name mapping; do not silently resolve them with `var_names_make_unique()` or arbitrary suffixes. Any approved conversion must save `tables/<secondary-task>/gene_id_to_gene_name_mapping.csv` and document whether `var_names` was changed or gene names were retained only as metadata/display labels.

Canonical audit code structure:

```python
import re
import pandas as pd

ENSEMBL_GENE_RE = re.compile(r"^ENS[A-Z]*G[0-9]+(?:\.[0-9]+)?$")

def audit_gene_identifiers(index, object_layer):
    values = pd.Index(index).astype(str)
    is_ensembl = values.to_series(index=range(len(values))).str.match(ENSEMBL_GENE_RE)
    fraction = float(is_ensembl.mean()) if len(values) else 0.0
    identifier_type = "gene_id" if fraction >= 0.90 else (
        "gene_name" if fraction <= 0.10 else "mixed_or_unknown"
    )
    return {
        "object_layer": object_layer,
        "var_names_type": identifier_type,
        "n_features": len(values),
        "n_unique": int(values.nunique()),
        "n_duplicates": int(values.duplicated().sum()),
        "ensembl_like_n": int(is_ensembl.sum()),
        "ensembl_like_fraction": fraction,
    }
```

## Optional Raw Inputs

Cell Ranger and velocyto are optional. If the user already provides count matrices, h5ad, loom, or ladata, skip raw preprocessing and validate the supplied files.

If FASTQ is supplied and counts are missing:

```text
cellranger count -> sample outs/filtered_feature_bc_matrix
```

If velocity is required and loom/ladata are missing:

```text
velocyto run10x -> sample velocyto loom
```

## Sample Naming Contract

Use stable sample IDs from the start of single-cell reading and merging:

```text
if sample folders are GEO/GSM accessions -> use the GSM folder name as sample_id
if samples do not have GSM-style names -> use the expression matrix directory name as sample_id
obs['sample'] = sample_id
```

Apply this rule in the single-cell read/merge implementation generated for the current run. The expression matrix directory means the folder that contains the 10x matrix, the Cell Ranger `outs/filtered_feature_bc_matrix`, or the user-supplied per-sample matrix. Do not rename samples to generic labels such as `sample1` or `sample2` unless the original directories are unavailable. If synthetic names are unavoidable, write a reversible mapping table:

```text
sample_id
source_matrix_path
source_directory_name
sample_label
```

Use the same `sample_id` later for h5ad metadata, Cell Ranger output directories, velocyto loom matching, spatial/bulk joins, and downstream figure labels unless a documented mapping table says otherwise.

## Required Metadata After Merge

```text
sample
series
status
original_barcode
```

Barcode rule:

```text
obs['original_barcode'] = raw source barcode
obs['sample'] = source sample ID
obs_names = globally unique merged cell ID, usually sample + separator + original_barcode
```

Keep a reversible mapping:

```text
merged_cell_id, sample, original_barcode, source_matrix_path
```

Metadata recovery:

```text
If the user does not provide sample metadata, recover it from primary sources before merging.
Use GEO/SRA/BioProject metadata, the original paper, supplementary tables, and dataset README files.
Determine disease/project-relevant sample status when this information is available, such as primary lesion, adjacent/normal-like tissue, metastatic lesion, thrombus, mixed sample, treated sample, control sample, or excluded out-of-scope sample. Use organ-specific labels only when the metadata support them.
Do not assign status only from a file name unless no better source exists; if inferred from a file name, mark it as inferred.
Save a sample metadata table with sample, series, status, tissue_site, disease_context, treatment_context, include/exclude decision, source_url_or_file, and notes.
```

## Integration Workflow

Implement or run:

```text
1. Read all samples.
2. Preserve raw counts in adata.raw or layers['counts'].
3. Merge samples with stable sample/series/status metadata using the sample merge contract below.
4. Save the just-merged object before QC as `h5ad/02-merge-metadata/adata_merge.h5ad`; this object should include the raw barcode column, so do not create a separate `adata_merge_raw_barcode.h5ad` by default.
5. QC cells and genes, including mitochondrial-content metrics and doublet detection/removal when raw counts support it.
6. Save the QC-completed object.
7. Run the post-QC integration workflow.
8. Also keep unintegrated coordinates/clusters when useful for diagnostics.
```

## Sample Merge Contract

Project-style merging is hierarchical:

```text
1. Within each dataset/series, merge that dataset's samples first.
2. Store each cell's source sample in obs['sample']; derive it from the expression matrix directory name when no GSM/sample accession is available.
3. Store each dataset name in obs['series'] during the cross-dataset merge.
4. Then merge the per-dataset AnnData objects into one atlas.
5. Use an inner gene join for both within-dataset and cross-dataset merges. This is the default and should change only if the user explicitly asks for a different join strategy.
6. Save the merged object before QC as `h5ad/02-merge-metadata/adata_merge.h5ad`, including `obs['original_barcode']`.
```

This mirrors the project pattern:

```text
sample-level AnnData objects -> per-dataset AnnData -> all-dataset AnnData
```

The final merged atlas should keep only shared genes under the default inner-join strategy. If the user requests a different join strategy, document why and report how many genes are retained.

Do not merge all samples from all datasets in one flat operation by default. The default contract is dataset-internal sample merge first, then cross-dataset merge, using inner gene joins at both stages.

## QC Contract

If the user does not specify QC parameters, use the default QC protocol below. The protocol order is fixed by this skill, but numeric thresholds are not hard-coded. The agent must inspect the current dataset's QC distributions, choose dataset-appropriate threshold values, record those values, then run the bundled QC executor. If the user provides a QC-completed h5ad, this step can be skipped after validating that doublets and low-quality cells were already handled.

Default project-style QC follows this order:

```text
1. Filter cells using an agent-selected initial cell-level minimum gene threshold.
2. Run Scrublet independently for every sample using batch_key = 'sample'; immediately exclude any whole sample that does not complete with valid Scrublet outputs.
3. For retained samples with valid Scrublet outputs, remove cells with obs['predicted_doublet'] == True.
4. Calculate mitochondrial metrics using human gene prefix 'MT-'.
5. Calculate ribosomal metrics using gene prefix 'RPS'.
6. Filter genes using an agent-selected minimum cell-count threshold.
7. Filter genes using an agent-selected minimum count threshold.
8. Recalculate QC metrics after gene filtering before final cell filtering.
9. Keep cells using agent-selected final cell-level QC thresholds.
10. Save the cleaned object as `h5ad/03-qc/adata_qc.h5ad`.
```

At minimum, compute and store these QC columns:

```text
n_genes_by_counts
total_counts
pct_counts_MT
total_counts_MT
pct_counts_RIBO
total_counts_RIBO
doublet_score
predicted_doublet
```

Mitochondrial genes:

```text
human gene symbols commonly start with MT-
mouse gene symbols commonly start with mt-
```

Use the workflow above as the default QC contract. The active threshold values are chosen by the agent from the current dataset's QC distributions unless the user supplies values. Report the selected thresholds, selection rationale, and how many cells/genes are removed per sample. If a supplied h5ad already has QC columns and the user trusts them, validate and reuse them instead of recomputing.

QC implementation scope:

```text
Use Scanpy/rapids-singlecell equivalents for cell filtering, gene filtering, QC metric calculation, Scrublet doublet scoring, and final cell filtering.
Parameter values should be selected for the current dataset or supplied by the user, not copied from an unavailable example workflow.
Record the active QC parameter range and removed cell/gene counts in the run report.
```

The QC executor must implement the default QC protocol and requires explicit threshold arguments. Before running it, the agent should compute/plot QC distributions and decide:

```text
min_genes_initial
min_cells_gene
min_counts_gene
min_genes_final
max_genes_final
max_pct_mt
```

Function/API range:

```text
Scanpy: sc.pp.filter_cells, sc.pp.filter_genes, sc.pp.calculate_qc_metrics, sc.pl.violin
rapids-singlecell when available: rsc.get.anndata_to_GPU, rsc.pp.scrublet, rsc.get.anndata_to_CPU
Project helper function names when present: calculate_percent_mt, calculate_percent_rps, subset_anndata
AnnData operations: obs/var flag creation, boolean cell filtering, write_h5ad
Legal placeholder call forms:
sc.pp.filter_cells(adata, min_genes=min_genes_threshold)
sc.pp.filter_genes(adata, min_cells=min_cells_threshold)
sc.pp.calculate_qc_metrics(adata, qc_vars=qc_var_names, inplace=inplace_flag)
rsc.pp.scrublet(adata, batch_key=sample_key)
sc.pl.violin(adata, keys=qc_metric_names, groupby=sample_key)
```

Default Scrublet sequence:

```python
failed_scrublet_samples = []
scrublet_failure_rows = []

for sample_id in adata.obs[sample_key].astype(str).unique():
    sample_mask = adata.obs[sample_key].astype(str) == sample_id
    sample_adata = adata[sample_mask].copy()
    input_names = sample_adata.obs_names.copy()
    try:
        rsc.get.anndata_to_GPU(sample_adata)
        rsc.pp.scrublet(sample_adata)
        rsc.get.anndata_to_CPU(sample_adata)
        required = {"doublet_score", "predicted_doublet"}
        if not required.issubset(sample_adata.obs.columns):
            raise ValueError("missing required Scrublet output columns")
        if not sample_adata.obs_names.equals(input_names):
            raise ValueError("Scrublet output cell IDs do not match input cell IDs")
        if sample_adata.obs["doublet_score"].isna().any():
            raise ValueError("doublet_score contains missing values")
        if sample_adata.obs["predicted_doublet"].isna().any():
            raise ValueError("predicted_doublet contains missing values")
    except Exception as exc:
        failed_scrublet_samples.append(str(sample_id))
        scrublet_failure_rows.append({
            "sample": str(sample_id),
            "n_cells_excluded": int(sample_mask.sum()),
            "failure_reason": str(exc),
            "action": "exclude_entire_sample",
        })
        continue
    adata.obs.loc[sample_adata.obs_names, "doublet_score"] = sample_adata.obs["doublet_score"]
    adata.obs.loc[sample_adata.obs_names, "predicted_doublet"] = sample_adata.obs["predicted_doublet"]

if failed_scrublet_samples:
    failed_mask = adata.obs[sample_key].astype(str).isin(failed_scrublet_samples)
    adata = adata[~failed_mask].copy()

if adata.n_obs == 0:
    raise RuntimeError("No samples remain after excluding failed Scrublet samples.")
if adata.obs[["doublet_score", "predicted_doublet"]].isna().any().any():
    raise RuntimeError("Retained samples must all have complete Scrublet outputs.")

adata[~adata.obs["predicted_doublet"].astype(bool)].write_h5ad(scrublet_filtered_h5ad)
adata = adata[~adata.obs["predicted_doublet"].astype(bool)].copy()
```

Scrublet must run as one single-sample AnnData subset at a time with `rapids-singlecell`/`rsc` by default. Keep the full merged AnnData on CPU, process one sample subset at a time with `rsc.pp.scrublet`, move the subset result back to CPU, and write each sample's `doublet_score` and `predicted_doublet` values back into the original merged `adata.obs` by cell barcode/index. Do not run default Scrublet on the full merged object, do not downsample cells, do not drop samples only to reduce memory use, and do not treat sample-wise Scrublet as a different doublet method. A sample passes Scrublet only when the run completes, input/output cell IDs match, and every cell has non-missing `doublet_score` and `predicted_doublet` values. Once a sample is classified as failed after the documented same-method backend attempt/fallback route, immediately remove every cell from that sample from the working AnnData; do not retain it with missing calls, treat missing calls as singlets, rerun Scrublet on the full merged object, or ask whether the failed sample should be retained. Keep all samples with valid Scrublet calls and then apply the normal `predicted_doublet` cell filter. If no samples remain, stop QC. Record the sample order, per-sample input/output cell counts, backend, fallback attempts, failed sample IDs, failure reasons, excluded sample IDs, and excluded cell counts in `tables/03-qc/02_scrublet_parameters.csv` and `tables/03-qc/02_scrublet_failed_sample_exclusions.csv`.

Default QC metric sequence, run once immediately before `subset_anndata(adata)`:

```python
def calculate_percent_mt(adata, pattern=mt_gene_prefix):
    adata.var[mt_qc_var] = adata.var_names.str.startswith(pattern)
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=[mt_qc_var],
        percent_top=percent_top_value,
        log1p=log1p_flag,
        inplace=inplace_flag,
    )
    return adata

def calculate_percent_rps(adata, pattern=ribo_gene_prefix):
    adata.var[ribo_qc_var] = adata.var_names.str.startswith(pattern)
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=[ribo_qc_var],
        percent_top=percent_top_value,
        log1p=log1p_flag,
        inplace=inplace_flag,
    )
    return adata

adata = calculate_percent_mt(adata)
adata = calculate_percent_rps(adata)
sc.pl.violin(
    adata,
    qc_metric_names,
    jitter=violin_jitter,
    multi_panel=multi_panel_flag,
)
adata = subset_anndata(adata)
```

The `rsc` calls require `rapids-singlecell`. If a usable GPU/CUDA driver is present but `rapids-singlecell` is unavailable, install a compatible `rapids-singlecell`/RAPIDS stack in the active or module-specific environment, verify `import rapids_singlecell as rsc`, then run the default single-sample `rsc` Scrublet workflow. If single-sample `rsc` Scrublet OOMs, release VRAM, record the OOM, and switch that sample to the equivalent `scanpy`/`sc` CPU Scrublet fallback without repairing or repeatedly retrying GPU. If the single-sample `rsc` workflow stalls, hangs, or fails for a non-OOM reason, diagnose and repair the `rsc`/GPU issue, clean stale processes only when safely identifiable, release memory, and retry the single-sample `rsc` workflow. Use the equivalent `scanpy`/`sc` CPU Scrublet workflow without user approval only if no usable GPU is present or if the single-sample `rsc` attempt OOMs after VRAM release. If the permitted single-sample route still fails or returns invalid/incomplete calls, classify that sample as failed and directly exclude the entire sample; do not stop the other samples' Scrublet runs and do not keep the failed sample in `adata_qc`. Keep the same output columns (`doublet_score`, `predicted_doublet`) and filtering rule for retained samples.

Doublet handling:

```text
default method = Scrublet
default execution = one single-sample AnnData subset at a time with rsc.pp.scrublet
write each sample's doublet_score and predicted_doublet back to the merged adata.obs
for each single-sample AnnData subset, directly remove the whole sample if Scrublet does not complete with valid calls after the permitted same-method attempt/fallback route
never retain a failed sample by filling missing predicted_doublet values with False
continue processing the remaining samples; stop only if no samples remain
store scores/calls in obs
default cleaned object removes predicted_doublet cells
keep an unfiltered object or a filter mask when possible
```

Do not replace Scrublet with another doublet method unless the user explicitly asks for it or supplies an external doublet table. If a non-Scrublet source is user-approved, label it clearly and do not mix doublet calls from different methods without documenting it.

## Post-QC Integration Contract

Project-style post-QC integration follows the fixed workflow order below. The method sequence is fixed; numeric parameters can be tuned when needed and should be recorded.

```text
1. Read `h5ad/03-qc/adata_qc.h5ad`.
2. If using rapids-singlecell, move to GPU with rsc.get.anndata_to_GPU(adata).
3. Normalize total counts; use the package default or a common project-style value unless the user provides a manual choice.
4. Log-transform with log1p.
5. Select/mark highly variable genes; use the package default or a common project-style value unless the user provides a manual choice. Keep the selected HVG count at or above 2000 unless the user approves a smaller value, and record the actual HVG count.
6. Store adata.raw after normalize/log/HVG marking and before HVG filtering.
7. Filter/subset to HVGs.
8. Regress out total_counts and pct_counts_MT, including mitochondrial-content regression through pct_counts_MT.
9. Scale the matrix; use the package default or a common project-style value unless the user provides a manual choice.
10. Run PCA; use the package default or a common project-style value unless the user provides a manual choice.
11. Run Harmony integration with key = "sample", basis = "X_pca", adjusted_basis = "X_pca_inte".
12. Save the Harmony output as `h5ad/04-integration-harmony/adata_harmony.h5ad`. This is the only required handoff from QC-to-Harmony into clustering.
```

Use GPU acceleration when the requested package or method itself has a valid GPU-accelerated path. For single-cell preprocessing/integration steps implemented with RAPIDS/rapids-singlecell, create a task-local preflight check first and save its report under the relevant secondary task table directory. If GPU/CUDA is usable but `rapids-singlecell` is missing, install a compatible RAPIDS/`rapids-singlecell` stack in the active or module-specific environment, verify the import and a minimal GPU transfer, rerun the preflight, and then use RAPIDS GPU functions. If the usable GPU path fails because of a CUDA tag/version mismatch, wrong `cu11`/`cu12` wheel, incompatible CuPy/RAPIDS package, missing CUDA runtime library, `libucx`/UCX error, or `cuCtxGetDevice`/CUDA context error, repair or recreate a compatible GPU environment automatically and retry before asking about CPU fallback. If the requested step has no GPU path, use its normal CPU path. If a usable GPU is present and the compatible `rsc`/GPU path runs out of VRAM, release VRAM, record the OOM, and switch that step directly to the equivalent CPU implementation without repairing or repeatedly retrying GPU. If a usable GPU is present and the compatible `rsc`/GPU path fails, stalls, hangs, or lacks a required operation for a non-OOM reason, diagnose and repair the `rsc`/GPU issue, clean stale processes only when safely identifiable, release memory, and retry the same `rsc`/GPU method. Do not switch to CPU just because the GPU run is stuck or slow. CPU `scanpy/sc` fallback is allowed without user approval when no usable GPU is present or when the GPU step OOMs after VRAM release; with a usable GPU, non-OOM CPU fallback requires documented repair attempts and explicit user approval. A CPU fallback for one failed GPU step does not imply CPU fallback for later steps: continue later GPU-capable steps with `rsc` when possible, record the per-step backend result, then rerun the final workflow from clean input if any fallback or partial failure occurred. If the first end-to-end execution runs through without fallback, keep it and do not rerun just for backend validation. The default statistical operation order should match the sequence above unless the user explicitly asks for a method change.

Clustering and embedding quality tuning is mandatory inside the fixed post-QC workflow, but clustering is a separate secondary task after Harmony. The clustering secondary task must read `h5ad/04-integration-harmony/adata_harmony.h5ad` as its source object. Do not rerun QC, normalize/log, HVG selection, regression, scaling, PCA, or Harmony inside each clustering candidate. If `adata_harmony.h5ad` is missing or invalid, stop the clustering step and fix or rerun `04-integration-harmony` first.

The tunable graph parameters are `n_neighbors` and `n_pcs`; Leiden `resolution` follows the default resolution-search rule. Do not accept a single blindly chosen clustering run as the final result. If the user explicitly specifies `n_neighbors`, `n_pcs`, and optionally `resolution`, run only the user-specified combination or resolution-search path and record it as user-specified. If the user does not specify graph/clustering parameters for Module 01 broad/all-cell clustering, run the fixed default grid: `n_pcs = 10, 15, 20, 25`, `n_neighbors = 10, 20, 30, 40`, and `resolution = 0.2, 0.3, 0.4, 0.5`. This default grid has 64 expected candidates. Do not only list candidate parameters in a plan; every candidate `n_neighbors`/`n_pcs` pair and every resolution in the default grid must be executed and must produce its own figure/table output directories and parameter record. The agent must not shrink the grid for speed, convenience, runtime, memory, or subjective judgment. If the full grid cannot run, stop the clustering task, document the blocker, and ask the user before reducing the grid. User-specified non-grid values must be used exactly and recorded as user-specified. After candidate generation, do not automatically choose a best candidate. Prepare the candidate outputs for manual user review; broad annotation starts only after the user explicitly chooses one completed candidate or has already provided the candidate parameters.

Before running the default Module 01 broad/all-cell clustering grid, write `tables/05-clustering-parameter-search/clustering_grid_manifest.csv` with one row per expected candidate. Required columns are `n_pcs`, `n_neighbors`, `resolution`, `candidate_label`, `expected_cluster_count_table`, `expected_parameter_table`, `expected_figure_dir`, `status`, `completed`, and `reason_if_skipped`. For the default grid, `expected_candidates` must equal 64. During execution, update candidate status as `planned`, `running`, `completed`, `failed`, or `skipped_user_approved`. Silent skipping is not allowed. Status values such as `skipped_for_speed`, `too_many_candidates`, or `not_needed` are not valid. Candidate manifest rows must not include an expected per-candidate h5ad path because Leiden parameter-search h5ad output is deferred until manual parameter selection.

For each `n_neighbors`/`n_pcs` graph candidate, start from a clean copy of `h5ad/04-integration-harmony/adata_harmony.h5ad`. A script may load the source object once as `adata_harmony = sc.read_h5ad(...)`. In the GPU/RAPIDS path, transfer this source template to GPU once with `rsc.get.anndata_to_GPU(adata_harmony)` before the candidate loop, then immediately before each `n_neighbors`/`n_pcs` graph candidate create a separate mutable object such as `adata_run = adata_harmony.copy()`. Run neighbors, UMAP, and the selected or fixed Leiden resolution values on `adata_run`; draw and save the candidate UMAP/diagnostic figures and lightweight tables; then delete `adata_run` and release GPU memory before the next graph candidate. In the CPU/Scanpy fallback path, use the same copy-before-mutation and delete-after-output pattern without GPU transfer. Do not run a second graph candidate on an AnnData object already mutated by a previous candidate's `.uns`, `.obsp`, `.obsm`, or `.obs` graph/UMAP/Leiden outputs. For Module 01 broad lineage or major cell-type clustering, use the fixed resolution values `0.2, 0.3, 0.4, 0.5` for every graph candidate unless the user provides another range. Do not carry over a resolution accepted for a previous `n_neighbors`/`n_pcs` pair; changing either value invalidates the previous resolution choice until the resolution search has been rerun for the new graph.

Save each executed clustering parameter-search candidate in matching lightweight directories named with the exact final values, for example `pcs-20_nn-30_res-0p5`: `figures/05-clustering-parameter-search/<pcs-nn-res>/`, `tables/05-clustering-parameter-search/<pcs-nn-res>/`, and `codes/05-clustering-parameter-search/<pcs-nn-res>/`. If a script evaluates several resolutions for the same `n_pcs`/`n_neighbors` graph and stores their figures/tables together, the shared directory must be named as an aggregate, for example `pcs-20_nn-30_res-all` or `pcs-20_nn-30_res-0p2-0p5`, not as a single resolution such as `pcs-20_nn-30_res-0p3`. A single-resolution directory must never contain outputs from other resolutions. Do not call `write_h5ad` for these candidates and do not create `h5ad/05-clustering-parameter-search/<pcs-nn-res>/`. Candidate parameter CSVs should record `source_h5ad = h5ad/04-integration-harmony/adata_harmony.h5ad`, `source_template_transferred_to_gpu_before_loop` when using rsc, `candidate_object_created_by_copy = True`, graph parameters, UMAP parameters, tested clustering resolution, clustering algorithm, raw cluster column name, cluster counts, output figure paths, cluster-count table path, code file, seed, backend/package, `candidate_h5ad_saved = False`, `candidate_object_deleted_after_outputs = True`, GPU memory release status when using rsc, and review notes. After candidate generation, write a candidate review manifest summarizing all completed candidates, cluster counts, and diagnostic figure paths. Do not write `tables/05-clustering-parameter-search/selected_clustering.csv` and do not save `h5ad/05-clustering-parameter-search/selected/adata_inte.h5ad` until the user manually chooses a completed candidate. Once the user chooses, rerun the selected graph/UMAP/clustering settings from a clean read of `h5ad/04-integration-harmony/adata_harmony.h5ad`; create a selected-run copy before mutation, save the required selected object at `h5ad/05-clustering-parameter-search/selected/adata_inte.h5ad`, and record that selected h5ad path plus the final raw clustering column in `selected_clustering.csv`. Broad annotation should read that selected h5ad path.

After all candidate jobs finish, write `tables/05-clustering-parameter-search/clustering_grid_completion_check.csv`. It must include `expected_candidates`, `completed_candidates`, `failed_candidates`, `skipped_candidates`, `missing_candidates`, `all_expected_lightweight_outputs_exist`, and `manual_selection_ready`. For the default grid, manual selection is ready only when `expected_candidates = 64`, `completed_candidates = 64`, `missing_candidates = 0`, and `all_expected_lightweight_outputs_exist = True`, unless a smaller grid was explicitly approved by the user and recorded in the manifest. `selected_clustering.csv`, when created after user choice, must select a candidate whose manifest row has `completed = True`.

Summarize cluster stability, batch/sample mixing, embedding diagnostics, and cluster-count behavior across candidate groups for manual review. Candidate-grid outputs are for parameter review and diagnostics, not the final object unless the user explicitly accepts a candidate run as final. The agent may flag obvious problems or suggest a shortlist, but must not declare a best candidate without user confirmation. Distinguish mandatory clustering/embedding quality tuning from resource-control tuning. Resource-control tuning should be justified by the observed resource blocker and safe cleanup/repair attempts, and should state whether the output remains the default full-data result or becomes a documented fallback.

Candidate graph/clustering jobs may be split into several independent scripts and run concurrently when candidate combinations are too many, provided memory and I/O limits are controlled. The preferred split granularity is one fixed `n_pcs` value per script, with that script covering the full `n_neighbors` and `resolution` grid for that `n_pcs`. For example, `codes/05-clustering-parameter-search/run_pcs10_nn10to40_res0p2to0p5.py` can handle `pcs=10`, `nn=10,20,30,40`, and `res=0.2,0.3,0.4,0.5`; another script can handle `pcs=20` with the same `nn/res` grid. The same `n_pcs` value may also be split across multiple scripts by non-overlapping `n_neighbors` subsets, for example one script for `pcs=10, nn=10,20` and another for `pcs=10, nn=30,40`. In every split script, each assigned `n_pcs`/`n_neighbors` pair must still run the full fixed resolution range `0.2,0.3,0.4,0.5`. Do not split in a way that drops a resolution or duplicates the same `pcs/nn/res` output. A split script may store all tested resolutions for one graph in a shared aggregate directory such as `pcs-10_nn-30_res-all`, but it must not store several resolutions under a single-resolution name such as `pcs-10_nn-30_res-0p3`. For the default grid, four fixed-`n_pcs` scripts or more granular fixed-`n_pcs`/`n_neighbors` subset scripts can cover all 64 candidates. GPU/rsc clustering scripts may also be launched concurrently, but only as bounded concurrency. Before launching script `i+1`, inspect current GPU processes and free VRAM with `nvidia-smi` or an equivalent API, account for the observed or estimated peak VRAM of each active clustering script plus a safety margin, and launch the next script only when it is unlikely to trigger OOM. If the per-script peak VRAM is unknown, start one GPU script first, monitor its peak usage, and use that measurement to set the concurrency limit before launching additional scripts. If there is not enough safe headroom for script `i+1`, queue it until an earlier script finishes and GPU memory is released. If concurrency causes OOM, reduce concurrency or queue the remaining scripts, release VRAM, record the OOM, and use CPU fallback for OOM candidates when needed; do not repair GPU solely for OOM and do not reduce the candidate grid. If concurrency causes stalls, CUDA context errors, or unstable half-written outputs for non-OOM reasons, reduce concurrency or queue the remaining scripts and follow the non-OOM GPU repair rule. Use explicit GPU assignment such as `CUDA_VISIBLE_DEVICES` or a scheduler when multiple GPUs are available, and record the GPU assignment, concurrency limit, free-VRAM check, safety margin, and launch order in the candidate parameter tables. Every script must read the same `h5ad/04-integration-harmony/adata_harmony.h5ad` and must write each `pcs/nn/res` combination to a unique single-resolution candidate label or to a clearly named aggregate resolution directory under `figures/05-clustering-parameter-search/`, `tables/05-clustering-parameter-search/`, and `codes/05-clustering-parameter-search/`, for example `pcs-10_nn-10_res-0p2`, `pcs-10_nn-40_res-0p5`, `pcs-20_nn-10_res-0p2`, or `pcs-25_nn-30_res-all`. Candidate scripts must not write per-candidate h5ad files. Never let two processes write the same figure directory, parameter table, cluster-count table, code output, `sc.settings.figdir`, or manual selected-output record. Record each process/script path and candidate parameters in its own parameter CSV, wait for all candidate jobs to finish, then prepare the completed outputs for manual user selection.

If Harmony does not converge or produces unstable correction, first increase `max_iter_harmony` as the same-method convergence adjustment, then rerun Harmony from the saved pre-Harmony/PCA input and overwrite only the intended Harmony-stage output for that rerun. Do not merely edit the parameter table. Record the original `max_iter_harmony`, the increased `max_iter_harmony`, convergence warnings/logs, elapsed time, output h5ad path, embedding key, and reason for accepting the rerun. Do not silently fall back to a failed or partially corrected embedding, and do not switch integration methods before the increased-`max_iter_harmony` Harmony rerun has been attempted or the user explicitly approves a method change.

Post-QC implementation scope:

```text
Use Scanpy/rapids-singlecell equivalents for normalize/log, HVG selection, raw preservation, HVG subsetting, regression, scaling, PCA, and Harmony integration in `04-integration-harmony`; use neighbor graph construction, UMAP, and Leiden clustering only in `05-clustering-parameter-search` after reading `h5ad/04-integration-harmony/adata_harmony.h5ad`.
Parameter values should use package defaults, common project-style values, or a user-approved plan, except for the fixed Module 01 broad/all-cell clustering grid documented above.
Record the active parameter range, embedding keys, candidate raw cluster columns, and manual selection status in the run report. Also save the validated `h5ad/04-integration-harmony/adata_harmony.h5ad` source path, per-candidate graph/UMAP/Leiden parameter CSVs under `tables/05-clustering-parameter-search/<pcs-nn-res>/`, explicit confirmation that no per-candidate h5ad was saved, and a selected-parameter file plus selected h5ad only after the user manually chooses a candidate.
The selected full-cell pre-Harmony parameters are the reference for Module 02 subtype reclustering. Record them explicitly enough that a subtype run can reuse the same normalization target sum, log transform, HVG parameters, raw-assignment point, HVG-subsetting rule, regression keys, scaling parameters, PCA settings, batch key, PCA basis, and Harmony input/output basis names.
```

Function/API range:

```text
RAPIDS/rapids-singlecell equivalents by default: normalize_total, log1p, highly_variable_genes, regress_out, scale, pca, harmony_integrate, neighbors, umap, leiden
AnnData operations: set raw after normalize/log/HVG marking, subset to HVGs, store obsm/obs outputs with clear integrated or unintegrated names
Legal placeholder call forms:
rsc.get.anndata_to_GPU(adata)
rsc.pp.normalize_total(adata, target_sum=target_sum_value)
rsc.pp.log1p(adata)
rsc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes_value)
adata.raw = adata
adata = adata[:, adata.var["highly_variable"]].copy()
rsc.pp.regress_out(adata, keys=regression_keys)
rsc.pp.scale(adata, max_value=scale_max_value)
rsc.tl.pca(adata, n_comps=pca_n_comps)
rsc.pp.harmony_integrate(adata, key=batch_key, basis=pca_basis, adjusted_basis=integrated_basis)
rsc.pp.neighbors(adata, n_neighbors=n_neighbors_value, n_pcs=neighbor_n_pcs, use_rep=integrated_basis)
rsc.tl.umap(adata)
rsc.tl.leiden(adata, resolution=leiden_resolution, key_added=cluster_key)
```

Optional pre-Harmony unintegrated diagnostic graph:

```text
Build an unintegrated neighbor graph, Leiden clustering, and UMAP using dataset-appropriate parameters when batch-effect diagnostics are needed.
Keep unintegrated output keys clearly labeled so they cannot be confused with Harmony-integrated outputs.
```

The unintegrated diagnostic graph above is optional and should not run by default. Use it only when the user asks to compare unintegrated versus integrated structure or when batch-effect diagnostics are needed.

Cell-count reduction is not allowed in Module 01. Do not downsample, subsample, randomly sample, sketch, or otherwise reduce cells for QC, integration, clustering, parameter tuning, visualization, or resource-control fallback. If the full dataset cannot run, document the blocker, close the process/kernel, release RAM/VRAM, and stop that step; do not produce a reduced-cell result.

## Optional scVI Integration

scVI integration is the documented backup integration route. It is used when the user explicitly requests scVI, when the user supplies a valid scVI-integrated object, or when the default Harmony workflow cannot complete and the user approves switching to the backup route. It is not the default integration route, and it should not silently replace the Harmony workflow above.

When scVI is requested, keep it as a separate secondary task named `04-integration-scvi`:

```text
1. Read `h5ad/03-qc/adata_qc.h5ad`; do not start from `adata_harmony.h5ad`, a scaled object, PCA output, UMAP output, or any already integrated object.
2. Use integer raw counts for scVI input. Prefer `adata.layers["counts"]`. If that layer is absent, use `.X` only after auditing that `.X` is raw count-like and recording the audit. If `.X` is normalized/log/scaled, stop and recover or regenerate a counts layer before training.
3. Preserve obs metadata: `sample`, `series`, `status`, `original_barcode`, and any dataset/source columns already present in the QC object.
4. Gene selection may subset genes but must not subset cells. If gene selection is used, default to the same HVG rule and `n_top_genes` used by Module 01 Harmony unless the user specifies otherwise. Record the before/after gene counts and keep the scVI count layer aligned to the selected genes.
5. Register AnnData with scvi-tools using `batch_key = "sample"` unless the user explicitly provides another batch column. Do not use `status`, `cell_type`, `leiden_coarse`, or downstream biological labels as batch covariates unless the user explicitly approves and the table records why.
6. Set fixed seeds before model setup/training: `SEED = 42`, `random.seed(SEED)`, `np.random.seed(SEED)`, `torch.manual_seed(SEED)`, and `scvi.settings.seed = SEED` when available.
7. Train `scvi.model.SCVI` with explicit or documented parameters. Default model parameters are `n_latent = 30`, `n_layers = 2`, and `gene_likelihood = "nb"` unless the user specifies otherwise. Use `max_epochs = None` to allow the scvi-tools default scheduler unless the user gives a fixed epoch count. Record whether `early_stopping` is supported/enabled.
8. Store the learned latent representation in `obsm["X_scVI"]`.
9. Save the scVI-integrated object as `h5ad/04-integration-scvi/adata_scvi.h5ad` and save the trained model under `h5ad/04-integration-scvi/scvi_model/` or record the external model path in `tables/04-integration-scvi/scvi_model_path.txt`.
10. Save scVI setup, training, latent/model, count audit, package-version, backend/GPU, elapsed-time, and seed records under `tables/04-integration-scvi/`.
11. Do not overwrite `h5ad/04-integration-harmony/adata_harmony.h5ad`. Do not call `X_scVI` Harmony, and do not use scVI merely because it is installed.
```

Canonical scVI code structure:

```python
import random
import numpy as np
import scanpy as sc
import torch
import scvi

SEED = 42
random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)
if hasattr(scvi, "settings"):
    scvi.settings.seed = SEED

adata = sc.read_h5ad("h5ad/03-qc/adata_qc.h5ad")
batch_key = "sample"
counts_layer = "counts"

if counts_layer not in adata.layers:
    # Use .X only after an input-count audit proves it is raw count-like.
    # Otherwise stop and regenerate a counts layer from the raw-count source.
    if not x_is_raw_count_like_after_audit:
        raise ValueError("scVI requires raw counts; recover/regenerate adata.layers['counts'] first.")
    adata.layers[counts_layer] = adata.X.copy()

sc.pp.highly_variable_genes(
    adata,
    layer=counts_layer,
    n_top_genes=n_top_genes_value,
    flavor=hvg_flavor_value,
    batch_key=batch_key if hvg_batch_key_enabled else None,
)
adata = adata[:, adata.var["highly_variable"]].copy()

scvi.model.SCVI.setup_anndata(
    adata,
    layer=counts_layer,
    batch_key=batch_key,
)
model = scvi.model.SCVI(
    adata,
    n_latent=30,
    n_layers=2,
    gene_likelihood="nb",
)
model.train(max_epochs=None)
adata.obsm["X_scVI"] = model.get_latent_representation()
model.save("h5ad/04-integration-scvi/scvi_model", overwrite=True)
adata.write_h5ad("h5ad/04-integration-scvi/adata_scvi.h5ad")
```

If scVI is selected as the integration source for clustering, the clustering task must read `h5ad/04-integration-scvi/adata_scvi.h5ad` instead of `h5ad/04-integration-harmony/adata_harmony.h5ad`, and neighbors must use `use_rep = "X_scVI"`. For scVI clustering candidates, vary `n_neighbors` and Leiden/Louvain `resolution`; do not treat Harmony `n_pcs` as PCA PCs. If a latent-dimension truncation is explicitly tested, record it as `latent_n_dims_used`, not as Harmony `n_pcs`, and keep the source representation `X_scVI` unchanged. Raw scVI clustering columns must start with `scvi_leiden_` or `scvi_louvain_`, for example `scvi_leiden_res0p5`. The selected final object may still be saved as `h5ad/05-clustering-parameter-search/selected/adata_inte.h5ad`, but `tables/05-clustering-parameter-search/selected_clustering.csv` must record `integration_source = scVI`, `source_h5ad = h5ad/04-integration-scvi/adata_scvi.h5ad`, and `use_rep = X_scVI`.

Use scvi-tools GPU acceleration through PyTorch/CUDA when available. If the scVI GPU environment is broken for a non-OOM reason, first try to repair it without changing the integration method. If the scVI GPU step OOMs, release VRAM, record the OOM, and continue with the equivalent scVI CPU run. If non-OOM GPU repair fails, document the reason and continue with scVI on CPU. Cell-count reduction is never an allowed scVI fallback.

## Long-Run Fallback

Run the default integration method first. The default method is the fixed Harmony workflow after QC. Use the default method on full data unless the user explicitly requests scVI or another documented route.

Before switching methods because the workflow cannot run, stop and document the failed step, resource state, and safe cleanup attempts. Do not change HVG count, PCA components, Harmony settings, neighbor parameters, or integration method to reduce memory or runtime without explicit user approval. Do not change cell count to reduce memory or runtime; cell-count reduction by downsampling, subsampling, sketching, or random sampling is not allowed. When another parameter is changed to reduce GPU memory, CPU memory, or runtime, label it as resource-control tuning rather than clustering-quality tuning.

If the fixed default integration step still cannot run after safe cleanup/repair attempts, stop that step cleanly, close the process/kernel, release RAM/VRAM, and ask before using any backup route. If the user approves a backup route, use this preference order:

```text
1. Full-data scVI integration.
```

This fallback option applies only after the default method has been attempted and the user approves using a backup route. scVI is the backup integration method, but it is not the default. Do not create reduced-cell fallback outputs.

Record the elapsed time before fallback, the failed/stalled step, whether GPU OOM or GPU repair failure occurred, the backup method selected, the user approval, and the reason for selecting it. Do not choose any backup method outside the documented options above unless the user explicitly approves it.

Global clustering QC:

```text
Apply this clustering QC standard to Module 01 discrete clustering steps, including unintegrated diagnostics, Harmony-based atlas clustering, and broad-lineage clustering. Module 02 should reuse this QC standard for lineage/subtype reclustering, but those subtype outputs must stay in the Module 02 analysis directory.
Use the "tight, clear, stable" rule: clusters should not be stringy, fragmented, or smeared.
Tune and validate clustering rather than accepting a resolution blindly.
Every final discrete clustering must come from either the user-specified parameter run or a documented candidate grid followed by manual user selection. If the user explicitly specifies `n_neighbors`, `n_pcs`, and optionally `resolution`, run that requested setting only through the active approved integration source and record it as user-specified. If the user does not specify these parameters for Module 01 broad/all-cell clustering, use the fixed default grid: `n_pcs = 10, 15, 20, 25`, `n_neighbors = 10, 20, 30, 40`, and `resolution = 0.2, 0.3, 0.4, 0.5`. The default grid must have 64 expected candidates, must be declared in `tables/05-clustering-parameter-search/clustering_grid_manifest.csv`, and must be verified in `tables/05-clustering-parameter-search/clustering_grid_completion_check.csv`. The agent must not run only a subset because the grid seems large. If the full grid cannot run, stop, document the failed or stalled step, and ask the user before reducing the grid. A grid that changes only clustering resolution while keeping `n_neighbors` and `n_pcs` fixed does not satisfy the default clustering QC requirement. Repeating the exact same parameter combination with the same fixed seed is not required because it should be deterministic; repeat identical parameters only when the user explicitly asks for seed sensitivity or reproducibility checks. After all candidates are generated, do not auto-select the best-supported parameters. Present or save the candidate review outputs, then wait for the user's manual choice before producing any final selected-clustering record, final cluster column handoff, or broad-annotation input. Once the user chooses parameters, rerun that selected setting from the clean approved integration object and save or redraw the selected result without deleting candidate-grid outputs. If old candidate-grid outputs are absent afterward and only the selected/final result remains, do not silently assume this is normal; inspect the selected-parameter record and final h5ad/provenance, record that the grid outputs are missing, and ask before rerunning only if the final selected result or its parameter record is missing or inconsistent. By default, run graph construction, UMAP, and clustering from `h5ad/04-integration-harmony/adata_harmony.h5ad`. If the user-approved scVI branch is the active source, run from `h5ad/04-integration-scvi/adata_scvi.h5ad` with `use_rep = "X_scVI"` and record `integration_source = scVI`. Do not rerun from `adata_qc` inside the clustering step. Never rerun QC for this parameter-search/final-clustering step unless the user explicitly changes QC thresholds or asks to redo QC.
For broad lineage or major cell-type clustering, use only the fixed default candidate range above unless the user provides another range. For lineage-specific or subtype reclustering in Module 02, the final resolution can be higher than the broad-lineage resolution if the extra clusters have stable UMAP/embedding structure; marker support should then be checked after clustering before assigning biological subtype labels.
Do not make the tested resolution range unnecessarily high. Stop increasing resolution once additional clusters mostly split existing stable groups, create tiny unstable clusters, or fragment continuous/bridge regions without clear structural support.
Do not default to overly high clustering resolution, overly high neighbor counts, or overly many PCs. For Module 01 broad/all-cell clustering, stay within the fixed default grid unless the user specifies another grid. When neighbor/PC parameters are varied, treat each `n_neighbors` and `n_pcs` pair as a separate graph candidate. In the rsc path for the default Harmony route, read `h5ad/04-integration-harmony/adata_harmony.h5ad`, transfer that Harmony template to GPU before the loop, copy it once per `n_neighbors`/`n_pcs` graph candidate, run the relevant clustering resolution values, save figures and lightweight tables, then delete the candidate object and release GPU memory. In the approved scVI route, read `h5ad/04-integration-scvi/adata_scvi.h5ad`, use `use_rep = "X_scVI"`, and record scVI latent-dimension handling separately from Harmony `n_pcs`. Use Leiden by default. For broad/all-cell clustering, treat a raw Leiden result with more than 100 clusters as abnormal unless the user explicitly pre-approves such a high-cluster analysis. Also treat clearly non-monotonic cluster-count behavior when lowering resolution or an output that fails the tight/clear/stable QC check in a way inconsistent with the graph and embedding as abnormal. When Leiden is abnormal, document the abnormal Leiden run and switch that candidate to Louvain rather than accepting the abnormal Leiden result. When switching to Louvain, every cluster column, candidate table, DEG directory, DEG filename, and downstream reference to the raw cluster label must use the Louvain name, such as `louvain_res0p3`, not `leiden_res0p3`. Save each candidate's outputs under a directory named with the algorithm when needed to avoid ambiguity, for example `pcs-20_nn-30_louvain_res-0p5`; do not overwrite or merge outputs across graph, resolution, or algorithm settings. Do not continue from another candidate object that already contains a different graph, UMAP, Leiden, Louvain, or integration source. Each such directory must include candidate parameter CSVs documenting the source h5ad, integration source, graph/clustering values, clustering algorithm, UMAP parameters, GPU-transfer behavior, copy-before-mutation behavior, delete-after-output behavior, and resolution-search path. Encode decimal resolutions with `p` instead of `.`, for example `res-0p5`.
Clusters should not be so many that obvious stable structures are fragmented during the clustering-parameter selection step.
Clusters should not be so few that major lineages or clearly distinct states are merged.
Clusters should not be dominated by individual samples, series, or processing batches unless this reflects a known biological or clinical difference.
For Harmony-based atlas clustering, default clustering-parameter comparison UMAPs should be colored by series, status, and the raw cluster column generated by the selected clustering algorithm; broad marker gene plots can be generated after the user manually selects a candidate for annotation validation. Use sample-colored UMAPs only as diagnostic outputs when checking batch/sample mixing or when the user explicitly asks for them.
For non-Harmony or lineage-specific clustering-parameter comparison, check the relevant embedding colored by series/status when available and the candidate cluster label; lineage/subtype marker genes can be inspected after candidate clusters are generated for annotation validation. Add sample-colored views only as diagnostics when sample effects are a concern.
Use structural embedding quality, cluster-count behavior, and sample/series mixing diagnostics to summarize whether each candidate looks acceptable for manual clustering-parameter review. Inspect marker expression after candidate clusters are generated for annotation and biological validation, not as the first-pass grid-review criterion.
Document all tested resolutions, tested neighbor/PC settings when varied, UMAP/embedding parameters, clustering algorithm, confirmation that every candidate started from `h5ad/04-integration-harmony/adata_harmony.h5ad`, confirmation that QC and Harmony were not rerun inside clustering candidates, confirmation that no per-candidate h5ad files were saved, candidate figure/table paths, candidate raw cluster column names, and cluster counts. If the user manually selects a candidate, append the selected resolution or equivalent clustering parameter, selected h5ad path, final source object, final cluster column name, and user-provided or documented selection note. Keep this clustering-quality record separate from resource-control changes such as reducing feature count, PCA dimensions, candidate grid size, or batch size to avoid memory/runtime failure. Do not reduce cell count for resource control.
Do not manually merge Leiden clusters as a hidden correction. If clusters are too fragmented, rerun or tune the graph/clustering parameters and document the final choice. It is acceptable for multiple Leiden clusters to be annotated as the same broad cell type or major lineage when marker evidence supports that annotation; keep the original Leiden cluster labels available.
If a cluster is stretched into a long string, thin filament, or bridge on UMAP, do not force it into a discrete subtype interpretation. Treat it as a possible continuous state, resolution/neighbor instability, or batch/sample gradient; inspect marker genes and consider lower resolution, adjusted neighbor parameters, diffusion map, or pseudotime analysis.
```

## Major Annotation

Major annotation belongs with the first broad/all-cell clustering step in this module. After the user manually selects a completed broad clustering candidate, annotate those broad clusters immediately using broad-lineage DEGs and canonical markers, and save the resulting broad labels before starting Module 02 lineage-specific subtype reclustering. Do not postpone the initial broad annotation until after subtype analysis, and do not mix broad-cluster annotation with the later subset-specific subtype clustering workflow.

The broad annotation workflow must preserve the active DEG export block and marker dotplot validation step. The broad score/rank QC workflow must preserve `score_genes(use_raw=True)`, exclude PCA variance-inspection-only code, and keep only cells whose `best_rank_type_global` maps to the same label as `leiden_coarse`.

Use this naming contract for the broad Leiden annotation pass:

```text
leiden_res<RES> = raw Leiden cluster label generated at resolution <RES>
louvain_res<RES> = raw Louvain cluster label generated at resolution <RES> when Louvain is used after abnormal Leiden fallback
leiden_coarse = initial broad cell_type label annotated from leiden_res<RES>
```

For broad/all-cell clustering, `leiden_coarse` is strictly a post-annotation
column name. It must not be used as the raw clustering column or as a synonym for
the selected Leiden/Louvain result. The raw clustering column must remain
`leiden_res<RES>` or `louvain_res<RES>` until DEG review and broad annotation
merge raw clusters into final broad labels.

During pure clustering and clustering-parameter search, do not assign biological names, lineage prefixes, or marker-gene labels. Candidate clustering outputs should keep only raw numeric clustering labels under `leiden_res<RES>` or `louvain_res<RES>` plus candidate parameter labels such as `pcs-20_nn-30_res-0p5`. Names in the form `<lineage_abbrev>_<marker_gene>` are allowed only after DEG review and annotation have started. If the workflow has not reached annotation yet, do not rename clusters into `leiden_coarse`, `cell_type`, or `cell_subtype`.

Encode decimal resolutions with `p`, for example `leiden_res0p5` or `louvain_res0p5` for `resolution=0.5`. Do not use `leiden_coarse` as the raw clustering column. Keep the selected raw cluster column, such as `leiden_res0p5` or `louvain_res0p5`, alongside the annotated `leiden_coarse` column.

If `leiden_res<RES>` contains more raw clusters than the final annotated biological categories, never copy `leiden_res<RES>` directly into `cell_type`, `cell_subtype`, or another biological annotation column. Raw Leiden clusters are technical grouping labels for DEG calculation and annotation review. Multiple `leiden_res<RES>` clusters may map to the same annotated `leiden_coarse` category when marker evidence supports the same initial broad `cell_type`.

Treat `leiden_coarse` as the merged broad label produced after annotating raw Leiden clusters. The raw clusters that were merged into one `leiden_coarse` label are not subtypes of that label. For example, if raw global clusters 10 and 11 both map to `leiden_coarse = 'Epithelial Cells'`, clusters 10 and 11 must not be called epithelial subtypes; they remain raw broad-clustering groups used only to create the broad annotation. Epithelial subtypes must be produced later by Module 02 after extracting cells with `leiden_coarse == 'Epithelial Cells'` and reclustering that subset.

For the first broad-lineage clustering/annotation pass, use full lineage names by default. Avoid terse labels such as `T`, `Mye`, or `Endo` as final broad-lineage names. Accepted exceptions are widely used names whose expanded form is awkward or less standard in figures, such as `pDC`. Abbreviated prefixes may still be used inside subtype labels and mapping tables.

The broad-lineage names below are examples, not a fixed answer key for every project. For a new organ, disease, species, platform, or atlas composition, determine the broad `leiden_coarse` categories from the actual top DEGs, canonical markers, organ biology, and user-provided annotation goals. Do not force missing lineages such as `pDC`, `Mast Cells`, or `B Cells` into a dataset when the marker evidence is absent. Add organ-specific broad labels when supported, and record the marker/DEG evidence for any added or omitted category.

Example broad lineages:

```text
Epithelial Cells
T Cells
NK Cells
Myeloid Cells
B Cells
Endothelial Cells
Stromal Cells
Mast Cells
pDC
```

After broad annotation exists, order `leiden_coarse` categories with epithelial first for plots, summaries, tables, and downstream handoff when the data contain epithelial cells. A recommended order is `Epithelial Cells`, `T Cells`, `NK Cells`, `Myeloid Cells`, `B Cells`, `Endothelial Cells`, `Stromal Cells`, `Mast Cells`, `pDC`, followed by any dataset-specific labels. This ordering rule does not apply before annotation, when only raw numeric Leiden clusters exist.

Use the broad-cluster annotation to create `leiden_coarse` before subtype reclustering. Treat `leiden_coarse` as the initial broad `cell_type`; if a downstream workflow requires an explicit `cell_type` column, initialize it from `leiden_coarse`. After subtype naming is finalized in Module 02 and projected back to the full atlas, regenerate or validate the final `cell_type` from the `cell_subtype` prefix. The subtype prefix can override the initial `leiden_coarse` broad label when they conflict. For example:

```text
leiden_coarse = NK Cells
cell_subtype = T_GZMK
final cell_type = T Cells
```

Recommended prefix-to-cell-type mapping for final validation:

```text
Epi -> Epithelial Cells
T -> T Cells
NK -> NK Cells
Mye -> Myeloid Cells
B -> B Cells
Endo -> Endothelial Cells
S -> Stromal Cells
Mast -> Mast Cells
pDC -> pDC
```

Broad-cluster DEG and marker-score QC rule:

```python
sc.tl.rank_genes_groups(adata, groupby=col, method='t-test', use_raw=True)
```

All DEG calculations must use `adata.raw` normalized/log expression explicitly
with `use_raw=True`. If `adata.raw` is absent at a DEG step, stop and fix the
upstream raw-normalized expression preservation before exporting DEGs.

Broad annotation must run two DEG rounds. Round 1 is before annotation, using the
manually selected raw cluster column. Round 2 is after annotation, after raw
clusters with the same broad identity have been merged into `leiden_coarse`.
Neither round is optional unless the user explicitly supplies trusted matching
DEG CSVs for that exact `groupby` column.

Round 1: use the manually selected raw Leiden cluster column first. For example,
if the user-selected candidate uses `resolution=0.5`, set `col =
'leiden_res0p5'`, compute raw-cluster DEGs, and save outputs as:

```text
tables/degs_leiden_res0p5_pcs<PCS>_nn<NN>_res<RES>/
{group}_degs_leiden_res0p5_pcs<PCS>_nn<NN>_res<RES>.csv
```

Every per-group DEG CSV must include the project-style canonical DEG columns
used by the reference analysis:

```text
gene, score, logfoldchanges, pvals, pvals_adj
```

If the chosen DEG method emits additional documented statistics, preserve them
after these canonical columns rather than dropping them. Do not replace these
canonical columns with unrelated names unless the user explicitly requests a
different DEG engine and column schema.

During annotation, DEG export is mandatory and must save one full-length CSV per
raw cluster or annotated group. Do not save top-only DEG files such as top50,
top100, top200, top300, or any `topXX` DEG CSV. Top-N subsets may be used only
in memory for marker review, naming, dotplots, or score-gene selection after the
full per-group DEG CSVs have been written.

Build each per-group CSV from `adata.uns['rank_genes_groups']`, for example:

```python
pd.DataFrame({
    'gene': result['names'][group],
    'score': result['scores'][group],
    'logfoldchanges': result['logfoldchanges'][group],
    'pvals': result['pvals'][group],
    'pvals_adj': result['pvals_adj'][group],
})
```

Use the raw-cluster DEGs together with canonical broad-lineage markers to annotate each `leiden_res<RES>` cluster into `leiden_coarse`. Do not annotate from cluster number alone. By default, review the first 50 rows from each already-saved full DEG CSV. If those rows do not provide enough interpretable marker evidence for the broad annotation level, expand the whole broad annotation level sequentially to the first 100, 200, and 300 rows from the same full CSVs. All raw Leiden clusters at this annotation level must use the same selected review depth. Record the selected review depth and why expansion was needed, but do not save separate top-only DEG CSVs or top-N DEG manifests.

Round 2: after `leiden_coarse` is assigned and raw clusters with the same broad
identity have been merged under the same `leiden_coarse` label, recompute DEGs
using `col = 'leiden_coarse'` and save outputs as:

```text
tables/degs_leiden_coarse_pcs<PCS>_nn<NN>_res<RES>/
{group}_degs_leiden_coarse_pcs<PCS>_nn<NN>_res<RES>.csv
```

For the broad annotation marker validation plot, use the annotated broad label column and raw-normalized expression. Use only the top 3 representative marker genes per broad cell category by default. The order of marker blocks in `var_names` must match the plotted `leiden_coarse` category order, and genes within each block should remain in the chosen marker order. For example, if the category order is `Epithelial Cells`, `T Cells`, `NK Cells`, then the dotplot marker dictionary should list epithelial markers first, T-cell markers second, and NK-cell markers third. Do not show long marker panels in the default broad dotplot unless the user asks for extended validation.

```python
sc.pl.dotplot(
    adata,
    var_names=markers,
    groupby="leiden_coarse",
    standard_scale="var",
    save="_all_leiden_coarse.pdf",
    use_raw=True,
)
```

If `adata.raw` is absent, stop and document why the raw-normalized expression layer was not preserved before drawing the broad `leiden_coarse` dotplot.

Use the `leiden_coarse` DEG tables to build broad-lineage marker-score evidence. This step must read the saved full DEG CSV files produced for `col = "leiden_coarse"`; do not build score/rank gene sets directly from preset marker lists, hard-coded lineage vocabularies, in-memory adata objects, or expected labels. Score/rank gene sets must come only from DEG tables for broad labels that were actually identified in the current data and are present in `adata.obs["leiden_coarse"]`. Do not compute score/rank columns from a preset or expected lineage list. For example, if the current annotation has no `B Cells`, there is no `B Cells` DEG table and the workflow must not calculate a B-cell `score_genes` column. For each observed `leiden_coarse` group, select the first 100 usable genes from that group's full `leiden_coarse` DEG table, restricted to genes present in the AnnData object. Use those genes with `sc.tl.score_genes` to create one score column per observed `leiden_coarse` label only. If labels contain spaces or punctuation, use stable sanitized score-column names and save the mapping between score columns and original `leiden_coarse` labels.

Rank the score columns per cell so the best-supported broad label can be recovered. Store both the score columns and a best-label column, for example:

```text
<leiden_coarse_label>_score
<leiden_coarse_label>_score_rank or <leiden_coarse_label>_score_rank_pct
best_rank_type_global
```

The selected/best rank label must correspond to the highest broad-lineage marker-score evidence for that cell. Compute rank percentiles for each broad-lineage score with larger score ranked better, for example `rank(ascending=False, pct=True)`. Therefore a smaller rank percentile means stronger expression/score evidence: 1 percent is better than 2 percent, and 2 percent is better than 3 percent. For each cell, set `best_rank_type_global` to the broad-lineage score column with the smallest rank percentile. For the broad post-annotation QC pass, keep a cell only when this best-ranked broad-lineage label maps exactly to the cell's current `leiden_coarse` identity. Example: if a cell is annotated as `B Cells` and the B-cell top-100 DEG score has rank percentile 1 percent while all other lineage scores have worse/larger rank percentiles such as 2 percent or 3 percent, keep it. If that same `B Cells` cell has a T-cell score rank percentile smaller than the B-cell score rank percentile, remove it from the default filtered object. Save both the unfiltered scored object for audit and the filtered object for downstream handoff. The score/rank consistency filter is the required default project QC branch whenever score/rank evidence is generated; downstream modules should use the filtered object unless the user explicitly cancels this filtering. Also save per-cell and summary tables separating consistent cells from inconsistent cells, plus kept/removed counts per `leiden_coarse` group. Do not silently overwrite the unfiltered atlas.

Expected major-annotation outputs:

```text
obs['leiden_res<RES>']
obs['leiden_coarse'] as the initial broad cell_type
obs['cell_type'] initialized from obs['leiden_coarse'] when required
obs['best_rank_type_global']
all broad *_score and *_score_rank or *_score_rank_pct columns
unfiltered broad-annotated h5ad
filtered broad-annotated h5ad after best-rank agreement QC
per-cell table marking best-rank/leiden_coarse consistent and inconsistent cells
separate consistent-cell and inconsistent-cell tables
major lineage UMAP
canonical marker dotplot
cell type composition tables/plots
tables/degs_<groupby>_pcs<PCS>_nn<NN>_res<RES>/
one full DEG CSV per raw cluster or annotated group; no topXX DEG CSVs
```

## Outputs

Expected handoff objects:

```text
h5ad/02-merge-metadata/adata_merge.h5ad
h5ad/03-qc/adata_qc.h5ad
h5ad/04-integration-harmony/adata_harmony.h5ad
h5ad/05-clustering-parameter-search/selected/adata_inte.h5ad required after user manual selection and clean selected rerun
tables/05-clustering-parameter-search/selected_clustering.csv only after user manual selection
h5ad/06-broad-annotation/adata_anno.h5ad
h5ad/07-score-rank-qc/adata_anno_score_genes_rank.h5ad
h5ad/07-score-rank-qc/adata_anno_score_genes_rank_consistent.h5ad
obs['leiden_coarse'] as the initial broad cell_type
obs['cell_type'] initialized from obs['leiden_coarse'] when required
unintegrated diagnostic UMAPs when requested or needed
integrated series/status UMAPs
sample-colored UMAPs only as diagnostic outputs, not default final figures
```

## Validation

- `obs_names` are unique.
- `original_barcode` is preserved.
- Raw counts are preserved.
- `sample`, `series`, and `status` exist.
- Mitochondrial-content metrics are present or explicitly marked as unavailable.
- Doublet scores/calls are present when doublet detection is requested, or the reason for skipping doublet detection is documented.
- Integration does not erase sample identity needed downstream.
- Batch-effect diagnostics include before/after views when possible.
- Broad `cell_type` annotation is saved with marker/DEG evidence before Module 02 subtype analysis starts.

## 02-Cell Subtype Integration Clustering

This block is the former Module 02, copied as the canonical instruction source
for the compact workflow. Inside this big skill, write this block's outputs under
`epi-cm-core-workflow/{codes,h5ad,tables,figures}/02-cell_subtype_integration_clustering/`
unless the user explicitly asks to use the original numbered module output tree.
Do not use the shorter summary as a substitute for these copied rules.

# 02-project-cell-annotation

Use this skill after Module 01 has produced broad all-cell clustering, `leiden_coarse` as the initial broad `cell_type`, and broad marker-score filtering outputs when available. Major all-cell annotation belongs in Module 01; this module handles lineage-specific subtype annotation, refinement, and projection back to the full atlas. The subtype clustering and annotation procedure is the same for every selected `leiden_coarse` broad label. If a example dataset exists, treat it as an example only; downstream branches must be chosen from the project biology and user goal rather than assumed.

For lineage-specific subclustering, use `adata_qc` as the expression/reclustering base. Use the Module 01 annotation/scored h5ad only to choose cell IDs and transfer broad labels by cell ID. By default, if `adata_anno_score_genes_rank_consistent.h5ad` exists, select lineage IDs from that file using the requested `leiden_coarse` label, then subset `adata_qc` by those IDs. If that consistent object does not exist and the user explicitly cancels or does not require score/rank filtering, use the final integrated annotated object such as `adata_inte.h5ad` only as the annotation/ID source, then still subset `adata_qc` for expression and reclustering. The score/rank filter is the required default project QC branch when rank evidence exists. The user can explicitly cancel it with `--rank-filter-mode never`, and that cancellation must be recorded.

## Non-Negotiable Subtype Annotation Gate

Module 02 subtype annotation is not complete when clustering or raw-cluster DEG
export finishes. Module 02 subtype annotation is complete only after all of
these artifacts exist for the selected lineage:

```text
selected clustered lineage h5ad from 04-subtype-selected-clustering
raw-cluster DEG CSVs for the selected clustering key
candidate subtype annotation table for user review or agent-authorized review
final raw-cluster-to-subtype mapping CSV
obs['cell_subtype'] written from that mapping
obs['functional_state'] or an equivalent state column written from that mapping
post-annotation DEG CSVs grouped by cell_subtype
final subtype UMAP colored by cell_subtype
final annotated lineage h5ad named h5ad/05-subtype-deg-annotation/<lineage>/adata_<cell_abbrev>.h5ad
```

Do not call a Module 02 subtype task finished if the latest h5ad still contains
only raw labels such as `leiden_res<RES>` or `louvain_res<RES>` and no biological
`cell_subtype` column. Do not copy a raw Leiden/Louvain cluster column into
`cell_subtype`. If the final mapping CSV is absent, incomplete, duplicated, or
not user-confirmed or explicitly agent-authorized, stop after writing the
candidate annotation table and ask for confirmation before saving the final
annotated h5ad.

The final mapping CSV must be one row per selected raw cluster and must cover
every selected raw cluster exactly once. It must contain at least:

```text
cluster
cell_subtype
functional_state
gene_selection_rationale
```

`cell_subtype` values must be unique within the current lineage. If two or more
raw clusters have the same top candidate subtype gene, assign that top-gene
label to the cluster with stronger evidence for that gene, then assign the other
cluster or clusters the next-ranked eligible positive DEG from their own DEG
tables. Evidence strength should be compared using adjusted P value when
available, then raw P value, then Scanpy score or absolute positive
logfoldchange. Record the conflict, the winner cluster, the skipped duplicate
gene, and the next-ranked gene chosen in `gene_selection_rationale`. Do not merge
raw clusters into one `cell_subtype` label during annotation.

If marker review shows that a raw subcluster resembles another lineage or
lineage subtype, do not move those cells into another lineage-specific h5ad, do
not edit the other lineage's h5ad, and do not change membership of any other
lineage object because a subcluster looks like that lineage. Annotate the
subcluster within the current lineage-specific object using the current lineage
prefix and the standard most-significant-positive-DEG naming rule, and document
the marker-supported alternate lineage in `functional_state` or
`gene_selection_rationale`. Move or restore cells to another lineage only if the
user explicitly asks for that operation.

## Pre-Execution Plan Requirement

Before executing code from this skill, write a concise method-and-result plan that the user can review and copy as the goal. Keep it result-oriented rather than overly procedural. Include only:

```text
analysis goal / expected result
method route to use
main inputs or provided intermediates
major code modules to run or skip
expected output figures/tables
key validation criterion
```

Do not start long-running analysis, dependency installation, or file-rewriting steps until this short plan has been stated. For simple inspection-only tasks, one or two sentences are enough.

If the user does not provide a manual choice for parameters, thresholds, method options, output naming, or optional branches, use the documented default settings in this skill and state that the default was used.

## Project Organization and Figure Output Contract

Treat each numbered module folder as its own output boundary with one shared four-directory layout: `figures/`, `tables/`, `codes/`, and `h5ad/`. Module output directory names must use the active project slug, for example `02-<project_slug>-cell-annotation/`; for BRCA use `02-brca-cell-annotation/`. The four top-level category directories may be created at module setup. Secondary-module/task directories inside those category directories must be created only when that category will receive at least one real output for that task. Do not pre-create empty task directories under `figures/`, `tables/`, `codes/`, or `h5ad/` just to mirror the layout. A directory creation command for a secondary task or candidate must be coupled to writing a real output file there; if the output is not generated, do not leave that task directory behind. Secondary modules are logical analysis units inside that numbered module, but they do not require a directory under every category. This is a write-location rule, not a read restriction: a module may read/reuse files and already generated outputs from other modules as inputs, but newly generated outputs for the current module must be written inside the current numbered module. Across the whole skill workflow, an agent must not delete, clear, overwrite, or move any existing output file or directory anywhere unless the user explicitly names the exact path and operation. During normal module execution, do not write, move, overwrite, clear, or delete files or directories inside any other numbered module output directory. Also do not delete, clear, overwrite, or move any existing output file or directory inside the current module unless the user explicitly names the exact path and operation. If a new result would conflict with an existing output, write to a new versioned path or stop and ask. Deleting any output directory is never part of a module run; it requires a separate explicit cleanup request naming the exact path. For example, Module 02 subtype outputs go under `02-<project_slug>-cell-annotation/{figures,tables,codes,h5ad}/03-subtype-clustering-grid/<lineage>/`, `04-subtype-selected-clustering/<lineage>/`, `05-subtype-deg-annotation/<lineage>/`, or `06-project-subtypes-to-full-adata/<lineage>/`, not under `01-<project_slug>-singlecell-integration/`, not directly as files under `02-<project_slug>-cell-annotation/`, and not under the top-level skill source directory.

Use stable numbered secondary-module/task names that describe the analysis step, lineage, method, or figure group. Reuse the same secondary-module/task name only under category directories that actually receive outputs from that analysis, so files stay aligned without creating empty placeholder task directories. For example:

Module 02 secondary tasks:

```text
01-lineage-selection
  choose eligible cells from Module 01 broad labels and optional score/rank consistency evidence
02-subtype-harmony
  subset adata_qc by eligible IDs, rerun normalize/log/HVG/raw/regress/scale/PCA/Harmony, then save the lineage Harmony h5ad
03-subtype-clustering-grid
  read the saved lineage Harmony h5ad, run the subtype clustering candidate grid, and save lightweight figures/tables/code only
04-subtype-selected-clustering
  read the saved lineage Harmony h5ad, rerun the user-selected graph/UMAP/clustering setting with GPU/rsc when available, then save the selected UMAP and selected clustered lineage h5ad
05-subtype-deg-annotation
  read the selected clustered lineage h5ad, calculate raw-cluster DEGs, create or consume the confirmed mapping CSV, write cell_subtype and functional_state, recompute final subtype DEGs, draw the final cell_subtype UMAP, and save the annotated lineage h5ad
06-project-subtypes-to-full-adata
  project annotated lineage/subtype labels back to the full atlas h5ad by stable cell IDs, save the projected full h5ad, draw required projection QC UMAPs, and report unmatched/duplicated cells
```

```text
02-<project_slug>-cell-annotation/
  codes/
    01-lineage-selection/
      epithelial/
        select_epithelial_ids.py
    02-subtype-harmony/
      epithelial/
        run_epithelial_harmony.py
    03-subtype-clustering-grid/
      epithelial/
        run_epi_subclustering_grid.py
    04-subtype-selected-clustering/
      epithelial/
        run_epi_selected_clustering.py
    05-subtype-deg-annotation/
      epithelial/
        run_epi_selected_subtype_annotation.py
    06-project-subtypes-to-full-adata/
      epithelial/
        project_epi_subtypes_to_full_adata.py
  h5ad/
    02-subtype-harmony/
      epithelial/
        adata_epithelial_harmony.h5ad
    04-subtype-selected-clustering/
      epithelial/
        adata_epi_selected_clustered.h5ad
    05-subtype-deg-annotation/
      epithelial/
        adata_epi.h5ad
    06-project-subtypes-to-full-adata/
      adata_anno_cellsubtype.h5ad
  figures/
    03-subtype-clustering-grid/
      epithelial/
        pcs-20_nn-30_res-all/
          umap_leiden_res0p4.pdf
    04-subtype-selected-clustering/
      epithelial/
        umap_leiden_res0p4_selected.pdf
    05-subtype-deg-annotation/
      epithelial/
        umap_cell_subtype.pdf
    06-project-subtypes-to-full-adata/
      umap_leiden_coarse_vs_projected_cell_type.pdf
      umap_projected_cell_subtype.pdf
      epithelial/
        umap_epithelial_cell_subtype_with_projected_subtype_palette.pdf
  tables/
    01-lineage-selection/
      epithelial/
        selected_cell_ids.csv
    02-subtype-harmony/
      epithelial/
        harmony_parameters.csv
    03-subtype-clustering-grid/
      epithelial/
        pcs-20_nn-30_res-all/
          cluster_counts.csv
    04-subtype-selected-clustering/
      epithelial/
        selected_clustering_parameters.csv
    05-subtype-deg-annotation/
      epithelial/
        epithelial_subtype_counts.csv
    06-project-subtypes-to-full-adata/
      epithelial/
        projection_match_report.csv
```

By default, save executable/reproducibility code under the current module's shared `codes/<secondary-module>/`, using ordered names such as `01_read_merge.ipynb`, `02_qc.py`, or `03_integrate.R`. Save AnnData-like objects under `h5ad/<secondary-module>/` as `.h5ad`, `.loom`, `.rds`, or equivalent files with stable names. Save corresponding figure files under `figures/<secondary-module>/` and use ordered names such as `01_umap.pdf` or `02_marker_dotplot.svg`. Save text-like and tabular outputs under `tables/<secondary-module>/`, such as CSV/TSV/XLSX/TXT/JSON/YAML logs, manifests, reports, mapping files, and parameter records. `figures/` should contain figure files only. `tables/` should contain text-like and tabular outputs only. `codes/` should contain executable/reproducibility code only. `h5ad/` should contain AnnData-like/intermediate object files only. Add `tables/<secondary-module>/readme.txt` documenting the input files, including any cross-module input/output files that were read, code order, h5ad/loom/rds objects, output figures/tables, and any skipped optional branches.

Do not write new h5ad, code, figures, or tables directly into the numbered module root. The numbered module root may contain the module `SKILL.md`, lightweight module-level index files, or manually curated high-level notes, but executable outputs should live under the shared four category directories. If a simple task has only one natural step, still use a small secondary-module/task name such as `01-main` or a task-specific name from the Module 02 secondary-task list, but create that task directory only under the category directories that receive real outputs.

If one analysis step outputs multiple files or figures, put that output set in the same named secondary-module subdirectory under `figures/`, `tables/`, `codes/`, or `h5ad/`, using the same analysis prefix when possible.

If an output already exists, do not rerun only to recreate it in the new layout. Do not move or delete existing outputs for layout cleanup unless the user explicitly names the exact path and operation. Prefer to leave existing outputs in place, copy them into the organized location only when provenance is recorded, then update the corresponding code paths so future runs write to the same organized location.

When a task creates a run, lineage, candidate parameter set, or method variant, create matching candidate subdirectories inside the active secondary-module/task directory only under category parents that receive outputs for that candidate. Lineage folders must sit under the secondary-module/task directory, not directly under `codes/`, `figures/`, `tables/`, or `h5ad/`. For example, use `codes/03-subtype-clustering-grid/epithelial/` and `figures/03-subtype-clustering-grid/epithelial/pcs-20_nn-30_res-all/`, not `codes/epithelial/` or `figures/epithelial/`. For multi-candidate or multi-condition runs, use matching candidate names under each relevant parent when needed; for example, create `figures/03-subtype-clustering-grid/epithelial/pcs-30_nn-15_res-0p8/` only if figures will be saved and create `tables/03-subtype-clustering-grid/epithelial/pcs-30_nn-15_res-0p8/` only if parameter or cluster-count tables will be saved. A directory name containing a single resolution, such as `pcs-25_nn-30_res-0p3`, must contain only outputs for that exact resolution. If one directory intentionally stores multiple resolutions for the same graph, do not name it after one resolution; use an explicit aggregate name such as `pcs-25_nn-30_res-all` or an explicit range such as `pcs-25_nn-30_res-0p2-0p4`. Inside such aggregate directories, every output filename and every table row must still include the exact algorithm and resolution, such as `louvain_res0p2`, `louvain_res0p3`, or `leiden_res0p4`. For Leiden clustering parameter searches in this module, candidate outputs are intentionally lightweight: do not create h5ad candidate directories and do not save per-candidate AnnData objects. Keep candidate code files under `codes/<secondary-module>/<lineage>/` with parameter-coded subdirectories when code is emitted. After the user selects final subtype-clustering parameters, subtask 04-subtype-selected-clustering must rerun the selected graph/UMAP/clustering settings from the saved lineage Harmony h5ad and save only the selected raw-cluster UMAP plus the selected clustered lineage h5ad under `04-subtype-selected-clustering/<lineage>/`. Subtask 05-subtype-deg-annotation then reads that selected clustered h5ad for DEG calculation and subtype annotation. Keep the old subtype grid figures/tables/code for other parameter candidates in place after the selected rerun is saved. Rerunning or redrawing the selected/final subtype result is not a cleanup request and must not delete grid-search outputs. Delete grid-search outputs only if the user gives a separate explicit cleanup request naming the exact path or output group to remove. Do not delete the selected clustered lineage h5ad, selected-parameter record, annotated lineage h5ad, subtype annotation tables, candidate manifests, or non-grid outputs during any cleanup.

Each analysis that produces an output should have corresponding source code under the current module's `codes/`. Acceptable code artifacts include `.ipynb`, `.py`, `.R`, and `.sh`, depending on the language actually used. Do not leave a figure, table, or exported result that can only be traced to manual GUI editing. If an analysis uses Python, keep the notebook and/or `.py` script that generates it; if it uses R, keep the `.R` script or R notebook; if both languages are used, keep both code artifacts under `codes/` with clear ordered prefixes. When converting notebooks to upload/download versions, keep the executable cells needed to reproduce the outputs and remove stale display output only when requested.

Each executed run should also create or update a parameter/provenance report under `tables/`, such as `tables/run_parameters.txt`, `tables/run_parameters.csv`, or a step-specific report in the same output subdirectory. The report should list the code file used, input files/objects, output files, exact parameters, random seeds, selected candidate/final settings, skipped steps, fallback decisions, and any user-approved method changes.

Do not substitute another analysis method, algorithm, statistical test, visualization strategy, database, or input layer without explicit user permission. If the specified method cannot run, stop that module, document the blocker in `tables/readme.txt`, and ask for confirmation before using any alternative. Any approved or documented method change should state why the original method was unsuitable or failed and why the replacement method is appropriate for the same analysis goal.

Parameter-integrity rule: if this skill specifies an exact parameter, threshold,
grid, iteration count, feature count, random seed, model, or output filename, do
not reduce, simplify, approximate, or replace it for speed, estimated runtime,
convenience, memory concerns, or an agent's judgment that a smaller run is
"enough". Runtime estimates are not authorization to change parameters. If the
specified run appears long, first verify the estimate from real local evidence
when possible, start or continue the specified run when resources are available,
or ask the user before changing anything. Any exploratory short run must be
explicitly requested or approved by the user, labeled exploratory, and must not
be used as the canonical result.

Missing-data rule: if a required input file, matrix, h5ad object, metadata column,
gene list, label column, or upstream result is missing, empty, unreadable, or
inconsistent, stop and ask the user for the correct input or permission to rerun
the upstream step. Do not randomly generate data, simulate fake matrices, create
placeholder labels, fabricate metadata, fill missing biological values with
invented numbers, or use toy/example data to complete a real analysis. This rule
applies to every module task, including subtype annotation, plotting,
correlation analysis, and validation. A tiny synthetic object may be used only
for isolated software smoke tests, and it must not be saved as an analysis
output or used to support any biological result.

When a task, notebook run, script run, or long interactive kernel finishes, promptly close the process/kernel/session and release CPU memory and GPU memory. Do not leave idle Python, R, Jupyter, CellChat, RAPIDS, PyTorch, TensorFlow, or CUDA processes holding RAM/VRAM after the requested work is complete.

After each module finishes, create or update `tables/package_versions.txt` describing the packages and tools used by that module. Include Python packages, R packages, command-line tools, CUDA/GPU libraries when relevant, interpreter/R version, environment name or path, and the code files that used them.

Install missing dependencies when they are required to execute the specified method or its approved acceleration path. This includes installing a compatible GPU-accelerated implementation when the method supports it and the machine has a usable GPU/CUDA driver, for example installing `rapids-singlecell`/RAPIDS to run Scanpy-style preprocessing through `rsc`. Dependency installation is allowed to make the requested method work; method substitution is not allowed without explicit user permission. For packages or methods that already provide GPU acceleration, enable and use the GPU-accelerated path after installing any missing compatible GPU packages and verifying imports/minimal execution. If the requested package/method has no GPU-accelerated implementation, use its normal CPU path. If an expected GPU path is installed but broken or incompatible for a non-OOM reason, including a CUDA/CUDA-tag mismatch such as `cu11` vs `cu12` wheels, CuPy/RAPIDS/PyTorch wheels incompatible with the visible driver, missing CUDA runtime libraries, `libucx`/UCX errors, or `cuCtxGetDevice`/CUDA context errors, first try to repair or reinstall a compatible GPU environment without changing the requested method. Choose a compatible wheel, channel, or uv environment automatically from `nvidia-smi`, Python version, platform, and package compatibility information; do not ask the user to choose the CUDA tag. Ask the user only before system-driver changes, OS package changes that require elevated privileges, deleting an existing environment, or replacing a working environment used by other analyses. If GPU runs out of memory, inspect active GPU processes, close stale or idle processes left by previous tasks/kernels when they can be safely identified, release VRAM, record the OOM, and switch that OOM step directly to the equivalent CPU implementation; do not repair, reinstall, or repeatedly retry GPU solely for OOM. Do not terminate unrelated active user processes unless the user explicitly approves. If compatible dependency/GPU installation or non-OOM repair fails, or no usable GPU is present, document the reason and continue with the normal CPU path for the same requested method.

This GPU backend rule applies to all GPU-capable code in every module. If the first full execution of the required task completes successfully with the planned backends, do not rerun only to validate the backend plan. If a GPU-accelerated step fails from GPU OOM, release VRAM, record the OOM, and treat CPU fallback for that step as pre-approved. If a GPU-accelerated step fails for a non-OOM reason and the user approves CPU fallback after documented repair/cleanup attempts, treat fallback at step granularity: use the CPU equivalent only for the failed step when needed, then continue later steps with GPU whenever those later steps have a valid GPU path. Do not mark the entire remaining workflow CPU-only because one GPU step failed. Record a backend capability table under the relevant `tables/<secondary-module>/` or `tables/<secondary-task>/` directory with one row per executed step, including `step`, `planned_backend`, `attempted_backend`, `status`, `error_summary`, `fallback_backend`, `clean_input_reloaded`, and `final_backend_for_rerun`. Also export `gpu_backend_capability_summary.csv` and, when useful, `gpu_backend_capability_summary.txt` in that same tables directory; these files must state which steps can use GPU, which steps must use CPU, and the reason for each CPU step. When any GPU failure, fallback, or partial object mutation occurs during this exploratory/profiling pass, finish the required task only to learn which steps can use GPU, then start a fresh Python/R process, reload the nearest clean upstream input h5ad/RDS/input files from disk, and rerun the whole required task once using the recorded final backend plan. In that final rerun, every step marked GPU-capable must use GPU, except steps with recorded GPU OOM or approved non-OOM CPU fallback, which should use the recorded CPU fallback so the final run avoids repeating known mid-run GPU failures while still using GPU wherever it works. Do not reuse in-memory AnnData/R objects, arrays, GPU buffers, fitted models, graphs, or partial metadata from the profiling pass. Do not present outputs from a partial-failure/profiling pass as canonical final outputs.

Python package management rule: use `uv` by default for Python dependencies, and save uv-managed environments under a dedicated subdirectory in the total analysis directory for the dataset/project. Create the layout `uv_envs/<category>/.venv` under the top-level analysis root, where `<category>` is a stable dependency category such as `main`, `rapids`, `velocity_cellrank`, `cellchat_liana`, or `survival`. Use `uv_envs/main/.venv` for the default shared Python stack, and create another category only when dependency compatibility requires it. Install with `uv pip install --python uv_envs/<category>/.venv/bin/python --index-url https://pypi.tuna.tsinghua.edu.cn/simple ...` or activate that category `.venv` before `uv pip install --index-url https://pypi.tuna.tsinghua.edu.cn/simple ...`. Do not create a separate per-module environment unless it is also a documented dependency category, do not create or reuse a global environment, do not put the run environment inside the skill source directory, and do not run bare `pip install` into the system/user Python unless `uv` is unavailable or the user explicitly requests it. Keep any category-level `pyproject.toml`, `uv.lock`, or requirements export inside `uv_envs/<category>/`, and record all environment paths, categories, package versions, and mirror fallbacks under the total analysis directory, usually `tables/package_versions.txt`. Prefer TUNA mirrors for package downloads: use the TUNA PyPI mirror for `uv pip install`, and use TUNA CRAN/Bioconductor mirrors for R packages when practical. If the TUNA mirror is unavailable, stale, or missing a required package, fall back to the official source only for the affected dependency and record the mirror fallback in `tables/package_versions.txt`. If `uv` is unavailable, install `uv` first when possible; otherwise document the fallback package manager in `tables/package_versions.txt`. R packages, Cell Ranger, velocyto, CUDA drivers, and system libraries are outside `uv` and should be installed with their appropriate manager while still recording versions.

Do not create a symlink in the current directory that points to an older output just to make it look renamed. If an output needs a new name, either regenerate it in the correct output directory or make a real copied file with documented provenance.

Python random seed rule: every Python script and notebook must define and use the default fixed random seed near the top of the file: `SEED = 42`, unless the user explicitly specifies another seed. Set `random.seed(SEED)` and `numpy.random.seed(SEED)` when those libraries are used, and set framework-specific seeds for stochastic packages when relevant, such as PyTorch, TensorFlow, scvi-tools, scVelo, veloVI, CellRank, scikit-learn, Scanpy, or UMAP. Pass `random_state=SEED`, `seed=SEED`, or the package-equivalent argument to every function that supports it, including PCA, neighbors/UMAP, Leiden/clustering, train/test splits, model fitting, bootstrapping, permutations, and plotting layouts when applicable. If a function has no seed argument or still has nondeterministic GPU behavior, document that limitation. Record the fixed seed in the module parameter log, per-candidate `harmony_params.txt` or equivalent, and `tables/package_versions.txt`; do not change the seed between candidate groups unless the user explicitly requests seed sensitivity testing.

Expression-layer rule: for expression-based plotting, marker visualization, dotplots, violin plots, gene scoring, or other Scanpy-compatible calculations, if the function exposes a `use_raw` argument and `adata.raw` is present, set `use_raw=True` by default unless the user explicitly requests a different layer or the method requires counts. Record any exception and the expression layer actually used.

Cell-count reduction rule: do not downsample or subsample cells for analysis, plotting, parameter tuning, or resource-control fallback unless the user explicitly requests or approves that specific reduced-cell run. If the full data cannot run, document the blocker and stop to ask instead of producing a reduced-cell result automatically. Any approved reduced-cell output must be labeled as reduced-cell/exploratory and must not be presented as the full-data result.

Default figure formats are PDF and/or SVG only. Do not create, save, convert, or request PNG files for final, intermediate, diagnostic, preview, thumbnail, or temporary figure outputs. If a tool defaults to PNG, override it to PDF/SVG or stop and ask; do not leave `.png` files under `figures/`, `tables/`, `codes/`, or `h5ad/`.

Global figure style for every module:
- Use the module-specific `Module Figure Style Contract` in this SKILL.md as the plotting-style source of truth. Do not rely on any material outside this SKILL.md to infer final figure style. Keep exactly one current canonical plotting route for each final figure, and label any other route as non-canonical or exploratory; do not duplicate final PDFs/SVGs for the same panel or file stem.
- Apply the style for the corresponding analysis type and module-specific section first; use these global rules only as baseline guardrails. Do not copy a figure style from an unrelated analysis type just because the file format or package is similar.
- Save figure outputs as PDF and/or SVG only; never create PNG previews, thumbnails, diagnostics, or temporary figure files.
- Keep text editable whenever the backend supports it. In Matplotlib set `pdf.fonttype = 42`, `ps.fonttype = 42`, `svg.fonttype = "none"`, and use a standard sans-serif font such as DejaVu Sans or Arial.
- Keep axis ticks visible on quantitative and categorical plots, including heatmaps, dotplots, barplots, forest plots, scatter plots, survival plots, and UMAP panels with axes. Do not call `axis("off")`, remove tick labels, or hide spines unless the plot type is a pure network/chord/graph layout where axes have no coordinate meaning.
- Use clean white backgrounds, black axis text, readable tick labels, and legends/colorbars outside or to the right when practical. Rotate dense x-axis labels rather than letting them overlap.
- Prefer the official plotting interface for the relevant package before manual low-level drawing. Use manual Matplotlib/ggplot2/Seaborn layout code only when the package interface cannot express the required final figure or when this skill gives explicit canonical manual code; record that reason in the code or parameter log.
- Size every figure element for the exported canvas: labels, legends, colorbars, risk tables, titles, p/q labels, arrows, node labels, and panel titles must be readable and must not collide. If any element overlaps or is clipped in the generated PDF/SVG, increase canvas size, margins, row/column spacing, legend placement, font size, label wrapping, or panel spacing and rerender before considering the figure complete.
- For UMAP-like embedding panels, keep panels square, use the official Scanpy/scVelo plotting interface when possible, keep framed axes, and use consistent palettes for the same labels across full-atlas and lineage-specific plots.
- For heatmaps, use a colorbar with visible ticks. Use centered diverging palettes for signed correlations/effects, sequential palettes for nonnegative scores or transition weights, and method-specific labels in titles/captions.
- For dotplots, keep gene and group orders explicit, preserve tick labels, and use standard scaling only when the method requires display scaling.
- For forest/survival plots, keep confidence intervals, hazard-ratio/reference lines, p/q labels, and axis ticks visible.
- For network/chord/directed-graph plots, use PDF/SVG, editable labels, clear legends/colorbars when edge weights are encoded, and document any intentional axis removal.


Global figure style rule: every figure produced by this skill must use a
`ticks` visual style unless the user explicitly requests another style or the
plotting package does not expose a comparable style option. In Python/Seaborn
code, call `sns.set_style("ticks")` or `sns.set_theme(style="ticks")` before
plotting. For Matplotlib-only and Scanpy figures, keep visible tick-style axes
when axes are shown and do not switch to `white`, `whitegrid`, or hidden-axis
styling as the default. For R/ggplot2 figures, use a ticked axis theme such as
`theme_classic()` unless the package-specific official plot does not support it.
Record any exception in the run-parameter table or a short code comment.

Scanpy plotting hard rule: For ordinary `sc.pl.*` outputs that support `save=...`, do not create Matplotlib axes, do not pass `ax=...`, and do not save with `fig.savefig` or `ax.figure.savefig`. Use `sc.settings.figdir` plus the Scanpy `save` argument. For multi-panel Scanpy plots, use the package interface such as `color=[...]`, `ncols`, `wspace`, or `standard_scale` instead of manual subplots. Manual `ax=...` is allowed only for documented special overlays or incompatible per-panel settings that Scanpy cannot express; record the reason in the run-parameter table or a code comment.


For both Python and R plotting, use the official plotting interface of the relevant package by default, such as Scanpy, Matplotlib, Seaborn, CellRank, scVelo, ggplot2, ComplexHeatmap, or CellChat plotting APIs. For Python-generated single-panel plots, use a default square canvas of 2.5 x 2.5 inches unless the user specifies another size or the plot type clearly requires more space, such as multi-panel layouts, heatmaps, wide dotplots, survival plots, or network/chord diagrams. For UMAP plots, use the default Scanpy save path and keep the call minimal. Set figure parameters with `sc.set_figure_params(figsize=(3, 3), dpi=150)` or `sc.settings.set_figure_params(figsize=(3, 3), dpi=150)`, set `sc.settings.figdir` to the target output directory, then call `sc.pl.umap(adata, color="leiden_coarse", save="_name.pdf")`. Keep the default Scanpy-style framed axes and outside legends; do not manually create `plt.subplots`, pass `ax=...`, or call `fig.savefig` for ordinary UMAPs. Use `ncols` only when `color` contains multiple objects, for example `sc.pl.umap(adata, color=["cell_subtype", "status", "cnv_score"], ncols=3, wspace=0.4, save="_celltype_status_cnvscore.pdf")`. Do not use `ncols` for a single-color UMAP. Keep `save` as a suffix/name handled by Scanpy rather than a full path. Manual axes and `fig.savefig` are reserved for special cases where different panels require incompatible per-panel parameters or post-processing that the official `color=[...]` interface cannot express, and the reason must be documented because manual saving can greatly increase PDF/SVG size. For ordinary UMAPs, specifically avoid `return_fig=True` followed by `ax.savefig(..., bbox_inches="tight")`; on large atlases this can inflate files from sub-MB Scanpy-saved outputs to multi-MB PDFs or very large SVGs. When the task is only to inspect, audit, or explain existing UMAP code/output size, do not modify the source code or rerun plotting unless the user explicitly asks for a fix or rerender. Keep all figure text editable as text whenever the plotting backend supports it; do not convert labels, legends, tick labels, titles, or annotations to outlines/paths unless the user explicitly requests it. For multi-panel figures, especially multi-panel UMAPs, verify that each UMAP panel remains square, legends/colorbars/titles do not overlap, and adjacent panels do not collide. It is acceptable to adjust the default figure width, height, `wspace`, `hspace`, legend font size, or margins to prevent overlap and preserve square UMAP panels. Do not add extra custom titles to UMAP panels; use the Scanpy default title derived from `color` unless the user explicitly asks for custom titles. Do not add automatic bitmap/raster conversion rules; let the user decide figure-size tradeoffs from the actual output files. Do not draw sample-colored UMAPs as default final figures; generate sample-colored UMAPs only when the user asks for them or when they are needed as integration/batch-mixing diagnostics, and label them as diagnostic outputs.

There are two QC layers:

```text
pre-annotation QC = remove low-quality cells/doublets before integration
post-annotation QC = required default marker-evidence filtering after broad labels are assigned when Module 01 score/rank evidence exists; user may explicitly cancel it
```

Post-annotation marker-score evidence is the required default project QC branch when Module 01 has generated score/rank evidence. Score/rank calculation saves a filtered companion object for downstream handoff, but it must not silently replace the unfiltered atlas, which remains an audit object. For lineage-specific detailed subclustering, use the filtered object or apply the same consistency rule in-script before selecting IDs. The user may explicitly cancel this filtering; if cancelled, document the cancellation and use coarse-label-only selection intentionally. Always record whether IDs came from a pre-filtered h5ad such as `adata_anno_score_genes_rank_consistent.h5ad` or were filtered in-script.

## Module Figure Style Contract

Use the following subtype-annotation figure styles unless the user explicitly
asks for a different style. Do not mention prior-run implementation provenance in generated reusable code, figure labels, captions, or readme files.

- Subtype clustering grid UMAPs: for each lineage and each candidate
  `pcs`/`nn` graph, put all tested resolutions for that graph in one
  clearly named PDF/SVG when practical. Use square framed Scanpy UMAP panels,
  stable palettes, and titles that include the exact algorithm/key such as
  `leiden_res0p5` or `louvain_res0p5`.
- Selected subtype clustering UMAPs: after the user chooses parameters, rerun
  from the saved lineage Harmony h5ad, plot the selected raw cluster key, and
  save the selected clustered h5ad. Do not delete grid-search figures unless
  the user separately asks for cleanup.
- Subtype DEG/annotation dotplots: use `use_raw=True` when available,
  `standard_scale='var'`, top marker genes in the reviewed order, and explicit
  subtype order. Keep gene and group tick labels visible.
- Projected full-atlas subtype UMAPs: plot `leiden_coarse` versus projected
  `cell_type` as one comparison figure, plot `cell_subtype` separately, and use
  the full-atlas palette consistently when drawing lineage-specific subtype
  h5ad objects.
## Inputs

```text
integrated h5ad with counts/log expression
sample/series/status metadata
Module 01 broad labels: obs['leiden_coarse'] as the initial broad cell_type, plus obs['cell_type'] if initialized from it
Module 01 marker-score/rank filtering columns when available
lineage-specific marker genes
lineage-specific subclustering results or marker tables
```

## Subtype Annotation

Subtype annotation is a second-stage workflow after major annotation. Apply the same subtype-reclustering workflow to every `leiden_coarse` broad label that needs subtype resolution; epithelial cells are only one example, not a special-case workflow. Extract the relevant broad lineage or confirmed cell population, rerun the same clustering workflow on that subset, calculate DEGs, and assign `cell_subtype` labels. Different `leiden_coarse` groups may use different marker sets, subtype names, and final resolutions, but the extraction, preprocessing, Harmony, candidate clustering, manual selection, DEG, annotation, and output-layout rules are the same. Perform or reuse lineage-specific subclustering and marker annotation. These lineage-specific outputs belong in the current Module 02 analysis directory, not in the Module 01 directory. Use Module 01 objects only as input/provenance sources for broad labels, score/rank evidence, and `adata_qc` IDs. Keep lineage-specific biological add-ons separate from this shared subtype step: only epithelial subtypes have the project-default extra epithelial-state branch for inferCNV projection, DPT/diffusion map/Monocle3, stemness scoring, epithelial pathway/marker analyses, and epithelial fate/velocity analyses; do not run those epithelial-only branches for every `leiden_coarse` group unless the user explicitly asks for an analogous analysis.

```text
obs['leiden_res<RES>'] or obs['louvain_res<RES>']
obs['cell_subtype']
pre-annotation cluster DEG tables
post-annotation subtype DEG tables
lineage-specific marker dotplots
subtype UMAPs
```

Any clustering performed in this module should follow the global clustering QC standard from Module 01: default clustering QC expects clusters to be tight, clear, and stable; not stringy, fragmented, or smeared; and not dominated by sample/series/batch without a documented biological reason. For lineage-specific or subtype-level reclustering, keep the workflow order consistent with the full-cell Module 01 pipeline through Harmony, neighbor graph construction, UMAP, clustering, DEG calculation, and annotation. The only extra first step is extracting target cells by the confirmed broad annotation, then returning to `adata_qc` for the expression matrix. Pre-Harmony operation order and numeric parameters should match the selected full-cell Module 01 settings by default, including `normalize_total` target, log transform, HVG method and count, raw-assignment point, HVG-subsetting rule, regression keys, scaling parameters, PCA settings, Harmony batch key, and Harmony input/output basis names. Subset-specific HVGs are still recomputed inside the extracted lineage, but the HVG rule and count should match Module 01 unless the user approves a change.

Subtype clustering should usually be more detailed than broad all-cell clustering. Do not accept a single coarse subtype run unless the user provides exact parameters. Generate a grid of candidate runs for the current lineage/subtype subset and let the user choose before final subtype annotation. Each broad cell lineage is an independent subtype-clustering task. Do not force all lineages to use the same final `n_pcs`, `n_neighbors`, `resolution`, UMAP parameters, or clustering algorithm. The set of lineage tasks must be derived from the actually observed and validated `leiden_coarse` categories in the current project, not from a hard-coded list. Examples such as T cells, epithelial cells, myeloid cells, endothelial cells, and B cells are examples only; run only the lineages present in the project annotation unless the user explicitly requests an additional subset.

For each subtype lineage, the default Leiden grid is fixed unless the user explicitly changes it: `n_pcs = 10, 15, 20, 25, 30, 35, 40, 45, 50`, `n_neighbors = 10, 15, 20, 25, 30, 35, 40, 45, 50`, and `resolution = 0.1, 0.2, ..., 1.5`. This means 81 graph/UMAP candidates and 1215 Leiden resolution outputs per lineage. Each `n_pcs` and `n_neighbors` pair is one graph candidate: build neighbors and UMAP once for that pair, then run every resolution from 0.1 through 1.5 on the same copied graph object. Every resolution must be represented in a UMAP output. For the same `n_pcs` and `n_neighbors` pair, all resolution-colored UMAP panels must be saved into one review figure file, preferably one multi-panel or multi-page PDF named with an aggregate token such as `pcs-25_nn-30_res-0p1-1p5`. Do not emit 15 separate UMAP PDFs for the same `n_pcs` and `n_neighbors` pair unless the user explicitly asks for split files. The companion cluster-count/parameter table for that graph must contain one row per resolution and include the exact raw cluster key, such as `leiden_res0p1` through `leiden_res1p5`.

Each observed lineage may still have different final selected parameters because its cell count, heterogeneity, and marker-supported structure differ. During this subtype grid-search/candidate phase, output only lightweight figures, tables, and code/provenance; do not save candidate h5ad files and do not create candidate directories under `h5ad/`. Save a selected clustered lineage h5ad only after the user selects a completed grid candidate or explicitly authorizes agent-led selection, then rerun that selected setting from the clean saved lineage Harmony h5ad in `04-subtype-selected-clustering`. This selected-clustering step must save the selected raw-cluster UMAP, selected clustering parameter record, and selected clustered h5ad before DEG or subtype annotation begins. Rerunning, redrawing, or saving the selected lineage/subtype result must preserve the grid-search figures/tables/code by default. If a future inspection sees only one selected/final lineage result and no visible grid candidate outputs, record that the grid outputs are missing; do not silently assume they were intentionally deleted. Check the selected-parameter record, selected clustered h5ad, annotated lineage h5ad, and provenance first, and request a rerun only if the selected result or required provenance is missing or inconsistent. Keep the default principle: enough resolution to separate stable marker-supported states, but not so much that it fragments continuous states, creates stringy/smeared groups, or splits clusters without DEG support. Use Leiden by default. If a raw Leiden result has more than 100 clusters, treat it as abnormal unless the user explicitly pre-approves such a high-cluster analysis; document the abnormal Leiden run and switch that candidate to Louvain. When switching to Louvain, every cluster column, candidate table, DEG directory, DEG filename, and downstream reference to the raw cluster label must use the Louvain name, such as `louvain_res0p3`, not `leiden_res0p3`. Document all candidate parameters, clustering algorithm, and the user-selected candidate when one is chosen. For subtype annotation, do not merge multiple raw subtype clusters into one `cell_subtype`. The final mapping must be one raw selected cluster to one final subtype label. If two or more selected raw clusters would receive the same subtype label only because their top positive DEG is identical, keep that top-gene label for the cluster with stronger evidence and assign the other cluster or clusters the next-ranked eligible positive DEG from their own DEG tables. If several clusters are genuinely indistinguishable even after reviewing ordered positive DEGs and marker evidence, flag this in the candidate table and ask the user whether to choose a less fragmented candidate or accept manual naming.

During subtype clustering and clustering-parameter search, do not assign biological subtype names yet. Candidate outputs should keep only raw numeric Leiden labels under `leiden_res<RES>` or Louvain labels under `louvain_res<RES>` plus candidate parameter labels. Prefix-plus-marker names such as `<lineage_prefix>_<GENESYMBOL>` are assigned only after DEG review and subtype annotation. If the workflow has not reached the annotation step, do not rename candidate clusters into `cell_subtype`.

For lineage-specific or subtype-level clustering parameter searches, write and run one independent script or script instance per broad cell lineage. A single subtype-clustering script run must process exactly one selected `leiden_coarse` lineage, such as one epithelial, myeloid, T-cell, endothelial, or B-cell subset; it must not loop over several lineages and write mixed outputs into one shared path. Several lineage-specific scripts may be launched at the same time when machine resources allow it, so multiple cell types can be processed concurrently, but each concurrent job must read its own eligible source cells and write only to its own lineage/candidate subdirectories under `figures/03-subtype-clustering-grid/<lineage>/`, `tables/03-subtype-clustering-grid/<lineage>/`, and `codes/03-subtype-clustering-grid/<lineage>/`; include the exact `pcs`, `nn`, algorithm, and resolution range in those paths. For the default subtype grid, each candidate directory should represent one `n_pcs` and `n_neighbors` graph, such as `pcs-25_nn-30_res-0p1-1p5`, and should contain one combined UMAP review figure covering all resolutions from 0.1 to 1.5 plus tables with one row per resolution. Do not store several resolutions under a single-resolution directory such as `pcs-25_nn-30_res-0p3`, and do not scatter same-graph resolutions across 15 separate figure files. Do not call `write_h5ad` for candidate clustering runs and do not create `h5ad/03-subtype-clustering-grid/<lineage>/<candidate>/` during parameter search. Subtype clustering is a separate secondary task after subtype Harmony, exactly like the full-cell Module 01 flow. The candidate-clustering task must read the saved subtype Harmony h5ad, such as `h5ad/<secondary-module>/<lineage>/adata_<lineage>_harmony.h5ad`, as its source object. Do not run neighbors, UMAP, Leiden, or Louvain on a just-created in-memory Harmony object before saving and reloading it for the clustering task. If the saved subtype Harmony h5ad is missing or invalid, stop the clustering step and fix or rerun the subtype Harmony step first. A clustering script may load the saved source object once as `adata_harmony = sc.read_h5ad(...)`. In the GPU/RAPIDS path, transfer this subtype Harmony template to GPU once before the candidate loop with `rsc.get.anndata_to_GPU(adata_harmony)`. Immediately before every `n_neighbors`/`n_pcs` graph candidate, create a separate mutable object such as `adata_run = adata_harmony.copy()`. Run neighbors and UMAP once on `adata_run`, then run Leiden for every default resolution from 0.1 to 1.5 on that same `adata_run`; save one combined UMAP review figure for all resolution columns, save lightweight tables, then delete `adata_run` and release GPU memory before the next graph candidate. In the CPU/Scanpy fallback path, use the same copy-before-mutation and delete-after-output pattern without GPU transfer. Do not run a second graph candidate on an AnnData object already mutated by a previous candidate's `.uns`, `.obsp`, `.obsm`, or `.obs` graph/UMAP/Leiden/Louvain outputs. Do not let parallel jobs share a figure directory, cluster-count table, run-parameter file, code output, or `sc.settings.figdir`. After all jobs finish, summarize candidate figures/tables for manual user review; do not write or accept the final `cell_subtype` annotation until the user chooses a completed subtype-clustering candidate or explicitly authorizes agent-led selection. After selection, subtask 04-subtype-selected-clustering must read the saved subtype Harmony h5ad, copy it before mutation, rerun the selected graph/UMAP/clustering parameters, save the selected raw-cluster UMAP, and save exactly one selected clustered lineage h5ad under `h5ad/04-subtype-selected-clustering/<lineage>/adata_<cell_abbrev>_selected_clustered.h5ad`. Subtask 04-subtype-selected-clustering must follow the same backend priority as the grid search: use GPU/RAPIDS/rsc for transfer, neighbors, UMAP, Leiden, and Louvain when available and recorded as GPU-capable; use CPU fallback only for documented GPU OOM or approved non-OOM GPU failure according to the global backend rule. Record the selected-run backend for each step in the selected-clustering parameter table. Subtask 05-subtype-deg-annotation must then read that selected clustered h5ad for raw-cluster DEG calculation, subtype annotation, post-annotation DEG calculation, subtype UMAP plotting, and annotated h5ad saving.

For each cell lineage or subtype-reclustering task, keep code files under `codes/<secondary-module>/<lineage>/`. Intermediate h5ad filenames should describe their state, for example `h5ad/02-subtype-harmony/epithelial/adata_epithelial_harmony.h5ad` and `h5ad/04-subtype-selected-clustering/epithelial/adata_epi_selected_clustered.h5ad`. Only the annotation-completed final subtype h5ad uses the short cell-abbreviation name exactly: `h5ad/05-subtype-deg-annotation/<lineage>/adata_<cell_abbrev>.h5ad`, where `<cell_abbrev>` is the user-confirmed or recorded cell abbreviation used for subtype prefixes, such as `epi`, `mye`, `t`, `b`, `endo`, or another project-specific abbreviation. Do not append `_annotated` or `_re` to this final annotated h5ad filename. Put lineage-specific figures and tables under `figures/<secondary-module>/<lineage>/` and `tables/<secondary-module>/<lineage>/`. Candidate parameter-search figures and tables may use parameter-coded subdirectories under `figures/03-subtype-clustering-grid/<lineage>/`, `tables/03-subtype-clustering-grid/<lineage>/`, and `codes/03-subtype-clustering-grid/<lineage>/`, but candidate parameter-search h5ad files are deferred until final user selection. Do not put lineage directories directly under `codes/`, `figures/`, `tables/`, or `h5ad/`.

Raw subtype-cluster columns such as `leiden_res<RES>` are technical clustering labels, not subtype labels. Do not directly copy `leiden_res<RES>` into `cell_subtype`. First use raw-cluster DEGs and marker evidence to assign marker-supported biological subtype labels or flag ambiguous groups for user review. Keep the raw `leiden_res<RES>` column separately for provenance.

For lineage-specific or subtype-level reclustering, the workflow is the same as the full-cell clustering workflow from Module 01, with one extra first step: extract the target cells. This rule applies identically to every observed and selected `leiden_coarse` value. Do not hard-code a fixed lineage set such as T, epithelial, myeloid, endothelial, and B cells. Instead, inspect the current project's validated `leiden_coarse` categories, report the observed lineages to the user, and run subtype reclustering only for the lineages that are present and selected. Use `adata_qc` as the expression/reclustering base. Use Module 01 broad annotation/scored h5ad objects only to identify target cell IDs and map confirmed annotation columns back by stable cell IDs. If a filtered annotation h5ad such as `adata_anno_score_genes_rank_consistent.h5ad` exists, use it by default to define the eligible IDs, then still subset `adata_qc` for the expression matrix. Then rerun the normal post-QC workflow on that subset: normalization/log handling as appropriate, HVG selection, scaling/regression if used, PCA, optional Harmony or other integration, neighbors, UMAP, Leiden, DEG calculation, and annotation. The subset must be based on `leiden_coarse` or another confirmed annotation label, not on raw global Leiden cluster numbers such as 1, 2, 3, or 4. Do not treat the multiple raw global clusters that map to one `leiden_coarse` label as subtypes. For example, if global raw clusters 10 and 11 both map to `leiden_coarse == 'Epithelial Cells'`, extract all cells with that `leiden_coarse` label and recluster the extracted epithelial subset; the same logic applies if several global raw clusters map to `T Cells`, `Myeloid Cells`, or any other observed `leiden_coarse` label. Do not call those broad raw clusters subtypes. Do not create the reclustering input by repeatedly subsetting already subsetted or integrated objects.

For subtype reclustering, the pre-Harmony operation order and default numeric parameters should match the full-cell Module 01 workflow. This includes normalization, log transform, HVG selection, raw-assignment point, HVG subsetting, regression if used, scaling, PCA, batch key, PCA basis, and Harmony input/output basis names. Recompute these steps on the subtype subset rather than reusing full-cell PCA/Harmony coordinates, unless the user explicitly asks for projection-only analysis. Subset-specific HVGs may differ because they are selected within the subset, but the HVG method and `n_top_genes` default to the Module 01 values. If Harmony does not converge for the subtype subset, increase `max_iter_harmony` as a same-method convergence adjustment and rerun Harmony for that lineage before changing integration methods. Do not merely edit the parameter table. Record the original and increased `max_iter_harmony`, the affected lineage, convergence warnings/logs, elapsed time, output h5ad path, and the final Harmony embedding name. Save the completed subtype Harmony template h5ad before starting the grid search, under the current module's `h5ad/<secondary-module>/<lineage>/` directory, using a state-aware name such as `adata_<lineage>_harmony.h5ad`. The subtype clustering/grid step must then read this saved Harmony h5ad as its source object, mirroring Module 01 where clustering reads `h5ad/04-integration-harmony/adata_harmony.h5ad`. Do not continue directly from the in-memory object produced by the Harmony step. After the saved subtype Harmony h5ad is read, run the fixed default subtype grid over `n_pcs = 10, 15, ..., 50`, `n_neighbors = 10, 15, ..., 50`, and Leiden `resolution = 0.1, 0.2, ..., 1.5` unless the user explicitly provides another grid. For every `n_neighbors`/`n_pcs` graph candidate, use the reloaded subtype Harmony object as the template; in the rsc path transfer the template to GPU before the loop, copy once per graph candidate, run all default resolution values with Leiden by default or Louvain after abnormal Leiden fallback, save one combined UMAP review figure plus lightweight tables, then delete the candidate copy and release memory. Do not save any candidate h5ad during this grid/candidate phase. After the user selects parameters, perform only the subtask 04-subtype-selected-clustering rerun from the saved subtype Harmony h5ad, using GPU/RAPIDS/rsc when available under the same backend rule as the grid search; save the selected UMAP and selected clustered h5ad with a state-aware selected-clustered filename, and stop before DEG/annotation. Then subtask 05-subtype-deg-annotation reads the selected clustered h5ad and performs DEG calculation, annotation, final subtype plotting, and annotated lineage h5ad saving as `adata_<cell_abbrev>.h5ad` unless the user explicitly asks to annotate multiple selected candidates.

Default DEG-assisted annotation workflow:

```text
1. Subtask 05-subtype-deg-annotation starts by reading the selected clustered lineage h5ad written by
   `04-subtype-selected-clustering`. This h5ad must already contain the
   unannotated raw cluster column and selected UMAP from the user-selected run.
   For a Leiden subtype/subset clustering pass, name the raw cluster column
   `leiden_res<RES>`, such as `leiden_res0p5`. For Louvain, use
   `louvain_res<RES>`.
2. Compute DEGs/marker genes for the unannotated clusters before assigning subtype
   names. Use Scanpy t-test for all DEG calculations unless the user explicitly
   requests another method:
   sc.tl.rank_genes_groups(adata, groupby=col, method='t-test', use_raw=True)
   All DEG calculations must use `adata.raw` normalized/log expression explicitly
   with `use_raw=True`. If `adata.raw` is absent at a DEG step, stop and fix the
   upstream raw-normalized expression preservation before exporting DEGs.
3. The DEG output directory and CSV filename must use the exact `groupby=col`
   column name used in `sc.tl.rank_genes_groups`. For example, with
   `col = 'leiden_res0p5'`, save outputs as:
   tables/05-subtype-deg-annotation/<lineage>/degs_leiden_res0p5_pcs<PCS>_nn<NN>_res<RES>/
   {group}_degs_leiden_res0p5_pcs<PCS>_nn<NN>_res<RES>.csv
4. Export one full-length CSV per group/cluster. Do not save top-only DEG files
   such as top50, top100, top200, top300, or any `topXX` DEG CSV. Do not write a
   combined DEG CSV by default; create one only if the user explicitly requests
   it. Every per-group CSV must include the project-style canonical DEG columns:
   gene, score, logfoldchanges, pvals, pvals_adj.
   If the chosen DEG method emits additional documented statistics, preserve
   them after these canonical columns rather than dropping them.
5. Build each per-group CSV from `adata.uns['rank_genes_groups']`, for example:
   pd.DataFrame({
       'gene': result['names'][group],
       'score': result['scores'][group],
       'logfoldchanges': result['logfoldchanges'][group],
       'pvals': result['pvals'][group],
       'pvals_adj': result['pvals_adj'][group],
   })
6. Use the already-saved full DEG CSVs together with canonical lineage/subtype
   markers to interpret subtype identity; do not annotate from cluster number
   alone. By default, review the first 50 rows from each full DEG CSV. If those
   rows do not provide enough interpretable marker evidence for that annotation
   level, expand the whole annotation level sequentially to the first 100, 200,
   and 300 rows from the same full CSVs. Different annotation levels may use
   different review depths, but all groups at the same annotation level must use
   the same review depth. Record the selected review depth for each annotation
   level and why expansion was needed, but do not save separate top-only DEG CSVs
   or top-N DEG manifests.
   Before writing any final subtype labels, save a candidate annotation table and
   either obtain user confirmation or explicitly record that the user authorized
   agent-led final labeling. Then save the final mapping CSV and use it to write
   `cell_subtype` and `functional_state` into `.obs`. Assign subtype labels from
   marker-supported biological identity and functional state rather than preserving
   Leiden labels as artificial subtype names. For subtype naming, use the most
   significant positive DEG for that raw cluster as the gene suffix by default.
   "Most significant" means the first gene in the exported DEG table after
   ordering by adjusted P value when available, then raw P value, then Scanpy
   score; if the exported table is already sorted by the DEG method, use its
   first positive-logfoldchange gene. The final subtype mapping must be one-to-one
   from each selected raw cluster to one unique `cell_subtype`; do not assign the
   same `cell_subtype` to several raw clusters inside the same lineage. If two or
   more raw clusters would receive the same prefix-plus-gene label because their
   top positive DEG is identical, give the top-gene label to the cluster with
   stronger evidence for that gene, then use the next-ranked eligible positive
   DEG from each losing cluster's own DEG table. Record this duplicate-top-gene
   resolution in the final mapping CSV and candidate table. If no confirmed
   or agent-authorized mapping CSV exists, Subtask 05-subtype-deg-annotation must stop here and must not
   save a final annotated h5ad.
7. After final subtype labels are written to `.obs`, recompute DEGs using the final annotation
   column, such as `cell_subtype`, again using:
   sc.tl.rank_genes_groups(adata, groupby=col, method='t-test', use_raw=True)
8. For post-annotation subtype DEGs with `col = 'cell_subtype'`, save outputs as:
   tables/05-subtype-deg-annotation/<lineage>/degs_cell_subtype_pcs<PCS>_nn<NN>_res<RES>/
   {group}_degs_cell_subtype_pcs<PCS>_nn<NN>_res<RES>.csv
   If the DEG `groupby` column is `cell_type` or another
   documented column, replace the directory and filename token with that exact
   column name. For subtype-level annotation or refinement, review the first 50
   rows from each saved full DEG CSV first. If marker evidence is insufficient for
   that subtype-analysis level, expand the whole level sequentially to the first
   100, 200, and 300 rows from the same full CSVs. All groups within the same
   subtype-analysis level must use the same review depth, even if one group
   appears annotatable with fewer genes. Record the selected review depth for that
   subtype-analysis level, but do not save separate top-only DEG CSVs or top-N DEG
   manifests.
9. Projection back to the full atlas is a separate secondary task,
   `06-project-subtypes-to-full-adata`. Do not silently write projected labels
   into a full h5ad during `05-subtype-deg-annotation`; save the annotated
   lineage h5ad first, then run task 06. If DEGs are later computed at the
   full-atlas broad-label level, use the exact groupby column such as
   `col = 'cell_type'` and save under the task that performs that DEG run.
10. Draw marker dotplots and subtype UMAPs from the final annotation labels, using the
    selected marker panel and/or most significant DEG used in each subtype name.
```

Before finalizing subtype names, ask the user whether they want to manually confirm labels or allow the agent to assign labels directly. Also proactively ask whether the user wants to provide manual broad-lineage abbreviations for subtype prefixes, such as a preferred abbreviation for myeloid cells. This is especially important for new datasets, weak marker separation, ambiguous transitional clusters, or when labels will be used in publication figures. If the user chooses agent-led annotation, still provide marker evidence, cluster sizes, proposed prefixes, and any uncertainty flags. If the user chooses manual confirmation, produce a candidate label table and wait for confirmation before writing final `cell_subtype` columns.

Default cluster/subtype naming convention:

```text
raw subtype-cluster column = leiden_res<RES>
if using Louvain, raw subtype-cluster column = louvain_res<RES>
subtype label = <selected_prefix>_<GENESYMBOL>
```

`<GENESYMBOL>` should be the most significant positive DEG for that raw cluster, written in official uppercase gene-symbol style when applicable. Do not hand-pick a more biologically attractive marker if a valid most-significant DEG exists. If two or more subclusters would receive the same `<selected_prefix>_<GENESYMBOL>` label because their top positive DEG is identical, keep the top-gene label for the cluster with stronger evidence for that gene and assign the other cluster or clusters the next-ranked eligible positive DEG from their own DEG tables. Compare evidence using adjusted P value when available, then raw P value, then Scanpy score or absolute positive logfoldchange. Record the duplicate-top-gene resolution in `gene_selection_rationale`; do not invent alternate marker names outside the DEG ranking and do not merge clusters into one subtype label. Use one gene in the subtype name; do not create combined gene labels such as `<prefix>_<GENE1>_<GENE2>`, `<prefix>_<GENE1>-with-<GENE2>`, `<prefix>_<GENE1>+<GENE2>`, or any equivalent two-gene naming form unless the user explicitly requests that naming style. `<selected_prefix>` should be chosen by this rule:

```text
0. Before writing subtype labels, ask whether the user wants to supply manual
   broad-lineage abbreviations. Use the user-provided abbreviation table when
   available and record it in the candidate label table.
1. If the subcluster belongs to the current broad lineage being analyzed, use the
   current broad-lineage prefix.
2. If the user does not provide a prefix for an observed broad lineage, derive a
   conservative default from the broad label. Common defaults include Epi for
   epithelial cells, T for T cells, B for B cells, Endo for endothelial cells,
   Mye for myeloid cells, NK for NK cells, Mast for mast cells, and Stromal or
   Str for stromal cells. The myeloid prefix should be abbreviated as Mye rather
   than using the full word Myeloid in subtype labels, unless the user overrides
   it.
3. If marker evidence shows that the subcluster resembles another lineage or
   lineage subtype, keep annotating it inside the current lineage-specific h5ad
   using the current lineage prefix and the standard most-significant-positive-DEG
   naming rule. Never modify another lineage-specific h5ad or move cells between
   lineage h5ad objects during subtype annotation just because a subcluster has
   marker evidence resembling another lineage. Move or restore cells to another
   lineage only if the user explicitly requests that operation.
5. Always report the marker evidence, functional-state rationale, prefix source,
   and naming rationale to the user before finalizing subtype names.
```

Do not use any special suffix rule for these cases. A subcluster that resembles
another lineage still receives a normal prefix-plus-DEG subtype label within the
current lineage-specific object unless the user explicitly requests
moving/restoring those cells to another lineage.

Examples of the label pattern:

```text
Epi_<GENESYMBOL>
Mye_<GENESYMBOL>
Endo_<GENESYMBOL>
T_<GENESYMBOL>
```

If a different subtype naming style is supplied by the user, keep it only after documenting the mapping to this convention.

After subtype labels are assigned, regenerate or validate `cell_type` from the subtype prefix unless the user explicitly provides a trusted `cell_type` column or custom mapping. If an off-lineage marker pattern is retained under the current broad-lineage prefix, add a note column such as `annotation_note`, `marker_evidence_lineage`, or `prefix_rationale` so the user can see what the marker evidence suggested.

For every subtype naming pass, provide a candidate table and a final mapping CSV.
The final mapping CSV used to write annotations must contain at least these
columns:

```text
cluster
cell_subtype
functional_state
gene_selection_rationale
```

Here `cluster` is the raw Leiden or Louvain cluster label from the selected
subtype-clustering run, `cell_subtype` is the final subtype name, and
`functional_state` is a concise biological state/function assigned from the same
marker review, such as metabolic, stress-response, mesenchymal-like,
proliferative, cytokine/chemokine, immune-like, off-lineage-like, or another
project-specific state. `gene_selection_rationale` must state that the gene in
`cell_subtype` is the most significant positive DEG for that raw cluster. The
`cell_subtype` values must be unique within the current lineage-specific
annotation table.

The candidate table should include:

```text
cluster
proposed cell_subtype
functional_state
selected_prefix
most_significant_deg_gene
gene_selection_rationale
top supporting markers
marker-supported lineage
prefix rationale
user_modified_marker when applicable
confidence or needs_user_confirmation
```

## Default Use of Module 01 Marker-Score Evidence

Use this step by default when Module 01 has already produced broad marker-score/rank evidence from `leiden_coarse` DEGs and the agent needs a cleaner lineage-specific input for subtype reclustering. Module 02 should consume these existing columns; do not recompute the broad `leiden_coarse` score/rank workflow here unless the Module 01 output is missing and the user explicitly asks to regenerate it.

```text
1. Validate that Module 01 columns are present, especially `leiden_coarse`.
2. Use `adata_qc` as the expression/reclustering base. Use the Module 01 broad
   annotation or scored h5ad to choose target IDs and transfer `leiden_coarse`
   plus score/rank evidence by exact cell ID when available.
3. If `adata_anno_score_genes_rank_consistent.h5ad` exists, use it by default
   as the annotation source. Select cells from it by `leiden_coarse`, then filter
   `adata_qc` by those selected IDs.
4. If the consistent object is absent and score/rank filtering is explicitly not
   required, use the final integrated annotated object such as `adata_inte.h5ad`
   only as the annotation/ID source; still subset `adata_qc` for expression.
5. If the input h5ad is a legacy unfiltered score/rank object, confirm that
   `best_rank_type_global` uses the same broad-label vocabulary as `leiden_coarse`,
   or load the saved Module 01 score-column mapping before applying any consistency
   rule.
6. If building a lineage-specific object for downstream subtype clustering,
   use marker-rank consistency filtering by default when score/rank evidence is
   available. Treat this as required project QC unless the user explicitly
   cancels it, and document whether it was used or cancelled.
```

Default project extraction rule when the consistent h5ad exists:

```text
annotation source = adata_anno_score_genes_rank_consistent.h5ad
select IDs where obs['leiden_coarse'] == '<target_leiden_coarse>'
subset adata_qc by those IDs
do not select by the raw global Leiden clusters that were merged into <target_leiden_coarse>
```

Coarse-label-only extraction rule when the user explicitly cancels score/rank filtering:

```text
annotation source = adata_inte.h5ad or another final integrated annotated h5ad
expression base = adata_qc.h5ad
keep cells for the requested subtype reclustering:
obs['leiden_coarse'] == '<target_leiden_coarse>'
subset adata_qc by those IDs
```

Legacy unfiltered input consistency rule:

```text
keep cells for the requested subtype reclustering:
obs['leiden_coarse'] == '<target_leiden_coarse>'
and mapped obs['best_rank_type_global'] == '<target_leiden_coarse>'

keep mast cells for mast subclustering:
obs['leiden_coarse'] == 'Mast Cells'
and mapped obs['best_rank_type_global'] == 'Mast Cells'
```

The same extraction and consistency rules must be used for every broad label selected from `leiden_coarse`; do not invent a different selection rule for one lineage unless the user explicitly requests and documents it.

The equality rule above is the default project QC rule. Do not replace it with no filtering, rank-percentile thresholds, score-margin thresholds, doublet-aware filtering, or manual retention of biologically plausible transitional cells unless the user explicitly requests that change and the run record documents it.

Legacy project objects may instead contain compact labels such as `Epi`, `Endo`, `Mast`, `Myeloid`, `NK`, `Stromal`, `T_cell`, and `B_cell` in `best_rank_type_global`. If these are present, map them to the coarse-label vocabulary before applying the consistency filter.

This is not the same as ordinary count/mitochondrial QC. Treat the Module 01 score/rank columns as annotation-confidence evidence, and treat any later cell removal as the default lineage-specific filtering step unless the user explicitly cancels it. Keep these artifacts when available:

```text
unfiltered annotated h5ad
filtered lineage-specific h5ad
obs['highest_score_type']
obs['best_rank_type_global']
all *_score, *_score_rank, and *_score_rank_pct columns when available
filter mask or a table with kept/removed counts per lineage if filtering is used
UMAP comparing assigned label versus marker-rank evidence
```

If the user supplies an already annotated h5ad and asks to skip this step, do not force it. If Module 01 score/rank evidence is inspected but no filtering is applied, report that no filtering was applied. If filtering is skipped entirely, continue with the supplied labels after validation. If filtering is applied before lineage/subtype clustering, report the rule, how many cells were removed, and whether downstream abundance-based analyses use filtered or unfiltered cells.

## 06 Project Subtypes To Full Adata

If a small/lineage-specific h5ad was annotated separately, project labels back to the full integrated atlas:

```text
1. Keep stable cell IDs between subset and full object.
2. Transfer subtype labels by exact obs_names match.
3. If obs_names changed, use sample + original_barcode mapping.
4. Do not overwrite existing full-object labels without saving the old column.
5. Record unmatched cells and duplicated mappings.
```

Recommended columns:

```text
cell_type
leiden_res<RES> or louvain_res<RES>
cell_subtype
functional_state
annotation_source
annotation_confidence or marker_evidence when available
```

Task 06 must read the final annotated lineage h5ad from
`h5ad/05-subtype-deg-annotation/<lineage>/adata_<cell_abbrev>.h5ad` and project
its labels back onto the same full-atlas annotation source used to select the
lineage cells for subtype analysis. If subtype reclustering used IDs from
`adata_anno_score_genes_rank_consistent.h5ad`, project labels back to that
consistent object. If subtype reclustering used a broad-label-only object
because the consistent object was absent or filtering was explicitly cancelled,
project labels back to that corresponding broad-annotated h5ad instead. Do not
mix a lineage subset derived from one full-atlas h5ad with a different
projection target unless the user explicitly requests it; record the source
h5ad, target h5ad, and reason. Then write outputs only under
`06-project-subtypes-to-full-adata`. Save the projected full h5ad as
`h5ad/06-project-subtypes-to-full-adata/adata_anno_cellsubtype.h5ad`,
a projection match report under
`tables/06-project-subtypes-to-full-adata/<lineage>/`, and code under
`codes/06-project-subtypes-to-full-adata/<lineage>/`. Projection QC figures are
required outputs, not optional. Save them under
`figures/06-project-subtypes-to-full-adata/` using the official Scanpy plotting
interface and the global Scanpy plotting rules. Draw one full-atlas comparison
UMAP containing both the original broad label `leiden_coarse` and the projected
broad label `cell_type`, named
`umap_leiden_coarse_vs_projected_cell_type.pdf` and, when possible, the matching
SVG. Use one shared broad-label palette for both panels so that unchanged broad
cell identities have the same colors in `leiden_coarse` and projected
`cell_type`. Draw the projected subtype label `cell_subtype` as its own separate
full-atlas UMAP, named `umap_projected_cell_subtype.pdf` and, when possible, the
matching SVG; do not combine `cell_subtype` into the broad-label comparison
figure. For each lineage-specific annotated h5ad that was projected, also draw
one lineage-level UMAP under
`figures/06-project-subtypes-to-full-adata/<lineage>/`, colored by
`cell_subtype`, using the same `cell_subtype` palette as the full-atlas
`umap_projected_cell_subtype` figure. The lineage-level `cell_subtype` plot is
intended to show the projected subtype colors on the corresponding lineage h5ad;
do not replace it with a broad `cell_type` or `leiden_coarse` plot. Record both
palette mappings, UMAP basis, source h5ad, target h5ad, and figure paths in the
projection match report or a companion figure manifest. Preserve the old
full-atlas columns by copying them to backup columns before writing new
projected labels, unless the user explicitly permits overwriting without backup.


## Validation

- Module 01 broad labels are present before subtype annotation starts.
- Subtype labels are internally consistent within each lineage.
- Any lineage-specific subclustering satisfies the global tight/clear/stable clustering QC standard.
- If post-annotation marker-score QC is used, filtered cells match both assigned label and marker-rank evidence, and the unfiltered atlas is still retained.
- Raw subtype-cluster columns follow `leiden_res<RES>` or `louvain_res<RES>` depending on the clustering algorithm.
- Subtype names follow `<selected_prefix>_<GENESYMBOL>` unless a documented mapping is provided.
- Final `cell_type` is consistent with the `cell_subtype` prefix or with a user-provided custom mapping.
- Projection/remapping preserves cell counts and reports unmatched cells.
- The final full h5ad contains the labels required by downstream skills: epithelial-state, CM-Epi, CellChat, and survival.

## 03-Epi-CM Discovery

This block is the former Module 04, copied as the canonical instruction source
for the compact workflow. Inside this big skill, write this block's outputs under
`epi-cm-core-workflow/{codes,h5ad,tables,figures}/03-epi-cm-discovery/`
unless the user explicitly asks to use the original numbered module output tree.
Do not use the shorter summary as a substitute for these copied rules.

# 04-project-cm-lineage-core

Use this skill for the core CM-lineage coupling analysis. If a example branch exists, treat its selected lineage and figure naming as examples only. For any project, set the selected lineage, status groups, CM naming convention, marker sets, and biological interpretation while keeping the same matrix/provenance contracts. This module can be entered independently from annotated h5ad, subtype abundance matrices, precomputed CM matrices, or figure-ready CSVs depending on what the user provides, but it is not logically isolated from the other modules. It should reuse the same `sample_id`, `status`, `cell_type`, `cell_subtype`, and subtype-prefix mapping contracts created or validated by the single-cell and annotation modules. Do not force upstream single-cell integration, annotation, spatial analysis, or survival analysis when valid intermediates are supplied.

Terminology:

```text
CM = a latent co-occurrence module/program learned from non-epithelial cell-subtype
     abundance patterns across samples by the balanced joint NMF workflow.

CM is not a cell type, not a cell subtype, and not a raw non-epithelial subtype
column. A non-epithelial subtype can be a node/loading contributor to one or more
CMs, but it should not itself be called a CM.

CM activity = sample-level W value from the NMF decomposition, saved as a
sample x CM activity matrix.

CM loading = subtype-level H value from the NMF decomposition, saved as a
non-epithelial subtype x CM or CM x non-epithelial subtype loading matrix
depending on the table orientation explicitly recorded in the filename.

Epi-CM association = association between epithelial subtype abundance and CM
activity across samples.
```

## Block 03 Sample-Status Preflight and Additive Single-Status Routes

Run this preflight after the sample-level metadata and `keep_for_cm` eligibility
table have been constructed, but before NMF, CM classification, association, or
plotting. Route from the final eligible sample set; if filtering removes an
entire status, the detected mode must reflect the remaining eligible samples.
This is an additive routing layer. The existing tumor-plus-normal-like balanced joint NMF, CM
classification, canonical naming, association, and plotting route below must
remain unchanged when both status groups are present.

The status mode must be detected from one row per biological sample, never from
cell counts or filenames alone. Use the trusted sample metadata `status` field.
If `status` is absent but the user explicitly states that every supplied sample
is tumor or every supplied sample is normal-like, create that single canonical
status for the supplied sample list and record `status_source` as
`user_declared_tumor_only` or `user_declared_normal_only`. Do not silently make
either assignment without that statement.

Required routing modes:

```text
tumor_normal:
  both canonical statuses tumor and normal-like are present
  -> run the existing Block 03 route below without changing it

tumor_only:
  every eligible biological sample is tumor
  -> run the additive tumor-only route defined in this section

normal_only:
  every eligible biological sample is normal-like
  -> run the additive normal-only route defined in this section

unsupported:
  unresolved/unknown labels, conflicting statuses within one sample, or more
  than two biological status classes without a user-approved map
  -> stop before NMF and request/repair the sample-status mapping
```

Before routing, write both of these preflight outputs under
`tables/01-cm-lineage-analysis/01_prepare_inputs_and_frequency_tables/`:

```text
cm_status_mode_detection.csv
cm_status_mode.json
```

`cm_status_mode_detection.csv` must contain one row per detected canonical
status and at least `detected_mode`, `canonical_status`, `n_samples`,
`n_total_samples`, `status_source`, `status_balance_applied`, and
`cm_classification_available`. The JSON must also record the raw-to-canonical
status mapping and the complete included sample list. Do not start K selection
until the sum of `n_samples` equals the number of unique eligible sample IDs.

Use this explicit detection pattern. Project-specific status aliases remain
replaceable, but the route decision and audit outputs are fixed:

```python
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd


def _status_token(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(value).strip().lower()).strip()


def detect_block03_status_mode(
    sample_status: pd.DataFrame,
    out_dir: Path,
    *,
    sample_col: str = "sample",
    status_col: str = "status",
    status_source: str = "sample_metadata",
    user_declared_tumor_only: bool = False,
    user_declared_normal_only: bool = False,
    status_aliases: dict[str, str] | None = None,
) -> tuple[pd.DataFrame, str]:
    """FIXED: detect Block 03 mode from unique biological samples."""
    frame = sample_status.copy()
    if sample_col not in frame.columns:
        if frame.index.name is None:
            raise ValueError(f"Missing {sample_col!r} column and named sample index.")
        frame = frame.reset_index().rename(columns={frame.index.name: sample_col})

    frame[sample_col] = frame[sample_col].astype(str)
    if frame[sample_col].duplicated().any():
        if status_col not in frame.columns:
            raise ValueError("Duplicate sample rows require an explicit status column.")
        conflicts = frame.groupby(sample_col)[status_col].nunique(dropna=False)
        conflicts = conflicts[conflicts > 1]
        if len(conflicts):
            raise ValueError(
                "Conflicting status labels within samples: "
                + ", ".join(conflicts.index.astype(str))
            )
        frame = frame.drop_duplicates(sample_col, keep="first")

    if status_col not in frame.columns or frame[status_col].isna().all():
        if user_declared_tumor_only == user_declared_normal_only:
            raise ValueError(
                "Missing sample status. Supply trusted metadata or explicitly declare "
                "exactly one of tumor-only or normal-only."
            )
        if user_declared_tumor_only:
            frame[status_col] = "tumor"
            status_source = "user_declared_tumor_only"
        else:
            frame[status_col] = "normal-like"
            status_source = "user_declared_normal_only"
    elif frame[status_col].isna().any():
        missing_samples = frame.loc[frame[status_col].isna(), sample_col].astype(str)
        raise ValueError("Missing status for samples: " + ", ".join(missing_samples))

    # REPLACEABLE: extend only from trusted project metadata.
    aliases = {
        "tumor": "tumor",
        "tumour": "tumor",
        "primary tumor": "tumor",
        "primary tumour": "tumor",
        "metastatic tumor": "tumor",
        "metastatic tumour": "tumor",
        "metastasis": "tumor",
        "normal": "normal-like",
        "normal like": "normal-like",
        "adjacent normal": "normal-like",
        "normal adjacent": "normal-like",
    }
    if status_aliases is not None:
        aliases.update({_status_token(k): v for k, v in status_aliases.items()})

    frame["raw_status"] = frame[status_col].astype(str)
    frame["canonical_status"] = frame["raw_status"].map(
        lambda value: aliases.get(_status_token(value))
    )
    unresolved = frame.loc[frame["canonical_status"].isna(), "raw_status"].unique()
    if len(unresolved):
        raise ValueError(
            "Unresolved sample status labels; provide a trusted status_aliases map: "
            + ", ".join(map(str, unresolved))
        )

    observed = set(frame["canonical_status"])
    if observed == {"tumor", "normal-like"}:
        mode = "tumor_normal"
        status_balance_applied = True
        cm_classification_available = True
    elif observed == {"tumor"}:
        mode = "tumor_only"
        status_balance_applied = False
        cm_classification_available = False
    elif observed == {"normal-like"}:
        mode = "normal_only"
        status_balance_applied = False
        cm_classification_available = False
    else:
        raise ValueError(f"Unsupported Block 03 status set: {sorted(observed)}")

    frame[status_col] = frame["canonical_status"]
    counts = (
        frame.groupby("canonical_status", observed=True)[sample_col]
        .nunique()
        .rename("n_samples")
        .reset_index()
    )
    counts.insert(0, "detected_mode", mode)
    counts["n_total_samples"] = frame[sample_col].nunique()
    counts["status_source"] = status_source
    counts["status_balance_applied"] = status_balance_applied
    counts["cm_classification_available"] = cm_classification_available

    out_dir.mkdir(parents=True, exist_ok=True)
    counts.to_csv(out_dir / "cm_status_mode_detection.csv", index=False)
    with (out_dir / "cm_status_mode.json").open("w", encoding="utf-8") as handle:
        json.dump(
            {
                "detected_mode": mode,
                "n_total_samples": int(frame[sample_col].nunique()),
                "status_source": status_source,
                "status_balance_applied": status_balance_applied,
                "cm_classification_available": cm_classification_available,
                "raw_to_canonical_status": dict(
                    frame[["raw_status", "canonical_status"]]
                    .drop_duplicates()
                    .itertuples(index=False, name=None)
                ),
                "included_samples": frame[sample_col].astype(str).tolist(),
            },
            handle,
            indent=2,
            ensure_ascii=False,
        )
    return frame, mode
```

For `tumor_normal`, continue into the original two-status code below exactly as
written. For `tumor_only` or `normal_only`, preserve the same frequency
construction, column-wise min-max input, K candidate range, seed grid, NMF
implementation, best-seed selection, and NNLS activity refit, but use equal row
weights and do not call `make_status_balanced_weights()`, `classify_raw_cms()`,
or `assign_canonical_cm_names()`.

Keep `CM1`, `CM2`, ..., `CMK` only as internal raw component IDs. Preserve the
original `s_CM* / n_CM* / t_CM*` naming schema by assigning the prefix that
matches the only observed status: tumor-only uses final IDs `t_CM1`, `t_CM2`,
..., `t_CMK` with `class = "tCM"`; normal-only uses final IDs `n_CM1`,
`n_CM2`, ..., `n_CMK` with `class = "normalCM"`. Use the same global
raw-component order as the original two-status route. Also set
`classification_available = False` and
`classification_basis = "single_status_presence_fallback"` with a reason
naming the absent comparison status. This records that the prefix follows the
only observed status rather than a tumor-vs-normal enrichment test. Write the
corresponding raw-to-final mapping in `raw_to_canonical_CM_mapping.csv`.

Reuse the original two-status route's file, matrix, column, and status-context
names for every output that is meaningful in a single-status cohort. Do not
invent parallel `tumor_only_*`, `normal_only_*`, or `single_status_*` names for
analysis matrices or figures. The only single-status-specific filename is the
skip audit `single_status_skipped_outputs.csv`.

```text
both single-status routes:
  W_df.csv
  H_df.csv
  activity_df_sample_by_CM.csv
  activity_df_CM_by_sample.csv
  loading_df_cell_subtype_by_CM.csv
  loading_df_cell_subtype_by_CM_fraction.csv
  joint_module_classification.csv
  raw_to_canonical_CM_mapping.csv

tumor_only reuses the existing tumor-context names:
  tumor_network_nodes_from_H_df.csv
  tumor_all_CM_nodeplot.pdf/.svg
  balanced_joint_cm_epi_cm_association_tumor_rho_matrix.csv
  balanced_joint_cm_epi_cm_association_tumor_q_matrix.csv

normal_only reuses the existing normal-like-context names:
  normal_like_network_nodes_from_H_df.csv
  normal_like_all_CM_nodeplot.pdf/.svg
  balanced_joint_cm_epi_cm_association_normal-like_rho_matrix.csv
  balanced_joint_cm_epi_cm_association_normal-like_q_matrix.csv
```

Keep the original `joint_module_classification.csv` schema. For a single-status
cohort, populate the corresponding existing mean column (`tumor_mean` or
`normal_like_mean`) and leave the absent-status mean, `delta`, `p`, and `q` as
`NA`; do not replace them with a new `single_status_mean_usage` column.

For either single-status route, comment out or bypass every analysis/plotting
call that compares statuses. Do not call a two-status function and catch its
error afterward. Keep only the outputs for the status that actually exists:

```text
KEEP IN BOTH SINGLE-STATUS ROUTES:
  K-selection diagnostics
  sample x CM activity and CM x subtype loading tables/heatmaps

KEEP IN tumor_only:
  tumor node tables, tumor-only nodeplots, and tumor node-correlation heatmaps
  tumor epithelial-subtype x CM associations, q-star heatmaps, and all scatters

KEEP IN normal_only:
  normal-like node tables, normal-only nodeplots, and normal node-correlation heatmaps
  normal epithelial-subtype x CM associations, q-star heatmaps, and all scatters

SKIP AS STRUCTURALLY INAPPLICABLE:
  status-balanced weighting
  sharedCM/normalCM/tCM classification
  tumor-vs-normal-like CM activity barplot
  every nodeplot, edge table, or correlation heatmap for the absent status
  every Epi-CM heatmap/scatter branch for the absent status
  every edge-origin shared-versus-status-specific classification
```

For the retained single-status association branch, compute BH-FDR across all
tested `epithelial subtype x CM` pairs within that one status and method
context. Do not use a per-CM or per-epithelial-subtype correction family. Write
`single_status_skipped_outputs.csv` with `output_or_step`,
`status = "skipped"`, and a reason naming the absent comparison status. These
skips are valid completion states and must not be reported as missing outputs.

## Pre-Execution Plan Requirement

Before executing code from this skill, write a concise method-and-result plan that the user can review and copy as the goal. Keep it result-oriented rather than overly procedural. Include only:

```text
analysis goal / expected result
method route to use
main inputs or provided intermediates
major code modules to run or skip
expected output figures/tables
key validation criterion
```

Do not start long-running analysis, dependency installation, or file-rewriting steps until this short plan has been stated. For simple inspection-only tasks, one or two sentences are enough.

If the user does not provide a manual choice for parameters, thresholds, method options, output naming, or optional branches, use the documented default settings in this skill and state that the default was used.

## Canonical Plotting Contract

For Module 03 and all later modules, the final plotting code and style rules written in this skill are the source of truth. Use the code, parameters, palettes, layout rules, and output names specified here or in module-local code files generated from this skill. Do not tell a future agent to inspect old notebooks or external project paths during execution. If both an analysis route and a later redraw/plot route are encoded for the same figure, the later redraw/plot implementation is canonical. Do not improvise alternate plot types, palettes, layouts, statistical labels, file formats, or single-pair shortcuts. If a required final figure lacks explicit plotting code or style rules in this skill, stop and ask for the skill to be updated before running; do not invent a new plotting route.

For CM analysis figures, later dedicated plotting/redraw implementations are
canonical over plotting embedded in the analysis code. In particular, use the
dedicated final plotting code for:
`joint_nmf_k_selection`, `activity_df_tumor_vs_normal_mean_sd_barplot`,
`w_df_activity_sample_activity_per_CM_standard_scale_col_with_series_clustermap`,
`tumor_centric_nodeplot_edge_origin`, top-node correlation heatmaps, and all
Epi-CM scatter plots. If an analysis function and a dedicated plotting function
would write the same figure stem, keep only the dedicated plotting output.

## Project Organization and Figure Output Contract

Treat each numbered module folder as its own output boundary with one shared four-directory layout: `figures/`, `tables/`, `codes/`, and `h5ad/`. Module output directory names must use the active project slug, for example `04-<project_slug>-cm-lineage-core/`; for BRCA use `04-brca-cm-lineage-core/`. The four top-level category directories may be created at module setup. Secondary-module/task directories inside those category directories must be created only when that category will receive at least one real output for that task. Do not pre-create empty task directories under `figures/`, `tables/`, `codes/`, or `h5ad/` just to mirror the layout. A directory creation command for a secondary task or candidate must be coupled to writing a real output file there; if the output is not generated, do not leave that task directory behind. Secondary modules are logical analysis units inside that numbered module, but they do not require a directory under every category. This is a write-location rule, not a read restriction: a module may read/reuse files and already generated outputs from other modules as inputs, but newly generated outputs for the current module must be written inside the current numbered module. Across the whole skill workflow, an agent must not delete, clear, overwrite, or move any existing output file or directory anywhere unless the user explicitly names the exact path and operation. During normal module execution, do not write, move, overwrite, clear, or delete files or directories inside any other numbered module output directory. Also do not delete, clear, overwrite, or move any existing output file or directory inside the current module unless the user explicitly names the exact path and operation. If a new result would conflict with an existing output, write to a new versioned path or stop and ask. Deleting any output directory is never part of a module run; it requires a separate explicit cleanup request naming the exact path. For example, Module 02 subtype outputs go under `02-<project_slug>-cell-annotation/{figures,tables,codes,h5ad}/01-epithelial-subclustering/`, not under `01-<project_slug>-singlecell-integration/`, not directly as files under `02-<project_slug>-cell-annotation/`, and not under the top-level skill source directory.

Use stable numbered secondary-module/task names that describe the analysis step, lineage, method, or figure group. Reuse the same secondary-module/task name only under category directories that actually receive outputs from that analysis, so files stay aligned without creating empty placeholder task directories. For example:

```text
02-<project_slug>-cell-annotation/
  codes/
    01-epithelial-subclustering/
      02_epithelial_subclustering.ipynb
      run_epi_subclustering.py
  h5ad/
    01-epithelial-subclustering/
      adata_epi.h5ad
  figures/
    01-epithelial-subclustering/
      umap_cell_subtype.pdf
  tables/
    01-epithelial-subclustering/
      epithelial_subtype_counts.csv
```

By default, save executable/reproducibility code under the current module's shared `codes/<secondary-module>/`, using ordered names such as `01_read_merge.ipynb`, `02_qc.py`, or `03_integrate.R`. Save AnnData-like objects under `h5ad/<secondary-module>/` as `.h5ad`, `.loom`, `.rds`, or equivalent files with stable names. Save corresponding figure files under `figures/<secondary-module>/` and use ordered names such as `01_umap.pdf` or `02_marker_dotplot.svg`. Save text-like and tabular outputs under `tables/<secondary-module>/`, such as CSV/TSV/XLSX/TXT/JSON/YAML logs, manifests, reports, mapping files, and parameter records. `figures/` should contain figure files only. `tables/` should contain text-like and tabular outputs only. `codes/` should contain executable/reproducibility code only. `h5ad/` should contain AnnData-like/intermediate object files only. Add `tables/<secondary-module>/readme.txt` documenting the input files, including any cross-module input/output files that were read, code order, h5ad/loom/rds objects, output figures/tables, and any skipped optional branches.

Do not write new h5ad, code, figures, or tables directly into the numbered module root. The numbered module root may contain the module `SKILL.md`, lightweight module-level index files, or manually curated high-level notes, but executable outputs should live under the shared four category directories. If a simple task has only one natural step, still use a small secondary-module/task name such as `01-main`, `01-qc`, or `01-epithelial-subclustering`, but create that task directory only under the category directories that receive real outputs.

If one analysis step outputs multiple files or figures that belong to one coherent output set, put that output set in the same named secondary-module subdirectory under `figures/`, `tables/`, `codes/`, or `h5ad/`, using the same analysis prefix when possible. If one code file produces multiple distinct output sets, create one additional child subdirectory per output set inside the relevant secondary-module directory instead of flattening unrelated outputs into one directory. For example, a single plotting script that writes association heatmaps, scatter plots, and node correlation heatmaps should write to `figures/02-cm-lineage-final-plotting/epi-cm-heatmaps/`, `figures/02-cm-lineage-final-plotting/epi-cm-scatterplots-spearman/`, and `figures/02-cm-lineage-final-plotting/node-correlation-heatmaps/`, with matching `tables/` child directories for manifests or source tables when produced. The code file itself may remain in `codes/02-cm-lineage-final-plotting/`, but its output paths must make the output-set boundary explicit.

If an output already exists, do not rerun only to recreate it in the new layout. Do not move or delete existing outputs for layout cleanup unless the user explicitly names the exact path and operation. Prefer to leave existing outputs in place, copy them into the organized location only when provenance is recorded, then update the corresponding code paths so future runs write to the same organized location.

When a task creates a run, lineage, candidate parameter set, or method variant, create matching candidate subdirectories inside the active secondary-module/task directory only under category parents that receive outputs for that candidate. For multi-candidate or multi-condition runs, use matching candidate names under each relevant parent when needed; for example, create `h5ad/01-integration-parameter-search/pcs-30_nn-15_res-0p8/` only if an AnnData-like object will be saved, and create `figures/01-integration-parameter-search/pcs-30_nn-15_res-0p8/` only if figures will be saved. Keep h5ad-like candidate objects under the module's shared `h5ad/` with secondary-module and parameter-coded subdirectories, for example `h5ad/01-integration-parameter-search/pcs-30_nn-15_res-0p8/adata_inte.h5ad`; keep candidate code files under the shared `codes/` with matching secondary-module and parameter-coded paths when code is emitted.

Each analysis that produces an output should have corresponding source code under the current module's `codes/`. Acceptable code artifacts include `.ipynb`, `.py`, `.R`, and `.sh`, depending on the language actually used. Do not leave a figure, table, or exported result that can only be traced to manual GUI editing. If an analysis uses Python, keep the notebook and/or `.py` script that generates it; if it uses R, keep the `.R` script or R notebook; if both languages are used, keep both code artifacts under `codes/` with clear ordered prefixes. When converting notebooks to upload/download versions, keep the executable cells needed to reproduce the outputs and remove stale display output only when requested.

Each executed run should also create or update a parameter/provenance report under `tables/`, such as `tables/run_parameters.txt`, `tables/run_parameters.csv`, or a step-specific report in the same output subdirectory. The report should list the code file used, input files/objects, output files, exact parameters, random seeds, selected candidate/final settings, skipped steps, fallback decisions, and any user-approved method changes.

Do not substitute another analysis method, algorithm, statistical test, visualization strategy, database, or input layer without explicit user permission. If the specified method cannot run, stop that module, document the blocker in `tables/readme.txt`, and ask for confirmation before using any alternative. Any approved or documented method change should state why the original method was unsuitable or failed and why the replacement method is appropriate for the same analysis goal.

When a task, notebook run, script run, or long interactive kernel finishes, promptly close the process/kernel/session and release CPU memory and GPU memory. Do not leave idle Python, R, Jupyter, CellChat, RAPIDS, PyTorch, TensorFlow, or CUDA processes holding RAM/VRAM after the requested work is complete.

After each module finishes, create or update `tables/package_versions.txt` describing the packages and tools used by that module. Include Python packages, R packages, command-line tools, CUDA/GPU libraries when relevant, interpreter/R version, environment name or path, and the code files that used them.

Install missing dependencies when they are required to execute the specified method or its approved acceleration path. This includes installing a compatible GPU-accelerated implementation when the method supports it and the machine has a usable GPU/CUDA driver, for example installing `rapids-singlecell`/RAPIDS to run Scanpy-style preprocessing through `rsc`. Dependency installation is allowed to make the requested method work; method substitution is not allowed without explicit user permission. For packages or methods that already provide GPU acceleration, enable and use the GPU-accelerated path after installing any missing compatible GPU packages and verifying imports/minimal execution. If the requested package/method has no GPU-accelerated implementation, use its normal CPU path. If an expected GPU path is installed but broken or incompatible for a non-OOM reason, including a CUDA/CUDA-tag mismatch such as `cu11` vs `cu12` wheels, CuPy/RAPIDS/PyTorch wheels incompatible with the visible driver, missing CUDA runtime libraries, `libucx`/UCX errors, or `cuCtxGetDevice`/CUDA context errors, first try to repair or reinstall a compatible GPU environment without changing the requested method. Choose a compatible wheel, channel, or uv environment automatically from `nvidia-smi`, Python version, platform, and package compatibility information; do not ask the user to choose the CUDA tag. Ask the user only before system-driver changes, OS package changes that require elevated privileges, deleting an existing environment, or replacing a working environment used by other analyses. If GPU runs out of memory, inspect active GPU processes, close stale or idle processes left by previous tasks/kernels when they can be safely identified, release VRAM, record the OOM, and switch that OOM step directly to the equivalent CPU implementation; do not repair, reinstall, or repeatedly retry GPU solely for OOM. Do not terminate unrelated active user processes unless the user explicitly approves. If compatible dependency/GPU installation or non-OOM repair fails, or no usable GPU is present, document the reason and continue with the normal CPU path for the same requested method.

This GPU backend rule applies to all GPU-capable code in every module. If the first full execution of the required task completes successfully with the planned backends, do not rerun only to validate the backend plan. If a GPU-accelerated step fails from GPU OOM, release VRAM, record the OOM, and treat CPU fallback for that step as pre-approved. If a GPU-accelerated step fails for a non-OOM reason and the user approves CPU fallback after documented repair/cleanup attempts, treat fallback at step granularity: use the CPU equivalent only for the failed step when needed, then continue later steps with GPU whenever those later steps have a valid GPU path. Do not mark the entire remaining workflow CPU-only because one GPU step failed. Record a backend capability table under the relevant `tables/<secondary-module>/` or `tables/<secondary-task>/` directory with one row per executed step, including `step`, `planned_backend`, `attempted_backend`, `status`, `error_summary`, `fallback_backend`, `clean_input_reloaded`, and `final_backend_for_rerun`. Also export `gpu_backend_capability_summary.csv` and, when useful, `gpu_backend_capability_summary.txt` in that same tables directory; these files must state which steps can use GPU, which steps must use CPU, and the reason for each CPU step. When any GPU failure, fallback, or partial object mutation occurs during this exploratory/profiling pass, finish the required task only to learn which steps can use GPU, then start a fresh Python/R process, reload the nearest clean upstream input h5ad/RDS/input files from disk, and rerun the whole required task once using the recorded final backend plan. In that final rerun, every step marked GPU-capable must use GPU, except steps with recorded GPU OOM or approved non-OOM CPU fallback, which should use the recorded CPU fallback so the final run avoids repeating known mid-run GPU failures while still using GPU wherever it works. Do not reuse in-memory AnnData/R objects, arrays, GPU buffers, fitted models, graphs, or partial metadata from the profiling pass. Do not present outputs from a partial-failure/profiling pass as canonical final outputs.

Python package management rule: use `uv` by default for Python dependencies, and save uv-managed environments under a dedicated subdirectory in the total analysis directory for the dataset/project. Create the layout `uv_envs/<category>/.venv` under the top-level analysis root, where `<category>` is a stable dependency category such as `main`, `rapids`, `velocity_cellrank`, `cellchat_liana`, or `survival`. Use `uv_envs/main/.venv` for the default shared Python stack, and create another category only when dependency compatibility requires it. Install with `uv pip install --python uv_envs/<category>/.venv/bin/python --index-url https://pypi.tuna.tsinghua.edu.cn/simple ...` or activate that category `.venv` before `uv pip install --index-url https://pypi.tuna.tsinghua.edu.cn/simple ...`. Do not create a separate per-module environment unless it is also a documented dependency category, do not create or reuse a global environment, do not put the run environment inside the skill source directory, and do not run bare `pip install` into the system/user Python unless `uv` is unavailable or the user explicitly requests it. Keep any category-level `pyproject.toml`, `uv.lock`, or requirements export inside `uv_envs/<category>/`, and record all environment paths, categories, package versions, and mirror fallbacks under the total analysis directory, usually `tables/package_versions.txt`. Prefer TUNA mirrors for package downloads: use the TUNA PyPI mirror for `uv pip install`, and use TUNA CRAN/Bioconductor mirrors for R packages when practical. If the TUNA mirror is unavailable, stale, or missing a required package, fall back to the official source only for the affected dependency and record the mirror fallback in `tables/package_versions.txt`. If `uv` is unavailable, install `uv` first when possible; otherwise document the fallback package manager in `tables/package_versions.txt`. R packages, Cell Ranger, velocyto, CUDA drivers, and system libraries are outside `uv` and should be installed with their appropriate manager while still recording versions.

Do not create a symlink in the current directory that points to an older output just to make it look renamed. If an output needs a new name, either regenerate it in the correct output directory or make a real copied file with documented provenance.

Python random seed rule: every Python script and notebook must define and use the default fixed random seed near the top of the file: `SEED = 42`, unless the user explicitly specifies another seed. Set `random.seed(SEED)` and `numpy.random.seed(SEED)` when those libraries are used, and set framework-specific seeds for stochastic packages when relevant, such as PyTorch, TensorFlow, scvi-tools, scVelo, veloVI, CellRank, scikit-learn, Scanpy, or UMAP. Pass `random_state=SEED`, `seed=SEED`, or the package-equivalent argument to every function that supports it, including PCA, neighbors/UMAP, Leiden/clustering, train/test splits, model fitting, bootstrapping, permutations, and plotting layouts when applicable. If a function has no seed argument or still has nondeterministic GPU behavior, document that limitation. Record the fixed seed in the module parameter log, per-candidate `harmony_params.txt` or equivalent, and `tables/package_versions.txt`; do not change the seed between candidate groups unless the user explicitly requests seed sensitivity testing.

Expression-layer rule: for expression-based plotting, marker visualization, dotplots, violin plots, gene scoring, or other Scanpy-compatible calculations, if the function exposes a `use_raw` argument and `adata.raw` is present, set `use_raw=True` by default unless the user explicitly requests a different layer or the method requires counts. Record any exception and the expression layer actually used.

Cell-count reduction rule: do not downsample or subsample cells for analysis, plotting, parameter tuning, or resource-control fallback unless the user explicitly requests or approves that specific reduced-cell run. If the full data cannot run, document the blocker and stop to ask instead of producing a reduced-cell result automatically. Any approved reduced-cell output must be labeled as reduced-cell/exploratory and must not be presented as the full-data result.

Default figure formats are PDF and/or SVG only. Do not create, save, convert, or request PNG files for final, intermediate, diagnostic, preview, thumbnail, or temporary figure outputs. If a tool defaults to PNG, override it to PDF/SVG or stop and ask; do not leave `.png` files under `figures/`, `tables/`, `codes/`, or `h5ad/`.

Global figure style for every module:
- Use the module-specific `Module Figure Style Contract` in this SKILL.md as the plotting-style source of truth. Do not rely on any material outside this SKILL.md to infer final figure style.
- Plot routing is code-dispatch based. Identify the figure family first, then call the corresponding canonical plotting code or package plotting route specified in this skill, analogous to using the correct `sc.pl.*` function for each Scanpy plot. Do not use a generic Matplotlib, Seaborn, ggplot2, or ad hoc plotting pattern for a figure family that has a dedicated canonical route. If no matching route exists in this skill for a required final figure, stop and ask for the skill to be updated before running; do not improvise a new plotting route.
- For Module 04 specifically, the canonical runnable code must be split into two numbered submodules: `01-cm-lineage-analysis` produces the CM-lineage analysis tables, and `02-cm-lineage-final-plotting` reads those tables to make final figures. The same canonical route should still produce the frequency matrices, balanced joint NMF outputs, CM classification, CM node/edge tables, Epi-CM association matrices, and the final CM nodeplots, heatmaps, and scatter plots. Keep final plotting logic inside the canonical submodule 02 instead of maintaining parallel implementations for the same figure.
- Keep exactly one current canonical plotting route for each final Module 04 figure, and label any other route as non-canonical or exploratory; do not duplicate final PDFs/SVGs for the same CM-lineage panel or file stem.
- Apply the style for the corresponding analysis type and module-specific section first; use these global rules only as baseline guardrails. Do not copy a figure style from an unrelated analysis type just because the file format or package is similar.
- Save figure outputs as PDF and/or SVG only; never create PNG previews, thumbnails, diagnostics, or temporary figure files.
- Keep text editable whenever the backend supports it. In Matplotlib set `pdf.fonttype = 42`, `ps.fonttype = 42`, `svg.fonttype = "none"`, and use a standard sans-serif font such as DejaVu Sans or Arial.
- Keep axis ticks visible on quantitative and categorical plots, including heatmaps, dotplots, barplots, forest plots, scatter plots, survival plots, and UMAP panels with axes. Do not call `axis("off")`, remove tick labels, or hide spines unless the plot type is a pure network/chord/graph layout where axes have no coordinate meaning.
- Use clean white backgrounds, black axis text, readable tick labels, and legends/colorbars outside or to the right when practical. Rotate dense x-axis labels rather than letting them overlap.
- Prefer the official plotting interface for the relevant package before manual low-level drawing. Use manual Matplotlib/ggplot2/Seaborn layout code only when the package interface cannot express the required final figure or when this skill gives explicit canonical manual code; record that reason in the code or parameter log.
- Size every figure element for the exported canvas: labels, legends, colorbars, risk tables, titles, p/q labels, arrows, node labels, and panel titles must be readable and must not collide. If any element overlaps or is clipped in the generated PDF/SVG, increase canvas size, margins, row/column spacing, legend placement, font size, label wrapping, or panel spacing and rerender before considering the figure complete.
- For UMAP-like embedding panels, keep panels square, use the official Scanpy/scVelo plotting interface when possible, keep framed axes, and use consistent palettes for the same labels across full-atlas and lineage-specific plots.
- For heatmaps, use a colorbar with visible ticks. Use centered diverging palettes for signed correlations/effects, sequential palettes for nonnegative scores or transition weights, and method-specific labels in titles/captions.
- For heatmaps with significance, annotate cells from corrected significance values only. Prefer a matching q-value / adjusted-p-value matrix or column. If only raw p values are available, first compute BH-FDR q values within the relevant analysis family and context, save those q values, and annotate from the q values. Do not annotate heatmaps from uncorrected raw p values. Cell text must be significance symbols only, not numeric matrix values. Use `ns` for non-significant tested cells, `*` for `0.01 <= q < 0.05`, `**` for `0.001 <= q < 0.01`, and `***` for `q < 0.001`. Missing/untested cells may be blank, but non-significant tested cells must be `ns`.
- For dotplots, keep gene and group orders explicit, preserve tick labels, and use standard scaling only when the method requires display scaling.
- For forest/survival plots, keep confidence intervals, hazard-ratio/reference lines, p/q labels, and axis ticks visible.
- For network/chord/directed-graph plots, use PDF/SVG, editable labels, clear legends/colorbars when edge weights are encoded, and document any intentional axis removal.


Scanpy plotting hard rule: For ordinary `sc.pl.*` outputs that support `save=...`, do not create Matplotlib axes, do not pass `ax=...`, and do not save with `fig.savefig` or `ax.figure.savefig`. Use `sc.settings.figdir` plus the Scanpy `save` argument. For multi-panel Scanpy plots, use the package interface such as `color=[...]`, `ncols`, `wspace`, or `standard_scale` instead of manual subplots. Manual `ax=...` is allowed only for documented special overlays or incompatible per-panel settings that Scanpy cannot express; record the reason in the run-parameter table or a code comment.


For both Python and R plotting, use the official plotting interface of the relevant package by default, such as Scanpy, Matplotlib, Seaborn, CellRank, scVelo, ggplot2, ComplexHeatmap, or CellChat plotting APIs. For Python-generated single-panel plots, use a default square canvas of 2.5 x 2.5 inches unless the user specifies another size or the plot type clearly requires more space, such as multi-panel layouts, heatmaps, wide dotplots, survival plots, or network/chord diagrams. For UMAP plots, use the default Scanpy save path and keep the call minimal. Set figure parameters with `sc.set_figure_params(figsize=(3, 3), dpi=150)` or `sc.settings.set_figure_params(figsize=(3, 3), dpi=150)`, set `sc.settings.figdir` to the target output directory, then call `sc.pl.umap(adata, color="leiden_coarse", save="_name.pdf")`. Keep the default Scanpy-style framed axes and outside legends; do not manually create `plt.subplots`, pass `ax=...`, or call `fig.savefig` for ordinary UMAPs. Use `ncols` only when `color` contains multiple objects, for example `sc.pl.umap(adata, color=["cell_subtype", "status", "cnv_score"], ncols=3, wspace=0.4, save="_celltype_status_cnvscore.pdf")`. Do not use `ncols` for a single-color UMAP. Keep `save` as a suffix/name handled by Scanpy rather than a full path. Manual axes and `fig.savefig` are reserved for special cases where different panels require incompatible per-panel parameters or post-processing that the official `color=[...]` interface cannot express, and the reason must be documented because manual saving can greatly increase PDF/SVG size. For ordinary UMAPs, specifically avoid `return_fig=True` followed by `ax.savefig(..., bbox_inches="tight")`; on large atlases this can inflate files from sub-MB Scanpy-saved outputs to multi-MB PDFs or very large SVGs. When the task is only to inspect, audit, or explain existing UMAP code/output size, do not modify the source code or rerun plotting unless the user explicitly asks for a fix or rerender. Keep all figure text editable as text whenever the plotting backend supports it; do not convert labels, legends, tick labels, titles, or annotations to outlines/paths unless the user explicitly requests it. For multi-panel figures, especially multi-panel UMAPs, verify that each UMAP panel remains square, legends/colorbars/titles do not overlap, and adjacent panels do not collide. It is acceptable to adjust the default figure width, height, `wspace`, `hspace`, legend font size, or margins to prevent overlap and preserve square UMAP panels. Do not add extra custom titles to UMAP panels; use the Scanpy default title derived from `color` unless the user explicitly asks for custom titles. Do not add automatic bitmap/raster conversion rules; let the user decide figure-size tradeoffs from the actual output files. Do not draw sample-colored UMAPs as default final figures; generate sample-colored UMAPs only when the user asks for them or when they are needed as integration/batch-mixing diagnostics, and label them as diagnostic outputs.

## Module Figure Style Contract

Use the following CM-lineage figure styles unless the user explicitly asks for a
different style. Do not mention implementation provenance from prior runs in generated reusable code, figure labels, captions, or readme files.

- CM nodeplots: draw compact network-style panels with stable node positions
  across conditions, edge colors or widths tied to signed/weighted association
  strength, and readable node labels.
  Save the node/edge tables used for each figure.
- CM activity and loading heatmaps: use clustered heatmap-style layouts when
  clustering is part of the analysis; otherwise preserve the supplied order.
  Never cluster the CM axis. Keep CM columns/rows in the canonical CM order
  from `joint_module_classification.csv`, `activity_df_sample_by_CM.csv`, or
  the matching source matrix. If clustering is used, it may apply only to the
  non-CM axis, such as samples or cell subtypes.
  For sample x CM activity heatmaps, hide sample labels on the heatmap axis and
  use exactly two sample annotation tracks by default: `Series` and `Status`.
  Do not replace these with sample-name tick labels or status-only annotations.
  Use clear row/column labels, visible colorbar ticks, and sequential palettes
  for nonnegative activity/loadings.
- Node-node correlation heatmaps: use a centered diverging palette on a fixed
  -1 to 1 scale unless the user explicitly requests another range. Keep row and
  column order identical when the matrix is square and symmetric.
- Epi-CM association heatmaps: use centered diverging colors for signed
  correlations/effects and overlay significance symbols from the matching
  q-value or adjusted-p-value table. If only raw p values exist, compute and
  save BH-FDR q values across all epithelial subtype x CM pairs within that
  method and status context before plotting. Never use raw p values directly for
  heatmap stars. Use `ns` for non-significant tested cells, `*` for
  `0.01 <= q < 0.05`, `**` for `0.001 <= q < 0.01`, and `***` for
  `q < 0.001`. Do not annotate heatmap cells with the current data value;
  the cell color and colorbar show the value, and the overlaid text shows only
  significance. Do not display non-significant tested entries as `NA`.
- When the optional Pearson branch is requested, Spearman and Pearson
  association figures must look parallel but remain
  method-labeled: Spearman uses `rho` labels, Pearson uses `r` labels, and each
  method reads only its own matrices.
- Scatter plots for CM-lineage associations: for Epi-CM association branches,
  generate one plot per epithelial subtype x CM pair, or an explicitly arranged
  multi-panel grid that still contains every pair. Show the fitted trend line
  only when that is part of the statistic, label the statistic and corrected
  q/adjusted-p value for the selected method, and avoid mixing Spearman and
  Pearson annotations. Scatter source CSVs must keep both raw `p_value` and
  corrected `q_value`; the figure annotation must use `q`, not raw `p`.

## Code Policy

```text
Fragile analysis and final-figure workflows must have explicit code in this
SKILL.md. Generated notebooks/scripts should be created from those code blocks
in the current task output directory when needed. Mark user/project-adjustable
values with REPLACEABLE comments and non-negotiable behavior with FIXED
comments. Do not rely only on prose such as "follow the old notebook" or
"use the same script" for canonical behavior.
```

The `CmEpi` package is intentionally not bundled in this module.

Bundled CM-Epi source-code resources:

```text
references/cm_epi_analysis_code_inventory.md
scripts/cm_epi_analysis_source/
```

Use the inventory to locate bundled code extracts from the CM-Epi analysis tree
when implementing or auditing Module 04. Notebook outputs and data objects are
not bundled; notebook code cells are extracted as Python scripts. Treat
`scripts/cm_epi_analysis_source/01-canonical-balanced-joint-nmf/` as the current
balanced joint CM workflow reference and
`scripts/cm_epi_analysis_source/02-final-redraw-and-plotting/` as the final
plotting/redraw reference. Treat directories beginning with `90-` or `91-` or
`92-` as legacy references only; do not use them as canonical unless the user
explicitly asks for a legacy/reference-guided branch.

## Canonical Workflow Specification

Module 04 must be split into two explicit submodules:

```text
01-cm-lineage-analysis
02-cm-lineage-final-plotting
```

Submodule 01 is the analysis submodule. It generates only canonical tables,
reports, and figure-source tables. Submodule 02 is the plotting submodule. It
reads only the canonical outputs from submodule 01 and generates final PDF/SVG figures
plus the figure manifest. Do not interleave final plotting code into submodule 01, and do
not recompute analysis results inside submodule 02. Each final figure must have exactly
one plotting block in submodule 02. Do not create multiple plotting blocks, notebooks, or
scripts that produce the same final CM-lineage figure.

Recommended canonical code organization:

```text
01-cm-lineage-analysis:
  01_prepare_inputs_and_frequency_tables
  02_balanced_joint_nmf
  03_cm_classification_nodes_and_edges
  04_epi_cm_association_spearman
  05_epi_cm_association_pearson_optional

02-cm-lineage-final-plotting:
  06_final_unified_plotting_and_manifest
```

The canonical workflow can be one notebook, one script, or a small ordered set
of scripts, but it must preserve the submodule boundary. Prefer the subtask
structure above: subtasks 01-05 belong to `01-cm-lineage-analysis` and produce
canonical analysis tables only; subtask 06 belongs to
`02-cm-lineage-final-plotting` and is the single unified final plotting step.
Do not draw final figures in subtasks 01-05. Subtask 06 reads the canonical
tables produced by subtasks 01-05 and writes all final Module 04 PDF/SVG figures
plus the figure manifest.

Within the Module 04 category-first layout, use the submodule name before the
task name. For example:

```text
04-<project_slug>-cm-lineage-core/
  codes/01-cm-lineage-analysis/01_prepare_inputs_and_frequency_tables.py
  codes/01-cm-lineage-analysis/02_balanced_joint_nmf.py
  tables/01-cm-lineage-analysis/01_prepare_inputs_and_frequency_tables/
  tables/01-cm-lineage-analysis/02_balanced_joint_nmf/
  codes/02-cm-lineage-final-plotting/06_final_unified_plotting_and_manifest.py
  figures/02-cm-lineage-final-plotting/
  tables/02-cm-lineage-final-plotting/figure_manifest.csv
```

Subtask responsibilities:

```text
01_prepare_inputs_and_frequency_tables:
  inputs: annotated obs or supplied frequency matrices
  tables: sample_status, epithelial counts/frequencies, non-epithelial counts/frequencies
  figures: none by default; optional QC summaries only when explicitly requested

02_balanced_joint_nmf:
  inputs: non-epithelial frequency matrix, sample_status
  tables: K-selection metrics, selected K, internal raw NMF W/H only if needed for provenance
  note: raw component names such as CM1..CMK are internal and must not be consumed by downstream final plotting or association
  figures: none

03_cm_classification_nodes_and_edges:
  inputs: W, H, activity, loading, sample_status
  tables: CM classification, raw_to_canonical_CM_mapping, canonical W/H/activity/loading tables, CM node tables, status-specific edge tables, node-node correlation matrices
  figures: none

04_epi_cm_association_spearman:
  inputs: epithelial frequency matrix, CM activity matrix, sample_status
  tables: all epithelial subtype x CM Spearman rho, P value, q value, significance labels, and scatter source tables
  figures: none

05_epi_cm_association_pearson_optional:
  inputs: epithelial frequency matrix, CM activity matrix, sample_status
  tables: all epithelial subtype x CM Pearson r, P value, q value, significance labels, and scatter source tables
  figures: none

06_final_unified_plotting_and_manifest:
  inputs: outputs from subtasks 01-05
  tables: final output manifest, package versions, parameter report, validation report
  figures: all final Module 04 PDF/SVG figures using the final plotting code
```

Canonical analysis flow:

```text
load annotated obs or supplied frequency matrices
validate sample IDs, status, subtype names, and matrix orientation
derive broad cell_type from cell_subtype prefix when needed
split epithelial and non-epithelial compartments
build sample x subtype count and frequency matrices
align sample_status, epithelial frequency, and non-epithelial frequency tables
apply sample eligibility filters before CM construction
apply column-wise min-max normalization to the non-epithelial frequency matrix
apply status-balanced sample weights for NMF when both status groups are present
run K-selection over the configured candidate K and seed grid
select final K from the K-selection report unless the user explicitly forces K
fit final NMF with the selected K on the weighted column-minmax non-epithelial matrix
refit sample CM activities by NNLS on the unweighted column-minmax matrix using the fixed H basis
classify modules as sharedCM, normalCM, or tCM
rename classified modules with the canonical s/n/t_CM<number> convention
write canonical W/H/activity/loading tables immediately after classification
derive CM node tables from top H loadings
derive status-specific CM edges and node-node correlation tables
run Epi-CM Spearman association from canonical CM activity as the default branch
run Epi-CM Pearson association only when requested, in a separate branch
write tables, figure inputs, parameter reports, package versions, and manifest
```

Canonical table outputs:

```text
non_epi_subtype_counts.csv
non_epi_subtype_frequency.csv
non_epi_subtype_frequency_column_minmax.csv
column_minmax_params.csv
group_balanced_sample_weights.csv
epi_subtype_counts.csv
epi_subtype_frequency.csv
sample_status.csv
sample_compartment_cell_counts.csv
sample_inclusion_exclusion.csv
joint_nmf_k_selection_metrics.csv
selected_module_k.json
raw_to_canonical_CM_mapping.csv
W_df.csv
H_df.csv
activity_df_sample_by_CM.csv
loading_df_cell_subtype_by_CM.csv
activity_df_CM_by_sample.csv
w_df_activity_sample_by_CM_raw.csv
w_df_activity_sample_by_CM_zscore.csv
w_df_activity_sample_by_CM_robust.csv
w_df_activity_sample_by_CM_standard_scale_col.csv
w_df_activity_sample_by_CM_standard_scale_col_with_series_plot_matrix.csv
w_df_activity_sample_by_CM_standard_scale_col_with_series_row_annotations.csv
h_df_loading_cell_subtype_by_CM_raw.csv
h_df_loading_cell_subtype_by_CM_zscore.csv
h_df_loading_cell_subtype_by_CM_robust.csv
h_df_loading_cell_subtype_by_CM_standard_scale_col.csv
loading_df_cell_subtype_by_CM_fraction.csv
balanced_joint_cm_subtype_loadings_raw_from_H_df.csv
balanced_joint_cm_subtype_loadings_fraction_from_H_df.csv
joint_CM_activity_tumor_vs_normal_mean_sd_summary.csv
joint_module_classification.csv
status_specific_nodeplot_edges.csv
balanced_joint_cm_reference_node_sets_after_edge_threshold.csv
tumor_network_nodes_from_H_df.csv
normal_like_network_nodes_from_H_df.csv
node_node_correlation_matrix.csv
node_node_correlation_q_matrix.csv
balanced_joint_cm_epi_cm_association_tumor_rho_matrix.csv
balanced_joint_cm_epi_cm_association_tumor_q_matrix.csv
balanced_joint_cm_epi_cm_association_normal-like_rho_matrix.csv
balanced_joint_cm_epi_cm_association_normal-like_q_matrix.csv
```

When Pearson is requested, add Pearson-specific tables with `pearson_r` and
`pearson_q` naming. Do not rename Spearman `rho` outputs into Pearson outputs.

Canonical plotting flow:

```text
run plotting only in 02-cm-lineage-final-plotting / 06_final_unified_plotting_and_manifest
read canonical figure input tables from subtasks 01-05
validate that required submodule 01 input tables exist before plotting
implement plotting with the final plotting code for each figure family
set PDF/SVG font options before creating figures
draw K-selection diagnostics from joint_nmf_k_selection_metrics.csv
draw CM activity heatmaps from the raw, z-score, robust, and column-min-max W/activity tables
draw CM loading heatmaps from the raw, z-score, robust, and column-min-max H/loading tables
draw CM activity tumor-vs-normal-like mean +/- SD barplot from joint_CM_activity_tumor_vs_normal_mean_sd_summary.csv
draw CM nodeplots from the canonical node tables and status_specific_nodeplot_edges.csv using the fixed nodeplot structure below
draw node-node correlation heatmaps from node-node correlation matrices
draw Epi-CM Spearman heatmaps from Spearman rho/q matrices only
draw all Epi-CM Spearman scatter plots from the same Spearman inputs
draw Epi-CM Pearson heatmaps/scatter plots only when Pearson tables exist, using all pairwise combinations
write one figure manifest listing every final PDF/SVG, direct input table, and plotting block
```

Canonical final plotting implementation:

```text
06_final_unified_plotting_and_manifest must contain the only final plotting code
for Module 04 and must live in submodule `02-cm-lineage-final-plotting`. It may
define shared helpers such as set_figure_style(), save_pdf_svg(), read_matrix(),
star_from_q(), and aligned_order(), then call one plotting function per final
figure family. Do not create separate export or figure-fix scripts for the same
final panel.
```

When adapting an existing project workflow into the canonical Module 04 code,
merge the final plotting logic into `06_final_unified_plotting_and_manifest`.
If a figure family has a later `plot*.ipynb`/`plot*.py` redraw in the project,
that redraw style and output stem replaces the original analysis-notebook
plotting block for the same figure. If no later redraw exists for a figure
family, keep the original analysis plotting block and add it to subtask 06.
Never keep both the original and redrawn plotting implementations as canonical
routes for the same final figure.

Use these concrete plotting rules in 06:

- Set Matplotlib/PDF/SVG font options before any figure is created:
  `pdf.fonttype = 42`, `ps.fonttype = 42`, `svg.fonttype = "none"`, white
  background, black tick/axis text, and editable labels.
- Save every final figure through one `save_pdf_svg(fig, stem)` helper that
  writes exactly `<stem>.pdf` and `<stem>.svg`; do not write PNG, previews, or
  duplicate alternative stems.
- K-selection diagnostics: plot candidate K on the x-axis with visible ticks;
  show selection score, explained fraction, reconstruction error, and stability
  in a compact 1 x 3 multi-panel diagnostic when those are the available final
  plotted metrics. Mark the selected K with a vertical dashed line or highlighted
  point. Source table:
  `joint_nmf_k_selection_metrics.csv`.
- Heatmap display scale definitions are fixed and must be recorded in the
  figure manifest: `raw` means the raw W activity or H loading value without
  display normalization; `standard_scale_col` means each CM column is
  independently min-max scaled to `[0, 1]`; `zscore` means each CM column is
  centered by its own mean and divided by its own standard deviation with
  `ddof=0`; `robust` means each CM column is centered by its own median and
  divided by its own IQR. These are display scales only; do not use
  display-scaled values for statistical tests or association calculations.
- CM activity heatmaps: rows are samples and columns are CM IDs. Preserve the
  supplied CM order and never cluster CM columns. Preserve sample order unless
  sample clustering is explicitly part of the analysis.
  Generate the same normalization/display variants as the analysis tables:
  raw, z-score, robust-scaled, and column min-max
  (`standard_scale_col`) W/activity heatmaps. Use sequential `viridis` palettes
  for raw/nonnegative activity and column min-max displays, and centered
  diverging palettes for signed z-score/robust displays. Keep visible colorbar
  ticks, readable annotations, and no hidden axes except that sample labels must
  be hidden on sample x CM activity heatmaps. Source
  tables include `activity_df_sample_by_CM.csv`,
  `w_df_activity_sample_by_CM_raw.csv`,
  `w_df_activity_sample_by_CM_zscore.csv`,
  `w_df_activity_sample_by_CM_robust.csv`, and
  `w_df_activity_sample_by_CM_standard_scale_col.csv`.
  The unannotated
  `w_df_activity_sample_activity_per_CM_standard_scale_col_clustermap` and the
  annotated
  `w_df_activity_sample_activity_per_CM_standard_scale_col_with_series_clustermap`
  must both be drawn when series/status annotations are available. The annotated
  column min-max activity heatmap must use the clustered style with exactly
  `Series` and `Status` row annotation tracks and write
  `w_df_activity_sample_activity_per_CM_standard_scale_col_with_series_clustermap`.
  The corresponding matrix and row annotation tables are
  `w_df_activity_sample_by_CM_standard_scale_col_with_series_plot_matrix.csv`
  and
  `w_df_activity_sample_by_CM_standard_scale_col_with_series_row_annotations.csv`.
- CM loading heatmaps: final plotting input must be converted to
  `cell_subtype x canonical CM` before plotting. Generate the same
  normalization/display variants as the analysis tables: raw, z-score,
  robust-scaled, and column min-max (`standard_scale_col`) heatmaps. Use
  sequential `viridis` palettes for raw/nonnegative loading and column min-max
  displays, and centered diverging palettes for signed z-score/robust displays.
  Never cluster the CM axis; keep CM IDs in canonical order even when subtype
  rows are clustered. Keep visible colorbar ticks and readable subtype labels. Source tables include
  `loading_df_cell_subtype_by_CM.csv`,
  `h_df_loading_cell_subtype_by_CM_raw.csv`,
  `h_df_loading_cell_subtype_by_CM_zscore.csv`,
  `h_df_loading_cell_subtype_by_CM_robust.csv`,
  and `h_df_loading_cell_subtype_by_CM_standard_scale_col.csv`.
  Despite their historical names, `H_df.csv` and
  `h_df_loading_cell_subtype_by_CM_{raw,zscore,robust,standard_scale_col}.csv`
  are `CM x cell_subtype` tables unless the matrix-orientation readme explicitly
  says otherwise; transpose them before plotting loading heatmaps. The
  `loading_df_cell_subtype_by_CM*.csv` tables are the canonical
  `cell_subtype x CM` orientation. Always validate and log the detected
  orientation instead of assuming from filename alone.
- CM activity mean +/- SD barplot: draw grouped bars for normal-like and tumor
  mean activity per CM with SD error bars, using
  `joint_CM_activity_tumor_vs_normal_mean_sd_summary.csv`. This figure is a
  required final output and must be written exactly as
  `activity_df_tumor_vs_normal_mean_sd_barplot.pdf` and
  `activity_df_tumor_vs_normal_mean_sd_barplot.svg`. Keep CM order from the
  table, color bars by status, keep y-axis ticks visible, rotate dense CM labels
  if needed, and mark/annotate MWU/BH q significance only from the summary
  table.
- CM nodeplots: use deterministic node coordinates from the final code
  (`SEED = 42`) and reuse the same coordinates across status panels so condition
  differences are visually comparable. Node size or label weight reflects
  loading rank or loading magnitude; edge
  width reflects edge weight; edge color distinguishes tumor-only,
  normal-like-only, and shared edges. Network/chord axes may be hidden only
  because graph layout coordinates have no axis meaning. Source tables:
  `tumor_network_nodes_from_H_df.csv`,
  `normal_like_network_nodes_from_H_df.csv`, and
  `status_specific_nodeplot_edges.csv`.
- Node-node correlation heatmap: use the top10 H-loading diagnostic nodes for
  each CM, regardless of whether those nodes survive edge filtering or appear in
  the nodeplot. Keep row and column order identical for square matrices, use a
  centered diverging palette fixed at -1 to 1, show colorbar ticks, and do not
  annotate cells with numeric correlation values. By default leave heatmap cells
  text-free; if the user explicitly requests significance labels, overlay
  q-value symbols only from the matching q table. Source tables:
  `node_node_correlation_matrix.csv` and
  `node_node_correlation_q_matrix.csv`.
- Spearman Epi-CM heatmaps: read only Spearman rho and Spearman q matrices. Use
  a centered diverging palette fixed at -1 to 1 unless the user explicitly
  requests another range. Overlay q-value symbols from the Spearman q matrix and
  label the statistic as `rho`. The canonical final heatmap stems use the
  `_heatmap_qstars` suffix for the separate tumor and normal-like association
  outputs only. Do not read Pearson tables in this plotting block.
- Pearson Epi-CM heatmaps: create only when Pearson tables exist or the user
  requests Pearson. Read only Pearson r and Pearson q matrices. Use parallel
  layout and colors to the Spearman figures, but label the statistic as `r`.
- Epi-CM scatter plots: read the same method-specific inputs used by the
  heatmap and generate the full epithelial subtype x CM Cartesian product for
  each eligible status context. Do not restrict scatter output to selected,
  representative, significant, top-ranked, or manuscript-highlighted pairs
  unless the user explicitly asks for an additional display subset after the
  full scatter set exists. Spearman scatter plots and Pearson scatter plots must
  be rendered in separate method-specific figure directories and must read only
  the matching method-specific correlation summary/table. Use one panel per pair
  by default; a planned multi-panel grid is allowed only if every pair in the
  full Cartesian product is still present. Keep ticks visible, plot points with
  moderate size and alpha, add a fitted trend line only when the method/config
  says to show it, and annotate the correct statistic and corrected q/adjusted-p
  value for the selected method. Do not annotate scatter plots with raw p values
  when q values are available or can be computed from the full pair set. Color scatter points by the epithelial subtype color used in
  `adata_epi`/the lineage-specific `cell_subtype` palette; if the palette is
  supplied as a CSV, read that exact mapping and validate that every epithelial
  subtype has a color. For the default Spearman branch, use the Spearman-only
	  scatter summary and labels: `Spearman rho=...`, `q=...`, `n=...`, y-axis fixed
	  to 0-1, one `sns.regplot` per pair, `sns.set_style("ticks")`, canvas 6 x 5
	  inches, point size 60, alpha 0.85, no marker edge, line width 2, and a white
	  rounded statistic box at the upper-left. The regression line must include the
	  default seaborn confidence band/shadow (`ci=95`); do not disable it with
	  `ci=None`. This is the canonical style for the
	  former `W_df_epi_frequency_scatter_plots_spearman_only` output, adapted to
	  PDF/SVG only. Do not compute or annotate Pearson in the same scatter output.
  Per-pair scatter stems should follow `scatter_<CM>_vs_<epi_subtype>` and write
  PDF/SVG only.

Use this explicit Epi-CM scatterplot implementation pattern in the final
plotting code. Items marked `REPLACEABLE` may be changed by user request or
project-specific style. Items marked `FIXED` define required behavior.

```python
from pathlib import Path
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr, spearmanr
from statsmodels.stats.multitest import multipletests

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"
sns.set_style("ticks")


def save_pdf_svg(fig, stem: Path):
    stem.parent.mkdir(parents=True, exist_ok=True)
    pdf = stem.with_suffix(".pdf")
    svg = stem.with_suffix(".svg")
    fig.savefig(pdf, bbox_inches="tight", dpi=300)
    fig.savefig(svg, bbox_inches="tight", dpi=300)
    plt.close(fig)
    return pdf, svg


def safe_corr(x: pd.Series, y: pd.Series, method: str):
    ok = x.notna() & y.notna()
    if ok.sum() < 3 or x[ok].nunique() < 2 or y[ok].nunique() < 2:
        return np.nan, np.nan
    if method == "spearman":
        return spearmanr(x[ok], y[ok])
    if method == "pearson":
        return pearsonr(x[ok], y[ok])
    raise ValueError(f"Unsupported method: {method}")


def safe_filename_part(x):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(x)).strip("_")


def fmt_num(x, digits=3):
    if pd.isna(x):
        return "NA"
    return f"{x:.{digits}f}"


def fmt_p(x):
    if pd.isna(x):
        return "NA"
    return f"{x:.1e}" if x < 0.001 else f"{x:.3f}"


def q_values_for_summary(summary: pd.DataFrame) -> pd.DataFrame:
    # FIXED: FDR is computed across all epithelial subtype x CM pairs within
    # each method and status context, not only plotted/highlighted pairs.
    out = summary.copy()
    out["q_value"] = np.nan
    for (method, status), idx in out.groupby(["method", "status"], sort=False).groups.items():
        p = out.loc[idx, "p_value"].to_numpy(float)
        valid = np.isfinite(p)
        q = np.full(p.shape, np.nan, dtype=float)
        if valid.any():
            q[valid] = multipletests(p[valid], method="fdr_bh")[1]
        out.loc[idx, "q_value"] = q
    return out


def load_cell_subtype_colors(color_csv: Path, required_subtypes: list[str]) -> dict[str, str]:
    # REPLACEABLE input path. FIXED behavior: validate full palette coverage.
    color_df = pd.read_csv(color_csv)
    if not {"cell_subtype", "color"}.issubset(color_df.columns):
        raise ValueError("Color CSV must contain columns: cell_subtype, color")
    colors = dict(zip(color_df["cell_subtype"].astype(str), color_df["color"].astype(str)))
    missing = sorted(set(map(str, required_subtypes)) - set(colors))
    if missing:
        raise ValueError(f"Missing epithelial subtype colors: {missing}")
    return colors


def plot_epi_cm_scatter_all_pairs(
    epi_freq: pd.DataFrame,
    cm_activity: pd.DataFrame,
    sample_status: pd.DataFrame,
    color_map: dict[str, str],
    out_dir: Path,
    method: str,
    status_order=("tumor", "normal-like"),  # REPLACEABLE status labels only.
    point_size=60,  # REPLACEABLE style only; canonical scatter uses 60.
    point_alpha=0.85,  # REPLACEABLE style only; canonical scatter uses 0.85.
):
    # FIXED: generate every epithelial subtype x CM pair for every status context.
    # Do not pass a selected/top/significant pair list into this function.
    stat_label = "rho" if method == "spearman" else "r"
    test_label = "Spearman" if method == "spearman" else "Pearson"
    out_dir = Path(out_dir)
    common_samples = epi_freq.index.intersection(cm_activity.index).intersection(sample_status.index)
    E = epi_freq.loc[common_samples].astype(float)
    C = cm_activity.loc[common_samples].astype(float)
    status = sample_status.loc[common_samples, "status"].astype(str)

    summary_rows = []
    plot_payloads = []
    figure_rows = []
    for status_value in status_order:
        status_samples = status.index[status.eq(status_value)]
        if len(status_samples) == 0:
            continue
        for cm in C.columns:
            for epi in E.columns:
                plot_df = pd.DataFrame(
                    {
                        "sample_id": status_samples,
                        "CM_score": C.loc[status_samples, cm].to_numpy(float),
                        "Epi_fraction": E.loc[status_samples, epi].to_numpy(float),
                    }
                ).dropna(subset=["CM_score", "Epi_fraction"])
                r, p = safe_corr(plot_df["CM_score"], plot_df["Epi_fraction"], method)
                summary_rows.append({
                    "method": method,
                    "status": status_value,
                    "CM": cm,
                    "epi_subtype": epi,
                    stat_label: r,
                    "p_value": p,
                    "n_samples": int(plot_df.shape[0]),
                    "epi_color": color_map[str(epi)],
                })
                plot_payloads.append((status_value, cm, epi, plot_df, r, p))

    summary = pd.DataFrame(summary_rows)
    if not summary.empty:
        summary = q_values_for_summary(summary)
    q_lookup = {
        (row.method, row.status, row.CM, row.epi_subtype): row.q_value
        for row in summary.itertuples(index=False)
    }

    # FIXED: scatter figure annotations use corrected q values, not raw p values.
    for status_value, cm, epi, plot_df, r, p in plot_payloads:
        status_dir = out_dir / method / status_value.replace("-", "_")
        q = q_lookup.get((method, status_value, cm, epi), np.nan)
        fig, ax = plt.subplots(figsize=(6, 5))  # REPLACEABLE canvas; canonical scatter uses 6 x 5.
        sns.regplot(
            data=plot_df,
            x="CM_score",
            y="Epi_fraction",
            ax=ax,
	            color=color_map[str(epi)],
	            scatter_kws={"s": point_size, "alpha": point_alpha, "edgecolor": "none"},
	            line_kws={"lw": 2},
	            ci=95,
	        )
        ax.set_xlabel(f"{cm} score")
        ax.set_ylabel(f"{epi} fraction")
        sample_label = "Tumor samples" if status_value == "tumor" else "Normal-like samples"
        ax.set_title(f"{sample_label}: {cm} vs {epi}")
        ax.set_ylim(0, 1)
        ax.text(
            0.02,
            0.98,
            f"{test_label} {stat_label}={fmt_num(r)}, q={fmt_p(q)}\nn={len(plot_df)}",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=10,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8),
        )
        sns.despine(ax=ax)
        fig.tight_layout()
        stem = status_dir / f"scatter_{safe_filename_part(cm)}_vs_{safe_filename_part(epi)}"
        pdf, svg = save_pdf_svg(fig, stem)
        figure_rows.extend([
            {"method": method, "status": status_value, "CM": cm,
             "epi_subtype": epi, "figure_file": str(pdf)},
            {"method": method, "status": status_value, "CM": cm,
             "epi_subtype": epi, "figure_file": str(svg)},
        ])

    out_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out_dir / f"epi_cm_{method}_scatter_all_pairs_summary.csv", index=False)
    pd.DataFrame(figure_rows).to_csv(out_dir / f"epi_cm_{method}_scatter_all_pairs_manifest.csv", index=False)
    return summary


# Required calls in the final plotting script:
# color_map = load_cell_subtype_colors(color_csv, list(epi_subtype_frequency.columns))
# plot_epi_cm_scatter_all_pairs(epi_freq, cm_activity, sample_status, color_map, out_dir, method="spearman")
# plot_epi_cm_scatter_all_pairs(epi_freq, cm_activity, sample_status, color_map, out_dir, method="pearson")
```

Per-figure code contract:

The final plotting script must contain code equivalent to the following
functions for each figure family. Do not replace these functions with vague
rules, external notebook references, or ad hoc plotting code. Items marked
`REPLACEABLE` are paths, labels, palettes, or project-specific ordering that
the user may change. Items marked `FIXED` define the canonical behavior.

```python
from pathlib import Path
import math
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.ticker import FormatStrFormatter

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["font.family"] = "DejaVu Sans"  # REPLACEABLE font only.
mpl.rcParams["axes.linewidth"] = 0.8
sns.set_style("ticks")


def save_pdf_svg(fig, stem: Path):
    # FIXED: final figures are PDF/SVG only.
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", dpi=300)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", dpi=300)
    plt.close(fig)


def cm_sort_key(cm: str) -> tuple[int, str]:
    # FIXED: canonical CM IDs are globally numbered, e.g. s_CM1, t_CM2, n_CM3.
    m = re.search(r"_CM(\d+)$", str(cm))
    return (int(m.group(1)) if m else 10**9, str(cm))


def q_stars(q: float, ns_label: str = "ns") -> str:
    # FIXED: heatmap annotations show significance only, never numeric values.
    # ns: q >= 0.05; *: 0.01 <= q < 0.05; **: 0.001 <= q < 0.01; ***: q < 0.001.
    if not np.isfinite(q):
        return ns_label
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return ns_label


def plot_joint_nmf_k_selection(metrics_csv: Path, out_dir: Path):
    # Inputs: joint_nmf_k_selection_metrics.csv.
    # Outputs: joint_nmf_k_selection.pdf/svg.
    metrics = pd.read_csv(metrics_csv).sort_values("k")
    k = metrics["k"].astype(int).to_numpy()
    selected = metrics.loc[metrics["selected"].astype(bool), "k"]
    selected_k = int(selected.iloc[0]) if len(selected) else None

    fig, axes = plt.subplots(1, 3, figsize=(10, 3), constrained_layout=True)
    panels = [
        ("best_balanced_explained_fraction", "Balanced fit", "Explained fraction", "#1f77b4", np.arange(0.50, 0.96, 0.10)),
        ("stability_matched_cosine", "Stability", "Matched cosine", "#2ca02c", np.arange(0.85, 1.01, 0.05)),
        ("selection_score", "Selection score", "Score", "#d62728", np.arange(0.50, 0.81, 0.10)),
    ]
    for ax, (col, title, ylabel, color, yticks) in zip(axes, panels):
        ax.plot(k, metrics[col].astype(float), marker="o", markersize=3.8, linewidth=1.8, color=color)
        ax.set_title(title, fontweight="bold")
        ax.set_xlabel("k")
        ax.set_ylabel(ylabel)
        ax.set_xticks(k)
        ax.set_xlim(k.min() - 0.5, k.max() + 0.5)
        ax.set_yticks(yticks)
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
        ax.tick_params(axis="x", labelsize=7, length=3)
        ax.tick_params(axis="y", labelsize=8, length=3)
        ax.grid(axis="y", color="#d9d9d9", linewidth=0.6, alpha=0.8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if selected_k is not None:
            ax.axvline(selected_k, color="black", linestyle="--", linewidth=1)
    save_pdf_svg(fig, Path(out_dir) / "joint_nmf_k_selection")


def plot_activity_mean_sd_barplot(summary_csv: Path, out_dir: Path):
    # Inputs: joint_CM_activity_tumor_vs_normal_mean_sd_summary.csv.
    # Outputs: activity_df_tumor_vs_normal_mean_sd_barplot.pdf/svg.
    df = pd.read_csv(summary_csv)
    cms = df["CM"].astype(str).tolist()
    statuses = ["normal-like", "tumor"]  # REPLACEABLE labels only if metadata differs.
    palette = {"normal-like": "#377EB8", "tumor": "#E41A1C"}  # REPLACEABLE colors.
    x = np.arange(len(cms))
    bar_width = 0.36
    offsets = {"normal-like": -bar_width / 2, "tumor": bar_width / 2}

    fig, ax = plt.subplots(figsize=(max(5, len(cms) * 0.72), 3.4))
    all_lows, all_highs = [], []
    for status in statuses:
        mean = df[f"{status}_mean"].astype(float).to_numpy()
        sd = df[f"{status}_sd"].astype(float).to_numpy()
        all_lows.append(mean - sd)
        all_highs.append(mean + sd)
        ax.bar(x + offsets[status], mean, width=bar_width, color=palette[status],
               label=status, edgecolor="white", linewidth=0.4)
        ax.errorbar(x + offsets[status], mean, yerr=sd, fmt="none",
                    ecolor="black", elinewidth=0.8, capsize=2, capthick=0.8, zorder=3)

    lower = float(np.nanmin(np.concatenate(all_lows)))
    upper = float(np.nanmax(np.concatenate(all_highs)))
    ymin = min(0.0, math.floor(lower / 2.0) * 2.0)
    ymax = max(10.0, math.ceil(upper / 2.0) * 2.0)
    ax.set_ylim(ymin, ymax)
    ax.set_yticks(np.arange(ymin, ymax + 0.1, 2.0))
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(cms, rotation=45, ha="right", fontsize=8)
    ax.tick_params(axis="y", labelsize=8, length=3)
    ax.tick_params(axis="x", length=3)
    ax.grid(axis="y", color="#d9d9d9", linewidth=0.6, alpha=0.8)
    ax.set_xlabel("")
    ax.set_ylabel("Mean CM activity +/- SD")
    ax.set_title("Joint CM activity in tumor vs normal-like samples", fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(title="", frameon=False)
    fig.tight_layout()
    save_pdf_svg(fig, Path(out_dir) / "activity_df_tumor_vs_normal_mean_sd_barplot")


def transform_cm_columns(df: pd.DataFrame, method: str) -> pd.DataFrame:
    # FIXED: transforms are per CM column, not whole-matrix/global transforms.
    # raw: unchanged raw W activity or H loading values.
    # standard_scale_col: each CM column independently min-max scaled to [0, 1].
    # zscore: each CM column centered by mean and divided by std(ddof=0).
    # robust: each CM column centered by median and divided by IQR.
    # These are display scales only; never feed them back into statistics.
    x = df.astype(float).copy()
    if method == "raw":
        return x
    if method == "standard_scale_col":
        mn = x.min(axis=0)
        rng = (x.max(axis=0) - mn).replace(0, np.nan)
        return x.sub(mn, axis=1).div(rng, axis=1).fillna(0.0).clip(0, 1)
    if method == "zscore":
        return x.sub(x.mean(axis=0), axis=1).div(x.std(axis=0, ddof=0).replace(0, np.nan), axis=1).fillna(0.0)
    if method == "robust":
        med = x.median(axis=0)
        iqr = (x.quantile(0.75, axis=0) - x.quantile(0.25, axis=0)).replace(0, np.nan)
        return x.sub(med, axis=1).div(iqr, axis=1).fillna(0.0)
    raise ValueError(f"Unknown transform method: {method}")


def heatmap_scale(df: pd.DataFrame, method: str, value_label: str):
    if method in {"zscore", "robust"}:
        values = df.to_numpy(float)
        values = values[np.isfinite(values)]
        vmax = float(np.nanpercentile(np.abs(values), 98)) if values.size else 1.0
        if not np.isfinite(vmax) or vmax == 0:
            vmax = 1.0
        return "vlag", -vmax, vmax, 0, f"{value_label} ({method})"
    if method == "standard_scale_col":
        return "viridis", 0, 1, None, f"Column min-max {value_label.lower()}"
    return "viridis", 0, None, None, value_label


def plot_cm_clustermap(df: pd.DataFrame, method: str, stem: Path, title: str,
                       value_label: str, row_colors: pd.DataFrame | None = None,
                       hide_sample_labels: bool = False):
    # FIXED: columns are canonical CMs and must not be clustered.
    plot_df = transform_cm_columns(df, method)
    plot_df = plot_df.loc[:, sorted(plot_df.columns.astype(str), key=cm_sort_key)]
    cmap, vmin, vmax, center, cbar_label = heatmap_scale(plot_df, method, value_label)
    figsize = (plot_df.shape[1] * 0.35 + 2.5, plot_df.shape[0] * 0.13 + 2.5)
    g = sns.clustermap(
        plot_df,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        center=center,
        cbar_kws={"label": cbar_label},
        figsize=figsize,
        row_colors=row_colors,
        colors_ratio=(0.06, 0.01),
        col_cluster=False,
        linewidths=0,
    )
    g.fig.suptitle(title, y=1.02)
    for spine in g.ax_heatmap.spines.values():
        spine.set_visible(False)
    if hide_sample_labels:
        # FIXED for sample x CM activity heatmaps: do not show sample labels.
        g.ax_heatmap.set_yticklabels([])
        g.ax_heatmap.set_yticks([])
        g.ax_heatmap.tick_params(left=False)
        g.ax_heatmap.set_ylabel("")
    g.fig.tight_layout(rect=[0, 0, 1, 0.98])
    save_pdf_svg(g.fig, stem)
    return plot_df


def plot_activity_heatmap_variants(activity_sample_by_cm: pd.DataFrame, sample_status: pd.DataFrame,
                                   series_by_sample: pd.Series, out_dir: Path, table_dir: Path):
    # Input orientation: sample x CM. Outputs four W/activity heatmap variants.
    status_palette = {"normal-like": "#377EB8", "tumor": "#E41A1C"}  # REPLACEABLE.
    series_values = sorted(series_by_sample.dropna().astype(str).unique())
    series_color_list = [
        "#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02",
        "#A6761D", "#666666", "#E7298A", "#00897B", "#7B3294",
    ]
    if len(series_values) > len(series_color_list):
        raise ValueError(f"Need more series colors: {len(series_values)} series values")
    series_palette = {series: series_color_list[i] for i, series in enumerate(series_values)}

    common = sample_status.index.intersection(activity_sample_by_cm.index).intersection(series_by_sample.index)
    df = activity_sample_by_cm.loc[common].copy()
    if df.empty:
        raise ValueError("No overlapping samples for sample x CM activity heatmap.")
    row_colors = pd.DataFrame(
        {
            # FIXED annotation columns for sample x CM heatmaps: Series and Status.
            "Series": series_by_sample.loc[common].astype(str).map(series_palette).fillna("#999999"),
            "Status": sample_status.loc[common, "status"].astype(str).map(status_palette).fillna("#999999"),
        },
        index=common,
    )
    titles = {
        "raw": "W_df activity: sample activity per CM (raw)",
        "standard_scale_col": "W_df activity: sample activity per CM (column min-max)",
        "zscore": "W_df activity: sample activity per CM (z-score)",
        "robust": "W_df activity: sample activity per CM (robust z-score)",
    }
    for method, title in titles.items():
        stem = Path(out_dir) / f"w_df_activity_sample_activity_per_CM_{method}_clustermap"
        plot_df = plot_cm_clustermap(
            df, method, stem, title, "Activity",
            row_colors=row_colors, hide_sample_labels=True,
        )
        plot_df.to_csv(Path(table_dir) / f"w_df_activity_sample_by_CM_{method}.csv")


def plot_loading_heatmap_variants(loading_cell_subtype_by_cm: pd.DataFrame, out_dir: Path, table_dir: Path):
    # Input orientation: cell_subtype x CM. Outputs four H/loading heatmap variants.
    titles = {
        "raw": "H_df loading: cell subtype weights per CM (raw)",
        "standard_scale_col": "H_df loading: cell subtype weights per CM (column min-max)",
        "zscore": "H_df loading: cell subtype weights per CM (z-score)",
        "robust": "H_df loading: cell subtype weights per CM (robust z-score)",
    }
    for method, title in titles.items():
        stem = Path(out_dir) / f"h_df_loading_cell_subtype_weights_per_CM_{method}_clustermap"
        plot_df = plot_cm_clustermap(loading_cell_subtype_by_cm, method, stem, title, "Loading")
        plot_df.to_csv(Path(table_dir) / f"h_df_loading_cell_subtype_by_CM_{method}.csv")


def plot_activity_minmax_with_series_clustermap(
    w_minmax: pd.DataFrame,
    sample_status: pd.DataFrame,
    series_by_sample: pd.Series,
    out_dir: Path,
    table_dir: Path,
):
    # Canonical dedicated redraw for the W min-max clustermap with status and
    # series row annotations. Use this instead of the older generic W heatmap
    # when `series_by_sample` is available.
    status_palette = {"normal-like": "#377EB8", "tumor": "#E41A1C"}
    series_values = sorted(series_by_sample.dropna().astype(str).unique())
    series_color_list = [
        "#1B9E77", "#D95F02", "#7570B3", "#66A61E", "#E6AB02",
        "#A6761D", "#666666", "#E7298A", "#00897B", "#7B3294",
    ]
    if len(series_values) > len(series_color_list):
        raise ValueError(f"Need more series colors: {len(series_values)} series values")
    series_palette = {series: series_color_list[i] for i, series in enumerate(series_values)}

    common = sample_status.index.intersection(w_minmax.index).intersection(series_by_sample.index)
    plot_df = w_minmax.loc[common].copy()
    status_series = sample_status.loc[common, "status"].astype(str)
    series_series = series_by_sample.loc[common].astype(str)
    if plot_df.empty:
        raise ValueError("No overlapping samples for W min-max with series clustermap.")

    row_colors = pd.DataFrame(
        {
            "Series": series_series.map(series_palette).fillna("#999999"),
            "Status": status_series.map(status_palette).fillna("#999999"),
        },
        index=plot_df.index,
    )
    figsize = (plot_df.shape[1] * 0.35 + 2.5, plot_df.shape[0] * 0.13 + 2.5)
    g = sns.clustermap(
        plot_df,
        cmap="viridis",
        vmin=0,
        vmax=1,
        center=None,
        cbar_kws={"label": "Column min-max activity"},
        col_cluster=False,
        figsize=figsize,
        row_colors=row_colors,
        colors_ratio=(0.06, 0.01),
    )
    g.fig.suptitle("W_df activity: sample activity per CM (column min-max, with series)", y=1.02)
    ax = g.ax_heatmap
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_yticklabels([])
    ax.set_yticks([])
    ax.tick_params(left=False)
    ax.set_ylabel("")

    status_handles = [mpl.patches.Patch(facecolor=color, label=label) for label, color in status_palette.items()]
    series_handles = [mpl.patches.Patch(facecolor=color, label=label) for label, color in series_palette.items()]
    status_legend = g.ax_col_dendrogram.legend(
        handles=status_handles,
        title="Status",
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        frameon=False,
    )
    g.ax_col_dendrogram.add_artist(status_legend)
    g.ax_col_dendrogram.legend(
        handles=series_handles,
        title="Series",
        loc="upper left",
        bbox_to_anchor=(1.01, 0.55),
        frameon=False,
    )
    g.fig.tight_layout(rect=[0, 0, 1, 0.98])
    stem = Path(out_dir) / "w_df_activity_sample_activity_per_CM_standard_scale_col_with_series_clustermap"
    save_pdf_svg(g.fig, stem)
    plot_df.to_csv(Path(table_dir) / "w_df_activity_sample_by_CM_standard_scale_col_with_series_plot_matrix.csv")
    pd.DataFrame({"Series": series_series, "Status": status_series}).to_csv(
        Path(table_dir) / "w_df_activity_sample_by_CM_standard_scale_col_with_series_row_annotations.csv"
    )


def plot_joint_module_top_subtype_heatmap(H_cm_by_subtype: pd.DataFrame, top_nodes: pd.DataFrame,
                                          out_dir: Path, top_n: int = 12):
    # Inputs: H_df.csv (CM x subtype) plus top subtype/node table.
    # Output: joint_module_top_subtype_heatmap.pdf/svg.
    selected = top_nodes.loc[top_nodes["rank"].astype(int) <= top_n, "cell_subtype"].drop_duplicates().tolist()
    h_frac = H_cm_by_subtype.div(H_cm_by_subtype.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    plot_df = h_frac.loc[:, [x for x in selected if x in h_frac.columns]]
    plot_df = plot_df.loc[sorted(plot_df.index.astype(str), key=cm_sort_key)]
    fig, ax = plt.subplots(figsize=(max(7.0, 0.18 * plot_df.shape[1]), max(3.5, 0.35 * plot_df.shape[0])))
    sns.heatmap(plot_df, cmap="Reds", ax=ax, linewidths=0.2, linecolor="white", cbar_kws={"label": "Loading fraction"})
    ax.set_xlabel("Cell subtype")
    ax.set_ylabel("CM")
    ax.set_title("Top subtype loadings", fontweight="bold")
    ax.tick_params(axis="x", rotation=90)
    fig.tight_layout()
    save_pdf_svg(fig, Path(out_dir) / "joint_module_top_subtype_heatmap")


def plot_qstar_heatmap(value: pd.DataFrame, q: pd.DataFrame, stem: Path, title: str,
                       cbar_label: str, cmap: str = "coolwarm", ns_label: str = "ns"):
    # FIXED: q stars come only from the matching q matrix/table.
    # FIXED: do not annotate numeric values in the heatmap cells.
    value = value.loc[:, sorted(value.columns.astype(str), key=cm_sort_key)]
    q = q.reindex(index=value.index, columns=value.columns)
    annot = q.map(lambda x: q_stars(float(x), ns_label=ns_label))
    vmax = np.nanmax(np.abs(value.to_numpy(dtype=float)))
    if not np.isfinite(vmax) or vmax == 0:
        vmax = 1.0
    fig, ax = plt.subplots(figsize=(max(7.5, 0.55 * value.shape[1] + 2.0),
                                    max(4.8, 0.42 * value.shape[0] + 1.8)))
    sns.heatmap(
        value.astype(float),
        cmap=cmap,
        center=0,
        vmin=-vmax,
        vmax=vmax,
        linewidths=0.25,
        linecolor="white",
        annot=annot,
        fmt="",
        annot_kws={"fontsize": 7, "color": "black"},
        cbar_kws={"label": cbar_label},
        ax=ax,
    )
    ax.set_title(title)
    ax.set_xlabel("CM")
    ax.set_ylabel("Epithelial subtype")
    ax.tick_params(axis="x", labelrotation=90)
    ax.tick_params(axis="y", labelrotation=0)
    fig.tight_layout()
    save_pdf_svg(fig, stem)


def plot_epi_cm_qstar_heatmaps(tumor_rho: pd.DataFrame, tumor_q: pd.DataFrame,
                               normal_rho: pd.DataFrame, normal_q: pd.DataFrame,
                               out_dir: Path, method: str = "spearman"):
    # Spearman uses rho labels; Pearson branch calls the same structure with r labels.
    stat = "Spearman rho" if method == "spearman" else "Pearson r"
    plot_qstar_heatmap(tumor_rho, tumor_q, Path(out_dir) / "balanced_joint_cm_epi_cm_association_tumor_heatmap_qstars",
                       "Balanced joint CM-Epi association (tumor)", stat)
    plot_qstar_heatmap(normal_rho, normal_q, Path(out_dir) / "balanced_joint_cm_epi_cm_association_normal-like_heatmap_qstars",
                       "Balanced joint CM-Epi association (normal-like)", stat)
```

Use these direct input and output mappings when wiring the code above into
`06_final_unified_plotting_and_manifest`:

Every figure in the direct input/output mapping below and every file listed in
`Canonical figure outputs` is mandatory for a complete Module 04 / Epi-CM
discovery plotting run. Do not skip a figure because it is inconvenient, has no
significant result, is empty after filtering, or was not requested by name. If a
required input table is missing, malformed, or empty, stop and report the exact
missing input instead of silently omitting the figure or substituting another
plot. A branch may be omitted only when the user explicitly disables that branch
before execution, and the omission must be recorded in the figure manifest.

```text
joint_nmf_k_selection_metrics.csv
  -> plot_joint_nmf_k_selection
  -> joint_nmf_k_selection.pdf/svg

joint_CM_activity_tumor_vs_normal_mean_sd_summary.csv
  -> plot_activity_mean_sd_barplot
  -> activity_df_tumor_vs_normal_mean_sd_barplot.pdf/svg

activity_df_sample_by_CM.csv or W_df.csv plus sample_status.csv plus sample-to-series mapping
  -> plot_activity_heatmap_variants
  -> w_df_activity_sample_activity_per_CM_{raw,zscore,robust,standard_scale_col}_clustermap.pdf/svg

w_df_activity_sample_by_CM_standard_scale_col.csv plus sample_status.csv plus sample-to-series mapping
  -> plot_activity_minmax_with_series_clustermap
  -> w_df_activity_sample_activity_per_CM_standard_scale_col_with_series_clustermap.pdf/svg
  -> w_df_activity_sample_by_CM_standard_scale_col_with_series_plot_matrix.csv
  -> w_df_activity_sample_by_CM_standard_scale_col_with_series_row_annotations.csv

loading_df_cell_subtype_by_CM.csv or H_df.csv transposed
  -> plot_loading_heatmap_variants
  -> h_df_loading_cell_subtype_weights_per_CM_{raw,zscore,robust,standard_scale_col}_clustermap.pdf/svg

H_df.csv plus joint_cm_cell_subtype_nodes_top*_from_H_df.csv
  -> plot_joint_module_top_subtype_heatmap
  -> joint_module_top_subtype_heatmap.pdf/svg

balanced_joint_cm_epi_cm_association_*_{rho|r}_matrix.csv and matching q matrices
  -> plot_epi_cm_qstar_heatmaps
  -> separate tumor and normal-like association heatmap_qstars PDF/SVG only
```

Nodeplots and all epithelial subtype x CM scatterplots have their own explicit
code blocks below because they are more error-prone; still call them from the
same `06_final_unified_plotting_and_manifest` script and record them in the same
figure manifest.
- The final manifest must include one row per PDF/SVG with columns:
  `figure_file`, `figure_family`, `plotting_function`, `direct_input_tables`,
  `method`, `output_format`, and `notes`.

Canonical figure outputs:

```text
joint_nmf_k_selection.pdf
joint_nmf_k_selection.svg
activity_df_tumor_vs_normal_mean_sd_barplot.pdf
activity_df_tumor_vs_normal_mean_sd_barplot.svg
w_df_activity_sample_activity_per_CM_raw_clustermap.pdf
w_df_activity_sample_activity_per_CM_raw_clustermap.svg
w_df_activity_sample_activity_per_CM_zscore_clustermap.pdf
w_df_activity_sample_activity_per_CM_zscore_clustermap.svg
w_df_activity_sample_activity_per_CM_robust_clustermap.pdf
w_df_activity_sample_activity_per_CM_robust_clustermap.svg
w_df_activity_sample_activity_per_CM_standard_scale_col_clustermap.pdf
w_df_activity_sample_activity_per_CM_standard_scale_col_clustermap.svg
w_df_activity_sample_activity_per_CM_standard_scale_col_with_series_clustermap.pdf
w_df_activity_sample_activity_per_CM_standard_scale_col_with_series_clustermap.svg
h_df_loading_cell_subtype_weights_per_CM_raw_clustermap.pdf
h_df_loading_cell_subtype_weights_per_CM_raw_clustermap.svg
h_df_loading_cell_subtype_weights_per_CM_zscore_clustermap.pdf
h_df_loading_cell_subtype_weights_per_CM_zscore_clustermap.svg
h_df_loading_cell_subtype_weights_per_CM_robust_clustermap.pdf
h_df_loading_cell_subtype_weights_per_CM_robust_clustermap.svg
h_df_loading_cell_subtype_weights_per_CM_standard_scale_col_clustermap.pdf
h_df_loading_cell_subtype_weights_per_CM_standard_scale_col_clustermap.svg
joint_module_top_subtype_heatmap.pdf
joint_module_top_subtype_heatmap.svg
normal_like_all_CM_nodeplot.pdf
normal_like_all_CM_nodeplot.svg
tumor_all_CM_nodeplot.pdf
tumor_all_CM_nodeplot.svg
normal_like_nodeplots_by_cm/<CM>_normal_like_nodeplot.pdf
normal_like_nodeplots_by_cm/<CM>_normal_like_nodeplot.svg
tumor_nodeplots_by_cm/<CM>_tumor_nodeplot.pdf
tumor_nodeplots_by_cm/<CM>_tumor_nodeplot.svg
tumor_centric_nodeplot_edge_origin.pdf
tumor_centric_nodeplot_edge_origin.svg
tumor_centric_nodeplots_by_cm/<CM>_tumor_centric_nodeplot_edge_origin.pdf
tumor_centric_nodeplots_by_cm/<CM>_tumor_centric_nodeplot_edge_origin.svg
tumor_centric_nodeplot_edge_origin_tumor_only_edge_colorbar.pdf
tumor_centric_nodeplot_edge_origin_tumor_only_edge_colorbar.svg
tumor_centric_nodeplot_edge_origin_shared_edge_colorbar.pdf
tumor_centric_nodeplot_edge_origin_shared_edge_colorbar.svg
tumor_centric_nodeplot_edge_origin_edge_class_legend.pdf
tumor_centric_nodeplot_edge_origin_edge_class_legend.svg
normal_like_top10_node_correlation_heatmap_no_edge_filter.pdf
normal_like_top10_node_correlation_heatmap_no_edge_filter.svg
tumor_top10_node_correlation_heatmap_no_edge_filter.pdf
tumor_top10_node_correlation_heatmap_no_edge_filter.svg
top10_node_correlation_heatmaps_no_edge_filter_by_cm/normal_like/<CM>_normal_like_top10_node_correlation_heatmap_no_edge_filter.pdf
top10_node_correlation_heatmaps_no_edge_filter_by_cm/normal_like/<CM>_normal_like_top10_node_correlation_heatmap_no_edge_filter.svg
top10_node_correlation_heatmaps_no_edge_filter_by_cm/tumor/<CM>_tumor_top10_node_correlation_heatmap_no_edge_filter.pdf
top10_node_correlation_heatmaps_no_edge_filter_by_cm/tumor/<CM>_tumor_top10_node_correlation_heatmap_no_edge_filter.svg
balanced_joint_cm_epi_cm_association_tumor_heatmap_qstars.pdf
balanced_joint_cm_epi_cm_association_tumor_heatmap_qstars.svg
balanced_joint_cm_epi_cm_association_normal-like_heatmap_qstars.pdf
balanced_joint_cm_epi_cm_association_normal-like_heatmap_qstars.svg
scatter_<CM>_vs_<epi_subtype>.pdf
scatter_<CM>_vs_<epi_subtype>.svg
```

Final plotting style is defined in `Module Figure Style Contract`. The
canonical workflow must implement that style directly. For each final figure,
keep one plotting block that contains the final style and writes the final
PDF/SVG output.

## Inputs

Minimum:

```text
sample x non-epithelial subtype abundance matrix
sample x epithelial subtype abundance matrix
sample metadata table with sample_id and status
```

Orientation:

```text
rows = samples
columns = cell subtypes
values = within-sample normalized frequencies or comparable abundance values
```

If starting from annotated h5ad, derive sample-by-subtype counts and frequencies first. Do not include epithelial columns in the non-epithelial CM matrix unless explicitly requested.

Optional entry points:

```text
annotated h5ad -> derive epithelial/non-epithelial abundance matrices
abundance matrices -> run NMF and Epi-CM associations
precomputed W/H CM matrices -> skip NMF and run nodeplots/Epi-CM associations
precomputed Epi-CM association and q matrices -> plot heatmaps only
figure-ready CSVs -> plot only
```

Before using supplied intermediates, validate orientation, sample IDs, subtype names, status labels, and whether values are raw proportions, min-max scaled, z-scored, or display-scaled.

## Cell Subtype to Broad Cell Type

If `cell_type` is missing or inconsistent, derive broad cell type from the `cell_subtype` prefix before building epithelial and non-epithelial matrices. Default mapping:

```text
Epi -> Epithelial Cells
T -> T Cells
NK -> NK Cells
Mye -> Myeloid Cells
B -> B Cells
Endo -> Endothelial Cells
S -> Stromal Cells
Mast -> Mast Cells
pDC -> pDC
```

Rules:

1. Parse the prefix before the first underscore in `cell_subtype`.
2. Map the prefix to `cell_type` using the table above unless the user provides a custom mapping.
3. If existing `cell_type` conflicts with the subtype prefix, prefer the subtype-prefix-derived value for CM-Epi abundance matrix construction, and write a conflict table.
4. Use `Epi`-derived epithelial cells only for the epithelial subtype abundance matrix.
5. Use all non-epithelial broad classes for the non-epithelial CM matrix unless the user explicitly excludes a class.

## Derive Frequency Matrices

This step is an input-preparation step for CM discovery, not the Epi-CM
association result. Do not describe raw non-epithelial `cell_subtype` columns as
CMs or "CM subtypes". A CM is a latent module/program produced later by the
balanced joint NMF step. Therefore the canonical Block 03 analysis must not stop
after sample x subtype frequency matrices and must not treat epithelial subtype
vs non-epithelial subtype Spearman correlation as the default Epi-CM discovery
output.

If starting from metadata or h5ad `obs`, compute:

```text
non_epi_subtype_counts.csv
non_epi_subtype_frequency.csv
epi_subtype_frequency.csv
sample_status.csv
sample_compartment_cell_counts.csv
sample_inclusion_exclusion.csv
```

Rules:

1. Count cells by `sample_id x cell_subtype`.
2. Split epithelial and non-epithelial subtypes using `cell_type` or a trusted subtype-prefix mapping.
3. Normalize within each sample so each compartment row sums to 1.
4. Before CM construction, keep only samples that pass the non-epithelial CM eligibility filters.
5. The canonical non-epithelial CM filter is `min_non_epi_cells_per_sample = 50`, applied to total non-epithelial cells per sample before building the CM input matrix.
6. Do not apply epithelial-cell-count filters to balanced joint NMF, K selection, W/H fitting, CM classification, or CM nodeplots. Epithelial-cell-count filters start only when computing Epi-CM association, Epi-CM scatterplots/heatmaps, or epithelial abundance high/low figures.
7. Write `sample_compartment_cell_counts.csv` with total epithelial and non-epithelial cells per sample before filtering.
8. Write `sample_inclusion_exclusion.csv` with one row per candidate sample, boolean keep flags, and explicit exclusion reasons. Use `keep_for_cm` for non-epithelial CM construction and `keep_for_epi_cm` for epithelial association. Do not merge those flags into one early filter.
9. Keep row sample IDs identical across non-epithelial and metadata tables after CM filtering. For Epi-CM association, intersect canonical CM activity samples with `keep_for_epi_cm` samples only at the association step.
10. Do not use display-scaled values for statistics.

## Balanced Joint NMF

Use non-epithelial subtype frequencies, then apply column-wise min-max
normalization across samples before NMF. Each non-epithelial subtype column is
scaled independently:

```text
frequency = sample x non-epithelial subtype frequency
V = column_minmax(frequency)
V ~= W @ H
W = sample x CM activity
H = CM x non-epithelial subtype loading
```

Column-wise min-max normalization is mandatory for the canonical balanced joint
NMF input:

```text
column_min[subtype] = min(frequency[:, subtype] over samples)
column_range[subtype] = max(frequency[:, subtype] over samples) - column_min[subtype]
V[:, subtype] = (frequency[:, subtype] - column_min[subtype]) / column_range[subtype]
if column_range[subtype] == 0, set that normalized subtype column to 0
```

Save both `non_epi_subtype_frequency_column_minmax.csv` and
`column_minmax_params.csv`. The parameter table must have one row per
non-epithelial subtype with at least `min` and `range` columns. Do not use raw
frequency, per-row normalization, whole-matrix scalar min-max, z-score
normalization, or L2 normalization as the NMF input unless the user explicitly
requests a method change.

Naming rule: code generated from this skill must use names such as
`column_minmax`, `column_min`, `column_range`, and `V_column_minmax` for this
operation. Do not name the helper, variables, or parameter files
`global_minmax`, because that name incorrectly suggests whole-matrix scalar
min-max normalization.

Required Python implementation pattern:

```python
def column_minmax_normalize(freq: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Column-wise min-max normalize sample x subtype frequencies.

    Rows are samples and columns are non-epithelial subtypes. Each subtype
    column is scaled independently over samples. This is not whole-matrix
    min-max normalization.
    """
    column_min = freq.min(axis=0)
    column_range = (freq.max(axis=0) - column_min).replace(0, np.nan)
    V_column_minmax = ((freq - column_min) / column_range).fillna(0.0).clip(lower=0.0)
    column_minmax_params = pd.DataFrame({
        "min": column_min,
        "range": column_range.fillna(0.0),
    })
    column_minmax_params.index.name = "cell_subtype"
    return V_column_minmax.astype(float), column_minmax_params


V_column_minmax, column_minmax_params = column_minmax_normalize(non_epi_frequency)
V_column_minmax.to_csv(out_dir / "non_epi_subtype_frequency_column_minmax.csv")
column_minmax_params.to_csv(out_dir / "column_minmax_params.csv")
```

Do not use `DataFrame.where()` with a column-indexed Series mask for this step.
That pattern can align the mask to rows instead of columns and silently turn the
whole normalized matrix into zeros/NaNs. The forbidden pattern is:

```python
# Forbidden for column-wise minmax.
V_norm = V_norm.where(column_range > 0, 0.0)
```

If using `where`, the mask must be explicitly aligned to columns, but the
preferred implementation is the division-by-NaN plus `fillna(0.0)` snippet
above. After normalization, verify that non-zero-variance columns remain before
running NMF:

```python
valid_columns = V_column_minmax.columns[V_column_minmax.var(axis=0) > 0]
if len(valid_columns) == 0:
    raise ValueError(
        "Column-minmax normalization produced zero valid non-epithelial subtype columns; "
        "check column alignment and the non-epithelial frequency input before NMF."
    )
```

Forbidden implementation pattern:

```python
v_min = freq.values.min()
v_max = freq.values.max()
V = (freq - v_min) / (v_max - v_min)
```

The forbidden pattern collapses all subtype columns onto one whole-matrix scale
and can change the selected CM rank. Do not use it for the canonical Module 04
balanced joint NMF.

Status balancing is required by default when both status groups are present:

```text
normal_group_total_weight = 0.5
tumor_group_total_weight = 0.5
weight(sample) = group_total_weight / number_of_samples_in_that_status_group
weight(sample) = weight(sample) / mean(weight over all included samples)
V_weighted = sqrt(weight(sample))[:, None] * V
```

Use the square root of sample weights when constructing `V_weighted` so the
Frobenius objective corresponds to weighted squared reconstruction error by
sample. The rescale-to-mean-1 step is mandatory and should be saved in
`group_balanced_sample_weights.csv`. This is not row normalization.

If there are more than two status groups, ask the user for the balancing plan.
If the user explicitly requests unweighted NMF, report that method change.

The canonical final CM model does not hard-code K. Run the K-selection step as
part of the analysis, use its selected K for the final NMF, and save the
K-selection diagnostic figure/table for provenance. Do not skip K-selection.

Default K-selection logic:

```text
candidate K range = 2 to 20 inclusive
actual max K = min(20, number of samples, number of non-epithelial subtypes)
random seeds = 0, 1, 2, 3, 4
input matrix = group-balanced weighted column-minmax non-epithelial frequency matrix
```

Canonical configuration pattern:

```python
# FIXED defaults unless the user explicitly changes them.
module_k = None
module_k_range = (2, 20)
rank_selection_seeds = (0, 1, 2, 3, 4)
random_state = 0

nmf_max_iter = 3000
nmf_tol = 1e-5
nmf_alpha_w = 0.0
nmf_alpha_h = 1e-3
nmf_l1_ratio = 0.1

normal_group_total_weight = 0.5
tumor_group_total_weight = 0.5
```

Implementation must derive the actual K list from `module_k_range` and the
active matrix dimensions:

```python
max_rank = min(module_k_range[1], V_weighted.shape[0], V_weighted.shape[1])
min_rank = max(2, module_k_range[0])
if max_rank < min_rank:
    raise ValueError(f"Invalid joint rank range after matrix size check: {min_rank}..{max_rank}")
candidate_ks = list(range(min_rank, max_rank + 1))
```

Do not hard-code `range(2, 11)`, `2..10`, three seeds, or alternate seed values
such as `[42, 123, 2024]` for the canonical branch. If a script reports K=10
while `actual max K` is greater than 10, treat that run as truncated and rerun
K-selection with the canonical K range.

Canonical balanced joint NMF analysis code pattern:

```python
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
import json
import math
import warnings

import numpy as np
import pandas as pd
from scipy.optimize import nnls
from scipy.optimize import linear_sum_assignment
from sklearn.decomposition import NMF
from sklearn.metrics.pairwise import cosine_similarity


@dataclass(frozen=True)
class CMNMFConfig:
    # REPLACEABLE input columns and labels.
    sample_col: str = "sample"
    status_col: str = "status"
    celltype_col: str = "cell_type"
    subtype_col: str = "cell_subtype"
    epi_label: str = "Epi"
    normal_status: str = "normal-like"
    tumor_status: str = "tumor"

    # FIXED defaults unless the user explicitly changes them.
    min_non_epi_cells_per_sample: int = 50
    forced_cm_k: int | None = None
    cm_k_range: tuple[int, int] = (2, 20)
    rank_selection_seeds: tuple[int, ...] = (0, 1, 2, 3, 4)
    nmf_max_iter: int = 3000
    nmf_tol: float = 1e-5
    nmf_alpha_w: float = 0.0
    nmf_alpha_h: float = 1e-3
    nmf_l1_ratio: float = 0.1
    normal_group_total_weight: float = 0.5
    tumor_group_total_weight: float = 0.5

    # FIXED CM class rules unless the user explicitly changes them.
    normal_specific_max_ratio: float = 0.5
    tumor_specific_min_ratio: float = 2.0
    min_active_fraction_for_specific: float = 0.05


def write_json(obj: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(obj, handle, indent=2, ensure_ascii=False)


def build_non_epi_and_epi_frequencies(obs: pd.DataFrame, cfg: CMNMFConfig):
    required = {cfg.sample_col, cfg.status_col, cfg.celltype_col, cfg.subtype_col}
    missing = required - set(obs.columns)
    if missing:
        raise ValueError(f"Missing required obs columns: {sorted(missing)}")

    clean = obs.dropna(subset=list(required)).copy()
    is_epi = clean[cfg.celltype_col].astype(str).eq(cfg.epi_label)
    non_epi = clean.loc[~is_epi].copy()
    epi = clean.loc[is_epi].copy()
    if non_epi.empty:
        raise ValueError("No non-epithelial cells remain for CM NMF.")

    non_epi_cells = non_epi.groupby(cfg.sample_col, observed=True).size().rename("non_epi_cells")
    keep_samples = non_epi_cells.index[non_epi_cells >= cfg.min_non_epi_cells_per_sample]
    non_epi = non_epi.loc[non_epi[cfg.sample_col].isin(keep_samples)].copy()
    if non_epi.empty:
        raise ValueError("No samples pass min_non_epi_cells_per_sample.")

    non_epi_counts = pd.crosstab(non_epi[cfg.sample_col], non_epi[cfg.subtype_col]).astype(float)
    non_epi_counts = non_epi_counts.sort_index(axis=0).sort_index(axis=1)
    non_epi_frequency = non_epi_counts.div(non_epi_counts.sum(axis=1), axis=0).fillna(0.0)

    status_counts = pd.crosstab(non_epi[cfg.sample_col], non_epi[cfg.status_col]).reindex(non_epi_frequency.index).fillna(0)
    sample_status = pd.DataFrame(index=non_epi_frequency.index)
    sample_status["status"] = status_counts.idxmax(axis=1)
    sample_status["status_majority_fraction"] = status_counts.max(axis=1) / status_counts.sum(axis=1)
    sample_status["non_epi_cells"] = non_epi_counts.sum(axis=1).astype(int)

    epi_counts = pd.crosstab(epi[cfg.sample_col], epi[cfg.subtype_col]).astype(float)
    epi_counts = epi_counts.reindex(non_epi_frequency.index).fillna(0.0)
    epi_frequency = epi_counts.div(epi_counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    epi_frequency = epi_frequency.sort_index(axis=0).sort_index(axis=1)
    return non_epi_counts, non_epi_frequency, epi_frequency, sample_status


def column_minmax_normalize(freq: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    # FIXED: each non-epithelial subtype column is scaled independently over samples.
    column_min = freq.min(axis=0)
    column_range = (freq.max(axis=0) - column_min).replace(0, np.nan)
    V_column_minmax = ((freq - column_min) / column_range).fillna(0.0).clip(lower=0.0)
    params = pd.DataFrame({"min": column_min, "range": column_range.fillna(0.0)})
    params.index.name = "cell_subtype"
    if (V_column_minmax.var(axis=0) > 0).sum() == 0:
        raise ValueError("Column-minmax matrix has no non-zero-variance subtype columns.")
    return V_column_minmax.astype(float), params


def make_status_balanced_weights(sample_status: pd.DataFrame, cfg: CMNMFConfig) -> pd.Series:
    status = sample_status["status"].astype(str)
    normal_mask = status.eq(cfg.normal_status)
    tumor_mask = status.eq(cfg.tumor_status)
    if normal_mask.sum() == 0 or tumor_mask.sum() == 0:
        raise ValueError("Need both normal-like and tumor samples for status-balanced CM NMF.")
    weights = pd.Series(0.0, index=sample_status.index, name="joint_nmf_row_weight")
    weights.loc[normal_mask] = cfg.normal_group_total_weight / normal_mask.sum()
    weights.loc[tumor_mask] = cfg.tumor_group_total_weight / tumor_mask.sum()
    return weights / weights.mean()


def make_nmf(n_components: int, seed: int, cfg: CMNMFConfig) -> NMF:
    return NMF(
        n_components=n_components,
        init="nndsvda",
        solver="cd",
        beta_loss="frobenius",
        random_state=seed,
        max_iter=cfg.nmf_max_iter,
        tol=cfg.nmf_tol,
        alpha_W=cfg.nmf_alpha_w,
        alpha_H=cfg.nmf_alpha_h,
        l1_ratio=cfg.nmf_l1_ratio,
    )


def row_l1_normalize(matrix: np.ndarray) -> np.ndarray:
    denom = matrix.sum(axis=1, keepdims=True)
    return np.divide(matrix, denom, out=np.zeros_like(matrix, dtype=float), where=denom != 0)


def matched_cosine_mean(H_a: np.ndarray, H_b: np.ndarray) -> float:
    sim = cosine_similarity(row_l1_normalize(H_a), row_l1_normalize(H_b))
    row_ind, col_ind = linear_sum_assignment(-sim)
    return float(sim[row_ind, col_ind].mean())


def stability_matched_cosine(H_by_seed: list[np.ndarray]) -> float:
    if len(H_by_seed) < 2:
        return float("nan")
    values = [
        matched_cosine_mean(H_by_seed[i], H_by_seed[j])
        for i in range(len(H_by_seed))
        for j in range(i + 1, len(H_by_seed))
    ]
    return float(np.mean(values))


def evaluate_cm_rank(V_weighted: np.ndarray, cfg: CMNMFConfig) -> tuple[pd.DataFrame, int]:
    # FIXED: derive candidate K dynamically from matrix dimensions and cm_k_range.
    max_rank = min(cfg.cm_k_range[1], V_weighted.shape[0], V_weighted.shape[1])
    min_rank = max(2, cfg.cm_k_range[0])
    if max_rank < min_rank:
        raise ValueError(f"Invalid CM K range after matrix size check: {min_rank}..{max_rank}")
    candidate_ks = list(range(min_rank, max_rank + 1))
    denom = float(np.square(V_weighted).sum())
    if denom == 0:
        raise ValueError("Weighted NMF input is all zero; check frequency/minmax inputs.")

    rows = []
    for k in candidate_ks:
        H_by_seed, seed_errors, seed_explained = [], [], []
        for seed in cfg.rank_selection_seeds:
            model = make_nmf(n_components=k, seed=seed, cfg=cfg)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=UserWarning)
                W = model.fit_transform(V_weighted)
            H = model.components_
            err = float(np.square(V_weighted - W @ H).sum())
            H_by_seed.append(H)
            seed_errors.append(err)
            seed_explained.append(1.0 - err / denom)
        best_i = int(np.argmin(seed_errors))
        rows.append({
            "k": k,
            "n_seeds": len(cfg.rank_selection_seeds),
            "mean_balanced_explained_fraction": float(np.nanmean(seed_explained)),
            "best_balanced_explained_fraction": float(np.nanmax(seed_explained)),
            "mean_reconstruction_error": float(np.mean(seed_errors)),
            "best_reconstruction_error": float(np.min(seed_errors)),
            "stability_matched_cosine": stability_matched_cosine(H_by_seed),
            "best_seed": int(cfg.rank_selection_seeds[best_i]),
        })

    metrics = pd.DataFrame(rows)
    metrics["selection_score"] = (
        metrics["best_balanced_explained_fraction"].fillna(0.0)
        + 0.05 * metrics["stability_matched_cosine"].fillna(0.0)
        - 0.01 * metrics["k"]
    )
    if cfg.forced_cm_k is not None:
        selected_k = int(cfg.forced_cm_k)
        if selected_k not in metrics["k"].tolist():
            raise ValueError(f"forced_cm_k={selected_k} outside evaluated ranks: {metrics['k'].tolist()}")
    else:
        selected_k = int(
            metrics.sort_values(
                ["selection_score", "best_balanced_explained_fraction"],
                ascending=False,
            ).iloc[0]["k"]
        )
    metrics["selected"] = metrics["k"].eq(selected_k)
    return metrics, selected_k


def fit_best_seed_nmf(V_weighted: np.ndarray, k: int, seeds: Iterable[int], cfg: CMNMFConfig):
    best_error = math.inf
    best_W = best_H = None
    for seed in seeds:
        model = make_nmf(n_components=k, seed=seed, cfg=cfg)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            W = model.fit_transform(V_weighted)
        H = model.components_
        err = float(np.square(V_weighted - W @ H).sum())
        if err < best_error:
            best_error, best_W, best_H = err, W, H
    if best_W is None or best_H is None:
        raise RuntimeError("Final NMF failed for all seeds.")
    return best_W, best_H


def refit_activity_by_nnls(H: np.ndarray, V_unweighted: np.ndarray) -> np.ndarray:
    # FIXED: canonical W/activity is NNLS refit on unweighted column-minmax V.
    design = H.T
    W = np.zeros((V_unweighted.shape[0], H.shape[0]), dtype=float)
    for i in range(V_unweighted.shape[0]):
        W[i, :] = nnls(design, V_unweighted[i, :], maxiter=design.shape[1] * 10)[0]
    return W


def classify_raw_cms(W: pd.DataFrame, sample_status: pd.DataFrame, cfg: CMNMFConfig) -> pd.DataFrame:
    status = sample_status["status"].astype(str)
    normal = W.loc[status.eq(cfg.normal_status)]
    tumor = W.loc[status.eq(cfg.tumor_status)]
    if normal.empty or tumor.empty:
        raise ValueError("Need both normal-like and tumor samples to classify CMs.")
    rows = []
    for raw_cm in W.columns:
        normal_mean = float(normal[raw_cm].mean())
        tumor_mean = float(tumor[raw_cm].mean())
        normal_active = float((normal[raw_cm] > 1e-8).mean())
        tumor_active = float((tumor[raw_cm] > 1e-8).mean())
        ratio = tumor_mean / normal_mean if normal_mean > 0 else np.inf
        if ratio >= cfg.tumor_specific_min_ratio and tumor_active >= cfg.min_active_fraction_for_specific:
            cm_class = "tCM"
        elif ratio <= cfg.normal_specific_max_ratio and normal_active >= cfg.min_active_fraction_for_specific:
            cm_class = "normalCM"
        else:
            cm_class = "sharedCM"
        rows.append({
            "raw_component": raw_cm,
            "raw_component_order": int(str(raw_cm).replace("CM", "")),
            "class": cm_class,
            "normal_mean_usage": normal_mean,
            "tumor_mean_usage": tumor_mean,
            "tumor_to_normal_mean_ratio": ratio,
            "normal_active_fraction": normal_active,
            "tumor_active_fraction": tumor_active,
        })
    return pd.DataFrame(rows).sort_values("raw_component_order", kind="stable")


def assign_canonical_cm_names(classification: pd.DataFrame) -> pd.DataFrame:
    prefix = {"sharedCM": "s", "normalCM": "n", "tCM": "t"}
    out = classification.sort_values("raw_component_order", kind="stable").copy()
    out["global_order"] = range(1, len(out) + 1)
    out["class_prefix"] = out["class"].map(prefix)
    if out["class_prefix"].isna().any():
        raise ValueError("Unknown CM class in classification table.")
    out["CM"] = out["class_prefix"] + "_CM" + out["global_order"].astype(str)
    return out


def run_balanced_joint_cm_nmf(non_epi_frequency: pd.DataFrame, sample_status: pd.DataFrame, out_dir: Path, cfg: CMNMFConfig):
    V_column_minmax, minmax_params = column_minmax_normalize(non_epi_frequency)
    sample_weights = make_status_balanced_weights(sample_status.loc[V_column_minmax.index], cfg)
    V = V_column_minmax.to_numpy(dtype=float)
    V_weighted = V * np.sqrt(sample_weights.to_numpy(dtype=float))[:, None]

    metrics, selected_k = evaluate_cm_rank(V_weighted, cfg)
    best_seed = int(metrics.loc[metrics["k"].eq(selected_k), "best_seed"].iloc[0])
    _, H_raw = fit_best_seed_nmf(V_weighted, selected_k, seeds=(best_seed,), cfg=cfg)
    W_raw = refit_activity_by_nnls(H_raw, V)

    raw_cms = [f"CM{i + 1}" for i in range(selected_k)]
    W_raw_df = pd.DataFrame(W_raw, index=V_column_minmax.index, columns=raw_cms)
    H_raw_df = pd.DataFrame(H_raw, index=raw_cms, columns=V_column_minmax.columns)
    classification = assign_canonical_cm_names(classify_raw_cms(W_raw_df, sample_status.loc[W_raw_df.index], cfg))
    raw_to_cm = dict(zip(classification["raw_component"], classification["CM"]))
    cm_order = classification["CM"].tolist()

    W_df = W_raw_df.rename(columns=raw_to_cm).loc[:, cm_order]
    H_df = H_raw_df.rename(index=raw_to_cm).loc[cm_order, :]
    loading_df = H_df.T.copy()
    loading_fraction = loading_df.div(loading_df.sum(axis=0).replace(0, np.nan), axis=1).fillna(0.0)
    activity_df = W_df.join(sample_status[["status", "non_epi_cells"]])

    out_dir.mkdir(parents=True, exist_ok=True)
    V_column_minmax.to_csv(out_dir / "non_epi_subtype_frequency_column_minmax.csv")
    minmax_params.to_csv(out_dir / "column_minmax_params.csv")
    sample_weights.to_frame().join(sample_status[["status", "non_epi_cells"]]).to_csv(out_dir / "group_balanced_sample_weights.csv")
    metrics.to_csv(out_dir / "joint_nmf_k_selection_metrics.csv", index=False)
    write_json({
        "selected_module_k": int(selected_k),
        "best_seed": best_seed,
        "candidate_Ks": metrics["k"].astype(int).tolist(),
        "rank_selection_seeds": list(cfg.rank_selection_seeds),
    }, out_dir / "selected_module_k.json")
    classification.to_csv(out_dir / "joint_module_classification.csv", index=False)
    classification[["raw_component", "raw_component_order", "CM", "class", "class_prefix", "global_order"]].to_csv(
        out_dir / "raw_to_canonical_CM_mapping.csv", index=False
    )
    W_df.to_csv(out_dir / "W_df.csv")
    H_df.to_csv(out_dir / "H_df.csv")
    loading_df.to_csv(out_dir / "loading_df_cell_subtype_by_CM.csv")
    loading_fraction.to_csv(out_dir / "loading_df_cell_subtype_by_CM_fraction.csv")
    activity_df.to_csv(out_dir / "activity_df_sample_by_CM.csv")
    return W_df, H_df, loading_df, activity_df, metrics, classification


def run_single_status_joint_cm_nmf(
    non_epi_frequency: pd.DataFrame,
    sample_status: pd.DataFrame,
    out_dir: Path,
    cfg: CMNMFConfig,
    detected_mode: str,
):
    """ADDITIVE route for a tumor-only or normal-like-only cohort.

    FIXED: keep the original NMF/K-selection/NNLS machinery, use equal sample
    weights, and map raw CM1..CMK to t_CM1..t_CMK for tumor-only or
    n_CM1..n_CMK for normal-only.
    """
    mode_to_status = {
        "tumor_only": cfg.tumor_status,
        "normal_only": cfg.normal_status,
    }
    if detected_mode not in mode_to_status:
        raise ValueError(
            "run_single_status_joint_cm_nmf requires tumor_only or normal_only; "
            f"received {detected_mode!r}"
        )
    present_status = mode_to_status[detected_mode]
    absent_status = (
        cfg.normal_status if detected_mode == "tumor_only" else cfg.tumor_status
    )
    if detected_mode == "tumor_only":
        final_prefix, final_class = "t", "tCM"
    else:
        final_prefix, final_class = "n", "normalCM"
    aligned_status = sample_status.loc[non_epi_frequency.index].copy()
    observed = set(aligned_status["status"].astype(str))
    if observed != {present_status}:
        raise ValueError(
            "Single-status NMF input does not match the detected route; "
            f"observed statuses={sorted(observed)}"
        )
    if "non_epi_cells" not in aligned_status.columns:
        raise ValueError("sample_status must contain non_epi_cells.")

    V_column_minmax, minmax_params = column_minmax_normalize(non_epi_frequency)
    V = V_column_minmax.to_numpy(dtype=float)
    sample_weights = pd.Series(
        1.0,
        index=V_column_minmax.index,
        name="joint_nmf_row_weight",
    )
    # Equal weights make V_weighted exactly equal to the original V.
    V_weighted = V

    metrics, selected_k = evaluate_cm_rank(V_weighted, cfg)
    best_seed = int(metrics.loc[metrics["k"].eq(selected_k), "best_seed"].iloc[0])
    _, H_raw = fit_best_seed_nmf(
        V_weighted,
        selected_k,
        seeds=(best_seed,),
        cfg=cfg,
    )
    W_raw = refit_activity_by_nnls(H_raw, V)

    raw_cms = [f"CM{i + 1}" for i in range(selected_k)]
    cm_order = [f"{final_prefix}_CM{i + 1}" for i in range(selected_k)]
    raw_to_cm = dict(zip(raw_cms, cm_order))
    W_raw_df = pd.DataFrame(W_raw, index=V_column_minmax.index, columns=raw_cms)
    H_raw_df = pd.DataFrame(H_raw, index=raw_cms, columns=V_column_minmax.columns)
    W_df = W_raw_df.rename(columns=raw_to_cm).loc[:, cm_order]
    H_df = H_raw_df.rename(index=raw_to_cm).loc[cm_order, :]
    loading_df = H_df.T.copy()
    loading_fraction = loading_df.div(
        loading_df.sum(axis=0).replace(0, np.nan), axis=1
    ).fillna(0.0)
    activity_df = W_df.join(aligned_status[["status", "non_epi_cells"]])

    classification = pd.DataFrame(
        {
            "raw_component": raw_cms,
            "raw_component_order": np.arange(1, selected_k + 1, dtype=int),
            "CM": cm_order,
            "class": final_class,
            "class_prefix": final_prefix,
            "global_order": np.arange(1, selected_k + 1, dtype=int),
            "classification_available": False,
            "classification_basis": "single_status_presence_fallback",
            "classification_reason": f"{absent_status} samples absent",
            "tumor_mean": np.nan,
            "normal_like_mean": np.nan,
            "delta": np.nan,
            "p": np.nan,
            "q": np.nan,
        }
    )
    present_mean_col = (
        "tumor_mean" if detected_mode == "tumor_only" else "normal_like_mean"
    )
    classification[present_mean_col] = [
        float(W_df[cm].mean()) for cm in cm_order
    ]

    out_dir.mkdir(parents=True, exist_ok=True)
    V_column_minmax.to_csv(out_dir / "non_epi_subtype_frequency_column_minmax.csv")
    minmax_params.to_csv(out_dir / "column_minmax_params.csv")
    weight_audit = sample_weights.to_frame().join(
        aligned_status[["status", "non_epi_cells"]]
    )
    weight_audit["status_balance_applied"] = False
    weight_audit["weight_reason"] = f"equal weights for {detected_mode} cohort"
    weight_audit.to_csv(out_dir / "group_balanced_sample_weights.csv")
    metrics.to_csv(out_dir / "joint_nmf_k_selection_metrics.csv", index=False)
    write_json(
        {
            "selected_module_k": int(selected_k),
            "best_seed": best_seed,
            "candidate_Ks": metrics["k"].astype(int).tolist(),
            "rank_selection_seeds": list(cfg.rank_selection_seeds),
            "detected_mode": detected_mode,
            "present_status": present_status,
            "absent_status": absent_status,
            "status_balance_applied": False,
            "cm_classification_available": False,
            "final_cm_naming": (
                "t_CM1..t_CMK in raw component order"
                if detected_mode == "tumor_only"
                else "n_CM1..n_CMK in raw component order"
            ),
            "classification_basis": "single_status_presence_fallback",
        },
        out_dir / "selected_module_k.json",
    )
    classification.to_csv(out_dir / "joint_module_classification.csv", index=False)
    classification[
        [
            "raw_component",
            "raw_component_order",
            "CM",
            "class",
            "class_prefix",
            "global_order",
            "classification_available",
            "classification_basis",
            "classification_reason",
        ]
    ].to_csv(out_dir / "raw_to_canonical_CM_mapping.csv", index=False)
    W_df.to_csv(out_dir / "W_df.csv")
    H_df.to_csv(out_dir / "H_df.csv")
    W_df.T.to_csv(out_dir / "activity_df_CM_by_sample.csv")
    loading_df.to_csv(out_dir / "loading_df_cell_subtype_by_CM.csv")
    loading_fraction.to_csv(out_dir / "loading_df_cell_subtype_by_CM_fraction.csv")
    activity_df.to_csv(out_dir / "activity_df_sample_by_CM.csv")

    skipped_steps = [
        "status-balanced weighting",
        "status-derived sharedCM/normalCM/tCM classification",
        "activity_df_tumor_vs_normal_mean_sd_barplot",
        "tumor_centric_nodeplot_edge_origin",
        f"{absent_status} nodeplots and correlation heatmaps",
        f"{absent_status} Epi-CM associations",
    ]
    skipped = pd.DataFrame({"output_or_step": skipped_steps})
    skipped["status"] = "skipped"
    skipped["reason"] = f"{absent_status} samples absent; structurally inapplicable"
    skipped.to_csv(out_dir / "single_status_skipped_outputs.csv", index=False)
    return W_df, H_df, loading_df, activity_df, metrics, classification


def run_cm_nmf_by_detected_mode(
    non_epi_frequency: pd.DataFrame,
    sample_status: pd.DataFrame,
    out_dir: Path,
    cfg: CMNMFConfig,
    detected_mode: str,
):
    """FIXED router; do not modify either analysis branch inside the router."""
    if detected_mode == "tumor_normal":
        # Preserve the original two-status route exactly.
        return run_balanced_joint_cm_nmf(
            non_epi_frequency, sample_status, out_dir, cfg
        )
    if detected_mode in {"tumor_only", "normal_only"}:
        return run_single_status_joint_cm_nmf(
            non_epi_frequency,
            sample_status,
            out_dir,
            cfg,
            detected_mode,
        )
    raise ValueError(f"Unsupported detected_mode: {detected_mode}")


def block03_status_output_policy(detected_mode: str) -> dict[str, object]:
    """FIXED dispatch for nodeplots, associations, and status comparisons."""
    if detected_mode == "tumor_normal":
        return {
            "status_contexts": ("tumor", "normal-like"),
            "run_status_comparisons": True,
        }
    if detected_mode == "tumor_only":
        return {
            "status_contexts": ("tumor",),
            "run_status_comparisons": False,
        }
    if detected_mode == "normal_only":
        return {
            "status_contexts": ("normal-like",),
            "run_status_comparisons": False,
        }
    raise ValueError(f"Unsupported detected_mode: {detected_mode}")


# FIXED generated analysis/plotting dispatch:
# policy = block03_status_output_policy(detected_mode)
# for status_context in policy["status_contexts"]:
#     run_node_tables_and_nodeplots_for_present_status(status_context)
#     run_epi_cm_association_and_plots_for_present_status(status_context)
# if policy["run_status_comparisons"]:
#     run_tumor_vs_normal_activity_barplot()
#     run_edge_origin_comparison()
# For tumor_only and normal_only, the two comparison calls above stay
# intentionally commented out/bypassed and are recorded in
# single_status_skipped_outputs.csv.
```

The code above intentionally uses `non_epi`, `CM`, `column_minmax`, `forced_cm_k`,
and `activity/loading` terminology. Do not reintroduce misleading names such as
`nonepi`, `global_minmax`, `module` for final CM labels, or raw `CM1..CMK`
outside the internal raw-component and mapping tables.

Before K-selection, validate and record the actual matrix dimensions used for
NMF:

```text
n_samples = number of eligible samples
n_subtypes_total = number of non-epithelial subtype columns before variance filtering
n_subtypes_valid = number of non-zero-variance columns after column-minmax
actual max K = min(20, n_samples, n_subtypes_valid)
```

If `n_subtypes_valid == 0`, or if `actual max K < 2`, stop with an explicit
error and do not write partial W/H/K-selection outputs. This indicates an input
or normalization bug, not a biological result. In particular, do not interpret
an empty or all-zero column-minmax matrix as evidence that the dataset has only
two CMs.

For each candidate K and seed, fit NMF on the weighted matrix with:

```text
init = "nndsvda"
solver = "cd"
beta_loss = "frobenius"
max_iter = 3000
tol = 1e-5
alpha_W = 0.0
alpha_H = 0.001
l1_ratio = 0.1
random_state = current seed
```

This uses Frobenius/L2 reconstruction loss on `V_weighted`. Do not apply
explicit L2 normalization to input rows, W columns, or H rows/columns unless the
user explicitly requests a method change. Sklearn NMF does not mean-center the
input matrix; do not mean-center because NMF requires nonnegative inputs. With
sklearn-style NMF regularization, `alpha_H = 0.001` and `l1_ratio = 0.1` means
the H/loading regularization contains 10% L1 and 90% L2 penalty; `alpha_W =
0.0` means no W/activity regularization. Do not add a separate extra L2 penalty
or change these values unless the user explicitly requests a method change.

For each K, calculate:

```text
reconstruction_error = sum((V_weighted - W @ H)^2)
balanced_explained_fraction = 1 - reconstruction_error / sum(V_weighted^2)
mean_balanced_explained_fraction across seeds
best_balanced_explained_fraction across seeds
mean_reconstruction_error across seeds
best_reconstruction_error across seeds
best_seed = seed with minimum reconstruction_error
stability_matched_cosine
```

Required reconstruction-error implementation:

```python
model = NMF(...)
W = model.fit_transform(V_weighted)
H = model.components_
reconstruction_error = float(np.square(V_weighted - W @ H).sum())
balanced_explained_fraction = 1.0 - reconstruction_error / float(np.square(V_weighted).sum())
```

Do not use `model.reconstruction_err_` directly as `reconstruction_error` in
the K-selection table or selection score. In scikit-learn, `reconstruction_err_`
is the Frobenius norm, equal to `sqrt(sum((V_weighted - W @ H)^2))`, while this
workflow requires the squared error `sum((V_weighted - W @ H)^2)`. Using the
unsquared norm falsely inflates explained fraction and can force selection
toward very small K.

`stability_matched_cosine` is computed by row-normalizing H components from each seed, comparing seed pairs by cosine similarity, matching components with Hungarian assignment, and averaging matched cosine similarities.

Select K with:

```text
selection_score =
  best_balanced_explained_fraction
  + 0.05 * stability_matched_cosine
  - 0.01 * K
```

Choose the K with the highest `selection_score`; break ties by highest
`best_balanced_explained_fraction`. Do not add an extra lower-K tie-breaker or
replace this score with reconstruction error alone. The final NMF must use the
K selected by this report, unless the user explicitly forces another
`module_k`. If the user forces K, still run and save the K-selection report when
feasible, then mark the forced K in `selected_module_k.json` and the run
report. Do not silently change the candidate K range, seeds, NMF parameters, or
scoring formula for speed.

Save K-selection outputs under `tables/`, for example:

```text
joint_nmf_k_selection_metrics.csv
selected_module_k.json
run_parameters.txt
```

K-selection diagnostics are generated by
`06_final_unified_plotting_and_manifest` from
`joint_nmf_k_selection_metrics.csv`, for example:

```text
joint_nmf_k_selection.pdf
joint_nmf_k_selection.svg
```

The final K is data-driven from the K-selection report and can be adjusted only
when the user explicitly supplies a fixed K or approves a documented
biological/statistical rationale after reviewing the K-selection report. Record
the final K, candidate K values, number of seeds, random seeds, metrics used,
selection score formula, best seed, and reason for accepting the final K. Fit
the final NMF with the selected K and its `best_seed` from
`joint_nmf_k_selection_metrics.csv` on `V_weighted`, then refit all samples by
NNLS on the unweighted column-minmax matrix `V` to the fixed H basis. The saved
canonical `W_df.csv` and downstream activity tables must come from this NNLS
refit, not from the weighted fit matrix.

Raw NMF component IDs (`CM1`, `CM2`, ..., or equivalent) are allowed only inside
the fitting function and in `raw_to_canonical_CM_mapping.csv`. Immediately after
the NNLS refit, classify the raw components as sharedCM, normalCM, or tCM and
assign final canonical IDs (`s_CM<number>`, `n_CM<number>`, `t_CM<number>`)
before writing any table that later steps may read. In other words, the first
public W/H/activity/loading tables should already use canonical CM IDs. If the
implementation saves raw W/H for debugging/provenance, place them in explicitly
named internal/provenance tables such as `W_df_raw_components.csv` and
`H_df_raw_components.csv`; never name raw-component tables `W_df.csv`,
`H_df.csv`, `activity_df_sample_by_CM.csv`, or
`loading_df_cell_subtype_by_CM.csv`.

Save matrix outputs with this orientation:

```text
W_df.csv = sample x CM
H_df.csv = CM x cell_subtype
loading_df_cell_subtype_by_CM.csv = cell_subtype x CM, equal to H_df transpose
loading_df_cell_subtype_by_CM_fraction.csv = cell_subtype x CM, column-normalized loading fractions
activity_df_CM_by_sample.csv = CM x sample, equal to W_df transpose
activity_df_sample_by_CM.csv = sample x CM plus sample status and nonepi cell counts
h_df_loading_cell_subtype_by_CM_raw.csv = CM x cell_subtype, equal to H_df raw values
h_df_loading_cell_subtype_by_CM_zscore.csv = CM x cell_subtype
h_df_loading_cell_subtype_by_CM_robust.csv = CM x cell_subtype
h_df_loading_cell_subtype_by_CM_standard_scale_col.csv = CM x cell_subtype
```

In this orientation contract, `CM` means final canonical CM ID, not raw NMF
component ID. Therefore `W_df.csv`, `H_df.csv`,
`activity_df_sample_by_CM.csv`, `activity_df_CM_by_sample.csv`,
`loading_df_cell_subtype_by_CM.csv`, and all `w_df_*`/`h_df_*` display tables
must already contain `s_CM*`, `n_CM*`, or `t_CM*` labels. Downstream steps should
not need to guess whether these tables are raw or renamed.

Also write a small matrix-orientation readme/table so downstream plotting and
association code cannot silently transpose these matrices incorrectly.

After CM classification, write every final W/H matrix with canonical names. Do
not leave raw NMF component names such as `CM1`, `CM2`, or `joint_01` in final
analysis, association, or plotting inputs. At minimum, write:

```text
raw_to_canonical_CM_mapping.csv
W_df.csv = sample x canonical CM
H_df.csv = canonical CM x cell_subtype
loading_df_cell_subtype_by_CM.csv = cell_subtype x canonical CM
h_df_loading_cell_subtype_by_CM_raw.csv = canonical CM x cell_subtype
h_df_loading_cell_subtype_by_CM_zscore.csv = canonical CM x cell_subtype
h_df_loading_cell_subtype_by_CM_robust.csv = canonical CM x cell_subtype
h_df_loading_cell_subtype_by_CM_standard_scale_col.csv = canonical CM x cell_subtype
loading_df_cell_subtype_by_CM_fraction.csv = cell_subtype x canonical CM
```

Do not require downstream code to choose between raw and `_canonical`/`_renamed`
versions. The canonical table names above are the final table names and must
already be canonicalized. If a legacy or partial run has only raw files, repair
that run by producing the canonical table names before executing association or
plotting code.

Use this helper pattern as a guard in final plotting code before any CM loading
heatmap or top-subtype heatmap. In a correct canonical run, it should confirm
the tables are already canonical. It is a fallback guard for legacy/partial
inputs, not a substitute for early canonical naming. Items marked `FIXED` are
required behavior.

```python
def load_raw_to_canonical_map(mapping_csv: Path) -> dict[str, str]:
    mapping = pd.read_csv(mapping_csv)
    if not {"raw_component", "CM"}.issubset(mapping.columns):
        raise ValueError("raw_to_canonical_CM_mapping.csv must contain raw_component and CM")
    return dict(zip(mapping["raw_component"].astype(str), mapping["CM"].astype(str)))


def canonical_cm_order(classification_csv: Path) -> list[str]:
    classification = pd.read_csv(classification_csv)
    if "CM" not in classification.columns:
        raise ValueError("joint_module_classification.csv must contain CM")
    return classification["CM"].astype(str).tolist()


def as_cell_subtype_by_canonical_cm(
    df: pd.DataFrame,
    raw_to_canonical: dict[str, str],
    cm_order: list[str],
    table_name: str,
) -> pd.DataFrame:
    # FIXED: accept either CM x cell_subtype or cell_subtype x CM input.
    # FIXED: rename raw CM IDs on whichever axis contains them.
    out = df.copy()
    out.index = out.index.astype(str)
    out.columns = out.columns.astype(str)
    canonical_set = set(cm_order)
    raw_set = set(raw_to_canonical)

    index_has_cm = bool((set(out.index) & canonical_set) or (set(out.index) & raw_set))
    columns_has_cm = bool((set(out.columns) & canonical_set) or (set(out.columns) & raw_set))
    if index_has_cm and columns_has_cm:
        raise ValueError(f"{table_name}: both axes look like CM axes; inspect orientation")
    if not index_has_cm and not columns_has_cm:
        raise ValueError(f"{table_name}: no CM axis found; cannot plot CM heatmap")

    if index_has_cm:
        out = out.rename(index=raw_to_canonical).T
    else:
        out = out.rename(columns=raw_to_canonical)

    keep_cm = [cm for cm in cm_order if cm in out.columns]
    if not keep_cm:
        raise ValueError(
            f"{table_name}: no canonical CM columns after mapping. "
            "Check raw_to_canonical_CM_mapping.csv and matrix orientation."
        )
    out = out.loc[:, keep_cm]
    if out.shape[0] == 0 or out.shape[1] == 0:
        raise ValueError(f"{table_name}: empty heatmap matrix after orientation/mapping")
    return out.astype(float)
```

For top-subtype heatmaps, derive `cm_order` from `joint_module_classification.csv`
and derive row/subtype selection from the node/top-node table after converting
the loading matrix with `as_cell_subtype_by_canonical_cm`. Subset only after
canonicalization. If any requested CM from `top_nodes` is absent after mapping,
raise an explicit error listing missing CMs instead of plotting an empty array.

## CM Classification

Classify modules as:

```text
sharedCM = active in tumor and normal-like samples
normalCM = enriched in normal-like samples
tCM = enriched in tumor samples
```

Canonical CM naming:

```text
s_CM<global_number> = sharedCM
n_CM<global_number> = normalCM
t_CM<global_number> = tCM
```

Use this naming convention for final CM IDs in all saved W/H/activity/loading
tables, node tables, nodeplots, association matrices, scatter-source tables,
figure labels, and manifests. Do not use `joint_<number>_sharedCM`,
`joint_<number>_normalCM`, `joint_<number>_tCM`, or free-form CM names as final
IDs. Keep a mapping table from the raw NMF component name/order to the canonical
final CM ID so provenance is not lost.

The raw-to-canonical mapping must be applied before final plotting. In
particular, `balanced_joint_cm_reference_node_sets_after_edge_threshold.csv`,
`joint_cm_cell_subtype_nodes_*_from_H_df.csv`, nodeplot inputs, top-subtype
heatmap inputs, and all CM activity/loading heatmap inputs must use the same
canonical CM IDs. If a heatmap input has raw columns/index such as `CM1..CMK`
while `top_nodes` or classification has `s_CM*`, `n_CM*`, or `t_CM*`, that is a
bug: rename the CM axis with `raw_to_canonical_CM_mapping.csv` before subsetting.
Never let an empty intersection between CM names create an empty heatmap.

Assign numbers deterministically after classification by using one single global
CM order across all modules, then adding the class prefix. Do not restart
numbering within each prefix class. The global number is the component's
position in the global ordered list of all selected CMs, preferably the raw NMF
component order or the numeric raw component ID order. For example, if raw
components in global order are `CM1`, `CM2`, `CM3`, and their classes are
sharedCM, tCM, and normalCM, the final names are `s_CM1`, `t_CM2`, and `n_CM3`.
Do not rename them to `s_CM1`, `t_CM1`, and `n_CM1`.

Use the canonical activity-ratio thresholds unless the user explicitly changes
them:

```text
normal_specific_max_ratio = 0.5
tumor_specific_min_ratio = 2.0
min_active_fraction_for_specific = 0.05
force_shared_modules = ()
force_tcm_modules = ()
force_normalcm_modules = ()
```

Classify from CM activity differences, effect sizes, active fractions, and
statistical comparisons. Forced module-class overrides are empty by default and
must remain empty unless the user explicitly supplies overrides.

Keep an explicit classification table:

```text
joint_module_classification.csv
columns: raw_component, raw_component_order, CM, class, class_prefix, global_order, tumor_mean, normal_like_mean, delta, p, q
```

In `joint_module_classification.csv`, `global_order` must match the number in
the final CM ID. A `class_order` column may be included for diagnostics only, but
it must never drive the final CM number.

Use this explicit implementation pattern for canonical CM naming immediately
after raw NMF/NNLS activity is available. This code is intentionally included in
the markdown so downstream agents do not invent per-prefix numbering.

```python
def assign_canonical_cm_names(classification: pd.DataFrame) -> pd.DataFrame:
    """Assign final CM names with one global numeric order and class prefix.

    Required input columns:
      raw_component: raw NMF component ID, for example CM1, CM2, ...
      class: sharedCM, normalCM, or tCM

    Required output:
      CM: s_CM<global_order>, n_CM<global_order>, or t_CM<global_order>
    """
    out = classification.copy()
    prefix_map = {"sharedCM": "s", "normalCM": "n", "tCM": "t"}
    missing_classes = sorted(set(out["class"].astype(str)) - set(prefix_map))
    if missing_classes:
        raise ValueError(f"Unknown CM classes: {missing_classes}")

    if "raw_component_order" not in out.columns:
        # FIXED: global order is one order across all CMs, not per class.
        # Prefer raw NMF order already present in the table; otherwise preserve row order.
        out["raw_component_order"] = range(1, len(out) + 1)

    out = out.sort_values("raw_component_order", kind="stable").copy()
    out["global_order"] = range(1, len(out) + 1)
    out["class_prefix"] = out["class"].map(prefix_map)
    out["CM"] = out["class_prefix"] + "_CM" + out["global_order"].astype(str)

    # FIXED validation: the number in CM must equal global_order.
    expected = out["class_prefix"] + "_CM" + out["global_order"].astype(str)
    if not out["CM"].equals(expected):
        raise AssertionError("Canonical CM naming does not match global_order")
    return out


def apply_canonical_cm_names(
    W_raw: pd.DataFrame,
    H_raw: pd.DataFrame,
    classification: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return canonical W/H plus the raw-to-canonical mapping table."""
    cls = assign_canonical_cm_names(classification)
    mapping = cls[["raw_component", "raw_component_order", "CM", "class",
                   "class_prefix", "global_order"]].copy()
    raw_to_cm = dict(zip(mapping["raw_component"].astype(str), mapping["CM"].astype(str)))

    missing_w = sorted(set(mapping["raw_component"]) - set(map(str, W_raw.columns)))
    missing_h = sorted(set(mapping["raw_component"]) - set(map(str, H_raw.index)))
    if missing_w or missing_h:
        raise ValueError(f"Missing raw components in W/H: W={missing_w}, H={missing_h}")

    W = W_raw.rename(columns=raw_to_cm)
    H = H_raw.rename(index=raw_to_cm)
    cm_order = mapping.sort_values("global_order")["CM"].tolist()
    W = W.loc[:, cm_order]
    H = H.loc[cm_order, :]
    return W, H, mapping
```

After this block, write `W_df.csv`, `H_df.csv`,
`loading_df_cell_subtype_by_CM.csv`, `activity_df_sample_by_CM.csv`, all
`w_df_*`/`h_df_*` display tables, node tables, Epi-CM association tables, and
figures from the canonical `W` and `H` outputs only.

## CM Nodes and Edges

Define CM nodes as selected high-loading non-epithelial subtypes in `H`, then
apply the nodeplot/edge-threshold selection used by the canonical workflow.
The canonical nodeplot membership table is
`balanced_joint_cm_reference_node_sets_after_edge_threshold.csv`. This table is
the direct CM-to-subtype-node definition used by nodeplots and downstream CM
mapping. It is not the top-10 diagnostic table and is not the all-subtype
loading table. This canonical table must contain exactly three columns:
`CM`, `reference_node_rank`, and `node`. Do not add provenance/helper columns
such as `loading_rank_in_H_top20` to this table. If a run needs to preserve the
original H-loading rank, write a separate diagnostic table instead.

Nodeplot graph semantics are fixed: each per-CM nodeplot is a subtype-subtype
correlation graph within one CM's retained subtype node set. Graph nodes are
cell subtypes only. Do not draw a bipartite graph with a central CM node and
subtype leaves, do not draw CM-to-subtype loading edges as the nodeplot, and do
not replace the per-CM subtype-subtype network with a single combined all-CM
graph. The CM label belongs in the title/filename, not as a graph node.

Use these canonical node and edge parameters unless the user explicitly changes
them:

```text
top_n_subtypes = 20
plot_top_n_subtypes = 12
top_n_nodes = 10
edge_r_threshold = 0.25
```

Compute node-node correlations across tumor and normal-like samples. Record the selected correlation method and use a -1 to 1 scale for signed correlation heatmaps.

Edge origins:

```text
tumor-only
normal-like-only
shared
```

Canonical nodeplot code structure:

```text
1. Build candidate CM nodes from H/loadings in descending rank.
2. For each CM, take top_n_nodes candidate nodes for edge screening.
3. In normal-like samples only, compute a node x node Pearson correlation matrix
   from the column-minmax non-epithelial subtype frequency table.
4. In tumor samples only, compute the same node x node Pearson correlation matrix.
5. A node is retained for the final reference node set if it has at least one
   edge with Pearson r >= edge_r_threshold in either normal-like or tumor.
6. Save the union-retained nodes, preserving H/loading order, to
   balanced_joint_cm_reference_node_sets_after_edge_threshold.csv.
   The saved CSV must have exactly `CM, reference_node_rank, node`; no extra
   columns are allowed in the canonical nodeplot membership table.
7. Recompute or subset the normal-like and tumor node tables/correlation
   matrices using only the retained nodes.
8. Save status_specific_nodeplot_edges.csv with all retained-node pairs for
   both contexts, not only passing edges.
9. Draw standard context nodeplots from context-specific correlations.
10. Draw the tumor-centric edge-origin nodeplot from tumor-passing edges only,
    with edge origin determined by whether the same pair also passes in
    normal-like.
11. Draw top-node correlation heatmaps from the top10 diagnostic node table,
    without applying the edge filter, so the heatmap remains a diagnostic view.
```

Required `status_specific_nodeplot_edges.csv` columns:

```text
context
CM
node_a
node_b
pearson_r
edge_pass_r_ge_0.25
```

For this table, `edge_pass_r_ge_0.25` means `pearson_r >= edge_r_threshold` in
that row's own `context`. Do not use a tumor correlation value to decide a
normal-like edge pass, and do not use a normal-like correlation value to draw
the tumor edge strength.
Do not save this as a wide table with columns such as `node_i`, `node_j`,
`r_normal_like`, and `r_tumor` unless a second canonical long-format table with
the required columns above is also written and used for plotting.

Required standard nodeplot figures:

```text
normal_like_all_CM_nodeplot.pdf
normal_like_all_CM_nodeplot.svg
tumor_all_CM_nodeplot.pdf
tumor_all_CM_nodeplot.svg
normal_like_nodeplots_by_cm/<CM>_normal_like_nodeplot.pdf
normal_like_nodeplots_by_cm/<CM>_normal_like_nodeplot.svg
tumor_nodeplots_by_cm/<CM>_tumor_nodeplot.pdf
tumor_nodeplots_by_cm/<CM>_tumor_nodeplot.svg
nodeplot_network_node_legend.pdf
nodeplot_network_node_legend.svg
nodeplot_network_edge_colorbar.pdf
nodeplot_network_edge_colorbar.svg
```

For `normal_like_all_CM_nodeplot`, draw only normal-like edges with
`edge_pass_r_ge_0.25 == True`, and encode line color/width from the
normal-like `pearson_r`. For `tumor_all_CM_nodeplot`, draw only tumor edges with
`edge_pass_r_ge_0.25 == True`, and encode line color/width from the tumor
`pearson_r`.
The all-CM nodeplots are overview figures only and must be multi-panel grids:
one subplot per CM, each subplot containing that CM's own subtype-subtype
network. Do not draw all CM nodes/edges in one shared graph coordinate system.
For any all-CM overview, compute the subplot grid automatically from the number
of CMs using the closest-to-square factor/grid rule below; do not hard-code
`ncols=2`, `ncols=3`, or `ncols=4`. Also write one single-panel nodeplot per CM
under the `normal_like_nodeplots_by_cm/` and `tumor_nodeplots_by_cm/` child
directories. These per-CM figures must use the same node set, node order, node
colors, edge threshold, edge color scale, and edge width scale as the all-CM
overview. Do not replace per-CM nodeplots with a single combined all-CM canvas.

For every per-CM standard nodeplot, read exactly that CM's retained nodes from
`balanced_joint_cm_reference_node_sets_after_edge_threshold.csv`, subset
`status_specific_nodeplot_edges.csv` to the same `CM` and context, and draw only
edges whose `edge_pass_r_ge_0.25` is true in that context. If a CM has no passing
edges after filtering, write a manifest row explaining that no panel was drawn
instead of silently drawing an unrelated graph.

Required tumor-centric edge-origin nodeplot figures:

```text
tumor_centric_nodeplot_edge_origin.pdf
tumor_centric_nodeplot_edge_origin.svg
tumor_centric_nodeplots_by_cm/<CM>_tumor_centric_nodeplot_edge_origin.pdf
tumor_centric_nodeplots_by_cm/<CM>_tumor_centric_nodeplot_edge_origin.svg
tumor_centric_nodeplot_edge_origin_tumor_only_edge_colorbar.pdf
tumor_centric_nodeplot_edge_origin_tumor_only_edge_colorbar.svg
tumor_centric_nodeplot_edge_origin_shared_edge_colorbar.pdf
tumor_centric_nodeplot_edge_origin_shared_edge_colorbar.svg
tumor_centric_nodeplot_edge_origin_edge_class_legend.pdf
tumor_centric_nodeplot_edge_origin_edge_class_legend.svg
```

For `tumor_centric_nodeplot_edge_origin`, first select only edges that pass the
threshold in tumor samples. Every drawn line uses the tumor-sample Pearson
correlation value for color intensity only. Edge width is fixed for all drawn
edges so tumor-only versus shared class remains visually clear. Classify each drawn line as
`shared` if the same unordered node pair also passes the threshold in
normal-like samples for the same CM; otherwise classify it as `tumor_only`.
Use two clearly separated edge palettes for these classes. Do not draw
normal-like-only edges in this tumor-centric figure.
Every nodeplot figure must embed its corresponding edge-correlation colorbar
outside the graph body. This applies to every single-CM figure and every all-CM
overview. A normal-like nodeplot embeds one normal-like edge colorbar; a tumor
nodeplot embeds one tumor edge colorbar; a tumor-centric nodeplot embeds two
separate colorbars when both tumor-only and shared edges are actually present.
For each single-CM figure, include only bars for edge classes drawn in that
figure; do not show an unused palette bar. For an all-CM overview, use the union
of edge classes actually drawn across its panels. One shared set of bars is
sufficient for an all-CM multi-panel overview because every panel uses the same
mapping, but standalone legend files do not substitute for the embedded bars.
Every embedded or standalone edge colorbar and edge-class legend must be
generated from the exact colormap objects and `Normalize` object passed to the
nodeplot drawing function. Place embedded bars in a reserved bottom margin and
increase the canvas when needed so they never cover nodes, labels, panel titles,
or neighboring panels. Do not recreate a similar palette for the legend and do
not hard-code independent legend colors. If the nodeplot palette or correlation
range changes, every embedded and standalone legend/colorbar must change
automatically in the same plotting call. The class-legend line samples must come
from the corresponding live colormap; continuous colorbars must use
`ScalarMappable` with the same live colormap and normalization used for the
plotted edge colors.
Detect the available nodeplot contexts from the actual edge table before
building the mode list. A tumor-only cohort must generate only tumor nodeplots
with the tumor edge bar. A normal-like-only cohort must generate only
normal-like nodeplots with the normal-like edge bar. Generate tumor-centric
edge-origin nodeplots, their two-class bars, and their standalone legend files
only when both tumor and normal-like contexts exist. Never create an empty
absent-status nodeplot or an unused absent-status colorbar.
The tumor-centric overview figure is not sufficient by itself. Also write one
single-panel tumor-centric edge-origin nodeplot per CM under
`tumor_centric_nodeplots_by_cm/`, using the same edge-origin rules as the
overview. Within each CM panel, place only that CM's retained nodes. Use a
deterministic circular layout per CM, preserving retained node/loading rank
order. Color nodes by subtype prefix/lineage using a stable prefix palette, not
by CM ID. Draw tumor-only edges with the tumor-only palette and shared edges
with the shared palette; edge color intensity must use the tumor sample Pearson
r, but edge width must be the same fixed value for all drawn tumor-centric
edges. Do not draw all edges in gray, do not combine nodes from
different CMs into one single graph coordinate system, and do not silently
ignore edge-origin classes.
The tumor-centric overview must also be a closest-to-square multi-panel grid,
one subplot per CM, using the same per-CM edge-origin drawing logic.

Use this explicit nodeplot implementation pattern in the final plotting code.
Items marked `REPLACEABLE` may be changed by user request or project-specific
style. Items marked `FIXED` define the required behavior and must not be changed
without explicit user approval.

```python
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D

SEED = 42  # FIXED unless the user explicitly requests another seed.
EDGE_THRESHOLD = 0.25  # REPLACEABLE only if the analysis threshold changes.

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"


def save_pdf_svg(fig, stem: Path):
    stem.parent.mkdir(parents=True, exist_ok=True)
    pdf = stem.with_suffix(".pdf")
    svg = stem.with_suffix(".svg")
    fig.savefig(pdf, bbox_inches="tight", dpi=300)
    fig.savefig(svg, bbox_inches="tight", dpi=300)
    plt.close(fig)
    return pdf, svg


def node_name_column(nodes: pd.DataFrame) -> str:
    for col in ("node", "cell_subtype"):
        if col in nodes.columns:
            return col
    raise ValueError("Node table must contain either 'node' or 'cell_subtype'.")


def edge_key(node_a: str, node_b: str) -> tuple[str, str]:
    return tuple(sorted((str(node_a), str(node_b))))


def closest_square_grid(n_panels: int) -> tuple[int, int]:
    # FIXED: choose a grid closest to square for all-CM overview panels.
    # Prefer rows <= columns. Allow empty trailing panels only when n is prime
    # or when an exact factorization is much less square.
    if n_panels <= 0:
        return 0, 0
    best = None
    for n_rows in range(1, int(np.ceil(np.sqrt(n_panels))) + 1):
        n_cols = int(np.ceil(n_panels / n_rows))
        empty = n_rows * n_cols - n_panels
        aspect_gap = abs(n_cols - n_rows)
        candidate = (aspect_gap, empty, n_rows, n_cols)
        if best is None or candidate < best:
            best = candidate
    _, _, n_rows, n_cols = best
    return n_rows, n_cols


def prefix_color_map(nodes: list[str]) -> dict[str, tuple[float, float, float]]:
    # REPLACEABLE palette, but colors must be stable and prefix/lineage-based.
    prefixes = sorted({str(node).split("_")[0] for node in nodes})
    palette = plt.get_cmap("tab20").colors
    return {prefix: palette[i % len(palette)] for i, prefix in enumerate(prefixes)}


def prepare_nodeplot_inputs(nodes: pd.DataFrame, edges: pd.DataFrame):
    # FIXED expected node table: one row per retained node per CM.
    node_col = node_name_column(nodes)
    required_node_cols = {"CM", node_col}
    missing_nodes = required_node_cols - set(nodes.columns)
    if missing_nodes:
        raise ValueError(f"Missing node columns: {sorted(missing_nodes)}")

    # FIXED expected edge table: long format, one context per row.
    required_edge_cols = {"context", "CM", "node_a", "node_b", "pearson_r"}
    missing_edges = required_edge_cols - set(edges.columns)
    if missing_edges:
        raise ValueError(f"Missing edge columns: {sorted(missing_edges)}")

    if "edge_pass_r_ge_0.25" not in edges.columns:
        edges = edges.copy()
        edges["edge_pass_r_ge_0.25"] = edges["pearson_r"].astype(float) >= EDGE_THRESHOLD

    rank_cols = [c for c in ("rank", "reference_node_rank") if c in nodes.columns]
    sort_cols = ["CM"] + rank_cols + [node_col]
    nodes = nodes.sort_values(sort_cols).copy()
    return nodes, edges, node_col


def circular_rank_layout(cm_nodes: list[str]) -> dict[str, tuple[float, float]]:
    # FIXED: deterministic circular layout within each CM, preserving input rank.
    n = len(cm_nodes)
    return {
        node: (float(np.cos(2 * np.pi * i / max(n, 1))),
               float(np.sin(2 * np.pi * i / max(n, 1))))
        for i, node in enumerate(cm_nodes)
    }


def draw_one_cm_nodeplot(
    ax,
    cm: str,
    cm_nodes: list[str],
    cm_edges: pd.DataFrame,
    node_colors: dict[str, tuple[float, float, float]],
    edge_norm: Normalize,
    tumor_only_cmap,
    shared_cmap,
    mode: str,
):
    # FIXED mode values:
    #   "normal_like": draw normal-like passing edges only.
    #   "tumor": draw tumor passing edges only.
    #   "tumor_centric": draw tumor passing edges and color by edge_origin.
    graph = nx.Graph()
    graph.add_nodes_from(cm_nodes)
    pos = circular_rank_layout(cm_nodes)

    if mode == "tumor_centric":
        plot_edges = cm_edges.loc[
            cm_edges["context"].eq("tumor")
            & cm_edges["edge_pass_r_ge_0.25"].astype(bool)
        ].copy()
    elif mode == "tumor":
        plot_edges = cm_edges.loc[
            cm_edges["context"].eq("tumor")
            & cm_edges["edge_pass_r_ge_0.25"].astype(bool)
        ].copy()
    elif mode == "normal_like":
        plot_edges = cm_edges.loc[
            cm_edges["context"].eq("normal-like")
            & cm_edges["edge_pass_r_ge_0.25"].astype(bool)
        ].copy()
    else:
        raise ValueError(f"Unsupported nodeplot mode: {mode}")

    normal_pass = {
        edge_key(row.node_a, row.node_b)
        for _, row in cm_edges.loc[
            cm_edges["context"].eq("normal-like")
            & cm_edges["edge_pass_r_ge_0.25"].astype(bool)
        ].iterrows()
    }

    for _, row in plot_edges.iterrows():
        a, b = str(row.node_a), str(row.node_b)
        if a not in pos or b not in pos or a == b:
            continue
        if mode == "tumor_centric":
            edge_class = "shared" if edge_key(a, b) in normal_pass else "tumor_only"
        elif mode == "tumor":
            edge_class = "tumor"
        else:
            edge_class = "normal_like"
        graph.add_edge(a, b, weight=float(row.pearson_r), edge_class=edge_class)

    for edge_class, cmap in (
        ("tumor_only", tumor_only_cmap),
        ("shared", shared_cmap),
        ("tumor", tumor_only_cmap),
        ("normal_like", shared_cmap),
    ):
        class_edges = [(a, b, d) for a, b, d in graph.edges(data=True)
                       if d["edge_class"] == edge_class]
        if not class_edges:
            continue
        nx.draw_networkx_edges(
            graph,
            pos,
            edgelist=[(a, b) for a, b, _ in class_edges],
            edge_color=[cmap(edge_norm(d["weight"])) for _, _, d in class_edges],
            width=2.6 if mode == "tumor_centric" else [2.0 + 3.0 * edge_norm(d["weight"]) for _, _, d in class_edges],
            alpha=0.88,  # REPLACEABLE style only.
            ax=ax,
        )

    nx.draw_networkx_nodes(
        graph,
        pos,
        nodelist=cm_nodes,
        node_color=[node_colors.get(str(n).split("_")[0], "#999999") for n in cm_nodes],
        node_size=1700,  # REPLACEABLE style only.
        linewidths=0.8,
        edgecolors="white",
        ax=ax,
    )
    nx.draw_networkx_labels(
        graph,
        pos,
        labels={node: node for node in cm_nodes},
        font_size=8,  # REPLACEABLE style only.
        font_color="black",
        ax=ax,
    )
    ax.set_title(cm, fontsize=12, fontweight="bold")
    ax.set_aspect("equal")
    ax.axis("off")
    return {str(d["edge_class"]) for _, _, d in graph.edges(data=True)}


def add_edge_colorbars(fig, mode: str, edge_norm: Normalize,
                       tumor_only_cmap, shared_cmap,
                       present_edge_classes: set[str]):
    # FIXED: every saved nodeplot embeds bars from the exact live cmap/norm
    # objects used to color its edges. Include only edge classes actually drawn.
    # Bars occupy a reserved bottom margin.
    if mode == "normal_like" and "normal_like" in present_edge_classes:
        specs = [(shared_cmap, "Normal-like edge Pearson r")]
    elif mode == "tumor" and "tumor" in present_edge_classes:
        specs = [(tumor_only_cmap, "Tumor edge Pearson r")]
    elif mode == "tumor_centric":
        specs = []
        if "tumor_only" in present_edge_classes:
            specs.append((tumor_only_cmap, "Tumor-only edge Pearson r"))
        if "shared" in present_edge_classes:
            specs.append((shared_cmap, "Shared edge Pearson r"))
    elif mode not in {"normal_like", "tumor", "tumor_centric"}:
        raise ValueError(f"Unsupported nodeplot mode: {mode}")
    else:
        specs = []

    if not specs:
        return

    left_margin = 0.14
    right_margin = 0.08
    gap = 0.06 if len(specs) > 1 else 0.0
    bar_width = (1.0 - left_margin - right_margin - gap * (len(specs) - 1)) / len(specs)
    for bar_i, (cmap, label) in enumerate(specs):
        left = left_margin + bar_i * (bar_width + gap)
        cax = fig.add_axes([left, 0.055, bar_width, 0.022])
        sm = ScalarMappable(norm=edge_norm, cmap=cmap)
        sm.set_array([])
        cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
        cb.set_label(label, fontsize=8)
        cb.ax.tick_params(labelsize=7, length=2)


def save_edge_legends(out_dir: Path, edge_norm: Normalize, tumor_only_cmap, shared_cmap):
    # FIXED: use the exact live cmap/norm objects used by draw_one_cm_nodeplot.
    # Never duplicate or hard-code separate legend colors.
    for stem, cmap, label in (
        ("tumor_centric_nodeplot_edge_origin_tumor_only_edge_colorbar",
         tumor_only_cmap, "Tumor-only edge correlation in tumor samples (r)"),
        ("tumor_centric_nodeplot_edge_origin_shared_edge_colorbar",
         shared_cmap, "Shared edge correlation in tumor samples (r)"),
    ):
        fig, ax = plt.subplots(figsize=(4.0, 0.45))
        sm = ScalarMappable(norm=edge_norm, cmap=cmap)
        sm.set_array([])
        cb = fig.colorbar(sm, cax=ax, orientation="horizontal")
        cb.set_label(label, fontsize=10)
        cb.ax.tick_params(labelsize=8)
        save_pdf_svg(fig, out_dir / stem)

    # The categorical swatches identify palette families. Sample them from the
    # same live colormaps at the upper end of the same plotted r normalization,
    # so any nodeplot palette change is inherited automatically.
    legend_r = float(edge_norm.vmax)
    handles = [
        Line2D(
            [0], [0],
            color=tumor_only_cmap(edge_norm(legend_r)),
            lw=3,
            alpha=0.88,
            label="Tumor only",
        ),
        Line2D(
            [0], [0],
            color=shared_cmap(edge_norm(legend_r)),
            lw=3,
            alpha=0.88,
            label="Shared with normal-like",
        ),
    ]
    fig, ax = plt.subplots(figsize=(3.8, 1.1))
    ax.axis("off")
    ax.legend(handles=handles, loc="center", frameon=False, ncol=2)
    save_pdf_svg(fig, out_dir / "tumor_centric_nodeplot_edge_origin_edge_class_legend")


def plot_nodeplots_by_cm(nodes: pd.DataFrame, edges: pd.DataFrame, out_dir: Path):
    # FIXED outputs: per-CM figures and all-CM overviews are all required.
    nodes, edges, node_col = prepare_nodeplot_inputs(nodes, edges)
    out_dir = Path(out_dir)
    all_nodes = nodes[node_col].astype(str).tolist()
    node_colors = prefix_color_map(all_nodes)
    tumor_only_cmap = LinearSegmentedColormap.from_list(
        "tumor_only_blue", ["#dbeafe", "#1d4ed8"]  # REPLACEABLE color palette.
    )
    shared_cmap = LinearSegmentedColormap.from_list(
        "shared_red", ["#fee2e2", "#b91c1c"]  # REPLACEABLE color palette.
    )
    edge_norm = Normalize(vmin=EDGE_THRESHOLD, vmax=1.0)

    available_contexts = set(edges["context"].dropna().astype(str))
    modes = {}
    if "normal-like" in available_contexts:
        modes["normal_like"] = (
            "normal_like_nodeplots_by_cm",
            "_normal_like_nodeplot",
        )
    if "tumor" in available_contexts:
        modes["tumor"] = ("tumor_nodeplots_by_cm", "_tumor_nodeplot")
    if {"tumor", "normal-like"}.issubset(available_contexts):
        modes["tumor_centric"] = (
            "tumor_centric_nodeplots_by_cm",
            "_tumor_centric_nodeplot_edge_origin",
        )
    if not modes:
        raise ValueError(
            "No supported nodeplot context found; expected 'tumor' and/or 'normal-like'."
        )

    written = []
    cm_order = nodes["CM"].drop_duplicates().tolist()
    for cm in cm_order:
        cm_nodes = nodes.loc[nodes["CM"].eq(cm), node_col].astype(str).tolist()
        cm_edges = edges.loc[edges["CM"].eq(cm)].copy()
        for mode, (subdir, suffix) in modes.items():
            canvas = (6.8, 4.8) if mode == "tumor_centric" else (4.8, 4.8)
            fig, ax = plt.subplots(figsize=canvas)  # REPLACEABLE canvas only.
            present_edge_classes = draw_one_cm_nodeplot(
                ax=ax,
                cm=cm,
                cm_nodes=cm_nodes,
                cm_edges=cm_edges,
                node_colors=node_colors,
                edge_norm=edge_norm,
                tumor_only_cmap=tumor_only_cmap,
                shared_cmap=shared_cmap,
                mode=mode,
            )
            fig.tight_layout(rect=[0, 0.16, 1, 1])
            add_edge_colorbars(
                fig,
                mode=mode,
                edge_norm=edge_norm,
                tumor_only_cmap=tumor_only_cmap,
                shared_cmap=shared_cmap,
                present_edge_classes=present_edge_classes,
            )
            written.extend(save_pdf_svg(fig, out_dir / subdir / f"{cm}{suffix}"))

    # FIXED: all-CM overviews are closest-to-square grids with one CM per panel.
    overview_stem_by_mode = {
        "normal_like": "normal_like_all_CM_nodeplot",
        "tumor": "tumor_all_CM_nodeplot",
        "tumor_centric": "tumor_centric_nodeplot_edge_origin",
    }
    overview_specs = {mode: overview_stem_by_mode[mode] for mode in modes}
    n_rows, n_cols = closest_square_grid(len(cm_order))
    for mode, stem in overview_specs.items():
        if n_rows == 0:
            continue
        fig, axes = plt.subplots(
            n_rows,
            n_cols,
            figsize=(n_cols * 4.4, n_rows * 4.2 + 0.9),  # Includes bar margin.
            squeeze=False,
        )
        flat_axes = axes.ravel()
        overview_edge_classes = set()
        for ax, cm in zip(flat_axes, cm_order):
            cm_nodes = nodes.loc[nodes["CM"].eq(cm), node_col].astype(str).tolist()
            cm_edges = edges.loc[edges["CM"].eq(cm)].copy()
            panel_edge_classes = draw_one_cm_nodeplot(
                ax=ax,
                cm=cm,
                cm_nodes=cm_nodes,
                cm_edges=cm_edges,
                node_colors=node_colors,
                edge_norm=edge_norm,
                tumor_only_cmap=tumor_only_cmap,
                shared_cmap=shared_cmap,
                mode=mode,
            )
            overview_edge_classes.update(panel_edge_classes)
        for ax in flat_axes[len(cm_order):]:
            ax.axis("off")
        fig.tight_layout(rect=[0, 0.11, 1, 1])
        add_edge_colorbars(
            fig,
            mode=mode,
            edge_norm=edge_norm,
            tumor_only_cmap=tumor_only_cmap,
            shared_cmap=shared_cmap,
            present_edge_classes=overview_edge_classes,
        )
        written.extend(save_pdf_svg(fig, out_dir / stem))

    if "tumor_centric" in modes:
        save_edge_legends(out_dir, edge_norm, tumor_only_cmap, shared_cmap)
    return written


# Required call in the final plotting script:
# nodes = pd.read_csv("tables/03-cm-classification-nodes-edges/tumor_network_nodes_from_H_df.csv")
# edges = pd.read_csv("tables/03-cm-classification-nodes-edges/status_specific_nodeplot_edges.csv")
# plot_nodeplots_by_cm(nodes, edges, Path("figures/02-cm-lineage-final-plotting/nodeplots"))
```

Required top-node diagnostic correlation heatmap figures:

```text
normal_like_top10_node_correlation_heatmap_no_edge_filter.pdf
normal_like_top10_node_correlation_heatmap_no_edge_filter.svg
tumor_top10_node_correlation_heatmap_no_edge_filter.pdf
tumor_top10_node_correlation_heatmap_no_edge_filter.svg
top10_node_correlation_heatmaps_no_edge_filter_by_cm/normal_like/<CM>_normal_like_top10_node_correlation_heatmap_no_edge_filter.pdf
top10_node_correlation_heatmaps_no_edge_filter_by_cm/normal_like/<CM>_normal_like_top10_node_correlation_heatmap_no_edge_filter.svg
top10_node_correlation_heatmaps_no_edge_filter_by_cm/tumor/<CM>_tumor_top10_node_correlation_heatmap_no_edge_filter.pdf
top10_node_correlation_heatmaps_no_edge_filter_by_cm/tumor/<CM>_tumor_top10_node_correlation_heatmap_no_edge_filter.svg
```

These heatmaps use the top10 loading diagnostic table directly, not the
edge-filtered membership table. For each CM, select exactly the top 10 subtypes
by H/loading rank from the diagnostic top-node table, even if some of those
subtypes are absent from the edge-threshold nodeplot. Use Pearson correlation,
identical row/column order, and a fixed centered -1 to 1 diverging scale. Do not
write numeric correlation values inside heatmap cells. The default top10
correlation heatmap has no cell text annotations; if significance labels are
explicitly requested, use `ns/*/**/***` symbols only, never numeric values.

The summary heatmap files are fixed multi-panel figures, not one merged
all-node matrix. In each status-specific summary (`tumor` or `normal_like`),
draw one separate square panel per CM. Each panel must contain only that CM's
own top10 H/loading diagnostic nodes and their within-status 10 x 10 Pearson
correlation matrix. Keep CM panels in canonical CM order and give all panels the
same `RdBu_r`, `vmin=-1`, `center=0`, `vmax=1` scale with one shared colorbar.
Never take the union of nodes from different CMs to draw one giant correlation
matrix, and never add cross-CM node correlations to these diagnostic figures.
In addition to the status-specific multi-panel summary, save every CM panel as
its own PDF and SVG under
`top10_node_correlation_heatmaps_no_edge_filter_by_cm/<status>/`. When both
tumor and normal-like samples exist, compute correlations independently within
each status, write separate `tumor/` and `normal_like/` panel directories, and
write both status-specific summary figures. Never pool tumor and normal-like
samples for these correlation matrices and never let one status overwrite the
other. For tumor-only or normal-only inputs, write only the corresponding
status outputs and explicitly skip the absent-status heatmaps.

Canonical top10 node-correlation heatmap code:

```python
def top10_nodes_by_cm(top_nodes: pd.DataFrame, cm: str, node_col: str = "cell_subtype") -> list[str]:
    # FIXED: use top10 H/loading diagnostic nodes, not edge-filtered nodeplot nodes.
    rank_col = "rank" if "rank" in top_nodes.columns else "loading_rank"
    sub = top_nodes.loc[top_nodes["CM"].astype(str).eq(str(cm))].copy()
    sub = sub.sort_values(rank_col, kind="stable").head(10)
    if node_col not in sub.columns and "node" in sub.columns:
        node_col = "node"
    nodes = sub[node_col].astype(str).tolist()
    if len(nodes) == 0:
        raise ValueError(f"No top10 diagnostic nodes found for {cm}")
    return nodes


def plot_top10_node_correlation_heatmap(
    corr: pd.DataFrame,
    top_nodes: pd.DataFrame,
    cm: str,
    stem: Path,
    title: str,
):
    nodes = top10_nodes_by_cm(top_nodes, cm)
    nodes = [n for n in nodes if n in corr.index and n in corr.columns]
    if len(nodes) < 2:
        raise ValueError(f"{cm}: fewer than 2 top10 nodes present in correlation matrix")
    plot_df = corr.loc[nodes, nodes].astype(float)
    fig, ax = plt.subplots(figsize=(max(4.2, 0.48 * len(nodes)), max(4.0, 0.48 * len(nodes))))
    sns.heatmap(
        plot_df,
        cmap="RdBu_r",
        center=0,
        vmin=-1,
        vmax=1,
        square=True,
        linewidths=0.2,
        linecolor="white",
        annot=False,  # FIXED: no numeric values in cells.
        cbar_kws={"label": "Pearson r"},
        ax=ax,
    )
    ax.set_title(title)
    ax.set_xticklabels(nodes, rotation=90, ha="center", fontsize=7)
    ax.set_yticklabels(nodes, rotation=0, fontsize=7)
    ax.tick_params(axis="both", length=0)
    fig.tight_layout()
    save_pdf_svg(fig, stem)


def plot_top10_node_correlation_summary(
    corr: pd.DataFrame,
    top_nodes: pd.DataFrame,
    cm_order: list[str],
    stem: Path,
    context_label: str,
    ncols: int = 4,
):
    # FIXED: one CM-specific matrix per panel; never build a union-node matrix.
    panel_data = []
    for cm in cm_order:
        nodes = top10_nodes_by_cm(top_nodes, cm)
        nodes = [n for n in nodes if n in corr.index and n in corr.columns]
        if len(nodes) < 2:
            raise ValueError(f"{cm}: fewer than 2 top10 nodes present in correlation matrix")
        panel_data.append((cm, nodes, corr.loc[nodes, nodes].astype(float)))

    if not panel_data:
        raise ValueError(f"No CM panels available for {context_label}")

    ncols = min(ncols, len(panel_data))
    nrows = int(np.ceil(len(panel_data) / ncols))
    panel_size = 4.0
    fig = plt.figure(figsize=(panel_size * ncols + 0.8, panel_size * nrows))
    grid = fig.add_gridspec(
        nrows,
        ncols + 1,
        width_ratios=[1.0] * ncols + [0.05],
        wspace=0.65,
        hspace=0.75,
    )
    cbar_ax = fig.add_subplot(grid[:, -1])

    for panel_i, (cm, nodes, plot_df) in enumerate(panel_data):
        row, col = divmod(panel_i, ncols)
        ax = fig.add_subplot(grid[row, col])
        sns.heatmap(
            plot_df,
            cmap="RdBu_r",
            center=0,
            vmin=-1,
            vmax=1,
            square=True,
            linewidths=0.2,
            linecolor="white",
            annot=False,
            cbar=(panel_i == 0),
            cbar_ax=cbar_ax if panel_i == 0 else None,
            cbar_kws={"label": "Pearson r"} if panel_i == 0 else None,
            ax=ax,
        )
        ax.set_title(str(cm))
        ax.set_xticklabels(nodes, rotation=90, ha="center", fontsize=6.5)
        ax.set_yticklabels(nodes, rotation=0, fontsize=6.5)
        ax.tick_params(axis="both", length=0)

    for empty_i in range(len(panel_data), nrows * ncols):
        row, col = divmod(empty_i, ncols)
        ax = fig.add_subplot(grid[row, col])
        ax.set_axis_off()

    fig.suptitle(f"{context_label}: top10 node correlations by CM", y=0.995)
    fig.subplots_adjust(top=0.94, bottom=0.08, left=0.06, right=0.95)
    save_pdf_svg(fig, stem)
```

Minimum node table columns:

```text
CM
cell_subtype
loading
rank
cell_lineage
status_context
```

Required nodeplot membership table:

```text
balanced_joint_cm_reference_node_sets_after_edge_threshold.csv
columns exactly: CM, reference_node_rank, node
one row per selected subtype node used in the CM nodeplot/reference CM definition
```

Do not substitute `joint_cm_cell_subtype_nodes_top10_from_H_df.csv`,
`joint_cm_cell_subtype_nodes_top20_from_H_df.csv`, or
`joint_cm_cell_subtype_nodes_all_from_H_df.csv` for this membership table. Those
tables are loading diagnostics. The nodeplot membership table should contain the
final selected node set per CM after the canonical node/edge selection, and may
have fewer or more than exactly 10 rows for a CM depending on the selection.
Do not include `loading_rank_in_H_top20`, raw loading values, lineage labels, or
other diagnostic columns in `balanced_joint_cm_reference_node_sets_after_edge_threshold.csv`.
If these values are useful, save them to a separate diagnostic CSV whose name
does not replace the canonical membership table.

## Epi-CM Association

Correlate epithelial subtype abundance with CM activity across samples.

This step can run only after canonical CM activity has been produced or supplied:
the required CM-side input is `activity_df_sample_by_CM.csv` or an equivalent
sample x CM activity matrix from the balanced joint NMF W matrix. Do not use
`non_epi_subtype_frequency.csv`, raw non-epithelial subtype frequencies, or raw
non-epithelial subtype columns as a substitute for CM activity. If W/CM activity
is missing, run or request the balanced joint NMF step first; do not relabel
non-epithelial subtypes as CMs.

Epithelial-cell-count eligibility is applied here, not during balanced joint
NMF. Read `sample_inclusion_exclusion.csv` and restrict association samples to
`keep_for_epi_cm == True` if that column exists; otherwise use the documented
epithelial threshold at this step and record it. The default threshold is
`min_epi_cells_per_sample = 1` unless the user supplies a stricter project
threshold. Do not refit K, W, H, CM classification, or nodeplots after applying
the epithelial association filter.

Run all pairwise combinations. For each method branch, compute the Cartesian
product of every epithelial subtype column and every canonical CM activity
column in the eligible sample set. Do not restrict the analysis to preselected,
significant, top-ranked, plotted, or claim-support pairs before computing the
correlation and q-value tables. Pair selection is allowed only after all-pair
tables have been written, for display subsets or manuscript-focused panels.

Keep Spearman and Pearson as separate analysis branches, not variants inside the
same output directory or figure set. The two branches may read the same
`epi_subtype_frequency.csv`, `activity_df_sample_by_CM.csv`, and
`sample_status.csv`, but they must have separate code files, table directories,
figure directories, and parameter/provenance records. Do not mix Spearman and
Pearson matrices in one heatmap/scatter directory, and do not reuse one method's
q values, stars, labels, or file names for the other method.

When both Spearman and Pearson are requested, their output inventory must be
parallel and complete. Spearman and Pearson must generate the same classes of
outputs, the same status contexts, the same epithelial subtype x CM pair set,
the same heatmap families, the same scatterplot families, and the same manifest
schema. The only expected differences are the statistical method, statistic
column names (`rho` for Spearman, `r` for Pearson), p/q values, method labels,
and method-specific file or directory stems. Do not let Pearson produce fewer
tables or figures than Spearman, and do not let Spearman produce fewer tables or
figures than Pearson, unless the user explicitly requests one method only.

Default branch:

```text
Spearman correlation
```

Optional separate branch:

```text
Pearson correlation
```

If Pearson is requested, create a Pearson-specific branch and recompute
correlation, p values, and FDR q values from the same input matrices. Use `r`
naming for Pearson and `rho` naming for Spearman. Do not rename Spearman
matrices as Pearson, and do not put Pearson outputs into the Spearman branch.
Use the same output structure and file inventory as the Spearman branch, with
only method-specific names changed.

Recommended secondary-task names:

```text
epi-cm-association-spearman
epi-cm-association-pearson
```

For each pair:

```text
association = corr(epithelial subtype abundance, CM activity)
p_value = method-specific correlation test p value
q_value = FDR-adjusted p value across all epithelial subtype x CM pairs in that analysis branch and sample context
```

If status-specific associations are part of the workflow, repeat the full
epithelial subtype x CM Cartesian product separately within each status context,
such as tumor and normal-like. Do not replace per-subtype associations with a
single total epithelial fraction unless the user explicitly requests that
alternative branch.

Significance labels:

```text
Use the thresholds from the plotting code/config or the user-provided plan.
State the active q-star thresholds in the figure legend, caption, or run report.
```

## Outputs

Expected tables:

```text
non_epi_subtype_frequency.csv
non_epi_subtype_frequency_column_minmax.csv
column_minmax_params.csv
group_balanced_sample_weights.csv
epi_subtype_frequency.csv
sample_status.csv
W_df.csv
H_df.csv
activity_df_sample_by_CM.csv
loading_df_cell_subtype_by_CM.csv
joint_CM_activity_tumor_vs_normal_mean_sd_summary.csv
joint_module_classification.csv
status_specific_nodeplot_edges.csv
balanced_joint_cm_reference_node_sets_after_edge_threshold.csv
tumor_network_nodes_from_H_df.csv
normal_like_network_nodes_from_H_df.csv
balanced_joint_cm_epi_cm_association_tumor_rho_matrix.csv
balanced_joint_cm_epi_cm_association_tumor_q_matrix.csv
```

Spearman association branch tables must live under a Spearman-specific task
directory and use Spearman-specific labels in plot titles/captions. Acceptable
Spearman stems may use `rho_matrix` naming or explicit
`spearman_rho_matrix` naming, but the branch directory must state `spearman`.
The Spearman branch table/figure manifest should be the template for the
Pearson branch inventory when Pearson is requested.

For Pearson branch tables:

```text
balanced_joint_cm_epi_cm_association_tumor_pearson_r_matrix.csv
balanced_joint_cm_epi_cm_association_tumor_pearson_q_matrix.csv
balanced_joint_cm_epi_cm_association_normal-like_pearson_r_matrix.csv
balanced_joint_cm_epi_cm_association_normal-like_pearson_q_matrix.csv
```

Pearson association branch tables must live under a Pearson-specific task
directory and use Pearson-specific labels in plot titles/captions.
The Pearson branch must mirror the Spearman branch output inventory one-to-one:
same contexts, same pair count, same matrix dimensions, same scatterplot count,
same heatmap count, and same manifest columns, with only method-specific stems
and statistic labels changed.

## Plots

Produce:

```text
CM nodeplot
CM activity heatmap
CM loading heatmap
Epi-CM association heatmap with q-value stars
all-pair Epi-CM scatter plots
NMF K-selection diagnostics
top-node correlation heatmaps
```

Plot requirements:

- Correlation heatmaps use a centered diverging palette.
- Epi-CM heatmaps must be generated inside the matching method branch: Spearman
  heatmaps from Spearman rho/q matrices only, Pearson heatmaps from Pearson r/q
  matrices only.
- Scatter plots must be generated inside the matching method branch and label the
  selected correlation method. Generate every epithelial subtype x CM pair for
  each status context, not only highlighted pairs. Do not place Pearson scatter
  plots in the Spearman-only scatter output directory or vice versa. Use the
  `adata_epi`/lineage `cell_subtype` color mapping for epithelial subtype point
  colors and validate palette coverage before plotting.
- Nodeplot edge style/color should distinguish tumor-only, normal-like-only, and shared edges when used.
- Node-node heatmaps use Pearson correlation on a -1 to 1 scale unless the user explicitly requests another method.

## Reference Code Snippets

These are canonical implementation snippets distilled from the project analysis
code. Treat `FIXED` behavior as required; only change `REPLACEABLE` paths,
column names, and user-specified parameters.

Build compartment frequency matrices:

```python
import pandas as pd

def build_frequency_tables(
    obs,
    epithelial_types=("Epithelial Cells",),
    min_non_epi_cells_per_sample=50,
    min_epi_cells_per_sample=1,  # REPLACEABLE; used only for Epi-CM association.
):
    required = {"sample_id", "status", "cell_type", "cell_subtype"}
    missing = required - set(obs.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")

    sample_status = (
        obs[["sample_id", "status"]]
        .drop_duplicates()
        .set_index("sample_id")
        .sort_index()
    )
    epi_mask = obs["cell_type"].isin(epithelial_types)

    def counts_and_freq(df):
        counts = pd.crosstab(df["sample_id"], df["cell_subtype"]).sort_index()
        freq = counts.div(counts.sum(axis=1), axis=0).fillna(0)
        return counts, freq

    non_epi_counts, non_epi_freq = counts_and_freq(obs.loc[~epi_mask])
    epi_counts, epi_freq = counts_and_freq(obs.loc[epi_mask])

    compartment_counts = pd.DataFrame(index=sample_status.index)
    compartment_counts["non_epi_cells"] = non_epi_counts.sum(axis=1).reindex(compartment_counts.index).fillna(0).astype(int)
    compartment_counts["epi_cells"] = epi_counts.sum(axis=1).reindex(compartment_counts.index).fillna(0).astype(int)

    inclusion = compartment_counts.copy()
    inclusion["has_status"] = inclusion.index.isin(sample_status.index)
    inclusion["has_non_epi"] = inclusion["non_epi_cells"] >= min_non_epi_cells_per_sample
    inclusion["has_epi_for_association"] = inclusion["epi_cells"] >= min_epi_cells_per_sample
    # FIXED: CM/NMF uses non-epithelial eligibility only.
    inclusion["keep_for_cm"] = inclusion["has_status"] & inclusion["has_non_epi"]
    # FIXED: epithelial filtering starts only for Epi-CM association/figures.
    inclusion["keep_for_epi_cm"] = inclusion["keep_for_cm"] & inclusion["has_epi_for_association"]

    def exclusion_reason(row):
        reasons = []
        if not row["has_status"]:
            reasons.append("missing_status")
        if not row["has_non_epi"]:
            reasons.append("non_epi_cells_below_threshold")
        if not row["has_epi_for_association"]:
            reasons.append("epi_cells_below_association_threshold")
        return ";".join(reasons) if reasons else "kept"

    inclusion["exclusion_reason"] = inclusion.apply(exclusion_reason, axis=1)
    keep_cm_samples = inclusion.index[inclusion["keep_for_cm"]].intersection(non_epi_counts.index).sort_values()
    keep_epi_cm_samples = inclusion.index[inclusion["keep_for_epi_cm"]].intersection(non_epi_counts.index).sort_values()
    if len(keep_cm_samples) == 0:
        raise ValueError("No samples pass CM sample eligibility filters.")

    return (
        sample_status.loc[keep_cm_samples],
        non_epi_counts.reindex(keep_cm_samples).fillna(0).astype(int),
        non_epi_freq.reindex(keep_cm_samples).fillna(0),
        epi_counts.reindex(keep_epi_cm_samples).fillna(0).astype(int),
        epi_freq.reindex(keep_epi_cm_samples).fillna(0),
        compartment_counts,
        inclusion,
    )
```

Compute Epi-CM correlations:

```python
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr
from statsmodels.stats.multitest import multipletests

def epi_cm_correlations(epi_freq, cm_activity, method="spearman"):
    common = epi_freq.index.intersection(cm_activity.index)
    E = epi_freq.loc[common].astype(float)
    C = cm_activity.loc[common].astype(float)
    corr = pd.DataFrame(index=E.columns, columns=C.columns, dtype=float)
    pval = pd.DataFrame(index=E.columns, columns=C.columns, dtype=float)
    test = spearmanr if method == "spearman" else pearsonr

    # Required: exhaustive epithelial subtype x CM Cartesian product.
    # Do not pre-filter to selected, significant, top-ranked, or plotted pairs.
    for epi in E.columns:
        for cm in C.columns:
            x, y = E[epi], C[cm]
            ok = x.notna() & y.notna()
            if ok.sum() < 3 or x[ok].nunique() < 2 or y[ok].nunique() < 2:
                corr.loc[epi, cm] = np.nan
                pval.loc[epi, cm] = np.nan
                continue
            r, p = test(x[ok], y[ok])
            corr.loc[epi, cm] = r
            pval.loc[epi, cm] = p

    flat = pval.to_numpy().ravel()
    valid = np.isfinite(flat)
    qflat = np.full(flat.shape, np.nan, dtype=float)
    qflat[valid] = multipletests(flat[valid], method="fdr_bh")[1]
    qval = pd.DataFrame(qflat.reshape(pval.shape), index=pval.index, columns=pval.columns)
    return corr, pval, qval
```

## Validation

- Confirm matrix orientation.
- Confirm sample IDs match across all tables.
- Confirm both tumor and normal-like groups exist for balanced joint analysis.
- Confirm column-minmax was performed by column alignment, not whole-matrix
  scaling and not `DataFrame.where()` with a misaligned column Series mask.
- Confirm the saved column-minmax matrix has at least one non-zero-variance
  non-epithelial subtype before NMF.
- Confirm `W = samples x K` and `H = K x non-epithelial subtypes`.
- Confirm `K` was selected from a non-empty K-selection table and that a small K
  is not caused by an empty/all-zero NMF input matrix.
- Confirm plot labels match the selected correlation method.
- Confirm q values are recomputed for the selected correlation method.
- When the optional Pearson branch is requested, confirm Spearman and Pearson
  outputs are in separate task directories and no plot or table combines the
  two methods unless the user explicitly requested a cross-method comparison.
- Confirm epithelial and non-epithelial subtype sets are not accidentally mixed.
- Do not use old HCM/reference-guided variants unless explicitly asked.

## 04-Spatial Validation Optional

This block is the former Module 05, copied as the canonical instruction source
for the compact workflow. Run it only when the user provides spatial inputs or
explicitly requests spatial validation. Inside this big skill, write this block's
outputs under
`epi-cm-core-workflow/{codes,h5ad,tables,figures}/04-spatial-validation-optional/`
unless the user explicitly asks to use the original numbered module output tree.
Do not use the shorter summary as a substitute for these copied rules.

# 05-project-spatial-validation

Use this skill for spatial validation of CM-lineage relationships. The default workflow is Tangram-based pseudobulk selected-lineage-CM spatial mapping. If a example branch exists, treat its selected lineage as an example only. For any project, set the lineage, spot-score features, marker profiles, status groups, and biological interpretation according to the user goal. The biological goal is to validate spatial colocalization or exclusion between selected lineage states and CM programs.

This module is not isolated from the single-cell, epithelial-state, or CM-Epi modules. Reuse the same epithelial subtype names, CM names, marker gene definitions, sample IDs, and subtype-prefix mapping. If spatial sample names differ from single-cell or CM sample names, write a mapping table before joining or plotting.

## Pre-Execution Plan Requirement

Before executing code from this skill, write a concise method-and-result plan that the user can review and copy as the goal. Keep it result-oriented rather than overly procedural. Include only:

```text
analysis goal / expected result
method route to use
main inputs or provided intermediates
major code modules to run or skip
expected output figures/tables
key validation criterion
```

Do not start long-running analysis, dependency installation, or file-rewriting steps until this short plan has been stated. For simple inspection-only tasks, one or two sentences are enough.

If the user does not provide a manual choice for parameters, thresholds, method options, output naming, or optional branches, use the documented default settings in this skill and state that the default was used.

## Code Inclusion Contract

For every fragile analysis or final-figure workflow in this module, keep the
canonical executable logic as explicit code in this `SKILL.md` or in a named
module-local `.py`, `.R`, or `.ipynb` file generated from the code block in this
skill. Do not rely only on prose such as "follow the old notebook" or "use the
same script". Code blocks must mark user/project-adjustable values with
`REPLACEABLE` comments and non-negotiable behavior with `FIXED` comments. When a
workflow has a final redrawn figure style, encode that final plotting code
directly and do not keep a second competing plotting implementation for the same
panel.

## Canonical Plotting Contract

For Module 03 and all later modules, the final plotting code and style rules written in this skill are the source of truth. Use the code, parameters, palettes, layout rules, and output names specified here or in module-local code files generated from this skill. Do not tell a future agent to inspect old notebooks or external project paths during execution. If both an analysis route and a later redraw/plot route are encoded for the same figure, the later redraw/plot implementation is canonical. Do not improvise alternate plot types, palettes, layouts, statistical labels, file formats, or single-pair shortcuts. If a required final figure lacks explicit plotting code or style rules in this skill, stop and ask for the skill to be updated before running; do not invent a new plotting route.

## Project Organization and Figure Output Contract

Treat each numbered module folder as its own output boundary with one shared four-directory layout: `figures/`, `tables/`, `codes/`, and `h5ad/`. Module output directory names must use the active project slug, for example `05-<project_slug>-spatial-validation/`; for BRCA use `05-brca-spatial-validation/`. The four top-level category directories may be created at module setup. Secondary-module/task directories inside those category directories must be created only when that category will receive at least one real output for that task. Do not pre-create empty task directories under `figures/`, `tables/`, `codes/`, or `h5ad/` just to mirror the layout. A directory creation command for a secondary task or candidate must be coupled to writing a real output file there; if the output is not generated, do not leave that task directory behind. Secondary modules are logical analysis units inside that numbered module, but they do not require a directory under every category. This is a write-location rule, not a read restriction: a module may read/reuse files and already generated outputs from other modules as inputs, but newly generated outputs for the current module must be written inside the current numbered module. Across the whole skill workflow, an agent must not delete, clear, overwrite, or move any existing output file or directory anywhere unless the user explicitly names the exact path and operation. During normal module execution, do not write, move, overwrite, clear, or delete files or directories inside any other numbered module output directory. Also do not delete, clear, overwrite, or move any existing output file or directory inside the current module unless the user explicitly names the exact path and operation. If a new result would conflict with an existing output, write to a new versioned path or stop and ask. Deleting any output directory is never part of a module run; it requires a separate explicit cleanup request naming the exact path. For example, Module 02 subtype outputs go under `02-<project_slug>-cell-annotation/{figures,tables,codes,h5ad}/01-epithelial-subclustering/`, not under `01-<project_slug>-singlecell-integration/`, not directly as files under `02-<project_slug>-cell-annotation/`, and not under the top-level skill source directory.

Use stable numbered secondary-module/task names that describe the analysis step, lineage, method, or figure group. Reuse the same secondary-module/task name only under category directories that actually receive outputs from that analysis, so files stay aligned without creating empty placeholder task directories. For example:

```text
02-<project_slug>-cell-annotation/
  codes/
    01-epithelial-subclustering/
      02_epithelial_subclustering.ipynb
      run_epi_subclustering.py
  h5ad/
    01-epithelial-subclustering/
      adata_epi.h5ad
  figures/
    01-epithelial-subclustering/
      umap_cell_subtype.pdf
  tables/
    01-epithelial-subclustering/
      epithelial_subtype_counts.csv
```

By default, save executable/reproducibility code under the current module's shared `codes/<secondary-module>/`, using ordered names such as `01_read_merge.ipynb`, `02_qc.py`, or `03_integrate.R`. Save AnnData-like objects under `h5ad/<secondary-module>/` as `.h5ad`, `.loom`, `.rds`, or equivalent files with stable names. Save corresponding figure files under `figures/<secondary-module>/` and use ordered names such as `01_umap.pdf` or `02_marker_dotplot.svg`. Save text-like and tabular outputs under `tables/<secondary-module>/`, such as CSV/TSV/XLSX/TXT/JSON/YAML logs, manifests, reports, mapping files, and parameter records. `figures/` should contain figure files only. `tables/` should contain text-like and tabular outputs only. `codes/` should contain executable/reproducibility code only. `h5ad/` should contain AnnData-like/intermediate object files only. Add `tables/<secondary-module>/readme.txt` documenting the input files, including any cross-module input/output files that were read, code order, h5ad/loom/rds objects, output figures/tables, and any skipped optional branches.

Do not write new h5ad, code, figures, or tables directly into the numbered module root. The numbered module root may contain the module `SKILL.md`, lightweight module-level index files, or manually curated high-level notes, but executable outputs should live under the shared four category directories. If a simple task has only one natural step, still use a small secondary-module/task name such as `01-main`, `01-qc`, or `01-epithelial-subclustering`, but create that task directory only under the category directories that receive real outputs.

If one analysis step outputs multiple files or figures, put that output set in the same named secondary-module subdirectory under `figures/`, `tables/`, `codes/`, or `h5ad/`, using the same analysis prefix when possible.

If an output already exists, do not rerun only to recreate it in the new layout. Do not move or delete existing outputs for layout cleanup unless the user explicitly names the exact path and operation. Prefer to leave existing outputs in place, copy them into the organized location only when provenance is recorded, then update the corresponding code paths so future runs write to the same organized location.

When a task creates a run, lineage, candidate parameter set, or method variant, create matching candidate subdirectories inside the active secondary-module/task directory only under category parents that receive outputs for that candidate. For multi-candidate or multi-condition runs, use matching candidate names under each relevant parent when needed; for example, create `h5ad/01-integration-parameter-search/pcs-30_nn-15_res-0p8/` only if an AnnData-like object will be saved, and create `figures/01-integration-parameter-search/pcs-30_nn-15_res-0p8/` only if figures will be saved. Keep h5ad-like candidate objects under the module's shared `h5ad/` with secondary-module and parameter-coded subdirectories, for example `h5ad/01-integration-parameter-search/pcs-30_nn-15_res-0p8/adata_inte.h5ad`; keep candidate code files under the shared `codes/` with matching secondary-module and parameter-coded paths when code is emitted.

Each analysis that produces an output should have corresponding source code under the current module's `codes/`. Acceptable code artifacts include `.ipynb`, `.py`, `.R`, and `.sh`, depending on the language actually used. Do not leave a figure, table, or exported result that can only be traced to manual GUI editing. If an analysis uses Python, keep the notebook and/or `.py` script that generates it; if it uses R, keep the `.R` script or R notebook; if both languages are used, keep both code artifacts under `codes/` with clear ordered prefixes. When converting notebooks to upload/download versions, keep the executable cells needed to reproduce the outputs and remove stale display output only when requested.

Each executed run should also create or update a parameter/provenance report under `tables/`, such as `tables/run_parameters.txt`, `tables/run_parameters.csv`, or a step-specific report in the same output subdirectory. The report should list the code file used, input files/objects, output files, exact parameters, random seeds, selected candidate/final settings, skipped steps, fallback decisions, and any user-approved method changes.

Do not substitute another analysis method, algorithm, statistical test, visualization strategy, database, or input layer without explicit user permission. If the specified method cannot run, stop that module, document the blocker in `tables/readme.txt`, and ask for confirmation before using any alternative. Any approved or documented method change should state why the original method was unsuitable or failed and why the replacement method is appropriate for the same analysis goal.

When a task, notebook run, script run, or long interactive kernel finishes, promptly close the process/kernel/session and release CPU memory and GPU memory. Do not leave idle Python, R, Jupyter, CellChat, RAPIDS, PyTorch, TensorFlow, or CUDA processes holding RAM/VRAM after the requested work is complete.

After each module finishes, create or update `tables/package_versions.txt` describing the packages and tools used by that module. Include Python packages, R packages, command-line tools, CUDA/GPU libraries when relevant, interpreter/R version, environment name or path, and the code files that used them.

Install missing dependencies when they are required to execute the specified method or its approved acceleration path. This includes installing a compatible GPU-accelerated implementation when the method supports it and the machine has a usable GPU/CUDA driver, for example installing `rapids-singlecell`/RAPIDS to run Scanpy-style preprocessing through `rsc`. Dependency installation is allowed to make the requested method work; method substitution is not allowed without explicit user permission. For packages or methods that already provide GPU acceleration, enable and use the GPU-accelerated path after installing any missing compatible GPU packages and verifying imports/minimal execution. If the requested package/method has no GPU-accelerated implementation, use its normal CPU path. If an expected GPU path is installed but broken or incompatible for a non-OOM reason, including a CUDA/CUDA-tag mismatch such as `cu11` vs `cu12` wheels, CuPy/RAPIDS/PyTorch wheels incompatible with the visible driver, missing CUDA runtime libraries, `libucx`/UCX errors, or `cuCtxGetDevice`/CUDA context errors, first try to repair or reinstall a compatible GPU environment without changing the requested method. Choose a compatible wheel, channel, or uv environment automatically from `nvidia-smi`, Python version, platform, and package compatibility information; do not ask the user to choose the CUDA tag. Ask the user only before system-driver changes, OS package changes that require elevated privileges, deleting an existing environment, or replacing a working environment used by other analyses. If GPU runs out of memory, inspect active GPU processes, close stale or idle processes left by previous tasks/kernels when they can be safely identified, release VRAM, record the OOM, and switch that OOM step directly to the equivalent CPU implementation; do not repair, reinstall, or repeatedly retry GPU solely for OOM. Do not terminate unrelated active user processes unless the user explicitly approves. If compatible dependency/GPU installation or non-OOM repair fails, or no usable GPU is present, document the reason and continue with the normal CPU path for the same requested method.

This GPU backend rule applies to all GPU-capable code in every module. If the first full execution of the required task completes successfully with the planned backends, do not rerun only to validate the backend plan. If a GPU-accelerated step fails from GPU OOM, release VRAM, record the OOM, and treat CPU fallback for that step as pre-approved. If a GPU-accelerated step fails for a non-OOM reason and the user approves CPU fallback after documented repair/cleanup attempts, treat fallback at step granularity: use the CPU equivalent only for the failed step when needed, then continue later steps with GPU whenever those later steps have a valid GPU path. Do not mark the entire remaining workflow CPU-only because one GPU step failed. Record a backend capability table under the relevant `tables/<secondary-module>/` or `tables/<secondary-task>/` directory with one row per executed step, including `step`, `planned_backend`, `attempted_backend`, `status`, `error_summary`, `fallback_backend`, `clean_input_reloaded`, and `final_backend_for_rerun`. Also export `gpu_backend_capability_summary.csv` and, when useful, `gpu_backend_capability_summary.txt` in that same tables directory; these files must state which steps can use GPU, which steps must use CPU, and the reason for each CPU step. When any GPU failure, fallback, or partial object mutation occurs during this exploratory/profiling pass, finish the required task only to learn which steps can use GPU, then start a fresh Python/R process, reload the nearest clean upstream input h5ad/RDS/input files from disk, and rerun the whole required task once using the recorded final backend plan. In that final rerun, every step marked GPU-capable must use GPU, except steps with recorded GPU OOM or approved non-OOM CPU fallback, which should use the recorded CPU fallback so the final run avoids repeating known mid-run GPU failures while still using GPU wherever it works. Do not reuse in-memory AnnData/R objects, arrays, GPU buffers, fitted models, graphs, or partial metadata from the profiling pass. Do not present outputs from a partial-failure/profiling pass as canonical final outputs.

Python package management rule: use `uv` by default for Python dependencies, and save uv-managed environments under a dedicated subdirectory in the total analysis directory for the dataset/project. Create the layout `uv_envs/<category>/.venv` under the top-level analysis root, where `<category>` is a stable dependency category such as `main`, `rapids`, `velocity_cellrank`, `cellchat_liana`, or `survival`. Use `uv_envs/main/.venv` for the default shared Python stack, and create another category only when dependency compatibility requires it. Install with `uv pip install --python uv_envs/<category>/.venv/bin/python --index-url https://pypi.tuna.tsinghua.edu.cn/simple ...` or activate that category `.venv` before `uv pip install --index-url https://pypi.tuna.tsinghua.edu.cn/simple ...`. Do not create a separate per-module environment unless it is also a documented dependency category, do not create or reuse a global environment, do not put the run environment inside the skill source directory, and do not run bare `pip install` into the system/user Python unless `uv` is unavailable or the user explicitly requests it. Keep any category-level `pyproject.toml`, `uv.lock`, or requirements export inside `uv_envs/<category>/`, and record all environment paths, categories, package versions, and mirror fallbacks under the total analysis directory, usually `tables/package_versions.txt`. Prefer TUNA mirrors for package downloads: use the TUNA PyPI mirror for `uv pip install`, and use TUNA CRAN/Bioconductor mirrors for R packages when practical. If the TUNA mirror is unavailable, stale, or missing a required package, fall back to the official source only for the affected dependency and record the mirror fallback in `tables/package_versions.txt`. If `uv` is unavailable, install `uv` first when possible; otherwise document the fallback package manager in `tables/package_versions.txt`. R packages, Cell Ranger, velocyto, CUDA drivers, and system libraries are outside `uv` and should be installed with their appropriate manager while still recording versions.

Do not create a symlink in the current directory that points to an older output just to make it look renamed. If an output needs a new name, either regenerate it in the correct output directory or make a real copied file with documented provenance.

Python random seed rule: every Python script and notebook must define and use the default fixed random seed near the top of the file: `SEED = 42`, unless the user explicitly specifies another seed. Set `random.seed(SEED)` and `numpy.random.seed(SEED)` when those libraries are used, and set framework-specific seeds for stochastic packages when relevant, such as PyTorch, TensorFlow, scvi-tools, scVelo, veloVI, CellRank, scikit-learn, Scanpy, or UMAP. Pass `random_state=SEED`, `seed=SEED`, or the package-equivalent argument to every function that supports it, including PCA, neighbors/UMAP, Leiden/clustering, train/test splits, model fitting, bootstrapping, permutations, and plotting layouts when applicable. If a function has no seed argument or still has nondeterministic GPU behavior, document that limitation. Record the fixed seed in the module parameter log, per-candidate `harmony_params.txt` or equivalent, and `tables/package_versions.txt`; do not change the seed between candidate groups unless the user explicitly requests seed sensitivity testing.

Expression-layer rule: for expression-based plotting, marker visualization, dotplots, violin plots, gene scoring, or other Scanpy-compatible calculations, if the function exposes a `use_raw` argument and `adata.raw` is present, set `use_raw=True` by default unless the user explicitly requests a different layer or the method requires counts. Record any exception and the expression layer actually used.

Cell-count reduction rule: do not downsample or subsample cells for analysis, plotting, parameter tuning, or resource-control fallback unless the user explicitly requests or approves that specific reduced-cell run. If the full data cannot run, document the blocker and stop to ask instead of producing a reduced-cell result automatically. Any approved reduced-cell output must be labeled as reduced-cell/exploratory and must not be presented as the full-data result.

Default figure formats are PDF and/or SVG only. Do not create, save, convert, or request PNG files for final, intermediate, diagnostic, preview, thumbnail, or temporary figure outputs. If a tool defaults to PNG, override it to PDF/SVG or stop and ask; do not leave `.png` files under `figures/`, `tables/`, `codes/`, or `h5ad/`.

Global figure style for every module:
- Use the module-specific `Module Figure Style Contract` in this SKILL.md as the plotting-style source of truth. Do not rely on any material outside this SKILL.md to infer final figure style. Keep exactly one current canonical plotting route for each final figure, and label any other route as non-canonical or exploratory; do not duplicate final PDFs/SVGs for the same panel or file stem.
- Plot routing is code-dispatch based. Identify the figure family first, then call the corresponding canonical plotting code or package plotting route specified in this skill, analogous to using the correct `sc.pl.*` function for each Scanpy plot. Do not use a generic Matplotlib, Seaborn, ggplot2, or ad hoc plotting pattern for a figure family that has a dedicated canonical route. If no matching route exists in this skill for a required final figure, stop and ask for the skill to be updated before running; do not improvise a new plotting route.
- Apply the style for the corresponding analysis type and module-specific section first; use these global rules only as baseline guardrails. Do not copy a figure style from an unrelated analysis type just because the file format or package is similar.
- Save figure outputs as PDF and/or SVG only; never create PNG previews, thumbnails, diagnostics, or temporary figure files.
- Keep text editable whenever the backend supports it. In Matplotlib set `pdf.fonttype = 42`, `ps.fonttype = 42`, `svg.fonttype = "none"`, and use a standard sans-serif font such as DejaVu Sans or Arial.
- Keep axis ticks visible on quantitative and categorical plots, including heatmaps, dotplots, barplots, forest plots, scatter plots, survival plots, and UMAP panels with axes. Do not call `axis("off")`, remove tick labels, or hide spines unless the plot type is a pure network/chord/graph layout where axes have no coordinate meaning.
- Use clean white backgrounds, black axis text, readable tick labels, and legends/colorbars outside or to the right when practical. Rotate dense x-axis labels rather than letting them overlap.
- Prefer the official plotting interface for the relevant package before manual low-level drawing. Use manual Matplotlib/ggplot2/Seaborn layout code only when the package interface cannot express the required final figure or when this skill gives explicit canonical manual code; record that reason in the code or parameter log.
- Size every figure element for the exported canvas: labels, legends, colorbars, risk tables, titles, p/q labels, arrows, node labels, and panel titles must be readable and must not collide. If any element overlaps or is clipped in the generated PDF/SVG, increase canvas size, margins, row/column spacing, legend placement, font size, label wrapping, or panel spacing and rerender before considering the figure complete.
- For UMAP-like embedding panels, keep panels square, use the official Scanpy/scVelo plotting interface when possible, keep framed axes, and use consistent palettes for the same labels across full-atlas and lineage-specific plots.
- For heatmaps, use a colorbar with visible ticks. Use centered diverging palettes for signed correlations/effects, sequential palettes for nonnegative scores or transition weights, and method-specific labels in titles/captions.
- For dotplots, keep gene and group orders explicit, preserve tick labels, and use standard scaling only when the method requires display scaling.
- For forest/survival plots, keep confidence intervals, hazard-ratio/reference lines, p/q labels, and axis ticks visible.
- For network/chord/directed-graph plots, use PDF/SVG, editable labels, clear legends/colorbars when edge weights are encoded, and document any intentional axis removal.


Scanpy plotting hard rule: For ordinary `sc.pl.*` outputs that support `save=...`, do not create Matplotlib axes, do not pass `ax=...`, and do not save with `fig.savefig` or `ax.figure.savefig`. Use `sc.settings.figdir` plus the Scanpy `save` argument. For multi-panel Scanpy plots, use the package interface such as `color=[...]`, `ncols`, `wspace`, or `standard_scale` instead of manual subplots. Manual `ax=...` is allowed only for documented special overlays or incompatible per-panel settings that Scanpy cannot express; record the reason in the run-parameter table or a code comment.


For both Python and R plotting, use the official plotting interface of the relevant package by default, such as Scanpy, Matplotlib, Seaborn, CellRank, scVelo, ggplot2, ComplexHeatmap, or CellChat plotting APIs. For Python-generated single-panel plots, use a default square canvas of 2.5 x 2.5 inches unless the user specifies another size or the plot type clearly requires more space, such as multi-panel layouts, heatmaps, wide dotplots, survival plots, or network/chord diagrams. For UMAP plots, use the default Scanpy save path and keep the call minimal. Set figure parameters with `sc.set_figure_params(figsize=(3, 3), dpi=150)` or `sc.settings.set_figure_params(figsize=(3, 3), dpi=150)`, set `sc.settings.figdir` to the target output directory, then call `sc.pl.umap(adata, color="leiden_coarse", save="_name.pdf")`. Keep the default Scanpy-style framed axes and outside legends; do not manually create `plt.subplots`, pass `ax=...`, or call `fig.savefig` for ordinary UMAPs. Use `ncols` only when `color` contains multiple objects, for example `sc.pl.umap(adata, color=["cell_subtype", "status", "cnv_score"], ncols=3, wspace=0.4, save="_celltype_status_cnvscore.pdf")`. Do not use `ncols` for a single-color UMAP. Keep `save` as a suffix/name handled by Scanpy rather than a full path. Manual axes and `fig.savefig` are reserved for special cases where different panels require incompatible per-panel parameters or post-processing that the official `color=[...]` interface cannot express, and the reason must be documented because manual saving can greatly increase PDF/SVG size. For ordinary UMAPs, specifically avoid `return_fig=True` followed by `ax.savefig(..., bbox_inches="tight")`; on large atlases this can inflate files from sub-MB Scanpy-saved outputs to multi-MB PDFs or very large SVGs. When the task is only to inspect, audit, or explain existing UMAP code/output size, do not modify the source code or rerun plotting unless the user explicitly asks for a fix or rerender. Keep all figure text editable as text whenever the plotting backend supports it; do not convert labels, legends, tick labels, titles, or annotations to outlines/paths unless the user explicitly requests it. For multi-panel figures, especially multi-panel UMAPs, verify that each UMAP panel remains square, legends/colorbars/titles do not overlap, and adjacent panels do not collide. It is acceptable to adjust the default figure width, height, `wspace`, `hspace`, legend font size, or margins to prevent overlap and preserve square UMAP panels. Do not add extra custom titles to UMAP panels; use the Scanpy default title derived from `color` unless the user explicitly asks for custom titles. Do not add automatic bitmap/raster conversion rules; let the user decide figure-size tradeoffs from the actual output files. Do not draw sample-colored UMAPs as default final figures; generate sample-colored UMAPs only when the user asks for them or when they are needed as integration/batch-mixing diagnostics, and label them as diagnostic outputs.

## Module Figure Style Contract

Use the following spatial-validation figure styles unless the user explicitly
asks for a different style. Do not mention prior-run implementation provenance in generated reusable code, figure labels, captions, or readme files.

- Single-feature spatial maps: use real tissue coordinates, default framed
  coordinate panels, no manual arrows or artificial axis labels, and consistent
  color limits within each feature family. For publication figures, prefer
  clear contrast and moderately saturated color scaling over washed-out maps,
  while recording the chosen `vmin`/`vmax` or percentile clipping.
- Per-sample spatial grids: keep every spot panel square or coordinate-faithful,
  keep sample titles compact, and preserve the same feature/color mapping across
  all samples in a grid.
- Four-quadrant CM-lineage spatial maps: use a fixed four-color legend for
  low-low, CM-high only, lineage-high only, and both-high spots; keep thresholds
  explicit in the table and figure metadata. Place the four-color quadrant
  legend at the upper-right of the corresponding quadrant panel, outside the
  plotted map body when needed, so it never covers spots. The legend must not
  overlap the map body, colorbars, panel titles, tick labels, or neighboring
  panels; if any overlap occurs, enlarge the canvas or right margin and
  rerender. Every spatial panel in the raw two-panel figure and percentile
  three-panel figure must be square; legends and colorbars must not distort
  panel aspect ratios. The figure should show the spatial map only, not extra
  explanatory text blocks.
- Fisher/Stouffer heatmaps: show signed effects with a centered diverging
  palette, visible colorbar ticks, q-value stars, and `ns` for non-significant
  entries when labels are displayed. Do not annotate heatmap cells with the
  underlying numeric values such as rho, odds ratio, Z score, p value, or q
  value. If annotation is enabled, the only text inside cells should be
  `ns`, `*`, `**`, or `***`. Keep 12-sample or per-sample panels separate when
  the task requests only one of them.
- Spatial forest/summary panels: keep confidence intervals or signed-Z/error
  intervals visible, include a zero/reference line, sort rows by the documented
  statistic, and keep axis ticks and labels visible. For the default
  per-sample CM-Epi forest, plot every CM-Epi pair present in the per-sample
  correlation table. Use a selected/core subset only when the user explicitly
  requests that subset.

Other spatial mapping methods can be used if the user requests them, but they are user-defined alternatives. Do not present them as the default.

## Method Integrity Rules

The canonical Module 05 workflow is Tangram-based. Runtime pressure, missing optional packages, or slow marker computation are not valid reasons to change the scientific method. A change recorded in `run_parameters.txt` is still a method change and is not allowed unless the user explicitly approved that change before execution.

Hard rules:

1. If Tangram is missing, install a compatible `tangram-sc`/PyTorch environment and run Tangram. Do not replace Tangram with `scipy.optimize.nnls`, non-negative least squares, linear regression, cosine matching, nearest-neighbor matching, `scanpy.tl.score_genes`, marker-score projection, average marker expression, AUCell/UCell/AddModuleScore-style scoring, or any other "Tangram-equivalent" method in the canonical branch.
2. If Tangram installation or execution cannot be fixed in the current environment, stop the canonical branch and report the blocker. An NNLS, marker-based `score_genes`, or other non-Tangram fallback may only be created as a clearly named alternative branch after explicit user approval. Such a fallback must never be labeled "Block 04 spatial validation completed", "canonical", "Tangram", or "method16"; label it `noncanonical_marker_score_exploratory` or another explicit user-approved alternative name.
3. Marker genes must come from Module 02's saved post-annotation subtype DEG CSVs grouped by the final biological `cell_subtype` labels. Use the annotation-completed DEG outputs, not raw-cluster DEG tables, preset marker lists, top-expression tables, or an in-memory AnnData object.
4. Do not recompute marker DEG inside the spatial-validation block, even when an input AnnData already contains `obs['cell_subtype']`. The presence of `cell_subtype` is not permission to call `sc.tl.rank_genes_groups` here. If the required final subtype DEG CSVs are absent, stop and request them. If the user explicitly requests DEG rebuilding, rerun Module 02's post-annotation DEG step first, save its full per-subtype CSV outputs, and only then restart spatial validation from those saved tables.
5. Do not replace missing or slow post-annotation DEG outputs with mean-expression ranking, expression-rate ranking, hard-coded markers, or DEG recalculated ad hoc in the spatial script. Use the existing annotation-completed DEG tables or stop and ask.
6. Spatial pair testing is per epithelial subtype and CM pair: `(Epi_subtype_fraction, CM_activity)`. By default, test every valid CM-Epi pair present in the spot-score table. Do not replace the epithelial subtype feature with `total_Epi_fraction` in the canonical pair tests.
7. `total_Epi_fraction` can be reported only as a separately named broad-epithelial sensitivity branch if the user asks for it. It must not overwrite or substitute the per-subtype CM-Epi Fisher/Stouffer or Spearman outputs.

Required Tangram preflight gate before any spatial mapping:

```python
# FIXED: canonical spatial validation requires real Tangram, not marker scores.
# Run this before computing spot scores. If it fails, install/repair Tangram or stop.
import importlib.util

if importlib.util.find_spec("tangram") is None:
    raise RuntimeError(
        "Tangram is required for canonical spatial validation. "
        "Install tangram-sc/PyTorch and rerun; do not substitute score_genes, "
        "marker-based scores, NNLS, regression, or nearest-neighbor matching."
    )

import tangram as tg
TANGRAM_AVAILABLE = True
```

Canonical Tangram call and parameters are fixed:

```python
# FIXED: this is the canonical Tangram mapping call for this workflow.
# `asc` is the subtype pseudobulk/reference AnnData, not raw single-cell rows.
# `asp` is the per-sample spatial AnnData normalized/log-transformed upstream.
tg.pp_adatas(asc, asp, genes=common)

mapper = tg.map_cells_to_space(
    asc,
    asp,
    mode="cells",
    device="cuda:0",
    num_epochs=350,
    learning_rate=0.05,
    random_state=42,
    verbose=False,
)

tg.project_cell_annotations(mapper, asp, annotation="cell_subtype")
abundance = asp.obsm["tangram_ct_pred"].copy()

# FIXED: EPIfrac is an epithelial-composition fraction, not the unmodified
# Tangram subtype abundance. Normalize only across epithelial subtype columns
# within each spatial observation so the retained EPIfrac columns sum to 1.
epi_subtypes = [
    subtype for subtype in abundance.columns
    if str(subtype).startswith("Epi_")
]
if not epi_subtypes:
    raise ValueError("No epithelial subtype columns found in tangram_ct_pred")
epi_abundance = abundance.loc[:, epi_subtypes].astype(float)
epi_total = epi_abundance.sum(axis=1)
epi_frac = epi_abundance.div(
    epi_total.replace(0, np.nan),
    axis=0,
).fillna(0.0)

# Audit invariant: rows with positive epithelial abundance must sum to 1;
# rows with zero total epithelial abundance remain all-zero after fillna(0).
positive_epi = epi_total > 0
if positive_epi.any() and not np.allclose(
    epi_frac.loc[positive_epi].sum(axis=1).to_numpy(),
    1.0,
    rtol=1e-6,
    atol=1e-8,
):
    raise ValueError("EPIfrac row-normalization failed")
```

`EPIfrac__<epithelial_subtype>` must always be written from `epi_frac` above.
Do not write the unmodified `tangram_ct_pred` epithelial columns under an
`EPIfrac__*` name. This is an epithelial-internal composition: its denominator
is the sum of all projected epithelial subtype abundances at the same spatial
observation, not all cell types and not a dataset-wide or sample-wide total.
Do not substitute a single `total_Epi_fraction` for these per-subtype columns.

Do not change these Tangram parameters to reduce runtime. If CUDA/PyTorch cannot
use `cuda:0`, repair the Tangram/PyTorch CUDA environment first. CPU execution
of the same Tangram call may be used only when no usable GPU is available or
after documented GPU repair/OOM handling under the global backend rule; changing
Tangram to marker scoring or another method is still forbidden.

The canonical per-sample spot-score CSV must be derived downstream of Tangram
mapping. If a code path calls `sc.tl.score_genes`, computes marker averages, or
uses a marker/module-score table to create `CMact__*` or `EPIfrac__*` columns,
that code path is noncanonical and must stop unless the user explicitly approved
an exploratory alternative branch.

## Inputs

Required:

```text
spatial h5ad per sample with real tissue coordinates
single-cell epithelial subtype pseudobulk/reference profiles
CM node/activity definitions
gene expression in spatial h5ad after per-sample normalize/log
```

Do not use matrix row/column indices as tissue coordinates if real spatial coordinates are available.

## Marker Source and Gene Subsetting

Default marker inputs:

```text
Module 02 post-annotation epithelial subtype DEG CSVs grouped by cell_subtype
Module 02 post-annotation non-epithelial subtype DEG CSVs grouped by cell_subtype
```

Canonical DEG source contract:

```text
tables/05-subtype-deg-annotation/<lineage>/degs_cell_subtype_pcs<PCS>_nn<NN>_res<RES>/
{cell_subtype}_degs_cell_subtype_pcs<PCS>_nn<NN>_res<RES>.csv
```

Use the equivalent annotation-completed output path only when the project uses
a different documented layout. Require one full DEG table per final
`cell_subtype`, with the final subtype name represented by the table group or
an explicit group column. Do not use pre-annotation raw-cluster DEG tables.
Do not calculate a replacement DEG table directly from a supplied
`cell_subtype` AnnData inside this block. Save the exact selected source path
for every subtype in `marker_genes_used.csv` as `source_deg_csv`.

Canonical reference-construction parameters:

```text
marker count per subtype = top 100 significant positive DEG genes
significant positive DEG rule = pvals_adj < 0.05 and logfoldchanges > 0
marker ordering = score descending only
reference cells = every available cell in each subtype
cell cap = none
random cell subsampling = none
reference expression = subtype mean from annotated single-cell adata.raw
```

Require numeric `pvals_adj`, `logfoldchanges`, and `score` columns in every
selected post-annotation subtype DEG table. Filter `pvals_adj < 0.05` and
`logfoldchanges > 0`, sort by `score` descending, keep source-table order for
exact score ties, drop duplicated gene names while retaining the highest-score
occurrence, and then take the first 100 genes. Use `pvals_adj` only as the
required significance filter; do not use `pvals`, `pvals_adj`, or q values as
sorting keys. If fewer than 100 significant positive unique genes remain, use
all of them and record the count. If `pvals_adj`, `score`, or `logfoldchanges`
is missing, stop and request a valid post-annotation DEG table instead of
inventing a fallback ordering.

These parameters are fixed defaults for the canonical branch. Do not restore a
top-12 marker limit, a 700-cell-per-subtype cap, or any other automatic
subsampling. Save `marker_genes_used.csv`, `reference_cells_used.csv`, and the
union/intersection gene list so the top-100 selection and all-cell usage are
auditable. `reference_cells_used.csv` must contain at least `subtype`,
`n_cells_total`, `n_cells_used`, and `all_cells_used`; require
`n_cells_used == n_cells_total` and `all_cells_used == True` for every subtype.

These marker tables must be the saved post-annotation subtype differential-expression results from Module 02. Use the union of selected marker genes and intersect it with genes present in the single-cell reference and spatial h5ad. Do not restrict the spatial sample to highly variable genes for this Tangram validation workflow. Do not generate marker tables from mean-expression ranking or by rerunning DEG on the spatial-block AnnData.

In other words:

```text
use DEG marker/intersected genes
do not use HVG-only spatial matrices
do not use top-expression-only marker substitutes
```

## Expression Contract

Single-cell reference side:

```text
use annotated single-cell adata.raw to build epithelial subtype pseudobulk/reference profiles
group cells by epithelial subtype
average expression per subtype over marker/intersected genes
use the AnnData only for subtype membership, cell counts, and expression means; never derive the marker DEG ranking from it in this block
do not use scaled X, PCA, integrated embeddings, or batch-corrected latent spaces for pseudobulk expression
if adata.raw is absent, use a documented normalized/log expression layer and record the fallback
```

For raw spatial h5ad:

```text
read each sample h5ad
normalize total per sample
log1p per sample
use adata[:, present_genes].X for marker expression
```

The `adata.raw` rule refers to the single-cell reference object. The spatial sample object should use normalized/log `.X` produced per sample.

Required Tangram spatial h5ad preprocessing pattern:

```python
def normalize_spatial_sample(adata: sc.AnnData, target_sum: float = 1e4) -> sc.AnnData:
    """Normalize/log a clean per-sample spatial h5ad before Tangram scoring."""
    adata.var_names_make_unique()
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    else:
        adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    return adata


def standardize_spatial_obs_names(adata: sc.AnnData, sample: str) -> sc.AnnData:
    if "barcode" not in adata.obs.columns:
        adata.obs["barcode"] = adata.obs_names.astype(str)
    adata.obs["sample"] = sample
    adata.obs_names = [f"{barcode}_{sample}" for barcode in adata.obs["barcode"].astype(str)]
    if not adata.obs_names.is_unique:
        raise ValueError(f"Non-unique obs_names after sample suffixing: {sample}")
    return adata


def read_spatial_sample(spatial_h5ad_dir: Path, sample: str, filename: str, genes: list[str]):
    adata = sc.read_h5ad(spatial_h5ad_dir / filename)
    adata = standardize_spatial_obs_names(adata, sample)
    adata = normalize_spatial_sample(adata)
    if "spatial" not in adata.obsm and {"array_col", "array_row"}.issubset(adata.obs.columns):
        adata.obsm["spatial"] = adata.obs[["array_col", "array_row"]].to_numpy(dtype=float)
    present_genes = [g for g in genes if g in adata.var_names]
    X_marker_log = adata[:, present_genes].X
    obs = adata.obs[["array_row", "array_col", "sample"]].copy()
    obs.index = adata.obs_names.astype(str)
    return adata, present_genes, X_marker_log, obs
```

This pattern is for clean/raw per-sample spatial h5ad input. Do not start from a
previously CM-assigned spatial h5ad for the canonical Tangram expression matrix
unless the user explicitly requests a reuse branch. If a provided h5ad has
already been normalized/logged, detect and record that fact; do not normalize it
twice.

## Main Output

For every sample:

```text
{sample}_tangram_pseudobulk_epi_cm_spot_scores.csv
```

This per-sample spot-score CSV is mandatory. It must contain the spot-level
epithelial fractions and CM activities needed for downstream statistics and
maps. Save this CSV immediately after Tangram mapping for each sample, before
any plotting, Spearman, Fisher, or Stouffer calculation. All downstream
canonical plots and statistics must reload this CSV and must not use an
in-memory Tangram object or a Tangram-augmented h5ad as the plotting/statistics
input. If an existing project file uses a legacy internal filename, treat it as
the same spot-score table but do not expose the internal method number in
captions or reusable skill instructions.

Canonical handoff:

```text
Tangram pseudobulk mapping -> per-sample spot-score CSV -> reload CSV -> plots/statistics
```

An h5ad with Tangram results written into `.obs`, `.obsm`, `.uns`, or `.layers`
may be saved as an optional intermediate/audit object when useful, but it is not
required for canonical plotting or statistics. Canonical plotting and
statistics must still reload the per-sample CSV. A spatial h5ad may be read as
the clean expression input before Tangram, as an optional saved audit object, or
as a coordinate fallback if the CSV is missing coordinates; the canonical CSV
should already contain the spot IDs and spatial coordinates needed for plotting.

The canonical spot-score values are the original numeric Tangram-derived
fractions/scores and CM activities. Do not replace these columns with min-max
scaled, robust-scaled, percentile, z-scored, clipped, rank-transformed, or
plot-normalized values. Any transformed value used for visualization must be
saved in a separate clearly named column, such as `<feature>_percentile` or
`<feature>_plot_scaled`, while preserving the original raw numeric column.

Minimum useful columns:

```text
sample
spot_id
barcode
array_row
array_col
spatial_x
spatial_y
EPIfrac__<epithelial_subtype> columns
CMact__<CM> columns
```

Canonical per-sample spot-score writer:

```python
# REPLACEABLE: obs contains spot metadata from the spatial h5ad, indexed by spot
# ID; cm_df has CM activity columns indexed by spot ID; epi_frac has epithelial
# subtype fraction columns indexed by spot ID.
# FIXED: write the CSV before any plotting/statistics; later code must reload it.
def save_spot_scores(sample, tables_dir, obs, cm_df, epi_frac):
    score = obs.copy()
    score.insert(0, "sample", sample)
    score.insert(1, "spot_id", score.index.astype(str))

    if "barcode" not in score.columns:
        score["barcode"] = score["spot_id"].astype(str)
    if "array_row" not in score.columns and "spatial_y" in score.columns:
        score["array_row"] = score["spatial_y"]
    if "array_col" not in score.columns and "spatial_x" in score.columns:
        score["array_col"] = score["spatial_x"]
    if "spatial_x" not in score.columns and "array_col" in score.columns:
        score["spatial_x"] = score["array_col"]
    if "spatial_y" not in score.columns and "array_row" in score.columns:
        score["spatial_y"] = score["array_row"]

    score = score.join(cm_df.add_prefix("CMact__"), how="left")
    score = score.join(epi_frac.add_prefix("EPIfrac__"), how="left")

    required = ["sample", "spot_id", "array_row", "array_col"]
    missing = [col for col in required if col not in score.columns]
    if missing:
        raise KeyError(f"{sample}: spot-score CSV missing required columns: {missing}")

    out = tables_dir / f"{sample}_tangram_pseudobulk_epi_cm_spot_scores.csv"
    score.to_csv(out, index=False)
    return out
```

## Canonical Sample Scopes

Run two canonical spatial-validation scopes unless the user explicitly requests
only one:

```text
all-samples
tumor-only
```

`all-samples` uses every spatial sample that passed input validation.
`tumor-only` uses only samples labeled as tumor in the user-provided sample
metadata or an explicitly documented sample-status mapping. Do not infer tumor
status from filename patterns, sample order, or project memory.

Tangram mapping and spot-score generation must run for all spatial samples that
passed input validation. Do not run Tangram only on tumor samples for the
canonical workflow, and do not create separate tumor-only Tangram mappings. The
`all-samples` and `tumor-only` scopes apply only after spot-score generation,
when building final statistics, cross-sample summaries, and Fisher/Stouffer
heatmaps. If no status mapping is available, still run Tangram for all valid
samples and complete all-samples statistics/heatmaps, then stop to ask for the
tumor/normal status table before creating tumor-only statistics/heatmaps.

Create and save a sample-scope table before statistics:

```text
sample
status
include_all_samples
include_tumor_only
exclusion_reason_all_samples
exclusion_reason_tumor_only
status_source
```

Spot-score generation must be performed once per valid sample. Only final
statistical summary tables and final heatmap figures must be materialized
separately by scope.
Use distinct output paths such as:

```text
tables/<secondary-task>/statistics/all-samples/
tables/<secondary-task>/statistics/tumor-only/
figures/<secondary-task>/heatmaps/all-samples/
figures/<secondary-task>/heatmaps/tumor-only/
```

Per-sample Tangram outputs and per-sample spatial maps should not be duplicated
solely for the tumor-only scope. Do not let tumor-only statistics/heatmaps
overwrite all-samples statistics/heatmaps, and do not present an all-samples
Fisher/Stouffer heatmap as the tumor-only result.

## Default Spatial Validation Logic

For every valid sample and every CM-Epi subtype pair, the canonical analysis is:

```text
1. Estimate spot-level EPIfrac for every valid epithelial subtype.
2. Estimate spot-level CMact for each CM program.
3. Save the original numeric CMact and EPIfrac columns to the per-sample spot-score CSV.
4. Reload the per-sample spot-score CSV and use those original numeric columns for Spearman correlation.
5. From the reloaded CSV, draw one original two-panel spatial figure per sample/pair: CMact and EPIfrac.
6. From the reloaded CSV, convert CMact and EPIfrac to within-sample percentile ranks.
7. From the reloaded CSV, draw one three-panel spatial figure per sample/pair: CMact percentile, EPIfrac percentile, and four-color percentile quadrant.
8. Run Fisher tests from the percentile quadrants per sample/pair.
9. Integrate Fisher directions across samples with signed Stouffer Z.
10. Draw exactly two default final Stouffer heatmaps: all-samples and tumor-only.
```

Higher percentile means higher relative abundance/activity/expression within
that sample. Do not use percentile, robust-scaled, clipped, z-scored, or
plot-scaled values for the raw Spearman correlation; use the original numeric
`CMact__<CM>` and `EPIfrac__<subtype>` columns.

Parallel execution is allowed only at the sample level. It is valid to run
independent workers where each worker owns one sample, reads that sample's
spatial h5ad/spot-score CSV, and writes only to that sample's output directory
and that sample's temporary/statistics table. After all sample workers finish,
merge per-sample statistics once and draw combined heatmaps once. Do not start
multiple full-workflow processes that all iterate over every sample, and do not
let two processes write the same sample, pair, PDF/SVG, or CSV path. If an
existing process is already rendering a sample, a new process must skip that
sample or stop before writing.

## Default Outputs

Default spatial figures/tables are:

```text
figures/<task>/raw_cmact_epifrac_by_sample/<sample>/<sample>__<CM>__<epi_subtype>__raw_cmact_epifrac.pdf
figures/<task>/raw_cmact_epifrac_by_sample/<sample>/<sample>__<CM>__<epi_subtype>__raw_cmact_epifrac.svg
figures/<task>/percentile_quadrant_by_sample/<sample>/<sample>__<CM>__<epi_subtype>__percentile_quadrant_0p5.pdf
figures/<task>/percentile_quadrant_by_sample/<sample>/<sample>__<CM>__<epi_subtype>__percentile_quadrant_0p5.svg
figures/<task>/stouffer_heatmaps/all_samples_stouffer_signedZ_qstars.pdf
figures/<task>/stouffer_heatmaps/all_samples_stouffer_signedZ_qstars.svg
figures/<task>/stouffer_heatmaps/tumor_only_stouffer_signedZ_qstars.pdf
figures/<task>/stouffer_heatmaps/tumor_only_stouffer_signedZ_qstars.svg
tables/<task>/per_sample_spearman.csv
tables/<task>/percentile_quadrant_fisher_per_sample.csv
tables/<task>/percentile_quadrant_fisher_sample_stouffer_all_samples.csv
tables/<task>/percentile_quadrant_fisher_sample_stouffer_tumor_only.csv
```

For per-sample spatial maps, the directory hierarchy is fixed as
`by_sample/<sample>/<sample>__<CM>__<epi_subtype>__...`. The sample directory
must be the first browsing level and the different CM-Epi pair figures must live
inside that sample directory. Do not write canonical spatial map outputs as
`per_pair_spot_maps/<pair>/<sample>...`, `by_pair/<pair>/<sample>...`, or any
other pair-first layout unless the user explicitly requests an additional
non-canonical browsing copy.

Do not generate extra default heatmap families, pooled-spot log2 odds-ratio
heatmaps, per-sample heatmap panels, bivariate overlay panels, split-spot
panels, or single-feature one-panel maps unless the user explicitly asks for
those optional branches. If a per-sample forest is requested, plot all CM-Epi
pairs by default and subset only when the user explicitly asks for selected
pairs.

Default pair scope is mandatory all CM-Epi pairs. For the canonical Module 05
branch, every valid sample must be crossed with every discovered CM and every
discovered epithelial subtype. Do not run only one representative CM, one
representative epithelial subtype, one sample, one user-looking pair, top pairs,
significant-only pairs, same-direction pairs, or selected/core pairs unless the
user explicitly asks for that exact subset as an additional non-canonical branch.
Every generated script must derive `cm_names`, `epi_names`, `samples`, and
`pairs` from the active input tables, save them as a run manifest, and iterate
through all rows of that manifest. If full all-pair plotting is too slow or
large, stop and ask; do not silently reduce the pair set.

Canonical all-pair loop skeleton:

```python
# REPLACEABLE: tables_dir and figures_dir.
# FIXED: default reads per-sample spot-score CSVs and iterates every
# sample x CM x epithelial subtype pair. Do not read Tangram-augmented h5ad
# objects for plotting/statistics.
import pandas as pd

def discover_cm_epi_pairs(score):
    cm_names = [c.split("__", 1)[1] for c in score.columns if c.startswith("CMact__")]
    epi_names = [c.split("__", 1)[1] for c in score.columns if c.startswith("EPIfrac__")]
    cm_names = sorted(set(cm_names))
    epi_names = sorted(set(epi_names))
    pairs = [{"CM": cm, "epi_subtype": epi} for cm in cm_names for epi in epi_names]
    if not pairs:
        raise ValueError("No CMact__/EPIfrac__ pairs found in the spot-score table")
    return cm_names, epi_names, pairs

score_table_paths = sorted(tables_dir.glob("*_tangram_pseudobulk_epi_cm_spot_scores.csv"))
if not score_table_paths:
    raise FileNotFoundError(f"No per-sample Tangram spot-score CSVs found in {tables_dir}")

manifest_rows = []
for score_path in score_table_paths:
    score = pd.read_csv(score_path)
    if "sample" not in score.columns:
        score["sample"] = infer_sample_from_path(score_path)
    if "array_col" not in score.columns and "spatial_x" in score.columns:
        score["array_col"] = score["spatial_x"]
    if "array_row" not in score.columns and "spatial_y" in score.columns:
        score["array_row"] = score["spatial_y"]
    sample_values = sorted(score["sample"].dropna().astype(str).unique().tolist())
    sample = sample_values[0] if len(sample_values) == 1 else infer_sample_from_path(score_path)
    cm_names, epi_names, pairs = discover_cm_epi_pairs(score)
    # FIXED: output directories are sample-first; pairs are files inside each sample.
    raw_sample_dir = figures_dir / "raw_cmact_epifrac_by_sample" / sample
    pct_sample_dir = figures_dir / "percentile_quadrant_by_sample" / sample
    raw_sample_dir.mkdir(parents=True, exist_ok=True)
    pct_sample_dir.mkdir(parents=True, exist_ok=True)

    for pair in pairs:
        cm = pair["CM"]
        epi = pair["epi_subtype"]
        manifest_rows.append({"sample": sample, "CM": cm, "epi_subtype": epi})

        # FIXED default outputs for every sample/pair:
        # 1. raw two-panel CMact/EPIfrac spatial figure
        # 2. percentile three-panel CM/Epi/quadrant spatial figure
        # 3. per-sample Spearman/Fisher rows for downstream Stouffer tables
        plot_raw_cmact_epifrac(score, sample, cm, epi, raw_sample_dir)
        plot_percentile_quadrant(score, sample, cm, epi, pct_sample_dir)
        compute_pair_statistics(score, sample, cm, epi)

pd.DataFrame(manifest_rows).to_csv(tables_dir / "all_sample_cm_epi_pair_manifest.csv", index=False)
```

The two final Stouffer heatmaps must then be built from the complete statistics
over this full manifest: one all-samples table and one tumor-only table. Missing
or invalid pair/sample rows may be excluded only with a saved reason column; do
not drop rows silently.

## Method16 Plotting Style Defaults

All Module 05 plots should follow the final non-old method16 plotting style
unless the user explicitly asks for a different style. The raw CMact/EPIfrac,
percentile, and four-quadrant spatial maps are CSV-driven: read
`{sample}_tangram_pseudobulk_epi_cm_spot_scores.csv`, derive all plot vectors
from that table, and do not use a Tangram-augmented h5ad for the plotted values.

```python
# FIXED baseline style for Module 05 figures.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "font.family": "DejaVu Sans",
})

# FIXED heatmap significance labels. Do not put numeric values in heatmap cells.
def q_label(q):
    if not np.isfinite(q):
        return "ns"
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return "ns"

# FIXED spatial four-quadrant colors.
CATEGORY_COLORS = {
    "both_low": "#d9d9d9",
    "epi_high_only": "#f59e0b",
    "cm_high_only": "#2b6cb0",
    "both_high": "#b91c1c",
}

# FIXED percentile maps.
CM_PERCENTILE_COLORS = ["#f2f2f2", "#b7d7ea", "#4c78a8", "#2b6cb0"]
EPI_PERCENTILE_COLORS = ["#f2f2f2", "#fde6b3", "#fbbf24", "#f59e0b"]

def robust_limits(values, q=(0.01, 0.99)):
    arr = pd.to_numeric(values, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna().to_numpy()
    if arr.size == 0:
        return 0.0, 1.0
    lo, hi = np.nanquantile(arr, q)
    if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
        lo, hi = float(np.nanmin(arr)), float(np.nanmax(arr))
    if lo == hi:
        hi = lo + 1e-9
    return float(lo), float(hi)

def percentile_rank(values):
    s = pd.to_numeric(values, errors="coerce")
    return s.rank(method="average", pct=True).fillna(0.0)

def classify_percentiles(cm_pct, epi_pct, cutoff=0.5):
    cm_high = pd.to_numeric(cm_pct, errors="coerce") >= cutoff
    epi_high = pd.to_numeric(epi_pct, errors="coerce") >= cutoff
    category = pd.Series("both_low", index=cm_pct.index, dtype="object")
    category.loc[epi_high & ~cm_high] = "epi_high_only"
    category.loc[cm_high & ~epi_high] = "cm_high_only"
    category.loc[cm_high & epi_high] = "both_high"
    return category

def set_spatial_axes(ax, x, y):
    ax.invert_yaxis()
    ax.set_box_aspect(1)
    ax.set_xlabel("array_col")
    ax.set_ylabel("array_row")

```

Use the plotting functions below verbatim for spatial pair figures. Do not
replace them with shorter scripts that use `figsize=(6, 3)`, `figsize=(9, 3)`,
`ax.set_aspect("equal")` without `ax.set_box_aspect(1)`, missing legends, or
PDF-only output. Those abbreviated versions violate the required figure
contract.

Canonical pair plotting functions:

```python
from pathlib import Path
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D

def save_pdf_svg(fig, stem: Path):
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", dpi=300)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", dpi=300)
    plt.close(fig)

def plot_raw_cmact_epifrac(score, sample, cm, epi, out_dir, size):
    cm_raw = pd.to_numeric(score[f"CMact__{cm}"], errors="coerce")
    epi_raw = pd.to_numeric(score[f"EPIfrac__{epi}"], errors="coerce")
    x = score["array_col"].to_numpy(dtype=float)
    y = score["array_row"].to_numpy(dtype=float)

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6), constrained_layout=True)
    sc_cm = axes[0].scatter(
        x, y, c=cm_raw, s=size, cmap="magma",
        norm=Normalize(*robust_limits(cm_raw)),
        marker="o", linewidths=0, rasterized=True,
    )
    set_spatial_axes(axes[0], x, y)
    axes[0].set_title(f"CMact: {cm}")
    fig.colorbar(sc_cm, ax=axes[0], fraction=0.046, pad=0.02).set_label("CMact")

    sc_epi = axes[1].scatter(
        x, y, c=epi_raw, s=size, cmap="viridis",
        norm=Normalize(*robust_limits(epi_raw)),
        marker="o", linewidths=0, rasterized=True,
    )
    set_spatial_axes(axes[1], x, y)
    axes[1].set_title(f"EPIfrac: {epi}")
    fig.colorbar(sc_epi, ax=axes[1], fraction=0.046, pad=0.02).set_label("EPIfrac")

    fig.suptitle(f"{sample}: {cm} | {epi}", y=1.04, fontsize=9, fontweight="bold")
    out_dir.mkdir(parents=True, exist_ok=True)
    save_pdf_svg(fig, out_dir / f"{sample}__{cm}__{epi}__raw_cmact_epifrac")

def plot_percentile_quadrant(score, sample, cm, epi, out_dir, size):
    cm_pct = percentile_rank(score[f"CMact__{cm}"])
    epi_pct = percentile_rank(score[f"EPIfrac__{epi}"])
    category = classify_percentiles(cm_pct, epi_pct, cutoff=0.5)
    x = score["array_col"].to_numpy(dtype=float)
    y = score["array_row"].to_numpy(dtype=float)

    fig, axes = plt.subplots(1, 3, figsize=(14.8, 4.8), constrained_layout=True)
    sc_cm = axes[0].scatter(
        x, y, c=cm_pct, cmap=CM_PERCENTILE_CMAP, vmin=0, vmax=1,
        s=size, marker="o", linewidths=0, rasterized=True,
    )
    set_spatial_axes(axes[0], x, y)
    axes[0].set_title(f"CM percentile: {cm}")
    fig.colorbar(sc_cm, ax=axes[0], fraction=0.046, pad=0.02).set_label("CMact percentile")

    sc_epi = axes[1].scatter(
        x, y, c=epi_pct, cmap=EPI_PERCENTILE_CMAP, vmin=0, vmax=1,
        s=size, marker="o", linewidths=0, rasterized=True,
    )
    set_spatial_axes(axes[1], x, y)
    axes[1].set_title(f"Epi percentile: {epi}")
    fig.colorbar(sc_epi, ax=axes[1], fraction=0.046, pad=0.02).set_label("EPIfrac percentile")

    for cat in CATEGORY_ORDER:
        mask = category.eq(cat).to_numpy()
        if mask.any():
            axes[2].scatter(
                x[mask], y[mask], c=CATEGORY_COLORS[cat], s=size,
                marker="o", linewidths=0, rasterized=True,
                label=CATEGORY_LABELS[cat],
            )
    set_spatial_axes(axes[2], x, y)
    axes[2].set_title("Percentile quadrant")
    handles = [
        Line2D([0], [0], marker="o", color="w", label=CATEGORY_LABELS[cat],
               markerfacecolor=CATEGORY_COLORS[cat], markersize=6)
        for cat in CATEGORY_ORDER
    ]
    axes[2].legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        frameon=True,
        fontsize=6,
        borderpad=0.35,
        borderaxespad=0.0,
    )
    fig.suptitle(f"{sample}: {cm} | {epi}", y=1.04, fontsize=9, fontweight="bold")
    out_dir.mkdir(parents=True, exist_ok=True)
    save_pdf_svg(fig, out_dir / f"{sample}__{cm}__{epi}__percentile_quadrant_0p5")
```

Raw CMact/EPIfrac two-panel spatial figure:

```python
# REPLACEABLE: score table, sample, CM, epithelial subtype, output paths.
# FIXED: one figure with exactly two panels: raw CMact and raw EPIfrac.
from matplotlib.colors import Normalize

cm_raw = pd.to_numeric(score[f"CMact__{cm}"], errors="coerce")
epi_raw = pd.to_numeric(score[f"EPIfrac__{epi}"], errors="coerce")
x = score["array_col"].to_numpy(dtype=float)
y = score["array_row"].to_numpy(dtype=float)

fig, axes = plt.subplots(1, 2, figsize=(6.4, 3.2), constrained_layout=True)
sc0 = axes[0].scatter(
    x, y, c=cm_raw, s=size, cmap="magma",
    norm=Normalize(*robust_limits(cm_raw)), linewidths=0, rasterized=False,
)
axes[0].invert_yaxis()
axes[0].set_box_aspect(1)
axes[0].set_xlabel("array_col")
axes[0].set_ylabel("array_row")
axes[0].set_title(f"CMact: {cm}")
fig.colorbar(sc0, ax=axes[0], fraction=0.046, pad=0.02).set_label("CMact")

sc1 = axes[1].scatter(
    x, y, c=epi_raw, s=size, cmap="viridis",
    norm=Normalize(*robust_limits(epi_raw)), linewidths=0, rasterized=False,
)
axes[1].invert_yaxis()
axes[1].set_box_aspect(1)
axes[1].set_xlabel("array_col")
axes[1].set_ylabel("array_row")
axes[1].set_title(f"EPIfrac: {epi}")
fig.colorbar(sc1, ax=axes[1], fraction=0.046, pad=0.02).set_label("EPIfrac")

fig.suptitle(f"{sample}: {cm} | {epi}", y=1.04, fontsize=9, fontweight="bold")
fig.savefig(pdf_path, bbox_inches="tight", dpi=300)
fig.savefig(svg_path, bbox_inches="tight", dpi=300)
plt.close(fig)
```

Percentile/quadrant three-panel spatial figure:

```python
# REPLACEABLE: score table, sample, CM, epithelial subtype, output paths.
# FIXED: one figure with exactly three panels: CM percentile, Epi percentile,
# and four-color percentile quadrant. Cutoff is 0.5 for both percentiles.
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

CM_PERCENTILE_CMAP = LinearSegmentedColormap.from_list(
    "cm_percentile_light_to_blue",
    ["#f2f2f2", "#b7d7ea", "#4c78a8", "#2b6cb0"],
)
EPI_PERCENTILE_CMAP = LinearSegmentedColormap.from_list(
    "epi_percentile_light_to_orange",
    ["#f2f2f2", "#fde6b3", "#fbbf24", "#f59e0b"],
)
CATEGORY_ORDER = ["both_low", "epi_high_only", "cm_high_only", "both_high"]
CATEGORY_LABELS = {
    "both_low": "CM<0.5, Epi<0.5",
    "epi_high_only": "Epi>=0.5 only",
    "cm_high_only": "CM>=0.5 only",
    "both_high": "Both>=0.5",
}
CATEGORY_COLORS = {
    "both_low": "#d9d9d9",
    "epi_high_only": "#f59e0b",
    "cm_high_only": "#2b6cb0",
    "both_high": "#b91c1c",
}

cm_pct = percentile_rank(score[f"CMact__{cm}"])
epi_pct = percentile_rank(score[f"EPIfrac__{epi}"])
category = classify_percentiles(cm_pct, epi_pct)

fig, axes = plt.subplots(1, 3, figsize=(14.8, 4.8), constrained_layout=True)
sc_cm = axes[0].scatter(
    x, y, c=cm_pct, cmap=CM_PERCENTILE_CMAP, vmin=0, vmax=1,
    s=size, marker="o", linewidths=0, rasterized=True,
)
set_spatial_axes(axes[0], x, y)
axes[0].set_box_aspect(1)
axes[0].set_title(f"CM percentile: {cm}")
fig.colorbar(sc_cm, ax=axes[0], fraction=0.046, pad=0.02).set_label("CMact percentile")

sc_epi = axes[1].scatter(
    x, y, c=epi_pct, cmap=EPI_PERCENTILE_CMAP, vmin=0, vmax=1,
    s=size, marker="o", linewidths=0, rasterized=True,
)
set_spatial_axes(axes[1], x, y)
axes[1].set_box_aspect(1)
axes[1].set_title(f"Epi percentile: {epi}")
fig.colorbar(sc_epi, ax=axes[1], fraction=0.046, pad=0.02).set_label("EPIfrac percentile")

for cat in CATEGORY_ORDER:
    mask = category.eq(cat).to_numpy()
    if mask.any():
        axes[2].scatter(
            x[mask], y[mask], c=CATEGORY_COLORS[cat], s=size,
            marker="o", linewidths=0, rasterized=True, label=CATEGORY_LABELS[cat],
        )
set_spatial_axes(axes[2], x, y)
axes[2].set_box_aspect(1)
axes[2].set_title("Percentile quadrant")
handles = [
    Line2D([0], [0], marker="o", color="w", label=CATEGORY_LABELS[cat],
           markerfacecolor=CATEGORY_COLORS[cat], markersize=6)
    for cat in CATEGORY_ORDER
]
axes[2].legend(
    handles=handles,
    loc="upper left",
    bbox_to_anchor=(1.02, 1.0),
    frameon=True,
    fontsize=6,
    borderpad=0.35,
    borderaxespad=0.0,
)
# REQUIRED: keep the four-color legend at the upper-right side of the quadrant
# panel and verify it does not overlap spots, colorbars, titles, tick labels,
# or adjacent panels.
# Every panel must remain square. If the legend overlaps or is clipped, increase
# figsize and/or bbox_to_anchor x and rerender; do not shrink or stretch panels.
fig.suptitle(f"{sample}: {cm} | {epi}", y=1.04, fontsize=9, fontweight="bold")
fig.savefig(pdf_path, bbox_inches="tight", dpi=300)
fig.savefig(svg_path, bbox_inches="tight", dpi=300)
plt.close(fig)
```

Stouffer signed-Z heatmaps:

```python
# REPLACEABLE: METHOD_DIR, TABLE_DIR, FIGURE_DIR, and the input CSV names if
# the active project uses non-default file stems.
# FIXED: sample-level Fisher/Stouffer heatmaps use combined_signed_z for color,
# BH q-values only for ns/*/**/*** text, no numeric cell values, coolwarm centered
# at 0, square cells, PDF and SVG output.
from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

METHOD_DIR = Path.cwd()
TABLE_DIR = METHOD_DIR / "tables"
FIGURE_DIR = METHOD_DIR / "figures" / "stouffer_heatmaps"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "text.usetex": False,
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.titlesize": 9,
})

def cm_sort_key(value: str) -> tuple[int, str]:
    match = re.search(r"joint_(\d+)", str(value))
    return (int(match.group(1)) if match else 999, str(value))

def q_label(q: float) -> str:
    if not np.isfinite(q):
        return "ns"
    if q < 0.001:
        return "***"
    if q < 0.01:
        return "**"
    if q < 0.05:
        return "*"
    return "ns"

def ordered_stouffer_matrix(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    z_value = df.pivot(index="epi_subtype", columns="CM", values="combined_signed_z")
    q_value = df.pivot(index="epi_subtype", columns="CM", values="q_value_bh")
    rows = sorted(z_value.index)
    cols = sorted(z_value.columns, key=cm_sort_key)
    return z_value.reindex(index=rows, columns=cols), q_value.reindex(index=rows, columns=cols)

def finite_symmetric_limit(values: pd.DataFrame) -> float:
    arr = values.to_numpy(dtype=float)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return 1.0
    vmax = float(np.nanmax(np.abs(finite)))
    return vmax if vmax > 0 else 1.0

def q_annot(value_matrix: pd.DataFrame, q_matrix: pd.DataFrame) -> pd.DataFrame:
    annot = pd.DataFrame("", index=value_matrix.index, columns=value_matrix.columns)
    for row in annot.index:
        for col in annot.columns:
            annot.loc[row, col] = q_label(q_matrix.loc[row, col])
    return annot

def plot_sample_stouffer_signed_z(input_csv: Path, output_stem: str, title: str) -> None:
    df = pd.read_csv(input_csv)
    z_matrix, q_matrix = ordered_stouffer_matrix(df)
    vmax = finite_symmetric_limit(z_matrix)
    fig_w = max(7.8, 0.52 * z_matrix.shape[1] + 2.6)
    fig_h = max(4.8, 0.42 * z_matrix.shape[0] + 1.9)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(
        z_matrix,
        cmap="coolwarm",
        center=0,
        vmin=-vmax,
        vmax=vmax,
        linewidths=0.25,
        linecolor="white",
        annot=q_annot(z_matrix, q_matrix),
        fmt="",
        annot_kws={"fontsize": 6.5, "color": "black", "linespacing": 0.9},
        cbar_kws={"label": "combined signed Z"},
        square=True,
        ax=ax,
    )
    ax.set_title(title)
    ax.set_xlabel("CM")
    ax.set_ylabel("Epithelial subtype")
    ax.tick_params(axis="x", labelrotation=90)
    ax.tick_params(axis="y", labelrotation=0)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / f"{output_stem}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(FIGURE_DIR / f"{output_stem}.svg", bbox_inches="tight", dpi=300)
    plt.close(fig)

# REQUIRED project-compatible outputs when reproducing the method16 branch:
# 14 samples:
plot_sample_stouffer_signed_z(
    TABLE_DIR / "method16_percentile_quadrant_fisher_sample_stouffer_14samples.csv",
    "method16_percentile_quadrant_fisher_sample_stouffer_14samples_signedZ_qstars",
    "Method16 percentile quadrant Fisher Stouffer (14 samples)",
)
# 12 tumor samples, excluding R_med and R_cor:
plot_sample_stouffer_signed_z(
    TABLE_DIR / "method16_percentile_quadrant_fisher_sample_stouffer_12samples_no_rmed_rcor.csv",
    "method16_percentile_quadrant_fisher_sample_stouffer_12samples_no_rmed_rcor_signedZ_qstars",
    "Method16 percentile quadrant Fisher Stouffer (12 samples, no rmed/rcor)",
)
```

## Validation

- Confirm `.X` is normalized/log expression before marker extraction.
- Confirm spot IDs and sample names match downstream CSVs.
- Confirm coordinates are real spatial coordinates.
- Confirm Fisher/Stouffer heatmaps and summary tables use the intended sample
  scope: all-samples versus tumor-only.
- Confirm tumor-only sample membership comes from an explicit sample-status
  table or user-approved mapping, not filename guessing.
- Do not mix outputs from different spatial mapping methods unless the user explicitly requests a method comparison.
- Do not expose internal method numbers in final captions, reusable skill text, or publication-facing output names.

## Completion Checklist

A compact Epi-CM run is complete only when:

```text
01-celltype_integration_clustering produced adata_merge, adata_qc, adata_harmony, and selected adata_inte
02-cell_subtype_integration_clustering produced adata_anno, adata_epi, and adata_anno_cellsubtype
03-epi-cm-discovery produced CM matrices, CM classification, all-pair Spearman tables, node-set CSVs, and final figures; require all-pair Pearson tables only when the user explicitly requested the optional Pearson branch
04-spatial-validation-optional is either completed with per-sample and combined outputs, or explicitly recorded as skipped because no spatial input was requested/provided
all executed blocks have code, readme/provenance, parameters, package versions, and PDF/SVG-only figures
```
