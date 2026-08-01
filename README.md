# Pan-cancer Epi-CM workflows

This repository contains the analysis code supporting the manuscript:

> **Coupled epithelial and microenvironmental states define spatial ecosystems in clear cell renal cell carcinoma**

It also provides a reusable agent skill for the core epithelial–cellular module
(Epi-CM) workflow. The skill covers single-cell preprocessing and annotation,
cell-subtype analysis, cellular-module discovery, Epi-CM association testing,
and optional spatial mapping and validation.

The manuscript and its biological conclusions focus on clear cell renal cell
carcinoma (ccRCC). The reusable workflow has additionally been technically
exercised on breast cancer (BRCA) and colorectal cancer (CRC) datasets to assess
cross-dataset portability. These applications should not be interpreted as
complete cross-cancer biological validation.

## Repository contents

```text
.
├── main_codes/                        # Manuscript analysis code organized by figure
│   ├── fig1/
│   ├── fig2/
│   ├── fig3/
│   ├── fig4_et_s8/
│   ├── fig5/
│   ├── fig6_et_s9/
│   ├── figs1/
│   ├── figs2/
│   ├── figs3/
│   ├── figs4/
│   ├── figs5/
│   ├── figs6/
│   └── figs7/
└── skills-project-epi-cm-core-workflow/ # Reusable Epi-CM agent skill
    ├── SKILL.md
    └── assets/
```

## Code organization

The analysis code is organized according to manuscript figures rather than the
original chronological analysis history.

- Files beginning with `fig` correspond to main figures.
- Files beginning with `figs` correspond to supplementary figures.
- `fig4_et_s8` contains analyses shared by Figure 4 and Supplementary Figure 8.
- `fig6_et_s9` contains analyses shared by Figure 6 and Supplementary Figure 9.
- The merged figure-level notebook or R script is the primary entry point in
  each directory. The additional component notebooks and scripts are retained
  to make the provenance of the merged code explicit.

Some figure-level workflows save and subsequently reload intermediate objects.
Therefore, the figure order is not necessarily a strict end-to-end execution
order. Required upstream objects are described in the corresponding code.

## Figure-to-code overview

| Code directory | Primary entry point | Main analysis content |
|---|---|---|
| `fig1` | `fig1-1.ipynb`, `fig1-2.ipynb` | Single-cell data merging, quality control, integration, broad annotation, score-based annotation support, and Figure 1 summaries |
| `fig2` | `fig2-1.ipynb`, `fig2-2.R` | Epithelial subtype analysis, differential expression, DPT analysis, DEG-expression visualization, and GO enrichment |
| `fig3` | `fig3-1.ipynb` | Group-balanced joint NMF, cellular-module definition, Epi-CM association analysis, and selected CM visualizations |
| `fig4_et_s8` | `fig4_et_s8-1.ipynb` | Tangram-based spatial mapping, spatial CM–epithelial association testing, Fisher/Stouffer integration, and spatial figures |
| `fig5` | `fig5-1.ipynb`, `fig5-2.ipynb` | Epi_VIM versus Epi_JUN pySCENIC analysis and selected CellChat analyses |
| `fig6_et_s9` | `fig6_s9-1.ipynb` | TCGA-KIRC and JAVELIN ssGSEA, Cox models, Kaplan–Meier analyses, and treatment-arm analyses |
| `figs1` | `figs1.ipynb` | Quality-control and integration diagnostics |
| `figs2` | `figs2-1.ipynb`, `figs2-2.ipynb`, `figs2-3.ipynb` | inferCNV, epithelial marker visualization, RNA velocity, veloVI, and CellRank analyses |
| `figs3` | `figs3-1.ipynb`, `figs3-2.R`, `figs3-3.ipynb` | DPT/stemness visualization, Monocle3 trajectories, and TCGA-KIRC epithelial survival analysis |
| `figs4` | `figs4-1.ipynb`, `figs4-2.ipynb` | Epi_CA9 reclustering, marker visualization, velocity, and CellRank analysis |
| `figs5` | `figs5-1.ipynb`, `figs5-2.ipynb`, `figs5-3.ipynb` | Non-epithelial subtype analysis, AutoGeneS/NuSVR deconvolution, CM activity summaries, and NMF diagnostics |
| `figs6` | `figs6-1.ipynb` | Tumor CM network outputs |
| `figs7` | `figs7-1.ipynb` | Normal-like CM network outputs |

## Figure-prefixed script guide

Only scripts whose file names begin with `fig` are documented below. These 23
public entry points comprise 21 merged Jupyter notebooks and two merged R
scripts. Inputs and outputs are summarized by data type rather than local file
path. The retained component files are source provenance and are not listed as
separate workflows.

### `fig1`

| Script | Main input | Main output |
|---|---|---|
| `fig1-1.ipynb` | Raw single-cell count matrices with sample annotations | Merged, QC-filtered, Harmony-integrated, broadly annotated and score-ranked AnnData objects; broad-cell DEG tables; QC violin and Leiden UMAP plots |
| `fig1-2.ipynb` | Score-ranked all-cell AnnData and lineage-specific subtype AnnData objects | Final all-cell subtype-annotated AnnData; annotation UMAPs, marker dot plot, cell-composition ring charts and box plot |

### `fig2`

| Script | Main input | Main output |
|---|---|---|
| `fig2-1.ipynb` | QC/annotated single-cell AnnData, epithelial cells, and spatial transcriptomic AnnData objects | Epithelial subtype DEG tables and AnnData; MIOX-rooted DPT/PAGA object; subtype, stemness, trajectory and diffusion-map plots; spatial gene/score maps |
| `fig2-2.R` | Epithelial-subtype DEG tables | GO enrichment result tables and selected epithelial-subtype GO dot plots |

### `fig3`

| Script | Main input | Main output |
|---|---|---|
| `fig3-1.ipynb` | Final all-cell subtype-annotated AnnData | Group-balanced joint-NMF matrices, CM classification, node/edge tables and all-pair Epi-CM Spearman results; association/loading heatmaps, tumor CM node plot, scatter plots and annotated CM activity clustermap |

### `fig4_et_s8`

| Script | Main input | Main output |
|---|---|---|
| `fig4_et_s8-1.ipynb` | Annotated single-cell reference, subtype DEG tables, CM loadings/statistics and spatial transcriptomic AnnData objects | Tangram reference and spot scores; per-sample and combined Spearman/Fisher/Stouffer statistics; CM–epithelial heatmaps, spatial pair maps and core-pair forest plot |

### `fig5`

| Script | Main input | Main output |
|---|---|---|
| `fig5-1.ipynb` | Raw-count all-cell AnnData and epithelial subtype annotations | Filtered loom inputs, pySCENIC/AUCell results and per-cell regulon comparison table; ranked Epi_VIM-versus-Epi_JUN regulon-difference plot |
| `fig5-2.ipynb` | Final subtype-annotated AnnData and CM-to-subtype mapping information | CellChat object, ligand–receptor result tables and selected bubble-plot data; communication circle plot and Epi_JUN-to-CM1 bubble plot |

### `fig6_et_s9`

| Script | Main input | Main output |
|---|---|---|
| `fig6_s9-1.ipynb` | Epithelial/CM DEG signatures, TCGA-KIRC and JAVELIN expression matrices, and clinical survival/treatment data | Epithelial and CM ssGSEA score matrices, Cox/Kaplan–Meier summary tables, Cox and CM-z forest plots, and TCGA-KIRC/JAVELIN survival and treatment-arm Kaplan–Meier figures |

### `figs1`

| Script | Main input | Main output |
|---|---|---|
| `figs1.ipynb` | QC-filtered and Harmony-integrated AnnData with sample status and study-series metadata | Normal-like/tumor composition box plots and study-series UMAP diagnostics before and after Harmony integration |

### `figs2`

| Script | Main input | Main output |
|---|---|---|
| `figs2-1.ipynb` | Raw-count/QC AnnData and annotated epithelial AnnData | Epithelial inferCNV result object; chromosome heatmap, all/tumor/normal-like CNV UMAP panels and epithelial marker dot plots |
| `figs2-2.ipynb` | Tumor-sample velocity objects and epithelial DPT/UMAP annotations | Merged-tumor veloVI/CellRank objects and transition matrices; directed-transition graph and matrix heatmap |
| `figs2-3.ipynb` | Processed epithelial velocity/DPT AnnData | Four-kernel CellRank object, state/fate/lineage-driver tables and Epi_JUN lineage-driver/trend gene heatmap |

### `figs3`

| Script | Main input | Main output |
|---|---|---|
| `figs3-1.ipynb` | Saved MIOX-rooted epithelial DPT AnnData | Stemness-score UMAP, stemness violin plot and epithelial diffusion-map plot |
| `figs3-2.R` | Six-subtype epithelial raw-count matrix, cell/gene metadata, UMAP coordinates and root cell | Monocle3 trajectory objects and pseudotime tables; subtype and pseudotime trajectory figures |
| `figs3-3.ipynb` | Epithelial DEG signatures, TCGA-KIRC expression matrix and overall-survival data | Epithelial-state ssGSEA scores, survival summaries and TCGA-KIRC Kaplan–Meier curves |

### `figs4`

| Script | Main input | Main output |
|---|---|---|
| `figs4-1.ipynb` | QC/annotated single-cell AnnData and Epi_CA9 cells | Epi_CA9 reclustered AnnData and DEG results; reclustering UMAP and representative-marker dot plot |
| `figs4-2.ipynb` | Eight tumor-sample loom/velocity datasets and Epi_CA9 annotations | Per-sample and merged Epi_CA9 velocity/CellRank objects and matrices; final single-direction CellRank graph |

### `figs5`

| Script | Main input | Main output |
|---|---|---|
| `figs5-1.ipynb` | QC and score-annotated all-cell AnnData | Lineage-specific subtype AnnData and DEG tables; publication-style subtype UMAP panels |
| `figs5-2.ipynb` | Subtype mean-expression reference and TCGA-KIRC bulk expression matrix | AutoGeneS/NuSVR subtype proportions and non-epithelial correlation matrices; proportion and Spearman-correlation clustermaps |
| `figs5-3.ipynb` | Joint-NMF loading/activity, sample-status, node-set and rank-selection tables | Normal-like CM node plot, tumor-versus-normal-like CM activity bar plot and joint-NMF rank-selection plot |

### `figs6` and `figs7`

| Script | Main input | Main output |
|---|---|---|
| `figs6/figs6-1.ipynb` | Tumor CM node/edge tables, subtype-frequency matrix and sample-status table | Tumor-centric CM node plot and tumor top-node correlation heatmap |
| `figs7/figs7-1.ipynb` | Normal-like top-node table, subtype-frequency matrix and sample-status table | Normal-like top-node correlation heatmap |

## Running the code

The code was assembled from the analysis environment used for the manuscript.
Before execution:

1. Inspect the first cells or configuration section of the relevant combined
   notebook or script.
2. Replace local absolute paths with paths appropriate for the target system.
3. Provide the required public or controlled-access input datasets and upstream
   intermediate objects.
4. Create isolated Python and R environments containing the package versions
   required by the selected workflow.
5. Run the analysis blocks before their corresponding plotting blocks.

The notebooks are provided with saved cell outputs removed. Large input data,
intermediate AnnData/loom/RDS objects, trained models, and generated figure
files are not necessarily included in the repository.

## Main software families

Different figure workflows use different subsets of the following tools:

- Python: Scanpy, AnnData, scVelo, veloVI/scvi-tools, CellRank, Tangram,
  infercnvpy, AutoGeneS, SciPy, statsmodels, pandas, seaborn, and matplotlib.
- R: Monocle3, CellChat, clusterProfiler, survival, survminer, and associated
  Bioconductor annotation packages.
- Regulatory-network analysis: pySCENIC and AUCell.

Exact dependencies should be taken from the relevant notebook/script and its
recorded analysis environment.

## Reproducibility notes

- Combined figure-level files are the recommended public entry points.
- Statistical source tables should be regenerated before figures when the
  upstream data or parameters change.
- Corrected p or q values, rather than raw p values, are used for significance
  labels where multiple testing is performed.
- Spearman is the default single-cell Epi-CM association method; Pearson is not
  required unless the optional Pearson branch is explicitly requested.
- The spatial workflow uses every available reference cell per subtype and does
  not apply automatic per-subtype cell caps or random subsampling.
- BRCA and CRC usage demonstrates technical workflow portability; interpretation
  of cancer-specific biological findings requires dedicated validation.

## Data availability

The repository primarily distributes analysis code and the reusable workflow
skill. Input data and large intermediate artifacts should be obtained from the
repositories or controlled-access resources described in the manuscript and
its data-availability statement. Users are responsible for complying with the
corresponding data-use terms.

## Citation

If you use this repository, please cite the associated manuscript:

> *Coupled epithelial and microenvironmental states define spatial ecosystems
> in clear cell renal cell carcinoma.*

Full author, journal, year, DOI, and preprint information can be added here when
available.

## License

No license is assigned in this draft README. Add an explicit repository license
before public release and verify that redistribution of all bundled code and
assets is compatible with the licenses of the upstream software and datasets.

## Reusable Epi-CM agent skill

The reusable agent workflow is provided in
`skills-project-epi-cm-core-workflow/SKILL.md`. To start the core workflow, provide the
agent with a raw single-cell count matrix; spatial inputs are needed only when
the optional spatial-validation block is requested.

The skill guides the agent through:

1. Single-cell data merging, metadata harmonization, quality control, doublet
   assessment, normalization, Harmony integration, clustering, and broad cell
   annotation.
2. Lineage-specific subclustering, subtype differential-expression analysis,
   subtype annotation, and projection of subtype labels to the complete object.
3. Construction of epithelial and non-epithelial subtype-frequency matrices,
   group-balanced joint NMF, cellular-module classification, node and edge
   derivation, and Epi-CM association analysis.
4. Spearman correlation as the required default Epi-CM association branch.
   Pearson correlation is a separate optional branch and is run only when
   explicitly requested.
5. Optional spatial validation using post-annotation subtype DEG markers,
   Tangram mapping, per-sample spatial statistics, cross-sample integration,
   multiple-testing correction, and publication-oriented visualization.
6. Explicit output contracts, provenance records, parameter reports, validation
   checks, and figure manifests.

For the canonical spatial marker-selection branch, the skill requires
`pvals_adj < 0.05` and `logfoldchanges > 0`, retains stable descending `score`
order, removes duplicated gene names, and selects up to 100 genes per subtype.
