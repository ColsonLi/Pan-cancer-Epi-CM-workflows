# Part 1: Build the six-subtype epithelial Monocle3 trajectory.
suppressPackageStartupMessages({
  library(monocle3)
  library(Matrix)
  library(ggplot2)
  library(ggrastr)
})

base_dir <- normalizePath(".")
input_dir <- file.path(base_dir, "inputs")
out_dir <- file.path(base_dir, "outputs")
fig_dir <- file.path(base_dir, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

plot_and_save <- function(p, filename, width = 5.5, height = 4.5) {
  ggsave(
    file.path(fig_dir, paste0(filename, ".pdf")),
    p,
    width = width,
    height = height,
    device = cairo_pdf
  )
  ggsave(
    file.path(fig_dir, paste0(filename, ".svg")),
    p,
    width = width,
    height = height
  )
}

matrix_path <- file.path(
  input_dir,
  "epi_6subtype_raw_allgenes_genes_by_cells.mtx"
)
cell_meta_path <- file.path(input_dir, "cell_metadata.tsv")
gene_meta_path <- file.path(input_dir, "gene_metadata.tsv")
umap_path <- file.path(input_dir, "umap_coordinates.tsv")
root_path <- file.path(input_dir, "root_cell.tsv")

expr <- readMM(matrix_path)
cell_meta <- read.delim(
  cell_meta_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
gene_meta <- read.delim(
  gene_meta_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
umap <- read.delim(
  umap_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
root_info <- read.delim(
  root_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

rownames(gene_meta) <- gene_meta$gene_id
rownames(cell_meta) <- cell_meta$cell
rownames(umap) <- umap$cell

if (!all(colnames(expr) == NULL)) {
  colnames(expr) <- NULL
}
rownames(expr) <- gene_meta$gene_id
colnames(expr) <- cell_meta$cell

if (!identical(rownames(cell_meta), rownames(umap))) {
  stop("Cell metadata and UMAP coordinates are not aligned")
}
if (!root_info$root_cell[1] %in% rownames(cell_meta)) {
  stop("Root cell is not present in cell metadata")
}

cell_meta$cell_subtype <- factor(
  cell_meta$cell_subtype,
  levels = c(
    "Epi_MIOX",
    "Epi_ALDOB",
    "Epi_GPX3",
    "Epi_CA9",
    "Epi_JUN",
    "Epi_VIM"
  )
)

cds <- new_cell_data_set(
  expr,
  cell_metadata = cell_meta,
  gene_metadata = gene_meta
)
cds <- preprocess_cds(cds, num_dim = 50, norm_method = "none")
reducedDims(cds)$UMAP <- as.matrix(umap[, c("UMAP_1", "UMAP_2")])
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = FALSE)
cds <- order_cells(cds, root_cells = root_info$root_cell[1])

pseudotime_df <- data.frame(
  cell = colnames(cds),
  monocle3_pseudotime = as.numeric(pseudotime(cds)),
  cell_subtype = colData(cds)$cell_subtype,
  status = colData(cds)$status,
  sample = colData(cds)$sample,
  series = colData(cds)$series,
  stemness_score = colData(cds)$stemness_score,
  stringsAsFactors = FALSE
)
write.csv(
  pseudotime_df,
  file.path(out_dir, "epi_6subtype_monocle3_pseudotime.csv"),
  row.names = FALSE
)

summary_df <- aggregate(
  monocle3_pseudotime ~ cell_subtype,
  data = pseudotime_df,
  FUN = function(x) {
    c(
      n = sum(is.finite(x)),
      mean = mean(x, na.rm = TRUE),
      median = median(x, na.rm = TRUE)
    )
  }
)
summary_out <- data.frame(
  cell_subtype = summary_df$cell_subtype,
  n = summary_df$monocle3_pseudotime[, "n"],
  mean = summary_df$monocle3_pseudotime[, "mean"],
  median = summary_df$monocle3_pseudotime[, "median"]
)
write.csv(
  summary_out,
  file.path(out_dir, "epi_6subtype_monocle3_pseudotime_by_subtype.csv"),
  row.names = FALSE
)

saveRDS(
  cds,
  file.path(out_dir, "epi_6subtype_monocle3_root_miox.cds.rds")
)
writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo.txt")
)

p_subtype <- plot_cells(
  cds,
  color_cells_by = "cell_subtype",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  show_trajectory_graph = TRUE
) + ggtitle("Monocle3 trajectory: epithelial subtypes")
plot_and_save(
  p_subtype,
  "epi_6subtype_monocle3_trajectory_cell_subtype"
)

p_pt <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  show_trajectory_graph = TRUE
) + ggtitle("Monocle3 pseudotime")
plot_and_save(p_pt, "epi_6subtype_monocle3_trajectory_pseudotime")

p_stem <- plot_cells(
  cds,
  color_cells_by = "stemness_score",
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  show_trajectory_graph = TRUE
) + ggtitle("Stemness score")
plot_and_save(p_stem, "epi_6subtype_monocle3_trajectory_stemness_score")

p_violin <- ggplot(
  pseudotime_df,
  aes(
    x = cell_subtype,
    y = monocle3_pseudotime,
    fill = cell_subtype
  )
) +
  geom_violin(scale = "width", linewidth = 0.2, na.rm = TRUE) +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    legend.position = "none"
  ) +
  xlab(NULL) +
  ylab("Monocle3 pseudotime")
plot_and_save(
  p_violin,
  "epi_6subtype_monocle3_pseudotime_violin",
  width = 4,
  height = 3.5
)

cat("root_cell", root_info$root_cell[1], "\n")
cat("cells", ncol(cds), "genes", nrow(cds), "\n")
print(summary_out)


# Part 2: Relearn the trajectory with ncenter = 300.
cds <- readRDS(
  file.path(out_dir, "epi_6subtype_monocle3_root_miox.cds.rds")
)
root_cell <- root_info$root_cell[1]

if (!root_cell %in% colnames(cds)) {
  stop("Root cell is not present in the Monocle3 CDS")
}

cds <- learn_graph(
  cds,
  use_partition = FALSE,
  learn_graph_control = list(ncenter = 300)
)
cds <- order_cells(cds, root_cells = root_cell)

saveRDS(
  cds,
  file.path(
    out_dir,
    "epi_6subtype_monocle3_root_miox_pruned_ncenter300.cds.rds"
  )
)

cat("root_cell", root_cell, "\n")
cat("cells", ncol(cds), "genes", nrow(cds), "\n")
cat(
  "principal_graph_nodes",
  igraph::vcount(principal_graph(cds)[["UMAP"]]),
  "principal_graph_edges",
  igraph::ecount(principal_graph(cds)[["UMAP"]]),
  "\n"
)


# Part 3: Export ncenter = 300 pseudotime and draw the basic plot.
cds <- readRDS(
  file.path(
    out_dir,
    "epi_6subtype_monocle3_root_miox_pruned_ncenter300.cds.rds"
  )
)

pseudotime_df <- as.data.frame(colData(cds))
pseudotime_df$cell <- rownames(pseudotime_df)
pseudotime_df$monocle3_pseudotime <- pseudotime(cds)
pseudotime_df <- pseudotime_df[, c(
  "cell",
  "cell_subtype",
  "sample",
  "status",
  "stemness_score",
  "monocle3_pseudotime"
)]
write.csv(
  pseudotime_df,
  file.path(
    out_dir,
    "epi_6subtype_monocle3_pseudotime_pruned_ncenter300.csv"
  ),
  row.names = FALSE
)

common_args <- list(
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  label_roots = FALSE,
  show_trajectory_graph = TRUE,
  trajectory_graph_segment_size = 0.6,
  cell_size = 0.08,
  cell_stroke = 0,
  alpha = 0.75
)

p_pt <- do.call(
  plot_cells,
  c(list(cds = cds, color_cells_by = "pseudotime"), common_args)
) + ggtitle("Monocle3 pseudotime: ncenter=300")
plot_and_save(
  p_pt,
  "epi_6subtype_monocle3_trajectory_pseudotime_pruned_ncenter300"
)


# Part 4: Draw the final rasterized subtype and pseudotime panels.
color_path <- file.path(
  input_dir,
  "adata_epi_cell_subtype_colors.tsv"
)
cell_subtype_colors <- NULL
if (file.exists(color_path)) {
  color_df <- read.delim(
    color_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  cell_subtype_colors <- setNames(
    color_df$color,
    color_df$cell_subtype
  )
}

common_args <- list(
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  label_roots = FALSE,
  show_trajectory_graph = TRUE,
  trajectory_graph_color = "black",
  trajectory_graph_segment_size = 0.6,
  cell_size = 0.3,
  cell_stroke = 0,
  alpha = 1
)

square_panel_theme <- theme(
  aspect.ratio = 1,
  plot.title = element_text(size = 10),
  axis.title = element_text(size = 9),
  axis.text = element_text(size = 8),
  legend.title = element_text(size = 9),
  legend.text = element_text(size = 8)
)

p_subtype <- do.call(
  plot_cells,
  c(list(cds = cds, color_cells_by = "cell_subtype"), common_args)
) +
  ggtitle("Monocle3 trajectory: ncenter=300") +
  square_panel_theme
if (!is.null(cell_subtype_colors)) {
  p_subtype <- p_subtype + scale_color_manual(values = cell_subtype_colors)
}
p_subtype <- rasterise(
  p_subtype,
  layers = "Point",
  dpi = 600,
  dev = "ragg"
)
plot_and_save(
  p_subtype,
  paste0(
    "epi_6subtype_monocle3_trajectory_cell_subtype_",
    "pruned_ncenter300_rasterdots_square_panel_larger_",
    "size0p3_alpha1"
  ),
  width = 5.2,
  height = 4.0
)

p_pt <- do.call(
  plot_cells,
  c(list(cds = cds, color_cells_by = "pseudotime"), common_args)
) +
  ggtitle("Monocle3 pseudotime: ncenter=300") +
  square_panel_theme
p_pt <- rasterise(
  p_pt,
  layers = "Point",
  dpi = 600,
  dev = "ragg"
)
plot_and_save(
  p_pt,
  paste0(
    "epi_6subtype_monocle3_trajectory_pseudotime_",
    "pruned_ncenter300_rasterdots_square_panel_larger_",
    "size0p3_alpha1"
  ),
  width = 4.6,
  height = 4.0
)
