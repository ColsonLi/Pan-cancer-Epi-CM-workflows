# Figure 2-2: epithelial subtype GO enrichment and selected-term dotplots

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(grid)
})

# Part 1: GO enrichment and per-subtype GO BP table export
setwd('/mnt/disk18t/data_disk18t/new_kirc_adata/cellsubtype_degs/')


read_and_filter_top1000 <- function(file){

  df <- read.csv(file, row.names = 1, check.names = FALSE)
  df_filtered <- df[df$logfoldchanges > 0.25 & df$pvals_adj < 0.05, ]
  syms <- rownames(df_filtered)[1:min(1000, nrow(df_filtered))]
  syms <- syms[!is.na(syms) & syms != ""]
  

  unique(syms)
}

syms_aldob <- read_and_filter_top1000('Epi_ALDOB_degs_epi.csv')
syms_aqp2  <- read_and_filter_top1000('Epi_AQP2_degs_epi.csv')
syms_ca12  <- read_and_filter_top1000('Epi_CA12_degs_epi.csv')
syms_ca9   <- read_and_filter_top1000('Epi_CA9_degs_epi.csv')
syms_gpx3  <- read_and_filter_top1000('Epi_GPX3_degs_epi.csv')
syms_jun   <- read_and_filter_top1000('Epi_JUN_degs_epi.csv')
syms_miox  <- read_and_filter_top1000('Epi_MIOX_degs_epi.csv')
syms_vim   <- read_and_filter_top1000('Epi_VIM_degs_epi.csv')


to_entrez <- function(symbols){
  bitr(symbols,
       fromType = "SYMBOL",
       toType   = "ENTREZID",
       OrgDb    = org.Hs.eg.db) |>
    dplyr::pull(ENTREZID) |>
    unique() |>
    as.character()
}

ent_aldob <- to_entrez(syms_aldob)
ent_aqp2  <- to_entrez(syms_aqp2)
ent_ca12  <- to_entrez(syms_ca12)
ent_ca9   <- to_entrez(syms_ca9)
ent_gpx3  <- to_entrez(syms_gpx3)
ent_jun   <- to_entrez(syms_jun)
ent_miox  <- to_entrez(syms_miox)
ent_vim   <- to_entrez(syms_vim)

run_go_all <- function(entrez_vec){
  enrichGO(gene          = entrez_vec,
           OrgDb         = org.Hs.eg.db,
           keyType       = "ENTREZID",
           ont           = "ALL",
           pAdjustMethod = "BH",
           pvalueCutoff  = 0.05,
           qvalueCutoff  = 0.05,
           readable      = TRUE)
}

ego_aldob <- run_go_all(ent_aldob)
ego_aqp2  <- run_go_all(ent_aqp2)
ego_ca12  <- run_go_all(ent_ca12)
ego_ca9   <- run_go_all(ent_ca9)
ego_gpx3  <- run_go_all(ent_gpx3)
ego_jun   <- run_go_all(ent_jun)
ego_miox  <- run_go_all(ent_miox)
ego_vim   <- run_go_all(ent_vim)


df_ALDOB <- as.data.frame(ego_aldob) |> mutate(Group = "Epi_ALDOB")
df_AQP2  <- as.data.frame(ego_aqp2)  |> mutate(Group = "Epi_AQP2")
df_CA12  <- as.data.frame(ego_ca12)  |> mutate(Group = "Epi_CA12")
df_CA9   <- as.data.frame(ego_ca9)   |> mutate(Group = "Epi_CA9")
df_GPX3  <- as.data.frame(ego_gpx3)  |> mutate(Group = "Epi_GPX3")
df_JUN   <- as.data.frame(ego_jun)   |> mutate(Group = "Epi_JUN")
df_MIOX  <- as.data.frame(ego_miox)  |> mutate(Group = "Epi_MIOX")
df_VIM   <- as.data.frame(ego_vim)   |> mutate(Group = "Epi_VIM")


go_df_all <- bind_rows(
  df_ALDOB, df_AQP2, df_CA12, df_CA9, 
  df_GPX3, df_JUN, df_MIOX, df_VIM
)


keep_cols <- c("ID","Description","ONTOLOGY","GeneRatio","Count",
               "pvalue","p.adjust","qvalue","geneID","Group")

for(g in unique(go_df_all$Group)){
  out <- go_df_all %>%

    dplyr::filter(Group == g, ONTOLOGY == "BP") %>% 
    dplyr::arrange(p.adjust) %>%
    dplyr::select(dplyr::all_of(intersect(keep_cols, names(go_df_all))))
  

  if(nrow(out) > 0){
    write.csv(out, paste0(as.character(g), "_GO_BP_group_20260405.csv"), row.names = FALSE)
  }
}

# Part 2: selected GO BP dotplots
setwd('/mnt/disk18t/data_disk18t/new_kirc_adata/gobp_filtered_top1000genes/')

df_MIOX <- read.csv('./Epi_MIOX_GO_BP_group_20260405.csv')
df_GPX3 <- read.csv('./Epi_GPX3_GO_BP_group_20260405.csv')
df_ALDOB <- read.csv('./Epi_ALDOB_GO_BP_group_20260405.csv')
df_CA9 <- read.csv('./Epi_CA9_GO_BP_group_20260405.csv')
df_JUN <- read.csv('./Epi_JUN_GO_BP_group_20260405.csv')
df_VIM <- read.csv('./Epi_VIM_GO_BP_group_20260405.csv')
df_AQP2 <- read.csv('./Epi_AQP2_GO_BP_group_20260405.csv')
df_CA12 <- read.csv('./Epi_CA12_GO_BP_group_20260405.csv')

go_df_all <- bind_rows(
  df_ALDOB, df_AQP2, df_CA12, df_CA9, 
  df_GPX3, df_JUN, df_MIOX, df_VIM
) |>

  mutate(Group = factor(Group, levels = c(
    "Epi_MIOX","Epi_GPX3","Epi_ALDOB", "Epi_CA9", 
     "Epi_JUN",  "Epi_VIM", "Epi_AQP2", "Epi_CA12"
  )))


groups_list <- list(
  Epi_ALDOB = c("GO:0034384","GO:0006520","GO:0006119","GO:0009636","GO:0006641"),
  Epi_AQP2  = c("GO:0050891","GO:0055078","GO:0072044","GO:0072073","GO:0031589"),
  Epi_CA12  = c("GO:1902600","GO:0006885","GO:0051452","GO:0072659","GO:0002181"),
  Epi_CA9   = c("GO:0001666","GO:0000280","GO:0007059","GO:0048002","GO:0007094"),
  Epi_GPX3  = c("GO:0045454","GO:0098754","GO:0044282","GO:0022904","GO:0042743"),
  Epi_JUN   = c("GO:0034976","GO:0097193","GO:0036293","GO:0043484","GO:0006986"),
  Epi_MIOX  = c("GO:1903825","GO:0003333","GO:0006006","GO:0170039","GO:1901293"),
  Epi_VIM   = c("GO:0090497","GO:0007015","GO:0030198","GO:0150117","GO:0001954") 
)


check_results <- lapply(names(groups_list), function(g) {
  selected_ids <- groups_list[[g]]

  existing_ids <- go_df_all$ID[go_df_all$Group == g]
  
  
  missing <- setdiff(selected_ids, existing_ids)
  
  if(length(missing) > 0) {
    message(paste0("⚠️ Group [", g, "] is missing IDs: ", paste(missing, collapse = ", ")))
  } else {
    message(paste0("✅ Group [", g, "] passed validation for all IDs!"))
  }
  return(missing)
})


go_ids <- c(
  "GO:1903825","GO:0003333","GO:0006006","GO:0170039","GO:1901293", ### Epi_MIOX includes subtype-specific terms
  "GO:0045454","GO:0098754","GO:0044282","GO:0022904","GO:0042743", ###Epi_GPX3
  "GO:0006520","GO:0006119","GO:0009636","GO:0006641","GO:0034384", ###Epi_ALDOB
  "GO:0001666","GO:0000280","GO:0007059","GO:0048002","GO:0007094", ###Epi_CA9
  "GO:0034976","GO:0097193","GO:0036293","GO:0043484","GO:0006986", ###Epi_JUN
  "GO:0090497","GO:0007015","GO:0030198","GO:0150117","GO:0001954",  ###Epi_VIM
  "GO:0050891","GO:0055078","GO:0072044","GO:0072073","GO:0031589", ###Epi_AQP2
  "GO:1902600","GO:0006885","GO:0051452","GO:0072659","GO:0002181" ###Epi_CA12
)

df_ALDOB_filtered <- df_ALDOB %>% filter(ID %in% go_ids)
df_AQP2_filtered  <- df_AQP2  %>% filter(ID %in% go_ids)
df_CA12_filtered  <- df_CA12  %>% filter(ID %in% go_ids)
df_CA9_filtered   <- df_CA9   %>% filter(ID %in% go_ids)
df_GPX3_filtered  <- df_GPX3  %>% filter(ID %in% go_ids)
df_JUN_filtered   <- df_JUN   %>% filter(ID %in% go_ids)
df_MIOX_filtered  <- df_MIOX  %>% filter(ID %in% go_ids)
df_VIM_filtered   <- df_VIM   %>% filter(ID %in% go_ids)

go_df_all <- bind_rows(
  df_ALDOB_filtered,
  df_AQP2_filtered,
  df_CA12_filtered,
  df_CA9_filtered,
  df_GPX3_filtered,
  df_JUN_filtered,
  df_MIOX_filtered,
  df_VIM_filtered
)

plot_df <- go_df_all %>%
  mutate(Group = factor(Group, levels = c(
    "Epi_MIOX", 
    "Epi_GPX3",
    "Epi_ALDOB", 
    "Epi_CA9", 
    "Epi_JUN", 
    "Epi_VIM",
    "Epi_AQP2", 
    "Epi_CA12"
  )))
plot_df2 <- plot_df %>%
  mutate(
    log10p = -log10(p.adjust),
    log10p_capped = squish(log10p, range = c(5, 20))
  )
p <- ggplot(plot_df2,
            aes(x = Group,
                y = Description,
                size = Count,
                colour = log10p_capped)) +
  geom_point() +
  scale_colour_gradient(
    low = "blue", high = "red",
    #limits = c(5, 20),
    oob = squish
  ) +
  scale_size(range = c(1.2, 5.5)) +  
  guides(
    colour = guide_colorbar(order = 1, barheight = unit(6, "cm")),
    size   = guide_legend(order = 2, override.aes = list(colour = "black"))
  ) +
  labs(
    title  = "GO Biological Process (Selected terms)",
    x = "Group", y = "GO Term",
    colour = "-log10(p.adjust)"
  ) +
  theme_bw(base_family = "Arial") +
  theme(
    plot.title  = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey85"),
    legend.title = element_text(size = 12, face = "bold", color = "black"),
    legend.text  = element_text(size = 10, color = "black"),
    legend.position = "right",
    legend.key.size = unit(1.0, "lines")
  )

print(p)


group_levels <- c(
  "Epi_MIOX",
  "Epi_GPX3",
  "Epi_ALDOB",
  "Epi_CA9",
  "Epi_JUN",
  "Epi_VIM",
  "Epi_AQP2",
  "Epi_CA12"
)


term_map <- bind_rows(lapply(names(groups_list), function(g) {
  data.frame(
    Group_owner = g,
    ID = groups_list[[g]],
    term_order = seq_along(groups_list[[g]]),
    stringsAsFactors = FALSE
  )
})) %>%
  mutate(Group_owner = factor(Group_owner, levels = group_levels))


plot_df <- go_df_all %>%
  mutate(Group = factor(Group, levels = group_levels)) %>%
  left_join(term_map, by = "ID") %>%
  mutate(
    term_key = paste(ID, Description, sep = " | ")
  )


term_levels <- term_map %>%
  left_join(
    plot_df %>% distinct(ID, Description, term_key),
    by = "ID"
  ) %>%
  arrange(Group_owner, term_order) %>%
  pull(term_key)

plot_df2 <- plot_df %>%
  mutate(
    term_key = factor(term_key, levels = rev(unique(term_levels))),
    log10p = -log10(p.adjust),
    log10p_capped = squish(log10p, range = c(2, 20))
  )

p <- ggplot(
  plot_df2,
  aes(
    x = Group,
    y = term_key,
    size = Count,
    colour = log10p_capped
  )
) +
  geom_point() +
  scale_y_discrete(
    labels = function(x) sub("^GO:[^|]+ \\| ", "", x)
  ) +
  scale_colour_gradient(
    low = "blue", high = "red",
    oob = squish
  ) +
  scale_size(range = c(1.2, 5.5)) +
  guides(
    colour = guide_colorbar(order = 1, barheight = unit(6, "cm")),
    size   = guide_legend(order = 2, override.aes = list(colour = "black"))
  ) +
  labs(
    title  = "GO Biological Process (Selected terms)",
    x = "Group",
    y = "GO Term",
    colour = "-log10(p.adjust)"
  ) +
  theme_bw(base_family = "Arial") +
  theme(
    plot.title  = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey85"),
    legend.title = element_text(size = 12, face = "bold", color = "black"),
    legend.text  = element_text(size = 10, color = "black"),
    legend.position = "right",
    legend.key.size = unit(1.0, "lines")
  )

print(p)
