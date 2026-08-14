
# Load Packages -----------------------------------------------------------

library(tidyverse)
library(ggpointdensity)
library(viridis)
library(patchwork)
library(DESeq2)
library(GSVA)
library(rstatix)
library(ComplexHeatmap)
library(circlize)

dir.create(
  file.path(
    "output_data",
    "plots",
    "ssGSEA_and_reversion_magnitude"
  ), recursive = T
)

dir.create(
  file.path(
    "output_data",
    "results",
    "ssGSEA_and_reversion_magnitude"
  ), recursive = T
)

cat(
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
         "Started 07_ssGSEA_and_reversion_magnitudes\n"),
  file = "PROGRESS_LOG.md",
  append = TRUE
)


# Load Data ---------------------------------------------------------------

deseq2_hgps_umk57_vs_dmso <- read.csv(
  file.path(
    "output_data",
    "results",
    "deseq2",
    "deseq2_hgps_umk57_vs_dmso.csv"
  )
)

deseq2_ds_umk57_vs_dmso <- read.csv(
  file.path(
    "output_data",
    "results",
    "deseq2",
    "deseq2_ds_umk57_vs_dmso.csv"
  )
)

deseq2_hgps_vs_healthy <- read.csv(
  file.path(
    "output_data",
    "results",
    "deseq2",
    "deseq2_hgps_vs_healthy.csv"
  )
)

deseq2_ds_vs_healthy <- read.csv(
  file.path(
    "output_data",
    "results",
    "deseq2",
    "deseq2_ds_vs_healthy.csv"
  )
)
deseq2_ds_vs_healthy |> head()


# Correlation Plot --------------------------------------------------------

## HGPS -------

df_hgps_reversion <- inner_join(
  deseq2_hgps_vs_healthy |> 
    dplyr::select(gene, lfc_disease = log2FoldChange, padj_disease = padj),
  deseq2_hgps_umk57_vs_dmso |> 
    dplyr::select(gene, lfc_drug = log2FoldChange, padj_drug = padj),
  by = "gene"
)

# Spearman correlation across all genes.

cor_hgps <- cor.test(df_hgps_reversion$lfc_disease, 
                     df_hgps_reversion$lfc_drug, 
                     method = "spearman", exact = T)

spearman_hgps <- cor_hgps$estimate

pvalue_hgps   <- cor_hgps$p.value

p_label <- if (cor_hgps$p.value < 0.0001) {
  "P < 0.0001"
} else {
  paste0("p = ", signif(cor_hgps$p.value, 3))
}

p1 <- ggplot(df_hgps_reversion,
             aes(x = lfc_disease,
                 y = lfc_drug)) +
  geom_pointdensity(adjust = 1, alpha = 0.5) +
  viridis::scale_color_viridis() +
  theme_bw(base_size = 14) +
  xlim(-6, 14) +
  ylim(-10, 6) +
  annotate("text", x = min(df_hgps_reversion$lfc_disease, na.rm = TRUE), y = -9, hjust = 0,
           label = paste0("Rho = ", round(spearman_hgps, 2))) +
  annotate("text", x = min(df_hgps_reversion$lfc_disease, na.rm = TRUE), y = -10, hjust = 0,
           label = p_label) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_smooth(method = "lm", color = "orange") +
  labs(x = "\nLog2FC: HGPS vs Healthy", y = "Log2FC: HGPS UMK57 vs DMSO\n",
       title = "LOG2FC Correlation - HGPS (all genes)")

p1

## DS -------
  
df_ds_reversion <- inner_join(
  deseq2_ds_vs_healthy |> 
    dplyr::select(gene, lfc_disease = log2FoldChange, padj_disease = padj),
  deseq2_ds_umk57_vs_dmso |> 
    dplyr::select(gene, lfc_drug = log2FoldChange, padj_drug = padj),
  by = "gene"
)

# Spearman correlation across all genes.

cor_ds <- cor.test(df_ds_reversion$lfc_disease, 
                   df_ds_reversion$lfc_drug, 
                   method = "spearman", exact = T)
spearman_ds <- cor_ds$estimate
pvalue_ds   <- cor_ds$p.value
p_label_ds <- if (cor_ds$p.value < 0.0001) {
  "P < 0.0001"
} else {
  paste0("p = ", signif(cor_ds$p.value, 3))
}

p2 <- ggplot(df_ds_reversion,
             aes(x = lfc_disease,
                 y = lfc_drug)) +
  geom_pointdensity(adjust = 1, alpha = 0.5) +
  viridis::scale_color_viridis() +
  theme_bw(base_size = 14) +
  annotate("text", x = min(df_ds_reversion$lfc_disease, na.rm = TRUE), y = -9, hjust = 0,
           label = paste0("Rho = ", round(spearman_ds, 2))) +
  annotate("text", x = min(df_ds_reversion$lfc_disease, na.rm = TRUE), y = -10, hjust = 0,
           label = p_label_ds) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_smooth(method = "lm", color = "orange") +
  labs(x = "\nLog2FC: DS vs Healthy", y = "Log2FC: DS UMK57 vs DMSO\n",
       title = "LOG2FC Correlation - DS (all genes)")
p2


p1 <- ggrastr::rasterise(p1, layers = "Point", dpi = 300)
p1

p2 <- ggrastr::rasterise(p2, layers = "Point", dpi = 300)
p2

p1 + p2

ggsave(
  "correlation_log2fc_plot.pdf",
  path = file.path(
    "output_data",
    "plots",
    "ssGSEA_and_reversion_magnitude"
  ),
  width = 12,
  height = 4
)


# ssGSEA Hallmarks --------------------------------------------------------

raw_counts <- read.csv(
  file.path("input_data", "counts_matrix.csv"),
  row.names = 1
)

metadata <- read.csv(
  file.path("input_data", "sample_info.csv"),
  row.names = 1
)

colnames(raw_counts) == row.names(metadata)


# Create DESeq Object:

dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = metadata,
  design = ~ Individual)

colnames(dds)

dds <- estimateSizeFactors(dds)

smallestGroupSize <- 4

keep <- rowSums(counts(dds) >=5) >= smallestGroupSize 

# Recomended in vignette.

dds <- dds[keep,]

dds

# VST Transformation:

vst_mat <- assay(vst(dds, blind = TRUE))

dim(vst_mat)
head(vst_mat[, 1:3])


# Hall GS:

Hall <- msigdbr::msigdbr(
  species = "Homo sapiens",
  collection = "H"
)
Hall

hall_t2g <- Hall %>% dplyr::distinct(gs_name, gene_symbol) %>% 
  as.data.frame()

# GSVA wants a named list: gene set name -> vector of gene symbols:

hall_list <- split(hall_t2g$gene_symbol, hall_t2g$gs_name)



# GSVA - ssGSEA -----------------------------------------------------------

ssgsea_par <- ssgseaParam(
  exprData    = vst_mat,
  geneSets    = hall_list,
  minSize     = 5,
  maxSize     = 500,
  normalize   = TRUE  
)

ssgsea_scores <- gsva(ssgsea_par)

dim(ssgsea_scores)   # 50 pathways x n_samples
ssgsea_scores[1:5, 1:5]  


ssgsea_df <- as.data.frame(ssgsea_scores) |> 
  rownames_to_column("pathway") |> 
  pivot_longer(-pathway, names_to = "sample", values_to = "es_score") |> 
  left_join(metadata |> rownames_to_column("sample"), by = "sample") |> 
  dplyr::mutate(group = paste0(Disease, "_", Treatment))

ssgsea_df$group <- factor(
  ssgsea_df$group,
  levels = c(
    "Healthy_DMSO",
    "HGPS_DMSO",
    "HGPS_UMK57",
    "DS_DMSO",
    "DS_UMK57"
  )
)

# Heatmap DS

col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_fun(seq(-2, 2))

metadata_correct_order <- metadata |> 
  dplyr::filter(Disease %in% c("Healthy", "DS")) |> 
  dplyr::arrange(
    factor(Disease, levels = c("Healthy", "DS")),
    Treatment
  ) |> 
  dplyr::mutate(Individual = case_when(
    Individual == "5y_DS" ~ "DS 5y",
    Individual == "14y_DS" ~ "DS 14y",
    Individual == "NEO" ~ "Healthy",
    .default = Individual
  ))

ha1 = HeatmapAnnotation(Individual = metadata_correct_order$Individual,
                    Treatment = metadata_correct_order$Treatment,
                    col = list(
                      Individual = c(
                        "Healthy" = "#636363",
                        "DS 5y" = "#95abbd",
                        "DS 14y" = "#95dbea"
                      ),
                      Treatment = c(
                        "DMSO" = "gray",
                        "UMK57" = "cyan"
                      )
                    ))

group_split <- metadata_correct_order$Disease      

group_split <- factor(group_split, levels = c("Healthy", "DS"))

ssgsea_scores_df <- as.data.frame(ssgsea_scores)

ssgsea_scores_df_filt <- ssgsea_scores_df[,row.names(metadata_correct_order)]

ssgsea_scores_scaled <- t(scale(t(ssgsea_scores_df_filt)))

correct_hall_names <- gsub("HALLMARK_", "", rownames(ssgsea_scores_scaled))

row.names(ssgsea_scores_scaled) <- correct_hall_names

ssgsea_scores_scaled <- as.data.frame(ssgsea_scores_scaled)[,row.names(metadata_correct_order)]


set.seed(123)

ht_ds <- Heatmap(
  ssgsea_scores_scaled,
  cluster_columns = F,
  cluster_rows = T,
  col = col_fun,
  column_split = group_split,
  top_annotation = ha1,
  row_km = 4,
  row_km_repeats = 100,
  border_gp = gpar(col = "black", lty = 1),
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 11),
  column_names_rot = 90,
  show_column_names = F,
  heatmap_legend_param = list(
    at = c(-2, 0, 2),
    labels = c("-2", "0", "2"),
    title = "Row Z-Score",
    border = "black",
    title_position = "leftcenter-rot"))
ht_ds

# Heatmap HGPS


metadata_correct_order <- metadata |> 
  dplyr::filter(Disease %in% c("Healthy", "HGPS")) |> 
  dplyr::arrange(
    factor(Disease, levels = c("Healthy", "HGPS")),
    Treatment
  ) |> 
  dplyr::mutate(Individual = case_when(
    Individual == "HGPS 169" ~ "HGPS 8y",
    Individual == "NEO" ~ "Healthy",
    .default = Individual
  ))

ha1 = HeatmapAnnotation(Individual = metadata_correct_order$Individual,
                        Treatment = metadata_correct_order$Treatment,
                        col = list(
                          Individual = c(
                            "Healthy" = "#636363",
                            "HGPS 8y" = "#812723"
                          ),
                          Treatment = c(
                            "DMSO" = "gray",
                            "UMK57" = "cyan"
                          )
                        ))

group_split <- metadata_correct_order$Individual      

group_split <- factor(group_split, levels = c("Healthy", "HGPS 8y"))

ssgsea_scores_df <- as.data.frame(ssgsea_scores)

ssgsea_scores_df_filt <- ssgsea_scores_df[,row.names(metadata_correct_order)]

ssgsea_scores_scaled <- t(scale(t(ssgsea_scores_df_filt)))

correct_hall_names <- gsub("HALLMARK_", "", rownames(ssgsea_scores_scaled))

row.names(ssgsea_scores_scaled) <- correct_hall_names

ssgsea_scores_scaled <- as.data.frame(ssgsea_scores_scaled)[,row.names(metadata_correct_order)]

set.seed(123)

ht_hgps <- Heatmap(
  ssgsea_scores_scaled,
  cluster_columns = F,
  cluster_rows = T,
  col = col_fun,
  column_split = group_split,
  top_annotation = ha1,
  row_km = 4,
  row_km_repeats = 100,
  border_gp = gpar(col = "black", lty = 1),
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 11),
  column_names_rot = 90,
  show_column_names = F,
  heatmap_legend_param = list(
    at = c(-2, 0, 2),
    labels = c("-2", "0", "2"),
    title = "Row Z-Score",
    border = "black",
    title_position = "leftcenter-rot"))
ht_hgps

graphics.off()


pdf(
  file.path("output_data",
            "plots",
            "ssGSEA_and_reversion_magnitude",
            "heatmap_ssGSEA_Hall_DS.pdf")
)

draw(ht_ds)

dev.off()




png(
  file.path("output_data",
            "plots",
            "ssGSEA_and_reversion_magnitude",
            "heatmap_ssGSEA_Hall_DS.png"), res = 600, width = 7, height = 9, units = "in"
)

draw(ht_ds)

dev.off()






pdf(
  file.path("output_data",
            "plots",
            "ssGSEA_and_reversion_magnitude",
            "heatmap_ssGSEA_Hall_HGPS.pdf")
)

draw(ht_hgps)

dev.off()



png(
  file.path("output_data",
            "plots",
            "ssGSEA_and_reversion_magnitude",
            "heatmap_ssGSEA_Hall_HGPS.png"), res = 600, width = 7, height = 9, units = "in"
)

draw(ht_hgps)

dev.off()


# Write Session Info and Log Outputs --------------------------------------


cat(
  paste0("\n## Step completed\n",
         "Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
         "Finished 07_ssGSEA_and_reversion_magnitude\n"),
  file = "PROGRESS_LOG.md",
  append = TRUE
)


writeLines(capture.output(sessionInfo()), 
           file.path("session_info_logs","07_ssGSEA_and_reversion_magnitud_sessionInfo.txt"))


