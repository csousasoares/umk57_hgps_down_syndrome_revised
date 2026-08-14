## Load Packages and Set Seed --------------------------------------------------

set.seed(123)

library(DESeq2)
library(org.Hs.eg.db)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(ggalt)
library(limma)
library(circlize)
library(ComplexHeatmap)
library(ggtree)
library(patchwork)
library(ggrastr)

log_text <- c(
  "# Project Progress Log",
  paste("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Steps performed: ",
  "Started 1_umk57_vs_dmso_rna_seq_analysis"
)

writeLines(log_text, "PROGRESS_LOG.md")


dir.create(file.path("output_data", "results", "deseq2"), recursive = T)
dir.create(file.path("output_data", "results", "gsea"), recursive = T)
dir.create(file.path("output_data", "plots"), recursive = T)


# Load Counts and Metadata ----------------------------------------------------

# Load correct sample info and counts for DS, HGPS and Healthy
# samples obtained from IonTorrent Suite:

counts_matrix_correct <- read.csv(file.path("input_data", "counts_matrix.csv"),
                                  row.names = 1)

sample_info_correct <- read.csv(file.path("input_data", "sample_info.csv"),
                                row.names = 1)


# PCA HGPS UMK57 vs DMSO + Healthy --------------------------------------------

sample_info_hgps_umk57_healthy <- sample_info_correct %>% 
  dplyr::filter(Disease %in% c("HGPS", "Healthy"))

counts_matrix_healthy_hgps_umk57 <- counts_matrix_correct[,row.names(
  sample_info_hgps_umk57_healthy)]

colnames(counts_matrix_healthy_hgps_umk57) == row.names(sample_info_hgps_umk57_healthy)

# Should be TRUE.

# Create DESeq Object:

dds_healthy_hgps_umk57 <- DESeqDataSetFromMatrix(
  countData = counts_matrix_healthy_hgps_umk57,
  colData = sample_info_hgps_umk57_healthy,
  design = ~ Individual) # Makes no difference for PCA.

colnames(dds_healthy_hgps_umk57)

## Design wont matter at this stage though:

smallestGroupSize <- nrow(sample_info_hgps_umk57_healthy) / 2 # ~ Half of samples

keep <- rowSums(counts(dds_healthy_hgps_umk57) >=5) >= smallestGroupSize 

# Recomended in vignette.

dds_healthy_hgps_umk57 <- dds_healthy_hgps_umk57[keep,]

dds_healthy_hgps_umk57

dds_healthy_hgps_umk57 <- DESeq(dds_healthy_hgps_umk57)

vsd <- vst(dds_healthy_hgps_umk57, blind = TRUE) # Blind to ignore conditions.

pcaData <- plotPCA(vsd, 
                   intgroup = c("Treatment"), 
                   returnData = TRUE, 
                   ntop = 1000) # Using Top 1000 features.

percentVar <- round(100 * attr(pcaData, "percentVar"))


pcaData <- pcaData |> dplyr::mutate(
  Individual = case_when(
    Individual == "NEO" ~ "Healthy",
    Individual == "HGPS 169" ~ "HGPS 8y",
    .default = Individual
  )
)

shapes <- c(
  "HGPS 8y" = 15,
  "Healthy" = 16,
  "DS 14y" = 17,
  "DS 5y" = 18
)

cols_pca <- c(
  "DMSO" = "gray30",
  "UMK57" = "cyan"
)


pca_1 <- ggplot(pcaData, aes(PC1, PC2, color=Treatment, 
                             label = Individual)) +
  geom_point(size=4, alpha = 0.5, aes(shape = Individual)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme_bw(base_size = 15) +
  xlim(-25, 60) +
  ylim(-15, 20) +
  scale_shape_manual(values = shapes) +
  scale_color_manual(values = cols_pca) +
  theme(panel.grid = element_blank()) +
  guides(fill = guide_legend(override.aes = list(linetype = 0)),
         color = guide_legend(override.aes = list(linetype = 0)),
         shape = guide_legend(override.aes = list(linetype = 0))) +
  labs(title = "Principal Component Analysis - HGPS")

pca_1


# PCA DS UMK57 vs DMSO + Healthy ----------------------------------------------

sample_info_ds_umk57_healthy <- sample_info_correct %>% 
  dplyr::filter(Disease %in% c("DS", "Healthy"))

counts_matrix_healthy_ds_umk57 <- counts_matrix_correct[,row.names(
  sample_info_ds_umk57_healthy)]

colnames(counts_matrix_healthy_ds_umk57) == row.names(sample_info_ds_umk57_healthy)

# Should be TRUE.

dds_healthy_ds_umk57 <- DESeqDataSetFromMatrix(
  countData = counts_matrix_healthy_ds_umk57,
  colData = sample_info_ds_umk57_healthy,
  design = ~ Individual) # Makes no difference for PCA.

colnames(dds_healthy_ds_umk57)

# Design wont matter at this stage though:

smallestGroupSize <- nrow(sample_info_ds_umk57_healthy) / 2 # Half of samples.

keep <- rowSums(counts(dds_healthy_ds_umk57) >=5) >= smallestGroupSize 

#Recomended in vignette.

dds_healthy_ds_umk57 <- dds_healthy_ds_umk57[keep,]

dds_healthy_ds_umk57

dds_healthy_ds_umk57 <- DESeq(dds_healthy_ds_umk57)


vsd <- vst(dds_healthy_ds_umk57, blind=TRUE) ## Blind to ignore conditions
pcaData <- plotPCA(vsd, intgroup=c("Treatment"), 
                   returnData=TRUE, ntop = 1000)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaData <- pcaData |> dplyr::mutate(
  Individual = case_when(
    Individual == "NEO" ~ "Healthy",
    Individual == "14y_DS" ~ "DS 14y",
    Individual == "5y_DS" ~ "DS 5y",
    .default = Individual
  )
)


pca_2 <- ggplot(pcaData, aes(PC1, PC2, color=Treatment, 
                             label = Individual)) +
  geom_point(size=4, alpha = 0.5, aes(shape = Individual)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw(base_size = 15) +
  xlim(-70, 50) +
  ylim(-35, 40) +
  scale_shape_manual(values = shapes) +
  scale_color_manual(values = cols_pca) +
  theme(panel.grid = element_blank()) +
  guides(fill = guide_legend(override.aes = list(linetype = 0)),
         color = guide_legend(override.aes = list(linetype = 0)),
         shape = guide_legend(override.aes = list(linetype = 0))) +
  labs(title = "Principal Component Analysis - DS")


pca_1 / pca_2

dir.create(file.path("output_data", "plots", "pca"),
           recursive = T)

ggsave(
  "pca_all_healthy_hgps_ds_umk57.pdf",
  path = file.path("output_data", "plots", "pca"),
  width = 7.5,
  height = 6.5, 
  create.dir = T
)


# DESeq2 HGPS UMK57 vs DMSO ---------------------------------------------------

hgps_of_interest <- sample_info_correct %>% 
  dplyr::filter(Disease == "HGPS") %>% 
  dplyr::pull(X)

length(hgps_of_interest) # 8 Samples

sample_info_hgps_of_interest <- sample_info_correct[hgps_of_interest,]

counts_matrix_correct_hgps_of_interest <- counts_matrix_correct[,hgps_of_interest]

row.names(sample_info_hgps_of_interest) == colnames(
  counts_matrix_correct_hgps_of_interest)

dds_hgps_of_interest <- DESeqDataSetFromMatrix(
  countData = counts_matrix_correct_hgps_of_interest,
  colData = sample_info_hgps_of_interest,
  design = ~ Treatment) 

smallestGroupSize <- nrow(sample_info_hgps_of_interest) / 2 # Half of samples.

keep <- rowSums(counts(dds_hgps_of_interest) >=5) >= smallestGroupSize 

# Filter very low expression genes, recommended in vignette:

dds_hgps_of_interest <- dds_hgps_of_interest[keep,]

dds_hgps_of_interest

dds_hgps_of_interest <- DESeq(dds_hgps_of_interest)

res_hgps_of_interest <- results(dds_hgps_of_interest, 
                                contrast = c("Treatment", "UMK57", "DMSO"), 
                                cooksCutoff = FALSE, 
                                independentFiltering = FALSE) 


df_hgps_of_interest <- as.data.frame(res_hgps_of_interest)


### HGPS UMK57 vs DMSO Volcano Plot -----------------------------------------


dir.create(file.path("output_data", "plots", "volcano_plots"),
           recursive = T)


df_hgps_of_interest$log10padj <- -log10(df_hgps_of_interest$padj)
df_hgps_of_interest$gene <- row.names(df_hgps_of_interest)


write.csv(
  df_hgps_of_interest, 
  file.path(
    "output_data",
    "results",
    "deseq2",
    "deseq2_hgps_umk57_vs_dmso.csv"
  )
)


df_hgps_of_interest <- df_hgps_of_interest %>% 
  dplyr::mutate(direction = case_when(
    log2FoldChange > 0 & padj < 0.05 ~ "Up",
    log2FoldChange < 0 & padj < 0.05 ~ "Down",
    .default = "No Change"
  )) %>% 
  dplyr::mutate(gene_label = case_when(
    direction %in% c("Up", "Down") ~ gene,
    .default = ""
  ))

cols_volcano <- c(
  "Up" = "red",
  "Down" = "blue",
  "No Change" = "gray"
)

n_up <- df_hgps_of_interest %>% 
  dplyr::filter(direction == "Up") %>% 
  nrow()

n_down <- df_hgps_of_interest %>% 
  dplyr::filter(direction == "Down") %>% 
  nrow()

n_unchanged <- df_hgps_of_interest %>% 
  dplyr::filter(direction == "No Change") %>% 
  nrow()


ggplot(df_hgps_of_interest, aes(
  x = log2FoldChange, y = log10padj, label = gene_label
)) +
  ggrastr::rasterise(geom_point(aes(color = direction), alpha = 0.5,
                                size = 0.75), dpi = 300) +
  theme_bw(base_size = 12) +
  xlim(-6.5, 5) +
  coord_fixed(0.25) +
  scale_color_manual(name = "Direction", values = cols_volcano, labels = c(
    paste0("Down: ", n_down), paste0("No Change: ", n_unchanged), 
    paste0("Up: ", n_up))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  labs(y = "-log10(padjust)\n", x = "Log2FC",
       title = "Volcano Plot - HGPS UMK57 vs DMSO")

ggsave("volcano_plot_hgps_umk57_dmso.pdf",
       path = file.path("output_data", "plots", "volcano_plots"),
       width = 7,
       height = 7)


# DESeq2 DS UMK57 vs DMSO  ----------------------------------------------------

sample_info_ds_of_interest <- sample_info_correct %>% 
  dplyr::filter(Individual %in% c("14y_DS", "5y_DS"))

counts_matrix_correct_ds_of_interest <- counts_matrix_correct[,row.names(
  sample_info_ds_of_interest)]

row.names(sample_info_ds_of_interest) == colnames(
  counts_matrix_correct_ds_of_interest)

dds_ds_of_interest <- DESeqDataSetFromMatrix(
  countData = counts_matrix_correct_ds_of_interest,
  colData = sample_info_ds_of_interest,
  design = ~ Treatment) 

smallestGroupSize <- nrow(sample_info_ds_of_interest) / 2

keep <- rowSums(counts(dds_ds_of_interest) >=5) >= smallestGroupSize 

# Recomended in vignette.

dds_ds_of_interest <- dds_ds_of_interest[keep,]

dds_ds_of_interest

dds_ds_of_interest <- DESeq(dds_ds_of_interest)
res_ds_of_interest <- results(dds_ds_of_interest, 
                              contrast = c("Treatment", "UMK57", "DMSO"), 
                              cooksCutoff = FALSE, 
                              independentFiltering = FALSE)


df_ds_of_interest <- as.data.frame(res_ds_of_interest)


### DS UMK57 vs DMSO Volcano Plot ----------------------------------------------

df_ds_of_interest$log10padj <- -log10(df_ds_of_interest$padj)
df_ds_of_interest$gene <- row.names(df_ds_of_interest)

write.csv(
  df_ds_of_interest, 
  file.path("output_data", "results", "deseq2", "deseq2_ds_umk57_vs_dmso.csv")
)


df_ds_of_interest <- df_ds_of_interest %>% 
  dplyr::mutate(direction = case_when(
    log2FoldChange > 0 & padj < 0.05 ~ "Up",
    log2FoldChange < 0 & padj < 0.05 ~ "Down",
    .default = "No Change"
  )) %>% 
  dplyr::mutate(gene_label = case_when(
    direction %in% c("Up", "Down") ~ gene,
    .default = ""
  ))

cols_volcano <- c(
  "Up" = "red",
  "Down" = "blue",
  "No Change" = "gray"
)

n_up <- df_ds_of_interest %>% 
  dplyr::filter(direction == "Up") %>% 
  nrow()

n_down <- df_ds_of_interest %>% 
  dplyr::filter(direction == "Down") %>% 
  nrow()

n_unchanged <- df_ds_of_interest %>% 
  dplyr::filter(direction == "No Change") %>% 
  nrow()


ggplot(df_ds_of_interest, aes(
  x = log2FoldChange, y = log10padj, label = gene_label
)) +
  ggrastr::rasterise(geom_point(aes(color = direction), alpha = 0.5,
                                size = 0.75), dpi = 300) +
  theme_bw(base_size = 12) +
  coord_fixed(0.25) +
  xlim(-6,12) +
  scale_color_manual(name = "Direction", values = cols_volcano, labels = c(
    paste0("Down: ", n_down), paste0("No Change: ", n_unchanged), paste0("Up: ", n_up))) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  labs(y = "-log10(padjust)\n", x = "Log2FC",
       title = "Volcano Plot - DS UMK57 vs DMSO")

ggsave("volcano_plot_ds_umk57_dmso.pdf",
       path = file.path("output_data", "plots", "volcano_plots"),
       width = 7,
       height = 7)




# GSEA HGPS UMK57 vs DMSO Hallmark --------------------------------------------

dir.create(file.path("output_data", "plots", "gsea"),
           recursive = T)


head(df_hgps_of_interest)

df_hgps_of_interest <- df_hgps_of_interest %>% 
  dplyr::arrange(desc(log2FoldChange))
head(df_hgps_of_interest)

hgps_umk57_ranked <- df_hgps_of_interest$log2FoldChange
names(hgps_umk57_ranked) <- df_hgps_of_interest$gene
head(hgps_umk57_ranked, 10)


saveRDS(
  hgps_umk57_ranked,
  file.path("output_data", "results", "gsea", "hgps_umk57_vs_dmso_ranked_list_for_GSEA.rds")
)

hall <- msigdbr::msigdbr(
  species = "Homo sapiens",
  collection = "H"
)

hall_t2g <- hall %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()

hgps_umk57_hall <- clusterProfiler::GSEA(
  hgps_umk57_ranked,
  exponent = 1,
  minGSSize = 0,
  maxGSSize = 1000,
  eps = 1e-50,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = hall_t2g,
  verbose = TRUE,
  seed = TRUE,
  by = "fgsea",
)

df_hgps_umk57_hall <- as.data.frame(hgps_umk57_hall)

write.csv(
  df_hgps_umk57_hall, 
  file.path("output_data", "results", "gsea", "hgps_umk57_dmso_GSEA_hallmarks.csv")
  )

df_hgps_umk57_hall_005 <- df_hgps_umk57_hall %>% 
  dplyr::filter(p.adjust < 0.05)

df_hgps_umk57_hall_005$Description <- gsub(
  "HALLMARK_", "", df_hgps_umk57_hall_005$Description
)

df_hgps_umk57_hall_005 <- df_hgps_umk57_hall_005 %>% 
  dplyr::mutate(direction = case_when(
    NES > 0 ~ "Up",
    NES < 0 ~ "Down"
  ))

cols <- c(
  "Up" = "red",
  "Down" = "blue"
)

ggplot(df_hgps_umk57_hall_005,
       aes(x = NES, y = forcats::fct_reorder(Description, NES))) +
  geom_segment(aes(xend = 0, yend = Description)) +
  geom_point(aes(fill = direction, size = p.adjust), shape = 21) +
  scale_fill_manual(values = cols, name = "Direction") +
  scale_size_continuous(range = c(8,3), transform = "log10", 
                        name = "padjust") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlim(-3,3) +
  labs(x = "NES", y = "Gene Sets\n",
       title = "HGPS UMK57 vs DMSO",
       subtitle = "Hallmark Gene Sets") +
  geom_vline(xintercept = 0, linetype = "dashed")


ggsave(
  "hallmarks_hgps_umk57_dmso_dotplot.pdf",
  path = file.path("output_data", "plots", "gsea"),
  width = 8,
  height = 8
)


# GSEA DS UMK57 vs DMSO Hallmark ----------------------------------------------

df_ds_of_interest$gene <- row.names(df_ds_of_interest)

df_ds_of_interest <- df_ds_of_interest %>% 
  dplyr::arrange(desc(log2FoldChange))

ds_of_interest_ranked <- df_ds_of_interest$log2FoldChange
names(ds_of_interest_ranked) <- df_ds_of_interest$gene
head(ds_of_interest_ranked)


saveRDS(
  ds_of_interest_ranked,
  file.path("output_data", "results", "gsea", "ds_umk57_vs_dmso_ranked_list_for_GSEA.rds")
)


Hall <- msigdbr::msigdbr(
  species = "Homo sapiens",
  collection = "H"
)

hall_t2g <- Hall %>% dplyr::distinct(gs_name, gene_symbol) %>% 
  as.data.frame()

Hall_ds_of_interest <- clusterProfiler::GSEA(
  ds_of_interest_ranked,
  exponent = 1,
  minGSSize = 0,
  maxGSSize = 1000,
  eps = 1e-50,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = hall_t2g,
  verbose = TRUE,
  seed = TRUE,
  by = "fgsea",
)

df_Hall_ds_of_interest <- as.data.frame(Hall_ds_of_interest)
head(df_Hall_ds_of_interest)

write.csv(
  df_Hall_ds_of_interest, 
  file.path("output_data", "results", "gsea", "ds_umk57_dmso_GSEA_hallmarks.csv")
  )


df_Hall_ds_of_interest_005 <- df_Hall_ds_of_interest %>% 
  dplyr::filter(p.adjust < 0.05)

df_Hall_ds_of_interest_005$Description <- gsub(
  "HALLMARK_", "", df_Hall_ds_of_interest_005$Description
)

df_Hall_ds_of_interest_005 <- df_Hall_ds_of_interest_005 %>% 
  dplyr::mutate(direction = case_when(
    NES > 0 ~ "Up",
    NES < 0 ~ "Down"
  ))

cols <- c(
  "Up" = "red",
  "Down" = "blue"
)

ggplot(df_Hall_ds_of_interest_005,
       aes(x = NES, y = forcats::fct_reorder(Description, NES))) +
  geom_segment(aes(xend = 0, yend = Description)) +
  geom_point(aes(fill = direction, size = p.adjust), shape = 21) +
  scale_fill_manual(values = cols, name = "Direction") +
  scale_size_continuous(range = c(8,3), transform = "log10", name = "padjust") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlim(-3,3) +
  labs(x = "NES", y = "Gene Sets\n",
       title = "DS UMK57 vs DMSO",
       subtitle = "Hallmark Gene Sets") +
  geom_vline(xintercept = 0, linetype = "dashed")


ggsave(
  "hallmarks_ds_umk57_dmso_dotplot.pdf",
  path = file.path("output_data", "plots", "gsea"),
  width = 8,
  height = 8
)


## GSEA HGPS UMK57 vs DMSO DNA Repair ------------------------------------------

mm_BP_sets <- msigdbr::msigdbr(
  species = "Homo sapiens")
mm_BP_sets

msigdbr_t2g = mm_BP_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% 
  as.data.frame()
msigdbr_t2g

msigdbr_t2g_dna_rep_terms <- msigdbr_t2g %>% 
  dplyr::filter(str_detect(gs_name, "REPAIR")) %>% 
  dplyr::group_by(gs_name) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::filter(n > 30) ## Get "Repair" terms across MSigDB with at least 30 genes

msigdbr_t2g_dna_rep_terms ## Check which ones pop up

## From 38 "Repair" sets, manually curated terms actually related to 
## global DNA repair (e.g. GOBP or Reactome DNA Repair)
## and specific DNA repair processes, such as DSB Repair,
## NER, BER, etc from well known collections:

misgdbr_t2g_filtered <- filter(msigdbr_t2g, gs_name %in% c("GOBP_BASE_EXCISION_REPAIR",
                                                           "GOBP_DNA_REPAIR",
                                                           "GOBP_DNA_SYNTHESIS_INVOLVED_IN_DNA_REPAIR",
                                                           "GOBP_DOUBLE_STRAND_BREAK_REPAIR",
                                                           "GOBP_DOUBLE_STRAND_BREAK_REPAIR_VIA_NONHOMOLOGOUS_END_JOINING",
                                                           "GOBP_INTERSTRAND_CROSS_LINK_REPAIR",
                                                           "GOBP_MISMATCH_REPAIR",
                                                           "GOBP_NUCLEOTIDE_EXCISION_REPAIR",
                                                           "GOBP_POSITIVE_REGULATION_OF_DNA_REPAIR",
                                                           "GOBP_POSITIVE_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR",
                                                           "GOBP_POSITIVE_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR_VIA_HOMOLOGOUS_RECOMBINATION",
                                                           "GOBP_POSTREPLICATION_REPAIR",
                                                           "GOBP_RECOMBINATIONAL_REPAIR",
                                                           "GOBP_REGULATION_OF_DNA_REPAIR",
                                                           "GOBP_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR",
                                                           "GOBP_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR_VIA_HOMOLOGOUS_RECOMBINATION",
                                                           "HALLMARK_DNA_REPAIR",
                                                           "KAUFFMANN_DNA_REPAIR_GENES",
                                                           "KEGG_BASE_EXCISION_REPAIR",
                                                           "KEGG_NUCLEOTIDE_EXCISION_REPAIR",
                                                           "REACTOME_BASE_EXCISION_REPAIR",
                                                           "REACTOME_BASE_EXCISION_REPAIR_AP_SITE_FORMATION",
                                                           "REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR",
                                                           "REACTOME_DNA_REPAIR",
                                                           "REACTOME_GLOBAL_GENOME_NUCLEOTIDE_EXCISION_REPAIR_GG_NER",
                                                           "REACTOME_HOMOLOGY_DIRECTED_REPAIR",
                                                           "REACTOME_NUCLEOTIDE_EXCISION_REPAIR",
                                                           "REACTOME_SUMOYLATION_OF_DNA_DAMAGE_RESPONSE_AND_REPAIR_PROTEINS",
                                                           "REACTOME_TRANSCRIPTION_COUPLED_NUCLEOTIDE_EXCISION_REPAIR_TC_NER",
                                                           "WP_BASE_EXCISION_REPAIR",
                                                           "WP_DNA_REPAIR_PATHWAYS_FULL_NETWORK",
                                                           "WP_NUCLEOTIDE_EXCISION_REPAIR"))
misgdbr_t2g_filtered


hgps_dna_repair <- clusterProfiler::GSEA(
  hgps_umk57_ranked,
  exponent = 1,
  minGSSize = 0,
  maxGSSize = 1000,
  eps = 1e-50,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = misgdbr_t2g_filtered,
  verbose = T,
  seed = T,
  by = "fgsea",
)

saveRDS(hgps_dna_repair, 
        file.path("output_data", "results", "gsea","hgps_dna_repair_gsea.rds")
        )

df_hgps_dna_repair <- as.data.frame(hgps_dna_repair)

write.csv(
  df_hgps_dna_repair,
  file.path("output_data", "results", "gsea", "df_hgps_umk57_vs_dmso_GSEA_DNA_repair.csv")
)

df_hgps_dna_repair_005 <- df_hgps_dna_repair %>% 
  dplyr::filter(p.adjust < 0.05)

df_hgps_dna_repair_005$Description <- gsub("_", " ", 
                                           df_hgps_dna_repair_005$Description)

df_hgps_dna_repair_005$Description <- stringr::str_wrap(df_hgps_dna_repair_005$Description,
                                                        width = 30)

df_hgps_dna_repair_005$Description

df_hgps_dna_repair_005 <- df_hgps_dna_repair_005 %>% 
  dplyr::mutate(direction = case_when(
    NES > 0 ~ "Up",
    NES < 0 ~ "Down"
  ))

ggplot(df_hgps_dna_repair_005,
       aes(x = NES, y = forcats::fct_reorder(Description, NES))) +
  geom_segment(aes(xend = 0, yend = Description)) +
  geom_point(aes(fill = direction, size = p.adjust), shape = 21) +
  scale_fill_manual(values = cols, name = "Direction") +
  scale_size_continuous(range = c(8,3), transform = "log10", name = "padjust") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlim(-3,3) +
  labs(x = "NES", y = "Gene Sets\n",
       title = "HGPS UMK57 vs DMSO",
       subtitle = "DNA Repair Gene Sets") +
  geom_vline(xintercept = 0, linetype = "dashed")


## GSEA DS UMK57 vs DMSO DNA Repair --------------------------------------------


t21_dna_repair <- clusterProfiler::GSEA(
  ds_of_interest_ranked,
  exponent = 1,
  minGSSize = 0,
  maxGSSize = 1000,
  eps = 1e-50,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = misgdbr_t2g_filtered,
  verbose = T,
  seed = T,
  by = "fgsea",
)

saveRDS(t21_dna_repair, 
        file.path("output_data", "results", "gsea", "ds_dna_repair_gsea.rds")
)


df_t21_dna_repair <- as.data.frame(t21_dna_repair)

write.csv(
  df_t21_dna_repair,
  file.path("output_data", "results", "gsea", "df_ds_umk57_vs_dmso_GSEA_DNA_repair.csv")
)

df_t21_dna_repair_005 <- df_t21_dna_repair %>% 
  dplyr::filter(p.adjust < 0.05)

df_t21_dna_repair_005$Description <- gsub("_", " ", 
                                          df_t21_dna_repair_005$Description)

df_t21_dna_repair_005$Description <- stringr::str_wrap(df_t21_dna_repair_005$Description,
                                                       width = 30)

df_t21_dna_repair_005$Description

df_t21_dna_repair_005 <- df_t21_dna_repair_005 %>% 
  dplyr::mutate(direction = case_when(
    NES > 0 ~ "Up",
    NES < 0 ~ "Down"
  ))

ggplot(df_t21_dna_repair_005,
       aes(x = NES, y = forcats::fct_reorder(Description, NES))) +
  geom_segment(aes(xend = 0, yend = Description)) +
  geom_point(aes(fill = direction, size = p.adjust), shape = 21) +
  scale_fill_manual(values = cols, name = "Direction") +
  scale_size_continuous(range = c(8,3), transform = "log10", name = "padjust") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlim(-3,3) +
  labs(x = "NES", y = "Gene Sets\n",
       title = "DS UMK57 vs DMSO",
       subtitle = "DNA Repair Gene Sets") +
  geom_vline(xintercept = 0, linetype = "dashed")


## GSEA DNA Repair Plot --------------------------------------------------------

df_hgps_dna_repair_005[1:5,1:5]
df_t21_dna_repair_005[1:5,1:5]


df_hgps_dna_repair_005$comparison <- "HGPS"
df_t21_dna_repair_005$comparison <- "DS"

dna_repair_convergent_plot <- rbind(
  df_hgps_dna_repair_005,
  df_t21_dna_repair_005
)

## DNA Repair DotPlot

new_dna_repair_terms <- data.frame(
  ID = c("GOBP_BASE_EXCISION_REPAIR",
         "GOBP_DNA_REPAIR",
         "GOBP_DNA_SYNTHESIS_INVOLVED_IN_DNA_REPAIR",
         "GOBP_DOUBLE_STRAND_BREAK_REPAIR",
         "GOBP_DOUBLE_STRAND_BREAK_REPAIR_VIA_NONHOMOLOGOUS_END_JOINING",
         "GOBP_INTERSTRAND_CROSS_LINK_REPAIR",
         "GOBP_MISMATCH_REPAIR",
         "GOBP_NUCLEOTIDE_EXCISION_REPAIR",
         "GOBP_POSITIVE_REGULATION_OF_DNA_REPAIR",
         "GOBP_POSITIVE_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR",
         "GOBP_POSITIVE_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR_VIA_HOMOLOGOUS_RECOMBINATION",
         "GOBP_POSTREPLICATION_REPAIR",
         "GOBP_RECOMBINATIONAL_REPAIR",
         "GOBP_REGULATION_OF_DNA_REPAIR",
         "GOBP_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR",
         "GOBP_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR_VIA_HOMOLOGOUS_RECOMBINATION",
         "HALLMARK_DNA_REPAIR",
         "KAUFFMANN_DNA_REPAIR_GENES",
         "KEGG_BASE_EXCISION_REPAIR",
         "KEGG_NUCLEOTIDE_EXCISION_REPAIR",
         "REACTOME_BASE_EXCISION_REPAIR",
         "REACTOME_BASE_EXCISION_REPAIR_AP_SITE_FORMATION",
         "REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR",
         "REACTOME_DNA_REPAIR",
         "REACTOME_GLOBAL_GENOME_NUCLEOTIDE_EXCISION_REPAIR_GG_NER",
         "REACTOME_HOMOLOGY_DIRECTED_REPAIR",
         "REACTOME_NUCLEOTIDE_EXCISION_REPAIR",
         "REACTOME_SUMOYLATION_OF_DNA_DAMAGE_RESPONSE_AND_REPAIR_PROTEINS",
         "REACTOME_TRANSCRIPTION_COUPLED_NUCLEOTIDE_EXCISION_REPAIR_TC_NER",
         "WP_BASE_EXCISION_REPAIR",
         "WP_DNA_REPAIR_PATHWAYS_FULL_NETWORK",
         "WP_NUCLEOTIDE_EXCISION_REPAIR"),
  "new" = c(1:32))


new_dna_repair_terms$new <- c(
  "(BP) Base Excision Repair",
  "(BP) DNA Repair",
  "(BP) DNA Synthesis Involved in DNA Repair",
  "(BP) DSB Repair",
  "(BP) DSB Repair via NHEJ",
  "(BP) Interstrand Crosslink Repair",
  "(BP) Mismatch Repair",
  "(BP) Nucleotide Excision Repair",
  "(BP) Positive Regulation of DNA Repair",
  "(BP)  Positive Regulation of DSB Repair",
  "(BP) Positive Regulation of DSB Repair via HR",
  "(BP) Postreplication Repair",
  "(BP) Recombinational Repair",
  "(BP) Regulation of DNA Repair",
  "(BP) Regulation of DSB Repair",
  "(BP) Regulation of DSB Repair via HR",
  "(H) DNA Repair",
  "Kauffman DNA Repair Genes",
  "(KEGG) Base Excision Repair",
  "(KEGG) Nucleotide Excision Repair (NER)",
  "(R) Base Excision Repair (BER)",
  "(R) BER - AP Site Formation",
  "(R) DNA DSB Repair",
  "(R) DNA Repair",
  "(R) GG-NER",
  "(R) Homology-Directed Repair",
  "(R) Nucleotide Excision Repair (NER)",
  "(R) Sumoylation of DNA Repair Proteins",
  "(R) TC-NER",
  "(WP) Base Excision Repair (BER)",
  "(WP) DNA Repair Pathways Full Network",
  "(WP) Nucleotide Excision Repair (NER)"
)


dna_repair_convergent_plot$comparison <- factor(
  dna_repair_convergent_plot$comparison,
  levels = c("HGPS",
             "DS")
)

dna_repair_convergent_plot <- dna_repair_convergent_plot %>% 
  left_join(new_dna_repair_terms, join_by(ID))

dna_repair_convergent_plot$comparison <- factor(
  dna_repair_convergent_plot$comparison,
  levels = c("HGPS", "DS")
)

dna_rep_dot <- ggplot(dna_repair_convergent_plot, aes(
  x = comparison, y = forcats::fct_reorder(new, setSize))) +
  geom_point(aes(fill = NES, size = p.adjust),
             shape = 21) +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0,
                       limits = c(-1, 3)) +
  scale_size_continuous(range = c(8, 2),
                        transform = "log10",
                        name = "padjust") +
  theme_bw(base_size = 14) +
  labs(x = "", y = "",
       title = "DNA Repair - UMK57 vs DMSO") +
  guides(fill = guide_colorbar(
    frame.colour = "black",
    frame.linewidth = 0.25,
    ticks.linewidth = 0.25,
    ticks.colour = "black"))
dna_rep_dot

ggsave(
  "dna_repair_GSEA_hgps_ds_umk57_dmso.pdf",
  path = file.path("output_data", "plots", "gsea"),
  width = 6,
  height = 6
)

## GSEA Different Metric -------------------------------------------------------

head(df_hgps_of_interest)

df_hgps_of_interest <- df_hgps_of_interest %>% 
  dplyr::mutate(sign = log2FoldChange/abs(log2FoldChange)) %>% 
  dplyr::mutate(metric = sign * (-log10(pvalue))) %>% 
  dplyr::arrange(desc(metric))

head(df_hgps_of_interest)

hgps_umk57_ranked <- df_hgps_of_interest$metric
names(hgps_umk57_ranked) <- df_hgps_of_interest$gene
head(hgps_umk57_ranked, 10)


hgps_umk57_hall_sign_pval <- clusterProfiler::GSEA(
  hgps_umk57_ranked,
  exponent = 1,
  minGSSize = 0,
  maxGSSize = 1000,
  eps = 1e-50,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = hall_t2g,
  verbose = TRUE,
  seed = TRUE,
  by = "fgsea",
)

df_hgps_umk57_hall_sign_pval <- as.data.frame(hgps_umk57_hall_sign_pval)

df_hgps_umk57_hall_sign_pval_005 <- df_hgps_umk57_hall_sign_pval %>% 
  dplyr::filter(p.adjust < 0.05)

df_hgps_umk57_hall_sign_pval_005$Description <- gsub(
  "HALLMARK_", "", df_hgps_umk57_hall_sign_pval_005$Description
)

df_hgps_umk57_hall_sign_pval_005 <- df_hgps_umk57_hall_sign_pval_005 %>% 
  dplyr::mutate(direction = case_when(
    NES > 0 ~ "Up",
    NES < 0 ~ "Down"
  ))

cols <- c(
  "Up" = "red",
  "Down" = "blue"
)

alex_hgps <- ggplot(df_hgps_umk57_hall_sign_pval_005,
                    aes(x = NES, y = forcats::fct_reorder(Description, NES))) +
  geom_segment(aes(xend = 0, yend = Description)) +
  geom_point(aes(fill = direction, size = p.adjust), shape = 21) +
  scale_fill_manual(values = cols, name = "Direction") +
  scale_size_continuous(range = c(8,3), transform = "log10", 
                        name = "padjust") +
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlim(-3,3) +
  labs(x = "NES", y = "Gene Sets\n",
       title = "HGPS UMK57 vs DMSO",
       subtitle = "Hallmark Gene Sets (sign * -log10(P-value))") +
  geom_vline(xintercept = 0, linetype = "dashed")
alex_hgps

## For DS UMK57 vs DMSO:

df_ds_of_interest <- df_ds_of_interest %>% 
  dplyr::mutate(sign = log2FoldChange/abs(log2FoldChange)) %>% 
  dplyr::mutate(metric = sign * (-log10(pvalue))) %>% 
  dplyr::arrange(desc(metric))

ds_of_interest_ranked <- df_ds_of_interest$metric
names(ds_of_interest_ranked) <- df_ds_of_interest$gene
head(ds_of_interest_ranked)

Hall_ds_of_interest_sign_pval <- clusterProfiler::GSEA(
  ds_of_interest_ranked,
  exponent = 1,
  minGSSize = 0,
  maxGSSize = 1000,
  eps = 1e-50,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = hall_t2g,
  verbose = TRUE,
  seed = TRUE,
  by = "fgsea",
)

df_Hall_ds_of_interest_sign_pval <- as.data.frame(
  Hall_ds_of_interest_sign_pval)

head(df_Hall_ds_of_interest_sign_pval)


df_Hall_ds_of_interest_sign_pval_005 <- df_Hall_ds_of_interest_sign_pval %>% 
  dplyr::filter(p.adjust < 0.05)

df_Hall_ds_of_interest_sign_pval_005$Description <- gsub(
  "HALLMARK_", "", df_Hall_ds_of_interest_sign_pval_005$Description
)

df_Hall_ds_of_interest_sign_pval_005 <- df_Hall_ds_of_interest_sign_pval_005 %>% 
  dplyr::mutate(direction = case_when(
    NES > 0 ~ "Up",
    NES < 0 ~ "Down"
  ))

cols <- c(
  "Up" = "red",
  "Down" = "blue"
)

alex_ds <- ggplot(df_Hall_ds_of_interest_sign_pval_005,
                  aes(x = NES, y = forcats::fct_reorder(Description, NES))) +
  geom_segment(aes(xend = 0, yend = Description)) +
  geom_point(aes(fill = direction, size = p.adjust), shape = 21) +
  scale_fill_manual(values = cols, name = "Direction") +
  scale_size_continuous(range = c(8,3), transform = "log10", name = "padjust") +
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlim(-3,3) +
  labs(x = "NES", y = "Gene Sets\n",
       title = "DS UMK57 vs DMSO",
       subtitle = "Hallmark Gene Sets (sign * -log10(P-value))") +
  geom_vline(xintercept = 0, linetype = "dashed")
alex_ds

alex_hgps + alex_ds

ggsave(
  "sign_pval_hallmarks_hgps_ds_umk57_dmso_padj_005_FOR_REBUTTAL.png",
  path = file.path("output_data", "plots", "gsea"),
  width = 14.5,
  height = 6.5
)

# Write Session Info and Log Outputs --------------------------------------

cat(
  paste0("\n## Step completed\n",
         "Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
         "Finished 1_umk57_vs_dmso_rna_seq_analysis\n"),
  file = "PROGRESS_LOG.md",
  append = TRUE
)

dir.create(
  file.path("session_info_logs")
)

writeLines(capture.output(sessionInfo()), 
           file.path("session_info_logs", 
                     "01_umk57_vs_dmso_rna_seq_analysis_sessionInfo.txt")
)
