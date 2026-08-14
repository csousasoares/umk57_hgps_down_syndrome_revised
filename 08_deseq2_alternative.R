## Load Packages and Set Seed --------------------------------------------------

set.seed(123)

library(DESeq2)
library(org.Hs.eg.db)
library(tidyverse)
library(patchwork)

cat(
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
         "Started 08_deseq2_alternative\n"),
  file = "PROGRESS_LOG.md",
  append = TRUE
)


dir.create(file.path("output_data", "results", "deseq2_alternative"), recursive = T)
dir.create(file.path("output_data", "results", "gsea_alternative"), recursive = T)
dir.create(file.path("output_data", "plots", "gsea_alternative"), recursive = T)


# Load Counts and Metadata ----------------------------------------------------

# Load correct sample info and counts for DS, HGPS and Healthy
# samples obtained from IonTorrent Suite:

counts_matrix_correct <- read.csv(file.path("input_data", "counts_matrix.csv"),
                                  row.names = 1)

sample_info_correct <- read.csv(file.path("input_data", "sample_info.csv"),
                                row.names = 1)



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
  design = ~ Exp_Batch_Correct + Treatment) 

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

df_hgps_of_interest$gene <- row.names(df_hgps_of_interest)

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
  design = ~ Individual + Treatment) 

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


df_ds_of_interest$gene <- row.names(df_ds_of_interest)


# GSEA HGPS UMK57 vs DMSO Hallmark ----------------------------------------

head(df_hgps_of_interest)

df_hgps_of_interest <- df_hgps_of_interest %>% 
  dplyr::arrange(desc(log2FoldChange))
head(df_hgps_of_interest)

hgps_umk57_ranked <- df_hgps_of_interest$log2FoldChange
names(hgps_umk57_ranked) <- df_hgps_of_interest$gene
head(hgps_umk57_ranked, 10)


saveRDS(
  hgps_umk57_ranked,
  file.path("output_data", "results", "gsea_alternative", "hgps_umk57_vs_dmso_ranked_list_for_GSEA.rds")
)

hall <- msigdbr::msigdbr(
  species = "Homo sapiens",
  collection = "H"
)

hall_t2g <- hall |> dplyr::distinct(gs_name, gene_symbol) |> as.data.frame()

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
  file.path("output_data", "results", "gsea_alternative", "hgps_umk57_dmso_GSEA_hallmarks.csv")
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

hall_hgps_plot <- ggplot(df_hgps_umk57_hall_005,
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
       title = "HGPS UMK57 vs DMSO (Batch Adjusted)",
       subtitle = "Hallmark Gene Sets") +
  geom_vline(xintercept = 0, linetype = "dashed")
hall_hgps_plot

ggsave(
  "hallmarks_hgps_umk57_dmso_dotplot.png",
  path = file.path("output_data", "plots", "gsea_alternative"),
  width = 7,
  height = 8, 
  dpi = 600
)


# GSEA DS UMK57 vs DMSO Hallmark ----------------------------------------------

df_ds_of_interest <- df_ds_of_interest %>% 
  dplyr::arrange(desc(log2FoldChange))

ds_of_interest_ranked <- df_ds_of_interest$log2FoldChange
names(ds_of_interest_ranked) <- df_ds_of_interest$gene
head(ds_of_interest_ranked)


saveRDS(
  ds_of_interest_ranked,
  file.path("output_data", "results", "gsea_alternative", "ds_umk57_vs_dmso_ranked_list_for_GSEA.rds")
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
  file.path("output_data", "results", "gsea_alternative", "ds_umk57_dmso_GSEA_hallmarks.csv")
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

hall_ds_plot <- ggplot(df_Hall_ds_of_interest_005,
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
       title = "DS UMK57 vs DMSO (Paired Analysis)",
       subtitle = "Hallmark Gene Sets") +
  geom_vline(xintercept = 0, linetype = "dashed")
hall_ds_plot

ggsave(
  "hallmarks_ds_umk57_dmso_dotplot.png",
  path = file.path("output_data", "plots", "gsea_alternative"),
  width = 7,
  height = 8, dpi = 600
)

hall_hgps_plot + hall_ds_plot

ggsave(
  "hallmarks_hgps_ds_umk57_dmso_dotplot_FOR_REBUTTAL.png",
  path = file.path("output_data", "plots", "gsea_alternative"),
  width = 12,
  height = 7, dpi = 600
)



# Write Session Info and Log Outputs --------------------------------------

cat(
  paste0("\n## Step completed\n",
         "Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
         "Finished 08_deseq2_alternative\n"),
  file = "PROGRESS_LOG.md",
  append = TRUE
)

dir.create(
  file.path("session_info_logs")
)

writeLines(capture.output(sessionInfo()), 
           file.path("session_info_logs", 
                     "08_deseq2_alternative_sessionInfo.txt")
)