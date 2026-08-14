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

cat(
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
         "Started 02_gsea_combined_healthy_controls\n"),
  file = "PROGRESS_LOG.md",
  append = TRUE
)


# Load Counts and Metadata ----------------------------------------------------

# Load correct sample info and counts for DS, HGPS and Healthy
# samples obtained from IonTorrent Suite:

counts_matrix_correct <- read.csv(file.path("input_data", "counts_matrix.csv"),
                                  row.names = 1)

sample_info_correct <- read.csv(file.path("input_data", "sample_info.csv"),
                                row.names = 1)





## DESeq2 HGPS vs Healthy ------------------------------------------------------

sample_info_correct_hgps_neo <- sample_info_correct %>% 
  dplyr:: filter(Disease %in% c("Healthy", "HGPS")) %>% 
  dplyr::filter(Treatment == "DMSO")
sample_info_correct_hgps_neo

counts_matrix_correct_hgps_neo <- counts_matrix_correct[,row.names(sample_info_correct_hgps_neo)]

colnames(counts_matrix_correct_hgps_neo) == rownames(sample_info_correct_hgps_neo)

# Should be TRUE.

dds_hgps_neo <- DESeqDataSetFromMatrix(
  countData = counts_matrix_correct_hgps_neo,
  colData = sample_info_correct_hgps_neo,
  design = ~ Disease)

# Design wont matter at this stage though.

smallestGroupSize <- nrow(sample_info_correct_hgps_neo) / 2 

# Half of samples.

keep <- rowSums(counts(dds_hgps_neo) >=5) >= smallestGroupSize 

# Recomended in vignette.

dds_hgps_neo <- dds_hgps_neo[keep,]

dds_hgps_neo

dds_hgps_neo <- DESeq(dds_hgps_neo)

res_hgps_neo <- results(dds_hgps_neo, contrast = 
                        c("Disease", "HGPS", "Healthy"), 
                        cooksCutoff = FALSE, 
                        independentFiltering = FALSE) 


df_hgps_neo <- as.data.frame(res_hgps_neo)

# GSEA HGPS vs Healthy Hallmark -----------------------------------------------

head(df_hgps_neo)

df_hgps_neo <- df_hgps_neo %>% 
  dplyr::arrange(desc(log2FoldChange))
head(df_hgps_neo)

df_hgps_neo$gene <- rownames(df_hgps_neo)

write.csv(
  df_hgps_neo,
  file.path("output_data", "results", "deseq2", "deseq2_hgps_vs_healthy.csv")
)

hgps_healthy_ranked <- df_hgps_neo$log2FoldChange
names(hgps_healthy_ranked) <- df_hgps_neo$gene
head(hgps_healthy_ranked, 10)


hall <- msigdbr::msigdbr(
  species = "Homo sapiens",
  collection = "H"
)

hall_t2g <- hall %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()

hgps_healthy_hall <- clusterProfiler::GSEA(
  hgps_healthy_ranked,
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

df_hgps_healthy_hall <- as.data.frame(hgps_healthy_hall)

write.csv(
  df_hgps_healthy_hall,
  file.path("output_data", "results", "gsea", 
            "hgps_vs_healthy_GSEA_hallmarks.csv")
)

df_hgps_healthy_hall_005 <- df_hgps_healthy_hall %>% 
  dplyr::filter(p.adjust < 0.05)

df_hgps_healthy_hall_005$Description <- gsub(
  "HALLMARK_", "", df_hgps_healthy_hall_005$Description
)

df_hgps_healthy_hall_005 <- df_hgps_healthy_hall_005 %>% 
  dplyr::mutate(direction = case_when(
    NES > 0 ~ "Up",
    NES < 0 ~ "Down"
  ))

cols <- c(
  "Up" = "red",
  "Down" = "blue"
)

ggplot(df_hgps_healthy_hall_005,
       aes(x = NES, y = forcats::fct_reorder(Description, NES))) +
  geom_segment(aes(xend = 0, yend = Description)) +
  geom_point(aes(fill = direction, size = p.adjust), shape = 21) +
  scale_fill_manual(values = cols, name = "Direction") +
  scale_size_continuous(range = c(8,3), transform = "log10", 
                        name = "padjust") +
  theme_bw(base_size = 13) +
  xlim(-4, 3) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  labs(x = "NES", y = "Gene Sets\n",
       title = "HGPS vs Healthy",
       subtitle = "Hallmark Gene Sets") +
  geom_vline(xintercept = 0, linetype = "dashed")




# DESeq2 DS vs Healthy --------------------------------------------------------

sample_info_correct_ds_neo <- sample_info_correct %>% 
  dplyr::filter(Treatment == "DMSO") %>% 
  dplyr::filter(Disease %in% c("Healthy", "DS"))

counts_matrix_correct_ds_neo <- counts_matrix_correct[,row.names(sample_info_correct_ds_neo)]


dds_ds_neo <- DESeqDataSetFromMatrix(
  countData = counts_matrix_correct_ds_neo,
  colData = sample_info_correct_ds_neo,
  design = ~ Disease) 

# Design wont matter at this stage though.

smallestGroupSize <- nrow(sample_info_correct_ds_neo)/2 

# Half of samples, can be changed later.

keep <- rowSums(counts(dds_ds_neo) >=5) >= smallestGroupSize 

# Recomended in vignette.

dds_ds_neo <- dds_ds_neo[keep,]

dds_ds_neo

dds_ds_neo <- DESeq(dds_ds_neo)

res_ds_neo <- results(dds_ds_neo, contrast = 
                      c("Disease", "DS", "Healthy"), 
                      cooksCutoff = FALSE, 
                      independentFiltering = FALSE) 


df_ds_neo <- as.data.frame(res_ds_neo)

# GSEA DS vs Healthy Hallmark -------------------------------------------------

head(df_ds_neo)

df_ds_neo$gene <- rownames(df_ds_neo)

df_ds_neo <- df_ds_neo %>% 
  dplyr::arrange(desc(log2FoldChange))
head(df_ds_neo)

write.csv(df_ds_neo, 
          file.path("output_data", "results", "deseq2", "deseq2_ds_vs_healthy.csv")
          )

ds_healthy_ranked <- df_ds_neo$log2FoldChange
names(ds_healthy_ranked) <- df_ds_neo$gene
head(ds_healthy_ranked, 10)


hall <- msigdbr::msigdbr(
  species = "Homo sapiens",
  collection = "H"
)

hall_t2g <- hall %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()

ds_healthy_hall <- clusterProfiler::GSEA(
  ds_healthy_ranked,
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

df_ds_healthy_hall <- as.data.frame(ds_healthy_hall)

write.csv(
  df_ds_healthy_hall,
  file.path("output_data", "results", "gsea", "ds_vs_healthy_GSEA_hallmarks.csv")
)

df_ds_healthy_hall_005 <- df_ds_healthy_hall %>% 
  dplyr::filter(p.adjust < 0.05)

df_ds_healthy_hall_005$Description <- gsub(
  "HALLMARK_", "", df_ds_healthy_hall_005$Description
)

df_ds_healthy_hall_005 <- df_ds_healthy_hall_005 %>% 
  dplyr::mutate(direction = case_when(
    NES > 0 ~ "Up",
    NES < 0 ~ "Down"
  ))

cols <- c(
  "Up" = "red",
  "Down" = "blue"
)

ggplot(df_ds_healthy_hall_005,
       aes(x = NES, y = forcats::fct_reorder(Description, NES))) +
  geom_segment(aes(xend = 0, yend = Description)) +
  geom_point(aes(fill = direction, size = p.adjust), shape = 21) +
  scale_fill_manual(values = cols, name = "Direction") +
  scale_size_continuous(range = c(8,3), transform = "log10", 
                        name = "padjust") +
  theme_bw(base_size = 13) +
  xlim(-4,3) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  labs(x = "NES", y = "Gene Sets\n",
       title = "DS vs Healthy",
       subtitle = "Hallmark Gene Sets") +
  geom_vline(xintercept = 0, linetype = "dashed")




## Combined Dotplot Hallmark ---------------------------------------------------

df_ds_healthy_hall[1:5,2:5]
df_hgps_healthy_hall[1:5,2:5]

df_hgps_umk57_hall <-  read.csv(
  file.path("output_data", "results", "gsea", "hgps_umk57_dmso_GSEA_hallmarks.csv"),
  row.names = 1)

df_Hall_ds_of_interest <- read.csv(
  file.path("output_data", "results", "gsea", "ds_umk57_dmso_GSEA_hallmarks.csv"),
  row.names = 1)

df_ds_healthy_hall$comparison <- "DS vs Healthy"
df_hgps_healthy_hall$comparison <- "HGPS vs Healthy"
df_hgps_umk57_hall$comparison <- "HGPS UMK57 vs DMSO"
df_Hall_ds_of_interest$comparison <- "DS UMK57 vs DMSO"

colnames(df_ds_healthy_hall) == colnames(df_Hall_ds_of_interest)

dodgeplot_ds <- rbind(df_ds_healthy_hall,
                      df_Hall_ds_of_interest)

dodgeplot_ds$Description <- gsub(
  "HALLMARK_",
  "",
  dodgeplot_ds$Description
)

dodgeplot_ds <- dodgeplot_ds %>% 
  dplyr::mutate(NES_NEO = case_when(
    comparison == "DS vs Healthy" ~ NES,
    .default = 0
  ))

dodgeplot_ds <- dodgeplot_ds %>% 
  dplyr::arrange(desc(NES_NEO))

dodgeplot_ds_005 <- dodgeplot_ds %>% 
  dplyr::filter(p.adjust < 0.05) %>% 
  dplyr::arrange(desc(NES_NEO))

unique(dodgeplot_ds_005$Description)

dodgeplot_hgps <- rbind(df_hgps_healthy_hall,
                        df_hgps_umk57_hall)

dodgeplot_hgps$Description <- gsub(
  "HALLMARK_",
  "",
  dodgeplot_hgps$Description
)

dodgeplot_hgps <- dodgeplot_hgps %>% 
  dplyr::mutate(NES_NEO = case_when(
    comparison == "HGPS vs Healthy" ~ NES,
    .default = 0
  ))

dodgeplot_hgps <- dodgeplot_hgps %>% 
  dplyr::arrange(desc(NES_NEO))

dodgeplot_hgps_005 <- dodgeplot_hgps %>% 
  dplyr::filter(p.adjust < 0.05) %>% 
  dplyr::arrange(desc(NES_NEO))

unique(dodgeplot_hgps_005$Description)


dodgeplot_hall_combined <- rbind(
  dodgeplot_ds_005,
  dodgeplot_hgps_005
)

dodgeplot_hall_combined$comparison <- factor(
  dodgeplot_hall_combined$comparison,
  levels = c("HGPS vs Healthy",
             "HGPS UMK57 vs DMSO",
             "DS vs Healthy",
             "DS UMK57 vs DMSO")
)

hall_factors <- dodgeplot_hall_combined %>% 
  filter(NES_NEO != 0) %>% 
  arrange(desc(NES_NEO)) %>% 
  pull(Description) %>% unique() %>% as.vector()

hall_factors_2 <- dodgeplot_hall_combined %>% 
  pull(Description) %>% unique() %>% as.vector()

hall_factors_2 <- setdiff(hall_factors_2, hall_factors)

hall_factors_final <- c(hall_factors_2, hall_factors)

dodgeplot_hall_combined$Description <- factor(
  dodgeplot_hall_combined$Description,
  levels = hall_factors_final
)

ggplot(dodgeplot_hall_combined, aes(
  x = comparison, 
  y = Description, 
  group = comparison)) +
  geom_point(shape = 21, aes(fill = NES, size = p.adjust)) +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0,
                       breaks = c(-3:2)) +
  scale_size_continuous(range = c(7, 2.5), transform = "log10",
                        name = "padjust") +
  theme_bw(base_size = 14) +
  labs(y = "Gene Sets\n", x = "",
       title = "Hallmark Gene Sets") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.grid.major = element_line(size = 0.25)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  guides(fill = guide_colorbar(
    frame.colour = "black",
    frame.linewidth = 0.25,
    ticks.linewidth = 0.25,
    ticks.colour = "black")) +
  geom_vline(xintercept = 2.5)

ggsave("GSEA_hallmarks_dotplot_all_comparisons.pdf",
       path = file.path("output_data", "plots", "gsea"),
       width = 7,
       height = 9)



## GSEA HGPS vs Healthy DNA Repair ---------------------------------------------

mm_BP_sets <- msigdbr::msigdbr(
  species = "Homo sapiens")
mm_BP_sets

msigdbr_t2g = mm_BP_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% 
  as.data.frame()
msigdbr_t2g

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


hgps__vs_healthy_dna_repair <- clusterProfiler::GSEA(
  hgps_healthy_ranked,
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
  nPerm = 20000 # Since it shows warning, to make sure all terms are included.
)

df_hgps__vs_healthy_dna_repair <- as.data.frame(hgps__vs_healthy_dna_repair)

write.csv(
  df_hgps__vs_healthy_dna_repair,
  file.path("output_data", "results", "gsea", "hgps_vs_healthy_GSEA_DNA_Repair.csv")
)

df_hgps__vs_healthy_dna_repair_005 <- df_hgps__vs_healthy_dna_repair %>% 
  dplyr::filter(p.adjust < 0.05)

df_hgps__vs_healthy_dna_repair_005$Description <- gsub("_", " ", df_hgps__vs_healthy_dna_repair_005$Description)

df_hgps__vs_healthy_dna_repair_005$Description <- stringr::str_wrap(df_hgps__vs_healthy_dna_repair_005$Description,
                                                                width = 30)

df_hgps__vs_healthy_dna_repair_005$Description

df_hgps__vs_healthy_dna_repair_005 <- df_hgps__vs_healthy_dna_repair_005 %>% 
  dplyr::mutate(direction = case_when(
    NES > 0 ~ "Up",
    NES < 0 ~ "Down"
  ))

cols = c(
  "Up" = "red",
  "Down" = "blue"
)

ggplot(df_hgps__vs_healthy_dna_repair_005,
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
       title = "HGPS vs Healthy",
       subtitle = "DNA Repair Gene Sets") +
  geom_vline(xintercept = 0, linetype = "dashed")


# GSEA DS vs Healthy DNA Repair -----------------------------------------------

ds__vs_healthy_dna_repair <- clusterProfiler::GSEA(
  ds_healthy_ranked,
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


df_ds__vs_healthy_dna_repair <- as.data.frame(ds__vs_healthy_dna_repair)

write.csv(
  df_ds__vs_healthy_dna_repair,
  file.path("output_data", "results", "gsea", "ds_vs_healthy_GSEA_DNA_Repair.csv")
)

df_ds__vs_healthy_dna_repair_005 <- df_ds__vs_healthy_dna_repair %>% 
  dplyr::filter(p.adjust < 0.05)

df_ds__vs_healthy_dna_repair_005$Description <- gsub("_", " ", df_ds__vs_healthy_dna_repair_005$Description)

df_ds__vs_healthy_dna_repair_005$Description <- stringr::str_wrap(df_ds__vs_healthy_dna_repair_005$Description,
                                                              width = 30)

df_ds__vs_healthy_dna_repair_005$Description

df_ds__vs_healthy_dna_repair_005 <- df_ds__vs_healthy_dna_repair_005 %>% 
  dplyr::mutate(direction = case_when(
    NES > 0 ~ "Up",
    NES < 0 ~ "Down"
  ))

ggplot(df_ds__vs_healthy_dna_repair_005,
       aes(x = NES, y = forcats::fct_reorder(Description, NES))) +
  geom_segment(aes(xend = 0, yend = Description)) +
  geom_point(aes(fill = direction, size = p.adjust), shape = 21) +
  scale_fill_manual(values = cols, name = "Direction") +
  scale_size_continuous(range = c(8,3), transform = "log10", name = "padjust") +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  xlim(-4,3) +
  labs(x = "NES", y = "Gene Sets\n",
       title = "DS vs Healthy",
       subtitle = "DNA Repair Gene Sets") +
  geom_vline(xintercept = 0, linetype = "dashed")



# Divergent Dotplots DNA Repair -----------------------------------------------

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


df_t21_dna_repair_005 <- read.csv(
  file.path("output_data", "results", "gsea", "df_ds_umk57_vs_dmso_GSEA_DNA_repair.csv"),
  row.names = 1
  ) |> 
  dplyr::mutate(direction = case_when(
    NES > 0 & p.adjust < 0.05 ~ "Up",
    NES < 0 & p.adjust < 0.05 ~ "Down",
    .default = "No Change"
  )) |> 
  dplyr::filter(p.adjust < 0.05) |> 
  dplyr::mutate(comparison = "DS UMK57 vs DMSO")

df_hgps_dna_repair_005 <- read.csv(
  file.path("output_data", "results", "gsea", "df_hgps_umk57_vs_dmso_GSEA_DNA_repair.csv"),
  row.names = 1
) |> 
  dplyr::mutate(direction = case_when(
    NES > 0 & p.adjust < 0.05 ~ "Up",
    NES < 0 & p.adjust < 0.05 ~ "Down",
    .default = "No Change"
  )) |> 
  dplyr::filter(p.adjust < 0.05) |> 
  dplyr::mutate(comparison = "HGPS UMK57 vs DMSO")

df_hgps__vs_healthy_dna_repair_005

df_ds__vs_healthy_dna_repair_005


df_hgps__vs_healthy_dna_repair_005$comparison <- "HGPS vs Healthy Control"
df_ds__vs_healthy_dna_repair_005$comparison <- "DS vs Healthy Control"

colnames(df_t21_dna_repair_005)
colnames(df_ds__vs_healthy_dna_repair_005)

ds_dna_repair_divergent <- rbind(
  df_t21_dna_repair_005,
  df_ds__vs_healthy_dna_repair_005
)

ds_dna_repair_divergent$comparison <- factor(
  ds_dna_repair_divergent$comparison,
  levels = c("DS vs Healthy Control",
             "DS UMK57 vs DMSO")
)

ds_dna_repair_divergent$Description <- str_wrap(
  ds_dna_repair_divergent$Description,
  width = 900) ## Make sure no wraping is made now.

df_ds__vs_healthy_dna_repair_005$Description <- str_wrap(
  df_ds__vs_healthy_dna_repair_005$Description,
  width = 900) ## Make sure no wraping is made now.

ds_dna_repair_divergent <- ds_dna_repair_divergent %>% 
  left_join(new_dna_repair_terms, join_by(ID))

ds_order <- ds_dna_repair_divergent %>% 
  dplyr::filter(comparison == "DS vs Healthy Control") %>% 
  dplyr::arrange(desc(NES)) %>% 
  dplyr::pull(new)

ds_dna_repair_divergent$new <- factor(
  ds_dna_repair_divergent$new,
  levels = ds_order
)

ds_healthy_dna_rep_plot <- ggplot(ds_dna_repair_divergent, aes(
  x = comparison, 
  y = new, 
  group = comparison)) +
  geom_point(shape = 21, aes(fill = NES, size = p.adjust)) +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0) +
  scale_size_continuous(range = c(6, 2.5), transform = "log10",
                        name = "padjust") +
  theme_bw(base_size = 14) +
  labs(y = "Gene Sets\n", x = "",
       title = "DS - DNA Repair Gene Sets") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90),
        panel.grid.major = element_line(size = 0.25)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  guides(fill = guide_colorbar(
    frame.colour = "black",
    frame.linewidth = 0.25,
    ticks.linewidth = 0.25,
    ticks.colour = "black"))

ds_healthy_dna_rep_plot

hgps_dna_repair_divergent <- rbind(
  df_hgps_dna_repair_005,
  df_hgps__vs_healthy_dna_repair_005
)

hgps_dna_repair_divergent$comparison <- factor(
  hgps_dna_repair_divergent$comparison,
  levels = c("HGPS vs Healthy Control",
             "HGPS UMK57 vs DMSO")
)

hgps_dna_repair_divergent$Description <- str_wrap(
  hgps_dna_repair_divergent$Description,
  width = 900) ## Make sure no wraping is made now.

df_hgps__vs_healthy_dna_repair_005$Description <- str_wrap(
  df_hgps__vs_healthy_dna_repair_005$Description,
  width = 900) ## Make sure no wraping is made now.

hgps_dna_repair_divergent <- hgps_dna_repair_divergent %>% 
  left_join(new_dna_repair_terms, join_by(ID))

hgps_order <- hgps_dna_repair_divergent %>% 
  dplyr::filter(comparison == "HGPS vs Healthy Control") %>% 
  dplyr::arrange(desc(NES)) %>% 
  dplyr::pull(new)

hgps_dna_repair_divergent$new <- factor(
  hgps_dna_repair_divergent$new,
  levels = hgps_order
)

all_dna_repair_divergent <- rbind(
  hgps_dna_repair_divergent,
  ds_dna_repair_divergent
)

all_dna_repair_divergent$comparison <- factor(
  all_dna_repair_divergent$comparison,
  levels = c("HGPS vs Healthy Control",
             "HGPS UMK57 vs DMSO",
             "DS vs Healthy Control",
             "DS UMK57 vs DMSO")
)

ggplot(all_dna_repair_divergent, aes(
  x = comparison, 
  y = forcats::fct_reorder(new, setSize), 
  group = comparison)) +
  geom_point(shape = 21, aes(fill = NES, size = p.adjust)) +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0) +
  scale_size_continuous(range = c(6, 2.5), transform = "log10",
                        name = "padjust") +
  theme_bw(base_size = 14) +
  labs(y = "Gene Sets\n", x = "") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        panel.grid.major = element_line(size = 0.25)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  guides(fill = guide_colorbar(
    frame.colour = "black",
    frame.linewidth = 0.25,
    ticks.linewidth = 0.25,
    ticks.colour = "black"))

ggsave("DNA_Repair_GSEA_dotplot_all_comparisons.pdf",
       path = file.path("output_data", "plots", "gsea"),
       width = 9,
       height = 10)



## GSEA DNA Repair Running Plots -----------------------------------------------

## DSB Repair Only:

hgps_dna_repair <- readRDS(
  file.path("output_data", "results", "gsea", "hgps_dna_repair_gsea.rds")
  )

t21_dna_repair <- readRDS(
  file.path("output_data", "results", "gsea", "ds_dna_repair_gsea.rds")
)

df_hgps_dna_repair <- as.data.frame(hgps_dna_repair)
df_t21_dna_repair <- as.data.frame(t21_dna_repair)


dsb_line <- match("GOBP_DOUBLE_STRAND_BREAK_REPAIR",
                  df_hgps_dna_repair$Description)

dna_rep_cols <- c(
  "HGPS UMK57 vs DMSO" = "#834746",
  "DS UMK57 vs DMSO" = "lightskyblue"
)

p1 <- enrichplot::gseaplot2(hgps_dna_repair, geneSetID = dsb_line,
                            title = "DSB Repair HGPS UMK57 vs DMSO",
                            base_size = 15,
                            color = "#834746",
                            rel_heights = c(1.5,0.2,0.5),
                            subplots = 1:2,
                            pvalue_table = F,
                            ES_geom = "line")
p1


df_hgps__vs_healthy_dna_repair

dsb_line <- match("GOBP_DOUBLE_STRAND_BREAK_REPAIR",
                  df_hgps__vs_healthy_dna_repair$Description)

p2 <- enrichplot::gseaplot2(hgps__vs_healthy_dna_repair, geneSetID = dsb_line,
                            title = "DSB Repair HGPS vs Healthy",
                            base_size = 15,
                            color = "#834746",
                            rel_heights = c(1.5,0.2,0.5),
                            subplots = 1:2,
                            pvalue_table = F,
                            ES_geom = "line")
p2


df_ds__vs_healthy_dna_repair

dsb_line <- match("GOBP_DOUBLE_STRAND_BREAK_REPAIR",
                  df_ds__vs_healthy_dna_repair$Description)

p3 <- enrichplot::gseaplot2(ds__vs_healthy_dna_repair, geneSetID = dsb_line,
                            title = "DSB Repair DS vs Healthy",
                            base_size = 15,
                            color = "lightskyblue",
                            rel_heights = c(1.5,0.2,0.5),
                            subplots = 1:2,
                            pvalue_table = F,
                            ES_geom = "line")
p3


df_t21_dna_repair

dsb_line <- match("GOBP_DOUBLE_STRAND_BREAK_REPAIR",
                  df_t21_dna_repair$Description)

p4 <- enrichplot::gseaplot2(t21_dna_repair, geneSetID = dsb_line,
                            title = "DSB Repair DS UMK57 vs DMSO",
                            base_size = 15,
                            color = "lightskyblue",
                            rel_heights = c(1.5,0.2,0.5),
                            subplots = 1:2,
                            pvalue_table = F,
                            ES_geom = "line")
p4


aplot::plot_list(p2, p1, p3, p4, nrow = 1) # This is the desired order.

ggsave("gseaplot2_gobp_dsb_repair_all_gsea_plots.pdf",
       path = file.path("output_data", "plots", "gsea"),
       width = 26,
       height = 5)


## GSEA SenMayo ----------------------------------------------------------------

## Ranked Lists of Interest

hgps_umk57_vs_dmso_ranked <- readRDS(
  file.path("output_data", "results", "gsea", "hgps_umk57_vs_dmso_ranked_list_for_GSEA.rds")
)

ds_umk57_vs_dmso_ranked <- readRDS(
  file.path("output_data", "results", "gsea", "ds_umk57_vs_dmso_ranked_list_for_GSEA.rds")
)


head(ds_healthy_ranked)
head(hgps_healthy_ranked)
head(ds_umk57_vs_dmso_ranked)
head(hgps_umk57_vs_dmso_ranked)

senmayo <- read.table(file.path("input_data", "senmayo.txt"),
                      header = T) %>% 
  dplyr::filter(gene_set == "SenMayo")

## DS vs H

ds_healthy_senmayo <- clusterProfiler::GSEA(
  ds_healthy_ranked,
  exponent = 1,
  minGSSize = 0,
  maxGSSize = 1000,
  eps = 1e-50,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = senmayo,
  verbose = TRUE,
  seed = TRUE,
  by = "fgsea",
)

## HGPS vs H

hgps_healthy_senmayo <- clusterProfiler::GSEA(
  hgps_healthy_ranked,
  exponent = 1,
  minGSSize = 0,
  maxGSSize = 1000,
  eps = 1e-50,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = senmayo,
  verbose = TRUE,
  seed = TRUE,
  by = "fgsea",
)


## HGPS UMK57 vs DMSO

hgps_umk57_dmso_senmayo <- clusterProfiler::GSEA(
  hgps_umk57_vs_dmso_ranked,
  exponent = 1,
  minGSSize = 0,
  maxGSSize = 1000,
  eps = 1e-50,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = senmayo,
  verbose = TRUE,
  seed = TRUE,
  by = "fgsea",
)

## DS UMK57 vs DMSO

ds_umk57_dmso_senmayo <- clusterProfiler::GSEA(
  ds_umk57_vs_dmso_ranked,
  exponent = 1,
  minGSSize = 0,
  maxGSSize = 1000,
  eps = 1e-50,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = senmayo,
  verbose = TRUE,
  seed = TRUE,
  by = "fgsea",
)

## GSEA SenMayo Running Plots --------------------------------------------------

p2 <- enrichplot::gseaplot2(hgps_umk57_dmso_senmayo, geneSetID = 1,
                            title = "SenMayo HGPS UMK57 vs DMSO",
                            base_size = 15,
                            color = "#834746",
                            rel_heights = c(1.5,0.2,0.5),
                            subplots = 1:2,
                            pvalue_table = F,
                            ES_geom = "line")
p2


p1 <- enrichplot::gseaplot2(hgps_healthy_senmayo, geneSetID = 1,
                            title = "SenMayo HGPS vs Healthy",
                            base_size = 15,
                            color = "#834746",
                            rel_heights = c(1.5,0.2,0.5),
                            subplots = 1:2,
                            pvalue_table = F,
                            ES_geom = "line")
p1

p3 <- enrichplot::gseaplot2(ds_healthy_senmayo, 
                            geneSetID = 1,
                            title = "SenMayo DS vs Healthy",
                            base_size = 15,
                            color = "lightskyblue",
                            rel_heights = c(1.5,0.2,0.5),
                            subplots = 1:2,
                            pvalue_table = F,
                            ES_geom = "line")
p3


p4 <- enrichplot::gseaplot2(ds_umk57_dmso_senmayo, geneSetID = 1,
                            title = "SenMayo DS UMK57 vs DMSO",
                            base_size = 15,
                            color = "lightskyblue",
                            rel_heights = c(1.5,0.2,0.5),
                            subplots = 1:2,
                            pvalue_table = F,
                            ES_geom = "line")
p4


aplot::plot_list(p1, p2, p3, p4, nrow = 1)

ggsave("gseaplot2_senmayo_all.pdf",
       path = file.path("output_data", "plots", "gsea"),
       width = 26,
       height = 5)


senmayo_results <- rbind(
  (ds_healthy_senmayo@result |> 
     dplyr::mutate(comparison = "DS_vs_Healthy")
   ),
  (ds_umk57_dmso_senmayo@result |> 
     dplyr::mutate(comparison = "DS_UMK57_vs_DMSO")
   ),
  (hgps_healthy_senmayo@result |> 
     dplyr::mutate(comparison = "HGPS_vs_Healthy")
   ),
  (hgps_umk57_dmso_senmayo@result |> dplyr::mutate(comparison = "HGPS_UMK57_vs_DMSO")
   )
)

write.csv(
  senmayo_results,
  file.path("output_data", "results", "gsea", "senmayo_results.csv")
)


## SenMayo DEGs in NEO + HGPS Full Overlap -------------------------------------


sample_info_hgps_umk57_neo <- sample_info_correct %>% 
  dplyr::filter(Disease %in% c("HGPS", "Healthy"))

counts_matrix_healthy_hgps_umk57 <- counts_matrix_correct[,row.names(
  sample_info_hgps_umk57_neo)]

colnames(counts_matrix_healthy_hgps_umk57) == row.names(sample_info_hgps_umk57_neo)

dds_healthy_hgps_umk57 <- DESeqDataSetFromMatrix(
  countData = counts_matrix_healthy_hgps_umk57,
  colData = sample_info_hgps_umk57_neo,
  design = ~ Individual) ## Makes no difference for PCA

colnames(dds_healthy_hgps_umk57)

## Design wont matter at this stage though.

smallestGroupSize <- 5 # Half of samples, can be changed later.

keep <- rowSums(counts(dds_healthy_hgps_umk57) >=5) >= smallestGroupSize 

# Recomended in vignette:

dds_healthy_hgps_umk57 <- dds_healthy_hgps_umk57[keep,]

dds_healthy_hgps_umk57

dds_healthy_hgps_umk57 <- DESeq(dds_healthy_hgps_umk57)

dds_healthy_hgps_umk57

sample_info_hgps_umk57_neo

vsd <- vst(dds_healthy_hgps_umk57, blind = TRUE)

degs_hgps_vs_neo <- df_hgps_neo %>% 
  dplyr::filter(padj < 0.05) %>% 
  dplyr::pull(gene)

degs_hgps <- read.csv(
  file.path("output_data", "results", "deseq2", "deseq2_hgps_umk57_vs_dmso.csv")
  ) |> 
  dplyr::filter(padj < 0.05) |> 
  dplyr::pull(gene)

senmayo_genes <- senmayo$gene_symbol

senamyo_degs_hgps_vs_neo <- intersect(
  degs_hgps_vs_neo,
  senmayo_genes
)

senamyo_degs_hgps_vs_neo <- intersect(
  senamyo_degs_hgps_vs_neo,
  degs_hgps
)

# Create heatmap with genes of interest:

vsd_df <- assay(vsd)
vsd_zscore <- scale(t(vsd_df))

hgps_correct_order <- sample_info_hgps_umk57_neo |> 
  dplyr::arrange(desc(Individual), Treatment) |> 
  dplyr::mutate(
    Individual = case_when(
      Individual == "HGPS 169" ~ "HGPS",
      Individual == "NEO" ~ "Healthy"
    )
  )

vsd_zscore_goi <- vsd_zscore %>% as.data.frame()
vsd_zscore_goi <- vsd_zscore_goi[row.names(hgps_correct_order),senamyo_degs_hgps_vs_neo]
vsd_zscore_goi <- na.omit(vsd_zscore_goi)


col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_fun(seq(-2, 2)) ## Necessary to assign values to colors

ha1 = rowAnnotation(Individual = hgps_correct_order$Individual,
                    Treatment = hgps_correct_order$Treatment,
                    col = list(
                      Individual = c(
                        "Healthy" = "#636363",
                        "HGPS" = "#831a10"
                      ),
                      Treatment = c(
                        "DMSO" = "gray",
                        "UMK57" = "cyan"
                      )
                    ))



group_split <- hgps_correct_order$Individual      

group_split <- factor(group_split, levels = c("Healthy", "HGPS"))


ht <- Heatmap(
  vsd_zscore_goi,
  cluster_columns = T,
  cluster_rows = F,
  show_row_names = F,
  row_split = group_split,
  column_km = 2,
  left_annotation = ha1,
  col = col_fun,
  border_gp = gpar(col = "black", lty = 1),
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 12),
  column_names_rot = 90,
  heatmap_legend_param = list(
    at = c(-2, 0, 2),
    labels = c("-2", "0", "2"),
    title = "Row Z-Score",
    border = "black",
    title_position = "leftcenter-rot"))
ht

dir.create(
  file.path("output_data", "plots", "heatmaps")
)

pdf(
  file.path("output_data", "plots", "heatmaps", "heatmap_senmayo_hgps.pdf"),
  width = 6,
  height = 3
)

draw(ht)

dev.off()


## SenMayo DEGs in NEO + DS Full Overlap ---------------------------------------

sample_info_ds_umk57_neo <- sample_info_correct %>% 
  dplyr::filter(Disease %in% c("DS", "Healthy"))

counts_matrix_healthy_ds_umk57 <- counts_matrix_correct[,row.names(
  sample_info_ds_umk57_neo)]

colnames(counts_matrix_healthy_ds_umk57) == row.names(sample_info_ds_umk57_neo)

dds_healthy_ds_umk57 <- DESeqDataSetFromMatrix(
  countData = counts_matrix_healthy_ds_umk57,
  colData = sample_info_ds_umk57_neo,
  design = ~ Individual) # Makes no difference for PCA.

colnames(dds_healthy_ds_umk57)

## Design wont matter at this stage though.

smallestGroupSize <- 5
keep <- rowSums(counts(dds_healthy_ds_umk57) >=5) >= smallestGroupSize 

# Recomended in vignette.

dds_healthy_ds_umk57 <- dds_healthy_ds_umk57[keep,]

dds_healthy_ds_umk57

dds_healthy_ds_umk57 <- DESeq(dds_healthy_ds_umk57)

dds_healthy_ds_umk57

sample_info_ds_umk57_neo

vsd <- vst(dds_healthy_ds_umk57, blind = TRUE) ## Blind to ignore conditions

degs_ds_vs_neo <- df_ds_neo %>% 
  dplyr::filter(padj < 0.05) %>% 
  dplyr::pull(gene)

degs_ds <- read.csv(
  file.path("output_data", "results", "deseq2", "deseq2_ds_umk57_vs_dmso.csv")
) |> 
  dplyr::filter(padj < 0.05) |> 
  dplyr::pull(gene)

senmayo_degs_ds_vs_neo <- intersect(
  degs_ds_vs_neo,
  senmayo_genes
)

senmayo_degs_ds_vs_neo <- intersect(
  senmayo_degs_ds_vs_neo,
  degs_ds
)

# Create heatmap with genes of interest:

vsd_df <- assay(vsd)
vsd_zscore <- scale(t(vsd_df)) %>% as.data.frame()

ds_correct_order <- sample_info_correct %>% 
  dplyr::filter(Disease %in% c("Healthy", "DS")) %>% 
  dplyr::arrange(Treatment, desc(Individual)) |> 
  dplyr::mutate(
    Individual = case_when(
      Disease == "DS" ~ "DS",
      Individual == "NEO" ~ "Healthy"
    )
  )
ds_correct_order

vsd_zscore <- vsd_zscore[row.names(ds_correct_order),]  

vsd_zscore_goi <- vsd_zscore %>% as.data.frame()
vsd_zscore_goi <- vsd_zscore_goi[,senmayo_degs_ds_vs_neo]


col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
col_fun(seq(-2, 2)) ## Necessary to assign values to colors



ha1 = rowAnnotation(Individual = ds_correct_order$Individual,
                    Treatment = ds_correct_order$Treatment,
                    col = list(
                      Individual = c(
                        "Healthy" = "#636363",
                        "DS" = "#95aace"
                      ),
                      Treatment = c(
                        "DMSO" = "gray",
                        "UMK57" = "cyan"
                      )
                    ))



group_split <- ds_correct_order$Disease      

group_split <- factor(group_split, levels = c("Healthy", "DS"))

ht_2 <- Heatmap(
  vsd_zscore_goi,
  cluster_columns = T,
  cluster_rows = F,
  show_row_names = F,
  row_split = group_split,
  column_km = 2,
  left_annotation = ha1,
  col = col_fun,
  border_gp = gpar(col = "black", lty = 1),
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 11),
  column_names_rot = 90,
  heatmap_legend_param = list(
    at = c(-2, 0, 2),
    labels = c("-2", "0", "2"),
    title = "Row Z-Score",
    border = "black",
    title_position = "leftcenter-rot"))
ht_2

pdf(
  file.path("output_data",
            "plots",
            "heatmaps",
            "heatmap_senmayo_ds.pdf"),
  width = 6,
  height = 3
)

draw(ht_2)

dev.off()



# Merged Heatmap:

ht + ht_2


pdf(
  file.path("output_data",
            "plots",
            "heatmaps",
            "heatmap_senmayo_hgps_and_ds.pdf"),
  width = 13,
  height = 3
)

draw(ht + ht_2)

dev.off()


# Write Session Info and Log Outputs --------------------------------------

cat(
  paste0("\n## Step completed\n",
         "Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
         "Finished 02_gsea_combined_healthy_controls\n"),
  file = "PROGRESS_LOG.md",
  append = TRUE
)

writeLines(capture.output(sessionInfo()), 
           file.path("session_info_logs",
                     "02_gsea_combined_healthy_controls_sessionInfo.txt"))

