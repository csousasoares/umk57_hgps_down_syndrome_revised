## Load Packages ---------------------------------------------------------------

set.seed(123)

library(tidyverse)
library(org.Hs.eg.db)
library(patchwork)

dir.create(
  file.path(
    "output_data",
    "plots",
    "age_correlation"
  )
)

cat(
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
         "Started 03_foxm1_kif2c_relation_analysis\n"),
  file = "PROGRESS_LOG.md",
  append = TRUE
)


# Get TPM From Aging HDFs -------------------------------------------------

aging_133_tpm <- read.csv(
  file.path("input_data", "aging_133_tpm.csv"),
  header = T, 
  sep = ";", 
  row.names = 1)

metadata <- read.csv(
  file.path("input_data", "sample_metadata_correct_tpm.csv"), 
   header = T, sep = ";")

metadata$Sample == colnames(aging_133_tpm)


# Log2(TPM + 1) Tranformation ---------------------------------------------

gene_tpm_log2 <- log2(aging_133_tpm + 1)
ages <- metadata$Age

correlations_log2 <- numeric(nrow(gene_tpm_log2))
pval_log2 <- numeric(nrow(gene_tpm_log2))


# Create loop to obtain Log2(TPM+1) correlation with age:

for (i in 1:nrow(gene_tpm_log2)) {
  counts_log <- as.numeric(gene_tpm_log2[i, ])
  test_log <- cor.test(counts_log, ages, method = "spearman")
  correlations_log2[i] <- test_log$estimate
  pval_log2[i] <- test_log$p.value
} # This may take a while...

result_df_log2 <- data.frame(
  Gene = row.names(gene_tpm_log2),
  Correlation = correlations_log2,
  P_value = pval_log2
)

result_df_log2 <- result_df_log2 %>% na.omit()

genes_2 <- result_df_log2[, "Gene"]

annots <- AnnotationDbi::select(
  org.Hs.eg.db, 
  keys = genes_2, 
  columns = "SYMBOL", 
  keytype = "ENTREZID")

result_df_log2_2 <- merge(result_df_log2, 
                          annots, 
                          by.x = "Gene", 
                          by.y = "ENTREZID")


# Now correct for FDR:

result_df_log2_2$fdr <- p.adjust(
  result_df_log2_2$P_value, 
  method = "BH")


dir.create(
  file.path("output_data", "results", "age_correlations")
)

write.csv(result_df_log2_2, 
          file.path("output_data", 
                    "results", 
                    "age_correlations",
                    "aging_log2tpm_spearman_corr.csv")
          )


# Get FOXM1 and KIF2C data:

goi <- c("FOXM1", "KIF2C")

annots <- AnnotationDbi::select(
  org.Hs.eg.db, 
  keys = goi, 
  columns = "ENTREZID", 
  keytype = "SYMBOL")

gene_tpm_log2_goi <- gene_tpm_log2[annots$ENTREZID,]

gene_tpm_log2_goi <- t(gene_tpm_log2_goi) %>% as.data.frame()

gene_tpm_log2_goi$sample <- row.names(gene_tpm_log2_goi)


foxm1_kif2c_source <- gene_tpm_log2_goi |> 
  dplyr::rename(Sample = sample) |> 
  dplyr::left_join(metadata, join_by("Sample"))

write.csv(
  foxm1_kif2c_source,
  file.path("output_data", 
            "results", 
            "age_correlations",
            "foxm1_kif2c_data_points_spearman_corr_with_aging.csv")
)

gene_tpm_log2_goi_2 <- merge(gene_tpm_log2_goi, 
                             metadata, 
                             by.x = "sample", 
                             by.y = "Sample")

colnames(gene_tpm_log2_goi_2)

gene_tpm_log2_goi_2 <- gene_tpm_log2_goi_2 |>  
  dplyr::rename("FOXM1" = "2305", "KIF2C" = "11004")


gene_tpm_log2_goi_2 <- gene_tpm_log2_goi_2 |> 
  pivot_longer(names_to = "Gene", 
               values_to = "log2TPM", 
               cols = c(FOXM1, KIF2C))

cols = c(
  "FOXM1" = "gray20",
  "KIF2C" = "darkorange"
)

result_df_log2_2 |> 
  dplyr::filter(SYMBOL %in% c("KIF2C", "FOXM1")) |> 
  head()

rho_aging_foxm1 <- round(result_df_log2_2 |> 
  dplyr::filter(SYMBOL == "FOXM1") |> 
  dplyr::pull(Correlation), 2)
rho_aging_foxm1

rho_aging_kif2c <- round(
  result_df_log2_2 |> 
  dplyr::filter(SYMBOL == "KIF2C") |> 
  dplyr::pull(Correlation), 2)
rho_aging_kif2c

padj_aging_foxm1 <- signif(result_df_log2_2 |> 
  dplyr::filter(SYMBOL == "FOXM1") |> 
  dplyr::pull(fdr), 3)
padj_aging_foxm1


padj_aging_kif2c <- result_df_log2_2 |> 
  dplyr::filter(SYMBOL == "KIF2C") |> 
  dplyr::pull(fdr) |> 
  formatC(format = "E", digits = 2) # Keeps 1.40 instead of 1.4
padj_aging_kif2c


p1 <- ggplot(
  gene_tpm_log2_goi_2, 
    aes(
      x = Age, 
      y = log2TPM, 
      colour = Gene)) + 
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(se = T, method = "lm") +
  xlim(0, 101) +
  scale_color_manual(values = cols) +
  theme_bw(base_size = 16) +
  labs(x = "Donor Age", y = "Log2(TPM + 1)",
       title = "Gene Expression vs Donor Age") +
  annotate("text", 
           label = paste0(
             "Rho = ", rho_aging_foxm1 ," (padjust = ", padj_aging_foxm1, ")"
             ), 
           x = 0.2, y = 2.2, hjust = 0, color = "gray20", size = 4) +
  annotate("text", 
           label = paste0(
             "Rho = ", rho_aging_kif2c ," (padjust = ", padj_aging_kif2c, ")"
           ), 
           x = 0.2, y = 1.7, hjust = 0, color = "darkorange", size = 4) 
p1

# Correlation between KIF2C and FOXM1 themselves:

gene_tpm_log2_goi <- gene_tpm_log2[annots$ENTREZID,]

gene_tpm_log2_goi <- t(gene_tpm_log2_goi) %>% as.data.frame()

gene_tpm_log2_goi$sample <- row.names(gene_tpm_log2_goi)

gene_tpm_log2_goi_2 <- merge(gene_tpm_log2_goi, 
                             metadata, 
                             by.x = "sample", 
                             by.y = "Sample")

colnames(gene_tpm_log2_goi_2)

gene_tpm_log2_goi_2 <- gene_tpm_log2_goi_2 |> 
  dplyr::rename("FOXM1" = "2305", "KIF2C" = "11004")



min(gene_tpm_log2_goi_2$FOXM1)
min(gene_tpm_log2_goi_2$KIF2C)



test_log <- cor.test(gene_tpm_log2_goi_2$FOXM1, 
                     gene_tpm_log2_goi_2$KIF2C, 
                     method = "spearman")

correl_kif_fox <- signif(test_log$estimate, 2)

pval_correl_kif_fox <- signif(test_log$p.value, 3)


p2 <- ggplot(
  gene_tpm_log2_goi_2, 
  aes(x = FOXM1, y = KIF2C)) + 
  geom_point(aes(color = Age), size = 2, alpha = 0.7) +
  scale_color_gradient2(low = "green", 
                        mid = "yellow",
                        high = "magenta", 
                        limits = c(0, 100), 
                        midpoint = 50) +
  geom_smooth(method = "lm") +
  xlim(2,8) +
  theme_bw(base_size = 16) +
  labs(x = "FOXM1 Log2(TPM + 1)", 
       y = "KIF2C Log2(TPM + 1)",
       title = "KIF2C vs FOXM1 Gene Expression") +
  annotate("text", 
           label = paste0(
             "Rho = ", 
             correl_kif_fox, 
             " (P = ",
             pval_correl_kif_fox,
             ")"
             ), 
           x = 4.2, 
           y = 1.6, 
           hjust = 0, 
           color = "gray20", 
           size = 4)
p2

p3 <- p1 + p2
p3

ggsave(
  "KIF2C_vs_FOXM1_age_plot.pdf",
  path = file.path("output_data", "plots", "age_correlation"),
  width = 13,
  height = 6
)

# Write Session Info and Log Outputs --------------------------------------


cat(
  paste0("\n## Step completed\n",
         "Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
         "Finished 03_foxm1_kif2c_relation_analysis\n"),
  file = "PROGRESS_LOG.md",
  append = TRUE
)


writeLines(capture.output(sessionInfo()), 
           file.path("session_info_logs",
           "03_foxm1_kif2c_relation_analysis_sessionInfo.txt"))
