
# Load Packages -----------------------------------------------------------

library(tidyverse)
library(rstatix)
library(writexl)

dir.create(
  file.path(
    "output_data",
    "plots",
    "tage"
  )
)

cat(
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
         "Started 06_tage_analysis\n"),
  file = "PROGRESS_LOG.md",
  append = TRUE
)


data <- read.csv(
  file.path("input_data","annotated_tAge_meta_taco_2026-05-28.csv"),
  sep = ";")

colnames(data)

# Chronological Age -------------------------------------------------------


mouse_chrono_elastic_net <- data |> 
  dplyr::select(X.1,
                Individual,
                Treatment,
                Disease,
                batch,
                X,
                geo_sample,
                Disease_Treatment,
                tAge.chronological.multi_mouse.scaled.EN_Human_years)

mouse_chrono_elastic_net$Disease_Treatment

mouse_chrono_elastic_net$Disease_Treatment <- factor(
  mouse_chrono_elastic_net$Disease_Treatment,
  levels = c(
    "Healthy_DMSO",
    "HGPS_DMSO",
    "HGPS_UMK57",
    "DS_DMSO",
    "DS_UMK57"
  )
)

my_comparisons <- c(
  paste0("Healthy_DMSO", "HGPS_DMSO"), 
  paste0("Healthy_DMSO", "DS_DMSO"), 
  paste0("HGPS_DMSO", "HGPS_UMK57"),
  paste0("DS_DMSO", "DS_UMK57")
)


# Performed Dunn Test (like in Graphpad Kruskall-Wallis with Dunn's test for indiv.
# multiple comparisons), but only keeping the comparisons of interest
# and correcting for multiple hypothesis testing
# using the Benjamini-Hochberg method.


dunn_results <- mouse_chrono_elastic_net %>%
  dunn_test(
    tAge.chronological.multi_mouse.scaled.EN_Human_years ~ Disease_Treatment, 
    p.adjust.method = "BH") |> 
  dplyr::mutate(p.adj = signif(p.adj, 4)) |> 
  dplyr::mutate(comparisons = paste0(group1, group2))

dunn_results_true <- dunn_results |> 
  dplyr::filter(comparisons %in% my_comparisons)

dunn_results_true$p.adj <- p.adjust(dunn_results_true$p, method = "BH")

dunn_results_true$p.adj <- round(dunn_results_true$p.adj, 4)


cols <- c(
  "Healthy_DMSO" = "#616362",
  "HGPS_DMSO" = "#812723",
  "HGPS_UMK57" = "#bd9091",
  "DS_DMSO" = "#95aacd",
  "DS_UMK57" = "#ccd7e9"
)

set.seed(123) # Fixed Jitter

ggplot(
  mouse_chrono_elastic_net,
  aes(
    x = Disease_Treatment,
    y = tAge.chronological.multi_mouse.scaled.EN_Human_years
  )
) +
  geom_boxplot(outlier.shape = NA,
               show.legend = F, staplewidth = 0.5, width = 0.5,
               aes(
                 fill = Disease_Treatment
               )) +
  geom_jitter(
    shape = 21,
    size = 3,
    fill = "gray30",
    width = 0.15,
    alpha = 0.5
  ) +
  theme_classic(base_size = 15) +
  ggpubr::stat_pvalue_manual(
    dunn_results_true,
    y.position = c(11, 13, 8, 15),
    label = "p.adj"
  ) +
  scale_fill_manual(
    values = cols
  ) +
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  ) +
  labs(
    title = "Chrono. Age (Mouse Multi Tissue, Scaled, Elastic Net)",
    x = "\nCondition",
    y = "Chronological tAge (Scaled)"
  )

ggsave(
  "tage_mouse_multitissue_elastic_net.pdf",
  path = file.path("output_data", "plots", "tage"),
  width = 6,
  height = 6
)



# Mortality Clock ---------------------------------------------------------

colnames(data)

mouse_mortality_elastic_net <- data |> 
  dplyr::select(X.1,
                Individual,
                Treatment,
                Disease,
                batch,
                X,
                geo_sample,
                Disease_Treatment,
                tAge.mortality.multi_mouse.scaled.EN)

mouse_mortality_elastic_net$Disease_Treatment

mouse_mortality_elastic_net$Disease_Treatment <- factor(
  mouse_mortality_elastic_net$Disease_Treatment,
  levels = c(
    "Healthy_DMSO",
    "HGPS_DMSO",
    "HGPS_UMK57",
    "DS_DMSO",
    "DS_UMK57"
  )
)

my_comparisons <- c(
  paste0("Healthy_DMSO", "HGPS_DMSO"), 
  paste0("Healthy_DMSO", "DS_DMSO"), 
  paste0("HGPS_DMSO", "HGPS_UMK57"),
  paste0("DS_DMSO", "DS_UMK57")
)

dunn_results <- mouse_mortality_elastic_net %>%
  dunn_test(
    tAge.mortality.multi_mouse.scaled.EN ~ Disease_Treatment, 
    p.adjust.method = "BH") |> 
  dplyr::mutate(p.adj = signif(p.adj, 4)) |> 
  dplyr::mutate(comparisons = paste0(group1, group2))

dunn_results_true <- dunn_results |> 
  dplyr::filter(comparisons %in% my_comparisons)

dunn_results_true$p.adj <- p.adjust(dunn_results_true$p, method = "BH")

dunn_results_true$p.adj <- round(dunn_results_true$p.adj, 4)

dna_rep_cols <- c(
  "Healthy_DMSO" = "#616362",
  "HGPS_DMSO" = "#812723",
  "HGPS_UMK57" = "#bd9091",
  "DS_DMSO" = "#95aacd",
  "DS_UMK57" = "#ccd7e9"
)

cols <- c(
  "Healthy_DMSO" = "#616362",
  "HGPS_DMSO" = "#812723",
  "HGPS_UMK57" = "#bd9091",
  "DS_DMSO" = "#95aacd",
  "DS_UMK57" = "#ccd7e9"
)

set.seed(123) # Fixed Jitter

ggplot(
  mouse_mortality_elastic_net,
  aes(
    x = Disease_Treatment,
    y = tAge.mortality.multi_mouse.scaled.EN
  )
) +
  geom_boxplot(outlier.shape = NA,
               show.legend = F, staplewidth = 0.5, width = 0.5,
               aes(
                 fill = Disease_Treatment
               )) +
  geom_jitter(
    shape = 21,
    size = 3,
    width = 0.15,
    fill = "gray30",
    alpha = 0.5
  ) +
  theme_classic(base_size = 15) +
  ggpubr::stat_pvalue_manual(
    dunn_results_true,
    y.position = c(1, 1.5, 0.5, 1),
    label = "p.adj"
  ) +
  scale_fill_manual(
    values = cols
  ) +
  theme(
    plot.title = element_text(size = 15, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  ) +
  labs(
    title = "Mortality (Mouse Multi Tissue, Scaled, Elastic Net)",
    x = "\nCondition",
    y = "Mortality tAge (Log HR)"
  )

ggsave(
  "tage_mouse_mortality_multitissue_elastic_net.pdf",
  path = file.path("output_data", "plots", "tage"),
  width = 6,
  height = 6
)


# Source Data -------------------------------------------------------------

chrono_source_data <- mouse_chrono_elastic_net |> 
  dplyr::select(
    Disease_Treatment, 
    tAge.chronological.multi_mouse.scaled.EN_Human_years) |> 
  group_by(Disease_Treatment) |> 
  mutate(row = row_number()) |> 
  pivot_wider(names_from = Disease_Treatment, values_from = tAge.chronological.multi_mouse.scaled.EN_Human_years) |> 
  dplyr::select(-row) |> 
  dplyr::select(
    Healthy_DMSO,
    HGPS_DMSO,
    HGPS_UMK57,
    DS_DMSO,
    DS_UMK57
  )

mortality_source_data <- mouse_mortality_elastic_net |> 
  dplyr::select(
    Disease_Treatment, 
    tAge.mortality.multi_mouse.scaled.EN) |> 
  group_by(Disease_Treatment) |> 
  mutate(row = row_number()) |> 
  pivot_wider(names_from = Disease_Treatment, 
              values_from = tAge.mortality.multi_mouse.scaled.EN) |> 
  dplyr::select(-row) |> 
  dplyr::select(
    Healthy_DMSO,
    HGPS_DMSO,
    HGPS_UMK57,
    DS_DMSO,
    DS_UMK57
  )

dir.create(
  file.path(
    "output_data",
    "results",
    "tage"
  )
)



writexl::write_xlsx(
  chrono_source_data,
  file.path("output_data","results", "tage", "chrono_source_data.xlsx"))     


writexl::write_xlsx(
  mortality_source_data,
  file.path("output_data","results", "tage", "mortality_source_data.xlsx"))   



# Write Session Info and Log Outputs --------------------------------------


cat(
  paste0("\n## Step completed\n",
         "Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
         "Finished 06_tage_analysis\n"),
  file = "PROGRESS_LOG.md",
  append = TRUE
)


writeLines(capture.output(sessionInfo()), 
           file.path("session_info_logs","06_tage_analysis_sessionInfo.txt"))
