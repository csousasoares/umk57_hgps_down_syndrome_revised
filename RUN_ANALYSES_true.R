# run_all_analyses.R

# If needed, use: BiocManager::install(version = "3.21")

renv::restore(transactional = TRUE, prompt = FALSE)

source("01_umk57_vs_dmso_rna_seq_analysis.R")
source("02_gsea_combined_healthy_controls.R")
source("03_foxm1_kif2c_relation_analysis.R")
source("04_survival_analysis.R")
source("05_gene_length_analysis.R")
source("06_tage_analysis.R")
source("07_ssGSEA_and_reversion_magnitudes.R")
source("08_deseq2_alternative.R")

