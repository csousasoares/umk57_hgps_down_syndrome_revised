# Small-molecule targeting of kinesin-13 KIF2C enhances double-strand break DNA repair and slows down aging

Repository reproducing multiple analyses in the manuscript and rebuttal of the paper "Small-molecule 
targeting of kinesin-13 KIF2C enhances double-strand break DNA repair and 
slows down aging". This includes the analysis pipeline for the 
UMK57 vs. DMSO RNA-seq experiments in Hutchinson-Gilford Progeria Syndrome 
(HGPS) and Down Syndrome (DS) fibroblasts, plus downstream survival, 
transcriptomic clocks, and gene-length correlation analyses.

Analyses can be reproduced by running `RUN_ANALYSES_true.R`

## Installation

Tested in **R 4.5.1** and Bioconductor **3.21**.

Requires renv to install all dependencies automatically:

```r
install.packages("renv")
```

All package versions are pinned in `renv.lock`. To restore the
environment, simply run `RUN_ANALYSES_true.R` or: 

```r
renv::restore(transactional = TRUE, prompt = FALSE)
```

This installs the exact versions of every package used below, including:

- **CRAN / tidyverse**: tidyverse, patchwork, viridis, ggpointdensity, ggrepel,
  ggrastr, ggpubr, ggalt, rstatix, circlize, survival, survminer, writexl
- **Bioconductor**: DESeq2, GSVA, limma, ComplexHeatmap, ggtree,
  AnnotationDbi, org.Hs.eg.db, GenomicFeatures
- **GitHub/other**: msigdbr, clusterProfiler

If `renv::restore()` is not available or fails, install these manually via
`install.packages()` / `BiocManager::install()` before proceeding.

## Directory Structure

```
project_root/
├── input_data/
│   ├── counts_matrix.csv
│   ├── sample_info.csv
│   ├── senmayo.txt
│   ├── aging_133_tpm.csv
│   ├── sample_metadata_correct_tpm.csv
│   ├── annotated_tAge_meta_taco_2026-05-28.csv
│   └── survival/
│       ├── survival_laki.csv
│       └── animal_sex_info.csv
├── 01_umk57_vs_dmso_rna_seq_analysis.R
├── 02_gsea_combined_healthy_controls.R
├── 03_foxm1_kif2c_relation_analysis.R
├── 04_survival_analysis.R
├── 05_gene_length_analysis.R
├── 06_tage_analysis.R
├── 07_ssGSEA_and_reversion_magnitudes.R
├── 08_deseq2_alternative.R
└── RUN_ANALYSES_true.R
```

`output_data/`, `session_info_logs/`, and `PROGRESS_LOG.md` are created
automatically when the scripts run and do not need to exist beforehand.
Therefore, they can be deleted safely since they are completely regenerated
by `RUN_ANALYSES_true.R`.

## Required Run Order

Scripts must be run in numeric order because later scripts read the CSV/RDS
outputs of earlier ones. `RUN_ANALYSES_true.R` runs them in the correct order
automatically:

| Order | Script | Depends on output of |
|---|---|---|
| 1 | `01_umk57_vs_dmso_rna_seq_analysis.R` | — |
| 2 | `02_gsea_combined_healthy_controls.R` | 01 |
| 3 | `03_foxm1_kif2c_relation_analysis.R` | — |
| 4 | `04_survival_analysis.R` | — |
| 5 | `05_gene_length_analysis.R` | 01, 02 |
| 6 | `06_tage_analysis.R` | — |
| 7 | `07_ssGSEA_and_reversion_magnitudes.R` | 01, 02 |
| 8 | `08_deseq2_alternative.R` | — |


## Expected Inputs

| File | Used by | Description |
|---|---|---|
| `input_data/counts_matrix.csv` | 01, 02, 07, 08 | Raw gene counts, genes × samples |
| `input_data/sample_info.csv` | 01, 02, 07, 08 | Sample metadata |
| `input_data/senmayo.txt` | 02 | SenMayo senescence gene set |
| `input_data/aging_133_tpm.csv` | 03 | TPM expression matrix for age-correlation analysis |
| `input_data/sample_metadata_correct_tpm.csv` | 03 | Sample/age metadata matching `aging_133_tpm.csv` columns |
| `input_data/annotated_tAge_meta_taco_2026-05-28.csv` | 06 | Transcriptomic-age (tAge) predictions per sample |
| `input_data/survival/survival_laki.csv` | 04 | LAKI mouse survival/death data |
| `input_data/survival/animal_sex_info.csv` | 04 | Animal sex metadata |

## Expected Outputs

All outputs are written under `output_data/`, organized by analysis:

```
output_data/
├── results/
│   ├── deseq2/                 # DESeq2 tables 
│   ├── gsea/                   # GSEA/Hallmark tables
│   ├── age_correlations/       # Spearman correlation tables for aging HDFs
│   ├── survival/               # Cox model summaries for survival analysis
│   ├── tage/                   # Dunn test tables for tAge clocks
│   ├── ssGSEA_and_reversion_magnitude/  (empty for the time being)
│   ├── deseq2_alternative/     # Batch-adjusted/paired DESeq2 tables
│   └── gsea_alternative/       # GSEA tables from DESeq2 alternative designs
└── plots/
    ├── age_correlation/        # Correlations using aging HDFs data
    ├── gene_length_analysis/   # Correlations of DGE with gene length
    ├── gsea/                   # GSEA analysis plots in the manuscript
    ├── gsea_alternative/       # GSEA analysis using alternative designs for rebuttal
    ├── heatmaps/               # Heatmap plots in the manuscript
    ├── pca/                    # PCA plots
    ├── ssGSEA_and_reversion_magnitude/    # ssGSEA heatmaps for rebuttal
    ├── survival/               # Survival analysis plots
    ├── tage/                   # Transcriptomic clock analysis plots
    ├── volcano_plots/          # Volcano plots


Additionally:

- `PROGRESS_LOG.md` — timestamped log of each script's start/completion,
  created fresh by script 01 and appended to by every subsequent script.
- `session_info_logs/*_sessionInfo.txt` — `sessionInfo()` output per script,
  for reproducibility records.

## Regenerating All Analyses (Single Command)

From a clean R session, with the working directory set to the project root:

```bash
Rscript RUN_ANALYSES_true.R
```

This will:
1. Restore the pinned package environment via `renv::restore()`
2. Run scripts 01–08 in order (dependencies satisfied automatically)

Or simply run the RUN_ANALYSES_true.R in `RStudio`

Rscript may need to be added to PATH before execution.

## Notes

- `set.seed()` is called on multiple scripts for reproducibility of
  stochastic steps (GSEA permutations, k-means clustering in heatmaps).
  Re-running an individual script in isolation reproduces its own results;
  re-running scripts out of order or skipping earlier scripts will cause
  downstream scripts to fail (missing input files) or silently use stale
  outputs from a previous run.
- Scripts assume a universal working directory layout (`file.path()` is used
  throughout).
