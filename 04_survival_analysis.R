## Load Packages ---------------------------------------------------------------

library(survival)
library(survminer)
library(tidyverse)
library(patchwork)


# Create Dir. -------------------------------------------------------------

dir.create(
  file.path(
    "output_data",
    "plots",
    "survival"
  ), recursive = T
)

dir.create(
  file.path(
    "output_data",
    "results",
    "survival"
  ), recursive = T
)

cat(
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
         "Started 04_survival_analysis\n"),
  file = "PROGRESS_LOG.md",
  append = TRUE
)

## LAKI Survival Data ----------------------------------------------------------

## Survival Analysis of LAKI mice treated with
## DMSO or UMK57, adjusted for sex.

## Load survival data:

df <- read.csv(
  file.path(
    "input_data",
    "survival",
    "survival_laki.csv"
    ), 
    header = TRUE, 
    sep = ";"
)

## Load information about sex of the animals:

sex_info <- read.csv(file.path(
  "input_data",
  "survival",
  "animal_sex_info.csv"), sep = ";")

length(unique(sex_info$Code))


## Ensures columns in df and sex_info2 are both named "Animal":

sex_info2 <- sex_info |> dplyr::distinct() |> dplyr::rename(Animal = Code)


## Merge survival and sex data:

df_sex <- merge(df, sex_info2, by.x = "Animal", by.y = "Animal")

nrow(df); nrow(df_sex)  # should match if no rows dropped


## Survival analysis inlcuding sex as covariate:

## coxph() fits a Cox regression model with Condition adjusted for Sez
## (DMSO vs UMK57) as the predictor. This model estimates 
## the hazard ratio (HR) — how much more (or less) likely a mouse 
## in one group is to die at any given moment compared to the other group. 
## A HR < 1 for UMK57 would suggest a protective effect.

cox_model <- coxph(
  Surv(time = Death, event = rep(1, nrow(df_sex))) ## 1 = Dead
  ~ Sex + Condition, data = df_sex) ## Adjusted for Sex

summary_results <- summary(cox_model) 

## Reveals p-value and effect size for each variable (sex or condition)

summary_results

## Confirm if proportional hazards change over time. 
## They shouldn't. Should be p > 0.05

cox.zph(cox_model) 

## No significant change over time, we can proceed.


## LogRank, Cox and Plot -------------------------------------------------------

## Get LogRank test results as comparison:

res <- survdiff(Surv(Death, rep(1, nrow(df_sex))) ~ Condition, data = df_sex)
res
pval_logrank <- 1 - pchisq(res$chisq, length(res$n) - 1)
pval_logrank


pval <- summary_results$coefficients
pval


pval_sex <- pval["SexM",5]
pval_condition <- pval["ConditionUMK57",5]

umk57_HR <- summary_results$coefficients["ConditionUMK57","exp(coef)"]

confidence_intervals <- summary_results[["conf.int"]]

umk57_HR_conf_lower_95 <- round(
  confidence_intervals["ConditionUMK57","lower .95"], 3)

umk57_HR_conf_upper_95 <- round(
  confidence_intervals["ConditionUMK57","upper .95"], 3)
  

## Create a plot, first survfit then ggsurvplot:

fit_long_rank <- survfit(Surv(Death, rep(1, nrow(df_sex))) ~ Condition, 
                         data = df_sex)

n_dmso <- nrow(df[df$Condition == "DMSO",])
n_umk57 <- nrow(df[df$Condition == "UMK57",])

p1 <- ggsurvplot(fit_long_rank, data = df_sex, pval = T, 
                 conf.int = F,
                 palette = c("#663e4e", "#d0c4ca"),
                 legend.labs = c(paste0("DMSO (n = ", n_dmso, ")"), 
                                 paste0("UMK57 (n = ", n_umk57, ")")))
p1

p1$plot + coord_cartesian(xlim = c(12, 20)) + 
  scale_x_continuous(breaks=seq(12, 20, 1)) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    labels = seq(0, 100, 25)
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_hline(aes(yintercept = 1, linetype = "WT (n = 13)"), 
             color = "gray50", linewidth = 0.8) +
  scale_linetype_manual(name = NULL, values = c("WT (n = 13)" = "solid")) + # Added Manually since all WT animals live well beyond 20w and were never part of the survival analysis to begin with.
  annotate(
    "text", 
    x = 12, 
    y = 0.25, 
    label = paste("Log-Rank Test: P = ", 
                  signif(pval_logrank, digits = 1)), hjust = 0) +
  annotate(
    "text", 
    x = 12, 
    y = 0.2, 
    label = paste("Cox Proportional Hazards (Condition): P = ", 
                  signif(pval_condition, digits = 2), "**"), hjust = 0) +
  annotate(
    "text", 
    x = 12, 
    y = 0.15, 
    label = paste("Cox Proportional Hazards (Sex): P = ", 
                  signif(pval_sex, digits = 2)), hjust = 0) +
  annotate(
    "text", 
    x = 12, 
    y = 0.10, 
    label = paste("Hazard Ratio UMK57 vs DMSO (HR) = ", 
                  signif(umk57_HR, digits = 2)), hjust = 0) +
  labs(y = "Survival Probability (%)", x = "Time (Weeks)") + 
  annotate(
    "text", 
    x = 12, 
    y = 0.05, 
    label = paste0("Hazard Ratio 95% CI = ", 
                   umk57_HR_conf_lower_95, " - ", umk57_HR_conf_upper_95),
    hjust = 0) +
  theme_bw(base_size = 14)

ggsave(
  "umk57_survival_cph_plot.pdf",
  path = file.path("output_data", "plots", "survival"),
  width = 12,
  height = 7,
  units = "in"
)


# Write Session Info and Log Outputs --------------------------------------


cat(
  paste0("\n## Step completed\n",
         "Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
         "Finished 04_survival_analysis\n"),
  file = "PROGRESS_LOG.md",
  append = TRUE
)


writeLines(capture.output(sessionInfo()), 
           file.path("session_info_logs","04_survival_analysis_sessionInfo.txt"))
