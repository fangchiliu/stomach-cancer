# ==============================================================================
# 05_survival_analysis.R — Kaplan-Meier survival by INHBA expression
#
# Splits TCGA-STAD tumors at the median INHBA level and compares overall
# survival. Reports the full follow-up curve and a 5-year censored version.
# ==============================================================================

source("00_setup.R")
theme_set_project()

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(SummarizedExperiment)
  library(dplyr)
})

GENE      <- "INHBA"
FIVE_YEAR <- 1825   # days

vst_matrix    <- readRDS(exp_path("TCGA_vst_matrix.rds"))
gene_metadata <- readRDS(exp_path("gene_metadata.rds"))
rna           <- readRDS(exp_path("stad_rna.rds"))
clinical      <- readRDS(exp_path("stad_clinical.rds"))

# ---- Survival endpoints -------------------------------------------------------
# TCGA records time-to-event in two columns depending on vital status.
clinical_surv <- clinical %>%
  mutate(
    time   = ifelse(vital_status == "Alive", days_to_last_follow_up, days_to_death),
    status = ifelse(vital_status == "Alive", 0, 1)
  ) %>%
  filter(!is.na(time), time > 0) %>%
  distinct(submitter_id, .keep_all = TRUE)

# ---- Gene expression in tumors only --------------------------------------------
gene_id   <- gene_metadata$gene_id[gene_metadata$gene_name == GENE]
gene_expr <- vst_matrix[gene_id, ]

tumor_barcodes <- colData(rna)$barcode[colData(rna)$sample_type == "Primary Tumor"]
gene_expr      <- gene_expr[names(gene_expr) %in% tumor_barcodes]

# The first 12 characters of a TCGA barcode identify the patient, which is the
# key the clinical table uses.
expr_df <- data.frame(
  submitter_id = substr(names(gene_expr), 1, 12),
  expression   = as.numeric(gene_expr)
) %>%
  distinct(submitter_id, .keep_all = TRUE)

# ---- Median split --------------------------------------------------------------
surv_df <- expr_df %>%
  inner_join(clinical_surv, by = "submitter_id") %>%
  mutate(
    expr_group = ifelse(expression > median(expression), "High", "Low"),
    expr_group = factor(expr_group, levels = c("Low", "High"))
  )

print(table(surv_df$expr_group))

# ---- Full follow-up ------------------------------------------------------------
fit_full <- survfit(Surv(time, status) ~ expr_group, data = surv_df)

km_full <- ggsurvplot(
  fit_full, data = surv_df,
  pval = TRUE, risk.table = TRUE, conf.int = FALSE,
  palette = c(COL_NORMAL, COL_TUMOR),
  legend.title = paste(GENE, "expression"), legend.labs = c("Low", "High"),
  xlab = "Time (days)", ylab = "Overall survival probability",
  title = sprintf("Overall survival by %s expression (TCGA-STAD)", GENE)
)

png(plot_path(sprintf("KM_%s_high_vs_low.png", GENE)),
    width = 2400, height = 2400, res = 300)
print(km_full)
dev.off()

# ---- 5-year administratively censored ------------------------------------------
# Late follow-up in TCGA is sparse and noisy, so the headline figure caps at 5y.
surv_df_5yr <- surv_df %>%
  mutate(
    status = ifelse(time > FIVE_YEAR, 0, status),
    time   = pmin(time, FIVE_YEAR)
  )

fit_5yr <- survfit(Surv(time, status) ~ expr_group, data = surv_df_5yr)

km_5yr <- ggsurvplot(
  fit_5yr, data = surv_df_5yr,
  pval = TRUE, pval.size = 5, pval.coord = c(50, 0.05),
  risk.table = TRUE, risk.table.height = 0.20,
  risk.table.y.text = TRUE, risk.table.y.text.col = TRUE,
  risk.table.col = "strata", fontsize = 4,
  conf.int = FALSE, censor.shape = "|", censor.size = 3, size = 1.1,
  palette = c(COL_NORMAL, COL_TUMOR),
  legend = "top", legend.title = paste(GENE, "expression"),
  legend.labs = c("Low", "High"),
  xlab = "Time (days)", ylab = "Overall survival probability",
  break.time.by = 365, xlim = c(0, FIVE_YEAR),
  surv.median.line = "hv",
  ggtheme = theme_minimal(base_size = 14) +
    theme(
      plot.title       = element_text(face = "bold", size = 16),
      plot.subtitle    = element_text(color = "grey40", size = 11.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey93", linewidth = 0.3),
      legend.title     = element_text(size = 12, face = "bold")
    )
)

km_5yr$plot <- km_5yr$plot +
  labs(
    title    = sprintf("High %s expression predicts worse survival in gastric cancer", GENE),
    subtitle = "TCGA-STAD · 5-year overall survival · median split · log-rank test"
  )

km_5yr$table <- km_5yr$table +
  labs(title = "Number at risk", y = NULL, x = "Time (days)") +
  theme(plot.title = element_text(size = 12, face = "bold"),
        panel.grid = element_blank())

png(plot_path(sprintf("KM_%s_high_vs_low_5yr.png", GENE)),
    width = 2400, height = 3000, res = 300)
print(km_5yr)
dev.off()

# ---- Log-rank statistic ---------------------------------------------------------
print(survdiff(Surv(time, status) ~ expr_group, data = surv_df_5yr))

write.csv(surv_df %>% dplyr::select(submitter_id, expression, expr_group, time, status),
          exp_path(sprintf("%s_survival_table.csv", GENE)), row.names = FALSE)

cat("Survival analysis complete.\n")
