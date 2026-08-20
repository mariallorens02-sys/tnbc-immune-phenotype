# =============================================================================
# 02_define_axes_phenotypes.R
# Compute immune and stromal scores and classify tumors into four 
# microenvironment phenotypes.
# =============================================================================

source("R/00_setup.R")
LOG_FILE <- "logs/02_define_axes_phenotypes.log"

log_msg("=== 02_define_axes_phenotypes.R START ===")

data_in  <- readRDS("data/processed/tcga_tnbc_inputs.rds")
exp_tnbc <- data_in$exp_tnbc
clinical <- data_in$clinical


# --- Compute scores -----------------------------------------

log_msg("Computing TIS score (Ayers 2017, 18 genes)...")
tis_score     <- compute_signature_score(exp_tnbc, SIGS$tis,     "TIS")

log_msg("Computing Stromal score (TGF-β/CAF, 10 genes)...")
stromal_score <- compute_signature_score(exp_tnbc, SIGS$stromal, "Stromal")

master_df <- data.frame(
  sample        = colnames(exp_tnbc),
  TIS_Score     = tis_score,
  Stromal_Score = stromal_score
) %>%
  left_join(clinical, by = "sample")


# --- Correlation between axes --------------------------------------

log_msg("Checking orthogonality of the two axes (Spearman)...")
cor_test <- cor.test(master_df$TIS_Score, master_df$Stromal_Score,
                     method = "spearman")

log_msg("Spearman rho = ", round(cor_test$estimate, 4),
        " | p = ", signif(cor_test$p.value, 3))


# --- Phenotype classification --------------------------------------

tis_med     <- median(master_df$TIS_Score,     na.rm = TRUE)
stromal_med <- median(master_df$Stromal_Score, na.rm = TRUE)

master_df <- master_df %>%
  mutate(
    Immune_Group  = ifelse(TIS_Score     >= tis_med,     "Hot",      "Cold"),
    Stromal_Group = ifelse(Stromal_Score >= stromal_med, "Excluded", "NonExcluded"),
    Phenotype     = factor(
      paste(Immune_Group, Stromal_Group, sep = "_"),
      levels = c("Hot_Excluded", "Hot_NonExcluded",
                 "Cold_Excluded", "Cold_NonExcluded")
    )
  )

stopifnot(all(table(master_df$Phenotype) > 0))

log_msg("Medians — TIS: ", round(tis_med, 3),
        " | Stromal: ", round(stromal_med, 3))


# --- Phenotype distribution ------------------------------------------

pheno_table <- master_df %>%
  count(Immune_Group, Stromal_Group, Phenotype, name = "N") %>%
  mutate(Pct = round(100 * N / sum(N), 1))

log_msg("--- Phenotype distribution ---")
print(pheno_table)

write.csv(pheno_table,
          "results/tables/table01_phenotype_distribution.csv",
          row.names = FALSE)


# --- Heatmap ----------------

all_sig_genes <- c(SIGS$tis, SIGS$stromal)
present_genes <- intersect(all_sig_genes, rownames(exp_tnbc))

sample_order <- master_df %>%
  arrange(Phenotype, desc(TIS_Score)) %>%
  pull(sample)

mat_hm <- exp_tnbc[present_genes, sample_order]
mat_z  <- t(scale(t(mat_hm)))  # per-gene z-score

ann_col <- data.frame(
  Phenotype     = master_df$Phenotype[match(sample_order, master_df$sample)],
  TIS_Score     = master_df$TIS_Score[match(sample_order, master_df$sample)],
  Stromal_Score = master_df$Stromal_Score[match(sample_order, master_df$sample)],
  row.names = sample_order
)

ann_row <- data.frame(
  Signature = ifelse(present_genes %in% SIGS$tis, "TIS (immune)", "Stromal (exclusion)"),
  row.names = present_genes
)

ann_colors <- list(
  Phenotype = PHENO_COLORS,
  Signature = c("TIS (immune)" = "#2980b9", "Stromal (exclusion)" = "#e67e22")
)

png("results/figures/fig01a_axes_heatmap_confirmatory.png",
    width = 4200, height = 2400, res = CFG$project$dpi)
pheatmap(
  mat_z,
  annotation_col    = ann_col,
  annotation_row    = ann_row,
  annotation_colors = ann_colors,
  cluster_cols      = FALSE,          # samples already ordered by phenotype
  cluster_rows      = TRUE,
  clustering_method = "ward.D2",
  show_colnames     = FALSE,
  color  = colorRampPalette(c("#3498db", "white", "#e74c3c"))(60),
  breaks = seq(-3, 3, length.out = 61),
  main   = paste0("Immune and stromal signatures",
                  ncol(exp_tnbc), ")\n",
                  "Spearman rho(TIS, Stromal) = ",
                  round(cor_test$estimate, 3))
)
dev.off()
log_msg("Figure saved: results/figures/fig01_axes_heatmap.png")


# --- Scatter plot -----------------

p_scatter <- ggplot(master_df,
                    aes(x = TIS_Score, y = Stromal_Score, color = Phenotype)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_vline(xintercept = tis_med,     linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = stromal_med, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = PHENO_COLORS) +
  annotate("text",
           x = min(master_df$TIS_Score) + 0.05,
           y = max(master_df$Stromal_Score) - 0.05,
           label = paste0("Spearman ρ = ",
                          round(cor_test$estimate, 3),
                          "\np = ", signif(cor_test$p.value, 2)),
           hjust = 0, vjust = 1, size = 3.5, color = "grey30") +
  labs(
    title    = "Two microenvironment axes in TNBC (TCGA)",
    x        = "TIS Score",
    y        = "Stromal Score",
    color    = "Phenotype"
  ) +
  theme(legend.position = "right")

ggsave("results/figures/fig02_phenotype_scatter.png",
       p_scatter, width = 7, height = 5, dpi = CFG$project$dpi)
log_msg("Figure saved: results/figures/fig02_phenotype_scatter.png")


# --- Save -------------------------------------------------------------------

saveRDS(
  list(
    master_df = master_df,
    exp_tnbc  = exp_tnbc,
    cor_test  = cor_test,
    tis_med   = tis_med,
    stromal_med = stromal_med
  ),
  file = "data/processed/tcga_phenotypes.rds"
)

log_msg("✓ Saved: data/processed/tcga_phenotypes.rds")
log_msg("=== 02_define_axes_phenotypes.R END ===")
