# =============================================================================
# 03_dea_hot_excluded_vs_nonexcluded.R
# Differential expression analysis (limma) between Hot_Excluded and
# Hot_NonExcluded tumors.
# =============================================================================

source("R/00_setup.R")
LOG_FILE <- "logs/03_dea_hot_excluded_vs_nonexcluded.log"

log_msg("=== 03_dea_hot_excluded_vs_nonexcluded.R START ===")

data_in   <- readRDS("data/processed/tcga_phenotypes.rds")
master_df <- data_in$master_df
exp_tnbc  <- data_in$exp_tnbc


# --- Select samples ------------------------------------------

meta_limma <- master_df %>%
  filter(Immune_Group == "Hot") %>%
  dplyr::select(sample, Phenotype, Stromal_Group, TIS_Score, Stromal_Score)

log_msg("Total Hot samples: ", nrow(meta_limma))
log_msg("  Hot_Excluded:    ", sum(meta_limma$Stromal_Group == "Excluded"))
log_msg("  Hot_NonExcluded: ", sum(meta_limma$Stromal_Group == "NonExcluded"))


# --- Prepare expression matrix -------------------------------------------

# Exclude genes from both signatures + additional markers
exp_dea <- exp_tnbc[
  !(rownames(exp_tnbc) %in% GENES_EXCLUDED_DEA),
  meta_limma$sample
]

n_excluded_present <- sum(GENES_EXCLUDED_DEA %in% rownames(exp_tnbc))
log_msg("Signature genes excluded from the analysis: ", n_excluded_present)
log_msg("Genes remaining for DEA: ", nrow(exp_dea),
        " | Samples: ", ncol(exp_dea))

stopifnot(sum(rownames(exp_dea) %in% GENES_EXCLUDED_DEA) == 0)


# --- Differential expression ----------------------------------------------------------

group <- factor(
  meta_limma$Stromal_Group,
  levels = c("NonExcluded", "Excluded")
)
design   <- model.matrix(~ 0 + group)
colnames(design) <- c("NonExcluded", "Excluded")

fit      <- lmFit(exp_dea, design)
contrast <- makeContrasts(Excluded_vs_NonExcluded = Excluded - NonExcluded,
                          levels = design)
fit2     <- contrasts.fit(fit, contrast)
fit2     <- eBayes(fit2)

dea_results <- topTable(fit2, coef = "Excluded_vs_NonExcluded",
                        number = Inf, sort.by = "P") %>%
  rownames_to_column(var = "Gene") %>%
  filter(!is.na(adj.P.Val))

log_msg("Genes analyzed: ", nrow(dea_results))
log_msg("Significant (adj.P.Val < 0.05, |logFC| > 1): ",
        sum(dea_results$adj.P.Val < 0.05 & abs(dea_results$logFC) > 1))


# --- Export results --------------------------------------------------------

write.csv(dea_results,
          "results/tables/table02_dea_full.csv",
          row.names = FALSE)

log_msg("Table saved in results/tables/")


# --- Volcano plot ---------------------------------------------------------

dea_plot <- dea_results %>%
  mutate(
    Sig = adj.P.Val < 0.05 & abs(logFC) > 1,
    Category = case_when(
      Sig & logFC > 1             ~ "Up in Excluded",
      Sig & logFC < -1            ~ "Down in Excluded",
      TRUE                        ~ "Not significant"
    )
  )

label_genes <- dea_plot %>%
  filter(Sig & (adj.P.Val < 1e-5 | abs(logFC) > 2.5))

cat_colors <- c(
  "Up in Excluded"              = "#c0392b",
  "Down in Excluded"            = "#2980b9",
  "Not significant"             = "grey80"
)

p_volcano <- ggplot(dea_plot,
                    aes(x = logFC, y = -log10(adj.P.Val), color = Category)) +
  geom_point(alpha = 0.55, size = 1.5) +
  geom_text_repel(data = label_genes,
                  aes(label = Gene),
                  size = 3, max.overlaps = 25,
                  segment.color = "grey60") +
  scale_color_manual(values = cat_colors) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "grey60", linewidth = 0.4) +
  labs(
    title    = "DEA: Hot_Excluded vs Hot_NonExcluded (TCGA-TNBC)",
    subtitle = paste0("Hot tumors (n=", nrow(meta_limma), ")"),
    x        = "Log2 Fold Change (Excluded vs NonExcluded)",
    y        = "-Log10 adj. p-value (BH)",
    color    = NULL
  ) +
  theme(legend.position = "bottom")

ggsave("results/figures/fig03_volcano_hot_excluded.png",
       p_volcano, width = 9, height = 7, dpi = CFG$project$dpi)
log_msg("Figure saved: results/figures/fig03_volcano_hot_excluded.png")


# --- Save ------------------------------------------------------------------

saveRDS(
  list(
    dea_results = dea_results,
    meta_limma  = meta_limma,
    fit2        = fit2
  ),
  file = "data/processed/tcga_dea.rds"
)

log_msg("✓ Saved: data/processed/tcga_dea.rds")
log_msg("=== 03_dea_hot_excluded_vs_nonexcluded.R END ===")
