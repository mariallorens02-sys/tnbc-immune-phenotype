# =============================================================================
# 05_metabric_validation.R
#
# Objectives
# 1. Reproduce the four phenotypes
# 2. Validate the TCGA stromal program
# 3. Perform an independent DEA in METABRIC and assess functional 
# concordance using Hallmark GSEA.
# 4. Common DEGs exploration.
# =============================================================================

source("R/00_setup.R")
LOG_FILE <- "logs/05_metabric_validation.log"

log_msg("=== 05_metabric_validation.R START ===")

tcga_dea   <- readRDS("data/processed/tcga_dea.rds")


# --- Load METABRIC ---------------------------------------------------------

log_msg("Loading METABRIC expression...")
mb_raw <- read.table(PATHS$metabric_expr,
                     header = TRUE, sep = "\t", check.names = FALSE)

# Resolve duplicates: keep the gene with the highest mean expression
row_means  <- rowMeans(mb_raw[, -c(1, 2)], na.rm = TRUE)
mb_ord     <- mb_raw[order(row_means, decreasing = TRUE), ]
mb_uniq    <- mb_ord[!duplicated(mb_ord$Hugo_Symbol), ]
rownames(mb_uniq) <- mb_uniq$Hugo_Symbol
mb_matrix  <- as.matrix(mb_uniq[, -c(1, 2)])
storage.mode(mb_matrix) <- "numeric"

log_msg("Unique genes: ", nrow(mb_matrix), " | Samples: ", ncol(mb_matrix))


# --- Filter TNBC -------------------------------------------------------------

log_msg("Loading METABRIC clinical data and filtering TNBC...")
mb_clin <- read.table(PATHS$metabric_clin,
                      header = TRUE, sep = "\t",
                      comment.char = "#", check.names = FALSE)
rownames(mb_clin) <- mb_clin[, 3]  # sample ID column

mb_tnbc_meta <- mb_clin[
  mb_clin$ER_STATUS  == "Negative" &
    mb_clin$PR_STATUS  == "Negative" &
    mb_clin$HER2_STATUS == "Negative", ]

log_msg("TNBC patients in METABRIC: ", nrow(mb_tnbc_meta))

common_mb   <- intersect(rownames(mb_tnbc_meta), colnames(mb_matrix))
mb_tnbc_exp <- mb_matrix[, common_mb]
log_msg("Samples with expression + clinical data: ", length(common_mb))


# --- Compute scores -----------------------------------------------------------

log_msg("Computing TIS score in METABRIC...")
tis_mb     <- compute_signature_score(mb_tnbc_exp, SIGS$tis,     "TIS_METABRIC")

log_msg("Computing Stromal score in METABRIC...")
stromal_mb <- compute_signature_score(mb_tnbc_exp, SIGS$stromal, "Stromal_METABRIC")

mb_df <- data.frame(
  sample        = common_mb,
  TIS_Score     = tis_mb,
  Stromal_Score = stromal_mb
)


# --- Classify phenotypes (METABRIC's own medians) ----------------------------

tis_med_mb     <- median(mb_df$TIS_Score,     na.rm = TRUE)
stromal_med_mb <- median(mb_df$Stromal_Score, na.rm = TRUE)

mb_df <- mb_df %>%
  mutate(
    Immune_Group  = ifelse(TIS_Score     >= tis_med_mb,     "Hot",      "Cold"),
    Stromal_Group = ifelse(Stromal_Score >= stromal_med_mb, "Excluded", "NonExcluded"),
    Phenotype     = factor(
      paste(Immune_Group, Stromal_Group, sep = "_"),
      levels = c("Hot_Excluded", "Hot_NonExcluded",
                 "Cold_Excluded", "Cold_NonExcluded")
    )
  )


# --- Phenotype distribution --------------------------------

mb_pheno_table <- mb_df %>%
  count(Phenotype, name = "N") %>%
  mutate(Pct = round(100 * N / sum(N), 1))

log_msg("--- Phenotype distribution in METABRIC ---")
print(mb_pheno_table)

write.csv(mb_pheno_table,
          "results/tables/table05_metabric_distribution.csv",
          row.names = FALSE)


# --- Figure: METABRIC phenotype scatter --------------------------------------

cor_mb <- cor.test(mb_df$TIS_Score, mb_df$Stromal_Score, method = "spearman")

p_mb_scatter <- ggplot(mb_df,
                       aes(x = TIS_Score, y = Stromal_Score, color = Phenotype)) +
  geom_point(alpha = 0.65, size = 2) +
  geom_vline(xintercept = tis_med_mb,     linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = stromal_med_mb, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = PHENO_COLORS) +
  annotate("text",
           x = min(mb_df$TIS_Score) + 0.05,
           y = max(mb_df$Stromal_Score) - 0.05,
           label = paste0("Spearman ρ = ", round(cor_mb$estimate, 3),
                          "\np = ", signif(cor_mb$p.value, 2)),
           hjust = 0, vjust = 1, size = 3.5, color = "grey30") +
  labs(
    title    = "Phenotypes in METABRIC-TNBC",
    subtitle = paste0("n=", nrow(mb_df), " | Illumina microarray"),
    x        = "TIS Score",
    y        = "Stromal Score",
    color    = "Phenotype"
  )

ggsave("results/figures/fig05_metabric_phenotypes.png",
       p_mb_scatter, width = 7, height = 5, dpi = CFG$project$dpi)
log_msg("Figure saved: results/figures/fig05_metabric_phenotypes.png")


# --- Program validation -----------------------

program_genes_tcga <- tcga_dea$dea_results %>%
  filter(adj.P.Val < 0.05, logFC > 1) %>%
  filter(!Gene %in% GENES_EXCLUDED_DEA) %>%
  pull(Gene)

log_msg("UP program genes in TCGA (adj.P<0.05, logFC>1, non-signature): ",
        length(program_genes_tcga))

program_in_mb <- intersect(program_genes_tcga, rownames(mb_tnbc_exp))
log_msg("Program genes present in METABRIC: ",
        length(program_in_mb), "/", length(program_genes_tcga))

if (length(program_in_mb) < 5) {
  log_msg("WARNING: <5 program genes available in METABRIC — ",
          "score may be unreliable", level = "WARN")
}


mb_hot_samples <- mb_df$sample[mb_df$Immune_Group == "Hot"]
program_score_mb <- compute_signature_score(
  mb_tnbc_exp[, mb_hot_samples, drop = FALSE],
  program_in_mb,
  "DEA_Program_METABRIC"
)

mb_hot_prog <- mb_df %>%
  filter(Immune_Group == "Hot") %>%
  mutate(Program_Score = program_score_mb[match(sample, names(program_score_mb))])

wilcox_prog <- wilcox.test(Program_Score ~ Stromal_Group,
                           data   = mb_hot_prog,
                           exact  = FALSE,
                           conf.int = TRUE)

log_msg("Wilcoxon Program_Score ~ Stromal_Group (METABRIC Hot):")
log_msg("  W = ", wilcox_prog$statistic,
        " | p = ", signif(wilcox_prog$p.value, 3))

if (!is.null(wilcox_prog$estimate)) {
  log_msg("  Location shift = ",
          round(wilcox_prog$estimate, 3))
}

prog_summary <- mb_hot_prog %>%
  group_by(Stromal_Group) %>%
  summarise(
    n      = n(),
    median = round(median(Program_Score, na.rm = TRUE), 3),
    IQR    = round(IQR(Program_Score,    na.rm = TRUE), 3),
    .groups = "drop"
  )
cat("\n--- Program Score by group (METABRIC Hot) ---\n")
print(prog_summary)


# --- Program score boxplot ---------------------------------------------------------

p_label <- if (wilcox_prog$p.value < 0.001) "p < 0.001" else
  paste0("p = ", round(wilcox_prog$p.value, 3))

p_program <- ggplot(mb_hot_prog,
                    aes(x = Stromal_Group, y = Program_Score,
                        fill = Stromal_Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.55, color = "grey30") +
  annotate("segment",
           x = 1, xend = 2,
           y = max(mb_hot_prog$Program_Score, na.rm = TRUE) * 1.05,
           yend = max(mb_hot_prog$Program_Score, na.rm = TRUE) * 1.05,
           color = "grey30") +
  annotate("text",
           x = 1.5,
           y = max(mb_hot_prog$Program_Score, na.rm = TRUE) * 1.10,
           label = p_label, size = 3.8, color = "grey20") +
  scale_fill_manual(
    values = c("Excluded"    = CFG$colors$hot_excluded,
               "NonExcluded" = CFG$colors$hot_non_excluded)
  ) +
  scale_x_discrete(labels = c("Excluded"    = "Hot_Excluded",
                              "NonExcluded" = "Hot_NonExcluded")) +
  labs(
    title    = "Program validation: TCGA DEA-UP score in METABRIC",
    subtitle = paste0("UP genes from the TCGA DEA (n=", length(program_in_mb),
                      " present) | Hot samples only | Wilcoxon"),
    x        = NULL,
    y        = "Program Score (mean z-score of TCGA-DEA UP genes)",
    fill     = NULL
  ) +
  theme(legend.position = "none")

ggsave("results/figures/fig07_metabric_program_score.png",
       p_program, width = 5, height = 5, dpi = CFG$project$dpi)
log_msg("Figure saved: results/figures/fig07_metabric_program_score.png")


# ---  GSEA (functional replication) ---------------------------------------------

log_msg("--- GSEA in METABRIC (functional replication) ---")

# --- (a) DEA Hot_Excluded vs Hot_NonExcluded in METABRIC ---------------------
# Same approach as the TCGA DEA (script 03), applied to METABRIC.

mb_hot_df <- mb_df %>% filter(Immune_Group == "Hot")
log_msg("METABRIC Hot: n = ", nrow(mb_hot_df),
        " | Excluded = ", sum(mb_hot_df$Stromal_Group == "Excluded"),
        " | NonExcluded = ", sum(mb_hot_df$Stromal_Group == "NonExcluded"))

# Expression matrix of Hot samples, excluding signature genes (circularity)
expr_mb_hot <- mb_tnbc_exp[, mb_hot_df$sample]
expr_mb_hot <- expr_mb_hot[!rownames(expr_mb_hot) %in% GENES_EXCLUDED_DEA, ]

grp_mb <- factor(mb_hot_df$Stromal_Group,
                 levels = c("NonExcluded", "Excluded"))
design_mb <- model.matrix(~ 0 + grp_mb)
colnames(design_mb) <- c("NonExcluded", "Excluded")

fit_mb  <- lmFit(expr_mb_hot, design_mb)
cont_mb <- makeContrasts(Excluded - NonExcluded, levels = design_mb)
fit2_mb <- eBayes(contrasts.fit(fit_mb, cont_mb))
dea_mb  <- topTable(fit2_mb, number = Inf, sort.by = "none")

log_msg("METABRIC DEA complete: ", nrow(dea_mb), " genes")

# Full METABRIC DEA table (Gene column added from rownames)
dea_mb_full <- data.frame(Gene = rownames(dea_mb), dea_mb, row.names = NULL)
write.csv(dea_mb_full, "results/tables/table_metabric_dea_full.csv",
          row.names = FALSE)


# --- (b) Hallmark GSEA on the moderated-t ranking -----------------------------

hallmark <- tryCatch({
  hm_df <- tryCatch(
    msigdbr(species = "Homo sapiens", collection = "H"),
    error = function(e) msigdbr(species = "Homo sapiens", category = "H")
  )
  split(hm_df$gene_symbol, hm_df$gs_name)
}, error = function(e) {
  log_msg("Could not load Hallmark via msigdbr: ", e$message, level = "WARN")
  NULL
})

if (!is.null(hallmark)) {
  ranks_mb <- dea_mb$t
  names(ranks_mb) <- rownames(dea_mb)
  ranks_mb <- sort(ranks_mb[is.finite(ranks_mb)], decreasing = TRUE)

  gsea_mb <- suppressWarnings(
    fgsea(hallmark, ranks_mb, minSize = 10, maxSize = 500)
  )
  gsea_mb <- gsea_mb[order(gsea_mb$NES, decreasing = TRUE), ]

  # Save full METABRIC GSEA
  gsea_mb_out <- as.data.frame(gsea_mb[, c("pathway","NES","pval","padj","size")])
  write.csv(gsea_mb_out,
            "results/tables/table10_metabric_gsea.csv",
            row.names = FALSE)
  log_msg("METABRIC GSEA saved: results/tables/table10_metabric_gsea.csv")
}

  # --- (c) Compare key pathways TCGA vs METABRIC -----------------------------

  key_pathways <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                    "HALLMARK_TGF_BETA_SIGNALING",
                    "HALLMARK_ANGIOGENESIS",
                    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                    "HALLMARK_E2F_TARGETS",
                    "HALLMARK_MYC_TARGETS_V1",
                    "HALLMARK_G2M_CHECKPOINT")

  # TCGA NES: taken from the TCGA GSEA results table on disk.
  tcga_gsea <- tryCatch(
    read.csv("results/tables/table04_gsea_results.csv", stringsAsFactors = FALSE),
    error = function(e) NULL
  )

  comp_rows <- lapply(key_pathways, function(pw) {
    nes_mb  <- gsea_mb$NES[gsea_mb$pathway == pw]
    padj_mb <- gsea_mb$padj[gsea_mb$pathway == pw]
    nes_tcga <- if (!is.null(tcga_gsea)) {
      tcga_gsea$NES[tcga_gsea$pathway == pw]
    } else NA
    data.frame(
      pathway   = sub("HALLMARK_", "", pw),
      NES_TCGA  = ifelse(length(nes_tcga) == 0, NA, round(nes_tcga, 2)),
      NES_METABRIC = ifelse(length(nes_mb) == 0, NA, round(nes_mb, 2)),
      padj_METABRIC = ifelse(length(padj_mb) == 0, NA, signif(padj_mb, 3)),
      Concordant  = ifelse(length(nes_mb) == 0 || length(nes_tcga) == 0, NA,
                            sign(nes_mb) == sign(nes_tcga))
    )
  })
  comp_tab <- do.call(rbind, comp_rows)

  write.csv(comp_tab,
            "results/tables/table11_gsea_tcga_vs_metabric.csv",
            row.names = FALSE)

  log_msg("--- Key pathway comparison TCGA vs METABRIC ---")
  for (i in seq_len(nrow(comp_tab))) {
    log_msg("  ", comp_tab$pathway[i],
            ": NES_TCGA = ", comp_tab$NES_TCGA[i],
            " | NES_METABRIC = ", comp_tab$NES_METABRIC[i],
            " | padj_MB = ", comp_tab$padj_METABRIC[i],
            " | concordant = ", comp_tab$Concordant[i])
  }
  n_conc <- sum(comp_tab$Concordant, na.rm = TRUE)
  log_msg("Key pathways concordant in direction: ", n_conc, "/",
          sum(!is.na(comp_tab$Concordant)))


  # --- (d) Figure: NES TCGA vs METABRIC on key pathways ------------------------

  comp_plot <- comp_tab[!is.na(comp_tab$NES_TCGA) &
                          !is.na(comp_tab$NES_METABRIC), ]

  if (nrow(comp_plot) > 0) {
    p_gsea_comp <- ggplot(comp_plot,
                        aes(x = NES_TCGA, y = NES_METABRIC, label = pathway)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey40") +
    geom_point(aes(color = NES_TCGA > 0), size = 3, alpha = 0.8) +
    geom_text_repel(size = 3, max.overlaps = 20) +
    scale_color_manual(values = c("TRUE" = "#c0392b", "FALSE" = "#27ae60"),
                       guide = "none") +
    labs(
      title    = "Functional replication: NES in TCGA vs METABRIC",
      subtitle = "Key pathways from the Hot_Excluded vs Hot_NonExcluded contrast | diagonal = concordance",
      x = "NES in TCGA (discovery)",
      y = "NES in METABRIC (validation)"
    )
    
    ggsave("results/figures/fig12_gsea_tcga_vs_metabric.png",
         p_gsea_comp, width = 7, height = 6, dpi = CFG$project$dpi)
    
    log_msg("Figure saved: results/figures/fig12_gsea_tcga_vs_metabric.png")
    
  } else {
    
    log_msg("No common pathways available for comparison.",
            level="WARN")
}

  
# --- Common DEGs exploration ----------------------------------------------------
log_msg("Common DEGs exploration...")
log_msg("Criteria: TCGA (adj.P<0.05 & |logFC|>1) ∩ METABRIC (adj.P<0.05) ",
          "with concordance in direction")

dea_tcga <- read.csv("results/tables/table02_dea_full.csv",
                     stringsAsFactors = FALSE)
write.csv(data.frame(Gene = rownames(dea_mb), dea_mb),
          'results/tables/table_metabric_dea_full.csv', row.names = FALSE)
dea_mb_path <- "results/tables/table_metabric_dea_full.csv"
dea_mb <- read.csv(dea_mb_path, stringsAsFactors = FALSE)

log_msg("DEA TCGA: ", nrow(dea_tcga), " genes | DEA METABRIC: ",
        nrow(dea_mb), " genes")

tcga_sig <- dea_tcga %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1) %>%
  select(Gene, logFC_TCGA = logFC, adjP_TCGA = adj.P.Val)

mb_sig <- dea_mb %>%
  filter(adj.P.Val < 0.05) %>%
  select(Gene, logFC_METABRIC = logFC, adjP_METABRIC = adj.P.Val)

log_msg("TCGA significant genes (adj.P<0.05, |logFC|>1): ", nrow(tcga_sig))
log_msg("METABRIC significant genes (adj.P<0.05): ", nrow(mb_sig))

common <- inner_join(tcga_sig, mb_sig, by = "Gene") %>%
  mutate(
    Concordant = sign(logFC_TCGA) == sign(logFC_METABRIC),
  ) %>%
  filter(Concordant) %>%
  arrange(desc(abs(logFC_TCGA)))

log_msg("Common DEGs (sig. in both and concordant in sign): ", nrow(common))
common_up   <- common %>% filter(logFC_TCGA > 0)
common_down <- common %>% filter(logFC_TCGA < 0)

log_msg("Common UP in Excluded: ", nrow(common_up))
log_msg("Common DOWN in Excluded: ", nrow(common_down))

write.csv(common_up,   "results/tables/table12_common_degs_up.csv",
          row.names = FALSE)
write.csv(common_down, "results/tables/table12b_common_degs_down.csv",
          row.names = FALSE)

dea_tcga_plot <- dea_tcga %>%
  mutate(
    Common    = Gene %in% common$Gene,
    Sig       = adj.P.Val < 0.05 & abs(logFC) > 1
  )

p_volcano_common <- ggplot(dea_tcga_plot,
                           aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(data = ~ filter(.x, !Sig), color = "grey80", size = 1, alpha = 0.5) +
  geom_point(data = ~ filter(.x, Sig & !Common), color = "grey50", size = 1.2, alpha = 0.6) +
  geom_point(data = ~ filter(.x, Common), color = "#2980b9", size = 1.8, alpha = 0.8) +
  geom_text_repel(data = ~ filter(.x, Common),
                  aes(label = Gene), size = 3, max.overlaps = 20,
                  color = "#c0392b") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60") +
  labs(
    title    = "Common DEGs in TCGA and METABRIC",
    x = "log2 Fold Change (Excluded vs NonExcluded, TCGA)",
    y = "-log10(adj. p-valor, TCGA)"
  )

ggsave("results/figures/fig13_common_degs_volcano_highlight.png",
       p_volcano_common, width = 8, height = 6, dpi = CFG$project$dpi)
log_msg("Figure saved in: results/figures/fig13_common_degs_volcano_highlight.png")

  
# --- Save -----------------------------------------------------------------------

saveRDS(
  list(
    mb_df           = mb_df,
    mb_hot_prog     = mb_hot_prog,
    program_in_mb   = program_in_mb,
    wilcox_prog     = wilcox_prog,
    prog_summary    = prog_summary,
    dea_mb          = dea_mb,
    gsea_mb         = if (exists("gsea_mb")) gsea_mb else NULL,
    comp_tab        = if (exists("comp_tab")) comp_tab else NULL
  ),
  "data/processed/metabric_phenotypes.rds"
)

log_msg("✓ Saved: data/processed/metabric_phenotypes.rds")
log_msg("=== 05_metabric_validation.R END ===")
