# =============================================================================
# 04_gsea_hallmark.R
# Hallmark GSEA using the ranked differential expression results.
# =============================================================================

source("R/00_setup.R")
LOG_FILE <- "logs/04_gsea_hallmark.log"

log_msg("=== 04_gsea_hallmark.R START ===")

data_in     <- readRDS("data/processed/tcga_dea.rds")
dea_results <- data_in$dea_results


# --- Build ranking ---------------------------------------------------------

ranks <- setNames(dea_results$t, dea_results$Gene)
ranks <- sort(ranks, decreasing = TRUE)

log_msg("Genes in ranking: ", length(ranks))


# --- Load Hallmark gene sets ------------------------------------------------

log_msg("Loading Hallmark gene sets (MSigDB)...")
h_df <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_list <- split(h_df$gene_symbol, h_df$gs_name)
log_msg("Gene sets available: ", length(hallmark_list))


# --- GSEA -------------------------------------------------------------------

set.seed(CFG$project$seed)
gsea_res <- fgsea(
  pathways = hallmark_list,
  stats    = ranks,
  minSize  = 15,
  maxSize  = 500,
)

log_msg("GSEA complete — gene sets evaluated: ", nrow(gsea_res))

gsea_clean <- as.data.frame(gsea_res) %>%
  filter(padj < 0.05) %>%
  mutate(
    pathway_label = gsub("^HALLMARK_", "", pathway) %>% gsub("_", " ", .),
  ) %>%
  arrange(desc(NES))

log_msg("Significant pathways (padj < 0.05): ", nrow(gsea_clean))

write.csv(
  gsea_clean %>% dplyr::select(-leadingEdge),
  "results/tables/table04_gsea_results.csv",
  row.names = FALSE
)


# --- Figure: NES barplot -----------------------------------------------------
p_gsea <- ggplot(gsea_clean,
                 aes(x = reorder(pathway_label, NES),
                     y = NES,
                     )) +
  geom_col(alpha = 0.85) +
  coord_flip() +
  labs(
    title    = "Hallmark GSEA",
    x        = NULL,
    y        = "Normalized Enrichment Score (NES)",
    fill     = NULL
  ) +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 9))

ggsave("results/figures/fig04_gsea_hallmark.png",
       p_gsea, width = 11, height = max(6, nrow(gsea_clean) * 0.35),
       dpi = CFG$project$dpi)
log_msg("Figure saved: results/figures/fig04_gsea_hallmark.png")


# --- Save ---------------------------------------------------------------------

saveRDS(
  list(
    gsea_res   = gsea_res,
    gsea_clean = gsea_clean,
  ),
  file = "data/processed/tcga_gsea.rds"
)

log_msg("✓ Saved: data/processed/tcga_gsea.rds")
log_msg("=== 04_gsea_hallmark.R END ===")
