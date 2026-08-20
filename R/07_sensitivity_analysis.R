# =============================================================================
# 07_sensitivity_analysis.R
# Sensitivity analysis of the classification threshold
# =============================================================================

source("R/00_setup.R")

suppressMessages({
  library(limma)
  library(fgsea)
  library(msigdbr)
  library(dplyr)
  library(ggplot2)
})

LOG_FILE <- "logs/07_sensitivity_analysis.log"

log_msg("=== 07_sensitivity_analysis.R START ===")
log_msg("Sensitivity analysis of classification thresholds")

# --- Load TCGA ---------------------------------------------------------------

pheno <- readRDS("data/processed/tcga_phenotypes.rds")

master_df <- pheno$master_df
exp_tnbc  <- pheno$exp_tnbc

if (!is.null(pheno$genes_excluded)) {
  genes_excluded <- pheno$genes_excluded
} else {
  t_markers <- c(
    "CD8A", "CD8B", "CD3D", "CD3E", "CD3G",
    "GZMA", "GZMB", "GZMK", "PRF1", "NKG7", "CD27", "IFNG"
  )
  
  genes_excluded <- unique(
    c(SIGS$tis, SIGS$stromal, t_markers)
  )
}

# --- Phenotype classification ------------------------------------------------

classify_phenotypes <- function(df, tis_low, tis_high,
                                str_low, str_high) {
  
  immune <- ifelse(
    df$TIS_Score >= tis_high, "Hot",
    ifelse(df$TIS_Score < tis_low, "Cold", NA)
  )
  
  stromal <- ifelse(
    df$Stromal_Score >= str_high, "Excluded",
    ifelse(df$Stromal_Score < str_low, "NonExcluded", NA)
  )
  
  ph <- ifelse(
    is.na(immune) | is.na(stromal),
    NA,
    paste(immune, stromal, sep = "_")
  )
  
  factor(
    ph,
    levels = c(
      "Hot_Excluded",
      "Hot_NonExcluded",
      "Cold_Excluded",
      "Cold_NonExcluded"
    )
  )
}

# --- Define threshold scenarios ----------------------------------------------

scenarios <- list(
  Reference = c(0.50, 0.50),
  P40       = c(0.40, 0.40),
  P60       = c(0.60, 0.60)
)

# --- Generate phenotype assignments ------------------------------------------

pheno_assignments <- list()

for (nm in names(scenarios)) {
  
  p <- scenarios[[nm]]
  
  tis_low  <- quantile(master_df$TIS_Score, p[1], na.rm = TRUE)
  tis_high <- quantile(master_df$TIS_Score, p[2], na.rm = TRUE)
  
  str_low  <- quantile(master_df$Stromal_Score, p[1], na.rm = TRUE)
  str_high <- quantile(master_df$Stromal_Score, p[2], na.rm = TRUE)
  
  pheno_assignments[[nm]] <- classify_phenotypes(
    master_df,
    tis_low,
    tis_high,
    str_low,
    str_high
  )
}

# --- Differential expression -------------------------------------------------

run_dea <- function(phenotype_vec) {
  
  keep <- phenotype_vec %in% c(
    "Hot_Excluded",
    "Hot_NonExcluded"
  )
  
  if (sum(keep) < 10) {
    return(NULL)
  }
  
  grp <- droplevels(
    factor(
      phenotype_vec[keep],
      levels = c("Hot_NonExcluded", "Hot_Excluded")
    )
  )
  
  if (length(levels(grp)) < 2) {
    return(NULL)
  }
  
  expr <- exp_tnbc[, keep]
  expr <- expr[!rownames(expr) %in% genes_excluded, ]
  
  design <- model.matrix(~ 0 + grp)
  colnames(design) <- c("NonExcluded", "Excluded")
  
  fit <- lmFit(expr, design)
  
  contrast <- makeContrasts(
    Excluded - NonExcluded,
    levels = design
  )
  
  fit2 <- eBayes(
    contrasts.fit(fit, contrast)
  )
  
  topTable(
    fit2,
    number = Inf,
    sort.by = "none"
  )
}

# --- Reference DEA -----------------------------------------------------------

dea_reference <- run_dea(
  as.character(pheno_assignments$Reference)
)

# --- Compare logFC across thresholds ----------------------------------------

logfc_results <- list()

for (nm in c("P40", "P60")) {
  
  dea <- run_dea(
    as.character(pheno_assignments[[nm]])
  )
  
  if (is.null(dea)) {
    log_msg(
      "DEA not evaluable for ", nm,
      " (insufficient sample size)",
      level = "WARN"
    )
    next
  }
  
  common_genes <- intersect(
    rownames(dea_reference),
    rownames(dea)
  )
  
  logfc_cor <- cor(
    dea_reference[common_genes, "logFC"],
    dea[common_genes, "logFC"],
    method = "pearson"
  )
  
  logfc_results[[nm]] <- data.frame(
    Scenario = nm,
    N_genes = length(common_genes),
    logFC_correlation = round(logfc_cor, 3)
  )
  
  log_msg(
    nm,
    " vs Reference: logFC correlation = ",
    round(logfc_cor, 3)
  )
}

logfc_results <- bind_rows(logfc_results)

write.csv(
  logfc_results,
  "results/tables/table08_sensitivity_logFC.csv",
  row.names = FALSE
)

# --- GSEA --------------------------------------------------------------------

# Hallmark gene sets

hallmark_df <- msigdbr(
  species = "Homo sapiens",
  collection = "H"
)

hallmark <- split(
  hallmark_df$gene_symbol,
  hallmark_df$gs_name
)

# Pathways relevant to the biological interpretation

key_pathways <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_G2M_CHECKPOINT"
)

gsea_results <- list()

for (nm in names(pheno_assignments)) {
  
  dea <- run_dea(
    as.character(pheno_assignments[[nm]])
  )
  
  if (is.null(dea)) {
    next
  }
  
  ranks <- dea$t
  names(ranks) <- rownames(dea)
  
  ranks <- sort(
    ranks[is.finite(ranks)],
    decreasing = TRUE
  )
  
  gsea <- suppressWarnings(
    fgsea(
      hallmark,
      ranks,
      minSize = 10,
      maxSize = 500
    )
  )
  
  gsea_key <- gsea %>%
    filter(pathway %in% key_pathways) %>%
    select(pathway, NES, padj) %>%
    mutate(Scenario = nm)
  
  gsea_results[[nm]] <- gsea_key
}

gsea_results <- bind_rows(gsea_results)

write.csv(
  gsea_results,
  "results/tables/table08_sensitivity_GSEA.csv",
  row.names = FALSE
)

# --- Plot logFC concordance --------------------------------------------------

plot_data <- list()

for (nm in c("P40", "P60")) {
  
  dea <- run_dea(
    as.character(pheno_assignments[[nm]])
  )
  
  if (is.null(dea)) {
    next
  }
  
  common_genes <- intersect(
    rownames(dea_reference),
    rownames(dea)
  )
  
  plot_data[[nm]] <- data.frame(
    Scenario = nm,
    logFC_reference = dea_reference[
      common_genes, "logFC"
    ],
    logFC_alternative = dea[
      common_genes, "logFC"
    ]
  )
}

plot_data <- bind_rows(plot_data)

p_logfc <- ggplot(
  plot_data,
  aes(
    x = logFC_reference,
    y = logFC_alternative
  )
) +
  geom_point(alpha = 0.25, size = 0.7) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  facet_wrap(~ Scenario) +
  labs(
    title = "Concordance of differential expression across thresholds",
    x = "logFC (reference, p50)",
    y = "logFC (alternative threshold)"
  ) +
  theme_classic()

ggsave(
  "results/figures/fig09_sensitivity_logFC_concordance.png",
  p_logfc,
  width = 9,
  height = 5,
  dpi = CFG$project$dpi
)

log_msg(
  "Figure saved: results/figures/fig09_sensitivity_logFC_concordance.png"
)

log_msg("=== 07_sensitivity_analysis.R END ===")
