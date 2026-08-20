# =============================================================================
# 06_ispy2_exploratory.R
# EXPLORATORY analysis: association of the four phenotypes with pCR in I-SPY2.
# =============================================================================

source("R/00_setup.R")
LOG_FILE <- "logs/06_ispy2_exploratory.log"

log_msg("=== 06_ispy2_exploratory.R START ===")
log_msg("NOTE: exploratory analysis — limited N, interpret with caution")
log_msg("APPROACH: main analysis restricted to the durvalumab+olaparib arm")


# --- Load I-SPY2 -----------------------------------------------------------

log_msg("Loading I-SPY2 data (GSE173839)...")
ispy2_expr <- read.table(PATHS$ispy2_expr,
                         header = TRUE, sep = "\t",
                         row.names = 1, check.names = FALSE)
ispy2_meta <- read.csv(PATHS$ispy2_meta)

meta_clean <- ispy2_meta %>%
  filter(pCR.status %in% c(0, 1)) %>%
  mutate(
    pCR         = factor(pCR.status, levels = c(0, 1),
                         labels = c("No_pCR", "pCR")),
    pCR_numeric = as.numeric(pCR.status),
    ResearchID  = as.character(ResearchID)
  )

colnames(ispy2_expr) <- as.character(colnames(ispy2_expr))
common_is <- intersect(meta_clean$ResearchID, colnames(ispy2_expr))
ispy2_clean <- ispy2_expr[, common_is]
meta_clean  <- meta_clean %>% filter(ResearchID %in% common_is)

log_msg("Expression dimensions: ", nrow(ispy2_clean), " genes x ",
        ncol(ispy2_clean), " samples")
log_msg("pCR available: ", nrow(meta_clean),
        " | Overall rate: ", round(mean(meta_clean$pCR_numeric), 2))
log_msg("Arms (all): ", paste(names(table(meta_clean$Arm)),
                                  table(meta_clean$Arm), sep = "=", collapse = " | "))


# --- Compute scores ---------------------------------------------------------
# Scores and phenotype cutoffs are computed using all I-SPY2 samples.
# The durvalumab arm is only used for the downstream pCR analysis.

log_msg("Computing TIS score in I-SPY2...")
tis_is     <- compute_signature_score(ispy2_clean, SIGS$tis,     "TIS_ISPY2")

log_msg("Computing Stromal score in I-SPY2...")
stromal_is <- compute_signature_score(ispy2_clean, SIGS$stromal, "Stromal_ISPY2")

is_df <- data.frame(
  ResearchID    = common_is,
  TIS_Score     = tis_is,
  Stromal_Score = stromal_is
) %>%
  left_join(meta_clean %>% dplyr::select(ResearchID, pCR, pCR_numeric, Arm),
            by = "ResearchID")


# --- Classify phenotypes (I-SPY2's own medians) -----------------------------

tis_med_is     <- median(is_df$TIS_Score,     na.rm = TRUE)
stromal_med_is <- median(is_df$Stromal_Score, na.rm = TRUE)

is_df <- is_df %>%
  mutate(
    Immune_Group  = ifelse(TIS_Score     >= tis_med_is,     "Hot",      "Cold"),
    Stromal_Group = ifelse(Stromal_Score >= stromal_med_is, "Excluded", "NonExcluded"),
    Phenotype     = factor(
      paste(Immune_Group, Stromal_Group, sep = "_"),
      levels = c("Hot_Excluded", "Hot_NonExcluded",
                 "Cold_Excluded", "Cold_NonExcluded")
    )
  )

log_msg("TIS-Stromal correlation in I-SPY2: rho = ",
        round(cor(is_df$TIS_Score, is_df$Stromal_Score, method = "spearman"), 3))

arm_levels <- unique(is_df$Arm)
log_msg("Arm labels found: ", paste(arm_levels, collapse = " | "))
durva_label <- arm_levels[grepl("durva", arm_levels, ignore.case = TRUE)]
if (length(durva_label) != 1)
  log_msg("WARNING: could not identify a single durvalumab arm (",
          paste(durva_label, collapse = ", "), ")", level = "WARN")

log_msg("Phenotype distribution in I-SPY2:")
print(is_df %>% count(Phenotype))

# --- Main analysis ---------------------------------

is_durva <- is_df %>% filter(Arm %in% durva_label)
log_msg("--- MAIN ANALYSIS: durvalumab arm (n = ", nrow(is_durva), ") ---")

pcr_durva <- is_durva %>%
  group_by(Phenotype) %>%
  summarise(
    n        = n(),
    n_pCR    = sum(pCR_numeric, na.rm = TRUE),
    pCR_rate = round(mean(pCR_numeric, na.rm = TRUE), 3),
    .groups  = "drop"
  ) %>%
  arrange(desc(pCR_rate))

log_msg("pCR rate by phenotype (durvalumab arm):")
print(pcr_durva)
write.csv(pcr_durva,
          "results/tables/table07_ispy2_pcr_durvalumab.csv",
          row.names = FALSE)

# Main test: Hot_Excluded vs Hot_NonExcluded in the durvalumab arm
is_durva_hot <- is_durva %>% filter(Immune_Group == "Hot")
log_msg("Hot subgroup in durvalumab arm: n = ", nrow(is_durva_hot))

if (nrow(is_durva_hot) >= 8 &&
    length(unique(is_durva_hot$Stromal_Group)) == 2 &&
    length(unique(is_durva_hot$pCR)) == 2) {
  tab_main <- table(is_durva_hot$Stromal_Group, is_durva_hot$pCR)
  fisher_main <- fisher.test(tab_main)
  log_msg("** MAIN TEST ** Hot_Excluded vs Hot_NonExcluded | durvalumab: ",
          "OR = ", round(fisher_main$estimate, 2),
          " | 95% CI: ", round(fisher_main$conf.int[1], 2), "-",
          round(fisher_main$conf.int[2], 2),
          " | p = ", signif(fisher_main$p.value, 3))
} else {
  fisher_main <- NULL
  log_msg("WARNING: Hot subgroup insufficient or lacking variation for Fisher's test",
          level = "WARN")
}


# --- Main figure --------------------

subtitle_main <- paste0("Durvalumab+olaparib arm | n=", nrow(is_durva),
                        " TNBC | Wilson 95% CI")
if (!is.null(fisher_main)) {
  subtitle_main <- paste0(
    subtitle_main,
    "\nHot_Excluded vs Hot_NonExcluded: OR=", round(fisher_main$estimate, 2),
    ", p=", signif(fisher_main$p.value, 3), " (exploratory)"
  )
}

p_durva <- ggplot(pcr_durva,
                  aes(x = reorder(Phenotype, pCR_rate),
                      y = pCR_rate,
                      fill = Phenotype)) +
  geom_col(alpha = 0.85, width = 0.6) +
  geom_errorbar(
    aes(
      ymin = pmax(0, (n_pCR + 1.92 - 1.96 * sqrt(n_pCR * (n - n_pCR) / n + 0.96)) /
                    (n + 3.84)),
      ymax = pmin(1, (n_pCR + 1.92 + 1.96 * sqrt(n_pCR * (n - n_pCR) / n + 0.96)) /
                    (n + 3.84))
    ),
    width = 0.2, color = "grey30"
  ) +
  geom_text(aes(label = paste0(round(pCR_rate * 100), "%\n(n=", n, ")")),
            vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = PHENO_COLORS) +
  scale_y_continuous(labels = scales::percent_format(),
                     limits = c(0, 0.85), expand = c(0, 0)) +
  coord_flip() +
  labs(
    title    = "pCR rate by phenotype — I-SPY2, immunotherapy arm",
    subtitle = subtitle_main,
    x        = NULL,
    y        = "pCR rate",
    fill     = NULL
  ) +
  theme(legend.position = "none")

ggsave("results/figures/fig07_ispy2_pcr_durvalumab.png",
       p_durva, width = 8, height = 5, dpi = CFG$project$dpi)
log_msg("MAIN figure saved: results/figures/fig07_ispy2_pcr_durvalumab.png")


# --- Score scatter (diagnostic, durvalumab arm only) -------------------------

p_is_scatter <- ggplot(is_durva,
                       aes(x = TIS_Score, y = Stromal_Score,
                           color = Phenotype, shape = pCR)) +
  geom_point(size = 2.5, alpha = 0.75) +
  geom_vline(xintercept = tis_med_is,     linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = stromal_med_is, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = PHENO_COLORS) +
  scale_shape_manual(values = c("No_pCR" = 1, "pCR" = 17)) +
  labs(
    title    = "Scores in I-SPY2 (durvalumab arm): two microenvironment axes",
    subtitle = "Triangle = pCR | Circle = No pCR",
    x        = "TIS Score",
    y        = "Stromal Score",
    color    = "Phenotype",
    shape    = "Response"
  )

ggsave("results/figures/fig08_ispy2_scores_scatter.png",
       p_is_scatter, width = 7, height = 5, dpi = CFG$project$dpi)
log_msg("Figure saved: results/figures/fig08_ispy2_scores_scatter.png")


# --- Save --------------------------------------------------------------------

saveRDS(
  list(
    is_df         = is_df,        
    is_durva      = is_durva,      
    pcr_durva     = pcr_durva,     
    fisher_main   = fisher_main
  ),
  file = "data/processed/ispy2_phenotypes.rds"
)

log_msg("✓ Saved: data/processed/ispy2_phenotypes.rds")
log_msg("REMINDER: exploratory analysis — state N limitations in the thesis")
log_msg("=== 06_ispy2_exploratory.R END ===")
