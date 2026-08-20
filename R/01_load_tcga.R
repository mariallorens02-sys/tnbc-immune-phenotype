# =============================================================================
# 01_load_tcga.R
# Load and prepare the TCGA-TNBC cohort for downstream analysis 
# =============================================================================

source("R/00_setup.R")
LOG_FILE <- "logs/01_load_tcga.log"

log_msg("=== 01_load_tcga.R START ===")


# --- Expression matrix ---------------------------------------------------

log_msg("Loading TCGA expression matrix...")
exp_raw <- read.table(PATHS$tcga_expr,
                      header = TRUE, sep = "\t", check.names = FALSE)

genes      <- exp_raw[, 1]
exp_matrix <- as.matrix(exp_raw[, -1])
rownames(exp_matrix) <- genes
storage.mode(exp_matrix) <- "numeric"

log_msg("Raw dimensions: ", nrow(exp_matrix), " genes x ",
        ncol(exp_matrix), " samples")


# --- Clinical data -----------------------------------

log_msg("Loading clinical data and filtering TNBC (ER-/PR-/HER2-)...")
clinical_raw <- read.table(PATHS$tcga_clinical,
                           header = TRUE, sep = "\t",
                           fill = TRUE, quote = "")

clinical <- clinical_raw %>%
  dplyr::select(
    sample,
    ER        = ER_Status_nature2012,
    PR        = PR_Status_nature2012,
    HER2      = HER2_Final_Status_nature2012,
    Age       = age_at_initial_pathologic_diagnosis,
    Stage     = AJCC_Stage_nature2012,
    Menopause = menopause_status,
    OS        = OS,
    OS.time   = OS.time,
    PFI       = PFI,
    PFI.time  = PFI.time
  ) %>%
  filter(ER == "Negative", PR == "Negative", HER2 == "Negative") %>%
  mutate(
    Menopause_Clean = case_when(
      grepl("Post",  Menopause) ~ "Post",
      grepl("Pre|Peri", Menopause) ~ "Pre/Peri",
      TRUE ~ NA_character_
    ),
    Stage_Simple = case_when(
      grepl("Stage I$|Stage IA|Stage IB", Stage) ~ "I",
      grepl("Stage II",  Stage) ~ "II",
      grepl("Stage III", Stage) ~ "III",
      grepl("Stage IV",  Stage) ~ "IV",
      TRUE ~ NA_character_
    )
  )

log_msg("TNBC patients identified: ", nrow(clinical))


# --- Match expression and clinical data  ----------------------------

common_samples <- intersect(clinical$sample, colnames(exp_matrix))
stopifnot(length(common_samples) > 0)
log_msg("Samples with expression + clinical data: ", length(common_samples))

exp_tnbc <- exp_matrix[, common_samples]
clinical  <- clinical %>% filter(sample %in% common_samples)


# --- Quality checks -------------------------------------------------

expr_range <- range(exp_tnbc, na.rm = TRUE)
log_msg("Expression range: [", round(expr_range[1], 2), ", ",
        round(expr_range[2], 2), "] — expected 0-20 if already log2")

if (expr_range[2] > 30) {
  log_msg("WARNING: range suggests linear scale; consider log2(x+1)",
          level = "WARN")
}


# --- Check presence of signature genes -----------------------------------

tis_present     <- intersect(SIGS$tis,     rownames(exp_tnbc))
stromal_present <- intersect(SIGS$stromal, rownames(exp_tnbc))

log_msg("TIS genes present: ", length(tis_present), "/", length(SIGS$tis))
log_msg("Stromal genes present: ",
        length(stromal_present), "/", length(SIGS$stromal))

if (length(tis_present) < 15)
  log_msg("WARNING: <15 TIS genes available — check data",
          level = "WARN")


# --- Save -----------------------------------------------------------------

saveRDS(
  list(
    exp_tnbc = exp_tnbc,
    clinical = clinical,
    samples  = common_samples
  ),
  file = "data/processed/tcga_tnbc_inputs.rds"
)

log_msg("✓ Saved: data/processed/tcga_tnbc_inputs.rds")
log_msg("  Final dimensions: ", nrow(exp_tnbc), " genes x ",
        ncol(exp_tnbc), " TNBC samples")
log_msg("=== 01_load_tcga.R END ===")
