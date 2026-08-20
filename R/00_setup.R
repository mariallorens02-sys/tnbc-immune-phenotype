# =============================================================================
# 00_setup.R
# Global setup: libraries, configuration, helper functions and gene signatures
#
# Run with source("R/00_setup.R") from the project root.
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(pheatmap)
  library(broom)
  library(rstatix)
  library(pROC)
  library(limma)
  library(fgsea)
  library(msigdbr)
  library(yaml)
})

# --- Configuration ----------------------------------------------------
CFG <- yaml::read_yaml("config/config.yml")
set.seed(CFG$project$seed)
theme_set(theme_bw(base_size = 12))

# --- Paths --------------------------------------------------------------------
PATHS <- CFG$paths

# --- Gene signatures -----------------------------------------------------------
SIGS <- list(
  tis     = CFG$signatures$tis_ayers$genes,          
  stromal = CFG$signatures$stromal_exclusion$genes    
)
GENES_EXCLUDED_DEA <- unique(c(
  CFG$dea_exclusions$tis_genes,
  CFG$dea_exclusions$stromal_genes,
  CFG$dea_exclusions$extra
))

# --- Color palette (four phenotypes) -----------------------------------------
PHENO_COLORS <- c(
  "Hot_Excluded"     = CFG$colors$hot_excluded,
  "Hot_NonExcluded"  = CFG$colors$hot_non_excluded,
  "Cold_Excluded"    = CFG$colors$cold_excluded,
  "Cold_NonExcluded" = CFG$colors$cold_non_excluded
)

# --- Logger -------------------------------------------------------------
LOG_FILE <- NULL  # set in each script

log_msg <- function(..., level = "INFO") {
  msg <- paste0("[", level, " ", format(Sys.time(), "%H:%M:%S"), "] ",
                paste(..., sep = ""))
  message(msg)
  if (!is.null(LOG_FILE)) {
    cat(msg, "\n", file = LOG_FILE, append = TRUE)
  }
}

# --- Signature score ---------------------------
compute_signature_score <- function(expr_matrix, genes, sig_name = "sig") {
  present <- intersect(genes, rownames(expr_matrix))
  missing <- setdiff(genes, rownames(expr_matrix))
  if (length(missing) > 0)
    log_msg("Missing genes in ", sig_name, ": ",
            paste(missing, collapse = ", "), level = "WARN")
  if (length(present) == 0)
    stop("No gene from signature '", sig_name, "' found in the matrix.")
  sub  <- expr_matrix[present, , drop = FALSE]
  zmat <- t(scale(t(sub)))
  colMeans(zmat, na.rm = TRUE)
}

message("✓ Setup complete — seed: ", CFG$project$seed,
        " | TIS: ", length(SIGS$tis), " genes",
        " | Stromal: ", length(SIGS$stromal), " genes")
