# TNBC Immune Phenotypes

Transcriptomic classification of triple-negative breast cancer (TNBC) into four tumor microenvironment phenotypes using two biological axes:

- **Immune activity:** Hot vs Cold (Tumor Inflammation Signature, TIS)
- **Stromal exclusion:** Excluded vs NonExcluded (TGF-β/CAF stromal signature)

Crossing both axes defines four phenotypes:

| | NonExcluded | Excluded |
|---|---|---|
| **Cold** | Cold_NonExcluded | Cold_Excluded |
| **Hot** | Hot_NonExcluded | Hot_Excluded |

The workflow includes phenotype classification, differential expression analysis, Hallmark GSEA, external validation in METABRIC, exploratory validation in I-SPY2, sensitivity analysis and generation of publication-ready figures and tables.

---

## Project structure

```
tnbc-immune-phenotype/
│
├── R/
│   ├── 00_setup.R
│   ├── 01_load_tcga.R
│   ├── 02_define_axes_phenotypes.R
│   ├── 03_dea_hot_excluded_vs_nonexcluded.R
│   ├── 04_gsea_hallmark.R
│   ├── 05_metabric_validation.R
│   ├── 06_ispy2_exploratory.R
│   └── 07_sensitivity_analysis.R
│
├── config/
│   └── config.yml
│
├── data/
│   ├── raw/
│   └── processed/
│
├── results/
│   ├── figures/
│   └── tables/
│
├── logs/
└── README.md
```

---

## Pipeline

Run the scripts in the following order:

```r
source("R/00_setup.R")
source("R/01_load_tcga.R")
source("R/02_define_axes_phenotypes.R")
source("R/03_dea_hot_excluded_vs_nonexcluded.R")
source("R/04_gsea_hallmark.R")
source("R/05_metabric_validation.R")
source("R/06_ispy2_exploratory.R")
source("R/07_sensitivity_analysis.R")
```

Each script:

- reads the output from the previous step (`.rds`);
- writes a log file;
- exports figures and tables;
- saves intermediate results in `data/processed/`.

---

## Data

The pipeline uses three public datasets:

- **TCGA-BRCA** – discovery cohort
- **METABRIC** – independent validation
- **I-SPY2** – exploratory validation

Raw data are **not included** in this repository.

Update the dataset locations in `config/config.yml` before running the analysis.

---

## Configuration

All user-defined settings are stored in `config/config.yml`, including:

- dataset paths;
- gene signatures;
- DEA exclusion genes;
- phenotype colours;
- project parameters.

---

## Output

The pipeline automatically generates:

- publication-ready figures (`results/figures/`);
- summary tables (`results/tables/`);
- intermediate R objects (`data/processed/`);
- execution logs (`logs/`).
