# E-MTAB-2990-Potato-Phosphate-Transcriptomics
E-MTAB-2990-Potato-Phosphate-Transcriptomics. First systematic transcriptomic analysis of E-MTAB-2990 — potato phosphate stress response across four cultivars | MSc Bioinformatics internship project

# Transcriptomic Analysis of Potato Phosphate Stress Response
### First systematic analysis of E-MTAB-2990 (ArrayExpress)

**MSc Bioinformatics Internship Project**  
University of Sharjah | Supervised by Dr. Reem

---

## Overview
This repository contains a five-stage R pipeline for the analysis of
a publicly archived Agilent one-colour microarray dataset (E-MTAB-2990),
examining how four potato cultivars (*Solanum tuberosum*) respond to
phosphate starvation in root tissue. This is the first formal systematic
analysis of this dataset.

**Dataset:** E-MTAB-2990 (ArrayExpress, EBI)  
**Platform:** Agilent A-MEXP-2272 (Potato 60K array, 54,317 probes)  
**Cultivars:** Maris Piper · Pentland Dell · Stirling · 12601  
**Samples:** 24 root tissue samples (4 cultivars × 2 conditions × 3 replicates)

---

## Key Findings
- **285 DEGs** in Maris Piper under low phosphate (73% upregulated)
- Top gene **PGSC0003DMT400069516** — an ABCB19-family auxin efflux
  transporter upregulated ~79-fold; conserved across 21 Solanaceae
  species and ~120 million years of evolution
- **100% upregulation** of all 11 acid phosphatases and all 6 PHT1
  phosphate transporters detected
- **555 Maris Piper-unique DEGs** vs 137 shared across all four cultivars
  (Stage 4, all-24-sample 8-group model — this model pools variance across
  cultivars, so its Maris Piper DEG count is larger than the 285 from Stage 2)

---

## Pipeline Structure

| Stage | Notebook | Description |
|-------|----------|-------------|
| 1 | `01_load_QC_FIXED.Rmd` | Data loading, QC, normalisation, filtering |
| 2 | `02_differential_expression_FINAL.Rmd` | limma DEG analysis (Analysis A & B) |
| 3 | `03_Functional_Enrichment_Final.Rmd` | GO enrichment via BioMart/clusterProfiler |
| 4 | `04_advanced_analysis_CLEAN.Rmd` | Cross-cultivar 8-group model, heatmaps, co-expression |
| 5 | `05_evolutionary_analysis.Rmd` | Phylogenetic tree across 21 Solanaceae species |

---

## Requirements

- R ≥ 4.1 (developed with a recent R 4.x release)
- Internet access for Stage 3 (Ensembl Plants BioMart queries)

| Source | Packages |
|--------|----------|
| CRAN | rmarkdown, knitr, ggplot2, dplyr, RColorBrewer, pheatmap, ggrepel, ggVennDiagram, ape, BiocManager |
| Bioconductor | limma, clusterProfiler, biomaRt |

Install everything missing with:

```r
source("install_packages.R")
```

or manually:

```r
install.packages(c("rmarkdown", "knitr", "ggplot2", "dplyr", "RColorBrewer",
                   "pheatmap", "ggrepel", "ggVennDiagram", "ape", "BiocManager"))
BiocManager::install(c("limma", "clusterProfiler", "biomaRt"))
```

Each notebook ends with `sessionInfo()`, so the rendered HTML records the
exact package versions used.

---

## Data setup

All paths in the notebooks are **relative to the repository root**. Before
running, create this layout (the `data/raw/` and `results/` folders are not
tracked by git):

```
E-MTAB-2990-Potato-Phosphate-Transcriptomics/
├── metadata/
│   └── E-MTAB-2990.sdrf.txt      # sample sheet from ArrayExpress
├── data/
│   ├── raw/                      # the 24 Agilent Feature Extraction .txt files
│   └── tree.nwk                  # (Stage 5 only) NCBI BLAST tree, Newick format
└── results/                      # created by Stage 1; all outputs go here
```

1. Open [E-MTAB-2990 on BioStudies/ArrayExpress](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-2990).
2. Download the SDRF file (`E-MTAB-2990.sdrf.txt`) into `metadata/`.
3. Download the raw data archive(s) and unzip the 24 `.txt` files into
   `data/raw/`. The file names must match the `Array Data File` column of
   the SDRF.
4. *(Stage 5 only)* BLAST the `PGSC0003DMT400069516` mRNA sequence
   (blastn, `refseq_rna`, organism = Viridiplantae) at NCBI, open
   **Distance tree of results**, download it in Newick format and save it
   as `data/tree.nwk`.

---

## How to run

The stages must be run **in order** — each one loads the `.RData` files
written by the previous stage(s) from `results/`.

```r
# from the repository root
for (f in c("01_load_QC_FIXED.Rmd",
            "02_differential_expression_FINAL.Rmd",
            "03_Functional_Enrichment_Final.Rmd",
            "04_advanced_analysis_CLEAN.Rmd",
            "05_evolutionary_analysis.Rmd")) {
  rmarkdown::render(f)
}
```

Or open a notebook in RStudio and click **Knit**. Stage 5 only needs
`data/tree.nwk` and can be run independently of Stages 1-4.

---

## Expected outputs

Each notebook renders to an `.html` report next to the `.Rmd`. Files written
to `results/`:

| Stage | Key outputs |
|-------|-------------|
| 1 | `normalised_data.RData` (filtered `EList` + `targets`), `MP_interarray_correlations.csv` |
| 2 | `DEG_results.RData`, `DEGs_A_with_outlier.csv`, `DEGs_B_without_outlier.csv`, `all_results_A.csv`, `all_results_B.csv`, `MP_PCA_scores.csv` |
| 3 | `Stage3_complete.RData`, `GO_enrichment_nominal_A.csv` / `_B.csv`, `Key_genes_annotated.csv`, `potato_ensembl_annotation.csv`, `full_array_annotation.csv`, `DEG_PGSC_IDs_A.csv` / `_B.csv`, GO dotplots, Venn and category plots (`.png`) |
| 4 | `Stage4_complete.RData`, per-cultivar DEG tables (`MP_`, `PD_`, `St_`, `12601_DEGs_stage4.csv`), `Common_response_all_cultivars.csv`, `MP_unique_DEGs.csv`, 10 figures (`Heatmap_top40.png`, `Volcano_plot.png`, `Venn_4cultivar.png`, …) |
| 5 | `phylogenetic_tree_final.png` |

---

## Data Availability
The raw microarray data used in this project are publicly available 
at ArrayExpress under accession number 
[E-MTAB-2990](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-2990).  
Data were originally deposited by Pete Hedley, James Hutton Institute (2014).  
No raw data files are hosted in this repository.

---

## Acknowledgements
Sincere thanks to **Dr. Reem** (University of Sharjah, Research Institute 
for Medical and Health Sciences) for her supervision, guidance, and 
support throughout this internship project.  
Thanks also to lab colleagues Fathima Nubla Latheef, Zoha Shaikh, 
and Omar Mohammad for their encouragement and support.

## 👩‍🔬 Author
**Halireena Rushdiha Mohomed**  
MSc Bioinformatics, University of Birmingham Dubai
Supervised by Dr. Reem
