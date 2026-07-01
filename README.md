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

## Tools & Packages
R · limma · ggplot2 · pheatmap · clusterProfiler · biomaRt · ape

---

## 👩‍🔬 Author
**Halireena Rushdiha Mohomed**  
MSc Bioinformatics, University of Sharjah  
Supervised by Dr. Reem
