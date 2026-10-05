# Transcriptomic Analysis of Potato Phosphate Stress Response (E-MTAB-2990)

**MSc Bioinformatics internship project** · Halireena Rushdiha Mohomed · supervised by Dr. Reem

A five-stage R Markdown pipeline that re-analyses a public Agilent microarray
dataset (ArrayExpress **E-MTAB-2990**) to find which genes potato roots switch
on or off when they are starved of phosphate, and whether that response
differs between four potato cultivars.

---

## Contents

1. [Project summary](#1-project-summary)
2. [Key findings](#2-key-findings)
3. [Pipeline at a glance](#3-pipeline-at-a-glance)
4. [Quick start](#4-quick-start)
5. [Get the data](#5-get-the-data)
6. [How to run](#6-how-to-run)
7. [Expected outputs](#7-expected-outputs)
8. [Glossary](#8-glossary)
9. [Troubleshooting](#9-troubleshooting)
10. [Known limitations and analysis caveats](#10-known-limitations-and-analysis-caveats)
11. [Data provenance and citation](#11-data-provenance-and-citation)
12. [Acknowledgements, author and licence](#12-acknowledgements-author-and-licence)

---

## 1. Project summary

**Biological question.** Phosphorus is essential for plant growth, but most
soil phosphate is locked in forms roots cannot take up, so crops rely on heavy
fertiliser use. When phosphate is scarce, plants switch on a "phosphate
starvation response", for example phosphate transporters and acid
phosphatases that release phosphate from organic matter. This project asks:

1. Which genes change expression in **Maris Piper** roots grown at low versus
   regular phosphate?
2. What biological functions do those genes have?
3. Is the response shared by all four cultivars, or specific to Maris Piper?

**Experimental design** (as recorded in the dataset's SDRF sample sheet):
root tissue, one-colour Agilent microarray, **24 arrays**.

| Cultivar | Low P (phosphate-starved) | Reg P (regular phosphate) | Arrays |
|---|:-:|:-:|:-:|
| Maris Piper | 3 | 3 | 6 |
| Pentland Dell | 3 | 3 | 6 |
| Stirling | 3 | 3 | 6 |
| 12601 (breeding line) | 3 | 3 | 6 |
| **Total** | **12** | **12** | **24** |

- **Platform:** Agilent potato gene-expression array, ArrayExpress design
  A-MEXP-2272 (one-colour, Cy3). Probes are named with PGSC potato transcript IDs.
- **Main comparison:** Low P vs Reg P within Maris Piper (Stage 2), then the
  same comparison in every cultivar (Stage 4).
- **One borderline outlier:** one Maris Piper Reg P array correlates slightly
  less well with its replicates. Stage 2 therefore runs the analysis with it
  (Analysis A, main result) and without it (Analysis B, sensitivity check).

---

## 2. Key findings

These numbers come from the author's rendered notebooks. This README does not
recompute them, so check them against your own HTML reports after running the
pipeline (see [section 10](#10-known-limitations-and-analysis-caveats)).

| Finding | Source | Figure to look at |
|---|---|---|
| **285 DEGs** in Maris Piper, Low P vs Reg P, about 73% up-regulated | Stage 2, Analysis A | volcano plots in the Stage 2 HTML report |
| The DEG list is largely the same when the outlier array is dropped (Analysis B: 227 DEGs) | Stage 2 A vs B comparison | `results/Venn_DEG_overlap_AB.png` |
| All 11 acid phosphatases and all 6 phosphate transporters (PHT1 family) found among the DEGs are **up-regulated**, the textbook starvation response | Stage 3 keyword annotation | `results/Heatmap_key_genes.png` |
| Root-hair development and other GO terms are enriched at nominal p < 0.05, but **none survive multiple-testing correction** | Stage 3 | `results/GO_dotplot_Analysis_A.png` |
| In a model fitted to all 24 arrays, **555 DEGs are unique to Maris Piper** and **137 are shared by all four cultivars** | Stage 4 | `results/Venn_4cultivar.png`, `results/Response_classification.png` |
| A strongly induced ABC transporter (ABCB family, `PGSC0003DMT400069516`) has close BLAST matches across the nightshade family (Solanaceae) | Stage 5 | `results/phylogenetic_tree_final.png` |

> Stage 4 fits one model to all 24 arrays, so its error estimate is pooled
> across cultivars and its Maris Piper DEG count is larger than Stage 2's 285.
> The two counts come from different models, so they do not contradict each other.

---

## 3. Pipeline at a glance

```
 ArrayExpress E-MTAB-2990                         NCBI BLAST (web)
 (SDRF + 24 raw .txt files)                              │
          │                                              ▼
          ▼                                         data/tree.nwk
 [1] 01_load_QC_FIXED ─────► normalised_data.RData       │
          │                        │                     ▼
          ▼                        │        [5] 05_evolutionary_analysis
 [2] 02_differential_expression_FINAL            (independent of 1-4)
          │ ─────► DEG_results.RData
          ▼
 [3] 03_Functional_Enrichment_Final   (+ internet: Ensembl Plants BioMart)
          │ ─────► Stage3_complete.RData, Key_genes_annotated.csv
          ▼
 [4] 04_advanced_analysis_CLEAN ─────► Stage4_complete.RData + figures
```

| Stage | Notebook | What it does | Reads | Writes (hand-off files in **bold**) | Approx. time |
|---|---|---|---|---|---|
| 1 | `01_load_QC_FIXED.Rmd` | Load raw arrays, QC plots, background correction, quantile normalisation, low-signal filtering, PCA, outlier check | `metadata/E-MTAB-2990.sdrf.txt`, `data/raw/*.txt` | **`results/normalised_data.RData`** | 1-3 min |
| 2 | `02_differential_expression_FINAL.Rmd` | limma DE for Maris Piper, Low P vs Reg P: Analysis A (6 arrays) and B (5 arrays) | `normalised_data.RData` | **`results/DEG_results.RData`** + CSV tables | < 1 min |
| 3 | `03_Functional_Enrichment_Final.Rmd` | Map probes to genes, GO enrichment (clusterProfiler), keyword-based gene-family annotation | `DEG_results.RData`, `normalised_data.RData`, one file in `data/raw/`, **internet** | **`results/Stage3_complete.RData`**, **`results/Key_genes_annotated.csv`** + GO tables and figures | 3-10 min (BioMart download) |
| 4 | `04_advanced_analysis_CLEAN.Rmd` | 8-group limma model (4 cultivars × 2 conditions), cross-cultivar comparison, heatmaps, Venn, co-expression | outputs of Stages 1-3 | `results/Stage4_complete.RData` + per-cultivar tables + 10 figures | 1-2 min |
| 5 | `05_evolutionary_analysis.Rmd` | Prune and plot the NCBI BLAST distance tree for the top ABC transporter | `data/tree.nwk` only | `results/phylogenetic_tree_final.png` | seconds |

The times are rough estimates for a laptop with ordinary broadband. The
author's notebooks quote 1-2 min for loading the raw files and 2-3 min for the
BioMart download. The one-off package installation (section 4) can take
10-30 minutes, mostly compiling Bioconductor dependencies.

Each notebook begins with a **"Stage N at a glance"** box that lists its
inputs, outputs and run time. If an input is missing, the notebook stops with
a message saying what to run first. Each notebook ends with `sessionInfo()`,
so the rendered HTML records the exact package versions.

---

## 4. Quick start

Requirements: **R ≥ 4.1**, [pandoc](https://pandoc.org/) (bundled with
RStudio), roughly **2 GB free disk space**, and internet access for the data
download and for Stage 3.

```bash
git clone https://github.com/halireena/E-MTAB-2990-Potato-Phosphate-Transcriptomics.git
cd E-MTAB-2990-Potato-Phosphate-Transcriptomics

Rscript install_packages.R   # 1. install R packages (once)
Rscript download_data.R      # 2. SDRF + 24 raw arrays -> metadata/ and data/raw/
#                              3. Stage 5 only: save data/tree.nwk (section 5, step 4)
Rscript -e 'for (f in sort(Sys.glob("0[1-4]_*.Rmd"))) rmarkdown::render(f, envir = new.env())'   # 4. Stages 1-4
```

| Source | Packages |
|--------|----------|
| CRAN | rmarkdown, knitr, ggplot2, dplyr, RColorBrewer, pheatmap, ggrepel, ggVennDiagram, ape, BiocManager |
| Bioconductor | limma, clusterProfiler, biomaRt |

The notebooks never install packages themselves. If a package is missing,
they stop and tell you to run `install_packages.R`.

---

## 5. Get the data

No data are stored in this repository. The notebooks expect these exact file
names and locations, relative to the repository root:

```
E-MTAB-2990-Potato-Phosphate-Transcriptomics/
├── 01_load_QC_FIXED.Rmd … 05_evolutionary_analysis.Rmd
├── install_packages.R
├── download_data.R
├── metadata/
│   ├── E-MTAB-2990.sdrf.txt     # REQUIRED: sample sheet (file name -> cultivar, condition)
│   └── E-MTAB-2990.idf.txt      # optional: experiment description, for reading
├── data/
│   ├── raw/                     # REQUIRED for Stages 1 and 3
│   │   ├── US10310381_253303310046_S01_GE1_107_Sep09_1_4.txt   # (example name)
│   │   └── …                    # 24 Agilent Feature Extraction .txt files, named exactly
│   │                            #   as in the SDRF "Array Data File" column, no sub-folders
│   └── tree.nwk                 # Stage 5 only (you create it, step 4 below)
└── results/                     # created by Stage 1; every output goes here
```

**Automatic (recommended):** run `Rscript download_data.R`. It downloads the
SDRF, reads the raw archive names from it, downloads and unzips the archives
into `data/raw/`, and checks that all 24 files listed in the SDRF are there.
It is safe to re-run. If EBI moves the files, set the environment variable
`EMTAB_BASE_URL` to the new folder, or download the files by hand.

**Manual:**

1. Open <https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-2990> and go to the **Files** section.
2. Download `E-MTAB-2990.sdrf.txt` (and, if you want it, `E-MTAB-2990.idf.txt`) into `metadata/`.
3. Download the **raw** data archive(s): `E-MTAB-2990.raw.1.zip` and any
   further `raw.N.zip`. Unzip the `.txt` files directly into `data/raw/`.
   Do not use "processed" files, because Stage 1 normalises from the raw
   Feature Extraction output itself.
4. *(Stage 5 only.)* Get the transcript sequence of `PGSC0003DMT400069516`
   (for example from Ensembl Plants or Spud DB). Run **blastn** at NCBI against
   `refseq_rna` with the organism limited to *Viridiplantae*, open **Distance
   tree of results**, download the tree in **Newick** format and save it as
   `data/tree.nwk`. BLAST databases keep growing, so your hits may differ from
   the author's.

What the files contain:

- **SDRF** (Sample and Data Relationship Format): a tab-separated table with
  one row per array. Stage 1 uses three of its columns: `Array Data File`,
  `Characteristics[cultivar]` and `Factor Value[growth condition]`
  (`Low P` / `Reg P`).
- **IDF** (Investigation Description Format): free-text description of the
  experiment, including who ran it and the protocols.
- **Raw `.txt` files**: one per array, written by Agilent Feature Extraction
  software. Each has a header block, then one row per probe with intensities
  and quality flags.

---

## 6. How to run

Run the stages **in order 1 → 2 → 3 → 4**. Each one loads `.RData` files
that earlier stages wrote to `results/`. Stage 5 is independent.

```r
# from the repository root, in R
for (f in c("01_load_QC_FIXED.Rmd",
            "02_differential_expression_FINAL.Rmd",
            "03_Functional_Enrichment_Final.Rmd",
            "04_advanced_analysis_CLEAN.Rmd",
            "05_evolutionary_analysis.Rmd")) {
  rmarkdown::render(f, envir = new.env())  # fresh environment per stage, so each
}                                          # stage really relies on its saved inputs
```

You can also open a notebook in RStudio and click **Knit**. Knit uses the
notebook's folder (the repository root) as the working directory, which is
what the relative paths expect.

Objects passed between stages:

| File | Written by | Loaded by | Objects inside |
|---|---|---|---|
| `results/normalised_data.RData` | Stage 1 | 2, 3, 4 | `norm_filtered` (limma `EList`: `$E` = log2 expression, `$genes` = probe annotation), `targets` (FileName, Cultivar, Condition) |
| `results/DEG_results.RData` | Stage 2 | 3, 4 | `degs_A`, `degs_B` (DEG tables), `results_A`, `results_B` (all probes), `mp_targets`, `mp_expr`, `mp_targets_B`, `mp_expr_B` |
| `results/Stage3_complete.RData` | Stage 3 | 4 | PGSC ID lists, `ego_A` / `ego_B` enrichment objects, `nom_A` / `nom_B`, `full_deg_anno`, `term2gene_full`, … |
| `results/Key_genes_annotated.csv` | Stage 3 | 4 | DEGs with BioMart descriptions and a functional `Category` |

The `.RData` files are git-ignored, so a fresh clone never contains them. You
must render the stages yourself.

---

## 7. Expected outputs

Each notebook renders to an `.html` report next to its `.Rmd`, which is the
easiest place to read the results. Files written to `results/`:

| Stage | Key outputs |
|-------|-------------|
| 1 | `normalised_data.RData`, `MP_interarray_correlations.csv` |
| 2 | `DEG_results.RData`, `DEGs_A_with_outlier.csv`, `DEGs_B_without_outlier.csv`, `all_results_A.csv`, `all_results_B.csv`, `MP_PCA_scores.csv` (and re-writes `MP_interarray_correlations.csv`) |
| 3 | `Stage3_complete.RData`, `Key_genes_annotated.csv`, `GO_enrichment_nominal_A.csv` / `_B.csv`, `potato_ensembl_annotation.csv`, `full_array_annotation.csv`, `DEG_PGSC_IDs_A.csv` / `_B.csv`, `GO_dotplot_Analysis_A.png` / `_B.png`, `Venn_DEG_overlap_AB.png`, `Functional_categories_barplot.png` |
| 4 | `Stage4_complete.RData`, `MP_` / `PD_` / `St_` / `12601_DEGs_stage4.csv`, `Common_response_all_cultivars.csv`, `MP_unique_DEGs.csv`, `Heatmap_top40.png`, `Cross_cultivar_heatmap.png`, `Volcano_plot.png`, `Top_DEGs_barplot.png`, `Venn_4cultivar.png`, `Response_classification.png`, `Heatmap_key_genes.png`, `Category_summary_plot.png`, `Coexpression_matrix.png`, `Integration_summary.png` |
| 5 | `phylogenetic_tree_final.png` |

**Reading a DEG table** (for example `DEGs_A_with_outlier.csv`): each row is
one probe. The first, unnamed column is the probe's row number in the filtered
matrix. `SystematicName` is the PGSC transcript ID, `logFC` is the log2 fold
change (Low P relative to Reg P), and `adj.P.Val` is the FDR-adjusted p-value.

---

## 8. Glossary

| Term | Meaning in this project |
|---|---|
| **Microarray** | A glass slide carrying tens of thousands of short DNA probes. Labelled RNA from a sample sticks to matching probes, and the brightness of each spot measures how much of that transcript was present. |
| **Agilent one-colour (single-channel)** | Each array is hybridised with **one** sample labelled with a single dye (Cy3, "green"). Two-colour arrays instead put two samples (Cy3 + Cy5) on one slide and measure their ratio. This dataset is one-colour, hence `green.only = TRUE` in Stage 1. |
| **Probe** | One spot on the array, designed against one transcript. Several probes can target the same gene. |
| **Control probes** | Built-in spots that are not potato genes (`ControlType` ≠ 0). Negative controls measure background noise. |
| **SDRF / IDF** | ArrayExpress metadata files (MAGE-TAB format). The SDRF maps each data file to its sample attributes, and the IDF describes the experiment. |
| **Background correction (normexp)** | Removes the non-specific glow under each spot without producing negative values. |
| **Quantile normalisation** | Gives every array the same intensity distribution, so that differences between arrays reflect biology rather than labelling or scanning. |
| **log2 / log2FC** | Expression is analysed on a log2 scale. A log2 fold change (logFC) of 1 means 2× higher in Low P, −1 means 2× lower, and 3 means 8× higher. |
| **PCA** | Principal component analysis. It compresses thousands of genes into 2-3 axes, so similar samples sit close together. Used here for QC and to spot outliers. |
| **limma** | Bioconductor package that fits a linear model to every probe. Its **empirical Bayes** step (`eBayes`) borrows variance information across probes, which matters with only 3 replicates per group. |
| **Design matrix** | A table telling limma which group each array belongs to. `~0 + group` gives one column per group. |
| **Contrast** | The comparison being tested, e.g. `Low.P - Reg.P`. A positive logFC means higher under phosphate starvation. |
| **p-value / adj.P.Val / FDR / padj** | With about 39,000 probes tested, roughly 5% would pass p < 0.05 by chance. The Benjamini-Hochberg (BH) method adjusts the p-values to control the **false discovery rate** (FDR): among genes called significant at adj.P.Val < 0.05, about 5% are expected to be false positives. DESeq2 calls the same quantity "padj". |
| **DEG** | Differentially expressed gene. Here it means adj.P.Val < 0.05 **and** \|logFC\| > 1, i.e. at least a 2-fold change. |
| **Analysis A vs B** | A uses all 6 Maris Piper arrays. B drops one borderline Reg P array as a robustness check. |
| **GO (Gene Ontology)** | A controlled vocabulary of gene functions with three branches: biological process, molecular function and cellular component. |
| **KEGG** | A database of metabolic and signalling pathways. Potato PGSC IDs could not be mapped to it, so there is no KEGG enrichment (see Stage 3). |
| **Enrichment (over-representation) analysis** | Tests whether a function appears among the DEGs more often than expected, given the **background (universe)**: the genes that could have been detected. It uses a hypergeometric test (`clusterProfiler::enricher`). |
| **PGSC IDs** | Potato Genome Sequencing Consortium identifiers from the DM v3.4 genome (2011). `PGSC0003DMT…` is a transcript, `PGSC0003DMG…` a gene. |
| **BioMart / Ensembl Plants** | The online database that Stage 3 queries for GO terms and gene descriptions of the PGSC IDs. |
| **BLAST distance tree** | A quick neighbour-joining tree that NCBI builds from BLAST pairwise alignments. It is a useful picture of sequence similarity, but not a formal phylogenetic analysis. |

---

## 9. Troubleshooting

| Problem | Fix |
|---|---|
| `install.packages("limma")`: *package 'limma' is not available* | limma is on **Bioconductor**, not CRAN. Use `BiocManager::install("limma")`, or simply run `install_packages.R`. |
| *Bioconductor version cannot be validated; no internet connection?* | BiocManager cannot reach bioconductor.org because you are offline, or a firewall or proxy is blocking it. Check `curl -I https://bioconductor.org`, set any proxy (`Sys.setenv(https_proxy = "http://host:port")`), then retry. |
| *Bioconductor version 'x' requires R version 'y'* | Each Bioconductor release works with one R version. Update R, or install the matching release: `BiocManager::install(version = "3.18")` for R 4.3, `"3.20"` for R 4.4, `"3.21"` for R 4.5. |
| Compilation errors on Linux (e.g. `xml2`, `curl`, `openssl`, `ragg`) | Install the system libraries first. On Debian/Ubuntu: `sudo apt install libxml2-dev libcurl4-openssl-dev libssl-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev`. On Windows and macOS the default binary packages normally need no compiler. |
| *lib is not writable* / *installation paths not writeable* | Accept R's offer to create a personal library, or set `R_LIBS_USER`. Do not run R as administrator just to install packages. Warnings that *base/recommended* packages cannot be updated are harmless. |
| `pandoc version 1.12.3 or higher is required` | Render from RStudio, which bundles pandoc, or install pandoc system-wide. |
| Stage 1: *SDRF not found* or *raw file(s) … are missing* | Run `Rscript download_data.R`, or compare your folders with the tree in section 5. File names must match the SDRF exactly, with no sub-folders inside `data/raw/`. |
| Stage 2/3/4: *Missing input file(s): results/…* | Render the earlier stages first, in order. |
| Stage 3: BioMart timeout, *Unable to query the Ensembl site*, HTTP 500/503 | Ensembl Plants is busy or down. Wait and re-knit, and check <https://plants.ensembl.org>. A VPN or institutional firewall can also block it. |
| Stage 3: GO results differ slightly from the README | Ensembl Plants annotation changes between releases (see section 10, point 7). Small differences are expected. |
| Stage 5: species "NOT found in tree" | Your BLAST hits differ from the author's. Edit `selected_species` to match the species listed in Step 1 of that notebook. |

---

## 10. Known limitations and analysis caveats

The statistical methods are documented here, not changed. Reviewers and
anyone reusing the pipeline should keep these points in mind:

1. **Enrichment background.** The GO universe is every PGSC-annotated gene on
   the array. The more conservative choice is only the genes that passed the
   Stage 1 expression filter. The current choice can inflate enrichment of
   root-expressed functions.
2. **No GO term survives BH correction.** GO results are reported at nominal
   p < 0.05 and should be read as exploratory.
3. **Control probes are kept** through filtering and limma, so they count
   towards the multiple-testing correction. Stage 3 drops them before
   annotation; the Stage 4 tables do not filter them explicitly.
4. **Filtering order.** The negative-control noise threshold is computed
   *after* a log2 > 5 pre-filter, so it is estimated only from the brighter
   negative controls.
5. **Probe vs gene counts.** DEGs are counted per probe. Several probes can
   map to one PGSC transcript, so gene-level counts can be lower.
6. **Hard-coded numbers in the text and figures.** Parts of the narrative,
   and the functional-category bar charts (Stage 3 `cat_summary`, Stage 4
   `category_summary`), are typed in rather than computed. Some of these
   statements also disagree with each other. The "top DEG" is given as:
   - an ethylene-responsive gene with logFC 5.61 (Stage 3 summary);
   - the ABC transporter at about 79-fold (earlier README);
   - logFC about 6.5, i.e. about 91-fold (a Stage 4 comment).

   Check these against your own rendered output.
7. **Genome versions.** The probes carry PGSC v3.4 transcript IDs, but
   Ensembl Plants now serves the DM v6.1 assembly. The `DMT → DMG` swap only
   matches genes whose PGSC IDs are still present in Ensembl, which lowers
   annotation coverage.
8. **KEGG attempts.** The KEGG work described in Stage 3 was done outside the
   notebooks and cannot be reproduced from this repository.
9. **Stage 5 is descriptive.** The tree is an NCBI BLAST distance tree for one
   query, not a phylogeny built from a curated multiple alignment with an
   evolutionary model. The code does not test claims about conservation
   across Solanaceae or divergence ages.
10. **"First systematic analysis".** This claim rests on the author's
    literature search and has not been verified here.
11. **Small n.** There are 3 replicates per group (2 for Reg P in Analysis B).
    limma's empirical Bayes helps, but statistical power is limited.

---

## 11. Data provenance and citation

- **Dataset:** ArrayExpress / BioStudies accession **E-MTAB-2990**, "Gene
  Expression Patterns in Potato Roots Associated with Greater Plant Growth at
  Low Phosphate Concentrations in the Rhizosphere". Submitted by Pete Hedley,
  James Hutton Institute, UK.
  <https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-2990>.
  The notebooks cite it as 2014 (deposition) and as 2017 in different places,
  so use the dates on the study page when citing.
- **Original publication:** none could be identified from this repository.
  If the study page's *Publication* field lists one, cite it too.
- **Archive:** Sarkans U. *et al.* (2021) From ArrayExpress to BioStudies.
  *Nucleic Acids Research* 49(D1):D1502-D1506.
- **Software and methods:** Ritchie M.E. *et al.* (2015) limma, *NAR* 43:e47;
  Smyth G.K. (2004) *Stat Appl Genet Mol Biol* 3:Article 3; Benjamini Y. &
  Hochberg Y. (1995) *JRSS-B* 57:289-300; Wu T. *et al.* (2021)
  clusterProfiler 4.0, *The Innovation* 2:100141; Durinck S. *et al.* (2009)
  biomaRt, *Nature Protocols* 4:1184-1191; Paradis E. & Schliep K. (2019)
  ape 5.0, *Bioinformatics* 35:526-528; Potato Genome Sequencing Consortium
  (2011) *Nature* 475:189-195.
- **This repository:** Mohomed H.R. *Transcriptomic analysis of potato
  phosphate stress response (E-MTAB-2990)*. GitHub:
  halireena/E-MTAB-2990-Potato-Phosphate-Transcriptomics.

No raw data are redistributed here.

---

## 12. Acknowledgements, author and licence

Sincere thanks to **Dr. Reem** (University of Sharjah, Research Institute for
Medical and Health Sciences) for her supervision, guidance and support
throughout this internship project. Thanks also to lab colleagues Fathima
Nubla Latheef, Zoha Shaikh and Omar Mohammad for their encouragement and
support.

**Author:** Halireena Rushdiha Mohomed, MSc Bioinformatics, University of
Birmingham Dubai (internship at the University of Sharjah, supervised by
Dr. Reem).

**Licence:** code under the MIT licence ([`LICENSE`](LICENSE)). The E-MTAB-2990
data are covered by the EMBL-EBI terms of use.

*About the file names:* the `_FIXED`, `_FINAL` and `_CLEAN` suffixes are
historical. The notebooks are simply Stages 1-5, run in numeric order.
