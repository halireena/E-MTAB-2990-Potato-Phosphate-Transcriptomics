# install_packages.R
# ---------------------------------------------------------------------------
# Installs every package used by the five Stage notebooks (01-05).
# Derived from the library() / requireNamespace() calls in the .Rmd files.
#
# Usage (from the repository root):
#   Rscript install_packages.R
# or, inside R:
#   source("install_packages.R")
#
# Only missing packages are installed; already-installed ones are left alone.
# ---------------------------------------------------------------------------

options(repos = c(CRAN = "https://cloud.r-project.org"))

# CRAN packages ------------------------------------------------------------
cran_pkgs <- c(
  "rmarkdown",      # rendering the notebooks
  "knitr",          # chunk engine used by rmarkdown
  "ggplot2",        # Stages 1-4: plotting
  "dplyr",          # Stages 2-4: data manipulation
  "RColorBrewer",   # Stages 1, 2, 4: colour palettes
  "pheatmap",       # Stages 1, 2, 4: heatmaps
  "ggrepel",        # Stage 4: non-overlapping volcano labels
  "ggVennDiagram",  # Stages 3, 4: Venn diagrams
  "ape",            # Stage 5: phylogenetic tree handling/plotting
  "BiocManager"     # installer for the Bioconductor packages below
)

# Bioconductor packages ----------------------------------------------------
bioc_pkgs <- c(
  "limma",            # Stages 1, 2, 4: Agilent import, normalisation, DE
  "clusterProfiler",  # Stage 3: GO enrichment (enricher())
  "biomaRt"           # Stage 3: Ensembl Plants annotation download
)

missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1),
                                  quietly = TRUE)]
if (length(missing_cran) > 0) {
  message("Installing CRAN packages: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran)
}

missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1),
                                  quietly = TRUE)]
if (length(missing_bioc) > 0) {
  message("Installing Bioconductor packages: ",
          paste(missing_bioc, collapse = ", "))
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

# Report -------------------------------------------------------------------
all_pkgs <- c(cran_pkgs, bioc_pkgs)
ok <- vapply(all_pkgs, requireNamespace, logical(1), quietly = TRUE)
for (p in all_pkgs) {
  cat(sprintf("  %-16s %s\n", p,
              if (ok[[p]]) as.character(packageVersion(p)) else "NOT INSTALLED"))
}
if (!all(ok)) {
  stop("Some packages failed to install: ",
       paste(all_pkgs[!ok], collapse = ", "),
       "\nSee README.md -> Troubleshooting for common Bioconductor install errors.")
}
cat("All packages available.\n")
