# download_data.R
# ---------------------------------------------------------------------------
# Downloads the E-MTAB-2990 files needed by Stages 1-3 and puts them where the
# notebooks expect them:
#
#   metadata/E-MTAB-2990.sdrf.txt   sample sheet (which file = which sample)
#   metadata/E-MTAB-2990.idf.txt    experiment description (optional, for reading)
#   data/raw/*.txt                  the 24 Agilent Feature Extraction files
#
# Usage (from the repository root):
#   Rscript download_data.R
#
# Safe to re-run: files that are already present are not downloaded again.
# Needs only base R (utils). Downloads several hundred MB; allow a few minutes.
#
# If this script fails (EBI changes its URLs from time to time), download the
# same files by hand from
#   https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-2990
# (Files tab) and put them in the folders above - see README -> Get the data.
# ---------------------------------------------------------------------------

options(timeout = max(3600, getOption("timeout")))  # large zips on slow links

acc      <- "E-MTAB-2990"
# Override with the environment variable EMTAB_BASE_URL if EBI moves the files
base_url <- Sys.getenv("EMTAB_BASE_URL",
  "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/990/E-MTAB-2990/Files/")
meta_dir <- "metadata"
raw_dir  <- file.path("data", "raw")
zip_dir  <- file.path("data", "zips")

if (!file.exists("01_load_QC_FIXED.Rmd")) {
  stop("Run this script from the repository root (the folder with the .Rmd files).")
}
for (d in c(meta_dir, raw_dir, zip_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

fetch <- function(file, dest_dir, required = TRUE) {
  dest <- file.path(dest_dir, file)
  if (file.exists(dest) && file.size(dest) > 0) {
    message("  already present: ", dest)
    return(TRUE)
  }
  message("  downloading ", file, " ...")
  ok <- tryCatch({
    utils::download.file(paste0(base_url, file), dest, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok) {
    unlink(dest)
    if (required) stop("Could not download ", paste0(base_url, file),
                       "\nDownload it manually (see README -> Get the data).")
  }
  ok
}

# 1. Metadata ---------------------------------------------------------------
message("Step 1/3: metadata")
invisible(fetch(paste0(acc, ".sdrf.txt"), meta_dir))
invisible(fetch(paste0(acc, ".idf.txt"),  meta_dir, required = FALSE))

sdrf <- read.delim(file.path(meta_dir, paste0(acc, ".sdrf.txt")), check.names = FALSE)
wanted <- unique(sdrf[["Array Data File"]])
message("  SDRF lists ", length(wanted), " raw data files")

# 2. Raw data archives ------------------------------------------------------
message("Step 2/3: raw data")
have_all <- function() all(file.exists(file.path(raw_dir, wanted)))

if (!have_all()) {
  # The SDRF usually names the archive(s) holding each raw file
  zip_col <- grep("ArrayExpress FTP file", names(sdrf), value = TRUE)
  zips <- if (length(zip_col)) unique(basename(sdrf[[zip_col[1]]])) else character(0)
  zips <- zips[nzchar(zips) & !is.na(zips)]
  # Otherwise try the standard ArrayExpress archive names raw.1.zip, raw.2.zip, ...
  if (length(zips) == 0) zips <- paste0(acc, ".raw.", 1:10, ".zip")

  for (z in zips) {
    if (have_all()) break
    got <- fetch(z, zip_dir, required = FALSE)
    if (!got) next
    message("  unzipping ", z)
    utils::unzip(file.path(zip_dir, z), exdir = raw_dir, junkpaths = TRUE)
  }
}

# 3. Check ------------------------------------------------------------------
message("Step 3/3: checking")
missing <- wanted[!file.exists(file.path(raw_dir, wanted))]
if (length(missing) > 0) {
  stop(length(missing), " of ", length(wanted), " raw files are still missing from ",
       raw_dir, ", e.g. ", paste(head(missing, 3), collapse = ", "),
       "\nDownload the raw archive(s) by hand (README -> Get the data).")
}
message("All ", length(wanted), " raw files present in ", raw_dir, ".")
message("You can delete ", zip_dir, " to save space. Next: render 01_load_QC_FIXED.Rmd")
