# ==============================================================
# 00_install.R
# One-time dependency installation for the ARDS epigenetics project.
# Run this BEFORE the analysis pipeline. Requires R >= 4.3 with
# BiocManager available (install.packages("BiocManager") if missing).
# ==============================================================

install_pkgs <- function(pkgs, from = c("bioc", "cran")) {
  from <- match.arg(from, several.ok = TRUE)
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      cat(sprintf(">> Installing %s ...\n", p))
      if ("bioc" %in% from) {
        tryCatch(
          BiocManager::install(p, ask = FALSE, update = FALSE),
          error = function(e) cat(sprintf("  [bioc] failed for %s: %s\n", p, e$message))
        )
      } else {
        install.packages(p, repos = "https://cloud.r-project.org")
      }
    } else {
      cat(sprintf("   %s already present.\n", p))
    }
  }
}

## CRAN packages
cran_pkgs <- c(
  "tidyverse",   # data wrangling / plotting
  "data.table",  # fast tables
  "survival",    # Cox models
  "survminer",   # survival plots
  "glmnet",      # penalised regression
  "pROC",        # ROC / AUC
  "cowplot",     # combined figures
  "here"         # portable paths
)

## Bioconductor packages
bioc_pkgs <- c(
  "BiocManager",
  "minfi",                  # IDAT -> beta (if processing raw)
  "ChAMP",                  # champ.load pipeline
  # NOTE: GEOquery deliberately NOT listed. 01_download.R uses base R
  # (download.file + gzfile) to fetch the GEO series matrix directly, so
  # GEOquery is unnecessary. On this machine its install is broken anyway
  # (library(GEOquery) segfaults). Install it only if you want the
  # optional fallback path in 01_download.R:
  #   BiocManager::install("GEOquery")
  "limma",                  # EWAS linear models
  "DSS",                    # smoothing / DMR
  "methylKit",              # alternative DMP
  "EpiDISH",                # cell-type deconvolution (RPC/Houseman)
  "meffil",                 # alternative deconvolution
  "wateRmelon",             # normalisation / QC
  "sva",                    # ComBat batch correction
  "ENmix",                  # additional QC
  "FlowSorted.Blood.450k",  # Houseman reference (450K)
  "FlowSorted.Blood.EPIC",  # Houseman reference (EPIC)
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "annotatr"                # genomic annotation
)

cat("=== Installing CRAN packages ===\n")
install_pkgs(cran_pkgs, from = "cran")

cat("=== Installing Bioconductor packages ===\n")
install_pkgs(bioc_pkgs, from = "bioc")

cat("\n=== Done. Verify with: \n")
cat('  sapply(c("GEOquery","limma","EpiDISH","FlowSorted.Blood.450k"), requireNamespace, quietly=TRUE)\n')
