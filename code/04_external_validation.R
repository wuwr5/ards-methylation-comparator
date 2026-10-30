# =====================================================================
# 04_external_validation.R
# External validation of the four-CpG "ARDS signature" in GSE199591
# (COVID-19, EPIC array, 123 Mild vs 64 Severe [ICU/death], 187 samples)
#
# Inputs (data/processed/):
#   GSE199591_COVID_Matrix_processed.txt.gz  -- normalized beta matrix
#   GSE199591_series_matrix.txt.gz           -- metadata (10 KB)
#
# Outputs (output/, prefix gse199591_):
#   metadata table, signature CpG stats, AUC/DeLong, deconvolution,
#   cell-adjusted signature stats
# =====================================================================

suppressMessages({
  library(limma)
  library(pROC)
  library(EpiDISH)
  library(data.table)
})

set.seed(42)
# ---- project root: env var ARDS_EPI_ROOT, else walk up from script dir ----
resolve_root <- function() {
  env <- Sys.getenv("ARDS_EPI_ROOT", "")
  if (nzchar(env)) return(normalizePath(env, mustWork = FALSE))
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  start <- if (length(f)) dirname(normalizePath(f, mustWork = FALSE)) else getwd()
  d <- normalizePath(start, mustWork = FALSE)
  while (!identical(d, dirname(d))) {
    if (file.exists(file.path(d, "R", "01_download.R")) || file.exists(file.path(d, "code", "01_download.R"))) return(d)
    d <- dirname(d)
  }
  d
}
# ---- project root: env var ARDS_EPI_ROOT, else walk up from script dir ----
resolve_root <- function() {
  env <- Sys.getenv("ARDS_EPI_ROOT", "")
  if (nzchar(env)) return(normalizePath(env, mustWork = FALSE))
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  start <- if (length(f)) dirname(normalizePath(f, mustWork = FALSE)) else getwd()
  d <- normalizePath(start, mustWork = FALSE)
  while (!identical(d, dirname(d))) {
    if (file.exists(file.path(d, "R", "01_download.R")) || file.exists(file.path(d, "code", "01_download.R"))) return(d)
    d <- dirname(d)
  }
  d
}
PROJ <- resolve_root()
DP    <- file.path(PROJ, "data/processed")
OUT   <- file.path(PROJ, "output")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

sig_cpgs <- c("cg03341377", "cg24310395", "cg07830557", "cg08418670")
log_msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")),
                             ..., "\n")

## ------------------------------------------------ 1. metadata
log_msg("parsing series matrix metadata ...")
sm <- gzfile(file.path(DP, "GSE199591_series_matrix.txt.gz"), "rt")
lines <- readLines(sm, n = 2000L)
close(sm)

grab <- function(tag) {
  ln <- grep(sprintf("!Series_%s\t", tag), lines, value = TRUE)
  if (length(ln)) trimws(gsub('"', '', strsplit(sub('^[^\t]+\t', '', ln[1]), '\t')[[1]])) else NULL
}
smp_ch  <- grep("!Sample_characteristics_ch1\t", lines, value = TRUE)
smp_geo <- gsub('"', '', strsplit(
  sub('^[^\t]+\t', '', grep("!Sample_geo_accession\t", lines, value = TRUE)[1]), '\t')[[1]])
smp_ttl <- gsub('"', '', strsplit(
  sub('^[^\t]+\t', '', grep("!Sample_title\t", lines, value = TRUE)[1]), '\t')[[1]])

pd <- data.frame(geo_accession = smp_geo, title = smp_ttl,
                 stringsAsFactors = FALSE)
for (cc in smp_ch) {
  vals <- gsub('"', '', strsplit(sub('^[^\t]+\t', '', cc), '\t')[[1]])
  key  <- sub(":.*$", "", vals[1])
  val  <- trimws(sub("^[^:]*:\\s*", "", vals))
  pd[[key]] <- val
}
# standardise keys
names(pd) <- make.names(names(pd))
ev_col <- grep("clinical|evolution|severe", names(pd), ignore.case = TRUE, value = TRUE)
stopifnot(length(ev_col) >= 1)
pd$severe <- as.integer(grepl("severe", pd[[ev_col[1]]], ignore.case = TRUE))
age_col <- grep("age", names(pd), ignore.case = TRUE, value = TRUE)
sex_col <- grep("sex", names(pd), ignore.case = TRUE, value = TRUE)
pd$age <- suppressWarnings(as.numeric(pd[[age_col[1]]]))
pd$male <- as.integer(grepl("^[Mm]$", pd[[sex_col[1]]]))
cat("samples:", nrow(pd), " severe:", sum(pd$severe == 1),
    " mild:", sum(pd$severe == 0), "\n")
cat("age non-missing:", sum(!is.na(pd$age)), " male:", sum(pd$male == 1), "\n")
fwrite(pd, file.path(OUT, "gse199591_pdata.csv"))

## ------------------------------------------------ 2. beta matrix
log_msg("streaming beta matrix (990 MB gz) ...")
beta_file <- list.files(DP, pattern = "Matrix_processed.*\\.txt\\.gz$",
                        full.names = TRUE)
stopifnot(length(beta_file) == 1)

need <- c(sig_cpgs)                      # rows to keep in full pass
first <- gzfile(beta_file[1], "rt"); hdr <- strsplit(readLines(first, n = 1L),
                                                     "\t")[[1]]; close(first)
cat("matrix columns:", length(hdr), "\n")

# read with data.table::fread; skip *_Detection_Pval columns
log_msg("fread on beta matrix (may take a few minutes) ...")
hdr_cn <- colnames(fread(beta_file[1], nrows = 0L, header = TRUE))
sel_cn <- c("ID_REF", setdiff(hdr_cn, c("ID_REF",
            grep("_Detection_Pval$", hdr_cn, value = TRUE))))
cat("total cols:", length(hdr_cn), " beta cols kept:", length(sel_cn) - 1, "\n")
bm <- fread(beta_file[1], header = TRUE, sep = "\t", showProgress = FALSE,
            stringsAsFactors = FALSE, data.table = TRUE, select = sel_cn)
setnames(bm, 1, "probe")
setkey(bm, probe)
log_msg("matrix dim: ", nrow(bm), " x ", ncol(bm) - 1)

# align columns to pdata by title/GSM
cn <- colnames(bm)[-1]
# try match: title pattern "DC1_01_AU_BIS01_B_MILD" or GSM ids
pd$key <- pd$title
map <- match(pd$key, cn)
if (any(is.na(map))) {                    # fallback: GSM accessions
  map <- match(pd$geo_accession, cn)
}
cat("column alignment: matched ", sum(!is.na(map)), "/", nrow(pd), "\n")
stopifnot(!any(is.na(map)))
bm2 <- as.matrix(bm[, map + 1L, with = FALSE])
rownames(bm2) <- bm$probe
colnames(bm2) <- pd$geo_accession
rm(bm); invisible(gc())

cat("beta range:", range(bm2), " n>1:", sum(bm2 > 1, na.rm = TRUE),
    " n<0:", sum(bm2 < 0, na.rm = TRUE), "\n")
bm2[bm2 < 0] <- 0; bm2[bm2 > 1] <- 1

sig_present <- sig_cpgs %in% rownames(bm2)
cat("signature CpGs present:", sum(sig_present), "/4:",
    paste(sig_cpgs[sig_present], collapse = ", "), "\n")
if (!all(sig_present)) {
  missing <- sig_cpgs[!sig_present]
  cat("MISSING (filtered by RnBeads pre-processing):",
      paste(missing, collapse = ", "), "\n")
}

## ------------------------------------------------ 3. EWAS severe vs mild
log_msg("limma Severe vs Mild (age+sex adjusted) ...")
keep <- !is.na(pd$severe)
bs   <- bm2[, pd$geo_accession[keep]]
pd2  <- pd[keep, ]
y    <- pd2$severe

X <- cbind(intercept = 1, severe = y)
if (sum(!is.na(pd2$age)) > 5 && length(unique(na.omit(pd2$age))) > 2)
  X <- cbind(X, age = ifelse(is.na(pd2$age), median(pd2$age, na.rm = TRUE), pd2$age))
if (length(unique(pd2$male)) > 1) X <- cbind(X, male = pd2$male)

M  <- qlogis(pmin(pmax(bs, 1e-4), 1 - 1e-4))
fit <- eBayes(lmFit(M, X))
tt  <- topTable(fit, coef = "severe", number = Inf, adjust.method = "BH")
tt$probe <- rownames(tt)
fwrite(tt, file.path(OUT, "gse199591_dmps_severe_vs_mild.csv"))
cat("DMPs FDR<0.05:", sum(tt$adj.P.Val < 0.05),
    " Bonferroni:", sum(tt$P.Value < 0.05/nrow(tt)), "\n")

## ------------------------------------------------ 4. signature CpG stats
log_msg("signature CpG table ...")
st <- data.frame()
for (g in intersect(sig_cpgs, rownames(bs))) {
  bt <- bs[g, ]; mb <- tapply(bt, y, mean)
  d  <- data.frame(probe = g,
                   beta_mild = mb["0"], beta_severe = mb["1"],
                   delta_beta = mb["1"] - mb["0"],
                   logFC = tt[g, "logFC"], adj_P = tt[g, "adj.P.Val"],
                   stringsAsFactors = FALSE)
  if (length(unique(y)) == 2) {
    roc <- roc(y, bt, quiet = TRUE)
    ci  <- as.numeric(ci.auc(roc, method = "delong"))
    d$auc <- as.numeric(auc(roc)); d$auc_lo <- ci[1]; d$auc_hi <- ci[3]
  }
  st <- rbind(st, d)
}
fwrite(st, file.path(OUT, "gse199591_signature_cpgs.csv"))
print(st, digits = 3)

## ------------------------------------------------ 5. deconvolution
log_msg("EpiDISH deconvolution (centDHSbloodDMC, RPC) ...")
data(centDHSbloodDMC.m)
common <- intersect(rownames(centDHSbloodDMC.m), rownames(bs))
cat("ref CpGs available on EPIC subset:", length(common), "/333\n")
ed <- epidish(bs[common, ], centDHSbloodDMC.m[common, ], method = "RPC")
frac <- ed$estF
frac <- cbind(data.frame(geo_accession = colnames(bs)), as.data.frame(frac))
fwrite(frac, file.path(OUT, "gse199591_cellprops.csv"))

## group means + tests
gm <- aggregate(. ~ grp,
                data = cbind(grp = ifelse(pd2$severe == 1, "Severe", "Mild"),
                             frac[, -1]),
                FUN = mean)
fwrite(gm, file.path(OUT, "gse199591_cellprops_group_means.csv"))
tests <- do.call(rbind, lapply(colnames(frac)[-1], function(ct) {
  wt <- wilcox.test(frac[, ct] ~ pd2$severe)
  data.frame(cell_type = ct, mild_mean = mean(frac[pd2$severe == 0, ct]),
             severe_mean = mean(frac[pd2$severe == 1, ct]),
             wilcox_p = wt$p.value)
}))
tests$p_bh <- p.adjust(tests$wilcox_p, "BH")
fwrite(tests, file.path(OUT, "gse199591_cellprops_group_tests.csv"))
print(tests, digits = 3)

## ------------------------------------------------ 6. cell-adjusted signature
log_msg("cell-adjusted re-modelling ...")
cells <- as.data.frame(scale(frac[, -1]))
# drop collinear component (Eosino ~ zero)
cells <- cells[, colMeans(abs(frac[, -1])) > 1e-4, drop = FALSE]
Xa <- cbind(X, cells[keep_seq <- seq_len(nrow(pd2)), , drop = FALSE])
fit2 <- eBayes(lmFit(M, as.matrix(Xa)))
tt2 <- topTable(fit2, coef = "severe", number = Inf, adjust.method = "BH")
tt2$probe <- rownames(tt2)
fwrite(tt2, file.path(OUT, "gse199591_dmps_severe_vs_mild_celladj.csv"))
cat("cell-adjusted DMPs FDR<0.05:", sum(tt2$adj.P.Val < 0.05), "\n")

st2 <- data.frame()
for (g in intersect(sig_cpgs, rownames(M))) {
  st2 <- rbind(st2, data.frame(probe = g,
                               delta_beta = mean(bs[g, y == 1]) - mean(bs[g, y == 0]),
                               adj_P_raw = tt[g, "adj.P.Val"],
                               adj_P_celladj = tt2[g, "adj.P.Val"]))
}
fwrite(st2, file.path(OUT, "gse199591_signature_celladj.csv"))
print(st2, digits = 3)

## ------------------------------------------------ 7. figures
log_msg("figures ...")
png(file.path(OUT, "fig_gse199591_cellprops.png"), width = 2400,
    height = 1800, res = 240)
par(mfrow = c(2, 4), mar = c(4, 4, 2.5, 1))
for (ct in colnames(frac)[-1]) {
  boxplot(frac[, ct] ~ factor(ifelse(pd2$severe == 1, "Severe", "Mild"),
                              levels = c("Mild", "Severe")),
          col = c("#9ecae1", "#fb6a4a"), main = ct, xlab = "", ylab = "proportion",
          las = 1, cex.main = 0.9, cex.lab = 0.8, cex.axis = 0.7)
}
dev.off()

png(file.path(OUT, "fig_gse199591_volcano.png"), width = 2000,
    height = 1800, res = 240)
plot(tt$logFC, -log10(tt$P.Value), pch = 16, cex = .35, col = "#888888",
     xlab = "log2 FC (Severe vs Mild)", ylab = "-log10 P",
     main = "GSE199591 EWAS: Severe vs Mild (EPIC)")
sig <- tt$adj.P.Val < 0.05
points(tt$logFC[sig], -log10(tt$P.Value)[sig], pch = 16, cex = .4,
       col = "#cb181d")
abline(h = -log10(0.05/nrow(tt)), lty = 2, col = "#333333")
if (all(sig_cpgs %in% rownames(tt))) {
  points(tt[sig_cpgs, "logFC"], -log10(tt[sig_cpgs, "P.Value"]),
         pch = 17, cex = 1.1, col = "#08519c")
  text(tt[sig_cpgs, "logFC"], -log10(tt[sig_cpgs, "P.Value"]),
       labels = sig_cpgs, pos = 4, cex = .55)
}
dev.off()

log_msg("DONE. Outputs in ", OUT)
