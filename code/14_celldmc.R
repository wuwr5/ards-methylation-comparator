# ==============================================================
# 14_celldmc.R  (round-5 reviewer response)
#
# Reviewer point: the 6-proportion adjustment is an over-adjustment /
# collinearity risk, and the "99.97% attenuation" may be an artefact.
# Orthogonal evidence requested: CellDMC (cell-type-specific differential
# methylation), which models the interaction instead of removing the
# composition.
#
# Model (identical to EpiDISH::CellDMC, fitted vectorised with limma):
#     beta_k ~ pheno + sum_j frac_j + pheno:sum_j frac_j + age + male + black
# The interaction pheno:frac_j tests whether the ARDS-vs-ICU difference in
# CpG k differs across the fraction of cell type j, i.e. whether the effect
# is cell-type specific. Pure composition confounding concentrates the
# signal in the cell type that actually shifts (neutrophils); a genuine
# cell-intrinsic signal survives inside at least one cell type.
#
# Why limma and not EpiDISH::CellDMC for the whole array:
#   * CellDMC loops per CpG in R; on Windows mc.cores must be 1 and the
#     measured cost is ~7.3 ms/CpG -> ~60 min for 485 577 CpGs, and it
#     returns BH-adjusted DMC calls per chunk rather than p-values.
#   * limma fits exactly the same per-CpG linear model, vectorised, with
#     global BH -> seconds, and gives p-values and effect sizes.
#   * Equivalence is not assumed: set RUN_EPIDISH=1 to additionally run
#     EpiDISH::CellDMC on a random subset and compare coefficients.
#
# Run:  Rscript code/14_celldmc.R            (limma only, ~1-2 min)
#       RUN_EPIDISH=1 Rscript code/14_celldmc.R
# Out:  results/celldmc_summary.csv
#       results/celldmc_signature.csv
#       results/celldmc_validation.csv     (only with RUN_EPIDISH=1)
# ==============================================================

resolve_root <- function() {
  env <- Sys.getenv("ARDS_EPI_ROOT", "")
  if (nzchar(env)) return(normalizePath(env, mustWork = FALSE))
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  start <- if (length(f)) dirname(normalizePath(f, mustWork = FALSE)) else getwd()
  d <- normalizePath(start, mustWork = FALSE)
  while (!identical(d, dirname(d))) {
    if (file.exists(file.path(d, "R", "01_download.R")) ||
        file.exists(file.path(d, "code", "01_download.R"))) return(d)
    d <- dirname(d)
  }
  d
}
ROOT <- resolve_root()
ca <- commandArgs(FALSE)
cf <- sub("^--file=", "", ca[grep("^--file=", ca)])
REPO_DIR <- if (length(cf)) dirname(dirname(normalizePath(cf[1], mustWork = FALSE))) else ROOT
RES <- file.path(REPO_DIR, "results"); dir.create(RES, showWarnings = FALSE, recursive = TRUE)
suppressMessages(library(limma))

pick <- function(f) {
  x <- c(file.path(ROOT, "data", f), file.path(REPO_DIR, "data", f),
         file.path(dirname(ROOT), "ards_epigenetics", "data", f))
  x <- x[file.exists(x)]; if (!length(x)) stop("not found: ", f); x[1]
}

cat("== loading ==\n")
beta  <- readRDS(pick("processed/gse67530_beta.rds"))
pdata <- readRDS(pick("processed/gse67530_pdata.rds"))
cellf <- c(file.path(RES, "cell_estimates_GSE67530.csv"),
           file.path(dirname(REPO_DIR), "ards_epigenetics", "output", "cell_estimates_GSE67530.csv"))
cells <- read.csv(cellf[file.exists(cellf)][1], row.names = 1, check.names = FALSE)
if ("geo_accession" %in% names(cells)) { rownames(cells) <- cells$geo_accession; cells$geo_accession <- NULL }

if (!"geo_accession" %in% names(pdata)) pdata$geo_accession <- rownames(pdata)
pdata$geo_accession <- as.character(pdata$geo_accession)
pv <- function(x) trimws(sub("^[^:]*:\\s*", "", as.character(x)))
kv <- list()
for (cc in grep("characteristics", names(pdata), ignore.case = TRUE, value = TRUE)) {
  v <- as.character(pdata[[cc]]); k <- trimws(sub(":.*$", "", v[!is.na(v)][1]))
  if (nzchar(k)) kv[[k]] <- pv(v)
}
grp <- ifelse(kv$ards == "1", "ARDS", ifelse(kv$ards == "0", "ICU_control", "healthy"))
CASE  <- Sys.getenv("CASE", "ARDS")
CTRL  <- Sys.getenv("CTRL", "ICU_control")
TAG   <- Sys.getenv("TAG", if (CTRL == "ICU_control") "ards_icu" else "ards_healthy")
keep <- which(grp %in% c(CASE, CTRL))
n <- length(keep)
cat(sprintf("%s vs %s [%s]: n=%d (%s %d / %s %d)\n", CASE, CTRL, TAG, n,
            CASE, sum(grp[keep] == CASE), CTRL, sum(grp[keep] == CTRL)))

PARTS <- c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro")
# NOTE ON RANK: a model with pheno + all 6 fractions + all 6 pheno:fraction
# terms is rank deficient by 1 even before adding covariates, because
# sum_j frac_j = 1 makes sum_j pheno:frac_j = pheno. (This is why the
# original six-covariate adjustment was aliased, and why EpiDISH::CellDMC
# silently drops a term internally.) We therefore use an explicit reference
# cell type. With Neutro as reference, frac_Neutro = 1 - sum(other), so:
#   pheno            = ARDS-vs-ICU effect *within neutrophils*
#   pheno:frac_j     = effect in cell type j MINUS effect in neutrophils
# which is the interpretable, full-rank version of the CellDMC model and
# answers the key question directly: does the difference survive inside
# neutrophils (the cell type that actually shifts), or is it composition?
REF <- "Neutro"
FRAC_USE <- setdiff(PARTS, REF)
frac <- as.matrix(cells[pdata$geo_accession[keep], PARTS])
pheno <- as.integer(grp[keep] == CASE)
age <- as.numeric(kv$age[keep]); male <- as.numeric(kv$male[keep])
black <- if (!is.null(kv$black)) as.numeric(kv$black[keep]) else
  as.integer(tolower(kv$race[keep]) %in% c("black", "black or african american"))
ok <- complete.cases(frac, age, male, black)
frac <- frac[ok, , drop = FALSE]; pheno <- pheno[ok]
age <- age[ok]; male <- male[ok]; black <- black[ok]

X <- cbind(intercept = 1, pheno = pheno, frac[, FRAC_USE, drop = FALSE],
           age = age, male = male, black = black)
for (j in seq_along(FRAC_USE)) X <- cbind(X, pheno * frac[, FRAC_USE[j]])
colnames(X) <- c("intercept", "pheno", paste0("frac_", FRAC_USE), "age", "male", "black",
                 paste0("pheno:frac_", FRAC_USE))
cat("design:", ncol(X), "params, rank", qr(X)$rank, "| n =", nrow(X),
    "| reference cell type:", REF, "\n")
stopifnot(qr(X)$rank == ncol(X))

samples <- pdata$geo_accession[keep][ok]
B <- beta[, samples, drop = FALSE]
M <- B
invisible(gc())
for (blk in split(seq_len(nrow(M)), ceiling(seq_len(nrow(M)) / 20000))) {
  b <- M[blk, , drop = FALSE]
  b[!is.na(b) & b <= 0] <- 1e-4; b[!is.na(b) & b >= 1] <- 1 - 1e-4
  M[blk, ] <- log2(b / (1 - b))
}
rm(B); invisible(gc())

fit <- eBayes(lmFit(M, X))
cat("fit done\n")

SIG <- c("cg03341377", "cg24310395", "cg07830557", "cg08418670")
summ <- list(); sig <- list()
for (ct in c(REF, FRAC_USE)) {
  co <- if (ct == REF) "pheno" else paste0("pheno:frac_", ct)
  tt <- topTable(fit, coef = co, number = Inf, sort.by = "none")
  summ[[ct]] <- data.frame(celltype = ct, term = co,
                           n_fdr05 = sum(tt$adj.P.Val < 0.05, na.rm = TRUE),
                           n_p1e3 = sum(tt$P.Value < 1e-3, na.rm = TRUE),
                           min_p = min(tt$P.Value, na.rm = TRUE))
  idx <- match(SIG, rownames(tt))
  sig[[ct]] <- data.frame(cpg = SIG, celltype = ct, term = co,
                          coef = tt$logFC[idx], t = tt$t[idx],
                          p = tt$P.Value[idx], adj.p = tt$adj.P.Val[idx],
                          stringsAsFactors = FALSE)
}
S <- do.call(rbind, summ); write.csv(S, file.path(RES, sprintf("celldmc_summary_%s.csv", TAG)), row.names = FALSE)
print(S, digits = 3)
SG <- do.call(rbind, sig)
write.csv(SG, file.path(RES, sprintf("celldmc_signature_%s.csv", TAG)), row.names = FALSE)
cat("\n=== signature CpG interactions (limma, global BH) ===\n")
print(reshape(SG[, c("cpg", "celltype", "coef", "adj.p")], idvar = "cpg", timevar = "celltype",
              direction = "wide"), digits = 3)

# ---- optional: verify against EpiDISH::CellDMC on a random subset ----
if (identical(Sys.getenv("RUN_EPIDISH"), "1")) {
  suppressMessages(library(EpiDISH))
  set.seed(20261103)
  idx <- sort(unique(c(sample.int(nrow(M), 3000), match(SIG, rownames(M)))))
  idx <- idx[!is.na(idx)]
  sub <- M[idx, , drop = FALSE]
  subb <- 2^sub / (1 + 2^sub)                      # CellDMC expects beta scale
  cov.mod <- data.frame(age = age, male = male, black = black)
  t0 <- Sys.time()
  r <- CellDMC(subb, pheno, frac[, FRAC_USE, drop = FALSE],
               adjPMethod = "none", cov.mod = cov.mod, mc.cores = 1)
  cat(sprintf("CellDMC on %d CpGs (%d cell types): %.1f s\n", nrow(subb),
              length(FRAC_USE), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  cat("CellDMC output elements:", paste(names(r), collapse = ", "), "\n")
  if (!is.null(r$coe)) cat("coe columns:", paste(colnames(r$coe), collapse = " | "), "\n")
  if (!is.null(r$dmct)) {
    cat("dmct columns:", paste(colnames(r$dmct), collapse = " | "), "\n")
    hit <- intersect(SIG, rownames(r$dmct))
    dmc <- data.frame(cpg = hit, r$dmct[hit, , drop = FALSE], check.names = FALSE)
    write.csv(dmc, file.path(RES, "celldmc_validation.csv"), row.names = FALSE)
    cat("\n=== CellDMC DMC calls for the signature CpGs ===\n"); print(dmc)
  }
}
cat("\nDONE 14\n")
