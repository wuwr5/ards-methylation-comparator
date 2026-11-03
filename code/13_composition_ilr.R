# ==============================================================
# 13_composition_ilr.R  (round-5 reviewer response)
#
# Question addressed: is the ~99.97% attenuation of the EWAS signal by
# "cell composition adjustment" an artefact of over-adjustment /
# collinearity, rather than genuine confounding?
#
# Design of the original analysis: six EpiDISH proportions entered as six
# separate covariates. But the six proportions sum to 1 by construction
# (Eosino is estimated as exactly 0 for every sample), so the design
# [1 | p1..p6] is rank deficient: condition number > 1e15, VIF = Inf.
# limma silently aliases one column, so the reported "adjusted" model was
# fitted on an arbitrary, non-reproducible 5-of-6 parameterisation.
#
# Fix (compositional data analysis, CoDA):
#   * drop Eosino (identically 0 -> carries no information)
#   * multiplicative zero replacement (0.65 x smallest positive value per
#     part), then closure, then isometric log-ratio (ilr) transform with
#     the Helmert orthonormal basis -> 5 full-rank, unconstrained
#     components, no aliasing
#   * CLR of the same closed composition is fitted as well; CLR (minus one
#     redundant component, since sum(clr) = 0) is a LINEAR reparameterisation
#     of ilr, so the two fits must give identical group coefficients.
#     Agreement is asserted numerically, not assumed.
#
# Three models are compared for each contrast:
#   unadj   : group + age + male + black
#   rawprop : + six raw proportions            (the original, aliased fit)
#   ilr     : + ilr components (5)             (the fix)
#   clr     : + clr components minus one (5)   (equivalence check)
#
# Outputs (repo results/):
#   composition_dmp_counts.csv
#   composition_signature_cpgs.csv
#   composition_model_agreement.csv
#   composition_ilr_components.csv
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

datadir <- function(f) {
  c(file.path(ROOT, "data", f),
    file.path(REPO_DIR, "data", f),
    file.path(dirname(ROOT), "ards_epigenetics", "data", f),
    file.path(dirname(REPO_DIR), "ards_epigenetics", "data", f))
}
pick <- function(f) { x <- datadir(f); x <- x[file.exists(x)]; if (!length(x)) stop("not found: ", f); x[1] }
outdirs <- unique(c(file.path(REPO_DIR, "results"),
                    file.path(ROOT, "output"),
                    file.path(dirname(ROOT), "ards_epigenetics", "output")))
outdirs <- outdirs[dir.exists(outdirs) | dir.exists(dirname(outdirs))]
RES <- file.path(REPO_DIR, "results")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

suppressMessages({ library(limma) })
set.seed(1)

CONFIG_COV <- c("age", "male", "black")
SIG <- c("cg03341377", "cg24310395", "cg07830557", "cg08418670")
PARTS <- c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro")

cat("ROOT     :", ROOT, "\n")
cat("REPO_DIR :", REPO_DIR, "\n")

# ---------------- data ----------------
beta  <- readRDS(pick("processed/gse67530_beta.rds"))
pdata <- readRDS(pick("processed/gse67530_pdata.rds"))
cellf <- c(file.path(RES, "cell_estimates_GSE67530.csv"),
           file.path(dirname(REPO_DIR), "ards_epigenetics", "output", "cell_estimates_GSE67530.csv"))
cellf <- cellf[file.exists(cellf)][1]
cells <- read.csv(cellf, row.names = 1, check.names = FALSE)
if ("geo_accession" %in% names(cells)) { rownames(cells) <- cells$geo_accession; cells$geo_accession <- NULL }
cat("beta:", dim(beta), "| cell est:", dim(cells), "\n")

if (!"geo_accession" %in% names(pdata)) pdata$geo_accession <- rownames(pdata)
pdata$geo_accession <- as.character(pdata$geo_accession)
pv <- function(x) trimws(sub("^[^:]*:\\s*", "", as.character(x)))
# pdata stores characteristics as "key: value" with keys only in the VALUES,
# not in the column names (columns are characteristics_ch1, _2, _3, _4).
chcols <- grep("characteristics", names(pdata), ignore.case = TRUE, value = TRUE)
kv <- list()
for (cc in chcols) {
  v <- as.character(pdata[[cc]])
  k <- trimws(sub(":.*$", "", v[!is.na(v)][1]))
  if (!nzchar(k) || is.na(k)) next
  kv[[k]] <- pv(v)
}
cat("characteristic keys:", paste(names(kv), collapse = ", "), "\n")
ards_raw <- trimws(kv$ards)
pdata$group <- factor(ifelse(ards_raw == "1", "ARDS",
                      ifelse(ards_raw == "0", "ICU_control",
                      ifelse(ards_raw %in% c("NA", "", "na"), "healthy", NA_character_))),
                      levels = c("healthy", "ICU_control", "ARDS"))
for (cv in CONFIG_COV) {
  if (!cv %in% names(kv)) { cat("  covariate missing:", cv, "\n"); next }
  pdata[[cv]] <- suppressWarnings(as.numeric(kv[[cv]]))
  cat(sprintf("  %s: n=%d range %.1f-%.1f\n", cv, sum(!is.na(pdata[[cv]])),
              min(pdata[[cv]], na.rm = TRUE), max(pdata[[cv]], na.rm = TRUE)))
}
print(table(pdata$group))
stopifnot(all(pdata$geo_accession %in% colnames(beta)))
ord <- match(pdata$geo_accession, colnames(beta))
cells <- cells[rownames(cells) %in% pdata$geo_accession, , drop = FALSE]
cells <- cells[pdata$geo_accession[pdata$geo_accession %in% rownames(cells)], , drop = FALSE]
stopifnot(identical(rownames(cells), as.character(pdata$geo_accession)))

X <- as.matrix(cells[, PARTS])
cat("sum of 6 parts: mean", round(mean(rowSums(X)), 6), "min", round(min(rowSums(X)), 6), "\n")

# ---------------- zero replacement + closure + ilr ----------------
zero_replace <- function(x, delta_factor = 0.65) {
  # multiplicative replacement, per part, then closure
  out <- x
  for (j in seq_len(ncol(x))) {
    v <- x[, j]
    pos <- v[v > 0]
    d <- if (length(pos)) delta_factor * min(pos) else 1e-6
    if (any(v == 0)) out[v == 0, j] <- d
  }
  out / rowSums(out)
}
close01 <- function(x) x / rowSums(x)

helmert_basis <- function(D) {
  V <- matrix(0, D, D - 1)
  for (j in 1:(D - 1)) {
    V[1:j, j] <- 1 / sqrt(j * (j + 1))
    V[(j + 1):D, j] <- -sqrt(j / (j + 1))
  }
  V
}

ilr_trans <- function(xc) {           # xc: closed composition, samples x D
  lg <- log(xc)
  clr <- lg - rowMeans(lg)
  clr %*% helmert_basis(ncol(xc))
}
clr_trans <- function(xc) { lg <- log(xc); lg - rowMeans(lg) }

Z <- zero_replace(X)
ILR <- ilr_trans(Z)
CLR <- clr_trans(Z)
colnames(ILR) <- paste0("ilr", 1:ncol(ILR))
colnames(CLR) <- paste0("clr_", PARTS)
cat("ilr dims:", dim(ILR), "| max |sum(clr)|:", max(abs(rowSums(CLR))), "\n")
cat("zero count per part:", apply(X == 0, 2, sum), "\n")
write.csv(data.frame(geo_accession = rownames(cells), cells[, PARTS], ILR, check.names = FALSE),
          file.path(RES, "composition_ilr_components.csv"), row.names = FALSE)

# ---------------- M values, in place (memory) ----------------
M <- beta
invisible(gc())
for (blk in split(seq_len(nrow(M)), ceiling(seq_len(nrow(M)) / 20000))) {
  b <- M[blk, , drop = FALSE]
  b[!is.na(b) & b <= 0]      <- 1e-4
  b[!is.na(b) & b >= 1]      <- 1 - 1e-4
  M[blk, ] <- log2(b / (1 - b))
}
rm(beta); invisible(gc())
cat("M ready:", dim(M), "\n")

# ---------------- model fitting ----------------
fit_one <- function(case, ctrl, extra, tag) {
  keep <- which(pdata$group %in% c(case, ctrl))
  pd2 <- pdata[keep, ]
  y <- as.integer(pd2$group == case)
  Xd <- cbind(intercept = 1, case_vs_ctrl = y)
  used <- character(0)
  for (cv in CONFIG_COV) {
    if (cv %in% names(pd2) && sum(!is.na(pd2[[cv]])) > 5) {
      Xd <- cbind(Xd, pd2[[cv]]); used <- c(used, cv)
    }
  }
  if (!is.null(extra)) { Xd <- cbind(Xd, extra[keep, , drop = FALSE]); used <- c(used, colnames(extra)) }
  colnames(Xd) <- c("intercept", "case_vs_ctrl", used)
  r <- qr(Xd)$rank
  msg <- if (r < ncol(Xd)) sprintf("RANK DEFICIENT %d/%d", r, ncol(Xd)) else "full rank"
  cat(sprintf("  [%s] %s vs %s n=%d design=%d rank=%d (%s)\n",
              tag, case, ctrl, length(keep), ncol(Xd), r, msg))
  Ms <- M[, keep, drop = FALSE]
  fit <- eBayes(lmFit(Ms, Xd))
  tt <- topTable(fit, coef = "case_vs_ctrl", number = Inf, sort.by = "none")
  out <- data.frame(probe = rownames(tt), logFC = tt$logFC, t = tt$t,
                    P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val,
                    stringsAsFactors = FALSE)
  # M-scale group difference for the signature CpGs (model-free); NOT a beta
  # difference -- beta-scale delta-beta is produced by 15_tost_prespecified.R
  sig <- out[out$probe %in% SIG, ]
  sig$delta_Mvalue <- sapply(sig$probe, function(p) {
    v <- as.numeric(M[match(p, rownames(M)), keep]); mean(v[y == 1]) - mean(v[y == 0])
  })
  attr(out, "sig") <- sig
  attr(out, "n") <- c(case = sum(y), ctrl = sum(!y))
  out
}

res <- list()
res[["ards_vs_icu_unadj"]]   <- fit_one("ARDS", "ICU_control", NULL, "unadj")
res[["ards_vs_icu_rawprop"]] <- fit_one("ARDS", "ICU_control", X, "rawprop")
res[["ards_vs_icu_ilr"]]     <- fit_one("ARDS", "ICU_control", ILR, "ilr")
res[["ards_vs_icu_clr"]]     <- fit_one("ARDS", "ICU_control", CLR[, -1, drop = FALSE], "clr")
res[["ards_vs_hlth_unadj"]]  <- fit_one("ARDS", "healthy", NULL, "unadj")
res[["ards_vs_hlth_rawprop"]]<- fit_one("ARDS", "healthy", X, "rawprop")
res[["ards_vs_hlth_ilr"]]    <- fit_one("ARDS", "healthy", ILR, "ilr")

# ---------------- summaries ----------------
cnt <- do.call(rbind, lapply(names(res), function(k) {
  d <- res[[k]]
  data.frame(model = k, n_probe = nrow(d),
             n_fdr05 = sum(d$adj.P.Val < 0.05, na.rm = TRUE),
             n_bonf = sum(d$P.Value < 0.05 / nrow(d), na.rm = TRUE),
             min_adjp = min(d$adj.P.Val, na.rm = TRUE),
             stringsAsFactors = FALSE)
}))
write.csv(cnt, file.path(RES, "composition_dmp_counts.csv"), row.names = FALSE)
print(cnt)

sigtab <- do.call(rbind, lapply(names(res), function(k) {
  s <- attr(res[[k]], "sig")
  if (is.null(s) || !nrow(s)) return(NULL)
  s$model <- k
  s[, c("model", "probe", "logFC", "t", "P.Value", "adj.P.Val", "delta_Mvalue")]
}))
write.csv(sigtab, file.path(RES, "composition_signature_cpgs.csv"), row.names = FALSE)

# attenuation relative to unadjusted, in |logFC| terms
att <- do.call(rbind, lapply(c("ards_vs_icu", "ards_vs_hlth"), function(cx) {
  u <- res[[paste0(cx, "_unadj")]]
  for (mv in c("rawprop", "ilr")) {
    m <- res[[paste0(cx, "_", mv)]]
    if (is.null(m)) next
    d <- merge(u[, c("probe", "logFC")], m[, c("probe", "logFC")], by = "probe", suffixes = c("_unadj", "_adj"))
    cat(sprintf("%s %-8s: median |logFC| unadj %.4f -> adj %.4f | attenuation %.2f%%\n",
                cx, mv, median(abs(d$logFC_unadj)), median(abs(d$logFC_adj)),
                100 * (1 - median(abs(d$logFC_adj)) / median(abs(d$logFC_unadj)))))
  }
  NULL
}))

# model agreement: correlation of group t-statistics
agree <- data.frame()
for (cx in c("ards_vs_icu", "ards_vs_hlth")) {
  u <- res[[paste0(cx, "_unadj")]]
  for (mv in c("rawprop", "ilr", "clr")) {
    m <- res[[paste0(cx, "_", mv)]]
    if (is.null(m)) next
    j <- merge(u[, c("probe", "t", "P.Value")], m[, c("probe", "t", "P.Value")], by = "probe", suffixes = c("_u", "_m"))
    ok <- is.finite(j$t_u) & is.finite(j$t_m)
    agree <- rbind(agree, data.frame(
      contrast = cx, model = mv,
      r_t = cor(j$t_u[ok], j$t_m[ok]),
      spearman_t = cor(j$t_u[ok], j$t_m[ok], method = "spearman"),
      median_abs_t_unadj = median(abs(j$t_u[ok])),
      median_abs_t_adj = median(abs(j$t_m[ok])),
      n_fdr_u = sum(j$P.Value_u < 0.05/nrow(j)), n_fdr_m = sum(j$P.Value_m < 0.05/nrow(j)),
      stringsAsFactors = FALSE))
  }
}
# CLR vs ILR must be identical (linear reparameterisation)
if (!is.null(res$ards_vs_icu_clr)) {
  a <- res$ards_vs_icu_ilr; b <- res$ards_vs_icu_clr
  j <- merge(a[, c("probe", "logFC")], b[, c("probe", "logFC")], by = "probe", suffixes = c("_ilr", "_clr"))
  cat(sprintf("CLR vs ILR: max |diff in group logFC| = %.3e\n", max(abs(j$logFC_ilr - j$logFC_clr))))
  agree$clr_ilr_maxdiff <- NA
  agree$clr_ilr_maxdiff[1] <- max(abs(j$logFC_ilr - j$logFC_clr))
}
write.csv(agree, file.path(RES, "composition_model_agreement.csv"), row.names = FALSE)
print(agree)
cat("\nDONE 13\n")
