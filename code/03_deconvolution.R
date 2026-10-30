# ==============================================================
# 03_deconvolution.R  (v2, 2026-10-26)
# Estimate blood cell-type proportions in GSE67530 and quantify
# how much of the EWAS signal is explained by leukocyte shifts.
#
# Reference strategy (NO FlowSorted.* required):
#   450K data -> EpiDISH built-in centDHSbloodDMC.m
#                (DHS-based 450K reference, 7 cell types:
#                 B, NK, CD4T, CD8T, Mono, Neutro, Eosino)
#   This is the reference EpiDISH itself recommends for 450K.
#   EPIC data (GSE199591) would need cent12CT.m; skipped until
#   that dataset is downloaded.
#
# Method: EpiDISH RPC (robust partial correlations).
#
# Outputs (output/):
#   cell_estimates_GSE67530.csv          samples x 7 proportions
#   deconv_group_tests_GSE67530.csv      KW + pairwise Wilcoxon
#   fig_cellprops_GSE67530.png           boxplots by group
#
# After this script, re-run R/02_ewas.R: it picks up
# cell_estimates_GSE67530.csv automatically and re-fits all three
# contrasts with cell proportions as covariates (*_celladj tags).
# ==============================================================

# ---- project root: env var ARDS_EPI_ROOT, else walk up from script dir ----
resolve_root <- function() {
  env <- Sys.getenv("ARDS_EPI_ROOT", "")
  if (nzchar(env)) return(normalizePath(env, mustWork = FALSE))
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  start <- if (length(f)) dirname(normalizePath(f, mustWork = FALSE)) else getwd()
  d <- normalizePath(start, mustWork = FALSE)
  while (!identical(d, dirname(d))) {
    if (file.exists(file.path(d, "R", "01_download.R"))) return(d)
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
    if (file.exists(file.path(d, "R", "01_download.R"))) return(d)
    d <- dirname(d)
  }
  d
}
BASE <- resolve_root()
OUT  <- file.path(BASE, "output")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

suppressMessages(library(EpiDISH))

# --------------------------------------------------------------
# 1. load data
# --------------------------------------------------------------
cat("== loading beta matrix (528 MB, be patient) ==\n")
beta  <- readRDS(file.path(BASE, "data", "processed", "gse67530_beta.rds"))
pdata <- readRDS(file.path(BASE, "data", "processed", "gse67530_pdata.rds"))

if (!"geo_accession" %in% names(pdata)) pdata$geo_accession <- rownames(pdata)
pdata$geo_accession <- as.character(pdata$geo_accession)
stopifnot(all(pdata$geo_accession %in% colnames(beta)))
beta <- beta[, pdata$geo_accession, drop = FALSE]

# ---- parse "key: value" characteristics (same convention as 02) ----
kv_parser <- function(pd) {
  ch <- grep("characteristics", names(pd), ignore.case = TRUE, value = TRUE)
  out <- data.frame(geo_accession = pd$geo_accession, stringsAsFactors = FALSE)
  for (cc in ch) {
    v   <- as.character(pd[[cc]])
    k   <- trimws(sub(":.*$", "", v))
    val <- trimws(sub("^[^:]*:\\s*", "", v))
    for (kk in unique(k)) {
      col <- ifelse(k == kk, val, NA_character_)
      if (!kk %in% names(out)) out[[kk]] <- col
      else out[[kk]] <- ifelse(is.na(col), out[[kk]], col)
    }
  }
  out
}
kv <- kv_parser(pdata)
cat("   characteristic keys:", paste(setdiff(names(kv), "geo_accession"),
                                      collapse = ", "), "\n")

ards_raw <- trimws(kv$ards)
pdata$group <- factor(
  ifelse(ards_raw == "1", "ARDS",
  ifelse(ards_raw == "0", "ICU_control",
  ifelse(ards_raw %in% c("NA", "", "na"), "healthy", NA_character_))),
  levels = c("healthy", "ICU_control", "ARDS"))
cat("   group counts:\n"); print(table(pdata$group, useNA = "ifany"))
for (cv in c("age", "male", "black")) {
  if (cv %in% names(kv))
    pdata[[cv]] <- suppressWarnings(as.numeric(kv[[cv]]))
}

# ---- clamp beta (same scale rule as 02; deposited values are beta) ----
beta[beta < 0] <- 0
beta[beta > 1] <- 1

# --------------------------------------------------------------
# 2. deconvolution: EpiDISH RPC with the 450K DHS reference
# --------------------------------------------------------------
data(centDHSbloodDMC.m)
ref_mat <- centDHSbloodDMC.m          # 333 CpG x 7 cell types
cell_names <- colnames(ref_mat)
cat(sprintf("\n== RPC deconvolution: reference %d CpG x %d cells (%s) ==\n",
            nrow(ref_mat), ncol(ref_mat), paste(cell_names, collapse = ", ")))

common <- intersect(rownames(ref_mat), rownames(beta))
cat(sprintf("   probes shared with the array: %d/%d\n",
            length(common), nrow(ref_mat)))
stopifnot(length(common) >= 200)      # guard: too few -> results meaningless

est <- epidish(beta.m = beta[common, , drop = FALSE],
               ref.m  = ref_mat[common, , drop = FALSE],
               method = "RPC")
props <- est$estF                     # samples x 7, rows sum to 1
cat(sprintf("   estimated for %d samples; row sums in [%.3f, %.3f]\n",
            nrow(props), min(rowSums(props)), max(rowSums(props))))

cell_est <- data.frame(geo_accession = rownames(props),
                       props, check.names = FALSE)
write.csv(cell_est, file.path(OUT, "cell_estimates_GSE67530.csv"),
          row.names = FALSE)

# --------------------------------------------------------------
# 3. group comparisons per cell type
# --------------------------------------------------------------
long <- data.frame(
  geo_accession = rep(rownames(props), times = ncol(props)),
  cell_type     = rep(cell_names, each = nrow(props)),
  prop          = as.vector(props),
  stringsAsFactors = FALSE)
long$group <- pdata$group[match(long$geo_accession, pdata$geo_accession)]

grp_mean <- aggregate(prop ~ cell_type + group, long, mean)
grp_mean$prop <- round(grp_mean$prop, 4)
cat("\n== mean estimated proportions by group ==\n")
print(grp_mean, row.names = FALSE)

tests <- data.frame()
for (ct in cell_names) {
  d   <- long[long$cell_type == ct & !is.na(long$group), ]
  kw  <- kruskal.test(prop ~ group, data = d)
  for (pr in list(c("ARDS", "healthy"), c("ICU_control", "healthy"),
                  c("ARDS", "ICU_control"))) {
    a <- d$prop[d$group == pr[1]]; b <- d$prop[d$group == pr[2]]
    w <- suppressWarnings(wilcox.test(a, b))
    tests <- rbind(tests, data.frame(
      cell_type   = ct,
      comparison  = paste(pr, collapse = "_vs_"),
      mean_case   = round(mean(a), 4),
      mean_ctrl   = round(mean(b), 4),
      median_case = round(median(a), 4),
      median_ctrl = round(median(b), 4),
      wilcox_p    = signif(w$p.value, 3),
      kruskal_p   = signif(kw$p.value, 3),
      stringsAsFactors = FALSE))
  }
}
tests$p_bh <- with(tests, ave(p.adjust(wilcox_p, method = "BH"),
                              cell_type, FUN = identity))
write.csv(tests, file.path(OUT, "deconv_group_tests_GSE67530.csv"),
          row.names = FALSE)

cat("\n== headline: Neutro (neutrophil) proportions ==\n")
neu <- tests[tests$cell_type == "Neutro", ]
print(neu[, c("comparison", "mean_case", "mean_ctrl", "wilcox_p", "p_bh")],
      row.names = FALSE)

# --------------------------------------------------------------
# 4. boxplot
# --------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressMessages(library(ggplot2))
  p <- ggplot(long, aes(x = group, y = prop, fill = group)) +
    geom_boxplot(outlier.size = 0.7, alpha = 0.75) +
    facet_wrap(~ cell_type, scales = "free_y", nrow = 2) +
    scale_fill_manual(values = c("healthy" = "#7fb3d5",
                                 "ICU_control" = "#f5b041",
                                 "ARDS" = "#e74c3c"), guide = "none") +
    labs(title = "GSE67530: estimated blood cell proportions (EpiDISH RPC, DHS 450K reference)",
         x = NULL, y = "estimated proportion") +
    theme_classic(base_size = 11) +
    theme(strip.text = element_text(face = "bold"))
  ggsave(file.path(OUT, "fig_cellprops_GSE67530.png"), p,
         width = 10, height = 5, dpi = 300)
  cat("\n>> wrote fig_cellprops_GSE67530.png\n")
}

cat("\n>> cell_estimates_GSE67530.csv written.\n")
cat(">> NEXT: re-run R/02_ewas.R (single process!) -> it will fit the three\n")
cat(">> contrasts again with cell proportions as covariates (*_celladj tags).\n")
