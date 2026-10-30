# ==============================================================
# 02_ewas.R
# EWAS on GSE67530 (Zhang 2019 / MYLK study, PMID 27543902; 450K, GPL13534).
#
# Cohort (from GEO !Sample_characteristics_ch1, verified 2026-10-25):
#   ards: 1  -> 39 ARDS cases
#   ards: 0  -> 75 non-ARDS ICU controls
#   ards: NA -> 30 healthy controls        (n = 144; header shows 143 + 1 tail)
# Covariates available in the same header: age, male, black (race).
#
# Contrasts (each: limma M-value model, DMP table, 4-CpG AUC + DeLong):
#   [primary]  ARDS vs healthy      -- replicates Zhang 2019 (39 vs 30)
#   [ards_icu] ARDS vs ICU control  -- isolates ARDS-specific signal (novel)
#   [icu_hlth] ICU control vs healthy -- critical-illness effect (context)
#
# Replication anchors (Zhang 2019, PMID 31623626): cg03341377, cg24310395,
# cg07830557, cg08418670 (reported single-CpG AUC up to 0.99, in-sample,
# no external validation -> expect optimism; we report 95% CI + DeLong).
#
# OPTIONAL: if output/cell_estimates_GSE67530.csv exists (from
# 03_deconvolution.R), the primary contrast is refit with blood cell
# proportions as covariates (ARDS blood is neutrophil-shifted).
#
# Dependencies: limma, pROC (ggplot2 optional for volcanoes).
# ==============================================================

options(stringsAsFactors = FALSE)

# ---- project root: ALWAYS resolve to an absolute path ----
find_project_root <- function() {
  env <- Sys.getenv("ARDS_EPI_ROOT", "")
  if (nzchar(env)) return(normalizePath(env, mustWork = FALSE))
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  start <- if (length(f)) dirname(normalizePath(f, mustWork = FALSE)) else getwd()
  d <- normalizePath(start, mustWork = FALSE)
  while (!identical(d, dirname(d))) {
    if (file.exists(file.path(d, "R", "02_ewas.R")) || file.exists(file.path(d, "code", "02_ewas.R"))) return(d)
    d <- dirname(d)
  }
  d
}
PROJECT_ROOT <- find_project_root()
here <- function(...) file.path(PROJECT_ROOT, ...)
cat(sprintf("   PROJECT_ROOT (absolute): %s\n", PROJECT_ROOT))

suppressMessages({
  library(limma)
  library(pROC)
})

CONFIG <- list(
  beta_path  = here("data", "processed", "gse67530_beta.rds"),
  pdata_path = here("data", "processed", "gse67530_pdata.rds"),
  out_dir    = here("output"),
  signature_cpgs = c("cg03341377", "cg24310395", "cg07830557", "cg08418670"),
  # covariates present in the GEO characteristics (age / male / black)
  covariates = c("age", "male", "black")
)
dir.create(CONFIG$out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- inputs must exist (fail loudly with absolute paths) ----
missing <- c(CONFIG$beta_path, CONFIG$pdata_path)[
  !file.exists(c(CONFIG$beta_path, CONFIG$pdata_path))]
if (length(missing)) {
  stop("Missing input file(s):\n  ", paste(missing, collapse = "\n  "),
       "\nRun source('R/01_download.R') first; GSE67530's series matrix is",
       "\n~568 MB (144 samples x 485,577 full-precision beta values), so use",
       "\nthe curl.exe one-liner in README section 9 to fetch it.")
}

beta  <- readRDS(CONFIG$beta_path)
pdata <- readRDS(CONFIG$pdata_path)
cat(sprintf("   beta: %d probes x %d samples | pdata: %d rows x %d cols\n",
            nrow(beta), ncol(beta), nrow(pdata), ncol(pdata)))
if (!"geo_accession" %in% names(pdata)) pdata$geo_accession <- rownames(pdata)
if (!all(pdata$geo_accession %in% colnames(beta)))
  stop("pdata$geo_accession does not match beta column names.")
beta <- beta[, pdata$geo_accession, drop = FALSE]

# ==============================================================
# parse the repeated "key: value" characteristics lines into columns
# ==============================================================
parse_keyvalue_chars <- function(pd) {
  ch <- grep("characteristics", names(pd), ignore.case = TRUE, value = TRUE)
  out <- data.frame(geo_accession = pd$geo_accession, stringsAsFactors = FALSE)
  for (cc in ch) {
    v <- as.character(pd[[cc]])
    key <- sub(":.*$", "", v)
    val <- sub("^[^:]*:\\s*", "", v)
    uk <- unique(key)
    for (k in uk) {
      col <- make.unique(c(names(out), k))[length(names(out)) + 1L]
      out[[col]] <- ifelse(key == k, val, NA_character_)
    }
  }
  out
}
kv <- parse_keyvalue_chars(pdata)
cat("   parsed characteristic keys:", paste(setdiff(names(kv), "geo_accession"),
                                            collapse = ", "), "\n")

# ---- group from the 'ards' key: 1 = ARDS, 0 = ICU control, NA = healthy ----
if (!"ards" %in% names(kv))
  stop("No 'ards' characteristic found. Inspect kv and code the group manually.")
ards_raw <- trimws(kv$ards)
pdata$group <- factor(
  ifelse(ards_raw == "1", "ARDS",
  ifelse(ards_raw == "0", "ICU_control",
  ifelse(ards_raw %in% c("NA", "", "na"), "healthy", NA_character_))),
  levels = c("healthy", "ICU_control", "ARDS"))
cat("   group counts:\n"); print(table(pdata$group, useNA = "ifany"))

# numeric covariates
for (cv in CONFIG$covariates) {
  if (cv %in% names(kv)) {
    pdata[[cv]] <- suppressWarnings(as.numeric(trimws(kv[[cv]])))
    cat(sprintf("   covariate %s: n=%d non-missing, range %.1f-%.1f\n",
                cv, sum(!is.na(pdata[[cv]])),
                min(pdata[[cv]], na.rm = TRUE), max(pdata[[cv]], na.rm = TRUE)))
  }
}

# ---- M-value vs beta detection ----
# GEO 450K series matrices almost always DEPOSIT BETA (0..1). A small
# fraction of values can fall slightly outside [0,1] (quantile-normalisation
# / detection-p-value artefacts at the extremes) WITHOUT being M-values.
# Only convert when the bulk of the data is clearly on an M scale
# (median near 0 AND either a large fraction or a wide range outside [0,1]).
# Empirically GSE67530 here: median ~0.64, only 0.37% outside [0,1]
# (range -0.41..1.12) -> unmistakably BETA; clamp boundary noise, do NOT convert.
rng <- range(beta, na.rm = TRUE)
frac_out <- mean(beta < 0 | beta > 1, na.rm = TRUE)
med <- median(beta, na.rm = TRUE)
looks_like_M <- (frac_out > 0.05) || (rng[1] < -3 || rng[2] > 3)
if (looks_like_M) {
  cat(sprintf("   %.1f%% values outside [0,1] (range %.2f..%.2f) -> M-values, converting to beta\n",
             100 * frac_out, rng[1], rng[2]))
  beta <- 2^(beta) / (1 + 2^(beta))
} else {
  cat(sprintf("   %.2f%% outside [0,1] (range %.4f..%.4f, median %.3f) -> BETA; clamping boundary noise to [0,1]\n",
             100 * frac_out, rng[1], rng[2], med))
  beta[beta < 0] <- 0
  beta[beta > 1] <- 1
  cat("   (deposited metadata mentions M-values; here the distribution is\n")
  cat("    unmistakably beta: median ~0.6, dense near 0 and 1)\n")
}

# ---- M-values for limma ----
beta_m <- beta
beta_m[beta_m <= 0] <- 1e-4
beta_m[beta_m >= 1] <- 1 - 1e-4
M <- log2(beta_m / (1 - beta_m))

# ==============================================================
# one contrast: case vs control, with optional extra covariates
# ==============================================================
run_contrast <- function(case, ctrl, tag, extra_cov = NULL) {
  cat(sprintf("\n========== contrast: %s vs %s  [%s] ==========\n", case, ctrl, tag))
  keep <- pdata$group %in% c(case, ctrl) & !is.na(pdata$group)
  pd2 <- pdata[keep, , drop = FALSE]
  M2  <- M[, keep, drop = FALSE]
  y   <- as.integer(pd2$group == case)            # 1 = case, 0 = control
  cat(sprintf("   n = %d (%s %d / %s %d)\n", sum(keep), case, sum(y),
              ctrl, sum(!y)))

  # NOTE: lmFit does NOT add an intercept automatically. Always include an
  # explicit intercept column so that `case_vs_ctrl` is estimated as the
  # CASE-minus-CONTROL contrast. Without the intercept, `case_vs_ctrl`
  # would be tested against a model that has no baseline level, making the
  # t-statistic test "mean methylation == 0" -> ~85% of probes falsely
  # "significant" (this was a bug fixed 2026-10-26).
  X <- cbind(intercept = rep(1, nrow(pd2)), case_vs_ctrl = as.integer(y))
  used_cov <- c()
  for (cv in c(CONFIG$covariates, extra_cov)) {
    if (!is.null(cv) && cv %in% names(pd2) && sum(!is.na(pd2[[cv]])) > 5) {
      X <- cbind(X, pd2[[cv]])
      used_cov <- c(used_cov, cv)
    }
  }
  colnames(X) <- c("intercept", "case_vs_ctrl", used_cov)
  cat("   design:", paste(colnames(X), collapse = " + "), "\n")

  fit <- lmFit(M2, X)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = "case_vs_ctrl", number = Inf, sort.by = "P")
  tt$probe <- rownames(tt)
  tt <- tt[, c("probe", setdiff(names(tt), "probe"))]
  n_fdr  <- sum(tt$adj.P.Val < 0.05, na.rm = TRUE)
  n_bonf <- sum(tt$P.Value < 0.05 / nrow(tt), na.rm = TRUE)
  cat(sprintf("   DMPs: %d | FDR<0.05: %d | Bonferroni: %d\n",
              nrow(tt), n_fdr, n_bonf))
  write.csv(tt, file.path(CONFIG$out_dir, sprintf("ewas_dmps_%s.csv", tag)),
            row.names = FALSE)

  # ---- signature CpGs: effect sizes + AUC (single + combined) ----
  sig <- CONFIG$signature_cpgs[CONFIG$signature_cpgs %in% rownames(beta)]
  cat(sprintf("   signature CpGs on platform: %d/%d (%s)\n",
              length(sig), length(CONFIG$signature_cpgs),
              paste(sig, collapse = ", ")))
  if (!length(sig)) return(invisible(NULL))

  b2 <- beta[, keep, drop = FALSE]
  sig_tab <- data.frame(
    probe        = sig,
    mean_control = rowMeans(b2[sig, !y, drop = FALSE]),
    mean_case    = rowMeans(b2[sig,  y, drop = FALSE]),
    stringsAsFactors = FALSE)
  sig_tab$delta_beta <- sig_tab$mean_case - sig_tab$mean_control
  sig_tab <- merge(sig_tab,
                   tt[tt$probe %in% sig,
                      c("probe", "logFC", "AveExpr", "P.Value", "adj.P.Val")],
                   by = "probe", all.x = TRUE)
  write.csv(sig_tab, file.path(CONFIG$out_dir,
                               sprintf("ewas_signature_cpgs_%s.csv", tag)),
            row.names = FALSE)
  cat("   signature CpG effects:\n"); print(sig_tab)

  yy <- as.numeric(y)
  roc_one <- function(x) {
    r <- roc(yy, x, quiet = TRUE)
    a <- as.numeric(auc(r)); ci <- as.numeric(ci(r, of = "auc"))
    list(roc = r, auc = a, lo = ci[1], hi = ci[3])
  }
  roc_rows <- data.frame(probe = sig, AUC = NA_real_, CI_lo = NA_real_,
                         CI_hi = NA_real_, direction = NA_character_)
  rocs <- list()
  for (i in seq_along(sig)) {
    rr <- roc_one(as.numeric(b2[sig[i], ]))
    rocs[[sig[i]]] <- rr$roc
    auc_rep <- max(rr$auc, 1 - rr$auc)   # report informative direction
    roc_rows$AUC[i] <- auc_rep
    roc_rows$CI_lo[i] <- if (rr$auc >= 0.5) rr$lo else 1 - rr$hi
    roc_rows$CI_hi[i] <- if (rr$auc >= 0.5) rr$hi else 1 - rr$lo
    roc_rows$direction[i] <- if (rr$auc >= 0.5) "hyper-in-case" else "hypo-in-case"
  }
  # combined score: mean of direction-signed beta (direction taken from this
  # same data -> in-sample optimistic; flagged in the output header)
  sgn <- ifelse(sig_tab$delta_beta[match(sig, sig_tab$probe)] >= 0, 1, -1)
  comb <- colMeans(b2[sig, , drop = FALSE] * sgn)
  rcomb <- roc(yy, comb, quiet = TRUE)
  ac <- as.numeric(auc(rcomb)); ci_c <- as.numeric(ci(rcomb, of = "auc"))
  roc_rows <- rbind(roc_rows, data.frame(
    probe = "combined_4cpgs_signed",
    AUC = max(ac, 1 - ac),
    CI_lo = if (ac >= 0.5) ci_c[1] else 1 - ci_c[3],
    CI_hi = if (ac >= 0.5) ci_c[3] else 1 - ci_c[1],
    direction = if (ac >= 0.5) "hyper-in-case" else "hypo-in-case"))
  write.csv(roc_rows, file.path(CONFIG$out_dir,
                                sprintf("ewas_signature_auc_%s.csv", tag)),
            row.names = FALSE)
  cat("   AUC:\n"); print(roc_rows)

  # ---- DeLong: combined vs each single CpG ----
  delong <- data.frame(comparison = character(), AUC_single = numeric(),
                       AUC_combined = numeric(), DeLong_p = numeric())
  for (i in seq_along(sig)) {
    pt <- tryCatch(roc.test(rocs[[sig[i]]], rcomb, paired = TRUE),
                   error = function(e) NULL)
    delong <- rbind(delong, data.frame(
      comparison = paste0("combined vs ", sig[i]),
      AUC_single = roc_rows$AUC[i],
      AUC_combined = max(ac, 1 - ac),
      DeLong_p = if (is.null(pt)) NA_real_ else as.numeric(pt$p.value)))
  }
  write.csv(delong, file.path(CONFIG$out_dir,
                              sprintf("ewas_delong_%s.csv", tag)),
            row.names = FALSE)
  cat("   DeLong (combined vs single):\n"); print(delong)

  # ---- volcano ----
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    suppressMessages(library(ggplot2))
    p <- ggplot(tt, aes(x = logFC, y = -log10(P.Value))) +
      geom_point(aes(color = adj.P.Val < 0.05), size = 0.6, alpha = 0.5) +
      scale_color_manual(values = c("grey40", "red"), guide = "none") +
      geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed",
                 color = "grey60") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                 color = "grey60") +
      labs(title = sprintf("GSE67530: %s vs %s (M-value EWAS)", case, ctrl),
           x = "log2 fold-change (M scale)", y = "-log10 P") +
      theme_classic(base_size = 12)
    ggsave(file.path(CONFIG$out_dir, sprintf("fig_volcano_%s.png", tag)), p,
           width = 6, height = 5, dpi = 300)
  }
  invisible(list(tt = tt, roc = roc_rows, delong = delong))
}

# ==============================================================
# MAIN
# ==============================================================
res_primary <- run_contrast("ARDS", "healthy", "primary_paper")
res_icu     <- run_contrast("ARDS", "ICU_control", "ards_vs_icu")
res_ctx     <- run_contrast("ICU_control", "healthy", "icu_vs_healthy")

# ---- OPTIONAL: cell-composition-adjusted primary contrast ----
cell_file <- file.path(CONFIG$out_dir, "cell_estimates_GSE67530.csv")
if (file.exists(cell_file)) {
  cat("\n>> cell estimates found -> adjusted primary contrast ...\n")
  cells <- read.csv(cell_file, row.names = 1, check.names = FALSE)
  if ("geo_accession" %in% names(cells)) rownames(cells) <- cells$geo_accession
  common <- intersect(rownames(cells), pdata$geo_accession)
  if (length(common) >= 10L) {
    pdata$geo_accession <- as.character(pdata$geo_accession)
    cells <- cells[as.character(common), , drop = FALSE]
    # attach cell proportions as extra covariates via a temporary pdata column
    cov_names <- setdiff(names(cells), c("geo_accession", names(pdata)))
    # 7 proportions sum to 1 -> collinear with the intercept; also Eosino is
    # ~0 for every sample (near-constant). Drop the last column as reference.
    if (length(cov_names) >= 2) {
      dropped <- cov_names[length(cov_names)]
      cov_names <- cov_names[-length(cov_names)]
      cat(sprintf("   cell covariates used: %s (reference/dropped: %s)\n",
                  paste(cov_names, collapse = ", "), dropped))
    }
    for (cn in cov_names) {
      pdata[[cn]] <- NA_real_
      pdata[[cn]][match(as.character(common), pdata$geo_accession)] <- cells[[cn]]
    }
    run_contrast("ARDS", "healthy", "primary_paper_celladj", extra_cov = cov_names)
    run_contrast("ARDS", "ICU_control", "ards_vs_icu_celladj", extra_cov = cov_names)
    run_contrast("ICU_control", "healthy", "icu_vs_healthy_celladj", extra_cov = cov_names)
  }
} else {
  cat("\n   (no cell_estimates_GSE67530.csv yet -> unadjusted only;\n")
  cat("    run 03_deconvolution.R first, then re-run this script)\n")
}

cat("\n>> EWAS done. Outputs (absolute paths) in:", CONFIG$out_dir, "\n")
cat(">> Next: source('R/03_deconvolution.R') (needs EpiDISH + FlowSorted refs)\n")
