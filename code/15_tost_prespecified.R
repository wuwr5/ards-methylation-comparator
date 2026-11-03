# ==============================================================
# 15_tost_prespecified.R  (round-5 reviewer response)
#
# Reviewer point: the +/-0.05 equivalence margin was derived from the
# observed data (half of the smallest case-vs-healthy delta-beta), so the
# equivalence claim is circular. Fix: (i) state the margin a priori from
# external, data-independent grounds; (ii) make the inference independent
# of the exact margin by reporting delta_min, the smallest margin at which
# TOST would still be passed.
#
# Pre-specified primary margin: delta = 0.05 (beta units).
# External anchors (both independent of this study's effect sizes):
#   (a) EWAS convention: |delta-beta| < 0.05 is conventionally treated as a
#       small / not biologically interpretable difference on 450K beta scale
#       (the same threshold used to grade effect magnitude in EWAS reports);
#   (b) array precision: Illumina 450K technical duplicates show within-pair
#       beta differences of roughly that order for CpGs with intermediate
#       methylation (Xu & Taylor, Epigenetics 2021, PMID 32749174: 128
#       technical duplicate pairs; whole-array duplicate R2 ~0.99 but
#       per-CpG reproducibility is far lower), i.e. a difference of ~0.05
#       is at the edge of what a single array pair can resolve reliably.
#
# delta_min is defined as the smallest margin for which the one-sided 5%
# TOST passes: delta > |d| + t(0.95, df) * se, i.e. delta_min = that value.
# Reporting delta_min removes the arbitrary constant from the inference:
# "differences larger than delta_min can be excluded" holds for any
# pre-specified margin >= delta_min, whatever its provenance.
#
# Run: Rscript code/15_tost_prespecified.R
# Out: results/tost_prespecified.csv, results/tost_margin_grid.csv
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
RES <- file.path(REPO_DIR, "results")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

sm_path <- c(file.path(ROOT, "data", "raw", "GSE67530_series_matrix.txt.gz"),
             file.path(dirname(ROOT), "ards_epigenetics", "data", "raw",
                       "GSE67530_series_matrix.txt.gz"))
sm_path <- sm_path[file.exists(sm_path)][1]
stopifnot(!is.na(sm_path))

DELTA0 <- 0.05
GRID <- c(0.01, 0.02, 0.025, 0.03, 0.04, 0.045, 0.05, 0.06, 0.08, 0.10)
SIG <- c("cg03341377", "cg24310395", "cg07830557", "cg08418670")

# ---- stream header (characteristics) then the 4 probe rows only ----
ch <- character(0); headers <- list()
con <- gzfile(sm_path, "rt")
repeat {
  l <- readLines(con, n = 1)
  if (length(l) == 0 || startsWith(l, "!series_matrix_table_begin")) break
  if (startsWith(l, "!Sample_characteristics_ch1")) ch <- c(ch, l)
  if (startsWith(l, "!Sample_geo_accession")) headers[["gsm"]] <- strsplit(l, "\t")[[1]][-1]
}
close(con)
ch <- gsub('"', "", ch)
key <- trimws(sub(":.*$", "", strsplit(ch[1], "\t")[[1]][2]))
chv <- lapply(ch, function(x) strsplit(x, "\t")[[1]][-1])
names(chv) <- sapply(chv, function(v) trimws(sub(":.*$", "", v[1])))
chv <- lapply(chv, function(v) trimws(sub("^[^:]*:\\s*", "", v)))

# probe rows via chunked scan
rows <- list()
con <- gzfile(sm_path, "rt")
repeat {
  l <- readLines(con, n = 2000)
  if (!length(l)) break
  id <- gsub('"', "", substr(l, 1, 12))          # '"cg03341377"' -> 'cg03341377'
  hit <- l[id %in% SIG]
  for (h in hit) {
    p <- strsplit(h, "\t")[[1]]
    pid <- gsub('"', "", p[1])
    if (pid %in% SIG) rows[[pid]] <- suppressWarnings(as.numeric(gsub('"', "", p[-1])))
  }
  if (length(rows) == length(SIG)) break
}
close(con)
stopifnot(all(SIG %in% names(rows)))

ards <- chv$ards
grp <- ifelse(ards == "1", "ARDS", ifelse(ards == "0", "ICU_control", "healthy"))
age <- suppressWarnings(as.numeric(chv$age))
male <- suppressWarnings(as.numeric(chv$male))
black <- if (!is.null(chv$black)) suppressWarnings(as.numeric(chv$black)) else
  as.integer(tolower(chv$race) %in% c("black", "black or african american"))
cat("n:", length(grp), "| groups:", paste(names(table(grp)), table(grp), collapse = " "), "\n")
cat("covariates: age", sum(!is.na(age)), "male", sum(!is.na(male)), "black", sum(!is.na(black)), "\n")
cat("characteristic keys:", paste(names(chv), collapse = ", "), "\n")

tost_row <- function(x, g, covars, delta) {
  d <- data.frame(x = x, g = g, covars)
  d <- d[stats::complete.cases(d), ]
  fit <- lm(x ~ ., data = d)
  co <- summary(fit)$coefficients
  est <- co["g", 1]; se <- co["g", 2]; df <- fit$df.residual
  tcrit <- qt(0.95, df)
  # TOST: H0 |delta_true| >= margin
  p_lo <- 1 - pt((delta - abs(est)) / se, df)      # upper bound of |effect|
  tost_p <- p_lo
  data.frame(delta_hat = est, se = se, df = df,
             ci90_lo = est - tcrit * se, ci90_hi = est + tcrit * se,
             delta_min = abs(est) + tcrit * se,
             tost_p = tost_p, pass = tost_p < 0.05)
}

out <- list(); grid_out <- list()
for (cx in list(c("ARDS", "ICU_control"), c("ARDS", "healthy"))) {
  case <- cx[1]; ctrl <- cx[2]
  keep <- grp %in% c(case, ctrl)
  g <- as.integer(grp[keep] == case)
  covars <- data.frame(age = age[keep], male = male[keep], black = black[keep])
  for (p in SIG) {
    x <- rows[[p]][keep]
    r <- tost_row(x, g, covars, DELTA0)
    r$contrast <- paste0(case, "_vs_", ctrl); r$probe <- p
    r$n_case <- sum(g == 1); r$n_ctrl <- sum(g == 0)
    out[[paste(cx[1], cx[2], p)]] <- r
    for (dl in GRID) {
      rr <- tost_row(x, g, covars, dl)
      grid_out[[paste(cx[1], cx[2], p, dl)]] <- data.frame(
        contrast = r$contrast, probe = p, delta = dl,
        delta_hat = rr$delta_hat, tost_p = rr$tost_p, pass = rr$pass)
    }
  }
}
O <- do.call(rbind, out); rownames(O) <- NULL
G <- do.call(rbind, grid_out); rownames(G) <- NULL
O <- O[, c("contrast", "probe", "n_case", "n_ctrl", "delta_hat", "se", "df",
           "ci90_lo", "ci90_hi", "delta_min", "tost_p", "pass")]
write.csv(O, file.path(RES, "tost_prespecified.csv"), row.names = FALSE)
write.csv(G, file.path(RES, "tost_margin_grid.csv"), row.names = FALSE)

cat("\n=== pre-specified margin delta =", DELTA0, "===\n"); print(O, digits = 4)
cat("\ndelta_min per contrast (max over the 4 CpGs) = the smallest margin at which ALL four pass:\n")
print(tapply(O$delta_min, O$contrast, max), digits = 4)
cat("\n=== grid pass counts ===\n")
print(with(G, tapply(pass, list(contrast, delta), sum)))
cat("\nDONE 15\n")
