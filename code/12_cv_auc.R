# Round 4 reviewer response: repeated 5-fold CV AUC of the four-CpG score
#
# Addresses the criticism that all AUCs are in-sample estimates.
# Design: 5-fold stratified CV, 50 repetitions; the score direction is
# re-determined within each training fold (no leakage), score = sum of
# direction-signed beta values of the four signature CpGs.
# Run: Rscript code/12_cv_auc.R   -> results/cv_auc_4cpg.csv

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
setwd(ROOT)
dir.create("results", showWarnings = FALSE)

suppressMessages(library(data.table))

cpgs <- c("cg03341377", "cg24310395", "cg07830557", "cg08418670")
sm_path <- c(file.path(ROOT, "data", "raw", "GSE67530_series_matrix.txt.gz"),
             Sys.getenv("GSE67530_SERIES_MATRIX", ""),
             file.path(dirname(ROOT), "ards_epigenetics", "data", "raw",
                       "GSE67530_series_matrix.txt.gz"))
sm_path <- sm_path[file.exists(sm_path)][1]
stopifnot(!is.na(sm_path))

# metadata: stream only the header lines (cheap)
ch <- character(0)
con <- gzfile(sm_path, "rt")
repeat {
  l <- readLines(con, n = 1)
  if (length(l) == 0 || startsWith(l, "!series_matrix_table_begin")) break
  if (startsWith(l, "!Sample_characteristics_ch1")) ch <- c(ch, l)
}
close(con)
ards_line <- grep("ards", ch, value = TRUE)[1]
ards <- gsub('"', "", strsplit(ards_line, "\t")[[1]][-1])
ards <- sub("^.*ards:\\s*", "", ards)

# expression table: stream in chunks, keep only the four CpG rows
hits <- list(); k <- 0
con <- gzfile(sm_path, "rt"); on.exit(close(con))
in_tab <- FALSE
repeat {
  lines <- readLines(con, n = 50000L)
  if (!length(lines)) break
  if (!in_tab) {
    i <- which(startsWith(lines, "!series_matrix_table_begin"))
    if (!length(i)) next
    lines <- lines[(i[1] + 1):length(lines)]
    in_tab <- TRUE
  }
  hit <- lines[startsWith(lines, '"') | startsWith(lines, "ID_REF")]
  hit <- hit[substr(hit, 2, 11) %in% cpgs | startsWith(hit, '"ID_REF')]
  if (length(hit)) { k <- k + 1; hits[[k]] <- hit }
  if (any(startsWith(lines, "!series_matrix_table_end"))) break
}
tab <- unlist(hits)
header <- strsplit(gsub('"', "", tab[startsWith(tab, '"ID_REF')])[1], "\t")[[1]]
stopifnot(identical(header[1], "ID_REF"))
rows <- strsplit(gsub('"', "", tab[!startsWith(tab, '"ID_REF')]), "\t")
names(rows) <- vapply(rows, function(r) r[1], "")
stopifnot(setequal(names(rows), cpgs))
mat <- do.call(rbind, lapply(rows, function(r) as.numeric(r[-1])))
beta <- t(mat)  # samples x 4 CpGs
stopifnot(nrow(beta) == length(ards))

auc_rank <- function(y, s) {
  r <- rank(s)
  n1 <- sum(y); n0 <- length(y) - n1
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

run_cv <- function(pos, neg, label) {
  keep <- grp %in% c(pos, neg)
  y <- as.integer(grp[keep] == pos)
  X <- beta[keep, , drop = FALSE]
  sgn <- sign(colMeans(X[y == 1, , drop = FALSE]) - colMeans(X[y == 0, , drop = FALSE]))
  in_sample <- auc_rank(y, as.vector(X %*% sgn))
  set.seed(42)
  aucs <- replicate(50, {
    fold <- rep(0, length(y))
    for (lv in 0:1) {
      ii <- which(y == lv)
      fold[ii[sample(length(ii))]] <- seq_along(ii) %% 5
    }
    s <- numeric(length(y))
    for (f in 0:4) {
      tr <- fold != f; te <- fold == f
      sg <- sign(colMeans(X[tr & y == 1, , drop = FALSE]) -
                 colMeans(X[tr & y == 0, , drop = FALSE]))
      sg[sg == 0] <- 1
      s[te] <- as.vector(X[te, , drop = FALSE] %*% sg)
    }
    auc_rank(y, s)
  })
  data.frame(contrast = label, n = length(y),
             auc_in_sample = round(in_sample, 3),
             auc_cv_mean = round(mean(aucs), 3),
             auc_cv_p2.5 = round(quantile(aucs, 0.025), 3),
             auc_cv_p97.5 = round(quantile(aucs, 0.975), 3))
}

grp <- ifelse(ards == "1", "ARDS", ifelse(ards == "0", "ICU", "healthy"))
out <- rbind(run_cv("ARDS", "healthy", "ARDS_vs_healthy"),
             run_cv("ICU", "healthy", "ICU_vs_healthy"),
             run_cv("ARDS", "ICU", "ARDS_vs_ICU"))
write.csv(out, file.path(ROOT, "results", "cv_auc_4cpg.csv"), row.names = FALSE)
print(out)
